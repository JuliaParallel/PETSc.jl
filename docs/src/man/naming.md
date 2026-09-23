# Naming and API conventions

These rules govern the **high-level interface**: everything reachable as `PETSc.foo`.
The low-level layer (`LibPETSc.*`) mirrors the PETSc C API name for name and is out of scope; it keeps its C names so that PETSc's own documentation stays usable.

The conventions take effect in v0.5 and break the existing API. Renamed functions keep a deprecation shim for one minor cycle.
A few changes are semantic and cannot be shimmed; those are listed in [§16](#16.-Breaking-changes-without-a-shim).

## 1. Scope and the tiebreak

| Layer | Namespace | Naming |
|---|---|---|
| Low level | `PETSc.LibPETSc` | Verbatim PETSc C names (`DMStagGetCorners`) |
| High level | `PETSc` | These conventions |

Every function outside `LibPETSc` follows these rules, including thin wrappers around a single C call.

Two goals pull against each other here: being idiomatic Julia, and staying recognisable to someone reading the PETSc manual. When they conflict:

> **Julia idiom decides structure. PETSc decides vocabulary.**

Dispatch, `!`, argument order, and extending `Base` follow Julia. Which *word* to use follows PETSc.
So `DMDAGetDof` becomes `ndofs(dm)`: Julian in shape, PETSc's term for the thing.

### 1.1 What these rules cover

The high-level surface is `names(PETSc; all = true)`, filtered to bindings whose methods are defined in PETSc's own source. `scripts/api_surface.jl` enumerates it, and is the only place a count of it should be read from: the number moves with every additive PR, and a figure written into prose here is wrong by the time anyone reads it.

Three subsets are treated differently:

**Internal helpers** are exempt: they may be renamed freely, with no shim.

A leading underscore stays legitimate for one thing: distinguishing an inner worker from the wrapper that shares its name, as `Base` does with `_growend!`. As of now, `_mul!`, `_unsafe_local_array`, `_local_arrays` and `_restore_local_arrays!` are the cases in this package.
That is disambiguation, not a visibility marker.

Because [§13](#13.-Exports) exports only types, "unexported" cannot by itself separate the public API from internals.
The register is `scripts/renames.jl`, a plain list of `old => new` pairs plus a set of names marked internal. Every binding must appear in one or the other, and `scripts/api_surface.jl --check` fails if a binding appears in neither.

The register is data rather than prose because five things have to agree with it: `src/deprecations.jl`, the `public` declaration ([§13.1](#13.1-The-public-keyword)), the rename table at the end of this document, the creator and destroyer name sets in `audit.jl`, and the `--check` comparison set. All five are generated from it. Kept by hand they drift, and the auditor's drift is silent: after a rename its patterns stop matching and it reports no leaks on leaking code.

`scripts/api_surface.jl --sweeps` covers the other half: every count this document states is derived rather than remembered.

It derives four lists: functions taking `petsclib` first alongside a dispatchable object ([§8](#8.-Argument-order)), functions returning a `NamedTuple` ([§12](#12.-Return-values)), exported names against §13's list, and readers returning a PETSc object without a `doc_borrowed` entry ([§3.3](#3.3-What-an-accessor-hands-back)). Done by hand the first three were all wrong, by one, by six, and by nine respectively, and each miscount survived a full draft.

**Macros** follow the same rules as functions: snake_case, and no prefix that repeats the module name. `PETSc.@petsc_residual_fn` says "petsc" twice.

```julia
@residual_fn      # was @petsc_residual_fn
@jacobian_fn      # was @petsc_jacobian_fn
@bd_fn            # was @petsc_bd_fn
@simple_fn        # was @petsc_simple_fn
```

**Base extensions** are not renamed. `Base.size`, `Base.getindex` and the rest keep their Base names by definition; [§11](#11.-Extending-Base-and-stdlib) governs which ones may exist at all.

## 2. Case

Default to `snake_case` for any multi-word name. Run words together only for names on the closed list below.

```julia
set_type!(ksp, :gmres)
ghost_corners(dm)
ownership_range(A)
local_indices(dm)
local_coordinate_array(dm)
project_function!(vec, dm, f)
```

**Closed run-together list**:

- Base-style predicates: `issimplex`, `issymmetric`, `ishermitian`, `isassembled`,
  `isinitialized`, `isfinalized`
- PETSc single tokens: `ndofs`, `nextra`, `l2diff`, `nullspace`, `seqaij`, `seqdense`
- Type accessors mirroring `Base.eltype`: `scalartype`, `inttype`, `realtype`

The default requires no judgement, and the exceptions are countable. That matters more than matching `Base`, which is itself split: `setindex!` and `set_zero_subnormals` both ship in `Base`, and the Julia style guide only asks for underscores "as necessary to improve readability".

Never use `camelCase`. `PascalCase` means "this is a type", without exception (see [§5](#5.-Types)).

## 3. Accessors and setters

Readers are nouns. Drop PETSc's `Get`:

```julia
corners(dm)      # was getcorners
info(dm)         # was getinfo
dm(ksp)          # was getDM
solution(ksp)    # was get_solution
type_name(ksp)   # was gettype
comm(v)          # was getcomm
```

Writers are `set_*!`:

```julia
set_dm!(snes, dm)          # was setDM!
set_type!(ksp, :gmres)     # was KSPSetType
set_values!(A, …)          # was setvalues!
set_name!(obj, name)       # was petsc_setname!
```

This targets PETSc's C `Get`/`Set` prefixes. Base's `get(collection, key, default)` idiom is unaffected, and a function following it may keep `get` in its name.

The reader for `XXXGetType` is `type_name`, not `type`. The setter stays `set_type!`, because `set_type!(ksp, :gmres)` reads well and nothing shadows it. The asymmetry is deliberate and buys two things:

- `type` is the conventional name for the argument of a setter, so a bare `type` reader is shadowed by its own writers. It had already happened three times in `ts.jl` before the rule was written: `set_type!`, `set_adapt_type!` and `set_problem_type!` each take a parameter called `type`, in the file that defines `type` twenty lines above them.
- What comes back is the *name* of a PETSc implementation, not a Julia type. `type_name(ksp) === :gmres` says that; `type(ksp) === :gmres` invites the reader to expect `typeof`.

### 3.1 What accessors return and setters take

PETSc has two kinds of named value, and they are represented differently.

**Enumerations** are fixed sets, and the wrapped layer already exposes 187 of them as
Julia enums. Use the enum:

```julia
DMStag(petsclib, comm, (DM_BOUNDARY_PERIODIC, DM_BOUNDARY_NONE), …)
```

**Type names** (`KSPType`, `MatType`, `DMType`, `PCType`, …) are not enums and cannot become ones: PETSc registers them at runtime as strings, so plugins add new ones the package has never heard of.
They are `Symbol` at the Julia API and `String` at the C boundary:

```julia
set_type!(ksp, :gmres)
type_name(ksp)                             # :gmres, not "gmres"
PetscMat(petsclib, m, n; type = :seqaij)
DMStag(…; stencil_type = :box)
```

Symbol costs nothing to compare, needs no fixed set (so plugin types keep working), and matches `parse_options`, which already returns a `NamedTuple` with Symbol keys.
The conversion to `String` happens once, where the call meets C.

#### Considered: singleton types instead of Symbols

A named type could be a Julia value carrying the name in its type, `PetscTypeName{:KSP, :gmres}()`, rather than a bare Symbol. The answer differs by direction.

**As a return value, no.** PETSc answers with a runtime string, so the parameter is only known at runtime and the return type infers as the `UnionAll`, making every downstream call dynamic. `Symbol` infers concretely. [§5.3](#5.3-DM-flavour-is-a-type,-not-a-string) already settles it: a property earns a type parameter only where a method dispatches on it, and nothing on the Julia side branches on `gmres` versus `cg`.

**As an argument, yes, and it is additive.** There the parameter is a literal, so inference is unaffected, and it buys tab completion, a typo caught as an `UndefVarError` at the call site, and a `MethodError` for `set_type!(ksp, KSPTypes.GAMG)` where `set_type!(ksp, :gamg)` compiles and fails at runtime. Adding it later is a new method, not a change to one, so v0.5 decides nothing. The `Symbol` overload stays either way, because plugins register names the package cannot know.

Scoped in `v05-progress.md` for after 0.5.0, alongside `PC`. One caution for whoever takes it: `LibPETSc` already binds `KSPType` and 59 names like it as C-side `Ptr{Cchar}` aliases, so the struct needs a different name. And `PETScDiffEq.jl`, the first external consumer of this API, had a free choice here and spelled these as plain `String`.

### 3.2 Readers whose names are also variable names

`dm`, `snes` and `ksp` are accessors and, at the same time, the conventional variable names for what they return.
Inside the package `dm = dm(ksp)` shadows the function and breaks every later call in that scope.
Users are unaffected, because [§13](#13.-Exports) exports none of them and a local binding cannot shadow a qualified `PETSc.dm`.

The names stay short, and the package pays for it with a convention rather than a longer spelling:

```julia
d = dm(ksp)            # not dm = dm(ksp)
s = snes(ts)           # as test/ts.jl already writes it
k = ksp(ts)
```

`type_name` ([§3](#3.-Accessors-and-setters)) is the one case where the reader was given a longer name instead, because the colliding word there is a *parameter* name on the matching setter rather than a local, and a convention cannot reach into a signature.

### 3.3 What an accessor hands back

A reader returning a PETSc object returns a **borrowed handle**: the object belongs to whatever it was asked of, no finalizer is attached, and calling `destroy!` on it is a user error that invalidates the owner's copy.

```julia
solution(ts)                # owned by ts
snes(ts), ksp(ts)           # owned by ts
dm(ksp)                     # owned by ksp
tolerances(ts).vatol        # owned by ts
```

Constructors ([§6](#6.-Constructors)) return owned objects, with a finalizer.

**The rule is enforced by the value, not only by the docstring.** `VecPtr(petsclib, ptr, own::Bool)` already carries an ownership flag and attaches a finalizer only when `own` is true. v0.5 puts that field on every high-level wrapper, has `destroy!` consult it, and has readers return `own = false` where constructors return `own = true`.

Part of that is a bug fix. v0.4's `destroy(m::AbstractPetscVec)` consults only `isdestroyable`, which tests the finalize state, the null pointer and the object's age, never `own`, so the no-op `VecPtr`'s docstring promises does not happen. `PETScDiffEq.jl` builds `VecPtr(pl, x_ptr, false)` in every callback on the strength of that sentence.

The docstring says so too, through a helper beside `doc_external`, so CI greps for a `doc_borrowed` call rather than for prose. A reader returning a PETSc object with no such entry fails, unless it is listed as owning its result: the constructors, and `narrow`.

`narrow` ([§5.4](#5.4-DMs-of-unknown-provenance)) is the same rule from the other side: it returns a second handle onto one PETSc object, and destroying either invalidates the other.

## 4. Object prefixes

Drop `dm_`, `mat_`, `vec_`, `plex_` when the first argument already carries the type.
Dispatch does the work:

```julia
project_function!(vec, dm, f)   # was dm_project_function!
set_nullspace!(A, ns)           # was mat_set_null_space!
distribute!(dm)                 # was plexdistribute!
```

Keep a prefix only for free functions, where no argument identifies the domain.
These take `petsclib` and `comm` rather than a PETSc object:

```julia
fe_create_default(petsclib, comm, dim, Nc, …)
fe_create_lagrange(petsclib, comm, dim, Nc, …)
mat_nullspace_create(petsclib, comm; has_const = true)
```

Check this against the argument list, not against the name. `vtk_save!` looks like a free function and is not one: it takes a `PetscVec` ([§8](#8.-Argument-order)).

### 4.1 When the subject is not the first argument

"Drop the prefix, dispatch does the work" assumes the first argument names the domain. Two rules in this document break that assumption, and where they do, the prefix comes back:

> **Keep a prefix when the object that identifies the operation is not the first argument**, either because a callback occupies that slot ([§8.1](#8.1-Callbacks-come-first)) or because the target is a sub-object reached inside the call.

```julia
set_snes_jacobian!(updateJ!, snes, J, PJ)   # callback first, so snes is not argument 1
set_adapt_type!(ts, :basic)                 # target is the TSAdapt inside, not the TS
snes_iterations(ts)                         # the count belongs to the sub-solver
ksp_iterations(ts)
snes_failures(ts)
```

For the sub-object cases the alternative is a reader and no prefix, `set_type!(adapt(ts), :basic)`. That is rejected for now because it makes a public `TSAdapt` for one setter. `PC` is the first sub-object that earned a type of its own, and it takes the reader form: `set_type!(pc(ksp), :jacobi)`, not `set_pc_type!(ksp, :jacobi)`. `set_adapt_type!` keeps its prefix until `TSAdapt` has a type.

One clarification, because this rule was previously justified on a collision that does not exist. `set_snes_jacobian!(updateJ!, snes, J, PJ)` and `set_jacobian!(ds, fieldI, fieldJ, g0, g1, g2, g3)` do not collide: one takes 3 or 4 arguments with an untyped first slot, the other 7 with a `PetscDS` first, and Julia tells them apart unaided. The prefix is there because `snes` sits in argument 2 and cannot carry the name, and because merging a callback registration with seven function pointers behind one generic is the error this document exists to prevent.

## 5. Types

### 5.1 One meaning for PascalCase

Before v0.5, every high-level constructor was a PascalCase *function* whose name differed from the *struct* it returned:

```
KSP()      -> PetscKSP      VecSeq()    -> PetscVec
SNES()     -> PetscSNES     MatSeqAIJ() -> PetscMat
Options()  -> PetscOptions  DMDA()      -> PetscDM
DMStag()   -> PetscDM       DMPlex()    -> PetscDM
```

PascalCase meant either a type or a factory returning something else, and the two never agreed.
`DMStag()` returning a plain `PetscDM` is the direct cause of the string-based dispatch described in [§5.3](#5.3-DM-flavour-is-a-type,-not-a-string).

The rule is now:

> **PascalCase always names a type. Construction always goes through the type.**

Most of what follows is a consequence of that one sentence.

### 5.2 Prefixes

Use PETSc's bare class name. Add a `Petsc` prefix only when the bare name is too generic to export safely:

```julia
# prefixed
PetscVec  PetscMat  PetscOptions  PetscDS
PetscDM   # low-level handle only, see §5.4

# bare
DMDA  DMStag  DMPlex
KSP  SNES  TS  Tao  IS  AO  PF
```

This renames `PetscKSP` to `KSP` and `PetscSNES` to `SNES`, which also removes their constructor/struct name mismatch. `Vec` and `Mat` stay prefixed: they are exported, and bare `Vec`/`Mat` would clash across the ecosystem.

### 5.3 DM flavour is a type, not a string

Before v0.5 there was a single concrete `PetscDM{PetscLib}`, and flavour was resolved at runtime by comparing a string returned from C:

```julia
# v0.4: type-unstable, the two branches return different NamedTuples
function getcorners(dm)
    type = gettype(dm)
    if type == "da";       return getcorners_dmda(dm)
    elseif type == "stag"; return getcorners_dmstag(dm)
    end
end
```

That forced the `_dmda`/`_dmstag` suffixes and made the return type unpredictable.
In v0.5 flavour and dimension are type parameters, so ordinary dispatch applies and each method has one concrete return type:

```julia
abstract type AbstractPetscDM{PetscLib} end

DMDA{PetscLib,N}   <: AbstractPetscDM{PetscLib}
DMStag{PetscLib,N} <: AbstractPetscDM{PetscLib}
DMPlex{PetscLib}   <: AbstractPetscDM{PetscLib}

corners(dm::DMDA{L,N})   where {L,N} = …   # (lower, upper, size)
corners(dm::DMStag{L,N}) where {L,N} = …   # (lower, upper, size, nextra)
```

Dimension earns a parameter on the same test as flavour: only where a method dispatches on it or a return type is shaped by it. That holds for DMDA and DMStag, whose corners and creation paths are written per dimension. It does not hold for DMPlex, where nothing in `dmplex.jl` dispatches on dimension and neither constructor could supply one honestly, since `dim` arrives as a runtime argument and `DMPlex(petsclib, comm)` leaves it unset until setup. A plex reports its dimension through `ndims(dm)`, which is what it is: a runtime property of the mesh.

There are therefore **no type suffixes** on function names. `getcorners_dmstag` becomes a method of `corners`, not a separate function.

Flavour becomes a type parameter only when a high-level method must dispatch on it.
DM qualifies, because `corners`, `local_indices` and coordinate handling all differ by flavour.
Vec and Mat do not: their flavour only affects construction, so it stays a keyword ([§9](#9.-Keyword-versus-positional)).

### 5.4 DMs of unknown provenance

`PetscDM{PetscLib}` remains as the low-level handle, since `LibPETSc.DMCreate` and friends must return something before the flavour is known. It is not part of the high-level API.

A DM obtained from another object (`dm(ksp)`, a `DMPlex` read from file) has its flavour and dimension known only at runtime. Such boundary functions query `DMGetType` and `DMGetDimension` and return the concrete type:

```julia
d = dm(ksp)   # ::Union{DMDA{L,1},…,DMStag{L,3},…}
```

That Union has more than nine members, well past `MAX_UNION_SPLITTING`, so inference gives up and the call is a dynamic dispatch.
One dispatch is cheap. What is not cheap is propagating an abstractly-typed DM through a hot loop, so narrow once behind a function barrier:

```julia
function work(ksp)
    d = dm(ksp)     # dynamic dispatch here, once
    _work(d)        # specialized on the concrete type
end
```

`narrow(dm)` performs that query. Accessors returning a DM that PETSc owns (`dm` on a KSP or a SNES) call it themselves, so the flavour is in the type by the time it reaches the caller. It is also the escape hatch for a handle from a low-level creator.

`narrow` returns a second handle onto the same PETSc object, borrowed in the sense of [§3.3](#3.3-What-an-accessor-hands-back). Flavours with no type of their own come back unchanged.

| v0.4 | v0.5 |
|---|---|
| (new) | `narrow` |

### 5.5 Abstract, callback and wrapper types

**Abstract types** are `Abstract` followed by whatever [§5.2](#5.2-Prefixes) gives the concrete name, so the prefix decision is made once:

```julia
AbstractPetscVec  AbstractPetscMat  AbstractPetscDM      # concrete is prefixed
AbstractKSP  AbstractPC  AbstractSNES  AbstractTS  AbstractIS  AbstractAO
```

This renames `AbstractPetscKSP` and `AbstractPetscSNES`, matching `PetscKSP` becoming `KSP`.
It also fixes `AbstractPETScMemBackend`, the one type spelling the prefix `PETSc` where every other type spells it `Petsc`.
That one is exported, so the inconsistency is user-visible.

**Callback types** take a `Fn` suffix, following the autowrapped layer's existing `DMDATSRHSFunctionLocalFn`. The current `Fn_` prefix puts an underscore inside a PascalCase name, which no rule here permits:

```julia
KSPComputeRHSFn        # was Fn_KSPComputeRHS
KSPComputeOperatorsFn  # was Fn_KSPComputeOperators
SNESSetFunctionFn      # was Fn_SNESSetFunction
SNESSetJacobianFn      # was Fn_SNESSetJacobian
```

**Wrapper and alias types** keep bare PETSc names. [§5.2](#5.2-Prefixes) asks whether *this name* is too generic to export, not whether its class is: `Mat` needs the prefix, `MatShell` does not.

```julia
MatShell  MatOp  MatPtr  VecPtr        # unchanged
```

`MatAT` was a `Union{PetscMat, Transpose, Adjoint}` alias whose name did not say so. It had no uses, so v0.5 deletes it rather than renaming it to `MatOrTranspose`: an alias nothing refers to is not worth a shim.

## 6. Constructors

Construction goes through the type, dispatching on argument types.
This replaces the `MatXxx`/`VecXxx` factory family, and removes the confusion between `MatCreateSeqAIJ` (converted a `SparseMatrixCSC`) and `MatSeqAIJ` (allocated from sizes), which were unrelated functions with near-identical names:

```julia
A = PetscMat(petsclib, S::SparseMatrixCSC)        # was MatCreateSeqAIJ
A = PetscMat(petsclib, m, n, nnz)                 # was MatSeqAIJ
A = PetscMat(petsclib, m, n; type = :dense)       # was MatSeqDense
A = PetscMat(petsclib, rowptr, colval, nzval)     # was MatSeqAIJWithArrays

v = PetscVec(petsclib, 10)                        # was VecSeq
v = PetscVec(petsclib, x::Vector)                 # was VecSeq
```

`petsclib` appears in constructors because no object exists yet to carry it.
It appears nowhere else ([§8](#8.-Argument-order)).

## 7. Mutation

A function gets `!` if and only if it mutates one of its arguments, **or** package-global state.
Configuring an opaque PETSc object counts as mutation:

```julia
corners(dm)              # reads
solution(ksp)            # reads
ndofs(dm)                # reads

destroy!(dm)             # was destroy
setup!(A)
assemble!(A)
ghost_update!(v)
set_type!(ksp, :gmres)
solve!(x, ksp, b)

set_library!(path)       # mutates package state, no argument mutated
unset_library!()
```

The global-state clause exists so that `set_library!` and `unset_library!` keep their bang, following `Random.seed!`.
It does not apply to `set_petsclib`, which despite its name mutates nothing: it builds and returns a library handle, and becomes a `PetscLibType` constructor ([§6](#6.-Constructors)).

`destroy` becomes `destroy!`, and finalizer registrations change with it (`finalizer(destroy!, v)`).

## 8. Argument order

Subject first, and mutated arguments before read-only ones, matching `mul!(C, A, B)` and `copyto!(dest, src)`:

```julia
solve!(x, ksp, b)                  # x is written
project_function!(vec, dm, f)      # vec is written
global_to_local!(lvec, dm, gvec)   # lvec is written

corners(dm)                        # readers: subject first
assemble!(A)
```

`petsclib` never leads a high-level call. The object carries `PetscLib` as a type parameter, so passing it again is redundant.
The exceptions are constructors ([§6](#6.-Constructors)) and the free functions of [§4](#4.-Object-prefixes).

**Eight** functions in v0.4 take `petsclib` first while also taking a dispatchable PETSc object, so the first argument is recoverable from the second and is dropped:

```
add_boundary!       add_natural_boundary!   dm_compute_l2diff
dm_project_field!   dm_project_function!    plex_set_snes_local_fem!
vtk_save!           dmda_star_fd_coloring
```

This list, and the two others this document states as counts ([§12](#12.-Return-values)'s `NamedTuple` returns and [§13](#13.-Exports)'s export diff), are produced by `scripts/api_surface.jl --sweeps` rather than by reading the source. The first version of this section said seven and missed `dmda_star_fd_coloring`, which is what a sweep asserted by hand is worth.

`vtk_save!` is in that list, which corrects [§4](#4.-Object-prefixes): it takes a `PetscVec`, so it is not a free function and does not keep its prefix.
It becomes `save_vtk!(vec, filename)`, with `comm` recovered from the vector.
`vtk_save_fields!` takes an iterable of vectors and becomes `save_vtk!(vecs, filename)` on the same name.

### 8.1 Callbacks come first

A function-valued argument goes first, ahead of the subject, so `do` syntax works.
This outranks subject-first, and follows `map(f, c)` and `open(f, path)`:

```julia
set_function!(snes, x) do f, x
    …
end

with_local_array!(v) do arr
    arr .= 1
end
```

PETSc callbacks are where `do` earns its keep, and `setfunction!`, `setjacobian!` and `withlocalarray!` already take the callback first in v0.4.

## 9. Keyword versus positional

Positional arguments are the ones the object cannot exist without.
Keyword arguments are the ones with a sensible default:

```julia
DMStag(petsclib, comm, boundary, dims, dof;
       stencil_width = 1, stencil_type = :box)

PetscMat(petsclib, m, n, nnz)             # required
PetscMat(petsclib, m, n; type = :dense)   # optional, defaults to :seqaij

distribute!(dm; overlap = 0)
```

Never take a keyword argument that changes the return type. `type = :dense` is acceptable because every variant returns `PetscMat`; a keyword selecting a DM flavour would not be, which is why flavour is a type parameter there ([§5.3](#5.3-DM-flavour-is-a-type,-not-a-string)).

## 10. Abbreviations

Follow PETSc's own abbreviations, so a user reading the PETSc manual guesses the Julia name correctly:

```julia
ndofs(dm)        # DMDAGetDof
l2diff(dm, …)    # DMComputeL2Diff
nextra           # DMStagGetCorners
seqaij, dense, mpi, ds, fe, is, pc, ksp, snes, ts
```

Do not invent abbreviations PETSc does not use, and do not expand ones it does.

## 11. Extending Base and stdlib

Add methods to `Base` and `LinearAlgebra` wherever the PETSc operation satisfies the existing contract, so that generic Julia code works on PETSc objects:

```julia
Base.size, Base.length, Base.ndims, Base.eltype, Base.axes
Base.similar, Base.fill!, Base.iterate, Base.copyto!
Base.getindex, Base.setindex!, Base.show
LinearAlgebra.norm, LinearAlgebra.mul!,
LinearAlgebra.issymmetric, LinearAlgebra.ishermitian
```

### 11.1 Blocklist

Do **not** extend these. The word matches but the contract does not, and a silent mismatch is worse than an unfamiliar name:

| Name | Why not | Use instead |
|---|---|---|
| `LinearAlgebra.nullspace` | Returns a basis matrix; PETSc's `MatNullSpace` is a solver hint object | `set_nullspace!`, `mat_nullspace_create` |
| `Base.values` | PETSc "values" are a setter concept (`MatSetValues`), not a view of contents | `entries` |
| `Base.view` | `PetscViewer` writes an object to a stream; it is not an array view | `petscview`, or `Base.show` |
| `Base.setfield!` | Core builtin that writes a struct field; PETSc's attaches a finite element to a DM | `set_field!` |
| `Base.empty` | Returns an empty collection of the same type; a DM is not a collection, and what v0.4 returned was a null handle rather than an empty DM | the type constructor, e.g. `DMStag{PetscLib, N}(C_NULL, petsclib.age)` |

Additions need the same three columns: the name, the contract mismatch, and the replacement.

`setfield!` is worth a note. v0.4 defines `PETSc.setfield!` with 8 methods as a function distinct from the core builtin, so inside the module the builtin is shadowed.
§2's snake_case rule resolves it without a special case: `set_field!` and `setfield!` are different identifiers.

## 12. Return values

Functions returning several related values return a `NamedTuple`, with field names from a fixed vocabulary:

| Field | Type | Meaning |
|---|---|---|
| `lower`, `upper` | `CartesianIndex{N}` | 1-based, inclusive |
| `size` | `NTuple{N,Int}` | Local extent |
| `center`, `vertex` | `NamedTuple` of `UnitRange` | Staggered index ranges, keyed `x`, `y`, `z` |
| `x`, `y`, `z` | `UnitRange` | Per-axis range inside `center`/`vertex` |
| `nextra` | `NTuple{N,Int}` | DMStag partial elements |

```julia
c = corners(dm2d)
c.lower                      # CartesianIndex{2}
(; lower, upper) = corners(dm)
```

This vocabulary is closed and covers **index-shaped returns only**: corners, ghost corners, and local and global index ranges. Extend it rather than inventing a synonym: a local extent is `size`, never `dims` or `extent`.

Every other `NamedTuple` return takes its field names from the PETSc parameters it reports, spelled in `snake_case`, one field per value:

```julia
tolerances(ts)     # (; atol, rtol, vatol, vrtol)      TSGetTolerances
info(dm)           # (; dim, global_size, procs, ndofs,
                   #    stencil_width, boundary_type, stencil_type)
library_info()     # (; source, path, scalar, int, real)
```

`info` is the only existing case that breaks the rule twice over: v0.4 returns the stencil width as both `s` and `stencil_width`, so `s` goes, and `dof` becomes `ndofs` per [§10](#10.-Abbreviations).

Twelve functions return a `NamedTuple` today and `library_info` becomes the thirteenth. Four use the closed vocabulary: `corners` and `ghost_corners` on both DM flavours, `local_indices` and `global_indices`. The rest take PETSc's parameter names: `info`, `tolerances`, `star_fd_coloring`, `parse_options`, `audit_file` and `check_wrappers_version`. That list is generated by `scripts/api_surface.jl --sweeps` ([§8](#8.-Argument-order)); the first hand sweep of this section reported six.

`center` and `vertex` are themselves `NamedTuple`s keyed by axis, not `NTuple`s. Under dimension-correct returns a 2D DM yields `(x = …, y = …)` with no `z`:

```julia
local_indices(dm3d).center     # (x = 2:9, y = 2:9, z = 2:9)
local_indices(dm2d).center     # (x = 2:9, y = 2:9)
```

Results are **dimension-correct**: a 2D DM returns `CartesianIndex{2}` and 2-tuples.
This is a semantic break, covered in [§16](#16.-Breaking-changes-without-a-shim).

This is a performance fix as much as a break. v0.4's `getcorners_dmda` builds a heap `Vector` and splats it into `CartesianIndex`, which hides the length from the compiler, so `lower` and `upper` infer as `Any` and there are three allocations per call. So the new methods must build their tuples with `ntuple(..., Val(N))` and index the axis names by the type parameter, never by splatting a collection. Written that way the names constant-fold and a 2D `local_indices(dm).center` infers concretely as `@NamedTuple{x::UnitRange{Int}, y::UnitRange{Int}}`.

Two returns stay abstract by design and are not covered: `corners` has three fields on a `DMDA` and four on a `DMStag`, so code generic over `AbstractPetscDM` sees a small union, and `dm(ksp)` returns a wide one ([§5.4](#5.4-DMs-of-unknown-provenance)).

### 12.1 Index base

> An index the high-level layer accepts or returns **into Julia data** is 1-based. `LibPETSc` is 0-based.

The qualifier matters, and v0.4 has counterexamples in both directions. The invariant covers positions in an array, and index ranges over one: corners, ghost corners, local and global indices, ownership ranges, `dof_slot`. It does not cover identifiers PETSc assigns and expects back unchanged.

| Kind | Rule | Examples |
|---|---|---|
| Index into Julia data | 1-based | `corners`, `ghost_corners`, `local_indices`, `global_indices`, `ownership_range`, `dof_slot` |
| Identifier PETSc issues | PETSc's own numbering | `set_residual!` field number, `add_boundary!` label value, `set_field!` field number |

A field number is not a position in anything the user indexes; it is the number `set_field!` gave that field, handed back to PETSc. Renumbering it would mean translating in both directions for no gain, and would silently disagree with every PETSc example.

**`ownership_range`.** v0.4 takes `base_one::Bool = true` as a *positional* argument, the one function whose index convention is a runtime choice. The default is already 1-based, so callers who never passed it see no change.

```julia
ownership_range(A)          # 1-based, always
ownership_range(A, false)   # warns in 0.5, MethodError in 0.6
```

The shim keeps the positional form, because that is the form v0.4 has. A keyword shim would compile and never fire.

**Bulk index arrays keep PETSc's base.** The scalar and slice interface is 1-based, because `Base.getindex` and `Base.setindex!` have no choice. A `Vector` of indices handed to C, or handed back from C, keeps the base C uses and says so in its name:

| Function | Base | Naming |
|---|---|---|
| `set_values!` rows and columns (`MatSetValues`) | 0-based | `rows_0b`, `cols_0b` |
| `set_values!` rows and columns (`MatSetValuesStencil`) | 0-based | `rows_0b`, `cols_0b`, holding `MatStencil` |
| `star_fd_coloring` returned index vectors | mixed | `row_coo_local_0b`, `perturb_cols_1b`, … |
| `A[i, j] = v`, `A[I::CartesianIndex, J] = v` | 1-based | Base's contract, converts internally |

Two things forced this rather than a flip to 1-based everywhere. `set_values!`'s second overload takes `MatStencil` structs whose layout must match PETSc's C struct field for field, and reinterpreting the contents of a C-layout struct would mean rebuilding the user's vector on a hot path, or leaving the two overloads on different bases. And the 1-based route already exists: `Base.setindex!` converts, including the `CartesianIndex` form that builds `MatStencil` for the user, which leaves `set_values!` as the explicit bulk route for code already holding PETSc indices.

`star_fd_coloring` returns both kinds at once, its `row_coo_local` going straight to `LibPETSc.MatSetPreallocationCOOLocal` while `perturb_cols` and `local_rows` index Julia arrays. `examples/ex19.jl` already renames the latter to `_1b` on receipt, by hand, which is this convention arriving on its own.

So nothing here is a semantic break, and `setvalues!` to `set_values!` is a plain rename with a forwarding shim. A name carrying `_0b` cannot quietly mislead anyone.

## 13. Exports

Only types and construction entry points are exported. Verbs and accessors stay qualified:

```julia
export LibPETSc
export DMDA, DMStag, DMPlex
export PetscVec, PetscMat, PetscOptions
export KSP, SNES, TS
export petsclibs
```

`petsclibs` is exported because every entry point takes one, and `LibPETSc` because the low-level layer is reached through it. Neither is a verb.

```julia
using PETSc
dm = DMStag(petsclib, comm, …)   # exported
c  = PETSc.corners(dm)           # qualified
PETSc.solve!(x, ksp, b)          # qualified
```

`initialize` and `finalize` are **not** exported: `Base.finalize` already exists, and shadowing it would be worse than typing `PETSc.finalize`.

Two reasons for staying narrow. `PETSc.solve!` has no `CommonSolve` dependency, so exporting it would make `using PETSc, LinearSolve` ambiguous for a common combination.
And adding an export later is non-breaking while removing one is not, so a narrow list keeps the option open.

### 13.1 The `public` keyword

Julia 1.11 added `public`, which marks a name as API without exporting it.
That is the distinction this package needs, because exporting only types ([§13](#13.-Exports)) leaves "unexported" unable to separate the API from internals.

v0.5 adopts it unconditionally: `Project.toml` requires Julia 1.12, so a version gate would be dead code.

```julia
public @bd_fn, @jacobian_fn, …, with_local_array!, wrap_local_array
```

The declaration is `src/public_names.jl`, generated from `scripts/renames.jl` ([§1.1](#1.1-What-these-rules-cover)), so it cannot drift from the register. A name cannot be both exported and `public`, so the generator subtracts §13's list. `names(PETSc)` then reports the API directly and `scripts/api_surface.jl --check` compares against it, rather than depending on searching this document for a substring.

This matters more than it looks, because [§3](#3.-Accessors-and-setters) creates short accessors (`dm`, `comm`, `info`, `ds`, `label`, `solution`) that a substring search cannot verify at all.

## 14. Errors

Argument problems raise standard Julia exceptions. `@assert` is reserved for invariants that cannot fail unless the package itself is wrong, and must never validate user input: it carries no useful message and is not guaranteed to run.

| Condition | Exception |
|---|---|
| Wrong argument value or type | `ArgumentError` |
| Mismatched sizes or shapes | `DimensionMismatch` |
| Library not initialized | `PetscNotInitialized` |
| Error returned by PETSc itself | `PetscError` (existing) |

Every high-level constructor that creates a PETSc object (`PetscVec`, `PetscMat`, `PetscOptions`, `KSP`, `SNES`, `TS` and the DM types) checks that its library is initialized before calling PETSc, so a missing `initialize` is reported as `PetscNotInitialized` and not as the `PetscError` PETSc would raise. `test/test_errors.jl` checks each of them. `PetscOptions` is included although PETSc's `PetscOptionsCreate` works before `PetscInitialize`: the object records the library's current `age`, `initialize` advances it, so an options database created first would already be stale and `destroy!` would skip it.

```julia
length(A) == prod(sz) ||
    throw(DimensionMismatch("array has \$(length(A)) entries, expected \$(prod(sz))"))

T === petsclib.PetscScalar ||
    throw(ArgumentError("scalar type \$T does not match the library's \$(petsclib.PetscScalar)"))

isinitialized(lib) || throw(PetscNotInitialized(lib))
```

This rule was written against the 43 `@assert` the high-level layer carried before #250, which converted three of the four kinds. Ten remained after it, and one after v0.5:

| Kind | Was | Now | Becomes |
|---|---|---|---|
| Size or length mismatch | 14 | 0 | `DimensionMismatch`, done in #250 |
| Library not initialized | 11 | 0 | `PetscNotInitialized`, done in #250 |
| Scalar or integer type mismatch | 8 | 0 | `ArgumentError`, done in #250 |
| DM flavour check (`gettype(dm) == "da"`) | 9 | 0 | deleted in v0.5; dispatch enforces it ([§5.3](#5.3-DM-flavour-is-a-type,-not-a-string)) |
| Startup invariant (`found_ref[] == PETSC_TRUE`) | 1 | 1 | unchanged, and correct: it cannot fail unless the package is wrong |

The fourth row was the point, and the only one left for v0.5. Those nine checks existed only because one DM type had to police itself at runtime: seven in `dmstag.jl`, one in `dmda.jl`, one in `dm.jl`. Giving DM flavour a type deleted them rather than converting them, which is why they survived #250.

The fifth row is what the rule is for. `@assert` stays where it guards a package invariant and never where it validates user input.

## 15. Docstrings and discoverability

Renaming thin wrappers costs discoverability: a user who knows `PetscFECopyQuadrature` cannot grep for it once it is `copy_quadrature!`.
Every high-level docstring must therefore carry an external link to the C function it wraps, using the existing `doc_external` helper:

```julia
"""
    corners(dm::DMStag)

…

# External Links
$(doc_external("DMSTAG/DMStagGetCorners"))
"""
```

CI fails on a high-level docstring with no `doc_external` entry.
The C-name-to-Julia-name index page is generated from those entries, so the lookup table maintains itself.

### 15.1 Say what goes wrong

Where a failure mode is known, the docstring states it. A link to the C page says what a function does; it rarely says what it does when misused, and PETSc is a library where misuse frequently succeeds quietly.

`PETScDiffEq.jl` sets the standard here, and it is worth matching:

> A wrong Jacobian is not caught here. Where the other implicit families fail to converge, this one reports success and returns a wrong answer, so check a hand-written `jac` against a finite-difference solve before trusting it.

> A mass matrix is rejected. PETSc's coupled-stage matrix assumes `dF/du̇ = I`, and with a non-identity mass matrix the answer drifts further from the true one as `dt` shrinks instead of failing, which is worse than an error.

Three cases in this package already qualify and should be written up as the rename touches them:

- `ksp(ts)` raises `PETSC_ERR_ARG_WRONG` unless the problem was declared `TS_LINEAR`
- `snes(ts)` on an explicit method creates an unused solver rather than reporting anything
- `type_name(ts)` answers `nothing` until `solve!` has applied the options

This is a documentation rule, not a CI rule: "is a failure mode known" cannot be checked mechanically. It belongs here because the alternative is that the knowledge stays in a consumer package, which is where it is today.

## 16. Breaking changes without a shim

A shim translates names, not semantics. These changes have no shim and must be read before upgrading:

| Change | Symptom |
|---|---|
| Dimension-correct returns ([§12](#12.-Return-values)) | `corners(dm2d).size[3]` returned `1`, now throws `BoundsError`. `c.lower` is `CartesianIndex{2}`, not `{3}` |
| `DMStag`/`DMDA`/`DMPlex` return concrete types ([§5.3](#5.3-DM-flavour-is-a-type,-not-a-string)) | Code annotated `::PetscDM` no longer matches. Use `AbstractPetscDM` |
| `dm(ksp)` returns a Union ([§5.4](#5.4-DMs-of-unknown-provenance)) | Type-unstable at the boundary; add a function barrier in hot code |
| Type names are `Symbol` ([§3.1](#3.1-What-accessors-return-and-setters-take)) | `type_name(ksp) == "gmres"` is now false; compare against `:gmres`. Setters still accept a `String` only through the v0.5 shim |
| `@assert` replaced by typed exceptions ([§14](#14.-Errors)) | Code catching `AssertionError` must catch `ArgumentError`, `DimensionMismatch` or `PetscNotInitialized` |
| Arguments reordered ([§8](#8.-Argument-order)) | `dm_project_function!` and `dm_project_field!` took the written vector **last**; it moves to first. `dm_global_to_local!`/`dm_local_to_global!` took the DM last; it moves after the written vector. A shim can forward these, but any call written positionally against the old order and passed through `invoke` or a function reference will not be caught |
| Subject-first callback setters removed ([§8.1](#8.1-Callbacks-come-first)) | v0.4 accepts both `setfunction!(snes, f!, v)` and `setfunction!(f!, snes, v)`, and the same for `setjacobian!`. Only the callback-first order survives |
| Twelve exported names lose their export ([§13](#13.-Exports)) | v0.4 exports nine functions — `audit_petsc_file`, `determine_memtype`, `dmda_star_fd_coloring`, `get_petsc_arrays`, `library_info`, `restore_petsc_arrays`, `set_library!`, `set_petsclib`, `unset_library!` — and three memory-backend types, `AbstractPetscMemBackend`, `AbstractPETScMemBackend` and `HostBackend`, the last of which is exported but never defined. It exports no other type, and in particular none of the ones a user constructs. §13 replaces that list wholesale, keeping only `LibPETSc`. A shim cannot help: the replacement is not exported either, so `using PETSc` code must qualify the call or import the name |
| Borrowed handles ([§3.3](#3.3-What-an-accessor-hands-back)) | No behaviour changed, but `destroy!` on the result of a reader was never correct and is now documented as an error. Code doing it was corrupting the owner's handle already |

## 17. Migration

Shims live in `src/deprecations.jl`, generated from `scripts/renames.jl` ([§1.1](#1.1-What-these-rules-cover)).
They are removed in v0.6.

### 17.1 The shims have to warn, and `@deprecate` does not

`Base.@deprecate` calls `Base.depwarn`, which checks `JLOptions().depwarn` and prints nothing unless Julia was started with `--depwarn=yes`. The default is `no`. A migration window that prints nothing is not a migration window: a user would see silence for a whole minor cycle, then a `MethodError` at v0.6 with no warning history to explain it.

So the shims use a local `@renamed old new` macro whose body wraps `@warn … maxlog = 1`, which prints whatever the flag is set to:

```julia
@renamed getcorners_dmda    corners
@renamed getghostcorners    ghost_corners
@renamed getDM              dm
@renamed gettype            type_name
@renamed destroy            destroy!
```

CI also runs the suite with `--depwarn=yes`, so `@test_deprecated` stays meaningful for anything still on `Base.@deprecate`.

### 17.2 Forwarding shims and erroring stubs

A shim forwards only when the old call and the new one mean the same thing.

| Kind | Shim | Cases in v0.5 |
|---|---|---|
| Name changed, meaning unchanged | forwards, warns once | the bulk of the rename table |
| Arguments reordered, types differ | forwards, warns once | `dm_project_function!`, `dm_project_field!`, `dm_global_to_local!`, `dm_local_to_global!` |
| Argument dropped | keeps the argument, warns when it is passed | `ownership_range(A, false)`, `set_type!(ksp, "gmres")` |
| Argument types unchanged, **meaning changed** | throws, naming the replacement | **none** |

The last row is policy with no instances, and that is the point of writing it down. Where the argument types are identical and only the meaning of the values changes, dispatch cannot tell the calls apart and a forwarding shim corrupts data in silence. [§12.1](#12.1-Index-base) was drafted with `set_values!` in that row and then kept the 0-based convention instead, precisely to avoid needing it. Any future change that would land here should be reconsidered first; if it survives, it throws rather than forwards.

### 17.3 Test coverage

The main test suite is converted to the new names, so CI exercises the API that ships.
`test/test_deprecations.jl` calls every shim and asserts the warning is emitted, and asserts a throw for any entry marked as an erroring stub. It is generated from `scripts/renames.jl` too, so a rename cannot land without its shim being covered.

### Rename table

#### `dm.jl`

| v0.4 | v0.5 |
|---|---|
| `destroy` | `destroy!` |
| `getinfo` | `info` (drops the duplicate `s` field, `dof` becomes `ndofs` and `mpi_proc_size` becomes `procs`, see [§12](#12.-Return-values)) |
| `getcorners`, `getcorners_dmda` | `corners` |
| `getghostcorners`, `getghostcorners_dmda` | `ghost_corners` |
| `dm_local_to_global`, `dm_local_to_global!` | `local_to_global`, `local_to_global!` |
| `dm_global_to_local`, `dm_global_to_local!` | `global_to_local`, `global_to_local!` |
| `setuniformcoordinates_dmda!` | `set_uniform_coordinates!` |
| `coordinatesDMLocalVec` | `local_coordinates` |
| `getlocalcoordinatearray` | `local_coordinate_array` |
| `MatAIJ` | `PetscMat` constructor |
| `DMGlobalVec` | `global_vec` (merged, see below) |
| `DMLocalVec` | `local_vec` (merged, see below) |
| `getdimension` | `Base.ndims` |
| `setfromoptions!` | `set_from_options!` |

`DMGlobalVec` in `dm.jl` and `dm_create_global_vec` in `dmplex.jl` both call `DMCreateGlobalVector`; the local pair duplicates likewise.
Each pair collapses to one name, so two shims point at each replacement.

#### `dmstag.jl`

| v0.4 | v0.5 |
|---|---|
| `getcorners_dmstag` | `corners` (method on `DMStag`) |
| `getghostcorners_dmstag` | `ghost_corners` (method on `DMStag`; it has no `nextra` field, which `DMStagGetGhostCorners` never reported although the v0.4 docstring promised it) |
| `local_indices_dmstag` | `local_indices` |
| `global_indices_dmstag` | `global_indices` |
| `setuniformcoordinates_stag!` | `set_uniform_coordinates!` |
| `DMStagDOF_Slot` | `dof_slot` |
| `to_petscint_tuple` | unchanged (internal) |

#### `dmda.jl`

| v0.4 | v0.5 |
|---|---|
| `reshapelocalarray` | `reshape_local_array` |
| `localinteriorlinearindex` | `local_interior_linear_index` |
| `dmda_star_fd_coloring` | `star_fd_coloring(da::DMDA{L,2})` (drops `petsclib` per [§8](#8.-Argument-order), loses its export per [§13](#13.-Exports), and every index field gains a `_0b` or `_1b` suffix per [§12.1](#12.1-Index-base)) |
| `ndofs` | unchanged (closed list) |

#### `dmplex.jl`

| v0.4 | v0.5 |
|---|---|
| `isplexsimplex` | `issimplex` |
| `plexdistribute!` | `distribute!` |
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
| `vtk_save!`, `vtk_save_fields!` | `save_vtk!` (drops `petsclib`, see [§8](#8.-Argument-order)) |
| `vtk_merge_tensor!` | unchanged (free function, no PETSc object) |
| `setfield!` | `set_field!` (resolves the `Base.setfield!` shadow) |
| `dmclone` | `clone` |
| `plex_set_snes_local_fem!` | `set_snes_local_fem!` |
| `snes_set_jacobian_null_space!` | `set_jacobian_nullspace!` |
| `fe_compose_constant_null_space!` | `compose_constant_nullspace!` |
| `add_boundary!`, `add_natural_boundary!` | unchanged |
| `create_split_boundary_labels!` | unchanged |
| `set_constants!`, `set_exact_solution!` | unchanged |
| `set_residual!`, `set_jacobian!`, `set_jacobian_preconditioner!` | unchanged |
| `@petsc_residual_fn`, `@petsc_jacobian_fn` | `@residual_fn`, `@jacobian_fn` |
| `@petsc_bd_fn`, `@petsc_simple_fn` | `@bd_fn`, `@simple_fn` |

#### `ksp.jl`, `snes.jl`

| v0.4 | v0.5 |
|---|---|
| `PetscKSP`, `PetscSNES` (types) | `KSP`, `SNES` |
| `KSP(…)`, `SNES(…)` (factories) | type constructors |
| `getDM` | `dm` |
| `setDM!` | `set_dm!` |
| `get_solution` | `solution` |
| `gettype` | `type_name` (see [§3](#3.-Accessors-and-setters)) |
| `destroy` | `destroy!` |
| `setcomputeoperators!` | `set_compute_operators!` |
| `setcomputerhs!` | `set_compute_rhs!` |
| `setfunction!` | `set_function!` |
| `setjacobian!` | `set_snes_jacobian!` (see note) |
| `setconvergencetest!` | `set_convergence_test!` (callback first, [§8.1](#8.1-Callbacks-come-first)) |
| `Fn_SNESSetConvergenceTest` | `SNESSetConvergenceTestFn` ([§5.5](#5.5-Abstract,-callback-and-wrapper-types)) |
| — | `narrow`, new ([§5.4](#5.4-DMs-of-unknown-provenance)) |
| — | `set_type!(obj, ::Symbol)` for Vec, Mat, KSP, SNES and DM, new: v0.4 had it on TS only ([§3.1](#3.1-What-accessors-return-and-setters-take)) |

`setjacobian!` and `dmplex.jl`'s `set_jacobian!` share an English word and nothing else. They stay separate names under [§4.1](#4.1-When-the-subject-is-not-the-first-argument), because [§8.1](#8.1-Callbacks-come-first) puts the callback in argument 1 and so leaves `snes` unable to carry the name.

`setfunction!` and `setjacobian!` each ship two argument orders in v0.4 (`src/snes.jl`), subject-first and callback-first. Only callback-first survives, which is a break with no shim: see [§16](#16.-Breaking-changes-without-a-shim).

#### `ts.jl`

No v0.4 column: PETSc.jl had no high-level TS, so every name here is new and written to these rules from the start.
They are listed because the rename table is the register of what is public, not only of what changed.

| v0.5 | Kind |
|---|---|
| `TS` | type, exported per [§13](#13.-Exports) |
| `set_rhs_function!`, `set_rhs_jacobian!` | callback setters, callback first per [§8.1](#8.1-Callbacks-come-first) |
| `set_ifunction!`, `set_ijacobian!` | callback setters, likewise |
| `set_monitor!` | callback setter, likewise |
| `set_user_ctx!`, `user_ctx` | writer and reader |
| `set_solution!` | writer |
| `set_time!`, `set_timestep!`, `set_max_time!`, `set_max_steps!` | writers |
| `current_time`, `timestep`, `max_time`, `max_steps` | noun readers per [§3](#3.-Accessors-and-setters) |
| `set_tolerances!`, `tolerances` | writer and reader. `tolerances` returns `(; atol, rtol, vatol, vrtol)`, PETSc parameter names per [§12](#12.-Return-values), and its two vectors are **borrowed** ([§3.3](#3.3-What-an-accessor-hands-back)) |
| `set_type!`, `type_name` | shared with `ksp.jl` and `snes.jl`, `Symbol` per [§3.1](#3.1-What-accessors-return-and-setters-take). `type_name` replaces `ts.jl`'s v0.4 spelling `type`, the one name in this file the document did change |
| `set_adapt_type!` | writer, `Symbol`. Prefixed per [§4.1](#4.1-When-the-subject-is-not-the-first-argument): the target is the `TSAdapt` fetched inside, not the `TS` |
| `set_problem_type!`, `set_exact_final_time!` | writers taking PETSc **enums** (`TSProblemType`, `TSExactFinalTimeOption`), per [§3.1](#3.1-What-accessors-return-and-setters-take)'s enum rule, not `Symbol` |
| `solve!`, `step!`, `interpolate!`, `reset!` | mutation per [§7](#7.-Mutation) |
| `converged_reason`, `solve_time`, `step_number` | noun readers |
| `snes_iterations`, `ksp_iterations`, `snes_failures` | noun readers, prefixed per [§4.1](#4.1-When-the-subject-is-not-the-first-argument): the counts belong to the sub-solver |
| `step_rejections` | noun reader, count, belongs to the `TS` itself |
| `snes`, `ksp`, `dm`, `solution` | readers returning **borrowed** handles ([§3.3](#3.3-What-an-accessor-hands-back)) |
| `destroy!` | mutation, shared |

`current_time` rather than `time`, because `Base.time` exists and [§11](#11.-Extending-Base-and-stdlib) does not allow shadowing it with something unrelated.

`type_name(ts)` returns `Union{Nothing, Symbol}`. `TSSetFromOptions` runs in `solve!` rather than in the constructor, so that a DM and the callbacks attached in between are visible to it, and until then PETSc reports no type at all. A reader that can decline to answer is a hole in [§3](#3.-Accessors-and-setters); it is recorded here rather than papered over with a default, because the alternative is to invent a type the object does not have.

#### `vec.jl`, `mat.jl`

| v0.4 | v0.5 |
|---|---|
| `VecSeq` | `PetscVec` constructor |
| `MatCreateSeqAIJ`, `MatSeqAIJ`, `MatSeqDense`, `MatSeqAIJWithArrays`, `MatAIJ` | `PetscMat` constructor |
| `VecPtr`, `MatPtr` | unchanged: they are the wrapper **types** for a handle PETSc owns ([§5.5](#5.5-Abstract,-callback-and-wrapper-types)), not spellings of a constructor |
| `unsafe_localarray`, `wrap_localarray` | `unsafe_local_array`, `wrap_local_array` |
| `acquire_petsc_local_array` | `acquire_local_array` |
| `release_petsc_local_array` | `release_local_array` |
| `get_petsc_arrays` | `local_arrays` |
| `restore_petsc_arrays` | `restore_local_arrays!` |
| `withlocalarray!` | `with_local_array!` |
| `ghostupdate!` | `ghost_update!` |
| `ghostupdatebegin!`, `ghostupdateend!` | `ghost_update_begin!`, `ghost_update_end!` |
| `ownershiprange` | `ownership_range` (drops the positional `base_one`, see [§12.1](#12.1-Index-base)) |
| `setvalues!` | `set_values!` (stays 0-based; `row0idxs`/`col0idxs` become `rows_0b`/`cols_0b`, see [§12.1](#12.1-Index-base)) |
| `addindex!` | `add_index!` |
| `determine_memtype` | `memtype` |
| `as_petsc_vec` | `PetscVec` constructor |
| `array_type`, `memtype_backend` | unchanged |
| `make_local_array` | unchanged (internal) |
| `assemble!`, `setup!` | unchanged |
| — | `owns`, new: it answers whether the wrapper owns its handle, which is what `destroy!` consults ([§3.3](#3.3-What-an-accessor-hands-back)) |
| `destroy` | `destroy!` |

#### `init.jl`

| v0.4 | v0.5 |
|---|---|
| `initialized` | `isinitialized` |
| `finalized` | `isfinalized` |
| `check_petsc_wrappers_version` | `check_wrappers_version` (returns a `NamedTuple`, see [§12](#12.-Return-values)) |
| `scalartype`, `inttype` | unchanged (closed list) |
| `initialize`, `finalize` | unchanged, not exported ([§13](#13.-Exports)) |
| `check_initialized` | unchanged (internal, see [§14](#14.-Errors)) |
| `tao_usable_after_reinitialize` | unchanged (a platform predicate, not a PETSc accessor) |
| `isdestroyable` | unchanged (internal) |

#### Internals

Not part of the API, so no shims.
They drop the `_` prefix, which was standing in for "internal" and is not how Julia expresses that:

| v0.4 | v0.5 |
|---|---|
| `_petsc_link` | `petsc_link` |
| `_petsc_subst` | `petsc_subst` |
| `_vtk_merge_one_tensor!` | `vtk_merge_one_tensor!` |
| `_build_petsc_options` | `build_petsc_options` |
| `_ensure_library_handle` | `ensure_library_handle` |
| `_ensure_mpi_initialized` | `ensure_mpi_initialized` |
| `_library_ptr` | `library_ptr` |
| `_post_initialize` | `post_initialize` |
| `_release_library_handle` | `release_library_handle` |
| `_doc_external`, `_lib_handles`, `_petsc_program_name` | prefix dropped likewise. `_doc_external` also stays as a `const` alias, because `src/autowrapped/` interpolates that spelling into several thousand docstrings |
| `_errorcode`, `_run_callback`, `_with_options` (`ts.jl`) | prefix dropped likewise |

The exceptions keep their underscore, because there it separates an inner worker from the wrapper of the same name rather than marking visibility:

| v0.4 | v0.5 |
|---|---|
| `_mul!` | unchanged (worker for `mul!`) |
| `_unsafe_localarray` | `_unsafe_local_array` (worker for `unsafe_local_array`) |
| `get_petsc_arrays_impl` | `_local_arrays` (worker for `local_arrays`) |
| `restore_petsc_arrays_impl` | `_restore_local_arrays!` |

#### Types ([§5.5](#5.5-Abstract,-callback-and-wrapper-types))

| v0.4 | v0.5 |
|---|---|
| `AbstractPetscKSP`, `AbstractPetscSNES` | `AbstractKSP`, `AbstractSNES` |
| `AbstractPETScMemBackend` | `AbstractPetscMemBackend` |
| `Fn_KSPComputeRHS`, `Fn_KSPComputeOperators` | `KSPComputeRHSFn`, `KSPComputeOperatorsFn` |
| `Fn_SNESSetFunction`, `Fn_SNESSetJacobian` | `SNESSetFunctionFn`, `SNESSetJacobianFn` |
| `TSSetRHSFunctionFn`, `TSSetRHSJacobianFn`, `TSSetIFunctionFn`, `TSSetIJacobianFn`, `TSMonitorSetFn` | unchanged (new, already §5.5) |
| `MatAT` | deleted (no uses; it would have been `MatOrTranspose`) |
| `MatShell`, `MatOp`, `MatPtr`, `VecPtr` | unchanged |
| `AbstractPetscDS`, `PetscDS` | unchanged |
| `DMStagGetIndices` | removed (already deprecated in v0.4) |

#### `options.jl`, `sys.jl`

| v0.4 | v0.5 |
|---|---|
| `Options` (factory) | `PetscOptions` constructor |
| `getcomm` | `comm` |
| `typedget` | `parse_option` |
| `parse_options` | unchanged |

`typedget` operates on a plain `NamedTuple`, not on `PetscOptions`, and coerces the looked-up value to the type of the supplied default.
`parse_option` names that, and pairs with `parse_options(args)` in the same file.

#### `init.jl`, `audit.jl`

| v0.4 | v0.5 |
|---|---|
| `set_petsclib` | `PetscLibType` constructor |
| `library_info` | returns a `NamedTuple`, printed via `show` |
| `audit_petsc_file` | `audit_file` |
| `set_library!`, `unset_library!` | unchanged ([§7](#7.-Mutation) global-state clause) |
| `audit_walk`, `audit_targets`, `audit_report`, `audit_creator`, `audit_destroyer`, `audit_callee`, `audit_argnames`, `audit_isbroadcast`, `audit_hasparseerror` | unchanged (internal) |

`library_info` printed a report and returned `nothing`, so its name promised data it never handed back.
It now returns the values and gets a `show` method, which keeps the REPL output and makes the data reachable from tests:

```julia
info = library_info()
info.scalar          # Float64

julia> library_info()
Source  : LocalPreferences.toml
Path    : /usr/lib/libpetsc.so
```

### Two consequences worth knowing

`audit_file` (formerly `audit_petsc_file`) regex-matches the API's own names to pair object creations against `destroy` calls, so the rename blinds it: `destroy` becomes `destroy!` and `MatSeqAIJ` becomes `PetscMat`, the patterns stop matching, and it reports no leaks on leaking code. It is rewritten in the same PR to walk the parsed AST against creator and destroyer name sets generated from `scripts/renames.jl` ([§1.1](#1.1-What-these-rules-cover)), which is what stops the next rename blinding it again.

The typed DM hierarchy ([§5.3](#5.3-DM-flavour-is-a-type,-not-a-string)) widens autowrapped signatures from `PetscDM{PetscLib}` to `AbstractPetscDM{PetscLib}`, roughly 3500 occurrences in `src/autowrapped/DM_wrappers.jl` alone. Mechanical over generated files, but the generator in `wrapping/` has to emit the wider type too, or the next regeneration undoes it.
