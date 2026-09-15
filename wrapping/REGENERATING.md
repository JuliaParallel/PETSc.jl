# Regenerating the wrappers

Read this before running `generatejuliabindings.jl` over a new PETSc release.

## The short version

The generator has not received the fixes that were applied to its output. Four
PRs fixed correctness bugs in `src/autowrapped/`, and none of them touched
`wrapping/`. Regenerating today reintroduces all of them, including 76
`*Destroy` wrappers that segfault.

| PR | What it fixed | In the generator? |
|---|---|---|
| #254 | 2 TS wrappers taking output parameters as inputs | no |
| #257 | ~25 wrappers returning PETSc-owned arrays sized from an undefined variable | no |
| #258 | 6 wrappers that nulled the handle PETSc returned | no |
| #259 | 76 `*Destroy` wrappers passing the handle by value, 82 bad `Ref` tests, string types, optional NULL arrays | no |
| #263 | input arguments take the abstract wrapper type | **yes** |

So the useful work here is not rerunning the script. It is moving these seven
rules into the script first, then rerunning it.

## 1. Output arguments are decided by star count, not by PETSc's docs

The generator infers that an argument is an output from its pointer depth and
type name. PETSc's own man pages already say, under `Input Parameters:` and
`Output Parameters:`, which is which, and `getAPI.py` exposes those sections.
Use them.

Getting this wrong has produced three distinct failures:

- an output parameter treated as an input, so the caller had to pass a value
  the wrapper then overwrote (`TSGetConvergedReason`, `TSGetTolerances`)
- a returned handle set to `C_NULL` instead of what PETSc wrote
  (`SNESGetSolution`, `SNESGetSolutionUpdate`, `DMAdaptorAdapt`,
  `TSMonitorEnvelopeGetBounds`, `TaoLineSearchGetStartingVector`,
  `VecNestGetSubVec`)
- a wrapper taking the object it should return, and nulling the caller's handle
  on the way (`TSGetSNES`, `TSGetKSP`), which needed hand-written bindings in
  `src/ts.jl`

A wrapper that returns a handle allocates its own `Ref`, calls, then wraps:

```julia
dm_ = Ref{CDM}()
@chk ccall((:DMCreate, $petsc_library), PetscErrorCode, (MPI_Comm, Ptr{CDM}), comm, dm_)
dm = PetscDM(dm_[], petsclib)
return dm
```

## 2. `*Destroy` takes a pointer to the handle, not the handle

PETSc's destroy functions take `T*` and null the caller's copy. The generator
emitted the handle by value, which compiles because `Ptr{T}` converts silently
to `Ptr{Ptr{T}}`, and then segfaults. 76 functions were affected.

Correct output accepts either form and always passes a `Ref`:

```julia
@for_petsc function ISLocalToGlobalMappingDestroy(
    petsclib::$UnionPetscLib,
    mapping::Union{ISLocalToGlobalMapping, Ref{ISLocalToGlobalMapping}},
)
    mapping_ = mapping isa Base.RefValue ? mapping : Ref{ISLocalToGlobalMapping}(mapping)
    @chk ccall(
        (:ISLocalToGlobalMappingDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{ISLocalToGlobalMapping},),
        mapping_,
    )
    return nothing
end
```

## 3. Test `isa Base.RefValue`, never `isa Ref`

`Ptr{T} <: Ref{T}` in Julia, so the generated `x isa Ref ? x : Ref{T}(x)` passes
a bare handle straight through for exactly the handle types it was meant to
wrap. All 82 occurrences now test `isa Base.RefValue`.

Check after regenerating:

```
grep -rn "isa Ref ?" src/autowrapped/     # must return nothing
```

## 4. A returned PETSc-owned array must be sized from something in scope

The generator emits

```julia
range = unsafe_wrap(Array, range_[], VecGetLocalSize(petsclib, x); own = false)
```

for any function returning a pointer to an array, with a literal `x` that is
usually not an argument of that function. Every such call raises
`UndefVarError`. #257 fixed the DMPlex, IS, PetscSection and PetscSF subset.

**126 functions still carry this bug**, mostly in `DM_wrappers.jl` (33) and
`Mat_wrappers.jl` (29). `PetscLayoutGetRanges` is a plain example: its arguments
are `(petsclib, map)` and the body asks for `VecGetLocalSize(petsclib, x)`.

The length has to come from the paired count output or the matching size query,
which PETSc documents alongside the array: `DMPlexGetCone` is sized by
`DMPlexGetConeSize`, `DMPlexGetSupport` by `DMPlexGetSupportSize`, and so on.
Where the generator cannot work it out, it is better to emit nothing and let the
function be written by hand than to emit a call that cannot run.

`own = false` is right: the memory belongs to PETSc, and the paired `Restore`
takes the count and the view back.

Find the rest:

```
grep -rn "VecGetLocalSize(petsclib, x)" src/autowrapped/
```

## 5. Type-name arguments need a `String` overload

PETSc registers implementation names at runtime, so `KSPType` and the 59 aliases
like it are C string types rather than enums. Dispatch happens before `ccall`
converts anything, so a Julia `String` does not match the generated signature and
`MatSetType(petsclib, mat, "mpiaij")` raises a `MethodError`.

These aliases are not spelled consistently, which matters if you generate the
overload from the type:

| Spelling | Count | Examples |
|---|---|---|
| `Cstring` | 2 | `MatType`, `VecType` |
| `Ptr{Cchar}` | 58 | `KSPType`, `DMType`, `SNESType`, `AOType`, `PCType`, … |

The convenience overloads taking `AbstractString` live in
`src/string_wrappers.jl` and `src/string_wrappers_extra.jl`, and are not
generated. Keep them, or teach the generator to emit one for every argument
whose type is one of the 60.

## 6. Optional array arguments accept `C_NULL`

PETSc lets several preallocation arguments be `NULL`.
`MatSeqAIJSetPreallocation`, `MatMPIAIJSetPreallocation` and
`MatXAIJSetPreallocation` must accept `C_NULL` for the `nnz` arrays, the way
`MatCreateAIJ` already did.

## 7. Input arguments take the abstract type

This one is already in the generator, via `abstract_arg_type`, and
`julia_function_doc_header` applies it to `str_in` only. Do not let a refactor
apply it to `str_out_doc`.

An input is annotated with the abstract supertype, so a borrowed handle
(`VecPtr`), a shell matrix, or a DM carrying its flavour in the type reaches the
same method:

```julia
function DMCopyDMKSP(petsclib::PetscLibType, dmsrc::AbstractPetscDM, dmdest::AbstractPetscDM) end
```

A return position keeps the concrete name, because the wrapper constructs that
object:

```
dm::PetscDM = DMCreate(petsclib::PetscLibType, comm::MPI_Comm)
```

Struct fields in `struct_wrappers.jl` (`JacActionCtx`, `DMDALocalInfo`,
`TSMonitorDMDARayCtx`) also stay concrete, since an abstract field type boxes
the value.

`test/wrapper_signatures.jl` enforces this and names the offending signatures.

## Files under `src/autowrapped/` that are not generator output

Do not overwrite these blindly:

- `petsc_library.jl` is `prologue.jl` after `move_prologue`, but see below
- `DMaddons_wrappers.jl` contains hand-written multi-line wrappers
  (`DMProjectFunction`, `DMComputeL2Diff`)
- `Vecs_wrappers.jl` documents `VecSetValues` with a bare signature expression
  rather than a `function … end` stub
- `PC_wrappers.jl` is deliberately excluded from the include list in
  `petsc_library.jl`, because `PC` in ccall signatures still needs fixing

## `prologue.jl` has drifted from `petsc_library.jl`

31 lines differ. The generated file carries five hand-edits that were never put
back into the prologue, so `move_prologue("prologue.jl")` would destroy them:

- `computerhs!`, `computeops!` and `opts` fields on `PetscKSP`
- `user_ctx` and `opts` fields on `PetscSNES`
- `Base.unsafe_convert(::Type{Ptr{CIS}}, v::AbstractIS)`, needed by
  `DMGetStratumIS` and friends
- `const DMLabel = Ptr{Cvoid}` instead of `mutable struct DMLabel end`
- the commented-out `include("PC_wrappers.jl")`

Reconcile the two before regenerating. `wrapping/local_types.jl` is a near
duplicate of `prologue.jl` that nothing reads, and can go.

## After regenerating

```
grep -rn "isa Ref ?" src/autowrapped/                    # expect nothing
grep -rn "VecGetLocalSize(petsclib, x)" src/autowrapped/ # expect fewer than 189
julia --project=. test/wrapper_signatures.jl             # abstract input types
julia --project=. -e 'using Pkg; Pkg.test()'             # full suite
```

Also worth running, since it catches a signature change that makes two methods
equally specific:

```julia
using Test, PETSc
length(detect_ambiguities(PETSc; recursive = true))   # 161 as of #263
```

The suite is the real check. It caught every one of the bugs above except the
126 still open in section 4, which no test reaches because nothing calls those
functions yet.
