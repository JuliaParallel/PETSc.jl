"""
    PETScTSAlgorithm

Abstract supertype for the PETSc TS time-integration algorithms exposed to
SciMLBase / OrdinaryDiffEq via the `PETScSciMLExt` extension.

Concrete subtypes carry just enough information to configure a PETSc `TS`
object: a top-level `TSSetType` argument, an optional subtype string, and an
algorithm-local `petsc_options::Vector{String}` of raw PETSc CLI tokens (e.g.
`["-snes_fd", "-ts_max_steps", "100"]`). The actual `solve` / `init` / `step!`
methods live in the extension and become available once `SciMLBase` is loaded.

Subtypes do **not** inherit from `SciMLBase.AbstractODEAlgorithm`, so they can
be defined here without a hard `SciMLBase` dependency. Dispatch into
SciMLBase's `solve` / `init` happens on the concrete type, which is sufficient.
"""
abstract type PETScTSAlgorithm end

"""
    TSRK(subtype::String[, petsc_options])

PETSc explicit Runge-Kutta integrator (`TSSetType(ts, "rk")`,
`TSRKSetType(ts, subtype)`). Suitable for non-stiff ODEs.

# Arguments
- `subtype`: PETSc RK subtype string, e.g. `"3bs"`, `"4"`, `"5dp"`, `"5bs"`.
  The full list comes from `TSRKType` (see PETSc docs / `-ts_rk_type`).
- `petsc_options`: optional `Vector{String}` of raw PETSc CLI tokens applied
  via `TSSetFromOptions` after the subtype is set.

# Example
```julia
using PETSc, OrdinaryDiffEq
prob = ODEProblem(f!, u0, tspan)
sol = solve(prob, PETSc.TSRK("3bs"); dt = 0.1)
```

# Common `petsc_options`
- `["-ts_adapt_type", "none"]` to disable PETSc's adaptive controller.
- `["-ts_max_steps", "1000"]` to cap the step count from PETSc's side.
- `["-ts_monitor"]` to print PETSc's per-step monitor (in addition to
  `save_everystep` on the SciML side).
"""
struct TSRK <: PETScTSAlgorithm
    subtype::String
    petsc_options::Vector{String}
end
TSRK(subtype::String) = TSRK(subtype, String[])
TSRK(subtype::String, petsc_options) =
    TSRK(subtype, String[String(s) for s in petsc_options])

"""
    TSRosW(subtype::String[, petsc_options])

PETSc Rosenbrock-W (linearly implicit) integrator
(`TSSetType(ts, "rosw")`, `TSRosWSetType(ts, subtype)`). Suitable for stiff
ODEs without a user-supplied analytic Jacobian — pass `"-snes_fd"` in
`petsc_options` to ask PETSc to compute Jacobians by finite differences.

# Arguments
- `subtype`: PETSc Rosenbrock-W subtype string, e.g. `"ra34pw2"`, `"rodas3"`,
  `"2m"`. See PETSc docs / `-ts_rosw_type` for the full list.
- `petsc_options`: optional `Vector{String}` of raw PETSc CLI tokens. For
  most stiff problems you will want at least `"-snes_fd"`.

# Example
```julia
sol = solve(prob, PETSc.TSRosW("ra34pw2", ["-snes_fd"]); dt = 1e-3)
```
"""
struct TSRosW <: PETScTSAlgorithm
    subtype::String
    petsc_options::Vector{String}
end
TSRosW(subtype::String) = TSRosW(subtype, String[])
TSRosW(subtype::String, petsc_options) =
    TSRosW(subtype, String[String(s) for s in petsc_options])

"""
    TSImplicit(subtype::String[, theta::Real][, petsc_options])

PETSc fully-implicit time integrator. The `subtype` selects the PETSc TS
type via `TSSetType(ts, subtype)`. Supported values:

| `subtype`  | PETSc method                | Notes                                                   |
|:-----------|:----------------------------|:--------------------------------------------------------|
| `"beuler"` | Backward Euler              | Order 1, L-stable                                       |
| `"cn"`     | Crank-Nicolson              | Order 2, A-stable; PETSc's endpoint-stage variant       |
| `"theta"`  | Theta method                | `theta` selects between BE (`1.0`) and CN (`0.5`)       |
| `"bdf"`    | Backward differentiation    | Set order via `petsc_options = ["-ts_bdf_order", "3"]`  |

`theta` defaults to `0.5` and is consulted only when `subtype == "theta"`.

Pass `"-snes_fd"` in `petsc_options` to let PETSc compute Jacobians via
finite differences when no user Jacobian is provided.

# Examples
```julia
solve(prob, PETSc.TSImplicit("beuler", ["-snes_fd"]); dt = 0.01)
solve(prob, PETSc.TSImplicit("theta", 0.5, ["-snes_fd"]); dt = 0.05)
solve(prob, PETSc.TSImplicit("bdf", ["-snes_fd", "-ts_bdf_order", "3"]); dt = 1e-3)
```
"""
struct TSImplicit <: PETScTSAlgorithm
    subtype::String
    theta::Float64
    petsc_options::Vector{String}
end
TSImplicit(subtype::String) = TSImplicit(subtype, 0.5, String[])
TSImplicit(subtype::String, theta::Real) = TSImplicit(subtype, Float64(theta), String[])
TSImplicit(subtype::String, petsc_options) =
    TSImplicit(subtype, 0.5, String[String(s) for s in petsc_options])
TSImplicit(subtype::String, theta::Real, petsc_options) =
    TSImplicit(subtype, Float64(theta), String[String(s) for s in petsc_options])

"""
    TSARKIMEX(subtype::String[, petsc_options])

PETSc Additive Runge-Kutta IMEX integrator (`TSSetType(ts, "arkimex")`,
`TSARKIMEXSetType(ts, subtype)`). Designed to integrate `SplitODEProblem`s
of the form `u' = f1(u,p,t) + f2(u,p,t)` where `f1` is the stiff/implicit
part and `f2` is the non-stiff/explicit part.

When a non-`SplitODEProblem` is passed, the full RHS is treated as the
implicit part and the explicit RHS is left at zero.

# Arguments
- `subtype`: PETSc ARKIMEX subtype string, e.g. `"2e"`, `"3"`, `"4"`,
  `"5"`. See PETSc docs / `-ts_arkimex_type` for the full list.
- `petsc_options`: optional `Vector{String}`. Typically include `"-snes_fd"`.

# Example
```julia
prob = SplitODEProblem(f1!, f2!, u0, tspan)   # f1 stiff/implicit, f2 explicit
sol = solve(prob, PETSc.TSARKIMEX("2e", ["-snes_fd"]); dt = 0.05)
```
"""
struct TSARKIMEX <: PETScTSAlgorithm
    subtype::String
    petsc_options::Vector{String}
end
TSARKIMEX(subtype::String) = TSARKIMEX(subtype, String[])
TSARKIMEX(subtype::String, petsc_options) =
    TSARKIMEX(subtype, String[String(s) for s in petsc_options])

"""
    TSGeneric(ts_type::String[, petsc_options]; explicit::Bool = false)

Pass-through algorithm that calls `TSSetType(ts, ts_type)` directly without
any subtype-specific configuration. Useful for PETSc TS types that do not
have a dedicated wrapper here yet (e.g. `"alpha"`, `"glle"`, `"glee"` for
implicit families; `"euler"`, `"ssp"` for explicit families).

By default `TSGeneric` registers an IFunction (residual `udot - f(u,p,t)`)
and sets `TS_NONLINEAR`, which is what implicit / Rosenbrock-style PETSc
TS types expect. Pass `explicit = true` to instead register the RHS via
`TSSetRHSFunction`, which is what explicit-only PETSc TS types like
`"euler"` and `"ssp"` require.

If you supply an `explicit = false` `TSGeneric` for a TS type that PETSc
classifies as explicit-only, `TSSetUp` / `TSStep` will fail with a clear
PETSc error indicating the residual is not consumed. In that case, retry
with `explicit = true`.

For the standard families prefer the dedicated wrappers
`TSRK` / `TSRosW` / `TSImplicit` / `TSARKIMEX`.

# Examples
```julia
# Implicit / Rosenbrock-style: default `explicit = false`.
solve(prob, PETSc.TSGeneric("alpha", ["-snes_fd"]); dt = 0.01)

# Explicit-only PETSc TS types: opt in via `explicit = true`.
solve(prob, PETSc.TSGeneric("euler"; explicit = true); dt = 0.01)
solve(prob, PETSc.TSGeneric("ssp"; explicit = true); dt = 0.01)
```
"""
struct TSGeneric <: PETScTSAlgorithm
    ts_type::String
    explicit::Bool
    petsc_options::Vector{String}
end
TSGeneric(ts_type::String; explicit::Bool = false) =
    TSGeneric(ts_type, explicit, String[])
TSGeneric(ts_type::String, petsc_options; explicit::Bool = false) =
    TSGeneric(ts_type, explicit, String[String(s) for s in petsc_options])
