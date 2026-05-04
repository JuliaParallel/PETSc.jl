# SciML Integration (ODE Time Stepping)

`PETSc.jl` ships a package extension (`PETScSciMLExt`) that lets you solve
in-place `ODEProblem`s with PETSc's TS time integrators through the standard
[SciML](https://sciml.ai/) interface.

The extension activates automatically when `SciMLBase` is loaded — no
`DiffEqBase` is required:

```julia
using PETSc, SciMLBase
```

`OrdinaryDiffEq` also works as the trigger since it re-exports `SciMLBase`,
but the extension itself only depends on `SciMLBase`.

## Quick start

```julia
using PETSc, SciMLBase

f!(du, u, p, t) = (du[1] = -u[1]; nothing)
prob = ODEProblem(f!, [1.0], (0.0, 1.0))

sol = solve(prob, PETSc.TSRK("3bs"); dt = 0.1)                               # explicit RK
sol = solve(prob, PETSc.TSRosW("ra34pw2", ["-snes_fd"]); dt = 1e-3)          # Rosenbrock-W
sol = solve(prob, PETSc.TSImplicit("bdf", ["-snes_fd"]); dt = 1e-3)          # BDF / theta / CN / BEuler

f1!(du, u, p, t) = (du[1] = -u[1]; nothing)   # stiff / implicit part
f2!(du, u, p, t) = (du[1] = cos(t); nothing)  # non-stiff / explicit part
prob_split = SplitODEProblem(f1!, f2!, [1.0], (0.0, 1.0))
sol = solve(prob_split, PETSc.TSARKIMEX("2e", ["-snes_fd"]); dt = 0.05)      # IMEX

integrator = init(prob, PETSc.TSRK("5dp"); dt = 0.05)                        # step!/solve! interface
step!(integrator)
sol = solve!(integrator)
```

## Algorithm types

All algorithm types are subtypes of `PETSc.PETScTSAlgorithm` and are
accessible from the top-level `PETSc` module.

### `PETSc.TSRK` — Explicit Runge-Kutta

```julia
PETSc.TSRK(subtype::String[, petsc_options])
```

Calls `TSSetType(ts, "rk")` and `TSRKSetType(ts, subtype)`. Suitable for
non-stiff ODEs.

Common subtypes: `"3bs"` (Bogacki-Shampine order 3), `"4"` (classical
4th-order RK), `"5dp"` (Dormand-Prince order 5), `"5bs"`. The full list
is available via `-ts_rk_type` in PETSc's option database.

### `PETSc.TSRosW` — Rosenbrock-W (linearly implicit)

```julia
PETSc.TSRosW(subtype::String[, petsc_options])
```

Calls `TSSetType(ts, "rosw")` and `TSRosWSetType(ts, subtype)`. Suitable
for stiff ODEs without a user-supplied analytic Jacobian — pass `"-snes_fd"`
in `petsc_options` to ask PETSc to compute Jacobians by finite differences.

Common subtypes: `"ra34pw2"`, `"rodas3"`, `"2m"`. See `-ts_rosw_type`.

### `PETSc.TSImplicit` — Fully implicit

```julia
PETSc.TSImplicit(subtype::String[, theta::Real][, petsc_options])
```

Selects a fully-implicit PETSc TS type via `TSSetType(ts, subtype)`.

| `subtype`  | PETSc method             | Notes                                                  |
|:-----------|:-------------------------|:-------------------------------------------------------|
| `"beuler"` | Backward Euler           | Order 1, L-stable                                      |
| `"cn"`     | Crank-Nicolson           | Order 2, A-stable                                      |
| `"theta"`  | Theta method             | `theta` selects between BE (`1.0`) and CN (`0.5`)      |
| `"bdf"`    | Backward differentiation | Set order via `["-ts_bdf_order", "3"]` in options      |

`theta` defaults to `0.5` and is consulted only when `subtype == "theta"`.
Pass `"-snes_fd"` in `petsc_options` when no analytic Jacobian is available.

### `PETSc.TSARKIMEX` — Additive Runge-Kutta IMEX

```julia
PETSc.TSARKIMEX(subtype::String[, petsc_options])
```

Calls `TSSetType(ts, "arkimex")` and `TSARKIMEXSetType(ts, subtype)`.
Designed for `SplitODEProblem`s of the form `u' = f1(u,p,t) + f2(u,p,t)`
where `f1` is the stiff/implicit part and `f2` is the non-stiff/explicit part.

When a plain `ODEProblem` is passed, the full RHS is treated as the implicit
part and the explicit part is left at zero.

Common subtypes: `"2e"`, `"3"`, `"4"`, `"5"`. See `-ts_arkimex_type`.

### `PETSc.TSGeneric` — Pass-through

```julia
PETSc.TSGeneric(ts_type::String[, petsc_options]; explicit::Bool = false)
```

Calls `TSSetType(ts, ts_type)` directly without any subtype-specific
configuration. Useful for PETSc TS types that do not have a dedicated wrapper
(e.g. `"alpha"`, `"glle"`, `"glee"`).

Pass `explicit = true` for PETSc TS types that register an RHS function
rather than an IFunction (e.g. `"euler"`, `"ssp"`).

## PETSc options

Per-solver PETSc command-line options are passed on the algorithm object as a
`Vector{String}` of raw tokens:

```julia
alg = PETSc.TSRK("3bs", ["-ts_monitor", "-ts_max_steps", "100"])
alg = PETSc.TSImplicit("bdf", ["-snes_fd", "-ts_bdf_order", "3"])
```

The tokens are applied via `TSSetFromOptions` after the TS type and subtype
are configured, so they can override any default.

## Supported `solve` keywords

| Keyword            | Meaning                                                               |
|:-------------------|:----------------------------------------------------------------------|
| `dt`               | Initial time step                                                     |
| `adaptive`         | `true` (default) enables PETSc's `TSAdapt` controller; `false` fixes dt |
| `dtmin`, `dtmax`   | Step-size bounds for the adaptive controller                          |
| `reltol`, `abstol` | Scalar tolerances forwarded to `TSSetTolerances`                      |
| `maxiters`         | Maximum number of accepted steps                                      |
| `saveat`           | Times at which to save the solution (vector or scalar spacing)        |
| `save_everystep`   | Save at every step endpoint (default: `true` when `saveat` is empty)  |
| `save_start`       | Save the initial state `u(t0)` (default: `true`)                      |
| `save_end`         | Save the final state `u(tf)` (default: `true`)                        |
| `save_on`          | Master switch: suppress all trajectory output when `false`            |
| `callback`         | `DiscreteCallback` or `CallbackSet` of discrete callbacks             |
| `initialize_save`  | Run post-`initialize!` save record (default: `true`)                  |
| `petsclib`         | Override the PETSc library instance to use                            |

Unrecognized keywords are rejected with a clear `ArgumentError` rather than
silently dropped.

## Callbacks

`DiscreteCallback`s and `terminate!` are supported. Their `initialize` and
`finalize` hooks are called at the start and end of the solve, matching the
standard SciML lifecycle.

`ContinuousCallback`s are **not** supported and are rejected with
`ArgumentError`. Wrap event detection through PETSc's `TSSetEventHandler`
directly if needed.

`tstops` is also **not** yet supported (rejected with `ArgumentError`) because
the extension drives `TSStep` directly and cannot currently guarantee exact
landing on requested times.

## Current limitations

- Only real-valued, in-place `ODEProblem`s (`f!(du, u, p, t)`) are supported.
- Only forward integration (`tspan[1] < tspan[2]`) is supported.
- The extension requires a PETSc library built with `PetscReal = Float64`.
- TSIRK/Gauss methods require an AIJ sparse Jacobian registered via
  `TSSetIJacobian`; the extension does not yet set one up, so those methods
  are not available through the SciML interface.
- `ContinuousCallback`s and `tstops` are not yet honoured (see above).

## Algorithm docstrings

```@docs
PETSc.TSRK
PETSc.TSRosW
PETSc.TSImplicit
PETSc.TSARKIMEX
PETSc.TSGeneric
```
