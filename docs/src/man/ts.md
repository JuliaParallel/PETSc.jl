# TS

The TS (Time Stepping) module integrates ordinary differential equations and differential algebraic equations in time, and is the layer PETSc builds its implicit, explicit and IMEX methods on.

The interface described here is the high-level one. It owns the `@cfunction` trampolines, keeps your callbacks rooted for as long as the stepper lives, and frees the object for you on a sequential communicator. The [low-level TS bindings](ts_lowlevel.md) remain available and unchanged.

## Overview

A problem is given to TS in one of two forms, or in both at once:

- **Explicit**, ``du/dt = G(t, u)``, through [`PETSc.set_rhs_function!`](@ref) and optionally [`PETSc.set_rhs_jacobian!`](@ref).
- **Implicit**, ``F(t, u, du/dt) = 0``, through [`PETSc.set_ifunction!`](@ref) and [`PETSc.set_ijacobian!`](@ref). This form also covers a DAE and a problem with a mass matrix.

Setting both is how an IMEX method splits a problem into the stiff part it treats implicitly and the rest it treats explicitly.

## Creating a stepper

```julia
ts = PETSc.TS(petsclib, MPI.COMM_WORLD)

# or with options, which are applied when you call solve!
ts = PETSc.TS(petsclib, MPI.COMM_WORLD;
    ts_type = "bdf",
    ts_adapt_type = "basic",
)
```

Options given here are held and applied inside [`PETSc.solve!`](@ref) rather than at construction, so a DM and the callbacks attached afterwards are in place before PETSc reads them. Until then the object has no type, and [`PETSc.type`](@ref) answers `nothing`.

On a communicator of size 1 the garbage collector calls [`PETSc.destroy!`](@ref). On a larger one, destruction is yours to do, since collection is asynchronous and the call is collective.

## An explicit problem

Solving ``du/dt = -u`` from ``u(0) = 1``:

```julia
ts = PETSc.TS(petsclib, MPI.COMM_SELF)
PETSc.set_type!(ts, :rk)

u = PETSc.VecSeq(petsclib, 1)
u[1] = 1.0
PETSc.assemble!(u)

PETSc.set_rhs_function!(ts) do F, ts, t, u
    PETSc.withlocalarray!((u, F); read = (true, false), write = (false, true)) do ua, Fa
        Fa[1] = -ua[1]
    end
    return 0
end

PETSc.set_time!(ts, 0.0)
PETSc.set_timestep!(ts, 0.01)
PETSc.set_max_time!(ts, 1.0)

PETSc.solve!(u, ts)
```

The callback comes first in the argument list, so `do` block syntax works. Pass it in the other order, `set_rhs_function!(ts, f!)`, when it is already a value.

## An implicit problem

The same equation written as ``F(t, u, u_t) = u_t + u = 0``, stepped with backward Euler. The Jacobian asked for is $\partial F/\partial u + \sigma\, \partial F/\partial u_t$, where PETSc supplies the shift ``\sigma``:

```julia
ts = PETSc.TS(petsclib, MPI.COMM_SELF)
PETSc.set_type!(ts, :beuler)

J = PETSc.MatSeqAIJ(petsclib, 1, 1, petsclib.PetscInt(1))

PETSc.set_ifunction!(ts) do F, ts, t, u, u_t
    PETSc.withlocalarray!(
        (u, u_t, F);
        read = (true, true, false),
        write = (false, false, true),
    ) do ua, uta, Fa
        Fa[1] = uta[1] + ua[1]
    end
    return 0
end

PETSc.set_ijacobian!(ts, J) do A, P, ts, t, u, u_t, shift
    A[1, 1] = shift + 1
    PETSc.assemble!(A)
    return 0
end
```

A Jacobian callback is always handed both the Jacobian `A` and the preconditioning matrix `P`. They are the same object unless you passed two, which is worth testing with `A.ptr == P.ptr` before filling `P` a second time.

## Callback signatures

| Setter | Callback |
| --- | --- |
| [`PETSc.set_rhs_function!`](@ref) | `f!(F, ts, t, u)` |
| [`PETSc.set_rhs_jacobian!`](@ref) | `updateJ!(A, P, ts, t, u)` |
| [`PETSc.set_ifunction!`](@ref) | `f!(F, ts, t, u, u_t)` |
| [`PETSc.set_ijacobian!`](@ref) | `updateJ!(A, P, ts, t, u, u_t, shift)` |
| [`PETSc.set_monitor!`](@ref) | `f(ts, step, t, u)` |

Each may return a PETSc error code; any other return value counts as success.

A callback that raises a Julia exception is reported and turned into a PETSc failure rather than being allowed to escape into C, where it would take the process down with it. The error is logged with its backtrace and [`PETSc.solve!`](@ref) then raises a `PetscError`.

## Passing your own data

Anything stored with [`PETSc.set_user_ctx!`](@ref) is handed back as a trailing argument, to whichever callbacks have a method that accepts one:

```julia
PETSc.set_user_ctx!(ts, (; viscosity = 1e-3))

PETSc.set_rhs_function!(ts) do F, ts, t, u, ctx
    # ctx.viscosity is available here
    return 0
end
```

The object is held on the Julia side by `ts`, so it stays alive without any pinning of your own.

## Watching a solve

```julia
PETSc.set_monitor!(ts) do ts, step, t, u
    @printf("step %3d  t = %.4f\n", step, t)
    return 0
end
```

The monitor runs once before the first step and once after each accepted one.
It does not replace the monitors PETSc installs from the options database, such as `-ts_monitor`.

## Solving, and reading the result

```julia
PETSc.solve!(u, ts)          # integrate, starting from u
PETSc.solve!(ts)             # or from the vector already set on ts

PETSc.converged_reason(ts)   # why it stopped
PETSc.solve_time(ts)         # the time actually reached
PETSc.step_number(ts)        # steps taken
PETSc.snes_iterations(ts)    # nonlinear iterations, summed over the steps
```

The final step lands exactly on [`PETSc.set_max_time!`](@ref) by default. 
This differs from PETSc, whose own default, `TS_EXACTFINALTIME_UNSPECIFIED`, integrates past the requested time without saying so; pass `exact_final_time` to the constructor or call [`PETSc.set_exact_final_time!`](@ref) to choose otherwise.

For finer control, [`PETSc.step!`](@ref) takes a single step and [`PETSc.interpolate!`](@ref) samples the solution inside the step just taken.

## Cleaning up

```julia
PETSc.destroy!(ts)
```

Safe to call more than once, and a no-op on a handle left over from a previous initialize/finalize cycle. `PETSc.destroy` also works, for consistency with the rest of the package.

## Functions

```@autodocs
Modules = [PETSc]
Pages   = ["ts.jl"]
```
