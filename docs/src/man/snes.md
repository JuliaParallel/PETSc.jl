# SNES

The SNES (Scalable Nonlinear Equations Solvers) module provides methods for solving nonlinear systems of the form `F(x) = 0`. It builds on KSP for the linear solves within Newton-like methods.

## Overview

SNES provides:
- **Newton methods**: Newton line search, Newton trust region
- **Quasi-Newton**: L-BFGS, Broyden
- **Nonlinear Richardson**: With various line search strategies
- **FAS multigrid**: Full Approximation Scheme for nonlinear problems
- **Composite solvers**: Combine multiple nonlinear solvers

## Creating a SNES Solver

```julia
# Basic creation
snes = SNES(petsclib, MPI.COMM_WORLD)

# With options
snes = SNES(petsclib, MPI.COMM_WORLD;
    snes_type = "newtonls",
    snes_rtol = 1e-8,
    snes_max_it = 50
)
```

## Setting the Nonlinear Function

Define the residual function `F(x)`. The callback comes first in the argument
list, so `do` block syntax works; v0.4's subject-first order
(`setfunction!(snes, f!, v)`) is gone and has no shim
([naming conventions](naming.md), §8.1):

```julia
function residual!(fx, snes, x)
    # Compute F(x) and store in fx
    fx[1] = x[1]^2 + x[2] - 1
    fx[2] = x[1] + x[2]^2 - 1
    return 0
end

PETSc.set_function!(residual!, snes, f_vec)
```

## Setting the Jacobian

Define the Jacobian `J = dF/dx`. `set_snes_jacobian!` rather than
`set_jacobian!`, because `set_jacobian!` belongs to `PetscDS` and the callback
occupies argument 1 here ([naming conventions](naming.md), §4.1):

```julia
function jacobian!(J, snes, x)
    # Fill Jacobian matrix
    J[1, 1] = 2*x[1]
    J[1, 2] = 1.0
    J[2, 1] = 1.0
    J[2, 2] = 2*x[2]
    PETSc.assemble!(J)
    return 0
end

PETSc.set_snes_jacobian!(jacobian!, snes, J, J)  # (J, P), P the preconditioner matrix
```

## Using a DM

For PDE problems, associate the SNES with a DM:

```julia
PETSc.set_dm!(snes, dm)

# Get the DM from SNES. `dm` is both the accessor and the usual variable name,
# so bind the result to something else (see the naming conventions, §3.2).
d = PETSc.dm(snes)
```

## Solving

```julia
# Solve with initial guess x
PETSc.solve!(x, snes)

# Get the solution vector. It is a **borrowed** handle: it belongs to `snes`,
# and `destroy!` on it is a no-op (naming conventions, §3.3)
sol = PETSc.solution(snes)

# What PETSc did
PETSc.converged_reason(snes)     # positive when it converged
PETSc.iteration_number(snes)     # Newton iterations
PETSc.ksp_iterations(snes)       # linear iterations, summed over them
PETSc.function_norm(snes)        # residual norm at the last iterate
PETSc.ksp(snes)                  # the linear solver, borrowed

PETSc.destroy!(snes)
```

An iterate the residual cannot be evaluated at (a negative density, say) is reported from inside the residual with `PETSc.set_function_domain_error!(snes)`. A few solvers then cut the step; the others stop and `converged_reason` says `SNES_DIVERGED_FUNCTION_DOMAIN`. Throwing instead stops the solve, and `solve!` rethrows the exception.

## Common Solver Options

### Nonlinear Solver Types (`snes_type`)
- `newtonls` - Newton with line search (default)
- `newtontr` - Newton with trust region
- `nrichardson` - Nonlinear Richardson
- `qn` - Quasi-Newton (L-BFGS)
- `fas` - Full Approximation Scheme multigrid

### Convergence Options
- `snes_rtol` - Relative tolerance
- `snes_atol` - Absolute tolerance
- `snes_stol` - Step tolerance
- `snes_max_it` - Maximum iterations
- `snes_monitor` - Print residual each iteration

### Line Search Options
- `snes_linesearch_type` - `bt` (backtracking), `basic`, `l2`, `cp`

## Example: Full Setup

```julia
petsclib = PETSc.petsclibs[1]
PETSc.initialize(petsclib)

snes = SNES(petsclib, MPI.COMM_WORLD;
    snes_monitor = true,
    ksp_type = "gmres",
    pc_type = "ilu"
)

PETSc.set_function!(residual!, snes, f)
PETSc.set_snes_jacobian!(jacobian!, snes, J, J)
PETSc.set_from_options!(snes)

PETSc.solve!(x, snes)
PETSc.destroy!(snes)
```

## Functions

```@autodocs
Modules = [PETSc]
Pages   = ["snes.jl"]
```
