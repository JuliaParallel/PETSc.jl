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

Define the residual function `F(x)`:

```julia
function residual!(fx, snes, x)
    # Compute F(x) and store in fx
    fx[1] = x[1]^2 + x[2] - 1
    fx[2] = x[1] + x[2]^2 - 1
    return 0
end

setfunction!(snes, residual!, f_vec)
```

## Setting the Jacobian

Define the Jacobian `J = dF/dx`:

```julia
function jacobian!(J, snes, x)
    # Fill Jacobian matrix
    J[1, 1] = 2*x[1]
    J[1, 2] = 1.0
    J[2, 1] = 1.0
    J[2, 2] = 2*x[2]
    assemble!(J)
    return 0
end

setjacobian!(snes, jacobian!, J, J)  # (J, P) where P is preconditioner matrix
```

## Using a DM

For PDE problems, associate the SNES with a DM:

```julia
setDM!(snes, dm)

# Get the DM from SNES
dm = getDM(snes)
```

## Solving

```julia
# Solve with initial guess x
solve!(x, snes)

# Get solution vector
sol = get_solution(snes)
```

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

setfunction!(snes, residual!, f)
setjacobian!(snes, jacobian!, J, J)
setfromoptions!(snes)

solve!(x, snes)
```

## Functions

```@autodocs
Modules = [PETSc]
Pages   = ["snes.jl"]
```
