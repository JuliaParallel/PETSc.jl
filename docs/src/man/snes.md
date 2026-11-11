# SNES

The SNES (Scalable Nonlinear Equations Solvers) module provides methods for solving nonlinear systems of the form `F(x) = 0`. It builds on KSP for the linear solves within Newton-like methods.

## Overview

SNES provides:
- **Newton methods**: Newton line search, Newton trust region
- **Quasi-Newton**: L-BFGS, Broyden
- **Nonlinear Richardson**: With various line search strategies
- **FAS multigrid**: Full Approximation Scheme for nonlinear problems
- **Composite solvers**: Combine multiple nonlinear solvers

## Example

Newton's method for x₁² + x₂ = 1, x₁ + x₂² = 1, started from (1, 1), finds the root x₁ = x₂ = (√5 − 1)/2:

```jldoctest snes
julia> using PETSc, MPI

julia> petsclib = PETSc.petsclibs[1];

julia> PETSc.initialize(petsclib)

julia> snes = SNES(petsclib, MPI.COMM_SELF; snes_rtol = 1e-12);

julia> f = PetscVec(petsclib, zeros(2));

julia> PETSc.set_function!(snes, f) do fx, snes, x
           fx[1] = x[1]^2 + x[2] - 1
           fx[2] = x[1] + x[2]^2 - 1
       end;

julia> J = PetscMat(petsclib, zeros(2, 2));

julia> PETSc.set_snes_jacobian!(snes, J) do J, snes, x
           J[1, 1] = 2x[1]; J[1, 2] = 1
           J[2, 1] = 1;     J[2, 2] = 2x[2]
           PETSc.assemble!(J)
       end;

julia> x = PetscVec(petsclib, [1.0, 1.0]);

julia> PETSc.solve!(x, snes);

julia> x[:] ≈ fill((√5 - 1) / 2, 2)
true

julia> foreach(PETSc.destroy!, (x, J, f, snes))
```

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

The options are applied once, at the first `setup!` or `solve!`, after the DM and the callbacks are attached. A setting made in code after that overrides them. The `options_at` keyword changes that; see [When the options apply](@ref) in the KSP page.

## Setting the Nonlinear Function

Define the residual function `F(x)`. The callback comes first in the argument
list, so `do` block syntax works:

```julia
function residual!(fx, snes, x)
    # Compute F(x) and store in fx
    fx[1] = x[1]^2 + x[2] - 1
    fx[2] = x[1] + x[2]^2 - 1
end

PETSc.set_function!(residual!, snes, f_vec)
```

## Setting the Jacobian

Define the Jacobian `J = dF/dx`. `set_snes_jacobian!` rather than
`set_jacobian!`, because `set_jacobian!` belongs to `PetscDS` and the callback
occupies argument 1 here:

```julia
function jacobian!(J, snes, x)
    # Fill Jacobian matrix
    J[1, 1] = 2*x[1]
    J[1, 2] = 1.0
    J[2, 1] = 1.0
    J[2, 2] = 2*x[2]
    PETSc.assemble!(J)
end

PETSc.set_snes_jacobian!(jacobian!, snes, J, J)  # (J, P), P the preconditioner matrix
```

## Using a DM

For PDE problems, associate the SNES with a DM:

```julia
PETSc.set_dm!(snes, dm)

# Get the DM from SNES. `dm` is both the accessor and the usual variable name,
# so bind the result to something else.
d = PETSc.dm(snes)
```

## Solving

```julia
# Solve with initial guess x
PETSc.solve!(x, snes)

# Get the solution vector. It is a **borrowed** handle: it belongs to `snes`,
# and `destroy!` on it is a no-op
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
- `snes_max_funcs` - Maximum residual evaluations
- `snes_monitor` - Print residual each iteration

In code, `PETSc.tolerances(snes)` returns `(; atol, rtol, stol, max_it, max_funcs)` and `PETSc.set_tolerances!(snes; rtol = 1e-10)` changes the ones given. Under the default `options_at = :first_solve`, a value set after the first solve stays; under `:each_solve` the constructor options win at every solve.

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
Pages   = ["snes/snes.jl"]
```
