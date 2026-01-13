# Tao (Optimization) - Low-level Interface

The Tao (Toolkit for Advanced Optimization) component provides methods for solving optimization problems, including unconstrained minimization, bound-constrained optimization, and constrained optimization.

## Overview

Tao supports various optimization algorithms:
- **Unconstrained**: Conjugate gradient (CG), limited-memory BFGS (BLMVM), Nelder-Mead
- **Bound-constrained**: Bound-constrained BFGS (BNCG, BNLS, BNTL), active-set methods
- **Constrained**: Augmented Lagrangian methods (ALMM), interior point methods
- **Least-squares**: Gauss-Newton methods, Levenberg-Marquardt
- **Complementarity**: Semismooth methods for complementarity problems
- **PDE-constrained**: Methods suitable for PDE-constrained optimization

The Tao object manages the optimization process, line searches, convergence monitoring, and provides a unified interface across different optimization algorithms.

## Basic Usage

```julia
using PETSc

# Initialize PETSc
petsclib = PETSc.getlib()

# Create a Tao object
tao = LibPETSc.TaoCreate(petsclib, MPI.COMM_SELF)

# Set the optimization algorithm (e.g., LMVM, BLMVM, NLS)
LibPETSc.TaoSetType(petsclib, tao, Base.unsafe_convert(Ptr{Int8}, "lmvm"))

# Set the objective function and gradient
# LibPETSc.TaoSetObjective(petsclib, tao, objective_function_ptr, C_NULL)
# LibPETSc.TaoSetGradient(petsclib, tao, C_NULL, gradient_function_ptr, C_NULL)

# For bound-constrained problems, set variable bounds
# LibPETSc.TaoSetVariableBounds(petsclib, tao, lower_bound_vec, upper_bound_vec)

# Set convergence tolerances
LibPETSc.TaoSetTolerances(petsclib, tao, 1e-8, 1e-8, 1e-8)

# Set maximum iterations
LibPETSc.TaoSetMaximumIterations(petsclib, tao, 1000)

# Set options from command line/options database
LibPETSc.TaoSetFromOptions(petsclib, tao)

# Set initial guess
# LibPETSc.TaoSetSolution(petsclib, tao, initial_vec)

# Solve the optimization problem
# LibPETSc.TaoSolve(petsclib, tao)

# Get solution information (values returned directly)
reason = Ref{LibPETSc.TaoConvergedReason}()
LibPETSc.TaoGetConvergedReason(petsclib, tao, reason)

iter = LibPETSc.TaoGetIterationNumber(petsclib, tao)

# Cleanup
LibPETSc.TaoDestroy(petsclib, tao)
```

## Common Workflow

1. **Create and configure Tao**: Use `TaoCreate`, `TaoSetType`
2. **Define objective**: `TaoSetObjective`, `TaoSetGradient`, optionally `TaoSetHessian`
3. **Set constraints** (if any): `TaoSetVariableBounds`, `TaoSetConstraints`
4. **Configure solver**: `TaoSetTolerances`, `TaoSetMaximumIterations`
5. **Set initial guess**: `TaoSetSolution`
6. **Solve**: `TaoSolve`
7. **Retrieve solution**: `TaoGetSolution`, `TaoGetConvergedReason`

## Optimization Algorithms

Available through `TaoSetType`:
- **Unconstrained**:
  - `TAOLMVM`: Limited-memory variable metric (quasi-Newton)
  - `TAOCG`: Conjugate gradient methods
  - `TAONM`: Nelder-Mead simplex method
  - `TAONLS`: Newton line search
  - `TAONTL`: Newton trust-region with line search

- **Bound-constrained**:
  - `TAOBLMVM`: Bound-constrained limited-memory variable metric
  - `TAOBNCG`: Bound-constrained conjugate gradient
  - `TAOBQNLS`: Bound-constrained quasi-Newton line search
  - `TAOBNTL`: Bound-constrained Newton trust-region
  - `TAOTRON`: Trust-region Newton method

- **Constrained**:
  - `TAOALMM`: Augmented Lagrangian multiplier method
  - `TAOIPM`: Interior point method
  - `TAOPDIPM`: Primal-dual interior point method

- **Least-squares**:
  - `TAOPOUNDERS`: POUNDERs model-based method
  - `TAOBRGN`: Bounded regularized Gauss-Newton

- **Complementarity**:
  - `TAOSSLS`: Semismooth least squares
  - `TAOASLS`: Active-set least squares

## Convergence Criteria

Tao monitors several convergence criteria:
- **Gradient tolerance**: `||∇f|| < gatol` or `||∇f||/||f|| < grtol`
- **Function tolerance**: `|f - f_prev| < fatol` or `|f - f_prev|/|f| < frtol`
- **Step tolerance**: `||x - x_prev|| < steptol`

Set using `TaoSetTolerances(tao, gatol, grtol, gttol)`.

## Hessian Options

For second-order methods:
- **Exact Hessian**: Provide via `TaoSetHessian`
- **Finite-difference approximation**: Use matrix-free approach
- **Quasi-Newton approximation**: LMVM methods build approximation automatically

## Function Reference

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/Tao_wrappers.jl"]
Order   = [:function]
```

## Tao Add-ons

Additional Tao utilities and helper functions:

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/Tao_addons_wrappers.jl"]
Order   = [:function]
```
