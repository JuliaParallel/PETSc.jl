# KSP

The KSP (Krylov Subspace Methods) module provides iterative linear solvers for systems of the form `Ax = b`. PETSc offers a wide variety of Krylov methods and preconditioners.

## Overview

KSP provides:
- **Krylov methods**: GMRES, CG, BiCGStab, and many more
- **Preconditioners**: Jacobi, ILU, multigrid, direct solvers, etc.
- **Runtime configuration**: Choose methods via command-line options
- **Convergence monitoring**: Built-in residual tracking

## Creating a KSP Solver

### From a Matrix

```julia
# Basic creation with default options
ksp = KSP(A)

# With preconditioner matrix P (for different preconditioning)
ksp = KSP(A, P)

# With options
ksp = KSP(A; 
    ksp_type = "gmres",
    pc_type = "ilu",
    ksp_rtol = 1e-8
)
```

### From a DM

```julia
# Create KSP associated with a DM (for multigrid, etc.)
ksp = KSP(dm; 
    ksp_type = "cg",
    pc_type = "mg"
)
```

### From a Sparse Matrix

```julia
# Directly from Julia SparseMatrixCSC
using SparseArrays
S = sprand(100, 100, 0.1) + 10I
ksp = KSP(petsclib, MPI.COMM_SELF, S)
```

## Solving

```julia
# Solve Ax = b
solve!(x, ksp, b)

# Or allocate solution vector
x = solve(ksp, b)
```

## Common Solver/Preconditioner Options

### Krylov Methods (`ksp_type`)
- `cg` - Conjugate Gradient (symmetric positive definite)
- `gmres` - Generalized Minimum Residual
- `bicgstab` - BiConjugate Gradient Stabilized
- `richardson` - Richardson iteration
- `preonly` - Apply preconditioner only (for direct solvers)

### Preconditioners (`pc_type`)
- `jacobi` - Diagonal scaling
- `ilu` - Incomplete LU factorization
- `lu` - Direct LU factorization
- `mg` - Geometric multigrid
- `gamg` - Algebraic multigrid
- `none` - No preconditioning

### Convergence Options
- `ksp_rtol` - Relative tolerance (default: 1e-5)
- `ksp_atol` - Absolute tolerance
- `ksp_max_it` - Maximum iterations
- `ksp_monitor` - Print residual each iteration

## Example: Multigrid Solver

```julia
ksp = KSP(dm;
    ksp_type = "cg",
    pc_type = "mg",
    pc_mg_levels = 4,
    pc_mg_galerkin = true,
    mg_levels_ksp_type = "richardson",
    mg_levels_pc_type = "jacobi",
    mg_coarse_pc_type = "lu"
)
```

## Functions

```@autodocs
Modules = [PETSc]
Pages   = ["ksp.jl"]
```
