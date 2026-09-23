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
# Basic creation with default options. The communicator comes from `A`
ksp = KSP(A)

# With the preconditioner construction matrix P
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
# Solve Ax = b, writing into x. The written vector comes first
PETSc.solve!(x, ksp, b)

# `ldiv!` and `\` are the LinearAlgebra spellings of the same call
using LinearAlgebra
ldiv!(x, ksp, b)
x = ksp \ b

# What PETSc did
PETSc.type_name(ksp)         # :gmres, a Symbol
PETSc.converged_reason(ksp)

PETSc.destroy!(ksp)
```

`KSP` is the type, not a factory function: `ksp isa PETSc.KSP` holds, and
`PETSc.dm(ksp)` hands back the DM it was built on as a **borrowed** handle —
`destroy!` on it is a no-op ([naming conventions](naming.md), §3.3).

Type names are `Symbol` at the Julia API: `PETSc.set_type!(ksp, :cg)` and
`PETSc.type_name(ksp) === :cg`. The `String` spelling warns in v0.5 and is a
`MethodError` in v0.6 (§3.1).

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

## The preconditioner

`PETSc.pc(ksp)` hands back the preconditioner as a borrowed `LibPETSc.PC`. It takes the same `set_type!`/`type_name` pair as the solver, and is how a split preconditioner gets its index sets, which options alone cannot supply:

```julia
p = PETSc.pc(ksp)                         # not pc = pc(ksp), see naming.md §3.2
PETSc.set_type!(p, :fieldsplit)
PETSc.set_fieldsplit_is!(p, "u", is_u)    # rows of the first split, 0-based
PETSc.set_fieldsplit_is!(p, "p", is_p)    # options prefix -fieldsplit_p_
PETSc.type_name(p)                        # :fieldsplit
```

A `LibPETSc.PC` is also what every low-level `PC*` function takes, so anything without a high-level verb is one call away: `LibPETSc.PCFieldSplitSetType(petsclib, p, LibPETSc.PC_COMPOSITE_SCHUR)`.

## Functions

```@autodocs
Modules = [PETSc]
Pages   = ["ksp.jl", "pc.jl"]
```
