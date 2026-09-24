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
PETSc.converged_reason(ksp)  # positive when it converged
PETSc.iteration_number(ksp)

PETSc.destroy!(ksp)
```

`KSP` is the type, not a factory function: `ksp isa PETSc.KSP` holds, and
`PETSc.dm(ksp)` hands back the DM it was built on as a **borrowed** handle —
`destroy!` on it is a no-op ([naming conventions](naming.md), §3.3).

Type names are `Symbol` at the Julia API: `PETSc.set_type!(ksp, :cg)` and
`PETSc.type_name(ksp) === :cg`. The `String` spelling warns in v0.5 and is a
`MethodError` in v0.6 (§3.1).

## A nested solver

An inner solve, such as the Schur complement solve inside a preconditioner, starts from a bare `KSP` and an assembled operator. An options prefix keeps its options apart from the outer solver's:

```julia
inner = LibPETSc.KSPCreate(petsclib, comm)
PETSc.set_operators!(inner, S)              # S also builds the preconditioner; pass P to differ
PETSc.set_options_prefix!(inner, "schur_")  # reads -schur_ksp_type, -schur_pc_type, ...
PETSc.set_from_options!(inner)              # apply them now rather than at the first solve!
PETSc.solve!(y, inner, r)
PETSc.options_prefix(inner)                 # "schur_"
```

A `KSP` built on a DM computes its operator, right-hand side and initial guess from it. `PETSc.set_dm_active!(ksp, false)` keeps the DM for its geometry (multigrid needs it) and takes the rest from `set_operators!` and `solve!`; `PETSc.set_dm_active!(ksp, :rhs, false)` switches off one part.

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

A `:shell` preconditioner is written in Julia. The action writes into its first argument, like every other callback here:

```julia
p = PETSc.pc(ksp)
PETSc.set_type!(p, :shell)
PETSc.set_shell_setup!(p) do p
    # rebuild whatever apply! needs, e.g. after the operator changed
end
PETSc.set_shell_apply!(p) do y, p, x
    PETSc.with_local_array!(y, x; read = (false, true), write = (true, false)) do ya, xa
        ya .= xa ./ diagonal                 # Jacobi, by hand
    end
end
```

The closures are kept with the PETSc preconditioner, not with the wrapper `p`, so `p` can be dropped. An exception thrown inside either one comes out of the solve that ran it, as itself: a `DomainError` in `apply!` makes `ksp \ b` throw that `DomainError` ([naming conventions](naming.md), §18.4).

A `LibPETSc.PC` is also what every low-level `PC*` function takes, so anything without a high-level verb is one call away: `LibPETSc.PCFieldSplitSetType(petsclib, p, LibPETSc.PC_COMPOSITE_SCHUR)`.

## Functions

```@autodocs
Modules = [PETSc]
Pages   = ["ksp.jl", "pc.jl"]
```
