# PC (Preconditioners) - Low-level Interface

The PC (Preconditioner) component provides methods for preconditioning linear systems to accelerate the convergence of iterative solvers. Preconditioners are essential for efficiently solving large sparse linear systems with KSP.

## Overview

Preconditioners transform the linear system $Ax = b$ into an equivalent system that is easier to solve iteratively. PETSc provides a wide variety of preconditioners:

- **Basic methods**: Jacobi, block Jacobi, SOR, ILU, ICC
- **Multigrid**: Algebraic multigrid (GAMG, HYPRE BoomerAMG), geometric multigrid (MG)
- **Domain decomposition**: Additive Schwarz (ASM), block Jacobi with local solves
- **Direct solvers**: LU, Cholesky (using external packages like MUMPS, SuperLU, UMFPACK)
- **Sparse approximations**: ILU(k), ICC(k), approximate inverses
- **Physics-based**: Field-split, composite preconditioners for coupled systems
- **Matrix-free**: Shell preconditioners for matrix-free implementations

The PC object can be used standalone or automatically by KSP during the solution process.

## Basic Usage

```julia
using PETSc, MPI

# Initialize MPI and PETSc
MPI.Init()
petsclib = PETSc.getlib()
PETSc.initialize(petsclib)

# Create a PC object
pc = Ref{LibPETSc.PC}()
LibPETSc.PCCreate(petsclib, LibPETSc.PETSC_COMM_SELF, pc)

# Set the preconditioner type
LibPETSc.PCSetType(petsclib, pc[], "ilu")  # String convenience wrapper

# Set the operator matrix
# LibPETSc.PCSetOperators(petsclib, pc[], A, A)

# Configure preconditioner-specific options
# For ILU: set fill level
# LibPETSc.PCFactorSetLevels(petsclib, pc[], 2)

# Set options from command line/options database
LibPETSc.PCSetFromOptions(petsclib, pc[])

# Set up the preconditioner
LibPETSc.PCSetUp(petsclib, pc[])

# Apply the preconditioner: y = P^{-1} x
# LibPETSc.PCApply(petsclib, pc[], x_vec, y_vec)

# Cleanup
LibPETSc.PCDestroy(petsclib, pc)

# Finalize PETSc and MPI
PETSc.finalize(petsclib)
MPI.Finalize()
```

## Integration with KSP

Preconditioners are typically used with KSP solvers:

```julia
# Create KSP and get its PC
ksp = Ref{LibPETSc.KSP}()
LibPETSc.KSPCreate(petsclib, LibPETSc.PETSC_COMM_SELF, ksp)

pc = Ref{LibPETSc.PC}()
LibPETSc.KSPGetPC(petsclib, ksp[], pc)

# Configure the preconditioner
LibPETSc.PCSetType(petsclib, pc[], "gamg")  # String convenience wrapper

# For GAMG, set additional options
LibPETSc.PCGAMGSetType(petsclib, pc[], LibPETSc.PCGAMGAGG)
LibPETSc.PCGAMGSetNlevels(petsclib, pc[], 10)
```

## Common Preconditioner Types

Available through `PCSetType`:

### Direct Solvers
- `PCLU`: LU factorization
- `PCCHOLESKY`: Cholesky factorization (symmetric positive definite)
- `PCQR`: QR factorization

### Incomplete Factorizations
- `PCILU`: Incomplete LU factorization
- `PCICC`: Incomplete Cholesky factorization
- Configure fill with `PCFactorSetLevels`

### Iterative/Relaxation Methods
- `PCJACOBI`: Jacobi (diagonal scaling)
- `PCSOR`: Successive over-relaxation
- `PCPBJACOBI`: Point-block Jacobi
- `PCBJACOBI`: Block Jacobi with local solves

### Multigrid Methods
- `PCMG`: Geometric multigrid
- `PCGAMG`: Algebraic multigrid (PETSc's built-in)
- `PCHYPRE`: HYPRE preconditioners (BoomerAMG, ParaSails, etc.)
- `PCML`: ML algebraic multigrid (Trilinos)

### Domain Decomposition
- `PCASM`: Additive Schwarz method
- `PCGASM`: Generalized additive Schwarz
- `PCBDDC`: Balancing domain decomposition by constraints

### Physics-Based
- `PCFIELDSPLIT`: Field-split for coupled systems
- `PCLSC`: Least-squares commutators (for saddle-point problems)
- `PCCOMPOSITE`: Combine multiple preconditioners

### Special Purpose
- `PCSHELL`: User-defined shell preconditioner
- `PCNONE`: No preconditioning (identity)
- `PCKSP`: Use another KSP as preconditioner
- `PCREDISTRIBUTE`: Redistribute matrix for better load balancing

## Multigrid Configuration

For geometric multigrid (PCMG):
```julia
LibPETSc.PCSetType(petsclib, pc[], LibPETSc.PCMG)
LibPETSc.PCMGSetLevels(petsclib, pc[], nlevels, C_NULL)
# Set up grid hierarchy, smoothers, coarse solver
```

For algebraic multigrid (PCGAMG):
```julia
LibPETSc.PCSetType(petsclib, pc[], LibPETSc.PCGAMG)
LibPETSc.PCGAMGSetType(petsclib, pc[], LibPETSc.PCGAMGAGG)  # Aggregation
LibPETSc.PCGAMGSetNSmooths(petsclib, pc[], 1)
LibPETSc.PCGAMGSetThreshold(petsclib, pc[], [0.0], 1)
```

## Field Split for Coupled Systems

For systems with multiple fields (e.g., velocity-pressure):
```julia
LibPETSc.PCSetType(petsclib, pc[], LibPETSc.PCFIELDSPLIT)
LibPETSc.PCFieldSplitSetType(petsclib, pc[], LibPETSc.PC_COMPOSITE_SCHUR)

# Define fields using index sets
# LibPETSc.PCFieldSplitSetIS(petsclib, pc[], "velocity", velocity_is)
# LibPETSc.PCFieldSplitSetIS(petsclib, pc[], "pressure", pressure_is)
```

## External Solver Packages

PETSc can use external direct solver packages:
- **MUMPS**: Parallel direct solver
- **SuperLU**: Sparse direct solver
- **UMFPACK**: Unsymmetric multifrontal solver
- **Pardiso**: Intel MKL parallel direct solver
- **STRUMPACK**: Structured sparse solver

Set via `PCFactorSetMatSolverType` after choosing `PCLU` or `PCCHOLESKY`.

## Performance Considerations

- **Jacobi/Block Jacobi**: Fast to apply, limited effectiveness
- **ILU/ICC**: Good for moderately difficult problems, configure fill with `PCFactorSetLevels`
- **AMG**: Excellent for elliptic PDEs, automatic coarse grid construction
- **Direct solvers**: Robust but memory-intensive for large 3D problems
- **Field split**: Essential for coupled multi-physics problems

## Function Reference

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PC_wrappers.jl"]
Order   = [:function]
```
