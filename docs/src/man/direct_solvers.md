# Direct solvers

`PETSc_jll` links several sparse direct solvers. `pc_type = "lu"` factors the matrix and `pc_factor_mat_solver_type` picks the package:

```julia
ksp = PETSc.KSP(A; ksp_type = "preonly", pc_type = "lu", pc_factor_mat_solver_type = "mumps")
```

Without `pc_factor_mat_solver_type` you get PETSc's own LU, which runs on one rank and doesn't pivot. For a saddle point, such as a Stokes system, that can give a wrong answer reported as a converged solve.

## Which solver to use

| Solver | `pc_factor_mat_solver_type` | Ranks | SPD | Saddle point |
|---|---|---|---|---|
| MUMPS | `mumps` | any | right | right |
| UMFPACK | `umfpack` | 1 | right | right |
| KLU | `klu` | 1 | right | right |
| CHOLMOD | `cholmod`, with `pc_type = "cholesky"` | 1 | right | n/a |
| SuperLU_DIST | `superlu_dist` | any | right | **wrong on one of two tests, reported as converged** |
| PETSc | `petsc` (the default) | 1 | right | **loses digits, reported as converged** |

On more than one rank, a solver that runs on one rank only stops with an error.

BoomerAMG (`pc_type = "hypre"`) isn't a direct solver, but it is in the same build: inside GMRES it solves the SPD test on 1 and 2 ranks.

For a saddle point, use MUMPS, or UMFPACK in serial. A fieldsplit that treats velocity and pressure separately avoids factoring the saddle point at all; [Stokes flow on a DMStag](dmstag_stokes.md) builds one.

## MUMPS on more than one rank

With MUMPS_jll 5.9.4 and earlier, a MUMPS solve on more than one rank on Linux crashes after the factorization with `double free or corruption`, or hangs. The libraries were loaded so that MUMPS didn't recognize MPICH's `MPI_IN_PLACE`. [Yggdrasil#15135](https://github.com/JuliaPackaging/Yggdrasil/pull/15135) fixes it in MUMPS_jll 5.9.5 and MPICH_jll 5.0.2+1. A serial MUMPS solve is not affected.

## Wrong answers reported as converged

A direct solver used with `ksp_type = "preonly"` always reports success: there is no iteration to check.

- **PETSc's LU** doesn't pivot. On the stiff saddle point below it gives a relative error of 10⁻⁸ and a relative residual ‖b − Ax‖/‖b‖ of 3·10⁻², where UMFPACK gives 10⁻¹³ and 10⁻⁷. On a 2D Stokes model at a viscosity contrast of 10⁶ we have seen velocities off by a factor of three, with SNES reporting convergence.
- **SuperLU_DIST** gives an error of 55 (serial) and 8 (2 ranks) on the milder saddle point below, with a residual of 2·10⁶. On the stiff saddle point it is right. We haven't found the cause.

GMRES around the factorization doesn't catch this: it stops on the preconditioned residual, which the inaccurate factorization keeps small, also with `ksp_norm_type = "unpreconditioned"`. FGMRES preconditions from the right and stops on the true residual. Around PETSc's LU it converges in 4 iterations to an error of 10⁻¹³; around SuperLU_DIST it doesn't converge, so the solve reports the failure instead of returning a wrong answer:

```julia
ksp = PETSc.KSP(A; ksp_type = "fgmres", ksp_rtol = 1e-10, pc_type = "lu")
```

Otherwise, compute ‖b − Ax‖/‖b‖ after the solve.

## How the table was made

[`scripts/direct_solvers.jl`](https://github.com/JuliaParallel/PETSc.jl/blob/main/scripts/direct_solvers.jl) runs each solver in its own `mpiexec` process, on 1 and 2 ranks, and compares the solution with Julia's sparse LU. "Right" means a relative error below 10⁻⁸. It solves three systems, each with `ksp_type = "preonly"` and with GMRES around the factorization:

- **SPD:** −∇·(η∇u) on a 24² grid, Dirichlet, η = 10⁴ in a central square and 1 elsewhere.
- **Stiff saddle point:** [A Bᵀ; B −εI] on a 24² grid, A the same operator for two velocity components with a contrast of 10⁶, B a divergence, ε = 10⁻¹².
- **Mild saddle point:** the same on a 32² grid with a contrast of 10⁴ and ε = 10⁻⁶.

The table is for Linux x86_64 with MPICH, PETSc_jll 3.25.4, and MUMPS_jll and MPICH_jll with the fix above. Other platforms and OpenMPI are untested. To check your own installation:

```
julia --project scripts/direct_solvers.jl
```

prints the full table, and `--check` fails if it differs from the one above.
