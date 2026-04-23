/*
  Laplacian in 3D with manufactured solution.

  PDE:  -Laplacian u = 12*pi^2 * sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z)
  on (0,1)^3 with Dirichlet BCs:
        u = sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z)  on all boundaries

  Exact solution: u = sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z)

  Matrix is scaled by h^3 (= Hx*Hy*Hz) so entries are O(h) rather than O(1/h^2),
  matching the scaling used in ex45.jl.

  Grid size is controlled by -Nx, -Ny, -Nz (matching ex45.jl convention).

  Output is formatted to match ex45.jl so that parse_scaling.jl works unchanged.
*/

static char help[] = "Solves 3D Laplacian (manufactured solution) using multigrid.\n\n";

#include <petscksp.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <math.h>

/* Context passed to callbacks */
typedef struct {
  PetscInt Nx, Ny, Nz;
} AppCtx;

extern PetscErrorCode ComputeMatrix(KSP, Mat, Mat, void *);
extern PetscErrorCode ComputeRHS(KSP, Vec, void *);
extern PetscErrorCode ComputeInitialGuess(KSP, Vec, void *);

static inline PetscReal exact_sol(PetscReal x, PetscReal y, PetscReal z)
{
  return sin(2.0 * PETSC_PI * x) * sin(2.0 * PETSC_PI * y) * sin(2.0 * PETSC_PI * z);
}

static inline PetscReal forcing(PetscReal x, PetscReal y, PetscReal z)
{
  return 12.0 * PETSC_PI * PETSC_PI *
         sin(2.0 * PETSC_PI * x) * sin(2.0 * PETSC_PI * y) * sin(2.0 * PETSC_PI * z);
}

int main(int argc, char **argv)
{
  KSP        ksp;
  DM         da;
  Vec        x, b, r;
  Mat        A;
  AppCtx     ctx;
  PetscReal  norm, solve_time, L2, maxerr;
  PetscInt   niter;
  PetscLogDouble t0, t1;

  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, (char *)0, help));

  /* Default grid size matches ex45.jl default */
  ctx.Nx = 7; ctx.Ny = 7; ctx.Nz = 7;
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-Nx", &ctx.Nx, NULL));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-Ny", &ctx.Ny, NULL));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-Nz", &ctx.Nz, NULL));

  /* Print grid info matching ex45.jl output */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Solving on %d × %d × %d grid\n",
                        (int)ctx.Nx, (int)ctx.Ny, (int)ctx.Nz));

  PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
  PetscCall(DMDACreate3d(PETSC_COMM_WORLD,
                         DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                         DMDA_STENCIL_STAR,
                         ctx.Nx, ctx.Ny, ctx.Nz,
                         PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                         1, 1, NULL, NULL, NULL, &da));
  PetscCall(DMSetFromOptions(da));
  PetscCall(DMSetUp(da));
  PetscCall(KSPSetDM(ksp, da));
  PetscCall(KSPSetComputeInitialGuess(ksp, ComputeInitialGuess, &ctx));
  PetscCall(KSPSetComputeRHS(ksp, ComputeRHS, &ctx));
  PetscCall(KSPSetComputeOperators(ksp, ComputeMatrix, &ctx));
  PetscCall(DMDestroy(&da));
  PetscCall(KSPSetFromOptions(ksp));

  /* Time the solve */
  PetscCall(PetscTime(&t0));
  PetscCall(KSPSolve(ksp, NULL, NULL));
  PetscCall(PetscTime(&t1));
  solve_time = (PetscReal)(t1 - t0);

  PetscCall(KSPGetSolution(ksp, &x));
  PetscCall(KSPGetRhs(ksp, &b));
  PetscCall(KSPGetIterationNumber(ksp, &niter));
  PetscCall(KSPGetResidualNorm(ksp, &norm));

  /* Compute explicit residual norm */
  PetscCall(VecDuplicate(b, &r));
  PetscCall(KSPGetOperators(ksp, &A, NULL));
  PetscCall(MatMult(A, x, r));
  PetscCall(VecAXPY(r, -1.0, b));
  PetscReal resnorm;
  PetscCall(VecNorm(r, NORM_2, &resnorm));

  /* Compute L2 and max error against exact solution */
  {
    DM          dm;
    DMDALocalInfo info;
    Vec         x_local;
    PetscScalar ***xu;
    PetscReal   local_err2 = 0.0, local_max = 0.0;
    PetscReal   Hx, Hy, Hz, vol;
    PetscInt    i, j, k;

    PetscCall(KSPGetDM(ksp, &dm));
    PetscCall(DMDAGetLocalInfo(dm, &info));
    Hx  = 1.0 / (PetscReal)(info.mx - 1);
    Hy  = 1.0 / (PetscReal)(info.my - 1);
    Hz  = 1.0 / (PetscReal)(info.mz - 1);
    vol = Hx * Hy * Hz;

    PetscCall(DMGetLocalVector(dm, &x_local));
    PetscCall(DMGlobalToLocalBegin(dm, x, INSERT_VALUES, x_local));
    PetscCall(DMGlobalToLocalEnd(dm, x, INSERT_VALUES, x_local));
    PetscCall(DMDAVecGetArrayRead(dm, x_local, &xu));

    for (k = info.zs; k < info.zs + info.zm; k++) {
      for (j = info.ys; j < info.ys + info.ym; j++) {
        for (i = info.xs; i < info.xs + info.xm; i++) {
          PetscReal xc   = i * Hx;
          PetscReal yc   = j * Hy;
          PetscReal zc   = k * Hz;
          PetscReal u_ex = exact_sol(xc, yc, zc);
          PetscReal err  = (PetscReal)xu[k][j][i] - u_ex;
          local_err2 += err * err * vol;
          if (PetscAbsReal(err) > local_max) local_max = PetscAbsReal(err);
        }
      }
    }

    PetscCall(DMDAVecRestoreArrayRead(dm, x_local, &xu));
    PetscCall(DMRestoreLocalVector(dm, &x_local));

    PetscCall(MPI_Allreduce(&local_err2, &L2,     1, MPIU_REAL, MPIU_SUM, PETSC_COMM_WORLD));
    PetscCall(MPI_Allreduce(&local_max,  &maxerr, 1, MPIU_REAL, MPIU_MAX, PETSC_COMM_WORLD));
    L2 = PetscSqrtReal(L2);
  }

  /* Output exactly matching ex45.jl format for parse_scaling.jl */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "L2 error: %e\n",       (double)L2));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Max error: %e\n",      (double)maxerr));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Residual norm %g\n",   (double)resnorm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Solve time: %f seconds\n", (double)solve_time));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Final residual norm %g\n", (double)norm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Fine grid resolution: %d × %d × %d\n",
                        (int)ctx.Nx, (int)ctx.Ny, (int)ctx.Nz));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Outer KSP iterations: %d\n", (int)niter));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Solve time: %f seconds\n",   (double)solve_time));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "L2 error: %e\n",             (double)L2));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Max error: %e\n",            (double)maxerr));

  PetscCall(VecDestroy(&r));
  PetscCall(KSPDestroy(&ksp));
  PetscCall(PetscFinalize());
  return 0;
}

PetscErrorCode ComputeInitialGuess(KSP ksp, Vec b, void *ctx)
{
  PetscFunctionBeginUser;
  PetscCall(VecSet(b, 0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode ComputeRHS(KSP ksp, Vec b, void *ctx)
{
  AppCtx            *user = (AppCtx *)ctx;
  DM                 dm;
  DMDALocalInfo      info;
  PetscScalar      ***barray;
  PetscInt           i, j, k;
  PetscReal          Hx, Hy, Hz, scale;

  PetscFunctionBeginUser;
  PetscCall(KSPGetDM(ksp, &dm));
  PetscCall(DMDAGetLocalInfo(dm, &info));

  Hx    = 1.0 / (PetscReal)(info.mx - 1);
  Hy    = 1.0 / (PetscReal)(info.my - 1);
  Hz    = 1.0 / (PetscReal)(info.mz - 1);
  scale = Hx * Hy * Hz;   /* same h^3 scaling as ex45.jl */

  PetscCall(DMDAVecGetArray(dm, b, &barray));
  for (k = info.zs; k < info.zs + info.zm; k++) {
    for (j = info.ys; j < info.ys + info.ym; j++) {
      for (i = info.xs; i < info.xs + info.xm; i++) {
        PetscReal x = i * Hx;
        PetscReal y = j * Hy;
        PetscReal z = k * Hz;
        if (i == 0 || j == 0 || k == 0 ||
            i == info.mx - 1 || j == info.my - 1 || k == info.mz - 1) {
          /* Dirichlet BC: identity row -> RHS = exact value */
          barray[k][j][i] = exact_sol(x, y, z);
        } else {
          /* Interior: scale forcing by h^3 */
          barray[k][j][i] = forcing(x, y, z) * scale;
        }
      }
    }
  }
  PetscCall(DMDAVecRestoreArray(dm, b, &barray));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode ComputeMatrix(KSP ksp, Mat jac, Mat B, void *ctx)
{
  DM            da;
  DMDALocalInfo info;
  PetscInt      i, j, k;
  PetscReal     Hx, Hy, Hz, scale;
  PetscReal     cx, cy, cz, diag;
  PetscScalar   v[7];
  MatStencil    row, col[7];

  PetscFunctionBeginUser;
  PetscCall(KSPGetDM(ksp, &da));
  PetscCall(DMDAGetLocalInfo(da, &info));

  Hx    = 1.0 / (PetscReal)(info.mx - 1);
  Hy    = 1.0 / (PetscReal)(info.my - 1);
  Hz    = 1.0 / (PetscReal)(info.mz - 1);
  scale = Hx * Hy * Hz;   /* h^3 scaling matching ex45.jl */

  /* Scaled FD coefficients: scale * (1/h^2) = h^3/h^2 = h */
  cx   = scale / (Hx * Hx);   /* = Hy*Hz/Hx */
  cy   = scale / (Hy * Hy);   /* = Hx*Hz/Hy */
  cz   = scale / (Hz * Hz);   /* = Hx*Hy/Hz */
  diag = 2.0 * (cx + cy + cz);

  for (k = info.zs; k < info.zs + info.zm; k++) {
    for (j = info.ys; j < info.ys + info.ym; j++) {
      for (i = info.xs; i < info.xs + info.xm; i++) {
        row.i = i; row.j = j; row.k = k;
        if (i == 0 || j == 0 || k == 0 ||
            i == info.mx - 1 || j == info.my - 1 || k == info.mz - 1) {
          /* Identity row for Dirichlet BC */
          v[0] = 1.0;
          PetscCall(MatSetValuesStencil(B, 1, &row, 1, &row, v, INSERT_VALUES));
        } else {
          v[0] = -cz; col[0].i = i;   col[0].j = j;   col[0].k = k-1;
          v[1] = -cy; col[1].i = i;   col[1].j = j-1; col[1].k = k;
          v[2] = -cx; col[2].i = i-1; col[2].j = j;   col[2].k = k;
          v[3] =  diag; col[3].i = i; col[3].j = j;   col[3].k = k;
          v[4] = -cx; col[4].i = i+1; col[4].j = j;   col[4].k = k;
          v[5] = -cy; col[5].i = i;   col[5].j = j+1; col[5].k = k;
          v[6] = -cz; col[6].i = i;   col[6].j = j;   col[6].k = k+1;
          PetscCall(MatSetValuesStencil(B, 1, &row, 7, col, v, INSERT_VALUES));
        }
      }
    }
  }
  PetscCall(MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY));
  PetscFunctionReturn(PETSC_SUCCESS);
}