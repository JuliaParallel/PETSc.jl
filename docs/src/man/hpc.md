# Running on HPC Systems

PETSc.jl can be used on HPC clusters in two main configurations: using precompiled binaries via MPITrampoline, or by pointing to a locally installed PETSc build.

## 1. Use precompiled binaries via MPITrampoline

[MPITrampoline](https://github.com/eschnett/MPITrampoline) is an MPI wrapper layer that lets MPI-linked binaries be redirected to any system MPI at runtime. The `PETSc_jll` binaries distributed with PETSc.jl are built against MPITrampoline, which means they can be used on clusters by simply configuring `MPI.jl` to use the system MPI.

#### Install MPITrampoline 

**Steps:**
1. Configure `MPI.jl` to use the system MPI (do this once per cluster):
```julia
using MPI
MPI.install_mpiexecjl()           # optional: installs mpiexecjl wrapper
```
Or set the environment variable before starting Julia:
```bash
export JULIA_MPI_BINARY=system
export JULIA_MPI_PATH=/path/to/system/mpi   # e.g. /usr/lib/openmpi
```

2. Use PETSc.jl as normal — no changes to your script needed:
```julia
using PETSc, MPI
petsclib = PETSc.getlib(; PetscScalar = Float64, PetscInt = Int64)
PETSc.initialize(petsclib)
# ...
```

3. Launch via the cluster's MPI:
```bash
mpiexec -n 128 julia --project myScript.jl
```
This is the easiest path for most clusters and requires no custom PETSc compilation.

## Option 2: Locally installed PETSc build

If you need a PETSc build with specific options (external packages, GPU support, custom BLAS, etc.), you can point PETSc.jl directly to your local installation.
The local PETSc must be compiled as a **shared library** (`--with-shared-libraries=1`) and linked against the **same MPI** that `MPI.jl` is configured to use.

### Via environment variables (recommended for batch scripts)

Set these before starting Julia:

```bash
export JULIA_PETSC_LIBRARY=/path/to/petsc/install/lib/libpetsc.so
export JULIA_PETSC_SCALAR=Float64    # Float32 | ComplexFloat64 | ComplexFloat32
export JULIA_PETSC_INT=Int64         # Int32
```

With `JULIA_PETSC_LIBRARY` set, `PETSc_jll` is **not** loaded or precompiled, which
avoids MPI incompatibility issues common on HPC clusters.

### Via `set_petsclib` (recommended for scripts)

```julia
using PETSc, MPI

petsclib = PETSc.set_petsclib(
    "/path/to/petsc/install/lib/libpetsc.so";
    PetscScalar = Float64,
    PetscInt    = Int64
)
PETSc.initialize(petsclib)
# ... your code ...
PETSc.finalize(petsclib)
```

To prevent `PETSc_jll` from being precompiled at all (which can fail on clusters
with incompatible MPI), set:
```bash
export JULIA_PETSC_SKIP_JLL=1
```

### Typical HPC job script

```bash
#!/bin/bash
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=128

module load PETSc/3.22.0-foss-2023b

export JULIA_PETSC_SKIP_JLL=1

srun julia --project ex45.jl \
  -N 257 \
  -ksp_type cg \
  -pc_type mg \
  -pc_mg_levels 5 \
  -mg_levels_ksp_type chebyshev \
  -mg_levels_pc_type sor \
  -ksp_rtol 1e-8 \
  -ksp_view \
  -pc_mg_log \
  -log_view
```

---

## Performance and Weak Scalability

This section summarises weak scalability results for the 3D Laplacian benchmark ([`ex45.jl`](https://github.com/JuliaParallel/PETSc.jl/blob/main/examples/ex45.jl)) using a CG solver with geometric multigrid preconditioning (`-pc_type mg`). The problem size is scaled proportionally with the number of MPI ranks so that the work per rank stays constant.

### Setup

- Problem: 3D Poisson equation on a uniform grid, solved with CG + geometric MG
- Grid size per rank: ~65³ degrees of freedom
- Solver: `-ksp_type cg -pc_type mg -pc_mg_levels 4 -ksp_rtol 1e-8`
- Machine: *[HPC system name, node type, interconnect]*

### a) PETSc.jl with precompiled `PETSc_jll` binaries

*Results to be added.*

| Ranks | Grid size | Solve time (s) | KSP iterations |
|------:|----------:|---------------:|---------------:|
|     1 |      65³  |              — |              — |
|     8 |     130³  |              — |              — |
|    64 |     260³  |              — |              — |
|   512 |     520³  |              — |              — |

### b) PETSc.jl linked against a local PETSc build (system MPI)

*Results to be added.*

| Ranks | Grid size | Solve time (s) | KSP iterations |
|------:|----------:|---------------:|---------------:|
|     1 |      65³  |              — |              — |
|     8 |     130³  |              — |              — |
|    64 |     260³  |              — |              — |
|   512 |     520³  |              — |              — |

### c) Native C build (`ex45.c`)

The equivalent PETSc C example compiled natively, serving as the baseline.

*Results to be added.*

| Ranks | Grid size | Solve time (s) | KSP iterations |
|------:|----------:|---------------:|---------------:|
|     1 |      65³  |              — |              — |
|     8 |     130³  |              — |              — |
|    64 |     260³  |              — |              — |
|   512 |     520³  |              — |              — |

### Observations

*To be filled in once results are available.*
