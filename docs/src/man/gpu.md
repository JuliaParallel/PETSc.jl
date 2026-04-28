# GPU Support (CUDA + KernelAbstractions)

Julia has outstanding support for GPUs as it compiles machine code for the particular devices. Importantly, all modern GPUs are supported, which implies that it is quite straightforward to write GPU kernels in Julia, for example using packages such as [KernelAbstractions](https://github.com/JuliaGPU/KernelAbstractions.jl).

PETSc has also added GPU support in recent years, and PETSc vector and matrix objects, along with many of the solvers, can be moved to the GPU.

GPU support in PETSc.jl requires a **locally built PETSc** with CUDA or HIP enabled — the precompiled `PETSc_jll` binaries do not include GPU support. See [Installation](@ref) for instructions on pointing PETSc.jl at a local library.
The examples below are given for CUDA. Doing this on AMD machines (HIP) will likely work the same but will require a specific extension to be added.

## Prerequisites

1. **PETSc built with CUDA** — configure with `--with-cuda=1` (and the matching `--with-cuda-dir`). Confirm with:
   ```bash
   grep -i cuda $PETSC_DIR/$PETSC_ARCH/include/petscconf.h
   # should show: #define PETSC_HAVE_CUDA 1
   ```

2. **CUDA.jl** in your Julia environment:
   ```julia
   pkg> add CUDA
   ```

3. **KernelAbstractions.jl** if you want to write portable GPU/CPU kernels:
   ```julia
   pkg> add KernelAbstractions
   ```

## How it works

When CUDA.jl is loaded alongside PETSc.jl, the `PETScCUDAExt` extension is activated automatically. It registers CUDA-aware implementations for the functions below via Julia's package extension mechanism — no extra configuration is needed.

PETSc manages where vector data lives (host or device). The extension inspects the `PetscMemType` returned by `VecGetArray*AndMemType` calls and either wraps the device pointer as a `CuArray` (zero-copy) or allocates a bounce buffer if the data needs to move between host and device.

## Public API

### `withlocalarray_device!`

Callback-based access to the underlying array of one or more Vecs, returning a `CuArray` when the data is on the device:

```julia
withlocalarray_device!(f!, vecs...; read=true, write=true)
```

```julia
using PETSc, CUDA, KernelAbstractions

withlocalarray_device!(my_vec; read=false, write=true) do arr
    # arr is a CuArray if the Vec lives on the GPU, plain Array otherwise
    fill!(arr, 42)
end
```

For multiple Vecs, pass keyword tuples to control read/write access per Vec:

```julia
withlocalarray_device!(
    (x_vec, f_vec);
    read  = (true,  false),
    write = (false, true),
) do x_arr, f_arr
    my_kernel!(backend)(f_arr, x_arr; ndrange = length(f_arr))
    KernelAbstractions.synchronize(backend)
end
```

### `get_petsc_arrays` / `restore_petsc_arrays`

Lower-level paired get/restore for the residual function pattern, where you need both a global output Vec and a local (ghost-padded) input Vec:

```julia
fx, lx, fx_arr, lx_arr, fx_bounce = get_petsc_arrays(petsclib, g_fx, l_x)
# launch kernel writing into fx, reading from lx
restore_petsc_arrays(petsclib, g_fx, l_x, fx, lx, fx_arr, lx_arr, fx_bounce)
```

- When both Vecs are on the GPU, `fx` and `lx` are zero-copy `CuArray` wrappers.
- When `l_x` is host-resident (e.g. on a coarser MG level), the data is copied host→device before the kernel and the result is copied device→host by `restore_petsc_arrays`.
- On a CPU-only path (CUDA.jl not loaded, or all Vecs on host), `fx`/`lx` are plain `Array`s with no copies.

## Writing portable kernels with KernelAbstractions

Select the backend at the top of your script based on the `useCUDA` flag:

```julia
using KernelAbstractions
using CUDA
import CUDA: CuArray, CuPtr, unsafe_wrap

const backend = CUDABackend()   # or CPU() for a CPU run
```

Write kernels with `@kernel` so the same code runs on both backends:

```julia
@kernel function my_kernel!(out, inp)
    i = @index(Global)
    out[i] = inp[i] * 2
end

# launch:
my_kernel!(backend, 256)(out_arr, inp_arr; ndrange = length(out_arr))
KernelAbstractions.synchronize(backend)
```

## Example

[`examples/ex19.jl`](https://github.com/JuliaParallel/PETSc.jl/blob/main/examples/ex19.jl) is a full 2D driven-cavity example (velocity–vorticity–temperature) that demonstrates:

- Switching between CPU and GPU with a single `useCUDA` flag.
- FD coloring-based Jacobian assembly running entirely on-device.
- `get_petsc_arrays` / `restore_petsc_arrays` in the residual callback.
- Multigrid preconditioning with coarser levels falling back to a CPU Jacobian.

To run it on a GPU:

1. Set `const useCUDA = true` near the top of `ex19.jl`.
2. Ensure your local PETSc build has CUDA support and is linked via `PETSc.set_library!`.
3. Launch:
   ```bash
   julia --project ex19.jl -da_grid_x 256 -da_grid_y 256 \
       -pc_type mg -pc_mg_levels 4 \
       -mg_levels_ksp_type chebyshev -mg_levels_pc_type jacobi \
       -snes_monitor -ksp_monitor
   ```

> [!NOTE]
> The `PetscInt` type in your `LocalPreferences.toml` must match the PETSc build.
> Check with `grep sizeof_PetscInt $PETSC_DIR/$PETSC_ARCH/include/petscconf.h`.


## Performance 

We have checked the performance of `examples/ex19.jl` by running it on GPU and on 1 or 32 CPU's of a Grace-Hopper 200 machine, using the following options:

```bash
$mpiexec -n 1 julia --project ex19.jl -snes_monitor -snes_converged_reason -snes_monitor   -pc_type mg  -mg_levels_ksp_type chebyshev -mg_levels_pc_type jacobi -ksp_monitor -log_view  -mg_levels_esteig_ksp_max_it 20  -mg_levels_esteig_ksp_max_it 20 -mg_levels_ksp_chebyshev_esteig 0,0.1,0,1.1 -mg_levels_ksp_max_it 3  -da_grid_x 513 -da_grid_y 513 -pc_mg_levels 6 
```
We performed tests in which we doubled the resolution in `x` and `y` while increasing the number of multigrid levels such that the coarse grid level has the same size in all cases.

The results of the full solve include booting up Julia and computing the coloring pattern (which doesn't scale well).

Below we report results for the inner solve itself (KSPSolve) on a GPU (with 1 MPI rank, which is a requirement of PETSc), and compare that with the results on CPU on 1 core and on 32 cores


**KSPSolve**:
| Resolution | GPU time (s) | GPU (GFlop/s) | CPU-1 time (s) | CPU-1 (GFlop/s) | CPU-32 time (s) | CPU-32 (GFlop/s) |
|---|---|---|---|---|---|---|
| 513²   |  0.116 | 144.3 |   4.337 |  4.1 |  0.297 | 61.0 |
| 1025²  |  0.299 | 249.5 |  19.10  |  3.9 |  1.196 | 61.7 |
| 2049²  |  1.118 | 295.5 |  89.76  |  3.6 |  5.757 | 56.6 |
| 4097²  |  4.540 | 312.4 | 422.2   |  3.3 | 28.30  | 49.4 |


**SNESSolve**:
| Resolution | GPU time (s) | GPU (GFlop/s) | CPU-1 time (s) | CPU-1 (GFlop/s) | CPU-32 time (s) | CPU-32 (GFlop/s) |
|---|---|---|---|---|---|---|
| 513²   |  3.645 |  6.7 |   8.041 | 3.2 |  2.118 | 12.3 |
| 1025²  |  5.127 | 20.5 |  32.75  | 3.2 |  3.698 | 28.2 |
| 2049²  | 11.37  | 39.7 | 144.1   | 3.1 | 10.57  | 42.3 |
| 4097²  | 36.88  | 47.7 | 658.3   | 2.9 | 43.58  | 43.2 |

From this it is clear that the `KSPSolve` itself is very efficient on the GPU (and clearly beats the CPU), but that there is quite some overhead when we compare it with `SNESSolve` where this difference is not so large anymore.


Lets have a look in detail on whethervthe example actually runs on the GPU:

***GPU Utilisation Evidence — ex19.jl on GH200 (SM90)***

All GPU runs use a CUDA-enabled PETSc build (cuSPARSE, cuBLAS) with KernelAbstractions.jl providing the residual kernel. The PETSc `GPU %F` column (fraction of flops executed on GPU) and the host↔device transfer logs provide direct evidence of efficient GPU utilisation.

*1. Flop fraction on GPU*

Nearly all floating point work executes on the GPU across all resolutions and all major solver phases:

| Event | GPU %F |
|---|---|
| KSPSolve | 100% |
| PCApply | 100% |
| PCSetUp | 98–99% |
| SNESSolve | 99–100% |

At the kernel level, all key GMRES and multigrid operations run fully on the GPU: `MatMult`, `MatResidual`, `VecMAXPY`, `VecMDot`, `VecAXPBYCZ`, `VecAYPX`, `VecPointwiseMult`, `VecNormalize` all report 100% GPU %F.

*2. Sparse matrix storage and transfers*

System matrices are assembled on the host and uploaded once per Newton iteration via `MatCUSPARSECopyTo`, then all SpMV operations run in cuSPARSE format. Only 3 GpuToCpu copies occur per solve (one per Newton iteration), confirming matrices remain GPU-resident throughout.

| Resolution | Matrix upload size (MB) |
|---|---|
| 513²  |    495 |
| 1025² |  1,990 |
| 2049² |  7,950 |
| 4097² | 31,800 |

*3. Vector transfers*

Vectors stay GPU-resident. The number of host↔device copies is small relative to the hundreds of thousands of vector operations performed:

| Resolution | CpuToGpu count | CpuToGpu (MB) | GpuToCpu count | GpuToCpu (MB) |
|---|---|---|---|---|
| 513²  | 37 |    109 | 24 |    42 |
| 1025² | 55 |    437 | 36 |   168 |
| 2049² | 64 |  1,750 | 42 |   672 |
| 4097² | 73 |  6,980 | 48 | 2,690 |

*4. Residual evaluation (KernelAbstractions)*

`MatFDColorApply` — which drives all residual evaluations for the finite-difference Jacobian — reports 0% GPU %F in PETSc's profiler. This is expected: the residual kernel is launched by Julia's CUDA.jl runtime (via KernelAbstractions) and its flops are invisible to PETSc's event system. GPU execution is confirmed indirectly by the GpuToCpu transfer pattern in `MatFDColorApply`: PETSc hands off perturbed vectors, the KA kernel evaluates the residual on the GPU, and the result is returned. On the GH200's unified memory architecture these transfers are intra-device and incur minimal latency.