# GPU Support (CUDA and HIP)

PETSc vectors, matrices and many of its solvers can live on a GPU. PETSc.jl gives you the device memory as a Julia array, so you write the residual or Jacobian kernel in Julia, for example with [KernelAbstractions](https://github.com/JuliaGPU/KernelAbstractions.jl), and the same code runs on NVIDIA GPUs, AMD GPUs and the CPU.

GPU support needs a **locally built PETSc** with CUDA or HIP enabled: the precompiled `PETSc_jll` binaries do not include it. See [Installation](@ref) for pointing PETSc.jl to a local library.

|  | NVIDIA | AMD |
|---|---|---|
| Julia package | CUDA.jl | AMDGPU.jl |
| Extension, loaded automatically | `PETScCUDAExt` | `PETScAMDGPUExt` |
| PETSc configure option | `--with-cuda=1` | `--with-hip=1` |
| PETSc vector and matrix types | `cuda`, `aijcusparse` | `hip`, `aijhipsparse` |
| Array type you get | `CuArray` | `ROCArray` |

The examples below use CUDA. For AMD, load AMDGPU instead and expect `ROCArray` where the text says `CuArray`.

## Setup

**1. Build PETSc with GPU support.** `scripts/build_petsc.sh` builds PETSc for `set_library!` and passes everything after `--` to PETSc's configure:

```bash
scripts/build_petsc.sh ~/petsc-gpu -- --with-cuda=1 --with-cuda-dir=/usr/local/cuda --with-cuda-arch=90   # your GPU's compute capability
scripts/build_petsc.sh ~/petsc-gpu -- --with-hip=1 --with-hip-dir=$ROCM_DIR --with-hip-arch=gfx1100   # your GPU's architecture
```

Confirm that the build has the GPU backend:

```bash
grep -i -e cuda -e hip $PETSC_DIR/$PETSC_ARCH/include/petscconf.h
# should show: #define PETSC_HAVE_CUDA 1  (or PETSC_HAVE_HIP 1)
```

Then point PETSc.jl to the library with `PETSc.set_library!`. Pass the `PetscInt` your build uses, which `grep sizeof_PetscInt $PETSC_DIR/$PETSC_ARCH/include/petscconf.h` shows.

**2. Use the same GPU runtime in Julia and in PETSc.**
- For HIP, AMDGPU.jl loads the system ROCm, so build PETSc against that ROCm.
- CUDA.jl downloads its own CUDA runtime by default. To use the toolkit PETSc was built with, run `CUDA.set_runtime_version!(local_toolkit = true)` once in your environment, then restart Julia.

**3. Add the Julia packages.**

```julia
pkg> add CUDA                # NVIDIA, or: add AMDGPU
pkg> add KernelAbstractions  # for portable kernels
```

**4. Check MPI.**
- A PETSc built with CUDA stops at startup with "your MPI is not GPU-aware" when MPI.jl uses its default MPI (MPICH_jll), which cannot pass device memory.
- With one rank, turn the check off with the option `-use_gpu_aware_mpi 0`, for example `PETSC_OPTIONS="-use_gpu_aware_mpi 0"` in the environment.
- With several ranks (one per GPU), use a GPU-aware system MPI instead; see [Running on HPC Systems](hpc.md).

## Putting data on the device

PETSc decides where a vector's data lives from its type. Ask for the device types when you create the object:

- **A DM:** pass the options as keywords, `DMStag(...; dm_vec_type = "cuda", dm_mat_type = "aijcusparse")`, or give `-dm_vec_type cuda -dm_mat_type aijcusparse` on the command line. The vectors and matrices the DM creates are then device objects.
- **A single Vec:** set the type before use, `PETSc.set_type!(v, :seqcuda)` (`:cuda` on several ranks). For HIP use `:seqhip` and `:hip`.

Everything else, including the solvers, follows the type of the objects. The `PetscMemType` that PETSc reports for each array tells PETSc.jl whether to hand you a `CuArray` or a plain `Vector`.

## Reading and writing the arrays

When CUDA.jl or AMDGPU.jl is loaded, the matching extension is active and no further setup is needed.

### `with_local_array!`

For one or more independent Vecs, `with_local_array!` gives you the local array of each inside a do-block. The array is a `CuArray` for a device Vec and a `Vector` for a host Vec, so the same code works on both. The device array wraps PETSc's memory without a copy.

```julia
using PETSc, CUDA, KernelAbstractions

@kernel function scale!(out, inp)
    i = @index(Global)
    out[i] = inp[i] * 2
end

PETSc.with_local_array!(out_vec, inp_vec; read = (true, true), write = (true, false)) do out, inp
    backend = KernelAbstractions.get_backend(out)
    scale!(backend, 256)(out, inp; ndrange = length(out))
    KernelAbstractions.synchronize(backend)
end
```

`get_backend` picks the CUDA, HIP or CPU backend from the array type, so the kernel runs on a host Vec as well.

### `with_local_arrays!`

In a residual callback you have an output Vec and an input Vec that PETSc may keep on different sides, for example on coarse multigrid levels or in PETSc's finite-difference coloring. `with_local_arrays!` checks the pair out so that you always get arrays on the same side:

```julia
PETSc.with_local_arrays!(g_fx, l_x) do fx, lx
    backend = KernelAbstractions.get_backend(fx)
    residual!(backend, 256)(fx, lx; ndrange = length(fx))
    KernelAbstractions.synchronize(backend)
end
```

- Both on the device: the arrays wrap the device memory.
- One on the host: that Vec is copied to the device and, for the output, copied back when the block returns.
- No GPU extension loaded: both are host arrays.

The device is synchronized before the arrays go back to PETSc, since the GPU packages launch kernels on their own streams.

For a staggered grid, `with_field_views!` hands out device views per field; see [DMStag](dmstag.md).

## Example

[`examples/dmda/ex19.jl`](https://github.com/JuliaParallel/PETSc.jl/blob/main/examples/dmda/ex19.jl) is a 2D driven-cavity problem (velocity, vorticity, temperature). It shows:

- switching between CPU, CUDA and HIP with a single `gpu` flag,
- a Jacobian by finite-difference coloring that runs entirely on the device,
- `with_local_arrays!` in the residual, which also covers the multigrid levels where one Vec is on the host,
- KernelAbstractions kernels that run on either vendor's GPUs or on the CPU,
- multigrid with the coarser levels falling back to a CPU Jacobian.

To run it on a GPU:

1. Set `const gpu = :cuda` (or `:hip`) near the top of `ex19.jl`.
2. Link a PETSc with CUDA (or HIP) support through `PETSc.set_library!`.
3. Launch:
   ```bash
   julia --project ex19.jl -da_grid_x 257 -da_grid_y 257 \
       -pc_type mg -pc_mg_levels 4 \
       -mg_levels_ksp_type chebyshev -mg_levels_pc_type jacobi \
       -snes_monitor -ksp_monitor
   ```

## Profiling

Wrap the solve with `CUDA.@profile`:

```julia
using CUDA
CUDA.@profile PETSc.solve!(x, snes)
```

This records a trace for NSight Systems (`nsys profile`) and NVTX. For a quick look at kernel timings, `CUDA.@profile external = true` launches the process under `nsys` for you.

PETSc's own report, `-log_view`, shows how much work ran on the device: its `GPU %F` column is the fraction of flops executed on the GPU, and the `CpuToGpu` and `GpuToCpu` lines count host-device copies.

## Performance

We timed `examples/dmda/ex19.jl` on one GPU of an NVIDIA GH200 (Grace Hopper, SM90) and on 1 and 32 cores of its CPU. The GPU runs use one MPI rank. The options were:

```bash
$mpiexec -n 1 julia --project ex19.jl -snes_monitor -snes_converged_reason -pc_type mg \
    -mg_levels_ksp_type chebyshev -mg_levels_pc_type jacobi -ksp_monitor -log_view \
    -mg_levels_esteig_ksp_max_it 20 -mg_levels_ksp_chebyshev_esteig 0,0.1,0,1.1 \
    -mg_levels_ksp_max_it 3 -da_grid_x 513 -da_grid_y 513 -pc_mg_levels 6
```

We doubled the resolution in `x` and `y` and added multigrid levels so that the coarse grid has the same size in every case. The full run also includes starting Julia and computing the coloring pattern, which doesn't scale well. The tables give PETSc's own timings of the `KSPSolve` and `SNESSolve` events.

**KSPSolve**

| Resolution | GPU (s) | GPU GFlop/s | 1 core (s) | 1 core GFlop/s | 32 cores (s) | 32 cores GFlop/s |
|---|---|---|---|---|---|---|
| 513²  |  0.116 | 144.3 |  4.337 |  4.1 |  0.297 | 61.0 |
| 1025² |  0.299 | 249.5 | 19.10  |  3.9 |  1.196 | 61.7 |
| 2049² |  1.118 | 295.5 | 89.76  |  3.6 |  5.757 | 56.6 |
| 4097² |  4.540 | 312.4 | 422.2  |  3.3 | 28.30  | 49.4 |

**SNESSolve**

| Resolution | GPU (s) | GPU GFlop/s | 1 core (s) | 1 core GFlop/s | 32 cores (s) | 32 cores GFlop/s |
|---|---|---|---|---|---|---|
| 513²  |  3.645 |  6.7 |  8.041 | 3.2 |  2.118 | 12.3 |
| 1025² |  5.127 | 20.5 | 32.75  | 3.2 |  3.698 | 28.2 |
| 2049² | 11.37  | 39.7 | 144.1  | 3.1 | 10.57  | 42.3 |
| 4097² | 36.88  | 47.7 | 658.3  | 2.9 | 43.58  | 43.2 |

`KSPSolve` is very efficient on the GPU and clearly beats the CPU. `SNESSolve` carries quite some overhead on top of it, so the gap there is much smaller.

### Is it running on the GPU?

All GPU runs use a CUDA-enabled PETSc build (cuSPARSE, cuBLAS) with KernelAbstractions.jl providing the residual kernel. PETSc's `-log_view` output answers the question directly.

**Flop fraction on the GPU.** Nearly all floating-point work runs on the GPU, at every resolution:

| Event | GPU %F |
|---|---|
| KSPSolve | 100% |
| PCApply | 100% |
| PCSetUp | 98–99% |
| SNESSolve | 99–100% |

`MatMult`, `MatResidual`, `VecMAXPY`, `VecMDot`, `VecAXPBYCZ`, `VecAYPX`, `VecPointwiseMult` and `VecNormalize` all report 100%.

**Matrix storage and transfers.** The system matrices are assembled on the host and uploaded once per Newton iteration with `MatCUSPARSECopyTo`. After that, all products run in cuSPARSE format. There are only 3 `GpuToCpu` copies per solve, one per Newton iteration, so the matrices stay on the GPU.

| Resolution | Matrix upload (MB) |
|---|---|
| 513²  |    495 |
| 1025² |  1,990 |
| 2049² |  7,950 |
| 4097² | 31,800 |

**Vector transfers.** Vectors stay on the GPU. The number of host-device copies is small next to the hundreds of thousands of vector operations:

| Resolution | CpuToGpu count | CpuToGpu (MB) | GpuToCpu count | GpuToCpu (MB) |
|---|---|---|---|---|
| 513²  | 37 |    109 | 24 |    42 |
| 1025² | 55 |    437 | 36 |   168 |
| 2049² | 64 |  1,750 | 42 |   672 |
| 4097² | 73 |  6,980 | 48 | 2,690 |

**The residual kernel.** `MatFDColorApply` reports 0% `GPU %F`. That is expected: CUDA.jl launches the residual kernel through KernelAbstractions, and its flops don't show up in PETSc's event system. The `GpuToCpu` transfers in `MatFDColorApply` show it runs on the GPU: PETSc hands over the perturbed vectors, the kernel evaluates the residual on the device, and the result goes back.
