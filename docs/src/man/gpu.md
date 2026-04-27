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
