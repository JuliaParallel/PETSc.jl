module PETScCUDAExt

using PETSc
using PETSc: LibPETSc, AbstractPetscVec
using PETSc.LibPETSc: PETSC_MEMTYPE_DEVICE
using CUDA

# ── CUDA memory backend ───────────────────────────────────────────────────────

struct CUDABackend <: PETSc.AbstractPETScMemBackend end

PETSc._memtype_backend(::Val{PETSC_MEMTYPE_DEVICE}) = CUDABackend()

# ── _wrap_localarray: device branch ──────────────────────────────────────────
#
# Called by `_unsafe_localarray_device` (defined in base vec.jl) when the Vec
# is device-resident.  Wraps the device pointer as a CuArray (zero-copy) and
# attaches a finalizer that calls the matching VecRestoreArray*AndMemType.
# `cpu_arr` (the raw Julia Vector wrapping the device pointer) is captured in
# the finalizer closure to keep it alive until the restore is done.

function PETSc._wrap_localarray(
    cpu_arr, ::CUDABackend, vec::AbstractPetscVec{PetscLib};
    read::Bool, write::Bool,
) where {PetscLib}
    T   = eltype(cpu_arr)
    n   = length(cpu_arr)
    ptr = reinterpret(CuPtr{T}, UInt(pointer(cpu_arr)))
    dev_arr = CUDA.unsafe_wrap(CuArray, ptr, n; own = false)
    finalizer(dev_arr) do _
        if write && read
            LibPETSc.VecRestoreArrayAndMemType(PetscLib, vec, cpu_arr)
        elseif write
            LibPETSc.VecRestoreArrayWriteAndMemType(PetscLib, vec, cpu_arr)
        else
            LibPETSc.VecRestoreArrayReadAndMemType(PetscLib, vec, cpu_arr)
        end
        return nothing
    end
    return dev_arr
end

# ── _get_petsc_arrays_impl: CUDA cases ───────────────────────────────────────
#
# Two methods cover all GPU sub-cases:
#
#   CUDABackend × CUDABackend → both Vecs on device: zero-copy wrap, no bounce
#   any other mix             → at least one Vec is host-resident:
#                               lx is wrapped zero-copy if on device, or
#                               copied H2D if on host;
#                               fx always gets a fresh scratch CuArray (bounce)
#                               so the kernel writes there and restore copies D2H.
#
# HostBackend × HostBackend is handled entirely in base (vec.jl) and never
# reaches these methods.

# Both Vecs on the device: zero-copy wrap, no scratch needed.
function PETSc._get_petsc_arrays_impl(
    petsclib, g_fx, l_x, ::Type{T}, fx_arr, lx_arr,
    ::CUDABackend, ::CUDABackend,
) where {T}
    fx = CUDA.unsafe_wrap(CuArray,
        reinterpret(CuPtr{T}, UInt(pointer(fx_arr))), length(fx_arr))
    lx = CUDA.unsafe_wrap(CuArray,
        reinterpret(CuPtr{T}, UInt(pointer(lx_arr))), length(lx_arr))
    return fx, lx, fx_arr, lx_arr, nothing
end

# At least one Vec is host-resident (e.g. MG coarser levels, FD-coloring path).
# Catch-all: less specific than (CUDABackend, CUDABackend), so Julia prefers
# the method above when both are on the device.
function PETSc._get_petsc_arrays_impl(
    petsclib, g_fx, l_x, ::Type{T}, fx_arr, lx_arr,
    fx_b::PETSc.AbstractPETScMemBackend, lx_b::PETSc.AbstractPETScMemBackend,
) where {T}
    lx_gpu = if lx_b isa CUDABackend
        CUDA.unsafe_wrap(CuArray,
            reinterpret(CuPtr{T}, UInt(pointer(lx_arr))), length(lx_arr))
    else
        tmp = CuArray{T}(undef, length(lx_arr))
        copyto!(tmp, lx_arr)    # H2D: send ghost input to GPU
        tmp
    end
    fx_gpu = CuArray{T}(undef, length(fx_arr))  # scratch buffer for residual
    return fx_gpu, lx_gpu, fx_arr, lx_arr, fx_gpu
end

# ── _restore_petsc_arrays_impl: CUDA ─────────────────────────────────────────
#
# When fx is a CuArray (returned by the GPU _get_petsc_arrays_impl above):
#   - if fx_bounce !== nothing, sync the device and copy the scratch D2H
#   - call VecRestoreArray*AndMemType on both raw PETSc arrays

function PETSc._restore_petsc_arrays_impl(
    petsclib, g_fx, l_x, fx::CuArray, lx, fx_arr, lx_arr, fx_bounce,
)
    if fx_bounce !== nothing
        CUDA.synchronize()
        copyto!(fx_arr, fx_bounce)  # D2H: copy residual back to host PETSc array
    end
    LibPETSc.VecRestoreArrayAndMemType(petsclib, g_fx, fx_arr)
    LibPETSc.VecRestoreArrayReadAndMemType(petsclib, l_x, lx_arr)
end

end # module PETScCUDAExt
