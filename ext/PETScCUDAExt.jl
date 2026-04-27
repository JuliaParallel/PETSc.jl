module PETScCUDAExt

using PETSc
using PETSc: LibPETSc, AbstractPetscVec
using PETSc.LibPETSc: PetscMemType, PETSC_MEMTYPE_HOST
using CUDA

# ── Internal: get one device-or-host array from a single Vec ─────────────────
#
# Uses VecGetArray{,Read,Write}AndMemType so PETSc tells us where the data is:
#   PETSC_MEMTYPE_HOST   → return a plain Julia Vector (no copy)
#   anything else (CUDA) → wrap the device pointer as a CuArray (no copy)
#
# The returned array has a finalizer that calls the matching VecRestoreArray*
# so the caller just needs to finalize it when done, exactly like unsafe_localarray.
#
function _unsafe_localarray_device(
    vec::AbstractPetscVec{PetscLib};
    read::Bool = true,
    write::Bool = true,
) where {PetscLib}

    if write && read
        cpu_arr, mtype = LibPETSc.VecGetArrayAndMemType(PetscLib, vec)
    elseif write
        cpu_arr, mtype = LibPETSc.VecGetArrayWriteAndMemType(PetscLib, vec)
    else
        cpu_arr, mtype = LibPETSc.VecGetArrayReadAndMemType(PetscLib, vec)
    end

    if mtype === PETSC_MEMTYPE_HOST
        # Data is on the host — attach a restore finalizer and return as-is.
        finalizer(cpu_arr) do a
            if write && read
                LibPETSc.VecRestoreArrayAndMemType(PetscLib, vec, a)
            elseif write
                LibPETSc.VecRestoreArrayWriteAndMemType(PetscLib, vec, a)
            else
                LibPETSc.VecRestoreArrayReadAndMemType(PetscLib, vec, a)
            end
            return nothing
        end
        return cpu_arr
    else
        # Data is on the GPU — wrap the device pointer as a CuArray.
        # cpu_arr holds the raw device pointer in a Julia Vector shell; we must
        # keep it alive (captured in the finalizer) so the pointer stays valid.
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
end

# ── Public hook: register CUDA implementation for withlocalarray_device! ──────
#
# We cannot override PETSc.withlocalarray_device! with the same signature
# during precompilation (Julia restriction).  Instead we register a closure
# in __init__ that the base method dispatches to when the hook is non-nothing.
#
function _cuda_withlocalarray_device_impl!(
    f!,
    vecs::NTuple{N};
    read::Union{Bool, NTuple{N, Bool}}  = true,
    write::Union{Bool, NTuple{N, Bool}} = true,
) where {N}
    read  isa NTuple{N, Bool} || (read  = ntuple(_ -> read,  N))
    write isa NTuple{N, Bool} || (write = ntuple(_ -> write, N))

    arrays = map(vecs, read, write) do v, r, w
        _unsafe_localarray_device(v; read = r, write = w)
    end

    val = f!(arrays...)

    map(Base.finalize, arrays)

    return val
end

# ── GPU-aware PETSc array helpers ─────────────────────────────────────────────
#
# CUDA implementations of PETSc.get_petsc_arrays / PETSc.restore_petsc_arrays.
# Registered as hooks in __init__ so the base-module functions dispatch here
# whenever CUDA.jl is loaded.
#
# Three sub-cases handled by get:
#   Both Vecs on GPU  → zero-copy CuArray wraps, fx_bounce = nothing
#   lx on GPU only    → wrap lx zero-copy; allocate GPU scratch for fx (bounce)
#   Both Vecs on CPU  → copy lx H2D; allocate GPU scratch for fx (bounce)
#
# restore then D2H-copies the bounce buffer (if any) before calling
# VecRestoreArray*AndMemType on both Vecs.

function _cuda_get_petsc_arrays_impl(petsclib, g_fx, l_x)
    T      = petsclib.PetscScalar
    fx_arr, fx_mtype = LibPETSc.VecGetArrayAndMemType(petsclib, g_fx)
    lx_arr, lx_mtype = LibPETSc.VecGetArrayReadAndMemType(petsclib, l_x)

    if fx_mtype == LibPETSc.PETSC_MEMTYPE_DEVICE &&
       lx_mtype == LibPETSc.PETSC_MEMTYPE_DEVICE
        # Both on GPU: zero-copy wrap, no bounce needed.
        fx = CUDA.unsafe_wrap(CuArray,
            reinterpret(CuPtr{T}, UInt(pointer(fx_arr))), length(fx_arr))
        lx = CUDA.unsafe_wrap(CuArray,
            reinterpret(CuPtr{T}, UInt(pointer(lx_arr))), length(lx_arr))
        return fx, lx, fx_arr, lx_arr, nothing
    else
        # At least one Vec is host-resident (e.g. freshly created coarser MG
        # level, or FD-coloring CPU path).  Wrap or copy lx to GPU as needed,
        # and allocate a GPU scratch buffer for fx so the kernel can write there;
        # restore_petsc_arrays copies it back D2H after the kernel.
        lx_gpu = if lx_mtype == LibPETSc.PETSC_MEMTYPE_DEVICE
            CUDA.unsafe_wrap(CuArray,
                reinterpret(CuPtr{T}, UInt(pointer(lx_arr))), length(lx_arr))
        else
            tmp = CuArray{T}(undef, length(lx_arr))
            copyto!(tmp, lx_arr)        # H2D: send ghost input to GPU
            tmp
        end
        fx_gpu = CuArray{T}(undef, length(fx_arr))
        return fx_gpu, lx_gpu, fx_arr, lx_arr, fx_gpu
    end
end

function _cuda_restore_petsc_arrays_impl(
    petsclib, g_fx, l_x, fx, lx, fx_arr, lx_arr, fx_bounce,
)
    if fx_bounce !== nothing
        # D2H: copy GPU residual result back to the host PETSc array.
        CUDA.synchronize()
        copyto!(fx_arr, fx_bounce)
    end
    LibPETSc.VecRestoreArrayAndMemType(petsclib, g_fx, fx_arr)
    LibPETSc.VecRestoreArrayReadAndMemType(petsclib, l_x, lx_arr)
end

function __init__()
    PETSc._withlocalarray_device_hook[] = _cuda_withlocalarray_device_impl!
    PETSc._get_petsc_arrays_hook[]      = _cuda_get_petsc_arrays_impl
    PETSc._restore_petsc_arrays_hook[]  = _cuda_restore_petsc_arrays_impl
end

end # module
