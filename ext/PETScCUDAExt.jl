module PETScCUDAExt

using PETSc
using PETSc: LibPETSc, AbstractPetscVec, AbstractPetscMemBackend
using PETSc.LibPETSc: PETSC_MEMTYPE_DEVICE
using CUDA

# ============================================================================
#   CUDA memory backend
# ============================================================================

struct CUDAMemBackend <: AbstractPetscMemBackend end

PETSc.memtype_backend(::Val{PETSC_MEMTYPE_DEVICE}) = CUDAMemBackend()
PETSc.array_type(::Val{PETSC_MEMTYPE_DEVICE}) = CuArray

# With the extension loaded, `memtype = :native` checkouts wrap device memory.
__init__() = push!(PETSc.device_memtypes, PETSC_MEMTYPE_DEVICE)

# CUDA.jl launches kernels on its own stream, PETSc reads on its own.
PETSc.synchronize_device(::CUDAMemBackend) = CUDA.synchronize()

# A CuArray over the device memory PETSc handed out. PETSc keeps ownership.
device_array(arr) = CUDA.unsafe_wrap(CuArray,
    reinterpret(CuPtr{eltype(arr)}, UInt(pointer(arr))), length(arr); own = false)

# ============================================================================
#   with_local_array!
# ============================================================================

PETSc.make_local_array(cpu_arr, ::CUDAMemBackend) = device_array(cpu_arr)

# The finalizer-based checkout, handed back when the CuArray is finalized.
# with_local_array! uses make_local_array and release_local_array instead.
function PETSc.wrap_local_array(
    cpu_arr, ::CUDAMemBackend, vec::AbstractPetscVec{PetscLib};
    read::Bool, write::Bool,
) where {PetscLib}
    dev_arr = device_array(cpu_arr)
    pv = PETSc.as_petsc_vec(vec)
    finalizer(dev_arr) do _
        if write && read
            LibPETSc.VecRestoreArrayAndMemType(PetscLib, pv, cpu_arr)
        elseif write
            LibPETSc.VecRestoreArrayWriteAndMemType(PetscLib, pv, cpu_arr)
        else
            LibPETSc.VecRestoreArrayReadAndMemType(PetscLib, pv, cpu_arr)
        end
        return nothing
    end
    return dev_arr
end

# ============================================================================
#   with_local_arrays!
# ============================================================================

# Both Vecs on the device: wrap both, nothing to copy.
function PETSc._local_arrays(
    petsclib, g_fx, l_x, ::Type, fx_arr, lx_arr, ::CUDAMemBackend, ::CUDAMemBackend,
)
    return device_array(fx_arr), device_array(lx_arr), fx_arr, lx_arr, nothing
end

# One Vec on the host (coarse multigrid levels, FD coloring): the kernel still
# runs on the device. A host `l_x` is copied over; a host `g_fx` gets a device
# copy (the bounce array) that the restore copies back.
function PETSc._local_arrays(
    petsclib, g_fx, l_x, ::Type, fx_arr, lx_arr,
    fx_b::Union{Nothing, CUDAMemBackend}, lx_b::Union{Nothing, CUDAMemBackend},
)
    lx = lx_b === nothing ? CuArray(lx_arr) : device_array(lx_arr)
    fx_bounce = fx_b === nothing ? CuArray(fx_arr) : nothing
    fx = fx_bounce === nothing ? device_array(fx_arr) : fx_bounce
    return fx, lx, fx_arr, lx_arr, fx_bounce
end

function PETSc._restore_local_arrays!(
    petsclib, g_fx, l_x, fx::CuArray, lx, fx_arr::Array, lx_arr::Array, fx_bounce,
)
    CUDA.synchronize()
    fx_bounce === nothing || copyto!(fx_arr, fx_bounce)
    LibPETSc.VecRestoreArrayAndMemType(petsclib, PETSc.as_petsc_vec(g_fx), fx_arr)
    LibPETSc.VecRestoreArrayReadAndMemType(petsclib, PETSc.as_petsc_vec(l_x), lx_arr)
    return nothing
end

end # module PETScCUDAExt
