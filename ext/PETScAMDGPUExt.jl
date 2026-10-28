module PETScAMDGPUExt

using PETSc
using PETSc: LibPETSc, AbstractPetscVec, AbstractPetscMemBackend
using PETSc.LibPETSc: PETSC_MEMTYPE_HIP
using AMDGPU

# ============================================================================
#   HIP memory backend
# ============================================================================

struct HIPMemBackend <: AbstractPetscMemBackend end

PETSc.memtype_backend(::Val{PETSC_MEMTYPE_HIP}) = HIPMemBackend()
PETSc.array_type(::Val{PETSC_MEMTYPE_HIP}) = ROCArray

# AMDGPU.jl launches kernels on its own stream, PETSc reads on its own.
PETSc.synchronize_device(::HIPMemBackend) = AMDGPU.synchronize()

# A ROCArray over the device memory PETSc handed out. PETSc keeps ownership.
device_array(arr) = unsafe_wrap(ROCArray, pointer(arr), length(arr); own = false)

# ============================================================================
#   with_local_array!
# ============================================================================

PETSc.make_local_array(cpu_arr, ::HIPMemBackend) = device_array(cpu_arr)

# ============================================================================
#   with_local_arrays!
# ============================================================================

# Both Vecs on the device: wrap both, nothing to copy.
function PETSc._local_arrays(
    petsclib, g_fx, l_x, ::Type, fx_arr, lx_arr, ::HIPMemBackend, ::HIPMemBackend,
)
    return device_array(fx_arr), device_array(lx_arr), fx_arr, lx_arr, nothing
end

# One Vec on the host (coarse multigrid levels, FD coloring): the kernel still
# runs on the device. A host `l_x` is copied over; a host `g_fx` gets a device
# copy (the bounce array) that the restore copies back.
function PETSc._local_arrays(
    petsclib, g_fx, l_x, ::Type, fx_arr, lx_arr,
    fx_b::Union{Nothing, HIPMemBackend}, lx_b::Union{Nothing, HIPMemBackend},
)
    lx = lx_b === nothing ? ROCArray(lx_arr) : device_array(lx_arr)
    fx_bounce = fx_b === nothing ? ROCArray(fx_arr) : nothing
    fx = fx_bounce === nothing ? device_array(fx_arr) : fx_bounce
    return fx, lx, fx_arr, lx_arr, fx_bounce
end

function PETSc._restore_local_arrays!(
    petsclib, g_fx, l_x, fx::ROCArray, lx, fx_arr::Array, lx_arr::Array, fx_bounce,
)
    AMDGPU.synchronize()
    fx_bounce === nothing || copyto!(fx_arr, fx_bounce)
    LibPETSc.VecRestoreArrayAndMemType(petsclib, PETSc.as_petsc_vec(g_fx), fx_arr)
    LibPETSc.VecRestoreArrayReadAndMemType(petsclib, PETSc.as_petsc_vec(l_x), lx_arr)
    return nothing
end

end # module PETScAMDGPUExt
