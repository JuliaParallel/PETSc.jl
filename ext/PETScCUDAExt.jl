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

function __init__()
    PETSc._withlocalarray_device_hook[] = _cuda_withlocalarray_device_impl!
end

end # module
