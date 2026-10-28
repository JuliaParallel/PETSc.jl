# test/gpu/device_arrays.jl
# PETScCUDAExt and PETScAMDGPUExt: Vecs in device memory come out of
# with_local_array! and with_local_arrays! as CuArrays or ROCArrays, DMStag
# field views wrap them, and what a kernel writes reaches PETSc.
#
# Not part of runtests.jl: it needs a GPU and a PETSc built with CUDA or HIP.
# Run it in an environment with PETSc.jl (dev'd), MPI.jl, KernelAbstractions.jl
# and CUDA.jl or AMDGPU.jl, whose LocalPreferences.toml points set_library! at
# that PETSc:
#
#   julia --project=ENV test/gpu/device_arrays.jl cuda    # or hip

using Test
using PETSc, MPI, KernelAbstractions
using PETSc: LibPETSc

const gpu = get(ARGS, 1, "")
# The device array type and PETSc's names for a sequential and a DM vector type
if gpu == "cuda"
    using CUDA
    const DeviceArray, seq_type, dm_vec_type = CuArray, :seqcuda, "cuda"
elseif gpu == "hip"
    using AMDGPU
    const DeviceArray, seq_type, dm_vec_type = ROCArray, :seqhip, "hip"
else
    error("pass the GPU backend as the argument: cuda or hip, got $gpu")
end

MPI.Initialized() || MPI.Init()

# Rp = ∂Vx/∂x + ∂Vy/∂y on the owned elements (unit spacing), offset by `lo` to
# the views' ghost-inclusive indices
@kernel function divergence!(Rp, Vx, Vy, lo)
    I = @index(Global, Cartesian) + lo
    Rp[I] = Vx[I + CartesianIndex(1, 0)] - Vx[I] + Vy[I + CartesianIndex(0, 1)] - Vy[I]
end

petsclib = PETSc.petsclibs[1]
PETSc.initialize(petsclib)
const T = petsclib.PetscScalar
const n = 100

# A sequential Vec of the given type (`seq_type` on the device, :seq on the host).
function vec_of_type(type, values)
    v = LibPETSc.VecCreate(petsclib, MPI.COMM_SELF)
    LibPETSc.VecSetSizes(petsclib, v, petsclib.PetscInt(n), petsclib.PetscInt(n))
    PETSc.set_type!(v, type)
    PETSc.with_local_array!(a -> copyto!(a, values), v)
    return v
end

@testset "device arrays ($gpu)" begin
    @testset "with_local_array! on a device Vec" begin
        x = rand(T, n)
        v = vec_of_type(seq_type, x)
        @test PETSc.memtype(v) === DeviceArray
        PETSc.with_local_array!(v) do a
            @test a isa DeviceArray{T}
            a .= 2 .* a
        end
        @test v[1:n] ≈ 2x
        @test PETSc.with_local_array!(sum, DeviceArray, v) ≈ sum(2x)
        @test_throws "expected array of type Vector" PETSc.with_local_array!(sum, Vector, v)
        # plain VecGetArray: PETSc copies to the host
        a = PETSc.unsafe_local_array(v; write = false)
        @test a isa Vector{T}
        @test a ≈ 2x
        finalize(a)
        PETSc.destroy!(v)
    end

    @testset "with_local_arrays!, $fx_type × $lx_type" for
            (fx_type, lx_type) in ((seq_type, seq_type), (:seq, seq_type), (seq_type, :seq))
        x = rand(T, n)
        g = vec_of_type(fx_type, ones(T, n))
        l = vec_of_type(lx_type, x)
        PETSc.with_local_arrays!(g, l) do fx, lx
            @test fx isa DeviceArray{T}
            @test lx isa DeviceArray{T}
            fx .+= 3 .* lx
        end
        @test g[1:n] ≈ 1 .+ 3x
        @test l[1:n] == x
        PETSc.destroy!(g)
        PETSc.destroy!(l)
    end

    # with_field_views! on a DMStag with device Vecs, against the same steps on the host
    @testset "DMStag field views on the device" begin
        none = ntuple(_ -> LibPETSc.DM_BOUNDARY_NONE, 2)
        function stag_result(vec_type)
            dm = PETSc.DMStag(petsclib, MPI.COMM_SELF, none, (8, 6), (0, 1, 1); dm_vec_type = vec_type)
            x = PETSc.local_vec(dm)
            r = PETSc.local_vec(dm)
            m = length(x)
            PETSc.with_local_array!(a -> copyto!(a, T.(1:m) .^ 2), x)
            PETSc.with_local_array!(a -> fill!(a, 0), r)
            flow = (PETSc.face_location(dm, 1) => 0, PETSc.face_location(dm, 2) => 0,
                    PETSc.element_location(dm) => 0)
            c = PETSc.corners(dm)
            lo = c.lower - CartesianIndex(1, 1)
            # OffsetArray broadcasts element by element whatever its parent, so
            # whole-view broadcasts go through `parent`; kernels take the views as they are
            PETSc.with_field_views!(dm, x, r; fields = flow, write = (false, true)) do (Vx, Vy, P), (Rx, Ry, Rp)
                vec_type == dm_vec_type && @test parent(parent(Rp)) isa DeviceArray{T}
                backend = get_backend(Rp)
                parent(Rx) .= 3
                divergence!(backend, 64)(Rp, Vx, Vy, lo; ndrange = Tuple(c.upper - c.lower) .+ 1)
                KernelAbstractions.synchronize(backend)
            end
            # the whole array, indexed [I..., slot]
            PETSc.with_field_views!(dm, r) do R
                parent(R) .*= 2
            end
            result = r[1:m]
            PETSc.destroy!(x)
            PETSc.destroy!(r)
            PETSc.destroy!(dm)
            return result
        end
        host = stag_result("standard")
        @test stag_result(dm_vec_type) == host
        @test count(!iszero, host) > 0
    end

    @testset "host Vecs stay on the host" begin
        g = vec_of_type(:seq, ones(T, n))
        l = vec_of_type(:seq, ones(T, n))
        PETSc.with_local_arrays!(g, l) do fx, lx
            @test fx isa Vector{T}
            @test lx isa Vector{T}
        end
        PETSc.destroy!(g)
        PETSc.destroy!(l)
    end
end

PETSc.finalize(petsclib)
