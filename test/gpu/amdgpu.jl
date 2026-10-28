# test/gpu/amdgpu.jl
# PETScAMDGPUExt: Vecs in HIP memory come out of with_local_array! and
# with_local_arrays! as ROCArrays, DMStag field views wrap them, and what a
# kernel writes reaches PETSc.
#
# Not part of runtests.jl: it needs an AMD GPU and a PETSc built with HIP. Run it
# in an environment with PETSc.jl (dev'd), AMDGPU.jl and MPI.jl whose
# LocalPreferences.toml points set_library! at that PETSc:
#
#   julia --project=ENV test/gpu/amdgpu.jl

using Test
using PETSc, MPI, AMDGPU, KernelAbstractions
using PETSc: LibPETSc

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

# A sequential Vec of the given type (:seqhip on the device, :seq on the host).
function vec_of_type(type, values)
    v = LibPETSc.VecCreate(petsclib, MPI.COMM_SELF)
    LibPETSc.VecSetSizes(petsclib, v, petsclib.PetscInt(n), petsclib.PetscInt(n))
    PETSc.set_type!(v, type)
    PETSc.with_local_array!(a -> copyto!(a, values), v)
    return v
end

@testset "PETScAMDGPUExt" begin
    @testset "with_local_array! on a HIP Vec" begin
        x = rand(T, n)
        v = vec_of_type(:seqhip, x)
        @test PETSc.memtype(v) === ROCArray
        PETSc.with_local_array!(v) do a
            @test a isa ROCArray{T}
            a .= 2 .* a
        end
        @test v[1:n] ≈ 2x
        @test PETSc.with_local_array!(sum, ROCArray, v) ≈ sum(2x)
        @test_throws "expected array of type Vector" PETSc.with_local_array!(sum, Vector, v)
        # plain VecGetArray: PETSc copies to the host
        a = PETSc.unsafe_local_array(v; write = false)
        @test a isa Vector{T}
        @test a ≈ 2x
        finalize(a)
        PETSc.destroy!(v)
    end

    @testset "with_local_arrays!, $fx_type × $lx_type" for
            (fx_type, lx_type) in ((:seqhip, :seqhip), (:seq, :seqhip), (:seqhip, :seq))
        x = rand(T, n)
        g = vec_of_type(fx_type, ones(T, n))
        l = vec_of_type(lx_type, x)
        PETSc.with_local_arrays!(g, l) do fx, lx
            @test fx isa ROCArray{T}
            @test lx isa ROCArray{T}
            fx .+= 3 .* lx
        end
        @test g[1:n] ≈ 1 .+ 3x
        @test l[1:n] == x
        PETSc.destroy!(g)
        PETSc.destroy!(l)
    end

    # with_field_views! on a DMStag with HIP Vecs, against the same steps on the host
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
                vec_type == "hip" && @test parent(parent(Rp)) isa ROCArray{T}
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
        @test stag_result("hip") == host
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
