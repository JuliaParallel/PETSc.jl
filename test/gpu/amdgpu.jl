# test/gpu/amdgpu.jl
# PETScAMDGPUExt: Vecs in HIP memory come out of with_local_array! and
# with_local_arrays! as ROCArrays, and what a kernel writes reaches PETSc.
#
# Not part of runtests.jl: it needs an AMD GPU and a PETSc built with HIP. Run it
# in an environment with PETSc.jl (dev'd), AMDGPU.jl and MPI.jl whose
# LocalPreferences.toml points set_library! at that PETSc:
#
#   julia --project=ENV test/gpu/amdgpu.jl

using Test
using PETSc, MPI, AMDGPU
using PETSc: LibPETSc

MPI.Initialized() || MPI.Init()

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
