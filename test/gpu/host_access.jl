# test/gpu/host_access.jl
# Host code on device Vecs with no GPU package loaded: with_local_array!,
# with_local_arrays! and with_field_views! give host arrays, PETSc copies the
# data between host and device, and the results match a host Vec.
#
# Not part of runtests.jl: it needs a PETSc built with CUDA or HIP. Run it in an
# environment with PETSc.jl (dev'd) and MPI.jl, whose LocalPreferences.toml points
# set_library! at that PETSc, and without loading CUDA.jl or AMDGPU.jl:
#
#   julia --project=ENV test/gpu/host_access.jl cuda    # or hip

using Test
using PETSc, MPI
using PETSc: LibPETSc
using LinearAlgebra: norm

const gpu = get(ARGS, 1, "")
# PETSc's names for a sequential and a DM vector type on the device
if gpu == "cuda"
    const seq_type, dm_vec_type = :seqcuda, "cuda"
elseif gpu == "hip"
    const seq_type, dm_vec_type = :seqhip, "hip"
else
    error("pass the GPU backend as the argument: cuda or hip, got $gpu")
end

MPI.Initialized() || MPI.Init()

petsclib = PETSc.petsclibs[1]
PETSc.initialize(petsclib)
const T = petsclib.PetscScalar
const n = 100

function vec_of_type(type, values)
    v = LibPETSc.VecCreate(petsclib, MPI.COMM_SELF)
    LibPETSc.VecSetSizes(petsclib, v, petsclib.PetscInt(n), petsclib.PetscInt(n))
    PETSc.set_type!(v, type)
    PETSc.with_local_array!(a -> copyto!(a, values), v)
    return v
end

# Rp = ∂Vx/∂x + ∂Vy/∂y on the owned elements (unit spacing), as a host loop
function divergence!(Rp, Vx, Vy, c)
    for I in c.lower:c.upper
        Rp[I] = Vx[I + CartesianIndex(1, 0)] - Vx[I] + Vy[I + CartesianIndex(0, 1)] - Vy[I]
    end
end

# with_field_views! on a DMStag of the given vector type, filled and read on the host
function stag_result(vec_type)
    none = ntuple(_ -> LibPETSc.DM_BOUNDARY_NONE, 2)
    dm = PETSc.DMStag(petsclib, MPI.COMM_SELF, none, (8, 6), (0, 1, 1); dm_vec_type = vec_type)
    x = PETSc.local_vec(dm)
    vec_type == "standard" || @test occursin(vec_type, String(PETSc.type_name(x)))
    r = PETSc.local_vec(dm)
    m = length(x)
    PETSc.with_local_array!(a -> copyto!(a, T.(1:m) .^ 2), x)
    PETSc.with_local_array!(a -> fill!(a, 0), r)
    flow = (PETSc.face_location(dm, 1) => 0, PETSc.face_location(dm, 2) => 0,
            PETSc.element_location(dm) => 0)
    c = PETSc.corners(dm)
    PETSc.with_field_views!(dm, x, r; fields = flow, write = (false, true)) do (Vx, Vy, P), (Rx, Ry, Rp)
        @test parent(parent(Rp)) isa Array{T}
        parent(Rx) .= 3
        divergence!(Rp, Vx, Vy, c)
    end
    # PETSc works on the device copy, so the host writes must have reached it
    LibPETSc.VecScale(petsclib, r, T(2))
    result = r[1:m]
    PETSc.destroy!(x)
    PETSc.destroy!(r)
    PETSc.destroy!(dm)
    return result
end

@testset "host access to device Vecs ($gpu, no GPU package)" begin
    @test isempty(PETSc.device_memtypes)

    @testset "with_local_array!" begin
        x = rand(T, n)
        v = vec_of_type(seq_type, x)
        @test PETSc.type_name(v) === seq_type
        # a device operation between two host accesses: each sees the other's data
        PETSc.with_local_array!(v) do a
            @test a isa Vector{T}
            a .= 2 .* a
        end
        LibPETSc.VecScale(petsclib, v, T(3))
        @test PETSc.with_local_array!(copy, v; write = false) ≈ 6x
        @test PETSc.with_local_array!(sum, v; memtype = :host) ≈ sum(6x)
        @test PETSc.memtype(v) === Vector
        PETSc.with_local_array!(a -> fill!(a, 1), v; read = false)
        @test norm(v) ≈ sqrt(n)
        PETSc.destroy!(v)
    end

    @testset "with_local_arrays!, $fx_type × $lx_type" for
            (fx_type, lx_type) in ((seq_type, seq_type), (:seq, seq_type), (seq_type, :seq))
        x = rand(T, n)
        g = vec_of_type(fx_type, ones(T, n))
        l = vec_of_type(lx_type, x)
        PETSc.with_local_arrays!(g, l) do fx, lx
            @test fx isa Vector{T}
            @test lx isa Vector{T}
            fx .+= 3 .* lx
        end
        @test g[1:n] ≈ 1 .+ 3x
        @test l[1:n] == x
        PETSc.destroy!(g)
        PETSc.destroy!(l)
    end

    @testset "DMStag field views" begin
        host = stag_result("standard")
        @test stag_result(dm_vec_type) == host
        @test count(!iszero, host) > 0
    end

    @testset "broadcasting" begin
        x = rand(T, n)
        v = vec_of_type(seq_type, x)
        v .= 2 .* v .+ 1
        @test v[1:n] ≈ 2 .* x .+ 1
        PETSc.destroy!(v)
    end
end

PETSc.finalize(petsclib)
