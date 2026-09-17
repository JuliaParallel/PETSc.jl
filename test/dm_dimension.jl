using Test
using PETSc, MPI
MPI.Initialized() || MPI.Init()

# Dimension-correct returns (docs/src/man/naming.md §12).
#
# The point of building these with `ntuple(..., Val(N))` from the DM's type
# parameter rather than by splatting a `Vector` is that the result infers
# concretely and costs nothing: v0.4's `getcorners_dmda` hid the length from the
# compiler, so `lower` and `upper` inferred as `Any` and every call allocated.
# These tests are what stops that coming back.

@testset "dimension-correct DM returns" begin
    comm = MPI.COMM_SELF
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)

    # A serial DMDA and DMStag of each dimension.
    dmda = ntuple(3) do N
        PETSc.DMDA(
            petsclib,
            comm,
            ntuple(_ -> PETSc.DM_BOUNDARY_NONE, N),
            ntuple(i -> 8 + i, N),
            2,
            1,
            PETSc.DMDA_STENCIL_STAR,
        )
    end
    dmstag = ntuple(3) do N
        PETSc.DMStag(
            petsclib,
            comm,
            ntuple(_ -> PETSc.DM_BOUNDARY_NONE, N),
            ntuple(i -> 8 + i, N),
            ntuple(_ -> 1, N + 1),
            1,
            :box,       # §3.1: the Symbol spelling of the stencil type
        )
    end

    @testset "shapes and keys" begin
        for N in 1:3
            da, stag = dmda[N], dmstag[N]

            for dm in (da, stag)
                c = PETSc.corners(dm)
                gc = PETSc.ghost_corners(dm)
                @test c.lower isa CartesianIndex{N}
                @test c.upper isa CartesianIndex{N}
                @test c.size isa NTuple{N, Int}
                @test gc.lower isa CartesianIndex{N}
                @test gc.size isa NTuple{N, Int}
                @test size(dm) isa NTuple{N, Int}
                @test size(dm) == ntuple(i -> 8 + i, N)
            end

            # `nextra` is a DMStag field, and only on `corners`:
            # `DMStagGetGhostCorners` does not report it.
            @test PETSc.corners(stag).nextra isa NTuple{N, Int}
            @test !hasproperty(PETSc.ghost_corners(stag), :nextra)
            @test !hasproperty(PETSc.corners(da), :nextra)

            # `center`/`vertex` are keyed by axis, up to the DM's dimension.
            axes_expected = (:x, :y, :z)[1:N]
            for idx in (PETSc.local_indices(stag), PETSc.global_indices(stag))
                @test keys(idx.center) === axes_expected
                @test keys(idx.vertex) === axes_expected
                @test all(r -> r isa UnitRange{Int}, values(idx.center))
            end

            # `info` reports the new field set, all of it `N`-dimensional.
            i = PETSc.info(da)
            @test keys(i) === (
                :dim,
                :global_size,
                :procs,
                :ndofs,
                :stencil_width,
                :boundary_type,
                :stencil_type,
            )
            @test i.dim == N
            @test i.ndofs == 2
            @test i.global_size isa NTuple{N, Int}
            @test i.procs isa NTuple{N, Int}
            @test length(i.boundary_type) == N
        end
    end

    @testset "inference" begin
        # The 2D DMDA is the case §12 calls out by name.
        da2, stag2 = dmda[2], dmstag[2]
        @test @inferred(PETSc.corners(da2)).lower isa CartesianIndex{2}
        @test @inferred(PETSc.ghost_corners(da2)).size isa NTuple{2, Int}
        @test @inferred(PETSc.corners(stag2)).nextra isa NTuple{2, Int}
        @test @inferred(PETSc.info(da2)).procs isa NTuple{2, Int}
        @test @inferred(size(da2)) isa NTuple{2, Int}
        @test @inferred(size(stag2)) isa NTuple{2, Int}

        # `center` infers as a concrete NamedTuple, with no `z` in 2D.
        @test @inferred(PETSc.local_indices(stag2)).center isa
              @NamedTuple{x::UnitRange{Int}, y::UnitRange{Int}}
        @test @inferred(PETSc.global_indices(stag2)).center isa
              @NamedTuple{x::UnitRange{Int}, y::UnitRange{Int}}
    end

    @testset "allocations" begin
        # A handful of words for the `Ref`s the generated wrappers use; v0.4
        # allocated a `Vector` per call on top of that.
        budget = 128
        for N in 2:3
            da, stag = dmda[N], dmstag[N]
            for f in (PETSc.corners, PETSc.ghost_corners, PETSc.info, size)
                f(da)   # warm up
                @test @allocated(f(da)) <= budget
            end
            for f in (
                PETSc.corners,
                PETSc.ghost_corners,
                PETSc.local_indices,
                PETSc.global_indices,
                size,
            )
                f(stag)
                @test @allocated(f(stag)) <= budget
            end
        end
    end

    for dm in (dmda..., dmstag...)
        PETSc.destroy!(dm)
    end
    PETSc.finalize(petsclib)
end
