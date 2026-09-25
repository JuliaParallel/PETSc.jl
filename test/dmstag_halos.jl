# test/dmstag_halos.jl
# local_to_local! (in place and into another vector), the product coordinate
# do-block, and matrices preallocated without their nonzero pattern.

using Test
using PETSc
using MPI

MPI.Initialized() || MPI.Init()

@testset "DMStag halos, coordinates and preallocation" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    comm = MPI.COMM_SELF
    LibPETSc = PETSc.LibPETSc
    none(n) = ntuple(_ -> PETSc.DM_BOUNDARY_NONE, n)

    @testset "local_to_local!" begin
        # 1D, periodic, one element dof: the local vector is [x4, x1, x2, x3, x4, x1]
        dm = PETSc.DMStag(petsclib, comm, (PETSc.DM_BOUNDARY_PERIODIC,), (4,), (0, 1), 1)
        g = PETSc.global_vec(dm)
        g[1:4] = [1.0, 2.0, 3.0, 4.0]
        PETSc.assemble!(g)
        l = PETSc.local_vec(dm)
        PETSc.global_to_local!(l, dm, g)
        @test l[:] == [4.0, 1.0, 2.0, 3.0, 4.0, 1.0]

        # change the first owned element; its periodic ghost follows in place
        l[2] = 10.0
        PETSc.assemble!(l)
        @test PETSc.local_to_local!(l, dm) === l
        @test l[6] == 10.0

        dst = PETSc.local_vec(dm)
        @test PETSc.local_to_local!(dst, dm, l, PETSc.INSERT_VALUES) === dst
        @test dst[1] == 4.0
        @test dst[6] == 10.0
        foreach(PETSc.destroy!, (dst, l, g, dm))
    end

    @testset "with_product_coordinates" begin
        dm = PETSc.DMStag(petsclib, comm, none(2), (3, 2), (1, 1, 1), 1)
        PETSc.set_uniform_coordinates!(dm, (0.0, 0.0), (3.0, 2.0))

        # slot 1 is the lower face, slot 2 the centre
        got = PETSc.with_product_coordinates(dm) do x, y
            (x[1, 1], x[1, 2], x[3, 2], y[2, 1], y[2, 2])
        end
        @test got == (0.0, 0.5, 2.5, 1.0, 1.5)

        # the arrays go back when the block throws, so the coordinates can be set again
        @test_throws ErrorException PETSc.with_product_coordinates((x, y) -> error("inside"), dm)
        PETSc.set_uniform_coordinates!(dm, (0.0, 0.0), (6.0, 4.0))
        @test PETSc.with_product_coordinates((x, y) -> x[1, 2], dm) == 1.0

        dm1 = PETSc.DMStag(petsclib, comm, none(1), (4,), (1, 1), 1)
        PETSc.set_uniform_coordinates!(dm1, (0.0,), (1.0,))
        @test PETSc.with_product_coordinates(x -> x[4, 2], dm1) == 0.875
        dm3 = PETSc.DMStag(petsclib, comm, none(3), (2, 2, 2), (0, 0, 0, 1), 1)
        PETSc.set_uniform_coordinates!(dm3, (0.0, 0.0, 0.0), (2.0, 2.0, 2.0))
        @test PETSc.with_product_coordinates((x, y, z) -> z[2, 1], dm3) == 1.0
        foreach(PETSc.destroy!, (dm3, dm1, dm))
    end

    @testset "set_matrix_preallocate_only!" begin
        dm = PETSc.DMStag(petsclib, comm, none(2), (3, 2), (0, 1, 1), 1)
        nonzeros(A) = LibPETSc.MatGetInfo(petsclib, A, LibPETSc.MAT_LOCAL).nz_used
        A = PETSc.PetscMat(dm)
        @test nonzeros(A) > 0                   # the stencil pattern, stored as zeros
        @test PETSc.set_matrix_preallocate_only!(dm, true) === dm
        B = PETSc.PetscMat(dm)
        @test nonzeros(B) == 0
        foreach(PETSc.destroy!, (B, A, dm))
    end

    PETSc.finalize(petsclib)
end
