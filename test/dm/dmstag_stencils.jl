# test/dm/dmstag_stencils.jl
# DMStag locations keyed by axis, stencils from 1-based indices, and the
# operations that take stencils: set_values!, zero_rows_local! and IS.

using Test
using PETSc
using MPI

MPI.Initialized() || MPI.Init()

@testset "DMStag locations and stencils" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    comm = MPI.COMM_SELF
    LibPETSc = PETSc.LibPETSc
    PetscInt = petsclib.PetscInt
    none(n) = ntuple(_ -> PETSc.DM_BOUNDARY_NONE, n)
    dm1 = PETSc.DMStag(petsclib, comm, none(1), (4,), (1, 1), 1)
    dm2 = PETSc.DMStag(petsclib, comm, none(2), (3, 2), (1, 1, 1), 1)
    dm3 = PETSc.DMStag(petsclib, comm, none(3), (2, 2, 2), (1, 1, 1, 1), 1)

    @testset "locations" begin
        @test PETSc.vertex_location(dm1) == LibPETSc.DMSTAG_LEFT
        @test PETSc.vertex_location(dm2) == LibPETSc.DMSTAG_DOWN_LEFT
        @test PETSc.vertex_location(dm3) == LibPETSc.DMSTAG_BACK_DOWN_LEFT
        @test [PETSc.face_location(dm3, a) for a in 1:3] ==
              [LibPETSc.DMSTAG_LEFT, LibPETSc.DMSTAG_DOWN, LibPETSc.DMSTAG_BACK]
        @test PETSc.face_location(dm2, ndims(dm2)) == LibPETSc.DMSTAG_DOWN
        @test PETSc.edge_location(dm3, 1, 2) == LibPETSc.DMSTAG_DOWN_LEFT
        @test PETSc.edge_location(dm3, 3, 1) == LibPETSc.DMSTAG_BACK_LEFT
        @test PETSc.edge_location(dm3, 2, 3) == LibPETSc.DMSTAG_BACK_DOWN
        @test PETSc.edge_location(dm2, 2, 1) == PETSc.vertex_location(dm2)
        @test PETSc.element_location(dm1) == LibPETSc.DMSTAG_ELEMENT

        @test_throws ArgumentError PETSc.face_location(dm2, 3)
        @test_throws ArgumentError PETSc.face_location(dm2, 0)
        @test_throws ArgumentError PETSc.edge_location(dm1, 1, 1)
        @test_throws ArgumentError PETSc.edge_location(dm3, 2, 2)
        @test_throws ArgumentError PETSc.edge_location(dm2, 1, 3)

        # on_lower_side against the names: LEFT, DOWN and BACK are the lower sides
        lower = ("LEFT", "DOWN", "BACK")
        for loc in instances(LibPETSc.DMStagStencilLocation)
            loc == LibPETSc.DMSTAG_NULL_LOCATION && continue
            name = string(loc)
            for axis in 1:3
                @test PETSc.on_lower_side(loc, axis) == occursin(lower[axis], name)
            end
        end
        @test PETSc.on_lower_side(PETSc.face_location(dm2, 2), 2)
        @test !any(a -> PETSc.on_lower_side(PETSc.element_location(dm3), a), 1:3)
        @test all(a -> PETSc.on_lower_side(PETSc.vertex_location(dm3), a), 1:3)
        @test_throws ArgumentError PETSc.on_lower_side(LibPETSc.DMSTAG_LEFT, 4)
        @test_throws ArgumentError PETSc.on_lower_side(LibPETSc.DMSTAG_NULL_LOCATION, 1)
    end

    @testset "stencil" begin
        # the struct matches PETSc's layout for the loaded integer width
        @test sizeof(LibPETSc.DMStagStencil) == (PetscInt === Int64 ? 40 : 20)
        @test sizeof(LibPETSc.MatStencil) == 4 * sizeof(PetscInt)

        s = PETSc.stencil(dm2, PETSc.face_location(dm2, 1), CartesianIndex(2, 1); dof = 1)
        @test (s.loc, s.i, s.j, s.k, s.c) == (LibPETSc.DMSTAG_LEFT, 1, 0, 0, 1)
        @test PETSc.stencil(dm3, PETSc.element_location(dm3), (1, 2, 2)) ==
              LibPETSc.DMStagStencil(LibPETSc.DMSTAG_ELEMENT, 0, 1, 1, 0)
        # ghost elements are not checked
        @test PETSc.stencil(dm2, PETSc.element_location(dm2), (0, 3)).i == -1

        build(dm) = PETSc.stencil(dm, PETSc.element_location(dm), (2, 1); dof = 0)
        @test @inferred(build(dm2)) isa LibPETSc.DMStagStencil
        @test @allocated(build(dm2)) == 0
    end

    @testset "set_values!, zero_rows_local! and IS" begin
        # element dofs only: element (i, j) has global row (i - 1) + 3(j - 1) + 1
        dm = PETSc.DMStag(petsclib, comm, none(2), (3, 2), (0, 0, 1), 1)
        e = PETSc.element_location(dm)
        at(I...) = PETSc.stencil(dm, e, I)
        A = PETSc.PetscMat(dm)

        # duplicate columns sum under ADD_VALUES; values need not be a Vector
        PETSc.set_values!(A, dm, [at(1, 1), at(2, 1)], [at(1, 1), at(1, 1), at(2, 1)],
            1.0:6.0, LibPETSc.ADD_VALUES)
        PETSc.assemble!(A)
        @test A[1, 1] == 3.0
        @test A[1, 2] == 3.0
        @test A[2, 1] == 9.0
        @test A[2, 2] == 6.0
        @test_throws DimensionMismatch PETSc.set_values!(A, dm, [at(1, 1)], [at(1, 1)], [1.0, 2.0])

        @test PETSc.zero_rows_local!(A, dm, view([at(2, 1), at(1, 1)], 1:1), 5) === A
        @test A[2, 1] == 0.0
        @test A[2, 2] == 5.0
        @test PETSc.zero_rows_local!(A, dm, LibPETSc.DMStagStencil[]) === A

        v = PETSc.global_vec(dm)
        @test PETSc.set_values!(v, dm, [at(3, 2)], [7.0]) === v
        PETSc.assemble!(v)
        @test v[6] == 7.0
        @test_throws DimensionMismatch PETSc.set_values!(v, dm, [at(3, 2)], Float64[])

        is = LibPETSc.IS(dm, e => 0)
        @test PETSc.owns(is)
        @test length(is) == 6
        @test length(is) isa Int
        foreach(PETSc.destroy!, (is, v, A, dm))

        # a split of several locations and components
        is = LibPETSc.IS(dm2, PETSc.face_location(dm2, 1) => 0, PETSc.face_location(dm2, 2) => 0)
        @test length(is) == 4 * 2 + 3 * 3
        PETSc.destroy!(is)
    end

    @testset "set_values! from reused buffers" begin
        dm = PETSc.DMStag(petsclib, comm, none(2), (3, 2), (0, 0, 1), 1)
        e = PETSc.element_location(dm)
        at(I...) = PETSc.stencil(dm, e, I)
        rows = [at(1, 1), at(2, 1), at(3, 2)]
        vals = [1.0, 2.0, 3.0, 4.0, 5.0]

        # a prefix view writes the same block as the Vector it stands for
        A = PETSc.PetscMat(dm)
        B = PETSc.PetscMat(dm)
        PETSc.set_values!(A, dm, view(rows, 1:2), view(rows, 1:2), view(vals, 1:4))
        PETSc.set_values!(B, dm, rows[1:2], rows[1:2], vals[1:4])
        PETSc.assemble!(A)
        PETSc.assemble!(B)
        @test [A[i, j] for i in 1:2, j in 1:2] == [B[i, j] for i in 1:2, j in 1:2] == [1.0 2.0; 3.0 4.0]

        # and passes the buffer without a copy; a view that is not a prefix is copied
        set!(M, n) = PETSc.set_values!(M, dm, view(rows, 1:n), view(rows, 1:n),
            view(vals, 1:(n * n)), LibPETSc.ADD_VALUES)
        set!(A, 1)
        @test @allocated(set!(A, 2)) == 0
        w = PETSc.global_vec(dm)
        setv!(x) = PETSc.set_values!(x, dm, view(rows, 1:2), view(vals, 1:2))
        setv!(w)
        @test @allocated(setv!(w)) == 0
        PETSc.set_values!(w, dm, view(rows, 2:3), view(vals, 4:5))
        PETSc.assemble!(w)
        @test w[1] == 1.0
        @test w[2] == 4.0
        @test w[6] == 5.0
        foreach(PETSc.destroy!, (w, B, A, dm))
    end

    foreach(PETSc.destroy!, (dm3, dm2, dm1))
    PETSc.finalize(petsclib)
end
