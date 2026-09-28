# test/dmstag_views.jl
# with_field_views!: views of DMStag local vectors by field, their types and
# allocations, and that every vector is handed back.

using Test
using PETSc
using MPI

MPI.Initialized() || MPI.Init()

@testset "with_field_views!" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    comm = MPI.COMM_SELF
    LibPETSc = PETSc.LibPETSc
    none(n) = ntuple(_ -> PETSc.DM_BOUNDARY_NONE, n)
    stag(n, dofs) = PETSc.DMStag(petsclib, comm, none(length(n)), n, dofs, 1)

    @testset "views address the points stencil names" begin
        dm = stag((3, 2), (1, 1, 1))
        flow = (PETSc.face_location(dm, 1) => 0, PETSc.face_location(dm, 2) => 0,
                PETSc.element_location(dm) => 0)
        x = PETSc.global_vec(dm)
        I = CartesianIndex(2, 1)
        PETSc.set_values!(x, dm, [PETSc.stencil(dm, flow[1].first, I),
                                  PETSc.stencil(dm, flow[3].first, I)], [5.0, 7.0])
        PETSc.assemble!(x)
        l = PETSc.local_vec(dm)
        PETSc.global_to_local!(l, dm, x)

        got = PETSc.with_field_views!(dm, l; fields = flow, write = false) do (Vx, Vy, P)
            (Vx[I], P[I], Vy[I], axes(P) == axes(Vx))
        end
        @test got == (5.0, 7.0, 0.0, true)

        # the whole array, indexed [I..., slot]
        slot = PETSc.dof_slot(dm, PETSc.element_location(dm), 0)
        @test PETSc.with_field_views!(A -> A[I, slot], dm, l) == 7.0

        # writing through a view changes the vector
        PETSc.with_field_views!(dm, l; fields = flow[3:3]) do (P,)
            P[I] = 9.0
        end
        @test PETSc.with_field_views!(A -> A[I, slot], dm, l; write = false) == 9.0
        foreach(PETSc.destroy!, (l, x, dm))
    end

    @testset "types, allocations, several vectors" begin
        function checkout(n)
            dm = stag((n, n), (1, 1, 1))
            flow = (PETSc.face_location(dm, 1) => 0, PETSc.face_location(dm, 2) => 0,
                    PETSc.element_location(dm) => 0)
            vs = ntuple(_ -> PETSc.local_vec(dm), 4)
            body(a, b, c, d) = a[3][1, 1] + d[1][2, 2]
            run() = PETSc.with_field_views!(body, dm, vs...; fields = flow,
                write = (false, true, true, false))
            run()
            allocations = @allocations run()
            concrete = PETSc.with_field_views!(dm, vs...; fields = flow) do args...
                all(t -> all(v -> isconcretetype(typeof(v)), t), args)
            end
            foreach(PETSc.destroy!, vs)
            PETSc.destroy!(dm)
            return allocations, concrete
        end
        small, concrete = checkout(8)
        large, _ = checkout(64)
        @test concrete
        @test small == large          # independent of the grid size
        @test small <= 4 * 20         # the budget: about 20 per vector

        for (n, dofs) in (((4,), (1, 1)), ((2, 2, 2), (1, 1, 1, 1)))
            dm = stag(n, dofs)
            l = PETSc.local_vec(dm)
            fields = (PETSc.vertex_location(dm) => 0, PETSc.element_location(dm) => 0)
            @test PETSc.with_field_views!(dm, l; fields) do (V, E)
                isconcretetype(typeof(V)) && ndims(E) == length(n)
            end
            foreach(PETSc.destroy!, (l, dm))
        end
    end

    @testset "errors hand the vectors back" begin
        dm = stag((3, 2), (0, 0, 1))
        l = PETSc.local_vec(dm)
        g = PETSc.global_vec(dm)
        @test_throws ErrorException PETSc.with_field_views!(A -> error("inside"), dm, l)
        # a vector still checked out could not be checked out for writing again
        @test PETSc.with_field_views!(A -> (A[1, 1, 1] = 1.0; true), dm, l)

        dm2 = stag((3, 3), (0, 0, 1))
        @test_throws ArgumentError PETSc.with_field_views!(A -> nothing, dm, l, g)
        l2 = PETSc.local_vec(dm2)
        @test_throws ArgumentError PETSc.with_field_views!(A -> nothing, dm, l2)
        foreach(PETSc.destroy!, (l2, dm2, g, l, dm))
    end

    PETSc.finalize(petsclib)
end
