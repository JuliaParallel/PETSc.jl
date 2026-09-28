# test/mpi_dmstag_halos.jl
# Run on several ranks: local_to_local! in place fills each rank's ghost
# points from the owned values its neighbours changed in their local vectors.

using Test
using MPI
MPI.Initialized() || MPI.Init()
using PETSc

@testset "local_to_local! across ranks" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    comm = MPI.COMM_WORLD
    nranks = MPI.Comm_size(comm)

    # 1D, one element dof, two elements per rank
    dm = PETSc.DMStag(petsclib, comm, (PETSc.DM_BOUNDARY_NONE,), (2nranks,), (0, 1), 1)
    owned = PETSc.corners(dm)
    ghosted = PETSc.ghost_corners(dm)
    element(k) = ghosted.lower[1] + k - 1       # global element of local entry k
    isowned(k) = owned.lower[1] <= element(k) <= owned.upper[1]

    g = PETSc.global_vec(dm)
    PETSc.with_local_array!(g) do a
        a .= owned.lower[1]:owned.upper[1]
    end
    l = PETSc.local_vec(dm)
    PETSc.global_to_local!(l, dm, g)

    # negate the owned entries only, then refresh the ghosts in place
    PETSc.with_local_array!(l) do a
        for k in eachindex(a)
            isowned(k) && (a[k] = -a[k])
        end
    end
    PETSc.local_to_local!(l, dm)
    PETSc.with_local_array!(l) do a
        # DMStag pads the last rank's ghost region past the grid; those entries stay unset
        @test all(a[k] == -element(k) for k in eachindex(a) if element(k) <= 2nranks)
    end

    foreach(PETSc.destroy!, (l, g, dm))
    PETSc.finalize(petsclib)
end
