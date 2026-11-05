using Test
using MPI
MPI.Initialized() || MPI.Init()
using PETSc

# On a periodic axis the STAR pattern wraps around. PETSc colors a periodic axis
# only when its two ends are on different ranks, so this needs several ranks.
@testset "star_fd_coloring on periodic axes" begin
    comm = MPI.COMM_WORLD
    nranks = MPI.Comm_size(comm)
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    P = PETSc.DM_BOUNDARY_PERIODIC
    N = PETSc.DM_BOUNDARY_NONE
    dof = 2
    procs = nranks == 4 ? (2, 2) : (nranks, 1)

    for (bx, by) in ((P, P), (P, N))
        da = PETSc.DMDA(petsclib, comm, (bx, by), (10, 10), dof;
                        stencil_width = 1, stencil_type = :star, processors = procs)
        c = PETSc.star_fd_coloring(da)
        # Each owned node couples to itself and its 4 neighbors, minus the ones
        # across a non-periodic edge
        @test sort!(reduce(vcat, c.coo_idxs_1b)) == 1:c.nnz_coo

        # The pattern the COO arrays describe, one per row
        A = PETSc.PetscMat(da)
        PETSc.LibPETSc.MatSetPreallocationCOOLocal(
            petsclib, A, petsclib.PetscInt(c.nnz_coo), c.row_coo_local_0b, c.col_coo_local_0b)
        PETSc.LibPETSc.MatSetValuesCOO(petsclib, A, ones(c.nnz_coo), PETSc.INSERT_VALUES)
        PETSc.assemble!(A)
        M = collect(A)   # the whole matrix, on every rank
        per_row = vec(count(!iszero, M; dims = 2))
        if by == P
            @test all(==(5dof), per_row)
        else
            # the rows of the bottom and top edge have no neighbor across it
            @test sort!(unique(per_row)) == [4dof, 5dof]
            @test count(==(4dof), per_row) == 2 * 10 * dof
        end
        @test (M .!= 0) == (M' .!= 0)
        foreach(PETSc.destroy!, (A, da))
    end

    da = PETSc.DMDA(petsclib, comm, (PETSc.DM_BOUNDARY_MIRROR, N), (10, 10), 1;
                    stencil_width = 1, stencil_type = :star, processors = procs)
    @test_throws "DM_BOUNDARY_MIRROR" PETSc.star_fd_coloring(da)
    PETSc.destroy!(da)
    PETSc.finalize(petsclib)
end
