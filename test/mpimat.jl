using Test
using MPI
MPI.Initialized() || MPI.Init()
using PETSc
using LinearAlgebra: mul!, norm

@testset "MatAIJ" begin
    comm = MPI.COMM_WORLD
    mpisize = MPI.Comm_size(comm)
    mpirank = MPI.Comm_rank(comm)

    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt

        loc_num_rows = 5
        loc_num_cols = 5
        diag_nonzeros = 3
        off_diag_non_zeros = 3

        mat = PETSc.MatAIJ(
            petsclib,
            comm,
            loc_num_rows,
            loc_num_cols,
            diag_nonzeros,
            off_diag_non_zeros,
        )

        # Get compatible vectors
        (right, left) = PETSc.createvecs(mat)

        # Fill the matrix and right vector
        row_rng = PETSc.ownershiprange(mat, false)
        for i in row_rng
            if i == 0
                vals = [-2, 1]
                row0idxs = [i]
                col0idxs = [i, i + 1]
            elseif i == mpisize * loc_num_rows - 1
                vals = [-2, 1]
                row0idxs = [i]
                col0idxs = [i, i - 1]
            else
                vals = [1, -2, 1]
                row0idxs = [i]
                col0idxs = [i - 1, i, i + 1]
            end
            PETSc.setvalues!(
                mat,
                PetscInt.(row0idxs),
                PetscInt.(col0idxs),
                PetscScalar.(vals),
            )
            right[i + 1] = i^3
        end
        PETSc.assemble!(mat)
        PETSc.assemble!(right)

        # Do matrix multiply and check result
        mul!(left, mat, right)
        for i in row_rng
            if i == 0
                v = -2 * i^3 + (i + 1)^3
            elseif i == mpisize * loc_num_rows - 1
                v = (i - 1)^3 - 2 * i^3
            else
                v = (i - 1)^3 - 2 * i^3 + (i + 1)^3
            end
            @test v == left[i + 1]
        end

        # Check the norm in parallel
        sz = loc_num_rows * mpisize
        exact_norm = sqrt(sz * 2^2 + 2 * (sz - 1))
        @test norm(mat) ≈ exact_norm

        PETSc.destroy(mat)
        PETSc.destroy(right)
        PETSc.destroy(left)

        PETSc.finalize(petsclib)
    end
end
nothing
