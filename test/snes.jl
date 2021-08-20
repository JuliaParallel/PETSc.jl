using Test
using PETSc
using MPI
MPI.Initialized() || MPI.Init()

@testset "SNES" begin
    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)

    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt

        snes = PETSc.SNES(
            petsclib,
            comm;
            ksp_rtol = 1e-4,
            pc_type = "none",
            ksp_monitor = false,
            snes_monitor = false,
            snes_converged_reason = false,
            ksp_converged_reason = false,
        )

        r = PETSc.VecSeqWithArray(petsclib, zeros(PetscScalar, 2))
        PETSc.setfunction!(snes, r) do fx, snes, x
            PETSc.withlocalarray!(
                fx,
                x;
                read = (false, true),
                write = (true, false),
            ) do fx, x
                fx[1] = x[1]^2 + x[1] * x[2] - 3
                fx[2] = x[1] * x[2] + x[2]^2 - 6
            end
            return 0
        end

        J = PETSc.MatSeqDense(petsclib, zeros(PetscScalar, (2, 2)))
        PETSc.setjacobian!(snes, J) do J, snes, x
            PETSc.withlocalarray!(x; write = false) do x
                J[1, 1] = 2x[1] + x[2]
                J[1, 2] = x[1]
                J[2, 1] = x[2]
                J[2, 2] = x[1] + 2x[2]
            end
            PETSc.assemble!(J)
            return 0
        end
        x = PETSc.VecSeqWithArray(petsclib, PetscScalar.([2, 3]))
        b = PETSc.VecSeqWithArray(petsclib, PetscScalar.([0, 0]))
        PETSc.solve!(x, snes, b)
        PETSc.withlocalarray!(x; read = true, write = false) do x
            @test x ≈ [1, 2] rtol = 1e-4
        end

        # cleanup
        PETSc.destroy(x)
        PETSc.destroy(b)
        PETSc.destroy(J)
        PETSc.destroy(snes);

        PETSc.finalize(petsclib)
    end
end
