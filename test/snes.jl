using Test
using PETSc
using MPI
MPI.Initialized() || MPI.Init()

@testset "SNES" begin
    comm = Sys.iswindows() ? LibPETSc.PETSC_COMM_SELF : MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)

    for petsclib in PETSc.petsclibs
        #@show petsclib
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt

        # Note: there are multiple ways to set the function and Jacobian
        # This is method 1 using withlocalarray! to access the local vector  
        # See below for other methods
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

        r = LibPETSc.VecCreateSeqWithArray(petsclib,comm, PetscInt(1), PetscInt(2), zeros(PetscScalar, 2))
        function fn!(cfx, snes, cx)
            PETSc.withlocalarray!(
                    cfx, cx;
                    read = (false, true),
                    write = (true, false),
            ) do fx, x
                fx[1] = x[1]^2 + x[1] * x[2] - PetscScalar(3)
                fx[2] = x[1] * x[2] + x[2]^2 - PetscScalar(6)
            end
            
            return PetscInt(0)
        end
        PETSc.setfunction!(snes, fn!, r)
        
       function jacobian!(J, snes, x)
            PETSc.withlocalarray!(x; write = false) do x
                J[1, 1] = 2x[1] + x[2]
                J[1, 2] = x[1]
                J[2, 1] = x[2]
                J[2, 2] = x[1] + 2x[2]
            end
            PETSc.assemble!(J)
            return PetscInt(0)
        end
        J = LibPETSc.MatCreateSeqDense(petsclib,comm, PetscInt(2), PetscInt(2), zeros(PetscScalar,4))
        PETSc.setjacobian!(jacobian!, snes, J)

        x = LibPETSc.VecCreateSeqWithArray(petsclib,comm, PetscInt(1), PetscInt(2), PetscScalar.([2, 3]))
        b = LibPETSc.VecCreateSeqWithArray(petsclib,comm, PetscInt(1), PetscInt(2), PetscScalar.([0, 0]))

        PETSc.solve!(x, snes, b)
        
        @test x[:] ≈ [1, 2] rtol = 1e-4

        # ----------------------------------------------------------------
        
        # Method 2 - use index notation in residual and jacobian functions
        snes2 = PETSc.SNES(
            petsclib,
            comm;
            ksp_rtol = 1e-4,
            pc_type = "none",
            ksp_monitor = false,
            snes_monitor = false,
            snes_converged_reason = false,
            ksp_converged_reason = false,
        )

        r2 = LibPETSc.VecCreateSeqWithArray(petsclib,comm, PetscInt(1), PetscInt(2), zeros(PetscScalar, 2))
       
     
        # use local indices - this may be slower (or allocate more); to be tested
        function fn2!(fx, snes, x)

            fx[1] = x[1]^2 + x[1] * x[2] - 3
            fx[2] = x[1] * x[2] + x[2]^2 - 6
  
            return PetscInt(0)
        end
        PETSc.setfunction!(snes2, fn2!, r2)

        function jacobian2!(J, snes2, x)
            J[1, 1] = 2x[1] + x[2]
            J[1, 2] = x[1]
            J[2, 1] = x[2]
            J[2, 2] = x[1] + 2x[2]

            PETSc.assemble!(J)
            return PetscInt(0)
        end
        J2 = LibPETSc.MatCreateSeqDense(petsclib,comm, PetscInt(2), PetscInt(2), zeros(PetscScalar,4))
        PETSc.setjacobian!(jacobian2!, snes2, J2)


        # 
        x2 = LibPETSc.VecCreateSeqWithArray(petsclib,comm, PetscInt(1), PetscInt(2), PetscScalar.([2, 3]))
        b2 = LibPETSc.VecCreateSeqWithArray(petsclib,comm, PetscInt(1), PetscInt(2), PetscScalar.([0, 0]))


        PETSc.solve!(x2, snes2, b2)
        
        
        @test x2[:] ≈ [1, 2] rtol = 1e-4
        # ----------------------------------------------------------------
        
        
        # Method 3 - use "do" to set residual and jacobian functions, for the ones of you that like this style
        snes3 = PETSc.SNES(
            petsclib,
            comm;
            ksp_rtol = 1e-4,
            pc_type = "none",
            ksp_monitor = false,
            snes_monitor = false,
            snes_converged_reason = false,
            ksp_converged_reason = false,
        )

        r3 = LibPETSc.VecCreateSeqWithArray(petsclib,comm, PetscInt(1), PetscInt(2), zeros(PetscScalar, 2))
        PETSc.setfunction!(snes3, r3) do fx, snes, x
            fx[1] = x[1]^2 + x[1] * x[2] - 3
            fx[2] = x[1] * x[2] + x[2]^2 - 6
            return PetscInt(0)
        end


        J3 = LibPETSc.MatCreateSeqDense(petsclib,comm, PetscInt(2), PetscInt(2), zeros(PetscScalar,4))
        PETSc.setjacobian!(snes3, J3) do J, snes, x
            J[1, 1] = 2x[1] + x[2]
            J[1, 2] = x[1]
            J[2, 1] = x[2]
            J[2, 2] = x[1] + 2x[2]

            PETSc.assemble!(J)
            return PetscInt(0)
        end

        # 
        x3 = LibPETSc.VecCreateSeqWithArray(petsclib,comm, PetscInt(1), PetscInt(2), PetscScalar.([2, 3]))
        b3 = LibPETSc.VecCreateSeqWithArray(petsclib,comm, PetscInt(1), PetscInt(2), PetscScalar.([0, 0]))


        PETSc.solve!(x3, snes3, b3)
        
        
        @test x3[:] ≈ [1, 2] rtol = 1e-4
        # ----------------------------------------------------------------
        

        # cleanup
        PETSc.destroy(x)
        PETSc.destroy(b)
        PETSc.destroy(r)
        PETSc.destroy(J)
        
        PETSc.destroy(x2)
        PETSc.destroy(b2)
        PETSc.destroy(r2)
        PETSc.destroy(J2)
     
        PETSc.destroy(x3)
        PETSc.destroy(b3)
        PETSc.destroy(r3)
        PETSc.destroy(J3)
     
        PETSc.destroy(snes) 
        PETSc.destroy(snes2) 
        PETSc.destroy(snes3) 

        PETSc.finalize(petsclib)
        
    end
end
