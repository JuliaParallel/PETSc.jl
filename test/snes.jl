using Test
using PETSc
using MPI
if !Sys.iswindows()
    MPI.Initialized() || MPI.Init()
end

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

        # setconvergencetest! — custom Julia-closure convergence test
        snes4 = PETSc.SNES(
            petsclib,
            comm;
            ksp_rtol = 1e-4,
            pc_type = "none",
            ksp_monitor = false,
            snes_monitor = false,
            snes_converged_reason = false,
            ksp_converged_reason = false,
        )
        r4 = LibPETSc.VecCreateSeqWithArray(petsclib, comm, PetscInt(1), PetscInt(2), zeros(PetscScalar, 2))
        PETSc.setfunction!(snes4, r4) do fx, snes, x
            PETSc.withlocalarray!(fx, x; read = (false, true), write = (true, false)) do fx, x
                fx[1] = x[1]^2 + x[1] * x[2] - PetscScalar(3)
                fx[2] = x[1] * x[2] + x[2]^2 - PetscScalar(6)
            end
            return PetscInt(0)
        end
        J4 = LibPETSc.MatCreateSeqDense(petsclib, comm, PetscInt(2), PetscInt(2), zeros(PetscScalar, 4))
        PETSc.setjacobian!(snes4, J4) do J, snes, x
            PETSc.withlocalarray!(x; write = false) do x
                J[1, 1] = 2x[1] + x[2]
                J[1, 2] = x[1]
                J[2, 1] = x[2]
                J[2, 2] = x[1] + 2x[2]
            end
            PETSc.assemble!(J)
            return PetscInt(0)
        end

        ntest_calls = Ref(0)
        seen_its = Int[]
        solupdate_ptrs_nonnull = Ref(true)
        solupdate_norms = Float64[]
        PETSc.setconvergencetest!(snes4) do snes, it, xnorm, gnorm, fnorm
            ntest_calls[] += 1
            push!(seen_its, it)

            # Regression check for a bug where SNESGetSolution/SNESGetSolutionUpdate's
            # wrappers discarded the C-returned vector pointer and set it to C_NULL instead
            # (`x.ptr = C_NULL` rather than `x.ptr = x_[]`), so any subsequent use of the
            # "returned" vector (e.g. VecGetArrayRead, as a GeoTech2D-style dU/dtol check
            # would do) operated on a null pointer and crashed. Exercise both accessors at
            # every iteration after the first (mirroring how a real dtol/step-size check
            # would use them) and confirm the vector is usable.
            if it > 0
                du = LibPETSc.SNESGetSolutionUpdate(petsclib, snes)
                solupdate_ptrs_nonnull[] &= (du.ptr != C_NULL)
                duarr = LibPETSc.VecGetArrayRead(petsclib, du)
                push!(solupdate_norms, maximum(abs, duarr))
                LibPETSc.VecRestoreArrayRead(petsclib, du, duarr)

                xsol = LibPETSc.SNESGetSolution(petsclib, snes)
                solupdate_ptrs_nonnull[] &= (xsol.ptr != C_NULL)
                xarr = LibPETSc.VecGetArrayRead(petsclib, xsol)
                LibPETSc.VecRestoreArrayRead(petsclib, xsol, xarr)
            end

            # a custom criterion in the same spirit as GeoTech2D's fres < rtol: converged
            # once the residual is small, otherwise keep iterating (never diverge here)
            return fnorm < 1e-6 ? LibPETSc.SNES_CONVERGED_FNORM_ABS : LibPETSc.SNES_CONVERGED_ITERATING
        end

        x4 = LibPETSc.VecCreateSeqWithArray(petsclib, comm, PetscInt(1), PetscInt(2), PetscScalar.([2, 3]))
        b4 = LibPETSc.VecCreateSeqWithArray(petsclib, comm, PetscInt(1), PetscInt(2), PetscScalar.([0, 0]))
        PETSc.solve!(x4, snes4, b4)

        @test x4[:] ≈ [1, 2] rtol = 1e-4
        @test ntest_calls[] > 0                 # the Julia closure was actually invoked
        @test issorted(seen_its)                # called once per iteration, in order
        @test maximum(seen_its) > 0             # the it>0 branch (SNESGetSolutionUpdate/Solution) ran
        @test solupdate_ptrs_nonnull[]           # neither accessor returned a null Vec
        @test !isempty(solupdate_norms) && all(isfinite, solupdate_norms)
        @test LibPETSc.SNESGetConvergedReason(petsclib, snes4) == LibPETSc.SNES_CONVERGED_FNORM_ABS

        # a convergence test that immediately reports divergence must produce that reason
        snes5 = PETSc.SNES(
            petsclib,
            comm;
            ksp_rtol = 1e-4,
            pc_type = "none",
            ksp_monitor = false,
            snes_monitor = false,
            snes_converged_reason = false,
            ksp_converged_reason = false,
        )
        r5 = LibPETSc.VecCreateSeqWithArray(petsclib, comm, PetscInt(1), PetscInt(2), zeros(PetscScalar, 2))
        PETSc.setfunction!(snes5, r5) do fx, snes, x
            PETSc.withlocalarray!(fx, x; read = (false, true), write = (true, false)) do fx, x
                fx[1] = x[1]^2 + x[1] * x[2] - PetscScalar(3)
                fx[2] = x[1] * x[2] + x[2]^2 - PetscScalar(6)
            end
            return PetscInt(0)
        end
        J5 = LibPETSc.MatCreateSeqDense(petsclib, comm, PetscInt(2), PetscInt(2), zeros(PetscScalar, 4))
        PETSc.setjacobian!(snes5, J5) do J, snes, x
            PETSc.withlocalarray!(x; write = false) do x
                J[1, 1] = 2x[1] + x[2]
                J[1, 2] = x[1]
                J[2, 1] = x[2]
                J[2, 2] = x[1] + 2x[2]
            end
            PETSc.assemble!(J)
            return PetscInt(0)
        end
        PETSc.setconvergencetest!(snes5) do snes, it, xnorm, gnorm, fnorm
            return LibPETSc.SNES_DIVERGED_LOCAL_MIN
        end
        x5 = LibPETSc.VecCreateSeqWithArray(petsclib, comm, PetscInt(1), PetscInt(2), PetscScalar.([2, 3]))
        b5 = LibPETSc.VecCreateSeqWithArray(petsclib, comm, PetscInt(1), PetscInt(2), PetscScalar.([0, 0]))
        PETSc.solve!(x5, snes5, b5)
        @test LibPETSc.SNESGetConvergedReason(petsclib, snes5) == LibPETSc.SNES_DIVERGED_LOCAL_MIN
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
     
        PETSc.destroy(x4)
        PETSc.destroy(b4)
        PETSc.destroy(r4)
        PETSc.destroy(J4)

        PETSc.destroy(x5)
        PETSc.destroy(b5)
        PETSc.destroy(r5)
        PETSc.destroy(J5)

        PETSc.destroy(snes)
        PETSc.destroy(snes2)
        PETSc.destroy(snes3)
        PETSc.destroy(snes4)
        PETSc.destroy(snes5)

        PETSc.finalize(petsclib)
        
    end
end
