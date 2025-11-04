using Test
using PETSc
using MPI
MPI.Initialized() || MPI.Init()

#@testset "SNES" begin
    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)

    #for petsclib in PETSc.petsclibs
        petsclib = PETSc.petsclibs[1]
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

        r = LibPETSc.VecCreateSeqWithArray(petsclib,comm, 1, 2, zeros(PetscScalar, 2))

        #=
        PETSc.setfunction!(snes, r) do fx, snes, x
            @show x, fx, snes
            PETSc.withlocalarray!(
                fx,
                x;
                read = (false, true),
                write = (true, false),
            ) do fx, x
            @show x, fx
                fx[1] = x[1]^2 + x[1] * x[2] - 3
                fx[2] = x[1] * x[2] + x[2]^2 - 6
            end
            return 0
        end
        =#

        function fn!(cfx, snes, cx)
            #@show cfx, snes, cx
           # @show cfx[:], cx[:]
            #=
            cfx[1] = cx[1]^2 + cx[1] * cx[2] - 3
            cfx[2] = cx[1] * cx[2] + cx[2]^2 - 6
            =#
            
            
            PETSc.withlocalarray!(
            cfx, cx;
            read = (false, true),
            write = (true, false),
            ) do fx, x
                fx[1] = x[1]^2 + x[1] * x[2] - 3
                fx[2] = x[1] * x[2] + x[2]^2 - 6
            #    @show x, fx
            end
            
            return 0
        end
        PETSc.setfunction!(snes, fn!, r)
        

        #J = PETSc.MatSeqDense(petsclib, zeros(PetscScalar, (2, 2)))
        J = LibPETSc.MatCreateSeqDense(petsclib,comm, 2, 2, zeros(4))

        #=
        PETSc.setjacobian!(snes, J) do J, snes, x
            @show J, snes
           # @show x
           #=
            PETSc.withlocalarray!(x; write = false) do x
                J[1, 1] = 2x[1] + x[2]
                J[1, 2] = x[1]
                J[2, 1] = x[2]
                J[2, 2] = x[1] + 2x[2]
            end
            PETSc.assemble!(J)
            =#
            return 0
        end
        =#

        function jacobian!(J, snes, cx)
            #@show J, snes, cx[:]
            J[1, 1] = 2x[1] + x[2]
            J[1, 2] = x[1]
            J[2, 1] = x[2]
            J[2, 2] = x[1] + 2x[2]

            PETSc.assemble!(J)
            return 0
        end
        PETSc.setjacobian!(jacobian!, snes, J)



        x = LibPETSc.VecCreateSeqWithArray(petsclib,comm, 1, 2, PetscScalar.([2, 3]))
        b = LibPETSc.VecCreateSeqWithArray(petsclib,comm, 1, 2, PetscScalar.([0, 0]))
       
      
        PETSc.solve!(x, snes, b)
        
        @test x[:] ≈ [1, 2] rtol = 1e-4
        
        #=
        
        PETSc.withlocalarray!(x; read = true, write = false) do x
            @show x
        end

        # cleanup
        PETSc.destroy(x)
        PETSc.destroy(b)
        PETSc.destroy(J)
        PETSc.destroy(snes);

        PETSc.finalize(petsclib)
        =#
    #end
#end
