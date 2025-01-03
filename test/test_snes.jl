using Test
using PETSc, MPI, LinearAlgebra, SparseArrays

MPI.Initialized() || MPI.Init()

@testset "SNES" begin

  comm    = MPI.COMM_WORLD
  mpirank = MPI.Comm_rank(comm)
  mpisize = MPI.Comm_size(comm)

   # structure with which we can pass data to the user-routines above
   mutable struct Data
    vec
    julia
  end

  for petsclib in PETSc.petsclibs
      PETSc.initialize(petsclib)
      PetscScalar = PETSc.scalartype(petsclib)
      PetscInt    = PETSc.inttype(petsclib)

      # This shows some examples of how SNES can be used. 
      # Note that you will receive pointers to PETSc vectors within your user routines.
      #       That is important for parallel simulations, where global residual
      #       and solution vectors are passed to the user-routines, but we work with
      #       local vectors

      function fn!(cfx, cx, user_ctx)
        # We could do Global->Local here on cfx/cx, provided a pointer to the local
        #  vector is available in user_ctx
        #x_in  = PETSc.unsafe_localarray(PetscScalar, cx;  write=false)   # read array
        #fx_in = PETSc.unsafe_localarray(PetscScalar, cfx; write=true)    # write array
      #
        #fx_in[1] = x_in[1]^2       +  x_in[1]*x_in[2] - 3
        #fx_in[2] = x_in[1]*x_in[2] +  x_in[2]^2 - 6
        #
        #Base.finalize(fx_in)
        #Base.finalize(x_in)
        PETSc.withlocalarray!(
          cfx,
          cx;
          read = (false, true),
          write = (true, false),
        ) do fx, x
            @show x
            fx[1] = x[1]^2 + x[1] * x[2] - 3
            fx[2] = x[1] * x[2] + x[2]^2 - 6
        end

      end

      #=
      function update_jac!(cx, J1, args...)
          #x_in  = PETSc.unsafe_localarray(PetscScalar, cx;  write=false)
          PETSc.withlocalarray!(cx; write = false) do x_in
            J1[1,1] = 2x_in[1] + x_in[2]
            J1[1,2] = x_in[1]
            J1[2,1] = x_in[2]
            J1[2,2] = x_in[1] + 2x_in[2] 
          end
          #Base.finalize(x_in)
          PETSc.assemble!(J1)          
      end
      =#
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
    
     

      julia_vec = 0;  # we want pointers to local vectors 

      # NOTE: for some reason, snes_monitor and ksp_monitor is not working if we specify that here.
      # To be sorted out why that is (converged_reason does work)
      S = PETSc.SNES(
                        petsclib,
                        comm;
                        ksp_rtol = 1e-4,
                        pc_type = "none",
                        ksp_monitor = false,
                        snes_monitor = false,
                        snes_converged_reason = false,
                        ksp_converged_reason = false,
                    )

      data        = Data(PetscScalar.([100;2]), 1)
      #S.user_ctx  = data;      # we can pack anything we need in this struct

      PJ = PETSc.MatSeqDense(petsclib,zeros(PetscScalar,(2,2)))
      #PETSc.setfunction!(S, fn!, PETSc.VecSeq(petsclib, 2))
      PETSc.setfunction!(fn!, S, PETSc.VecSeq(petsclib, 2))
      PETSc.setjacobian!(S, update_jac!, PJ, PJ)
 
      x  = PETSc.VecSeqWithArray(petsclib,PetscScalar.([2.0, 3.0]));
      b  = PETSc.VecSeqWithArray(petsclib,PetscScalar.([0.0, 0.0]));
      PETSc.solve!(x, S, b)

      sol = PETSc.unsafe_localarray(PetscScalar, x.ptr)
      @test sol ≈ [1.0,2.0] rtol=1e-4

      # cleanup
      PETSc.destroy(x);
      PETSc.destroy(b);
      PETSc.destroy(PJ);
      #PETSc.destroy(S);

      PETSc.finalize(petsclib)
    end
end
