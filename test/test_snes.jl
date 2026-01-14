using Test
using PETSc, MPI, LinearAlgebra, SparseArrays

MPI.Initialized() || MPI.Init()

@testset "test_snes" begin

  comm = Sys.iswindows() ? LibPETSc.PETSC_COMM_SELF : MPI.COMM_WORLD
  mpirank = MPI.Comm_rank(comm)
  mpisize = MPI.Comm_size(comm)

   # structure with which we can pass data to the user-routines above
   mutable struct Data
    vec
    julia
  end

  for petsclib in PETSc.petsclibs
      #petsclib = PETSc.petsclibs[1]
      PETSc.initialize(petsclib)
      PetscScalar = PETSc.scalartype(petsclib)
      PetscInt    = PETSc.inttype(petsclib)

      # This shows some examples of how SNES can be used. 
      # Note that you will receive pointers to PETSc vectors within your user routines.
      #       That is important for parallel simulations, where global residual
      #       and solution vectors are passed to the user-routines, but we work with
      #       local vectors

      function fn!(fx, snes, x, user_ctx)
        fx[1] = x[1]^2 + x[1] * x[2] - 3
        fx[2] = x[1] * x[2] + x[2]^2 - 6
        
        return PetscInt(0) # success code of PETSc
      end


      function update_jac!(J1, snes, x, user_ctx)
        J1[1,1] = 2x[1] + x[2]
        J1[1,2] = x[1]
        J1[2,1] = x[2]
        J1[2,2] = x[1] + 2x[2]
        PETSc.assemble!(J1)
        return PetscInt(0)
      end

     

      julia_vec = 0;  # we want pointers to local vectors 

      # NOTE: for some reason, snes_monitor and ksp_monitor is not working if we specify that here.
      # To be sorted out why that is (converged_reason does work)
      S = PETSc.SNES(petsclib,comm;  
                      ksp_rtol=1e-4, 
                      pc_type="none", 
                      ksp_monitor=false, 
                      snes_monitor=false, 
                      snes_converged_reason=false, 
                      ksp_converged_reason=false)

      data        = Data(PetscScalar.([100;2]), 1)
      S.user_ctx  = data;      # we can pack anything we need in this struct

      #PJ = PETSc.MatSeqDense(zeros(PetscScalar,(2,2)))
      PJ = LibPETSc.MatCreateSeqDense(petsclib,comm, PetscInt(2), PetscInt(2), zeros(PetscScalar,4))

      r = LibPETSc.VecCreateSeqWithArray(petsclib,comm, PetscInt(1), PetscInt(2), zeros(PetscScalar, 2))
      PETSc.setfunction!(S, fn!, r)
      PETSc.setjacobian!(S, update_jac!, PJ, PJ)
 
      x = LibPETSc.VecCreateSeqWithArray(petsclib,comm, PetscInt(1), PetscInt(2), PetscScalar.([2, 3]))
      b = LibPETSc.VecCreateSeqWithArray(petsclib,comm, PetscInt(1), PetscInt(2), PetscScalar.([0, 0]))
      PETSc.solve!(x, S, b)

      #sol = PETSc.unsafe_localarray(PetscScalar, x.ptr)
      sol = x[:]
      @test sol ≈ [1.0,2.0] rtol=1e-4

      # cleanup - destroy SNES first, then the vectors/matrices it references
      PETSc.destroy(S);
      PETSc.destroy(x);
      PETSc.destroy(b);
      PETSc.destroy(r);
      PETSc.destroy(PJ);

      PETSc.finalize(petsclib)
      GC.gc()  # Force garbage collection to clean up any remaining objects
    end
end
