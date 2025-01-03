using Test
using PETSc, MPI, LinearAlgebra, SparseArrays
#PETSc.initialize()

petsclib = PETSc.petsclibs[1]
PETSc.initialize(petsclib)

@testset "Tests" begin
  m,n = 20,20
  x = randn(n)
  V = PETSc.VecSeqWithArray(petsclib,x)


  @test norm(x) ≈ norm(V) rtol=10eps()

  S = sprand(m,n,0.1) + I
  M = PETSc.MatSeqAIJ(petsclib,S)

  @test norm(S) ≈ norm(M) rtol=10eps()

  w = M*x
  @test w ≈ S*x

  ksp = PETSc.KSP(M; ksp_rtol=1e-8, pc_type="jacobi", ksp_monitor=false)
  #PETSc.settolerances!(ksp; rtol=1e-8)

  # TO BE FIXED
  ##@test PETSc.gettype(ksp) == "gmres" # default

  # TO BE FIXED
  #=
  pc = PETSc.PC(ksp)
  @test PETSc.nrefs(pc) == 2
  @test PETSc.gettype(pc) == "jacobi"

  # create an extra handle, check ref count is incremented
  pc_extra = PETSc.PC(ksp);
  @test PETSc.nrefs(pc) == 3
  # destroy extra handle, check ptr is set to null, ref count is decremented
  PETSc.destroy(pc_extra)
  @test pc_extra.ptr == C_NULL
  @test PETSc.nrefs(pc) == 2

  # safe to call destroy on null pointer
  PETSc.destroy(pc_extra)

  # set new pc, check ref counts are modified by ksp
  pc_new = PETSc.PC{Float64}(MPI.COMM_SELF);
  @test PETSc.nrefs(pc_new) == 1
  PETSc.setpc!(ksp, pc_new)
  @test PETSc.nrefs(pc_new) == 2
  @test PETSc.nrefs(pc) == 1

  PETSc.setpc!(ksp, pc)
  @test PETSc.nrefs(pc_new) == 1
  @test PETSc.nrefs(pc) == 2
  =#


  y = ksp \ w
  @test S*y ≈ w  rtol=1e-8


  w = M'*x
  @test w ≈ S'*x

  # to be fixed
#  y = ksp' \ w
#  @test S'*y ≈ w rtol=1e-8

  f!(y,x) = y .= 2 .*x

  M = PETSc.MatShell(petsclib,f!,MPI.COMM_SELF, 10,10)

  x = rand(10)

  @test M*x ≈ 2x
  @test PETSc.KSP(M) \ x ≈ x/2


  function F!(cfx, snes, cx)
    fx = PETSc.unsafe_localarray(cfx; read = false, write = true)
    x = PETSc.unsafe_localarray(cx;  read = true,  write = false)
 
    fx[1] = x[1]^2 + x[1] * x[2] - 3
    fx[2] = x[1] * x[2] + x[2]^2 - 6

    finalize(fx)
    finalize(x)
  end

  J = zeros(2,2)
  PJ = PETSc.MatSeqDense(petsclib,J)
  
  function updateJ!(J, snes, x)
    PETSc.unsafe_localarray(x; read = true, write = false)

    #PETSc.withlocalarray!(cx; write = false) do x
      J[1, 1] = 2x[1] + x[2]
      J[1, 2] = x[1]
      J[2, 1] = x[2]
      J[2, 2] = x[1] + 2x[2]
    #end
    PETSc.assemble!(J)
    Base.finalize(x)
    
   
  end

  S = PETSc.SNES(petsclib,MPI.COMM_SELF; ksp_rtol=1e-4, pc_type="none")
  r = PETSc.VecSeq(petsclib, 2)
  #PETSc.setfunction!(S, F!, S,r)
  
  PETSc.setfunction!(S, r) do cfx, snes, cx
      F!(cfx, snes, cx)
      return 0
  end


  PETSc.setjacobian!(S, PJ) do J, S, x
    PETSc.withlocalarray!(x; write = false) do x
        J[1, 1] = 2x[1] + x[2]
        J[1, 2] = x[1]
        J[2, 1] = x[2]
        J[2, 2] = x[1] + 2x[2]
    end
    PETSc.assemble!(J)
    return 0
  end

  #PETSc.setjacobian!(S, updateJ!, PJ, PJ)

  a = PETSc.VecSeqWithArray(petsclib,[2.0,3.0])
  @test PETSc.solve!(a, S) ≈ [1.0,2.0] rtol=1e-4
end
