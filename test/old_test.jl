using Test
using PETSc, MPI, LinearAlgebra, SparseArrays
#PETSc.initialize()
PETSc.initialize(PETSc.petsclibs[1])

@testset "Tests" begin
  m,n = 20,20
  x = randn(n)
  V = PETSc.VecSeq(x)


  @test norm(x) ≈ norm(V) rtol=10eps()

  S = sprand(m,n,0.1) + I
  M = PETSc.MatSeqAIJ(S)

  @test norm(S) ≈ norm(M) rtol=10eps()

  w = M*x
  @test w ≈ S*x

  ksp = PETSc.KSP(M; ksp_rtol=1e-8, pc_type="jacobi", ksp_monitor=false)
  #PETSc.settolerances!(ksp; rtol=1e-8)

  @test PETSc.gettype(ksp) == "gmres" # default

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


  y = ksp \ w
  @test S*y ≈ w  rtol=1e-8


  w = M'*x
  @test w ≈ S'*x

  y = ksp' \ w
  @test S'*y ≈ w rtol=1e-8

  f!(y,x) = y .= 2 .*x

  M = PETSc.MatShell{Float64}(f!,10,10)

  x = rand(10)

  @test M*x ≈ 2x
  @test PETSc.KSP(M) \ x ≈ x/2


  function F!(cfx, cx, a)
    x = PETSc.unsafe_localarray(Float64,cx,  write=false, read=true)
    fx = PETSc.unsafe_localarray(Float64,cfx,write=true, read=true)
    fx[1] = x[1]^2 + x[1]*x[2] - 3
    fx[2] = x[1]*x[2] + x[2]^2 - 6
    Base.finalize(x)
    Base.finalize(fx)
  end

  J = zeros(2,2)
  PJ = PETSc.MatSeqDense(J)
  function updateJ!(cx, args...)
    x = PETSc.unsafe_localarray(Float64,cx, write=false, read=true)
    J[1,1] = 2x[1] + x[2]
    J[1,2] = x[1]
    J[2,1] = x[2]
    J[2,2] = x[1] + 2x[2]
    Base.finalize(x)
  end

  S = PETSc.SNES{Float64}(PETSc.petsclibs[1],MPI.COMM_SELF; ksp_rtol=1e-4, pc_type="none")
  PETSc.setfunction!(S, F!, PETSc.VecSeq(zeros(2)))
  PETSc.setjacobian!(S, updateJ!, PJ, PJ)
  a = PETSc.VecSeq([2.0,3.0])
  @test PETSc.solve!(a, S) ≈ [1.0,2.0] rtol=1e-4
end
