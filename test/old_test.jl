using Test
using PETSc, MPI, LinearAlgebra, SparseArrays
#PETSc.initialize()

petsclib = PETSc.petsclibs[1]
PETSc.initialize(petsclib)

MPI.Initialized() || MPI.Init()
# Windows PETSc binaries are built without MPI support, use PETSC_COMM_SELF instead
comm = Sys.iswindows() ? LibPETSc.PETSC_COMM_SELF : MPI.COMM_WORLD

@testset "Old Tests" begin
  m,n = 20,20
  x = randn(n)
  #V = PETSc.VecSeqWithArray(petsclib,x)
  
  PetscScalar     = petsclib.PetscScalar
  PetscInt        = petsclib.PetscInt

  bs = PetscInt(1)
  V = LibPETSc.VecCreateSeqWithArray(petsclib,comm, bs, PetscInt(length(x)), x)

  @test norm(x) ≈ norm(V) rtol=10eps()

  S = sprand(m,n,0.1) + I
  #M = PETSc.MatSeqAIJ(petsclib,S)
  M = PETSc.MatCreateSeqAIJ(petsclib, comm, S)
       

  @test norm(S) ≈ norm(M) rtol=10eps()
  vec_x = LibPETSc.VecCreateSeqWithArray(petsclib, LibPETSc.PETSC_COMM_SELF, PetscInt(1), PetscInt(length(x)), copy(x))

  w = M*vec_x
  @test w[:] ≈ S*x

  ksp = PETSc.KSP(M; ksp_rtol=1e-8, pc_type="jacobi", ksp_monitor=false)
  #PETSc.settolerances!(ksp; rtol=1e-8)

  @test PETSc.type(ksp) == "gmres" # default

  y = ksp \ w
  @test S*y[:] ≈ w[:]  rtol=1e-6


  #w = M'*x
  #@test w ≈ S'*x

  # to be fixed
#  y = ksp' \ w
#  @test S'*y ≈ w rtol=1e-8

  f!(y,x) = y .= 2 .*x

  M = PETSc.MatShell(petsclib,f!, Sys.iswindows() ? LibPETSc.PETSC_COMM_SELF : MPI.COMM_SELF, 10,10)

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

    return 0
  end

  J = zeros(2,2)
  PJ = LibPETSc.MatCreateSeqDense(petsclib, comm, PetscInt(size(J,1)), PetscInt(size(J,2)), J[:])
        

  
  function updateJ!(J, snes, x)
    PETSc.unsafe_localarray(x; read = true, write = false)

    J[1, 1] = 2x[1] + x[2]
    J[1, 2] = x[1]
    J[2, 1] = x[2]
    J[2, 2] = x[1] + 2x[2]

    PETSc.assemble!(J)
    Base.finalize(x)
    
    return 0  # needs to be zero!
  end

  S = PETSc.SNES(petsclib, Sys.iswindows() ? LibPETSc.PETSC_COMM_SELF : MPI.COMM_SELF; ksp_rtol=1e-4, pc_type="none")
  r = LibPETSc.VecCreateSeqWithArray(petsclib, LibPETSc.PETSC_COMM_SELF, PetscInt(1), PetscInt(2), zeros(PetscScalar, 2))

  
  # You can do this to set the callback functions.
  # Please be aware that the functions above MUST return 0
  PETSc.setfunction!(F!, S, r)
  PETSc.setjacobian!(updateJ!, S, PJ, PJ)

  # The alternative is to do this (tested below):
  #PETSc.setfunction!(S, r) do cfx, snes, cx
  #    F!(cfx, snes, cx)
  #    return 0
  #end

  #PETSc.setjacobian!(S, PJ) do J, S, x
  #  updateJ!(J, S, x)
  #  return 0
  #end

  #a = PETSc.VecSeqWithArray(petsclib,[2.0,3.0])
  a = LibPETSc.VecCreateSeqWithArray(petsclib, LibPETSc.PETSC_COMM_SELF, PetscInt(1), PetscInt(2), PetscScalar.([2.0,3.0]))

  sol = PETSc.solve!(a, S)
  @test sol[:] ≈ [1.0,2.0] rtol=1e-4


  # Test the alternative is to do this:
  PETSc.setfunction!(S, r) do cfx, snes, cx
      F!(cfx, snes, cx)
      return 0
  end

  PETSc.setjacobian!(S, PJ) do J, S, x
    updateJ!(J, S, x)
    return 0
  end
  sol1 = PETSc.solve!(a, S)
  @test sol1[:] ≈ [1.0,2.0] rtol=1e-4

  PETSc.destroy(S)
  PETSc.destroy(M)
  PETSc.destroy(ksp)
  PETSc.destroy(V)
  PETSc.destroy(r)

end

PETSc.finalize(petsclib)