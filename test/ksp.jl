@testset "Testing KSP" begin
  A = PETSc.Mat(ST, 3,3)
  A_julia = zeros(ST,3,3)
  for i=1:3
    for j=1:3
      pos = (i-1)*3 + j
      val = RC(Complex(pos, pos))
      A[i, j] = val
      A_julia[i,j] = val
     end
  end

  A[3,2] = RC(Complex(6,6))  # make A non singular
  A_julia[3,2] = RC(Complex(6,6))

  b = PETSc.Vec(ST, 3, PETSc.C.VECMPI)
  b_julia = zeros(ST, 3)
  b[3] = RC(Complex(1,1))
  b_julia[3] = RC(Complex(1,1))

  @testset "ksp GMRES solves" begin
      kspg = PETSc.KSP(A, ksp_monitor="")
      x = kspg\b
      x_julia = A_julia\b_julia

      @test x ≈ x_julia

      println("ksp info:\n",petscview(kspg))
      
      pc   = PETSc.PC(ST,comm=comm(kspg),pc_type="jacobi")
      PETSc.chk(PETSc.C.PCSetOperators(pc.p,A.p,A.p))
      kspg = PETSc.KSP(pc, ksp_monitor="")
      x = kspg\b
      x_julia = A_julia\b_julia

      @test x ≈ x_julia
      println("ksp info:\n",petscview(kspg))
      println("pc info:\n",petscview(pc))
  end

  @testset "ksp BCGS solves" begin
      kspb = PETSc.KSP(A, ksp_type="bcgs", ksp_monitor="")
      x = kspb\b
      x_julia = A_julia\b_julia
      @test x ≈ x_julia
      println("ksp info:\n",petscview(kspb))
  end
end
