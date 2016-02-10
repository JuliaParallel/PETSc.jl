@testset "KSP{$ST}" begin
  A = PETSc.Mat(ST, 3,3)
  println("A = ", A)
  println("typeof(A) = ", typeof(A))
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
      kspg = PETSc.KSP(A)
      x = kspg\b
      x_julia = A_julia\b_julia
      @test x ≈ x_julia

      pc   = PETSc.PC(ST,comm=comm(kspg),pc_type="jacobi")
      PETSc.chk(PETSc.C.PCSetOperators(pc.p,A.p,A.p))
      kspg = PETSc.KSP(pc)
      x = kspg\b
      x_julia = A_julia\b_julia

      @test x ≈ x_julia
      # fixme: test petscview without side effects, e.g.
      # via PetscViewerStringOpen
  end

  @testset "ksp BCGS solves" begin
      kspb = PETSc.KSP(A, ksp_type="bcgs")
      x = kspb\b
      x_julia = A_julia\b_julia
      @test x ≈ x_julia
  end
end
