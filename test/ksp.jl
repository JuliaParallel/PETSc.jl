facts("\nTesting KSP") do
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

  kspg = PETSc.KSP(A, ksp_monitor="")

  println("performing ksp GMRES solve")
  x = kspg\b
  println("finished ksp solve")
  x_julia = A_julia\b_julia

  for i=1:3
    @fact x[i] --> roughly(x_julia[i])
  end

  println("A = ", A)
  println("b = ", b)
  println(" x = ", x)
  println("ksp info:\n",petscview(kspg))
  
  pc   = PETSc.PC(ST,comm=comm(kspg),pc_type="jacobi")
  PETSc.chk(PETSc.C.PCSetOperators(pc.p,A.p,A.p))
  kspg = PETSc.KSP(pc, ksp_monitor="")
  println("performing ksp GMRES solve with Jacobi preconditioner")
  x = kspg\b
  println("finished ksp solve")
  x_julia = A_julia\b_julia

  for i=1:3
    @fact x[i] --> roughly(x_julia[i])
  end

  println("A = ", A)
  println("b = ", b)
  println(" x = ", x)
  println("ksp info:\n",petscview(kspg))
  println("pc info:\n",petscview(pc))

  kspb = PETSc.KSP(A, ksp_type="bcgs", ksp_monitor="")
  println("performing ksp BCGS solve")
  x = kspb\b
  println("finished ksp solve")
  x_julia = A_julia\b_julia

  for i=1:3
    @fact x[i] --> roughly(x_julia[i])
  end

  println("A = ", A)
  println("b = ", b)
  println(" x = ", x)
  println("ksp info:\n",petscview(kspb))
end
