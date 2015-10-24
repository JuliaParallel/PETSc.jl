facts("\ntesting KSP") do
  A = PETSc.Mat(ST, 3,3)
  A_julia = zeros(ST,3, 3)
  for i=1:3
    for j=1:3
      pos = (i-1)*3 + j
      val = RC(Complex(pos, pos))
      A[i, j] = val
      A_julia[i,j] = val
     end
  end

  A[3,2] = RC(Complex(6, 6))  # make A non singular
  A_julia[3,2] = RC(Complex(6,6))

  PETSc.AssemblyBegin(A, PETSc.C.MAT_FINAL_ASSEMBLY)

  PETSc.AssemblyEnd(A, PETSc.C.MAT_FINAL_ASSEMBLY)

  b = PETSc.Vec(ST, 3, PETSc.C.VECMPI)
  b_julia = zeros(ST, 3)
  b[3] = RC(Complex(1, 1))
  b_julia[3] = RC(Complex(1,1))

  ksp = PETSc.KSP(A)

  x = PETSc.solve(ksp, b)
  x_julia = A_julia\b_julia

  for i=1:3
    @fact x[i] => roughly(x_julia[i])
  end

  println("A = ", A)
  println("b = ", b)
  println(" x = ", x)
end
