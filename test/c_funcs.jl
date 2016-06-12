@testset "C functions {$ST}" begin

 v_ptr = PETSc.C.VecCreate(ST)
 b = Vec{ST, :mpi}(v_ptr)
 resize!(b, mlocal=sys_size)
 global_indices = localpart(b) - 1  # zero based
  for i=1:sys_size
    idxm = [global_indices[i] ]
    val = [ rhs[i] ]
    PETSc.C.SetValues(b.p, idxm, val)
  end

  PETSc.C.AssemblyBegin(b.p)
  PETSc.C.AssemblyEnd(b.p)

   for i=1:sys_size
    idxm = [global_indices[i] ]
    val = ST[ 0.0 ]
    PETSc.C.GetValues(b.p, idxm, val)
    @test val[1] == rhs[i]
  end

  # matrix
  A = Mat(ST, mlocal=sys_size, nlocal=sys_size)
  row_range, col_range = localranges(A)

  for i=1:sys_size
    for j=1:sys_size
      row_idx = PETSc.C.PetscInt[ row_range[i] - 1 ]
      col_idx = PETSc.C.PetscInt[ row_range[j] - 1 ]
      val = ST[ A_julia[i, j] ]
      PETSc.C.SetValues(A.p, row_idx, col_idx, val)
    end
  end

  PETSc.C.AssemblyBegin(A.p)
  PETSc.C.AssemblyEnd(A.p)

  for i=1:sys_size
    for j=1:sys_size
      row_idx = PETSc.C.PetscInt[ row_range[i] - 1 ]
      col_idx = PETSc.C.PetscInt[ row_range[j] - 1 ]
      val = ST[ 0.0 ]
      PETSc.C.GetValues(A.p, row_idx, col_idx, val)
      @test val[1] == A_julia[i, j]
    end
  end

  idx = collect(row_range - 1)
  vals = collect(1:(sys_size*sys_size))
  vals2 = convert(Array{ST, 1}, vals)
  PETSc.C.SetValues(A.p, idx, idx, vals2)
  PETSc.C.AssemblyBegin(A.p)
  PETSc.C.AssemblyEnd(A.p)


  for i=1:sys_size
    for j=1:sys_size
      row_idx = PETSc.C.PetscInt[ row_range[i] - 1 ]
      col_idx = PETSc.C.PetscInt[ row_range[j] - 1 ]
      val = ST[ 0.0 ]
      PETSc.C.GetValues(A.p, row_idx, col_idx, val)
      @test val[1] == (i + (j-1)*sys_size)
    end
  end


  
  # block matrix
  B = Mat(ST, mlocal=sys_size, nlocal=sys_size, bs=sys_size, mtype=PETSc.C.MATMPIBAIJ)
  idx = PETSc.C.PetscInt[comm_rank]
  PETSc.C.SetValues(B.p, idx, idx, vals2)
  PETSc.C.AssemblyBegin(B.p)
  PETSc.C.AssemblyEnd(B.p)
  for i=1:sys_size
    for j=1:sys_size
      row_idx = PETSc.C.PetscInt[ row_range[i] - 1 ]
      col_idx = PETSc.C.PetscInt[ row_range[j] - 1 ]
      val = ST[ 0.0 ]
      PETSc.C.GetValues(A.p, row_idx, col_idx, val)
      @test val[1] == (i + (j-1)*sys_size)
    end
  end


  # shell matrix
  if ST == Float64
    ctx = (1, 2, 3)
    ctx_ptr = pointer_from_objref(ctx)
    c_ptr = PETSc.C.MatCreateShell(sys_size, sys_size, PETSc.C.PETSC_DETERMINE, PETSc.C.PETSC_DETERMINE, ctx_ptr)
    C = Mat{ST, PETSc.C.MATSHELL}(c_ptr)
   
     
    f_ptr = cfunction(mymult, PETSc.C.PetscErrorCode, (PETSc.C.Mat{ST}, PETSc.C.Vec{ST}, PETSc.C.Vec{ST}))
    PETSc.C.MatShellSetOperation(C.p, PETSc.C.MATOP_MULT, f_ptr)

    ctx_ret = PETSc.C.MatShellGetContext(C.p)
    @test unsafe_pointer_to_objref(ctx_ret) == ctx

    x = Vec(ST[1:sys_size;])
    xlocal = LocalVector(x)
    b = Vec(zeros(ST, sys_size))
    *(C, x, b)
    gc()  # avoid a finalizer problem
    blocal = LocalVector(b)
    for i=1:length(blocal)
      @test blocal[i] == ST(i*i)
    end

    LocalVectorRestore(blocal)
  end
end 
