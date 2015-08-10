facts("\n   ---testing KSP solvers---") do
# create vectors and matrices


  b = PetscVec(comm);
  PetscVecSetType(b, "mpi");
  PetscVecSetSizes(b,sys_size, PetscInt(comm_size*sys_size));

  low, high = PetscVecGetOwnershipRange(b)
  b_global_indices = Array(low:PetscInt(high - 1))
 

  x = PetscVec(comm);
  PetscVecSetType(x,"mpi");
  PetscVecSetSizes(x,sys_size, PetscInt(comm_size*sys_size));

  low, high = PetscVecGetOwnershipRange(x)
  x_global_indices = Array(low:PetscInt(high - 1))
 

  for i=1:sys_size
    idxm = [ b_global_indices[i] ]   # index
    val = [ rhs[i] ]  # value
    PetscVecSetValues(b, idxm, val, PETSC_INSERT_VALUES)
  end

  PetscVecAssemblyBegin(b)
  PetscVecAssemblyEnd(b)



  A = PetscMat(comm)
  PetscMatSetType(A, "mpiaij")

  PetscMatSetSizes(A,sys_size,sys_size,PetscInt(comm_size*sys_size),PetscInt(comm_size*sys_size));
  PetscSetUp(A);

  low, high = PetscMatGetOwnershipRange(A)
  mat_global_indices = Array(low:PetscInt(high - 1))
 

  for i=1:sys_size
    for j = 1:sys_size
      idxm = [ mat_global_indices[i] ]  # row index
      idxn = [ mat_global_indices[j] ]  # column index
      PetscMatSetValues(A,idxm, idxn, [A_julia[i,j]],PETSC_INSERT_VALUES);
    end
  end

  PetscMatAssemblyBegin(A,PETSC_MAT_FINAL_ASSEMBLY);
  PetscMatAssemblyEnd(A,PETSC_MAT_FINAL_ASSEMBLY);



# perform solve
ksp = KSP(comm)
KSPSetOperators(ksp, A, A)
KSPSetFromOptions(ksp)
KSPSetUp(ksp)
KSPSolve(ksp, b, x)
reason = KSPGetConvergedReason(ksp)

println("KSP convergence reason = ", reason)
@fact reason => greater_than(0)  # convergence

PetscView(ksp)

# copy solution back to Julia
x_copy = zeros(PetscScalar, sys_size)
#idx = Array(0:2)
#idx = zeros(PetscInt, sys_size)
#for i=1:sys_size
#  idx[i] = i-1
#end

PetscVecGetValues(x, sys_size, x_global_indices, x_copy)
println("x_copy = ", x_copy)
println("x_julia = ", x_julia)
for i=1:sys_size
    @fact x_copy[i] => roughly(x_julia[i], atol=1e-14)
end

#PetscView(x, 0)

PetscDestroy(ksp)


println("   \n--- Testing LGMRES --- ")
# test using non default KSP method

# perform solve
ksp = KSP(comm)
KSPSetOperators(ksp, A, A)
KSPSetFromOptions(ksp)
KSPSetTolerances(ksp, rtol, abstol, dtol, maxits)
KSPSetInitialGuessNonzero(ksp, PetscBool(true))
KSPSetType(ksp, PETSc.KSPLGMRES)
KSPSetUp(ksp)
KSPSolve(ksp, b, x)
reason = KSPGetConvergedReason(ksp)
ksptype = KSPGetType(ksp)
println("finished calling KSPGetType")
println("typeof(ksptype) = ", typeof(ksptype))
println("ksptype = ", ksptype)

@fact ksptype => PETSc.KSPLGMRES


rtol_ret, abstol_ret, dtol_ret, maxits_ret = KSPGetTolerances(ksp)

@fact rtol_ret => roughly(rtol)
@fact abstol_ret => roughly(abstol)
@fact dtol_ret => roughly(dtol)
@fact maxits_ret => maxits
@fact KSPGetInitialGuessNonzero(ksp) => true
println("KSP convergence reason = ", reason)
@fact reason => greater_than(0)  # convergence

rnorm = KSPGetResidualNorm(ksp)
@fact rnorm => less_than(abstol)
PetscView(ksp)

# copy solution back to Julia
x_copy = zeros(PetscScalar, sys_size)
#idx = Array(0:2)
#idx = zeros(PetscInt, sys_size)
#for i=1:sys_size
#  idx[i] = i-1
#end

PetscVecGetValues(x, sys_size, x_global_indices, x_copy)
println("x_copy = ", x_copy)
println("x_julia = ", x_julia)
for i=1:sys_size
    @fact x_copy[i] => roughly(x_julia[i], atol=1e-14)
end





PetscDestroy(x)
PetscDestroy(b)
PetscDestroy(A)
PetscDestroy(ksp)
end
