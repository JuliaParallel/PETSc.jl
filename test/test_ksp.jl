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

# test PC
println("   \n---Testing PC Options ---")
pc = KSPGetPC(ksp)

PCSetType(pc, PETSc.PCBJACOBI)
pctype = PCGetType(pc)
println("pctype = ", pctype)
@fact pctype => PETSc.PCBJACOBI
#PCFactorSetUseInPlace(pc, PetscBool(true))

### Do some KSP setup

KSPSetFromOptions(ksp)
KSPSetTolerances(ksp, rtol, abstol, dtol, maxits)
KSPSetInitialGuessNonzero(ksp, PetscBool(true))
KSPSetType(ksp, PETSc.KSPLGMRES)
KSPSetUp(ksp)

### More PC testing
println("getting sub KSP objects")
n_local, first_local, ksp_arr = PCBJacobiGetSubKSP(pc)
ksp2 = ksp_arr[1]

println("ksp2 = ", ksp2)
println("typeof(ksp2) = ", typeof(ksp2))
println("getting sub KSP PC")
pc2 = KSPGetPC(ksp2)
println("finsihed getting sub KSP PC")

println("setting ReusePreconditioner")
PCSetReusePreconditioner(pc2, PetscBool(true))
@fact PCGetReusePreconditioner(pc2) => true
println("finished setting ReusePreconditioner")

PCFactorSetAllowDiagonalFill(pc2, PetscBool(true))
@fact PCFactorGetAllowDiagonalFill(pc2) => true

PCFactorSetLevels(pc2, PetscInt(1))
@fact PCFactorGetLevels(pc2) => 1  # should be pc2

PCSetReusePreconditioner(pc2, PetscBool(true))
@fact PCGetReusePreconditioner(pc2) => true

PCFactorSetFill(pc, PetscReal(7.0))

#=
PCJacobiSetType(pc, PETSc.PC_JACOBI_ROWMAX)
@fact PCJacobiGetType(pc) => PETSc.PC_JACOBI_ROWMAX
=#








#=
tmp = PCFactorGetUseInPlace(pc)
println("tmp = ", tmp)
@fact PCFactorGetUseInPlace(pc) => true
=#


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
