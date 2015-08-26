import MPI
using mod_PETSc
MPI.Init()

comm = MPI.COMM_WORLD
comm_size = MPI.Comm_size(comm)
comm_rank = MPI.Comm_rank(comm)

petsc_args =  ["-ksp_monitor","-malloc","-malloc_debug","-malloc_dump"]

mod_PETSc.PetscInitialize(petsc_args)
sys_size = 3
sys = [1.0 2 3; 4 5 7; 7 8 9]
rhs = [1., 2, 3]

syst = sys.'

# create b
tmp = Array(mod_PETSc.Vec{Float64}, 1)
mod_PETSc.VecCreate( comm, tmp)
b = tmp[1]
mod_PETSc.VecSetType(b, mod_PETSc.VECMPI)
mod_PETSc.VecSetSizes( b, sys_size, sys_size)

# create x
tmp = Array(mod_PETSc.Vec{Float64}, 1)
mod_PETSc.VecCreate( comm, tmp)
x = tmp[1]
mod_PETSc.VecSetType(x, mod_PETSc.VECMPI)
mod_PETSc.VecSetSizes( x, sys_size, sys_size)


# populate b
idx = Int32[0:(sys_size - 1);]
mod_PETSc.VecSetValues(b, sys_size, idx, rhs, mod_PETSc.INSERT_VALUES)

mod_PETSc.VecAssemblyBegin(b)
mod_PETSc.VecAssemblyEnd(b)

#=
# create PetscViewer
tmp = Array(PetscViewer,1)
PetscViewerCreate(comm, tmp)
viewer = tmp[1]
=#

# verify b is populated
vals_ret = zeros(sys_size)
mod_PETSc.VecGetValues(b, sys_size, idx, vals_ret)

println("vals_ret = ", vals_ret) 


tmp = Array(mod_PETSc.Mat{Float64}, 1)
mod_PETSc.MatCreate(comm, tmp)
A = tmp[1]

mod_PETSc.MatSetType(A, mod_PETSc.MATMPIAIJ)
mod_PETSc.MatSetSizes(A, sys_size, sys_size, sys_size, sys_size)
mod_PETSc.MatSetUp(A)

mod_PETSc.MatSetValues(A, sys_size, idx, sys_size, idx, sys, mod_PETSc.INSERT_VALUES)

println("finished setting values")

mod_PETSc.MatAssemblyBegin(A, mod_PETSc.MAT_FINAL_ASSEMBLY)
mod_PETSc.MatAssemblyEnd(A, mod_PETSc.MAT_FINAL_ASSEMBLY)

A_ret = zeros(sys_size, sys_size)

mod_PETSc.MatGetValues(A, sys_size, idx, sys_size, idx, A_ret)

println("A_ret = ", A_ret)


# do KSP Solve
tmp = Array(mod_PETSc.KSP{Float64}, 1)
mod_PETSc.KSPCreate(comm, tmp)
ksp = tmp[1]
mod_PETSc.KSPSetOperators(ksp, A, A)
mod_PETSc.KSPSetFromOptions(ksp)
mod_PETSc.KSPSetUp(ksp)
mod_PETSc.KSPSolve(ksp, b, x)

tmp = Array(mod_PETSc.KSPConvergedReason, 1)
mod_PETSc.KSPGetConvergedReason(ksp, tmp)
reason = tmp[1]
println("ksp convergence reason = ", reason)

x_ret = zeros(sys_size)

mod_PETSc.VecGetValues(x, sys_size, idx, x_ret)

println("x_ret = ", x_ret)

x_julia = syst\rhs
println("x_julia = ", x_julia)


