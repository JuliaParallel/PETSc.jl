using PETSc

vec = PETSc.Vec(Float64, 5)

vec[3] = 7
display(vec)

mat = PETSc.Mat(Float64, 3, 4)

mat[2,3] = 7

println("mat[2,3] = ", mat[2,3])

display(mat)


PETSc.VecDestroy(vec)
PETSc.MatDestroy(mat)

PETSc.C.PetscFinalize(Float64)


