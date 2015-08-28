using PETSc

vec = PETSc.Vec(Float64, 5)

vec[3] = 7
display(vec)
print("\n")

mat = PETSc.Mat(Float64, 3, 4)

mat[2,3] = 7

println("mat[2,3] = ", mat[2,3])

display(mat)
print("\n")

PETSc.VecDestroy(vec)
PETSc.MatDestroy(mat)

PETSc.C.PetscFinalize(Float64)


