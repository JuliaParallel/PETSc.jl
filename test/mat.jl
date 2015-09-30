facts("--- Testing Matrix Functions ---") do

mat = PETSc.Mat(Float64, 3, 4)

m,n = size(mat)

@fact m => 3
@fact n => 4

m,n = sizelocal(mat)

@fact m => 3
@fact n => 4

mtype = PETSc.gettype(mat)
@fact mtype => PETSc.C.MATMPIAIJ


# test set/get index
mat[1,1] = 3
mat[1,2] = 5
PETSc.assemble(mat)
val_ret = mat[1,1]
@fact val_ret => roughly(3.0)
mat[1,2] = 4
PETSc.assemble(mat)
println("inserted an additional value after final assembly")

mat2 = similar(mat)

m2, n2 = size(mat2)
@fact m2 => m
@fact n2 => n
@fact mat2[1,1] => not(mat[1,1])
mat3 = similar(mat, Float64, 4, 4)

m3, n3 = size(mat3)

@fact m3 => 4
@fact n3 => 4
@fact mat2[1,1] => not(mat[1,1])


mat4 = copy(mat)
@fact mat4[1,1] => roughly(mat[1, 1])

println("getting matrix info")
info = PETSc.getinfo(mat4)
bs = info.block_size
println("bs = ", bs)
@fact bs => 1


assembled = isassembled(mat4)
@fact assembled => true

println("mat = ", mat)
#println("mat3 = ", mat3)
println("size(mat) = ", size(mat))
println("inserting value into mat")
println("inserting value into mat3")
mat3.assembling = false
mat3[1,2] = 3
mat3[1,3] = 5
mat3.assembling = false
PETSc.assemble(mat3)

@fact mat3[1,2] => 3
@fact mat3[1,3] => 5

mat5 = Mat( Float64, 3, 3)

function increasing_diag()
  println("inserting values into mat5")

  (m,n) = size(mat5)
  println("m,n = ", m, ", ", n)
  dim = min(m, n)
  println("dim = ", dim)

  for i=1:dim
    println("i = ", i)
    mat5[i,i] = i
  end
end

assemble(increasing_diag, mat5)

(m,n) = size(mat5)
dim = min(m, n)

for i=1:dim
  @fact mat5[i,i] => roughly(i)
end

mat6 = PETSc.Mat(Float64, 3, 3)

vals = rand(3, 2)
idx = Array(1:3)
idy = Array(1:2)

mat6[idx, idy] = vals
assemble(mat6)
mat6j = zeros(3,3)
mat6j[1:3, 1:2] = vals
for i=1:3
  for j=1:3
    @fact mat6[i,j] => mat6j[i,j]
  end
end 

vals_ret = mat6[idx, idy]
@fact vals_ret => roughly(vals)

mat7 = PETSc.Mat(Float64, 3, 3)
vals = rand(3)
mat7[1, idx] = vals
assemble(mat7)
mat7j = zeros(3,3)
mat7j[1, idx] = vals
for i=1:3
  for j=1:3
    @fact mat7[i,j] => mat7j[i,j]
  end
end 

vals_ret = mat7[1, idx]
println("mat7 = ", mat7)
println("vals = ", vals)
println("vals_ret = ", vals_ret)
for i=1:3
  @fact vals_ret[i] => roughly(vals[i])
end

mat8 = PETSc.Mat(Float64, 3, 3)
vals = rand(3)
mat8[idx, 1] = vals
assemble(mat8)
mat8j = zeros(3,3)
mat8j[idx, 1] = vals


println("mat8 = ", mat8)
for i=1:3
  for j=1:3
    @fact mat8[i,j] => mat8j[i,j]
  end
end 

vals_ret = mat8[idx, 1]
@fact vals_ret => roughly(vals, atol= 1e-13)


mat9 = PETSc.Mat(Float64, 3, 3)
mat9[idx, idy] = 3
assemble(mat9)
mat9j = zeros(3,3)
mat9j[1:3, 1:2] = 3
for i=1:3
  for j=1:3
    @fact mat9[i,j] => mat9j[i,j]
  end
end 

mat10 = PETSc.Mat(Float64, 3, 3)
mat10[idx, 1] = 3
assemble(mat10)
mat10j = zeros(3,3)
mat10j[1:3, 1] = 3
for i=1:3
  for j=1:3
    @fact mat10[i,j] => mat10j[i,j]
  end
end 


mat11 = PETSc.Mat(Float64, 3, 3)
mat11[1, idy] = 3
assemble(mat11)
mat11j = zeros(3,3)
mat11j[1, 1:2] = 3
for i=1:3
  for j=1:3
    @fact mat11[i,j] => mat11j[i,j]
  end
end 






# test ranges and colon
mat12 = PETSc.Mat(Float64, 3, 3)
mat12[1:3, 1:2] = 3
assemble(mat12)
mat12j = zeros(3,3)
mat12j[1:3, 1:2] = 3
for i=1:3
  for j=1:3
    @fact mat12[i,j] => mat12j[i,j]
  end
end

mat13 = PETSc.Mat(Float64, 3, 3)
mat13[:, idy] = 3
assemble(mat13)
mat13j = zeros(3,3)
mat13j[:, 1:2] = 3
for i=1:3
  for j=1:3
    @fact mat13[i,j] => mat13j[i,j]
  end
end 


# test conversion of values to a new type

vals = [1, 2, 3]
mat14 = PETSc.Mat(Float64, 3, 3)
mat14[:, 1] = 3
assemble(mat14)
mat14j = zeros(3,3)
mat14j[:, 1] = 3
for i=1:3
  for j=1:3
    @fact mat14[i,j] => mat14j[i,j]
  end
end 

mat15 = PETSc.Mat(Float64, 3,3)
fill!(mat15, 1.0)
assemble(mat15)

for i=1:3
  for j=1:3
    @fact mat15[i,j] => roughly(1.0)
  end
end 


mat15jd = full(mat15)
for i=1:3
  for j=1:3
    @fact mat15[i,j] => roughly(mat15jd[i,j])
  end
end 

# it appears get/set values don't work on a transposed matrix
mat15_t = PETSc.MatTranspose(mat15)
#=
#assemble(mat15_t)
for i=1:3
  for j=1:3
    println("mat15_t[j,i] = ", mat15_t[j, i])
    @fact mat15_t[j, i] => roughly(mat15[i,j])
  end
end 
=#

mat16 = PETSc.Mat(Float64, 3, 3)
mat17 = PETSc.Mat(Float64, 3, 3)
mat16j = zeros(3, 3)
mat17j = zeros(3, 3)
vec1 = PETSc.Vec(Float64, 3)
vec1j = zeros(3)
cnt = 1
for i=1:3
  for j=1:3
    mat16[i,j] = cnt
    mat16j[i,j] = cnt
    mat17[i,j] = cnt + 9
    mat17j[i,j] = cnt + 9
    cnt += 1
  end
  vec1[i] = i
  vec1j[i] = i
end

assemble(mat16)
assemble(mat17)

vec2 = mat16*vec1
vec2j = mat16j*vec1j

for i=1:3
  @fact vec2[i] => vec2j[i]
end

mat18 = 2*mat16
assemble(mat18)
mat18j = 2*mat16j

println("mat18 = ", mat18)
println("mat18j = ", mat18j)
for i=1:3
  for j=1:3
    @fact mat18[i,j] => roughly(mat18j[i,j])
  end
end

#=
println("mat16 = ", mat16)
println("mat17 = ", mat17)
mat19 = mat16*mat17
mat19j = mat16j*mat17j
for i=1:3
  for j=1:3
    @fact mat19[i,j] => roughly(mat19j[i,j])
  end
end
=#

end
