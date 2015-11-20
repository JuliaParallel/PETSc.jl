facts("--- Testing Matrix Functions ---") do

mat = PETSc.Mat(ST, 3, 4)

@fact size(mat) --> (3,4)

@fact sizelocal(mat) --> (3,4)
@fact lengthlocal(mat) --> 12

mtype = PETSc.gettype(mat)
@fact mtype --> PETSc.C.MATMPIAIJ


# test set/get index
vt1 = RC(complex(3., 3.))
vt2 = RC(complex(5., 5.))
mat[1,1] = vt1
mat[1,2] = vt2
PETSc.assemble(mat)
val_ret = mat[1,1]
@fact val_ret --> roughly(vt1)
vt1 = RC(complex(4., 4.))
mat[1,2] = vt1
PETSc.assemble(mat)
println("inserted an additional value after final assembly")

mat2 = similar(mat)

@fact size(mat2) --> (3,4)
@fact mat2[1,1] --> not(mat[1,1])
mat3 = similar(mat, ST, 4, 4)

@fact size(mat3) --> (4,4)
@fact mat2[1,1] --> not(mat[1,1])
@fact_throws ArgumentError resize!(mat2)
@fact_throws ArgumentError resize!(mat2,5,mlocal=2)

mat4 = copy(mat)
@fact mat4[1,1] --> roughly(mat[1, 1])

mat5 = copy(mat)
@fact conj(conj(mat5))[1,1] --> roughly(mat[1,1])

println("getting matrix info")
@fact PETSc.getinfo(mat4).block_size --> 1

@fact isassembled(mat4) --> true

println("mat = ", mat)
println("size(mat) = ", size(mat))
println("inserting value into mat")
println("inserting value into mat3")
mat3.assembling = false
vt1 = RC(complex(3., 3.))
vt2 = RC(complex(5., 5.))
mat3[1,2] = vt1
mat3[1,3] = vt2
mat3.assembling = false
PETSc.assemble(mat3)

@fact mat3[1,2] --> vt1
@fact mat3[1,3] --> vt2

mat5 = Mat( ST, 3, 3)

function increasing_diag()
  println("inserting values into mat5")

  (m,n) = size(mat5)
  println("m,n = ", m, ", ", n)
  dim = min(m, n)
  println("dim = ", dim)

  for i=1:dim
    println("i = ", i)
    i_float = Float64(i)
    mat5[i,i] = RC(complex(i_float, i_float))
  end
end

assemble(increasing_diag, mat5)

(m,n) = size(mat5)
dim = min(m, n)

for i=1:dim
  i_float = Float64(i)
  @fact mat5[i,i] --> roughly(RC(complex(i_float, i_float)))
end

mat6 = PETSc.Mat(ST, 3, 3)

vals = RC(complex(rand(3, 2), rand(3,2)))
idx = Array(1:3)
idy = Array(1:2)

mat6[idx, idy] = vals
assemble(mat6)
mat6j = zeros(ST, 3,3)
mat6j[1:3, 1:2] = vals
for i=1:3, j=1:3
  @fact mat6[i,j] --> mat6j[i,j]
end

@fact mat6[idx,idy] --> roughly(vals)

mat7 = PETSc.Mat(ST, 3, 3)
vals = RC( complex(rand(3), rand(3)))
mat7[1, idx] = vals
assemble(mat7)
mat7j = zeros(ST, 3,3)
mat7j[1, idx] = vals
for i=1:3, j=1:3
  @fact mat7[i,j] --> mat7j[i,j]
end

vals_ret = mat7[1, idx]
println("mat7 = ", mat7)
println("vals = ", vals)
println("vals_ret = ", vals_ret)
for i=1:3
  @fact vals_ret[i] --> roughly(vals[i])
end

mat8 = PETSc.Mat(ST, 3, 3)
vals = RC(complex(rand(3), rand(3)))
mat8[idx, 1] = vals
assemble(mat8)
mat8j = zeros(ST, 3,3)
mat8j[idx, 1] = vals

println("mat8 = ", mat8)
for i=1:3, j=1:3
  @fact mat8[i,j] --> mat8j[i,j]
end

vals_ret = mat8[idx, 1]
@fact vals_ret --> roughly(vals, atol= 1e-13)

mat9 = PETSc.Mat(ST, 3, 3)
vt = RC(complex(3., 3.))
mat9[idx, idy] = vt
assemble(mat9)
mat9j = zeros(ST, 3,3)
mat9j[1:3, 1:2] = vt
for i=1:3, j=1:3
  @fact mat9[i,j] --> mat9j[i,j]
end

mat10 = PETSc.Mat(ST, 3, 3)
mat10[idx, 1] = vt
assemble(mat10)
mat10j = zeros(ST, 3,3)
mat10j[1:3, 1] = vt
for i=1:3, j=1:3
  @fact mat10[i,j] --> mat10j[i,j]
end

mat11 = PETSc.Mat(ST, 3, 3)
mat11[1, idy] = vt
assemble(mat11)
mat11j = zeros(ST, 3,3)
mat11j[1, 1:2] = vt
for i=1:3, j=1:3
  @fact mat11[i,j] --> mat11j[i,j]
end

context("test ranges and colon") do
    mat12 = PETSc.Mat(ST, 3, 3)
    mat12[1:3, 1:2] = vt
    assemble(mat12)
    mat12j = zeros(ST, 3,3)
    mat12j[1:3, 1:2] = vt
    for i=1:3, j=1:3
      @fact mat12[i,j] --> mat12j[i,j]
    end

    mat13 = PETSc.Mat(ST, 3, 3)
    mat13[:, idy] = vt
    assemble(mat13)
    mat13j = zeros(ST, 3,3)
    mat13j[:, 1:2] = vt
    for i=1:3, j=1:3
      @fact mat13[i,j] --> mat13j[i,j]
    end
end

context("test conversion of values to a new type") do

    vals = [1, 2, 3]
    mat14 = PETSc.Mat(ST, 3, 3)
    mat14[:, 1] = vt
    assemble(mat14)
    mat14j = zeros(ST, 3,3)
    mat14j[:, 1] = vt
    for i=1:3, j=1:3
       @fact mat14[i,j] --> mat14j[i,j]
    end

    vt = RC(complex(1.,1))
    mat15 = PETSc.Mat(ST, 3,3)
    fill!(mat15, vt)
    assemble(mat15)

    for i=1:3, j=1:3
      @fact mat15[i,j] --> roughly(vt)
    end

    mat15jd = full(mat15)
    for i=1:3, j=1:3
      @fact mat15[i,j] --> roughly(mat15jd[i,j])
    end

    # it appears get/set values don't work on a transposed matrix
    mat15_t = PETSc.MatTranspose(mat15)
    #=
    #assemble(mat15_t)
    for i=1:3
      for j=1:3
        println("mat15_t[j,i] = ", mat15_t[j, i])
        @fact mat15_t[j, i] --> roughly(mat15[i,j])
      end
    end
    =#

    mat16 = PETSc.Mat(ST, 3, 3)
    mat17 = PETSc.Mat(ST, 3, 3)
    mat16j = zeros(ST, 3, 3)
    mat17j = zeros(ST, 3, 3)
    vec1 = PETSc.Vec(ST, 3)
    vec1j = zeros(ST, 3)
    cnt = 1
    for i=1:3
      for j=1:3
        cnt_f = RC(complex(Float64(cnt), Float64(cnt)))
        cnt_f2 = RC(complex(Float64(cnt + 9), Float64(cnt + 9)))

        mat16[i,j] = cnt_f
        mat16j[i,j] = cnt_f
        mat17[i,j] = cnt_f2
        mat17j[i,j] = cnt_f2
        cnt += 1
      end
      vec1[i] = RC(complex(Float64(i), i))
      vec1j[i] = RC(complex(Float64(i), i))
    end

    assemble(mat16)
    assemble(mat17)

    vec2 = mat16*vec1
    vec2j = mat16j*vec1j

    for i=1:3
      @fact vec2[i] --> vec2j[i]
    end

    mat18 = 2*mat16
    assemble(mat18)
    mat18j = 2*mat16j

    println("mat18 = ", mat18)
    println("mat18j = ", mat18j)
    for i=1:3, j=1:3
      @fact mat18[i,j] --> roughly(mat18j[i,j])
    end

    println("mat16 = ", mat16)
    println("mat17 = ", mat17)
    mat19 = mat16*mat17
    mat19j = mat16j*mat17j
    for i=1:3, j=1:3
      @fact mat19[i,j] --> roughly(mat19j[i,j])
    end

    mat20= mat16+mat17
    mat20j = mat16j+mat17j
    for i=1:3, j=1:3
      @fact mat20[i,j] --> roughly(mat20j[i,j])
    end

    mat21= mat16-mat17
    mat21j = mat16j-mat17j
    for i=1:3, j=1:3
      @fact mat21[i,j] --> roughly(mat21j[i,j])
    end
    
    mat22 = mat16/2
    mat22j = mat16j/2

    for i=1:3
      for j=1:3
        @fact mat22[i,j] --> roughly(mat22j[i,j])
      end
    end

    mat23 = 2\mat16
    mat23j = 2\mat16j

    for i=1:3
      for j=1:3
        @fact mat23[i,j] --> roughly(mat23j[i,j])
      end
    end

    mat24  = -mat16
    mat24j = -mat16j
    for i=1:3
      for j=1:3
        @fact mat24[i,j] --> roughly(mat24j[i,j])
      end
    end
end

end
