function make_mat(dims=(3,4))
    mat = PETSc.Mat(ST, dims...)
    vt1 = RC(complex(3., 3.))
    vt2 = RC(complex(5., 5.))
    mat[1,1] = vt1
    mat[1,2] = vt2
    PETSc.assemble(mat)
    return mat
end

@testset "Testing Matrix Functions" begin

  @testset "Preallocator" begin
    for mt in [PETSc.C.MATMPIAIJ,PETSc.C.MATMPIBAIJ,PETSc.C.MATMPISBAIJ]
      @test_throws ArgumentError PETSc.Mat(ST, 3, 4, mtype=mt, nnz=collect(1:10))
      @test_throws ArgumentError PETSc.Mat(ST, 3, 4, mtype=mt, onnz=collect(1:10))
    end
    for mt in [PETSc.C.MATBLOCKMAT,PETSc.C.MATSEQAIJ,PETSc.C.MATSEQBAIJ,PETSc.C.MATSEQSBAIJ]
      @test_throws ArgumentError PETSc.Mat(ST, 3, 4, mtype=mt, nnz=collect(1:10))
    end
  end

  vt1 = RC(complex(3., 3.))
  vt2 = RC(complex(5., 5.))
  @testset "Utility functions" begin
    mat = make_mat()

    @test size(mat) == (3,4)
    @test sizelocal(mat) == (3,4)
    @test lengthlocal(mat) == 12
    val_ret = mat[1,1]
    vt1 = RC(complex(3., 3.))
    @test val_ret ≈ vt1
    vt1 = RC(complex(4., 4.))
    mat[1,2] = vt1
    PETSc.assemble(mat)

    mtype = PETSc.gettype(mat)
    @test mtype == PETSc.C.MATMPIAIJ
    # test set/get index
    vt1 = RC(complex(3., 3.))
    vt2 = RC(complex(5., 5.))
    mat[1,1] = vt1
    mat[1,2] = vt2
    PETSc.assemble(mat)
    val_ret = mat[1,1]
    @test val_ret ≈ vt1
    vt1 = RC(complex(4., 4.))
    mat[1,2] = vt1
    PETSc.assemble(mat)

    #test nnz
    @test nnz(mat) == 2
  end
  @testset "real and imag" begin
    mat  = make_mat()
    rmat = PETSc.Mat(ST, 3, 4)
    rmat[1,1] = RC(complex(3., 0.))
    rmat[1,2] = RC(complex(5., 0.))
    PETSc.assemble(rmat)
    @test real(mat)[1,1] == rmat[1,1]
    @test real(mat)[1,2] == rmat[1,2]
    if ST <: Complex
        @test imag(mat)[1,1] == rmat[1,1]
        @test imag(mat)[1,2] == rmat[1,2]
    end
  end
  @testset "diag and trace" begin
    dmat = similar(make_mat(),3,3)
    dmat[1,1] = vt1
    assemble(dmat)
    d = diag(dmat)
    @test d[1] == vt1
    @test trace(dmat) == vt1
  end
  @testset "similar and resize" begin
    mat = similar(make_mat())
    @test size(mat) == (3,4)
    @test mat[1,1] != vt1
    @test size(similar(mat, ST)) == (3,4)
    @test size(similar(mat, 4, 4)) == (4,4)
    @test size(similar(mat, (4, 4))) == (4,4)
    @test mat[1,1] != make_mat()[1,1]
    @test_throws ArgumentError resize!(mat)
    @test_throws ArgumentError resize!(mat,5,mlocal=2)
  end
  @testset "copy and conj" begin
    mat = similar(make_mat((4,4)))
    @test size(mat) == (4,4)
    mat2 = copy(mat)
    @test mat2[1,1] ≈ mat[1, 1]
    @test conj(conj(mat))[1,1] ≈ mat[1,1]
  end
  @testset "getting Mat info, inserting and assembling" begin
    mat = make_mat()
    @test PETSc.getinfo(mat).block_size == 1
    @test isassembled(mat)

    mat2 = similar(mat, ST, 4, 4)
    mat2.assembling = false
    mat2[1,2] = vt1
    mat2[1,3] = vt2
    mat2.assembling = false
    PETSc.assemble(mat2)

    @test mat2[1,2] == vt1
    @test mat2[1,3] == vt2
    mat3 = Mat( ST, 3, 3)

    function increasing_diag()
      (m,n) = size(mat3)
      dim = min(m, n)
      for i=1:dim
        i_float = Float64(i)
        mat3[i,i] = RC(complex(i_float, i_float))
      end
    end

    assemble(increasing_diag, mat3)

    (m,n) = size(mat3)
    dim = min(m, n)

    for i=1:dim
      i_float = Float64(i)
      @test mat3[i,i] ≈ RC(complex(i_float, i_float))
    end
  end
  @testset "transpose and transpose!" begin
    mat = make_mat((3,3))
    ctmat = copy(mat)
    @test transpose!(transpose!(copy(ctmat))) == mat
    @test transpose(transpose(ctmat)) == mat
  end
  vt = RC(complex(3., 3.))
  @testset "array indexing" begin 
    vals = RC(complex(rand(3, 2), rand(3,2)))
    idx = Array(1:3)
    idy = Array(1:2)
    @testset "sub indexing" begin
      mat = PETSc.Mat(ST, 3, 3)
      mat[idx, idy] = vals
      assemble(mat)
      matj = zeros(ST, 3,3)
      matj[1:3, 1:2] = vals
      @test mat == matj
      @test mat[idx,idy] ≈ vals
    end
    @testset "y indexing" begin
      mat = PETSc.Mat(ST, 3, 3)
      vals = RC( complex(rand(3), rand(3)))
      mat[1, idx] = vals
      assemble(mat)
      matj = zeros(ST, 3,3)
      matj[1, idx] = vals
      @test mat == matj

      vals_ret = mat[1, idx]
      @test vals_ret.' ≈ vals[idx] 
    end
    @testset "x indexing" begin
      mat = PETSc.Mat(ST, 3, 3)
      vals = RC(complex(rand(3), rand(3)))
      mat[idx, 1] = vals
      assemble(mat)
      matj = zeros(ST, 3,3)
      matj[idx, 1] = vals
      @test mat == matj
      vals_ret = mat[idx, 1]
      @test vals_ret ≈ vals
    end
    @testset "x,y set and fetch" begin 
      mat = PETSc.Mat(ST, 3, 3)
      mat[idx, idy] = vt
      assemble(mat)
      matj = zeros(ST, 3,3)
      matj[1:3, 1:2] = vt
      @test mat == matj
    end
    @testset "x set and fetch" begin 
      mat = PETSc.Mat(ST, 3, 3)
      mat[idx, 1] = vt
      assemble(mat)
      matj = zeros(ST, 3,3)
      matj[1:3, 1] = vt
      @test mat == mat
    end
    @testset "y set and fetch" begin 
      mat = PETSc.Mat(ST, 3, 3)
      mat[1, idy] = vt
      assemble(mat)
      matj = zeros(ST, 3,3)
      matj[1, 1:2] = vt
      @test mat == matj
    end
  end
  @testset "test ranges and colon" begin
    idy = Array(1:2)
    @testset "submatrix" begin
      mat = PETSc.Mat(ST, 3, 3)
      mat[1:3, 1:2] = vt
      assemble(mat)
      matj = zeros(ST, 3,3)
      matj[1:3, 1:2] = vt
      @test mat == matj
    end
    @testset "on an axis with range" begin
      mat = PETSc.Mat(ST, 3, 3)
      mat[:, idy] = vt
      assemble(mat)
      matj = zeros(ST, 3,3)
      matj[:, 1:2] = vt
      @test mat == matj
    end
    @testset "on a column" begin
      vals = [1, 2, 3]
      mat = PETSc.Mat(ST, 3, 3)
      mat[:, 1] = vt
      assemble(mat)
      matj = zeros(ST, 3,3)
      matj[:, 1] = vt
      @test mat == matj
    end
  end
    
  @testset "full and fill" begin
    vt = RC(complex(1.,1.))
    mat = PETSc.Mat(ST, 3, 3)
    fill!(mat, vt)
    assemble(mat)
    @test mat == fill(vt,(3,3))
    matjd = full(mat)
    @test mat == matjd
  end
  @testset "test conversion of values to a new type" begin
    mata = PETSc.Mat(ST, 3, 3)
    matb = PETSc.Mat(ST, 3, 3)
    mataj = zeros(ST, 3, 3)
    matbj = zeros(ST, 3, 3)
    vec = PETSc.Vec(ST, 3)
    vecj = zeros(ST, 3)
    cnt = 1
    for i=1:3
      for j=1:3
        cnt_f = RC(complex(Float64(cnt), Float64(cnt)))
        cnt_f2 = RC(complex(Float64(cnt + 9), Float64(cnt + 9)))

        mata[i,j] = cnt_f
        mataj[i,j] = cnt_f
        matb[i,j] = cnt_f2
        matbj[i,j] = cnt_f2
        cnt += 1
      end
      vec[i] = RC(complex(Float64(i), i))
      vecj[i] = RC(complex(Float64(i), i))
    end

    assemble(mata)
    assemble(matb)

    @testset "matrix-vector product" begin
      result = mata*vec
      resultj = mataj*vecj
      @test result == resultj
    
      result = mata.'*vec
      resultj = mataj.'*vecj
      @test result == resultj
    
      result = mata'*vec
      resultj = mataj'*vecj
      @test result == resultj
    end
    @testset "binary matrix operations" begin
      result = mata + matb
      assemble(result)
      resultj = mataj + matbj
      @test result == resultj
      result = mata - matb
      assemble(result)
      resultj = mataj - matbj
      @test result == resultj
      result = mata * matb
      assemble(result)
      resultj = mataj * matbj
      @test result == resultj
    end
    @testset "matrix operations with numbers" begin
      result = 2*mata
      assemble(result)
      resultj = 2*mataj

      @test result == resultj
      result = mata/2
      assemble(result)
      resultj = mataj/2
      @test result == resultj

      result = 2.\mata
      assemble(result)
      resultj = 2.\mataj
      @test result == resultj

      result  = -mata
      resultj = -mataj
      @test result == resultj
    end 
  end

  @testset "Testing {c}transpose mults" begin
    mat1  = PETSc.Mat(ST,3,3,mtype=PETSc.C.MATSEQAIJ)
    mat2  = PETSc.Mat(ST,3,3,mtype=PETSc.C.MATSEQAIJ)
    mat1j = zeros(ST,3,3)
    mat2j = zeros(ST,3,3)
    vt1 = RC(complex(3., 0.))
    vt2 = RC(complex(5., 5.))
    vt3 = RC(complex(4., 4.))
    vt4 = RC(complex(10., 0.))
    vt5 = RC(complex(1., 0.))
    mat1[1,1] = vt1
    mat1[1,2] = vt2
    mat1[2,1] = conj(vt2)
    mat1[1,3] = vt3
    mat1[3,1] = conj(vt3)
    mat1[2,2] = vt4
    mat1[3,3] = vt5
    mat1j[1,1] = vt1
    mat1j[1,2] = vt2
    mat1j[2,1] = conj(vt2)
    mat1j[1,3] = vt3
    mat1j[3,1] = conj(vt3)
    mat1j[2,2] = vt4
    mat1j[3,3] = vt5
    mat2[1,1] = vt1
    mat2[1,2] = vt2
    mat2[2,1] = vt2
    mat2[1,3] = vt3
    mat2[3,1] = vt3
    mat2[2,2] = vt4
    mat2[3,3] = vt5
    mat2j[1,1] = vt1
    mat2j[1,2] = vt2
    mat2j[2,1] = vt2
    mat2j[1,3] = vt3
    mat2j[3,1] = vt3
    mat2j[2,2] = vt4
    mat2j[3,3] = vt5
    PETSc.assemble(mat1)
    PETSc.assemble(mat2)
    @test ishermitian(mat1)
    @test issym(mat2)

    mat3 = mat1.'*mat2
    mat3j = mat1j.'*mat2j
    assemble(mat3)
    @test mat3 == mat3j

    mat4 = mat1*mat2.'
    mat4j = mat1j*mat2j.'
    assemble(mat4)
    @test mat4 == mat4j
  end

end
