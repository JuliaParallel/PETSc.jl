# create Vec
@testset "Vec{$ST}" begin
  vtype = PETSc.C.VECMPI
  vec = PETSc.Vec(ST, vtype)
  resize!(vec, 4)
  @test_throws ArgumentError resize!(vec)
  len_ret = length(vec)

  @test length(vec) == 4
  @test size(vec) == (4,)
  @test lengthlocal(vec) == 4
  @test sizelocal(vec) == (4,)
  @test PETSc.gettype(vec) == PETSc.C.VECMPI


  vt = complex(2.,2)  # use vt to hold temporary values
  vec[1] = RC(vt)
  val_ret = vec[1]
  @test vec[1] == RC(vt)

  vec2 = similar(vec,ST)
  PETSc.AssemblyBegin(vec2)
  PETSc.AssemblyEnd(vec2)
  @test isassembled(vec2)
  val2_ret = vec2[1]

  @test val2_ret != val_ret

  if gettype(vec2) == PETSc.C.VECSEQ
    lv2 = localpart(vec2)
    @test lv2 == vec2
  end

  vec_tmp = Vec([1., 2, 3])
  @test PETSc.isfinalized(vec_tmp) == false
  PETSc.PetscDestroy(vec_tmp)
  @test PETSc.isfinalized(vec_tmp) == true

  vec3 = similar(vec, ST, 5)
  @test length(vec3) == 5

  vec4 = copy(vec)
  @test vec4 ≈ vec

  idx = [1,3, 4]
  vt = RC(complex(2.,2))
  vec4[idx] = vt
  vals_ret = vec4[idx]
  @test vals_ret == fill(vt,length(idx))

  vt = RC(complex(3.,3))
  fill!(vec4, vt)

  @test vec4 ≈ fill(vt,length(vec4))

  vt = RC(complex( 4.,4))
  vec4[1:2] = vt

  @test vec4[1:2] == [vt, vt]

  vals = [RC(complex(1,1.)), RC(complex(3.,3)), RC(complex(4., 3))]
  vec4[idx] = vals

  @test vec4[idx] == vals


  vec5 = Vec(Float64, 4)
  varr = LocalArray(vec5)
  @test length(vec5) == 4
  vec5j = [1., 2, 3, 4]
  for i=1:length(vec5)  varr[i] = vec5j[i] end
  
  @test varr[1] == vec5j[1]
  @test varr == vec5j

  LocalArrayRestore(varr)

  @test vec5 == vec5j

  varr = LocalArrayRead(vec5)
  for i=1:length(vec5) @test varr[i] ==  vec5[i] end
  LocalArrayRestore(varr)
  

  # test mlocal constructor
  vec5 = Vec(ST, mlocal=3)
  @test length(vec5) == 3

  @testset "testing logical indexing" begin
      logicals = Array(Bool, length(vec4))
      for i=eachindex(logicals)
        logicals[i] = false
      end
      logicals[2] = true

      vt = RC(complex(5,5.))
      vec4[logicals] = vt

      @test vec4[2] ≈ vt
      @test vec4[1] != vt

      vt = RC(complex(rand(), rand()))
      vals = [vt]
      vec4[logicals] = vals
      @test vec4[2] ≈ vals[1]
      @test vec4[1] != vals[1]

      # reset vec4
      vec4_j = zeros(ST, length(vec4))
      for i=1:length(vec4)
        vec4[i] = RC(complex(Float64(-i), Float64(-i)))
        vec4_j[i] = RC(complex(Float64(-i), Float64(-i)))
      end
  end
  @testset "testing math functions" begin

    @testset "testin chop" begin
       jvec = RC([complex(1.0, 1.0), complex(2.0, 2.0), complex(3.0, 3.0)])
       pvec = Vec(jvec)
       chop!(pvec, RT(1.5))
       jvec[1] = 0.0
       @test pvec ≈ jvec
    end 

    vec4_j = zeros(ST, length(vec4))
    for i=1:length(vec4)
      vec4[i] = RC(complex(Float64(-i), Float64(-i)))
      vec4_j[i] = RC(complex(Float64(-i), Float64(-i)))
    end
    @testset "testing abs" begin
      vec4_j = abs(vec4_j)
      absv4  = abs(vec4)
      abs!(vec4)
      if VERSION >= v"0.5.0-dev+0"
          @test real(vec4) ≈ vec4_j
          @test real(absv4) ≈ vec4_j
          @test imag(vec4) ≈ zeros(vec4_j)
          @test imag(absv4) ≈ zeros(vec4_j)
      else
          @test vec4 == vec4_j
          @test absv4 == vec4_j
      end
    end
    @testset "testing exp" begin
      vec4_j = exp(vec4_j)
      exp!(vec4)
      @test vec4 ≈ vec4_j
    end
    @testset "testing log" begin
      vec4_j = log(vec4_j)
      log!(vec4)
      @test vec4 ≈ vec4_j
    end
    onevec = PETSc.Vec(ST, vtype)
    resize!(onevec, 4)
    PETSc.AssemblyBegin(onevec)
    PETSc.AssemblyEnd(onevec)
    for i=1:length(onevec)
        onevec[i] = one(ST)
    end

    @testset "testing norm" begin
      @test_throws ArgumentError norm(onevec,3)
      @test norm(onevec,Inf) == 1
      normvec = copy(onevec)
      PETSc.normalize!(normvec)
      @test norm(normvec,2) == one(ST)
    end
    if ST <: Real
      @testset "testing max and min" begin
        maxvec = copy(onevec)
        maxvec[1] = ST(2)
        @test maximum(maxvec) == 2
        @test findmax(maxvec) == (2.0,1)
        minvec = copy(onevec)
        minvec[1] = ST(0)
        @test minimum(minvec) == 0
        @test findmin(minvec) == (0.0,1)
      end
    end

    @testset "testing pointwise max, min, /" begin
      div1vec = 2*copy(onevec)
      div2vec = 4*copy(onevec)
      @test max(div1vec,div2vec) == div2vec
      @test min(div1vec,div2vec) == div1vec
      @test div1vec .* div2vec == 8*onevec
      @test div2vec ./ div1vec == div1vec
    end
    @testset "testing scale! and negation" begin
      scalevec = scale!(copy(onevec),2)
      @test scalevec == fill(2,length(onevec))
      minusvec = -onevec
      @test minusvec == -onevec
    end

    @testset "testing sum, +, -, *, and /" begin
      @test sum(onevec) == length(onevec)
      multvec = copy(onevec)
      multvec = multvec * 2 * 3 * 4
      @test multvec == 24*onevec
      multvec = copy(onevec)
      multvec = 2 .* multvec
      @test multvec == 2*onevec
      divvec = copy(onevec)
      divvec = divvec * 2 * 3
      divvec = divvec ./ 2
      @test divvec == 3*onevec
      divvec = 3 .\ divvec
      @test divvec == onevec

      divvec = 2*copy(onevec)
      divvec = 2 ./ divvec
      @test divvec == onevec
      addvec = copy(onevec)
      addvec = addvec + 2
      addvec = addvec - 2
      @test addvec == onevec
      addvec = copy(onevec)
      addvec = 2 - addvec
      addvec = 2 + addvec
      @test addvec == 3*onevec
    end
  end

  @testset "testing dot product" begin
    val = dot(vec4, vec)
    val_j = dot(vec4, vec)
    @test val == val_j
  end
  # make copies of vecs 1 2 4

  @testset "testing level 1 Blas" begin

    vecj = zeros(ST, length(vec))
    vec2j = zeros(ST, length(vec))
    vec4j = zeros(ST, length(vec))

    for i=1:length(vec)
      vecj[i] = vec[i]
      vec2j[i] = vec2[i]
      vec4j[i] = vec4[i]
    end

    @testset "testing axpy" begin
      vt = RC(complex(2.,2))
      axpy!(vt, vec, vec2)
      vec2j = vt*vecj + vec2j
      @test vec2j == vec2

      @testset "testing 4 argument axpy" begin
        axpy!(vt, vec, vec2, vec4)
        vec4j = vt*vecj + vec2j
        @test vec2j == vec2
      end

      @testset "testing aypx" begin
        aypx!(vec, vt, vec2)
        vec2j = vt*vec2j + vec
        @test vec2j == vec2
      end

      vt2 = RC(complex(3.,3))
      vt3 = RC(complex(4.,4))
      @testset "testing axpby" begin
        axpby!(vt, vec, vt2, vec2)
        vec2j = vt*vecj + vt2*vec2j
        @test vec2j == vec2

        axpbypcz!(vt, vec, vt2, vec2, vt3, vec4)
        vec4j = vt*vecj + vt2*vec2j + vt3*vec4j
        @test vec4j == vec4
      end

      vecs = Array(typeof(vec), 2)
      vecs[1] = vec
      vecs[2] = vec2
      alphas = [vt2, vt3]
      axpy!(vec4, alphas, vecs)
      vec4j = vec4j + vt2*vecj + vt3*vec2j
      @test vec4j == vec4
    end
    @testset "testing .*, ./, .^" begin
      vec5 = Vec(ST, 3, vtype=PETSc.C.VECMPI)
      vec6 = similar(vec5)
      vec5j = zeros(ST, 3)
      vec6j = zeros(ST, 3)

      for i=1:3
        i_float = Float64(i)

        vec5[i] = RC(complex(i_float, i_float))
        vec6[i] = RC(complex(i_float+3, i_float+3))
        vec5j[i] = RC(complex(i_float, i_float))
        vec6j[i] = RC(complex(i_float +3, i_float+3))
      end

      vec7 = vec5.*vec6
      vec7j = vec5j.*vec6j
      @test vec7 ≈ vec7j

      vec8 = vec5./vec6
      vec8j = vec5j./vec6j
      @test vec8 ≈ vec8j

      vec9 = vec5.^3
      vec9j = vec5j.^3
      @test vec9 ≈ vec9j

      vec10 = vec5 + vec6
      vec10j = vec5j + vec6j
      @test vec10 ≈ vec10j

      vec11 = vec5 - vec6
      vec11j = vec5j - vec6j
      @test vec11 ≈ vec11j
    end
    @testset "test unconjugated dot product" begin
      x = Vec(ST, 2)
      y = Vec(ST, 2)
      copy!(y, [1, 1])
      if ST <: Complex
          copy!(x, [1, im])
          @test (x'*y)[1] == 1-im
          @test (x.'*y)[1] == 1+im
      else
          copy!(x, [2, 3])
          @test (x'*y)[1] == 5
          @test (x.'*y)[1] == 5
      end
    end
  end
  let x = rand(ST, 7)
    @test Vec(x) == x
  end

  @testset "map" begin
    x = rand(3)
    y = Vec(x)
    map!(sin, x)
    map!(sin, y)
    @test x ≈ y
    println("finished first test")
    x2 = map(sin, x)
    y2 = map(sin, y)
    @test x2 ≈ y2
    println("finished second test")
#=
    function myfunc(a, b)
      return a + b
    end

    x3 = copy(x2)
    y3 = copy(y2)
    map!(myfunc, x3, x2, x)
    map!(myfunc, y3, y2, y)
    @test x3 ≈ y3
    println("finished third test")
=#
  end
end
