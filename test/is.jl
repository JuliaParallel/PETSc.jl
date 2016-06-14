# Test index sets and vector gather/scatter

@testset "IS{$ST}" begin
  let i = IS(ST, 10:-2:1), j = sort(i)
    @test !issorted(i)
    @test issorted(j)
    @test i == j
    @test i ∪ j == i
    @test Vector{Int}(i) == [10:-2:1;]
    @test Vector{Int}(j) == [2:2:10;]
    @test length(i) == 5
    @test lengthlocal(i) == 5
    @test minimum(i) == 2
    @test maximum(i) == 10
    @test Set(setdiff(i, IS(ST, 1:4))) == Set(6:2:10)
    @test Set(union(i, IS(ST, [1,2,3,4]))) == Set(2:2:10) ∪ Set(1:4)
  end

  let x = Vec(ST[1,17,24,2]), y = Vec(ST, 2)
    @test scatter!(x, 2:3, y, 1:2) == [17,24]
  end
  
  let x = Vec(ST[1,17,24,2]), y = Vec(ST, 2)
    VS = VecScatter(x, IS(ST, 2:3, comm=comm(x)), y, IS(ST, 1:2, comm=comm(y)))
    @test scatter!(copy(VS),x,y) == [17,24]
  end

  let
    idx = 2:4
    is = IS(Float64, idx)
    bs = 2
    set_blocksize(is, bs)
    @test get_blocksize(is) == bs
  end
end
