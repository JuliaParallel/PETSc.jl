# Test index sets and vector gather/scatter

facts("\n --- Testing IS Functions ---") do
  let i = IS(ST, 10:-2:1), j = sort(i)
    @fact issorted(i) --> false
    @fact issorted(j) --> true
    @fact i == j --> true
    @fact i ∪ j == i --> true
    @fact Vector{Int}(i) --> [10:-2:1;]
    @fact Vector{Int}(j) --> [2:2:10;]
    @fact length(i) --> 5
    @fact lengthlocal(i) --> 5
    @fact minimum(i) --> 2
    @fact maximum(i) --> 10
    @fact Set(setdiff(i, IS(ST, 1:4))) --> Set(6:2:10)
    @fact Set(union(i, IS(ST, [1,2,3,4]))) --> Set(2:2:10) ∪ Set(1:4)
  end

  let x = Vec(ST[1,17,24,2]), y = Vec(ST, 2)
    @fact scatter!(x, 2:3, y, 1:2) --> [17,24]
  end
end
