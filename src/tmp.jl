type Vec
  a::Ptr{Void}
end

function (./)(a::Number, x::Vec)

    chk(C.VecReciprocal(x.p))
    if a != 1.0
        scale!(x, a)
    end
    x
end


