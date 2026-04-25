function LibPETSc.TSSetType(petsclib, ts, type::AbstractString)
    s = String(type)
    GC.@preserve s LibPETSc.TSSetType(petsclib, ts, Base.unsafe_convert(Ptr{Cchar}, s))
end

function LibPETSc.TaoSetType(petsclib, tao, type::AbstractString)
    s = String(type)
    GC.@preserve s LibPETSc.TaoSetType(petsclib, tao, Base.unsafe_convert(Ptr{Cchar}, s))
end
