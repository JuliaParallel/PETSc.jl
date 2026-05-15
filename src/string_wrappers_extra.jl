"""
    TSSetType(petsclib, ts, type::AbstractString)

Set the TS time-stepping type. Accepts any `AbstractString`.

# External Links
$(_doc_external("TS/TSSetType"))
"""
function LibPETSc.TSSetType(petsclib, ts, type::AbstractString)
    s = String(type)
    GC.@preserve s LibPETSc.TSSetType(petsclib, ts, Base.unsafe_convert(Ptr{Cchar}, s))
    return nothing
end

"""
    TaoSetType(petsclib, tao, type::AbstractString)

Set the Tao optimization solver type. Accepts any `AbstractString`.

# External Links
$(_doc_external("Tao/TaoSetType"))
"""
function LibPETSc.TaoSetType(petsclib, tao, type::AbstractString)
    s = String(type)
    GC.@preserve s LibPETSc.TaoSetType(petsclib, tao, Base.unsafe_convert(Ptr{Cchar}, s))
    return nothing
end
