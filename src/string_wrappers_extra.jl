"""
    TSSetType(petsclib, ts, type::String)

Convenience wrapper for setting TS (time-stepping) type using a Julia string.

# Example
```julia
ts = LibPETSc.TSCreate(petsclib, LibPETSc.PETSC_COMM_SELF)
LibPETSc.TSSetType(petsclib, ts, "bdf")
```
"""
function LibPETSc.TSSetType(petsclib::LibPETSc.PetscLibType, ts, type::String)
    c_str = Vector{UInt8}(type * "\0")
    ptr = Base.unsafe_convert(Ptr{Int8}, pointer(c_str))
    LibPETSc.TSSetType(petsclib, ts, ptr)
    return nothing
end

"""
    TaoSetType(petsclib, tao, type::String)

Convenience wrapper for setting Tao solver type using a Julia string.

# Example
```julia
tao = LibPETSc.TaoCreate(petsclib)
LibPETSc.TaoSetType(petsclib, tao, "lmvm")
```
"""
function LibPETSc.TaoSetType(petsclib::LibPETSc.PetscLibType, tao, type::String)
    c_str = Vector{UInt8}(type * "\0")
    ptr = Base.unsafe_convert(Ptr{Int8}, pointer(c_str))
    LibPETSc.TaoSetType(petsclib, tao, ptr)
    return nothing
end
