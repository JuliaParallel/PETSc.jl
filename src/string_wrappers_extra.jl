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

"""
    TSRKSetType(petsclib, ts, subtype::String)

Convenience wrapper for selecting a Runge-Kutta subtype on an explicit RK
`TS` using a Julia string (e.g. `"3bs"`, `"5dp"`).

# Example
```julia
ts = LibPETSc.TSCreate(petsclib, LibPETSc.PETSC_COMM_SELF)
LibPETSc.TSSetType(petsclib, ts, "rk")
LibPETSc.TSRKSetType(petsclib, ts, "3bs")
```
"""
function LibPETSc.TSRKSetType(petsclib::LibPETSc.PetscLibType, ts, subtype::String)
    c_str = Vector{UInt8}(subtype * "\0")
    ptr = Base.unsafe_convert(LibPETSc.TSRKType, pointer(c_str))
    LibPETSc.TSRKSetType(petsclib, ts, ptr)
    return nothing
end

"""
    TSRosWSetType(petsclib, ts, subtype::String)

Convenience wrapper for selecting a Rosenbrock-W subtype on a `TS` of type
`"rosw"` using a Julia string (e.g. `"ra34pw2"`, `"rodas3"`).

# Example
```julia
ts = LibPETSc.TSCreate(petsclib, LibPETSc.PETSC_COMM_SELF)
LibPETSc.TSSetType(petsclib, ts, "rosw")
LibPETSc.TSRosWSetType(petsclib, ts, "ra34pw2")
```
"""
function LibPETSc.TSRosWSetType(petsclib::LibPETSc.PetscLibType, ts, subtype::String)
    c_str = Vector{UInt8}(subtype * "\0")
    ptr = Base.unsafe_convert(LibPETSc.TSRosWType, pointer(c_str))
    LibPETSc.TSRosWSetType(petsclib, ts, ptr)
    return nothing
end

"""
    TSARKIMEXSetType(petsclib, ts, subtype::String)

Convenience wrapper for selecting an ARK IMEX subtype on a `TS` of type
`"arkimex"` using a Julia string (e.g. `"2e"`, `"3"`).

# Example
```julia
ts = LibPETSc.TSCreate(petsclib, LibPETSc.PETSC_COMM_SELF)
LibPETSc.TSSetType(petsclib, ts, "arkimex")
LibPETSc.TSARKIMEXSetType(petsclib, ts, "2e")
```
"""
function LibPETSc.TSARKIMEXSetType(petsclib::LibPETSc.PetscLibType, ts, subtype::String)
    c_str = Vector{UInt8}(subtype * "\0")
    ptr = Base.unsafe_convert(LibPETSc.TSARKIMEXType, pointer(c_str))
    LibPETSc.TSARKIMEXSetType(petsclib, ts, ptr)
    return nothing
end
