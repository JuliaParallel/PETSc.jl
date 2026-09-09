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
