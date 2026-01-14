# Convenience wrappers for PETSc SetType functions that accept Julia strings
# instead of C string pointers
#
# These wrappers are defined in the parent PETSc module and delegate to
# LibPETSc functions with automatic string-to-pointer conversion.

"""
    MatSetType(petsclib, mat, type::String)

Convenience wrapper for setting matrix type using a Julia string.

# Example
```julia
mat = LibPETSc.MatCreate(petsclib, LibPETSc.PETSC_COMM_SELF)
LibPETSc.MatSetType(petsclib, mat, "seqaij")
```
"""
function LibPETSc.MatSetType(petsclib::LibPETSc.PetscLibType, mat, type::String)
    c_str = Vector{UInt8}(type * "\0")
    ptr = Base.unsafe_convert(Ptr{Int8}, pointer(c_str))
    LibPETSc.MatSetType(petsclib, mat, ptr)
    return nothing
end

"""
    VecSetType(petsclib, vec, type::String)

Convenience wrapper for setting vector type using a Julia string.

# Example
```julia
vec = LibPETSc.VecCreate(petsclib, LibPETSc.PETSC_COMM_SELF)
LibPETSc.VecSetType(petsclib, vec, "seq")
```
"""
function LibPETSc.VecSetType(petsclib::LibPETSc.PetscLibType, vec, type::String)
    c_str = Vector{UInt8}(type * "\0")
    ptr = Base.unsafe_convert(Ptr{Int8}, pointer(c_str))
    LibPETSc.VecSetType(petsclib, vec, ptr)
    return nothing
end

"""
    KSPSetType(petsclib, ksp, type::String)

Convenience wrapper for setting KSP solver type using a Julia string.

# Example
```julia
ksp = LibPETSc.KSPCreate(petsclib, LibPETSc.PETSC_COMM_SELF)
LibPETSc.KSPSetType(petsclib, ksp, "gmres")
```
"""
function LibPETSc.KSPSetType(petsclib::LibPETSc.PetscLibType, ksp, type::String)
    c_str = Vector{UInt8}(type * "\0")
    ptr = Base.unsafe_convert(Ptr{Int8}, pointer(c_str))
    LibPETSc.KSPSetType(petsclib, ksp, ptr)
    return nothing
end

"""
    SNESSetType(petsclib, snes, type::String)

Convenience wrapper for setting SNES solver type using a Julia string.

# Example
```julia
snes = LibPETSc.SNESCreate(petsclib, LibPETSc.PETSC_COMM_SELF)
LibPETSc.SNESSetType(petsclib, snes, "newtonls")
```
"""
function LibPETSc.SNESSetType(petsclib::LibPETSc.PetscLibType, snes, type::String)
    c_str = Vector{UInt8}(type * "\0")
    ptr = Base.unsafe_convert(Ptr{Int8}, pointer(c_str))
    LibPETSc.SNESSetType(petsclib, snes, ptr)
    return nothing
end

"""
    DMSetType(petsclib, dm, type::String)

Convenience wrapper for setting DM type using a Julia string.

# Example
```julia
dm = LibPETSc.DMCreate(petsclib, LibPETSc.PETSC_COMM_SELF)
LibPETSc.DMSetType(petsclib, dm, "da")
```
"""
function LibPETSc.DMSetType(petsclib::LibPETSc.PetscLibType, dm, type::String)
    c_str = Vector{UInt8}(type * "\0")
    ptr = Base.unsafe_convert(Ptr{Int8}, pointer(c_str))
    LibPETSc.DMSetType(petsclib, dm, ptr)
    return nothing
end
