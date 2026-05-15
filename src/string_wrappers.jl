# Convenience overloads for PETSc Set*Type functions.
# Each accepts AbstractString and converts to the Ptr{Cchar} the C API expects.
# GC.@preserve keeps the String alive across the ccall inside the LibPETSc wrapper.

"""
    MatSetType(petsclib, mat, type::AbstractString)

Set the matrix type. Accepts any `AbstractString`.

# External Links
$(_doc_external("Mat/MatSetType"))
"""
function LibPETSc.MatSetType(petsclib, mat, type::AbstractString)
    s = String(type)
    GC.@preserve s LibPETSc.MatSetType(petsclib, mat, Base.unsafe_convert(Ptr{Cchar}, s))
    return nothing
end

"""
    VecSetType(petsclib, vec, type::AbstractString)

Set the vector type. Accepts any `AbstractString`.

# External Links
$(_doc_external("Vec/VecSetType"))
"""
function LibPETSc.VecSetType(petsclib, vec, type::AbstractString)
    s = String(type)
    GC.@preserve s LibPETSc.VecSetType(petsclib, vec, Base.unsafe_convert(Ptr{Cchar}, s))
    return nothing
end

"""
    KSPSetType(petsclib, ksp, type::AbstractString)

Set the KSP solver type. Accepts any `AbstractString`.

# External Links
$(_doc_external("KSP/KSPSetType"))
"""
function LibPETSc.KSPSetType(petsclib, ksp, type::AbstractString)
    s = String(type)
    GC.@preserve s LibPETSc.KSPSetType(petsclib, ksp, Base.unsafe_convert(Ptr{Cchar}, s))
    return nothing
end

"""
    SNESSetType(petsclib, snes, type::AbstractString)

Set the SNES nonlinear solver type. Accepts any `AbstractString`.

# External Links
$(_doc_external("SNES/SNESSetType"))
"""
function LibPETSc.SNESSetType(petsclib, snes, type::AbstractString)
    s = String(type)
    GC.@preserve s LibPETSc.SNESSetType(petsclib, snes, Base.unsafe_convert(Ptr{Cchar}, s))
    return nothing
end

"""
    DMSetType(petsclib, dm, type::AbstractString)

Set the DM type. Accepts any `AbstractString`.

# External Links
$(_doc_external("DM/DMSetType"))
"""
function LibPETSc.DMSetType(petsclib, dm, type::AbstractString)
    s = String(type)
    GC.@preserve s LibPETSc.DMSetType(petsclib, dm, Base.unsafe_convert(Ptr{Cchar}, s))
    return nothing
end

# DMSetVecType and DMSetMatType accept AbstractString directly (VecType/MatType = Cstring).
