# Convenience overloads for PETSc Set*Type functions.
# Each accepts AbstractString and converts to the Ptr{Cchar} the C API expects.
# GC.@preserve keeps the String alive across the ccall inside the LibPETSc wrapper.

function LibPETSc.MatSetType(petsclib, mat, type::AbstractString)
    s = String(type)
    GC.@preserve s LibPETSc.MatSetType(petsclib, mat, Base.unsafe_convert(Ptr{Cchar}, s))
end

function LibPETSc.VecSetType(petsclib, vec, type::AbstractString)
    s = String(type)
    GC.@preserve s LibPETSc.VecSetType(petsclib, vec, Base.unsafe_convert(Ptr{Cchar}, s))
end

function LibPETSc.KSPSetType(petsclib, ksp, type::AbstractString)
    s = String(type)
    GC.@preserve s LibPETSc.KSPSetType(petsclib, ksp, Base.unsafe_convert(Ptr{Cchar}, s))
end

function LibPETSc.SNESSetType(petsclib, snes, type::AbstractString)
    s = String(type)
    GC.@preserve s LibPETSc.SNESSetType(petsclib, snes, Base.unsafe_convert(Ptr{Cchar}, s))
end

function LibPETSc.DMSetType(petsclib, dm, type::AbstractString)
    s = String(type)
    GC.@preserve s LibPETSc.DMSetType(petsclib, dm, Base.unsafe_convert(Ptr{Cchar}, s))
end

# DMSetVecType and DMSetMatType accept AbstractString directly (VecType/MatType = Cstring).
