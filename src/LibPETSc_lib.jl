"""
    PetscLibType{PetscScalar, PetscInt}(petsc_library)

A container for specific PETSc libraries.

All other containers for PETSc objects should be typed on this to ensure that
dispatch is correct.
"""
mutable struct PetscLibType{PetscScalar, PetscInt, LibType}
    petsc_library::LibType
    age::Int
end
function PetscLibType{ST, IT}(petsc_library) where {ST, IT}
    LT = typeof(petsc_library)
    return PetscLibType{ST, IT, LT}(petsc_library, 0)
end
const UnionPetscLibType = Union{PetscLibType, Type{<:PetscLibType}}

"""
    getlib(; PetscScalar = Float64, PetscInt = Int64)

Return the `PetscLibType` with the associated parameters
"""
function getlib(; PetscScalar = Float64, PetscInt = Int64)
    return PetscLibType{PetscScalar, PetscInt}()
end

function Base.getproperty(petsclib::UnionPetscLibType, name::Symbol)
    if name == :PetscScalar
        return scalartype(petsclib)
    elseif name == :PetscReal
        return realtype(petsclib)
    elseif name == :PetscInt
        return inttype(petsclib)
    else
        return getfield(petsclib, name)
    end
end

"""
    scalartype(petsclib::PetscLibType)

return the scalar type for the associated `petsclib`
"""
scalartype(::PetscLibType{ST}) where {ST} = ST
scalartype(::Type{PetscLib}) where {PetscLib <: PetscLibType{ST}} where {ST} =
    ST

"""
    realtype(petsclib::PetscLibType)

return the real type for the associated `petsclib`
"""
realtype(::PetscLibType{ST}) where {ST} = real(ST)
realtype(::Type{PetscLib}) where {PetscLib <: PetscLibType{ST}} where {ST} =
    real(ST)

"""
    inttype(petsclib::PetscLibType)

return the int type for the associated `petsclib`
"""
inttype(::PetscLibType{ST, IT}) where {ST, IT} = IT
inttype(
    ::Type{PetscLib},
) where {PetscLib <: PetscLibType{ST, IT}} where {ST, IT} = IT

function initialize(libhdl::Ptr{Cvoid})
    PetscInitializeNoArguments_ptr = dlsym(libhdl, :PetscInitializeNoArguments)
    @chk ccall(PetscInitializeNoArguments_ptr, PetscErrorCode, ())
end

const petsclibs = map(libs) do (lib, PetscScalar, PetscInt)
    return PetscLibType{PetscScalar, PetscInt}(lib[1])
end

macro for_petsc(expr)
    quote
        for petsclib in petsclibs
            # String for the library
            petsc_library = petsclib.petsc_library

            # types we dispatch on
            PetscLib = typeof(petsclib)
            UnionPetscLib = Union{PetscLib, Type{PetscLib}}

            PetscScalar = scalartype(petsclib)
            PetscReal = realtype(petsclib)
            PetscInt = inttype(petsclib)
            PetscComplex = complex(PetscReal)

            @eval esc($expr)
        end
    end
end

@for_petsc begin
    getlib(::Type{$PetscLib}) = $petsclib
    PetscLibType{$PetscScalar, $PetscInt}() = $petsclib
end
