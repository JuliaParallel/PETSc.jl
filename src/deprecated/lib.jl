using Libdl
"""
    PetscLibType{PetscScalar, PetscInt}(petsc_library)

A container for specific PETSc libraries.

All other containers for PETSc objects should be typed on this to ensure that
dispatch is correct.
"""
struct PetscLibType{PetscScalar, PetscInt, LibType}
    petsc_library::LibType
end
function PetscLibType{ST, IT}(petsc_library) where {ST, IT}
    LT = typeof(petsc_library)
    return PetscLibType{ST, IT, LT}(petsc_library)
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

const petsclibs = map(libs) do lib
    libhdl = dlopen(lib...)

    # initialize petsc
    PetscInitializeNoArguments_ptr =
        dlsym(libhdl, :PetscInitializeNoArguments)
    @chk ccall(PetscInitializeNoArguments_ptr, PetscErrorCode, ())

    PETSC_REAL = DataTypeFromString(libhdl, "Real")
    PETSC_SCALAR = DataTypeFromString(libhdl, "Scalar")
    PETSC_INT_SIZE = PetscDataTypeGetSize(libhdl, PETSC_INT)

    PetscReal =
        PETSC_REAL == PETSC_DOUBLE ? Cdouble :
        PETSC_REAL == PETSC_FLOAT ? Cfloat :
        error("PETSC_REAL = $PETSC_REAL not supported.")

    PetscScalar =
        PETSC_SCALAR == PETSC_REAL ? PetscReal :
        PETSC_SCALAR == PETSC_COMPLEX ? Complex{PetscReal} :
        error("PETSC_SCALAR = $PETSC_SCALAR not supported.")

    PetscInt =
        PETSC_INT_SIZE == 4 ? Int32 :
        PETSC_INT_SIZE == 8 ? Int64 :
        error("PETSC_INT_SIZE = $PETSC_INT_SIZE not supported.")

    PetscFinalize_ptr = dlsym(libhdl, :PetscFinalize)
    @chk ccall(PetscFinalize_ptr, PetscErrorCode, ())

    # TODO: PetscBLASInt, PetscMPIInt ?
    return PetscLibType{PetscScalar, PetscInt}(lib[1])
end

# New macro is really to track the update
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

            @eval esc($expr)
        end
    end
end

# TODO: Remove this macro once we have switch over to @for_petsc macro
macro for_libpetsc(expr)
    quote
        for petsclib in petsclibs
            # String for the library
            libpetsc = petsclib.petsc_library

            # types we dispatch on
            PetscLib = typeof(petsclib)
            UnionPetscLib = Union{PetscLib, Type{PetscLib}}

            PetscScalar = scalartype(petsclib)
            PetscReal = realtype(petsclib)
            PetscInt = inttype(petsclib)

            @eval esc($expr)
        end
    end
end

# TODO: These should be removed!
const scalar_types = map(x -> scalartype(x), petsclibs)
@assert length(scalar_types) == length(unique(scalar_types))
@for_libpetsc inttype(::Type{$PetscScalar}) = $PetscInt
