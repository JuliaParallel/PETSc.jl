# define common PETSc constants
# this excludes configurable constants (e.g. PetscScalar) which are set in lib.jl

const PetscErrorCode = Cint

struct PetscError <: Exception
    code::PetscErrorCode
end

macro chk(expr)
    :((errcode = $(esc(expr))) == 0 || throw(PetscError(errcode)))
end

@enum PetscBool::UInt32 begin
    PETSC_FALSE = 0
    PETSC_TRUE = 1
end

@enum PetscDataType::UInt32 begin
    PETSC_DATATYPE_UNKNOWN = 0
    PETSC_DOUBLE = 1
    PETSC_COMPLEX = 2
    PETSC_LONG = 3
    PETSC_SHORT = 4
    PETSC_FLOAT = 5
    PETSC_CHAR = 6
    PETSC_BIT_LOGICAL = 7
    PETSC_ENUM = 8
    PETSC_BOOL = 9
    PETSC___FLOAT128 = 10
    PETSC_OBJECT = 11
    PETSC_FUNCTION = 12
    PETSC_STRING = 13
    PETSC___FP16 = 14
    PETSC_STRUCT = 15
    PETSC_INT = 16
    PETSC_INT64 = 17
end
