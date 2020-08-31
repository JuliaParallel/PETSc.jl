# define common PETSc constants
# this excludes configurable constants (e.g. PetscScalar) which are set in lib.jl

const PetscErrorCode      = Cint

struct PetscError <: Exception
    code::PetscErrorCode
end

macro chk(expr)
    :((errcode = $(esc(expr))) == 0 || throw(PetscError(errcode)))
end

const PETSC_DEFAULT = -2


@enum PetscBool PETSC_FALSE PETSC_TRUE

@enum PetscDataType begin
    PETSC_DATATYPE_UNKNOWN  = 0
    PETSC_DOUBLE            = 1
    PETSC_COMPLEX           = 2
    PETSC_LONG              = 3
    PETSC_SHORT             = 4
    PETSC_FLOAT             = 5
    PETSC_CHAR              = 6
    PETSC_BIT_LOGICAL       = 7
    PETSC_ENUM              = 8
    PETSC_BOOL              = 9
    PETSC_FLOAT128          = 10
    PETSC_OBJECT            = 11
    PETSC_FUNCTION          = 12
    PETSC_STRING            = 13
    PETSC___FP16            = 14
    PETSC_STRUCT            = 15
    PETSC_INT               = 16
    PETSC_INT64             = 17
end

const Petsc64bitInt       = Int64
const PetscLogDouble      = Cdouble

@enum InsertMode NOT_SET_VALUES INSERT_VALUES ADD_VALUES MAX_VALUES MIN_VALUES INSERT_ALL_VALUES ADD_ALL_VALUES INSERT_BC_VALUES ADD_BC_VALUES



@enum NormType begin
    NORM_1 = 0
    NORM_2 = 1
    NORM_FROBENIUS = 2
    NORM_INFINITY = 3
    NORM_1_AND_2 = 4
end

@enum MatAssemblyType MAT_FLUSH_ASSEMBLY=1 MAT_FINAL_ASSEMBLY=0

@enum MatFactorType begin
    MAT_FACTOR_NONE     = 0
    MAT_FACTOR_LU       = 1
    MAT_FACTOR_CHOLESKY = 2
    MAT_FACTOR_ILU      = 3
    MAT_FACTOR_ICC      = 4
    MAT_FACTOR_ILUDT    = 5
end
