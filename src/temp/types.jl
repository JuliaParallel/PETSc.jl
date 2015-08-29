# Simple Petsc datatypes and constants

typealias PetscDataType Cint # enum in petscsys.h 
const PETSC_INT = convert(Cint, 0)
const PETSC_DOUBLE = convert(Cint, 1)
const PETSC_COMPLEX = convert(Cint, 2) # actually PetscScalar
const PETSC_LONG = convert(Cint, 3)
const PETSC_SHORT = convert(Cint, 4)
const PETSC_FLOAT = convert(Cint, 5)
const PETSC_CHAR = convert(Cint, 6)
const PETSC_BIT_LOGICAL = convert(Cint, 7)
const PETSC_ENUM = convert(Cint, 8)
const PETSC_BOOL = convert(Cint, 9)
const PETSC___FLOAT128 = convert(Cint, 10)
const PETSC_OBJECT = convert(Cint, 11) 
const PETSC_FUNCTION = convert(Cint, 12)

function petsc_sizeof(t::PetscDataType)
    s = Array(Csize_t, 1)
    chk(ccall((:PetscDataTypeGetSize,petsc), PetscErrorCode,
              (PetscDataType, Ptr{Csize_t}), t, s))
    s[1]
end
macro petsc_type(name, T)
    DT = symbol(string("PETSC_", uppercase(string(name)[6:end])))
    sz = convert(Int, petsc_sizeof(eval(DT))) * 8
    quote
        typealias $(esc(name)) $(symbol(string(T, sz)))
    end
end

export PetscInt, PetscScalar, PetscReal

@petsc_type PetscInt Int
@petsc_type PetscLong Int
@petsc_type PetscShort Int
@petsc_type PetscEnum Int
@petsc_type PetscBool Int
@petsc_type PetscDouble Float
@petsc_type PetscFloat Float

# Not sure how to get sizeof(PetscReal).  For now, assume double precision.
const scalartype = try
    petsc_sizeof(PETSC_COMPLEX) == 2*sizeof(Float64) ? 
    Complex128 : Complex64
catch
    Float64 # TODO: how to detect single precision?
end
if scalartype == Complex128
    typealias PetscScalar Complex128
    typealias PetscReal Float64
elseif scalartype == Complex64
    typealias PetscScalar Complex64
    typealias PetscReal Float64
else # scalartype == Float64
    typealias PetscScalar Float64
    typealias PetscReal Float64
end

typealias InsertMode PetscEnum # enum from petscsys.h
const DEFAULT_INSERTMODE = convert(PetscEnum,0) # only used internally
const INSERT_VALUES = convert(PetscEnum,1)
const ADD_VALUES = convert(PetscEnum,2)
const MAX_VALUES = convert(PetscEnum,3)
insertmode(x) = x.insertmode==DEFAULT_INSERTMODE ? INSERT_VALUES : x.insertmode

typealias NormType PetscEnum # enums from petscvec.h
const NORM_1 = convert(PetscEnum, 0)
const NORM_2 = convert(PetscEnum, 1)
const NORM_FROBENIUS = convert(PetscEnum, 2)
const NORM_INFINITY = convert(PetscEnum, 3)
const NORM_1_AND_2 = convert(PetscEnum, 4)

typealias MatOption PetscEnum # enum from petscmat.h
const MAT_OPTION_MIN = convert(PetscEnum, -8)
const MAT_NEW_NONZERO_LOCATION_ERR = convert(PetscEnum, -7)
const MAT_NO_OFF_PROC_ZERO_ROWS = convert(PetscEnum, -6)
const MAT_NO_OFF_PROC_ENTRIES = convert(PetscEnum, -5)
const MAT_UNUSED_NONZERO_LOCATION_ERR = convert(PetscEnum, -4)
const MAT_NEW_NONZERO_ALLOCATION_ERR = convert(PetscEnum, -3)
const MAT_ROW_ORIENTED = convert(PetscEnum, -2)
const MAT_NEW_NONZERO_LOCATIONS = convert(PetscEnum, -1)
const MAT_SYMMETRIC = convert(PetscEnum, 1)
const MAT_STRUCTURALLY_SYMMETRIC = convert(PetscEnum, 2)
const MAT_NEW_DIAGONALS = convert(PetscEnum, 3)
const MAT_IGNORE_OFF_PROC_ENTRIES = convert(PetscEnum, 4)
const MAT_USE_HASH_TABLE = convert(PetscEnum, 5)
const MAT_KEEP_NONZERO_PATTERN = convert(PetscEnum, 6)
const MAT_IGNORE_ZERO_ENTRIES = convert(PetscEnum, 7)
const MAT_USE_INODES = convert(PetscEnum, 8)
const MAT_HERMITIAN = convert(PetscEnum, 9)
const MAT_SYMMETRY_ETERNAL = convert(PetscEnum, 10)
const MAT_CHECK_COMPRESSED_ROW = convert(PetscEnum, 11)
const MAT_IGNORE_LOWER_TRIANGULAR = convert(PetscEnum, 12)
const MAT_ERROR_LOWER_TRIANGULAR = convert(PetscEnum, 13)
const MAT_GETROW_UPPERTRIANGULAR = convert(PetscEnum, 14)
const MAT_SPD = convert(PetscEnum, 15)
const MAT_OPTION_MAX = convert(PetscEnum, 16)

const PETSC_DECIDE = convert(Cint, -1) # from petscsys.h

typealias MatDuplicateOption PetscEnum # enum from petscmat.h
const MAT_DO_NOT_COPY_VALUES = convert(PetscEnum, 0)
const MAT_COPY_VALUES = convert(PetscEnum, 1)

immutable MatInfo # mirror of struct in petscmat.h
    block_size::Float64
    nz_allocated::Float64
    nz_used::Float64
    nz_unneeded::Float64
    memory::Float64
    assemblies::Float64
    mallocs::Float64
    fill_ratio_given::Float64
    fill_ratio_needed::Float64
    factor_mallocs::Float64
end
typealias MatInfoType PetscEnum # from petscmat.h 
const MAT_LOCAL=convert(PetscEnum, 1)
const MAT_GLOBAL_MAX=convert(PetscEnum, 2)
const MAT_GLOBAL_SUM=convert(PetscEnum, 3)

typealias MatAssemblyType PetscEnum # enum in petscmat.h
const MAT_FLUSH_ASSEMBLY = convert(PetscEnum,1)
const MAT_FINAL_ASSEMBLY = convert(PetscEnum,0)

typealias MatReuse PetscEnum # from petscmat.h
const MAT_INITIAL_MATRIX = convert(PetscEnum, 0)
const MAT_REUSE_MATRIX = convert(PetscEnum, 1)
const MAT_IGNORE_MATRIX = convert(PetscEnum, 2)
