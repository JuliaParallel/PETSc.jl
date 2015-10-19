# Automatically generated using Clang.jl wrap_c, version 0.0.0

using Compat

#skipping undefined const PETSC_FUNCTION_NAME = PETSC_FUNCTION_NAME_C
#skipping undefined const PETSC_RESTRICT = PETSC_C_RESTRICT
#skipping undefined const PETSC_STATIC_INLINE = PETSC_C_STATIC_INLINE
#skipping undefined const PETSC_VISIBILITY_PUBLIC = PETSC_DLLIMPORT
#= # Skipping MacroDefinition: PETSC_EXTERN extern PETSC_VISIBILITY_PUBLIC =#
#= # Skipping MacroDefinition: PETSC_INTERN extern PETSC_VISIBILITY_INTERNAL =#
const PETSC_VERSION_RELEASE = 1
const PETSC_VERSION_MAJOR = 3
const PETSC_VERSION_MINOR = 6
const PETSC_VERSION_SUBMINOR = 0
const PETSC_VERSION_PATCH = 0
const PETSC_RELEASE_DATE = symbol("Jun, 9, 2015")
const PETSC_VERSION_DATE = symbol("Jun, 09, 2015")
const PETSC_VERSION_GIT = symbol("v3.6")
const PETSC_VERSION_DATE_GIT = symbol("2015-06-09 16:15:46 -0500")

#= # Skipping MacroDefinition: PETSC_VERSION_ ( MAJOR , MINOR , SUBMINOR ) ( ( PETSC_VERSION_MAJOR == ( MAJOR ) ) && ( PETSC_VERSION_MINOR == ( MINOR ) ) && ( PETSC_VERSION_SUBMINOR == ( SUBMINOR ) ) && ( PETSC_VERSION_RELEASE == 1 ) ) =#
#= # Skipping MacroDefinition: PETSC_VERSION_LT ( MAJOR , MINOR , SUBMINOR ) ( PETSC_VERSION_RELEASE == 1 && ( PETSC_VERSION_MAJOR < ( MAJOR ) || ( PETSC_VERSION_MAJOR == ( MAJOR ) && ( PETSC_VERSION_MINOR < ( MINOR ) || ( PETSC_VERSION_MINOR == ( MINOR ) && ( PETSC_VERSION_SUBMINOR < ( SUBMINOR ) ) ) ) ) ) ) =#
#= # Skipping MacroDefinition: PETSC_VERSION_LE ( MAJOR , MINOR , SUBMINOR ) ( PETSC_VERSION_LT ( MAJOR , MINOR , SUBMINOR ) || PETSC_VERSION_ ( MAJOR , MINOR , SUBMINOR ) ) =#
#= # Skipping MacroDefinition: PETSC_VERSION_GT ( MAJOR , MINOR , SUBMINOR ) ( 0 == PETSC_VERSION_LE ( MAJOR , MINOR , SUBMINOR ) ) =#
#= # Skipping MacroDefinition: PETSC_VERSION_GE ( MAJOR , MINOR , SUBMINOR ) ( 0 == PETSC_VERSION_LT ( MAJOR , MINOR , SUBMINOR ) ) =#
const PETSC_AUTHOR_INFO = symbol("       The PETSc Team\n    petsc-maint@mcs.anl.gov\n http://www.mcs.anl.gov/petsc/\n")
const MPICH_SKIP_MPICXX = 1
const OMPI_SKIP_MPICXX = 1

#= # Skipping MacroDefinition: PetscAttrMPIPointerWithType ( bufno , typeno ) __attribute__ ( ( pointer_with_type_tag ( MPI , bufno , typeno ) ) ) =#
#= # Skipping MacroDefinition: PetscAttrMPITypeTag ( type ) __attribute__ ( ( type_tag_for_datatype ( MPI , type ) ) ) =#
#= # Skipping MacroDefinition: PetscAttrMPITypeTagLayoutCompatible ( type ) __attribute__ ( ( type_tag_for_datatype ( MPI , type , layout_compatible ) ) ) =#
#skipping undefined const MPIU_INT = MPI_INT
#skipping undefined const MPIU_INT64 = MPI_LONG_LONG_INT
#skipping undefined const MPIU_SIZE_T = MPI_UNSIGNED
#= # Skipping MacroDefinition: PetscUnlikely ( cond ) ( cond ) =#
#= # Skipping MacroDefinition: PetscLikely ( cond ) ( cond ) =#
#= # Skipping MacroDefinition: PetscExpPassiveScalar ( a ) PetscExpScalar ( ) =#
#skipping undefined const MPIU_REAL = MPI_FLOAT
#= # Skipping MacroDefinition: PetscSqrtReal ( a ) sqrt ( a ) =#
#= # Skipping MacroDefinition: PetscExpReal ( a ) exp ( a ) =#
#= # Skipping MacroDefinition: PetscLogReal ( a ) log ( a ) =#
#= # Skipping MacroDefinition: PetscLog10Real ( a ) log10 ( a ) =#
#= # Skipping MacroDefinition: PetscSinReal ( a ) sin ( a ) =#
#= # Skipping MacroDefinition: PetscCosReal ( a ) cos ( a ) =#
#= # Skipping MacroDefinition: PetscTanReal ( a ) tan ( a ) =#
#= # Skipping MacroDefinition: PetscAsinReal ( a ) asin ( a ) =#
#= # Skipping MacroDefinition: PetscAcosReal ( a ) acos ( a ) =#
#= # Skipping MacroDefinition: PetscAtanReal ( a ) atan ( a ) =#
#= # Skipping MacroDefinition: PetscAtan2Real ( a , b ) atan2 ( a , b ) =#
#= # Skipping MacroDefinition: PetscSinhReal ( a ) sinh ( a ) =#
#= # Skipping MacroDefinition: PetscCoshReal ( a ) cosh ( a ) =#
#= # Skipping MacroDefinition: PetscTanhReal ( a ) tanh ( a ) =#
#= # Skipping MacroDefinition: PetscPowReal ( a , b ) pow ( a , b ) =#
#= # Skipping MacroDefinition: PetscCeilReal ( a ) ceil ( a ) =#
#= # Skipping MacroDefinition: PetscFloorReal ( a ) floor ( a ) =#
#= # Skipping MacroDefinition: PetscFmodReal ( a , b ) fmod ( a , b ) =#
#= # Skipping MacroDefinition: PetscTGamma ( a ) tgammaf ( a ) =#
#skipping undefined const MPIU_SCALAR = MPIU_REAL
#= # Skipping MacroDefinition: PetscRealPart ( a ) ( a ) =#
#= # Skipping MacroDefinition: PetscImaginaryPart ( a ) ( ( PetscReal ) 0. ) =#
#= # Skipping MacroDefinition: PetscConj ( a ) ( a ) =#
#= # Skipping MacroDefinition: PetscSqrtScalar ( a ) sqrt ( a ) =#
#= # Skipping MacroDefinition: PetscPowScalar ( a , b ) pow ( a , b ) =#
#= # Skipping MacroDefinition: PetscExpScalar ( a ) exp ( a ) =#
#= # Skipping MacroDefinition: PetscLogScalar ( a ) log ( a ) =#
#= # Skipping MacroDefinition: PetscSinScalar ( a ) sin ( a ) =#
#= # Skipping MacroDefinition: PetscCosScalar ( a ) cos ( a ) =#
#= # Skipping MacroDefinition: PetscAsinScalar ( a ) asin ( a ) =#
#= # Skipping MacroDefinition: PetscAcosScalar ( a ) acos ( a ) =#
#= # Skipping MacroDefinition: PetscTanScalar ( a ) tan ( a ) =#
#= # Skipping MacroDefinition: PetscSinhScalar ( a ) sinh ( a ) =#
#= # Skipping MacroDefinition: PetscCoshScalar ( a ) cosh ( a ) =#
#= # Skipping MacroDefinition: PetscTanhScalar ( a ) tanh ( a ) =#
#= # Skipping MacroDefinition: PetscSign ( a ) ( ( ( a ) >= 0 ) ? ( ( a ) == 0 ? 0 : 1 ) : - 1 ) =#
#= # Skipping MacroDefinition: PetscAbs ( a ) ( ( ( a ) >= 0 ) ? ( a ) : - ( a ) ) =#
#= # Skipping MacroDefinition: PetscMin ( a , b ) ( ( ( a ) < ( b ) ) ? ( a ) : ( b ) ) =#
#= # Skipping MacroDefinition: PetscMax ( a , b ) ( ( ( a ) < ( b ) ) ? ( b ) : ( a ) ) =#
#= # Skipping MacroDefinition: PetscClipInterval ( x , a , b ) ( PetscMax ( ( a ) , PetscMin ( ( x ) , ( b ) ) ) ) =#
#= # Skipping MacroDefinition: PetscAbsInt ( a ) ( ( ( a ) < 0 ) ? - ( a ) : ( a ) ) =#
#= # Skipping MacroDefinition: PetscAbsReal ( a ) ( ( ( a ) < 0 ) ? - ( a ) : ( a ) ) =#
#= # Skipping MacroDefinition: PetscSqr ( a ) ( ( a ) * ( a ) ) =#
const PETSC_PI = 3.141592653589793
const PETSC_MAX_INT = 2147483647
const PETSC_MIN_INT = -PETSC_MAX_INT - 1
const PETSC_MAX_REAL = 3.4028234663852886e38
const PETSC_MIN_REAL = -PETSC_MAX_REAL
const PETSC_MACHINE_EPSILON = 1.1920929e-7
const PETSC_SQRT_MACHINE_EPSILON = 0.000345266983
const PETSC_SMALL = 1.0e-5
const PETSC_INFINITY = PETSC_MAX_REAL / 4.0
const PETSC_NINFINITY = -PETSC_INFINITY

# excluding lhs of  typealias PetscReal Cfloat
const PassiveReal = Float32

# excluding lhs of  typealias PetscScalar PetscReal
const PassiveScalar = Float32

#skipping undefined const MPIU_MATSCALAR = MPIU_SCALAR
#skipping undefined const MPIU_2INT = MPI_2INT
const PETSC_NULL = C_NULL
const PETSC_IGNORE = C_NULL
const PETSC_DECIDE = -1
const PETSC_DETERMINE = PETSC_DECIDE
const PETSC_DEFAULT = -2
const PETSC_COMM_SELF = MPI_COMM_SELF

#= # Skipping MacroDefinition: PetscMalloc ( a , b ) ( ( * PetscTrMalloc ) ( ( a ) , __LINE__ , PETSC_FUNCTION_NAME , __FILE__ , ( void * * ) ( b ) ) ) =#
#= # Skipping MacroDefinition: PetscAddrAlign ( a ) ( void * ) ( ( ( ( PETSC_UINTPTR_T ) ( a ) ) + ( PETSC_MEMALIGN - 1 ) ) & ~ ( PETSC_MEMALIGN - 1 ) ) =#
#= # Skipping MacroDefinition: PetscMalloc1 ( m1 , r1 ) PetscMalloc ( ( m1 ) * sizeof ( * * ( r1 ) ) , r1 ) =#
#= # Skipping MacroDefinition: PetscCalloc1 ( m1 , r1 ) ( PetscMalloc1 ( ( m1 ) , r1 ) || PetscMemzero ( * ( r1 ) , ( m1 ) * sizeof ( * * ( r1 ) ) ) ) =#
#= # Skipping MacroDefinition: PetscMalloc2 ( m1 , r1 , m2 , r2 ) ( PetscMalloc1 ( ( m1 ) , ( r1 ) ) || PetscMalloc1 ( ( m2 ) , ( r2 ) ) ) =#
#= # Skipping MacroDefinition: PetscCalloc2 ( m1 , r1 , m2 , r2 ) ( PetscMalloc2 ( ( m1 ) , ( r1 ) , ( m2 ) , ( r2 ) ) || PetscMemzero ( * ( r1 ) , ( m1 ) * sizeof ( * * ( r1 ) ) ) || PetscMemzero ( * ( r2 ) , ( m2 ) * sizeof ( * * ( r2 ) ) ) ) =#
#= # Skipping MacroDefinition: PetscMalloc3 ( m1 , r1 , m2 , r2 , m3 , r3 ) ( PetscMalloc1 ( ( m1 ) , ( r1 ) ) || PetscMalloc1 ( ( m2 ) , ( r2 ) ) || PetscMalloc1 ( ( m3 ) , ( r3 ) ) ) =#
#= # Skipping MacroDefinition: PetscCalloc3 ( m1 , r1 , m2 , r2 , m3 , r3 ) ( PetscMalloc3 ( ( m1 ) , ( r1 ) , ( m2 ) , ( r2 ) , ( m3 ) , ( r3 ) ) || PetscMemzero ( * ( r1 ) , ( m1 ) * sizeof ( * * ( r1 ) ) ) || PetscMemzero ( * ( r2 ) , ( m2 ) * sizeof ( * * ( r2 ) ) ) || PetscMemzero ( * ( r3 ) , ( m3 ) * sizeof ( * * ( r3 ) ) ) ) =#
#= # Skipping MacroDefinition: PetscMalloc4 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 ) ( PetscMalloc1 ( ( m1 ) , ( r1 ) ) || PetscMalloc1 ( ( m2 ) , ( r2 ) ) || PetscMalloc1 ( ( m3 ) , ( r3 ) ) || PetscMalloc1 ( ( m4 ) , ( r4 ) ) ) =#
#= # Skipping MacroDefinition: PetscCalloc4 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 ) ( PetscMalloc4 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 ) || PetscMemzero ( * ( r1 ) , ( m1 ) * sizeof ( * * ( r1 ) ) ) || PetscMemzero ( * ( r2 ) , ( m2 ) * sizeof ( * * ( r2 ) ) ) || PetscMemzero ( * ( r3 ) , ( m3 ) * sizeof ( * * ( r3 ) ) ) || PetscMemzero ( * ( r4 ) , ( m4 ) * sizeof ( * * ( r4 ) ) ) ) =#
#= # Skipping MacroDefinition: PetscMalloc5 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 , m5 , r5 ) ( PetscMalloc1 ( ( m1 ) , ( r1 ) ) || PetscMalloc1 ( ( m2 ) , ( r2 ) ) || PetscMalloc1 ( ( m3 ) , ( r3 ) ) || PetscMalloc1 ( ( m4 ) , ( r4 ) ) || PetscMalloc1 ( ( m5 ) , ( r5 ) ) ) =#
#= # Skipping MacroDefinition: PetscCalloc5 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 , m5 , r5 ) ( PetscMalloc5 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 , m5 , r5 ) || PetscMemzero ( * ( r1 ) , ( m1 ) * sizeof ( * * ( r1 ) ) ) || PetscMemzero ( * ( r2 ) , ( m2 ) * sizeof ( * * ( r2 ) ) ) || PetscMemzero ( * ( r3 ) , ( m3 ) * sizeof ( * * ( r3 ) ) ) || PetscMemzero ( * ( r4 ) , ( m4 ) * sizeof ( * * ( r4 ) ) ) || PetscMemzero ( * ( r5 ) , ( m5 ) * sizeof ( * * ( r5 ) ) ) ) =#
#= # Skipping MacroDefinition: PetscMalloc6 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 , m5 , r5 , m6 , r6 ) ( PetscMalloc1 ( ( m1 ) , ( r1 ) ) || PetscMalloc1 ( ( m2 ) , ( r2 ) ) || PetscMalloc1 ( ( m3 ) , ( r3 ) ) || PetscMalloc1 ( ( m4 ) , ( r4 ) ) || PetscMalloc1 ( ( m5 ) , ( r5 ) ) || PetscMalloc1 ( ( m6 ) , ( r6 ) ) ) =#
#= # Skipping MacroDefinition: PetscCalloc6 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 , m5 , r5 , m6 , r6 ) ( PetscMalloc6 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 , m5 , r5 , m6 , r6 ) || PetscMemzero ( * ( r1 ) , ( m1 ) * sizeof ( * * ( r1 ) ) ) || PetscMemzero ( * ( r2 ) , ( m2 ) * sizeof ( * * ( r2 ) ) ) || PetscMemzero ( * ( r3 ) , ( m3 ) * sizeof ( * * ( r3 ) ) ) || PetscMemzero ( * ( r4 ) , ( m4 ) * sizeof ( * * ( r4 ) ) ) || PetscMemzero ( * ( r5 ) , ( m5 ) * sizeof ( * * ( r5 ) ) ) || PetscMemzero ( * ( r6 ) , ( m6 ) * sizeof ( * * ( r6 ) ) ) ) =#
#= # Skipping MacroDefinition: PetscMalloc7 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 , m5 , r5 , m6 , r6 , m7 , r7 ) ( PetscMalloc1 ( ( m1 ) , ( r1 ) ) || PetscMalloc1 ( ( m2 ) , ( r2 ) ) || PetscMalloc1 ( ( m3 ) , ( r3 ) ) || PetscMalloc1 ( ( m4 ) , ( r4 ) ) || PetscMalloc1 ( ( m5 ) , ( r5 ) ) || PetscMalloc1 ( ( m6 ) , ( r6 ) ) || PetscMalloc1 ( ( m7 ) , ( r7 ) ) ) =#
#= # Skipping MacroDefinition: PetscCalloc7 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 , m5 , r5 , m6 , r6 , m7 , r7 ) ( PetscMalloc7 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 , m5 , r5 , m6 , r6 , m7 , r7 ) || PetscMemzero ( * ( r1 ) , ( m1 ) * sizeof ( * * ( r1 ) ) ) || PetscMemzero ( * ( r2 ) , ( m2 ) * sizeof ( * * ( r2 ) ) ) || PetscMemzero ( * ( r3 ) , ( m3 ) * sizeof ( * * ( r3 ) ) ) || PetscMemzero ( * ( r4 ) , ( m4 ) * sizeof ( * * ( r4 ) ) ) || PetscMemzero ( * ( r5 ) , ( m5 ) * sizeof ( * * ( r5 ) ) ) || PetscMemzero ( * ( r6 ) , ( m6 ) * sizeof ( * * ( r6 ) ) ) || PetscMemzero ( * ( r7 ) , ( m7 ) * sizeof ( * * ( r7 ) ) ) ) =#
#= # Skipping MacroDefinition: PetscNew ( b ) PetscCalloc1 ( 1 , ( b ) ) =#
#= # Skipping MacroDefinition: PetscNewLog ( o , b ) ( PetscNew ( ( b ) ) || PetscLogObjectMemory ( ( PetscObject ) o , sizeof ( * * ( b ) ) ) ) =#
#= # Skipping MacroDefinition: PetscFree ( a ) ( ( * PetscTrFree ) ( ( void * ) ( a ) , __LINE__ , PETSC_FUNCTION_NAME , __FILE__ ) || ( ( a ) = 0 , 0 ) ) =#
#= # Skipping MacroDefinition: PetscFreeVoid ( a ) ( ( * PetscTrFree ) ( ( a ) , __LINE__ , PETSC_FUNCTION_NAME , __FILE__ ) , ( a ) = 0 ) =#
#= # Skipping MacroDefinition: PetscFree2 ( m1 , m2 ) ( PetscFree ( m2 ) || PetscFree ( m1 ) ) =#
#= # Skipping MacroDefinition: PetscFree3 ( m1 , m2 , m3 ) ( PetscFree ( m3 ) || PetscFree ( m2 ) || PetscFree ( m1 ) ) =#
#= # Skipping MacroDefinition: PetscFree4 ( m1 , m2 , m3 , m4 ) ( PetscFree ( m4 ) || PetscFree ( m3 ) || PetscFree ( m2 ) || PetscFree ( m1 ) ) =#
#= # Skipping MacroDefinition: PetscFree5 ( m1 , m2 , m3 , m4 , m5 ) ( PetscFree ( m5 ) || PetscFree ( m4 ) || PetscFree ( m3 ) || PetscFree ( m2 ) || PetscFree ( m1 ) ) =#
#= # Skipping MacroDefinition: PetscFree6 ( m1 , m2 , m3 , m4 , m5 , m6 ) ( PetscFree ( m6 ) || PetscFree ( m5 ) || PetscFree ( m4 ) || PetscFree ( m3 ) || PetscFree ( m2 ) || PetscFree ( m1 ) ) =#
#= # Skipping MacroDefinition: PetscFree7 ( m1 , m2 , m3 , m4 , m5 , m6 , m7 ) ( PetscFree ( m7 ) || PetscFree ( m6 ) || PetscFree ( m5 ) || PetscFree ( m4 ) || PetscFree ( m3 ) || PetscFree ( m2 ) || PetscFree ( m1 ) ) =#
#skipping undefined const MPIU_PETSCLOGDOUBLE = MPI_DOUBLE
#= # begin enum PetscDataType =#
typealias PetscDataType UInt32

const PETSC_INT = (UInt32)(0)
const PETSC_DOUBLE = (UInt32)(1)
const PETSC_COMPLEX = (UInt32)(2)
const PETSC_LONG = (UInt32)(3)
const PETSC_SHORT = (UInt32)(4)
const PETSC_FLOAT = (UInt32)(5)
const PETSC_CHAR = (UInt32)(6)
const PETSC_BIT_LOGICAL = (UInt32)(7)
const PETSC_ENUM = (UInt32)(8)
const PETSC_BOOL = (UInt32)(9)
const PETSC___FLOAT128 = (UInt32)(10)
const PETSC_OBJECT = (UInt32)(11)
const PETSC_FUNCTION = (UInt32)(12)
const PETSC_STRING = (UInt32)(12)

#= # end enum PetscDataType =#
const PETSC_SCALAR = PETSC_FLOAT
const PETSC_REAL = PETSC_FLOAT
const PETSC_FORTRANADDR = PETSC_LONG

#skipping undefined const MPIU_SUM = MPI_SUM
#skipping undefined const MPIU_MAX = MPI_MAX
#skipping undefined const MPIU_MIN = MPI_MIN
const PETSC_ERR_MIN_VALUE = 54
const PETSC_ERR_MEM = 55
const PETSC_ERR_SUP = 56
const PETSC_ERR_SUP_SYS = 57
const PETSC_ERR_ORDER = 58
const PETSC_ERR_SIG = 59
const PETSC_ERR_FP = 72
const PETSC_ERR_COR = 74
const PETSC_ERR_LIB = 76
const PETSC_ERR_PLIB = 77
const PETSC_ERR_MEMC = 78
const PETSC_ERR_CONV_FAILED = 82
const PETSC_ERR_USER = 83
const PETSC_ERR_SYS = 88
const PETSC_ERR_POINTER = 70
const PETSC_ERR_ARG_SIZ = 60
const PETSC_ERR_ARG_IDN = 61
const PETSC_ERR_ARG_WRONG = 62
const PETSC_ERR_ARG_CORRUPT = 64
const PETSC_ERR_ARG_OUTOFRANGE = 63
const PETSC_ERR_ARG_BADPTR = 68
const PETSC_ERR_ARG_NOTSAMETYPE = 69
const PETSC_ERR_ARG_NOTSAMECOMM = 80
const PETSC_ERR_ARG_WRONGSTATE = 73
const PETSC_ERR_ARG_TYPENOTSET = 89
const PETSC_ERR_ARG_INCOMP = 75
const PETSC_ERR_ARG_NULL = 85
const PETSC_ERR_ARG_UNKNOWN_TYPE = 86
const PETSC_ERR_FILE_OPEN = 65
const PETSC_ERR_FILE_READ = 66
const PETSC_ERR_FILE_WRITE = 67
const PETSC_ERR_FILE_UNEXPECTED = 79
const PETSC_ERR_MAT_LU_ZRPVT = 71
const PETSC_ERR_MAT_CH_ZRPVT = 81
const PETSC_ERR_INT_OVERFLOW = 84
const PETSC_ERR_FLOP_COUNT = 90
const PETSC_ERR_NOT_CONVERGED = 91
const PETSC_ERR_MISSING_FACTOR = 92
const PETSC_ERR_MAX_VALUE = 93

#= # Skipping MacroDefinition: PetscStringizeArg ( a ) # a =#
#= # Skipping MacroDefinition: PetscStringize ( a ) PetscStringizeArg ( a ) =#
#= # Skipping MacroDefinition: SETERRQ ( c , n , s ) =#
#= # Skipping MacroDefinition: SETERRQ1 ( c , n , s , a1 ) =#
#= # Skipping MacroDefinition: SETERRQ2 ( c , n , s , a1 , a2 ) =#
#= # Skipping MacroDefinition: SETERRQ3 ( c , n , s , a1 , a2 , a3 ) =#
#= # Skipping MacroDefinition: SETERRQ4 ( c , n , s , a1 , a2 , a3 , a4 ) =#
#= # Skipping MacroDefinition: SETERRQ5 ( c , n , s , a1 , a2 , a3 , a4 , a5 ) =#
#= # Skipping MacroDefinition: SETERRQ6 ( c , n , s , a1 , a2 , a3 , a4 , a5 , a6 ) =#
#= # Skipping MacroDefinition: SETERRQ7 ( c , n , s , a1 , a2 , a3 , a4 , a5 , a6 , a7 ) =#
#= # Skipping MacroDefinition: SETERRQ8 ( c , n , s , a1 , a2 , a3 , a4 , a5 , a6 , a7 , a8 ) =#
#= # Skipping MacroDefinition: SETERRABORT ( comm , n , s ) =#
#= # Skipping MacroDefinition: CHKERRQ ( n ) ; =#
#= # Skipping MacroDefinition: CHKERRABORT ( comm , n ) ; =#
#= # Skipping MacroDefinition: CHKERRCONTINUE ( n ) ; =#
const PETSCSTACKSIZE = 64

#= # Skipping MacroDefinition: PetscStackPushNoCheck ( funct , petsc_routine , hot ) do { } while ( 0 ) =#
#= # Skipping MacroDefinition: PetscStackPopNoCheck do { } while ( 0 ) =#
#= # Skipping MacroDefinition: PetscFunctionReturn ( a ) return ( a ) =#
#= # Skipping MacroDefinition: PetscFunctionReturnVoid ( ) return =#
#skipping undefined const PetscStackPop = CHKMEMQ
#= # Skipping MacroDefinition: PetscStackPush ( f ) CHKMEMQ =#
#= # Skipping MacroDefinition: PetscStackCall ( name , routine ) do { PetscStackPush ( name ) ; routine ; PetscStackPop ; } while ( 0 ) =#
#= # Skipping MacroDefinition: PetscStackCallStandard ( func , args ) do { PetscStackPush ( # func ) ; ierr = func args ; PetscStackPop ; if ( ierr ) SETERRQ1 ( PETSC_COMM_SELF , PETSC_ERR_LIB , "Error in %s()" , # func ) ; } while ( 0 ) =#
const PETSC_SMALLEST_CLASSID = 1211211

#= # Skipping MacroDefinition: PetscObjectComposeFunction ( a , b , d ) PetscObjectComposeFunction_Private ( a , b , ( PetscVoidFunction ) ( d ) ) =#
#= # Skipping MacroDefinition: PetscOptionsBegin ( comm , prefix , mess , sec ) 0 ; do { PetscOptions PetscOptionsObjectBase ; PetscOptions * PetscOptionsObject = & PetscOptionsObjectBase ; PetscMemzero ( PetscOptionsObject , sizeof ( PetscOptions ) ) ; for ( PetscOptionsObject -> count = ( PetscOptionsPublish ? - 1 : 1 ) ; PetscOptionsObject -> count < 2 ; PetscOptionsObject -> count ++ ) { PetscErrorCode _5_ierr = PetscOptionsBegin_Private ( PetscOptionsObject , comm , prefix , mess , sec ) ; CHKERRQ ( _5_ierr ) ; =#
#= # Skipping MacroDefinition: PetscObjectOptionsBegin ( obj ) 0 ; do { PetscOptions PetscOptionsObjectBase ; PetscOptions * PetscOptionsObject = & PetscOptionsObjectBase ; for ( PetscOptionsObject -> count = ( PetscOptionsPublish ? - 1 : 1 ) ; PetscOptionsObject -> count < 2 ; PetscOptionsObject -> count ++ ) { PetscErrorCode _5_ierr = PetscObjectOptionsBegin_Private ( PetscOptionsObject , obj ) ; CHKERRQ ( _5_ierr ) ; =#
#= # Skipping MacroDefinition: PetscOptionsEnd ( ) _5_ierr = PetscOptionsEnd_Private ( PetscOptionsObject ) ; CHKERRQ ( _5_ierr ) ; } } while ( 0 ) =#
#= # Skipping MacroDefinition: PetscOptionsTail ( ) 0 ; { if ( PetscOptionsObject -> count != 1 ) PetscFunctionReturn ( 0 ) ; } =#
#= # Skipping MacroDefinition: PetscOptionsEnum ( a , b , c , d , e , f , g ) PetscOptionsEnum_Private ( PetscOptionsObject , a , b , c , d , e , f , g ) =#
#= # Skipping MacroDefinition: PetscOptionsInt ( a , b , c , d , e , f ) PetscOptionsInt_Private ( PetscOptionsObject , a , b , c , d , e , f ) =#
#= # Skipping MacroDefinition: PetscOptionsReal ( a , b , c , d , e , f ) PetscOptionsReal_Private ( PetscOptionsObject , a , b , c , d , e , f ) =#
#= # Skipping MacroDefinition: PetscOptionsScalar ( a , b , c , d , e , f ) PetscOptionsScalar_Private ( PetscOptionsObject , a , b , c , d , e , f ) =#
#= # Skipping MacroDefinition: PetscOptionsName ( a , b , c , d ) PetscOptionsName_Private ( PetscOptionsObject , a , b , c , d ) =#
#= # Skipping MacroDefinition: PetscOptionsString ( a , b , c , d , e , f , g ) PetscOptionsString_Private ( PetscOptionsObject , a , b , c , d , e , f , g ) =#
#= # Skipping MacroDefinition: PetscOptionsBool ( a , b , c , d , e , f ) PetscOptionsBool_Private ( PetscOptionsObject , a , b , c , d , e , f ) =#
#= # Skipping MacroDefinition: PetscOptionsBoolGroupBegin ( a , b , c , d ) PetscOptionsBoolGroupBegin_Private ( PetscOptionsObject , a , b , c , d ) =#
#= # Skipping MacroDefinition: PetscOptionsBoolGroup ( a , b , c , d ) PetscOptionsBoolGroup_Private ( PetscOptionsObject , a , b , c , d ) =#
#= # Skipping MacroDefinition: PetscOptionsBoolGroupEnd ( a , b , c , d ) PetscOptionsBoolGroupEnd_Private ( PetscOptionsObject , a , b , c , d ) =#
#= # Skipping MacroDefinition: PetscOptionsFList ( a , b , c , d , e , f , g , h ) PetscOptionsFList_Private ( PetscOptionsObject , a , b , c , d , e , f , g , h ) =#
#= # Skipping MacroDefinition: PetscOptionsEList ( a , b , c , d , e , f , g , h ) PetscOptionsEList_Private ( PetscOptionsObject , a , b , c , d , e , f , g , h ) =#
#= # Skipping MacroDefinition: PetscOptionsRealArray ( a , b , c , d , e , f ) PetscOptionsRealArray_Private ( PetscOptionsObject , a , b , c , d , e , f ) =#
#= # Skipping MacroDefinition: PetscOptionsScalarArray ( a , b , c , d , e , f ) PetscOptionsScalarArray_Private ( PetscOptionsObject , a , b , c , d , e , f ) =#
#= # Skipping MacroDefinition: PetscOptionsIntArray ( a , b , c , d , e , f ) PetscOptionsIntArray_Private ( PetscOptionsObject , a , b , c , d , e , f ) =#
#= # Skipping MacroDefinition: PetscOptionsStringArray ( a , b , c , d , e , f ) PetscOptionsStringArray_Private ( PetscOptionsObject , a , b , c , d , e , f ) =#
#= # Skipping MacroDefinition: PetscOptionsBoolArray ( a , b , c , d , e , f ) PetscOptionsBoolArray_Private ( PetscOptionsObject , a , b , c , d , e , f ) =#
#= # Skipping MacroDefinition: PetscOptionsEnumArray ( a , b , c , d , e , f , g ) PetscOptionsEnumArray_Private ( PetscOptionsObject , a , b , c , d , e , f , g ) =#
#= # Skipping MacroDefinition: PetscObjectQueryFunction ( obj , name , fptr ) PetscObjectQueryFunction_Private ( ( obj ) , ( name ) , ( PetscVoidFunction * ) ( fptr ) ) =#
#= # Skipping MacroDefinition: PetscSAWsBlock ( ) 0 =#
#= # Skipping MacroDefinition: PetscObjectSAWsViewOff ( obj ) 0 =#
#= # Skipping MacroDefinition: PetscObjectSAWsSetBlock ( obj , flg ) 0 =#
#= # Skipping MacroDefinition: PetscObjectSAWsBlock ( obj ) 0 =#
#= # Skipping MacroDefinition: PetscObjectSAWsGrantAccess ( obj ) 0 =#
#= # Skipping MacroDefinition: PetscObjectSAWsTakeAccess ( obj ) 0 =#
#= # Skipping MacroDefinition: PetscStackViewSAWs ( ) 0 =#
#= # Skipping MacroDefinition: PetscStackSAWsViewOff ( ) 0 =#
#= # Skipping MacroDefinition: PetscFunctionListAdd ( list , name , fptr ) PetscFunctionListAdd_Private ( ( list ) , ( name ) , ( PetscVoidFunction ) ( fptr ) ) =#
#= # Skipping MacroDefinition: PetscFunctionListFind ( list , name , fptr ) PetscFunctionListFind_Private ( ( list ) , ( name ) , ( PetscVoidFunction * ) ( fptr ) ) =#
#= # Skipping MacroDefinition: PetscNot ( a ) ( ( a ) ? PETSC_FALSE : PETSC_TRUE ) =#
const PETSC_EVENT = 1311311

#= # Skipping MacroDefinition: PetscInfo ( A , S ) 0 =#
#= # Skipping MacroDefinition: PetscInfo1 ( A , S , a1 ) 0 =#
#= # Skipping MacroDefinition: PetscInfo2 ( A , S , a1 , a2 ) 0 =#
#= # Skipping MacroDefinition: PetscInfo3 ( A , S , a1 , a2 , a3 ) 0 =#
#= # Skipping MacroDefinition: PetscInfo4 ( A , S , a1 , a2 , a3 , a4 ) 0 =#
#= # Skipping MacroDefinition: PetscInfo5 ( A , S , a1 , a2 , a3 , a4 , a5 ) 0 =#
#= # Skipping MacroDefinition: PetscInfo6 ( A , S , a1 , a2 , a3 , a4 , a5 , a6 ) 0 =#
#= # Skipping MacroDefinition: PetscInfo7 ( A , S , a1 , a2 , a3 , a4 , a5 , a6 , a7 ) 0 =#
#= # Skipping MacroDefinition: PetscLogFlops ( n ) 0 =#
#= # Skipping MacroDefinition: PetscLogEventActivate ( a ) 0 =#
#= # Skipping MacroDefinition: PetscLogEventDeactivate ( a ) 0 =#
#= # Skipping MacroDefinition: PetscLogEventActivateClass ( a ) 0 =#
#= # Skipping MacroDefinition: PetscLogEventDeactivateClass ( a ) 0 =#
#= # Skipping MacroDefinition: PetscLogEventSetActiveAll ( a , b ) 0 =#
const PetscLogPLB = 0
const PetscLogPLE = 0
const PetscLogPHC = 0
const PetscLogPHD = 0

#= # Skipping MacroDefinition: PetscGetFlops ( a ) ( * ( a ) = 0.0 , 0 ) =#
#= # Skipping MacroDefinition: PetscLogEventBegin ( e , o1 , o2 , o3 , o4 ) 0 =#
#= # Skipping MacroDefinition: PetscLogEventEnd ( e , o1 , o2 , o3 , o4 ) 0 =#
#= # Skipping MacroDefinition: PetscLogEventBarrierBegin ( e , o1 , o2 , o3 , o4 , cm ) 0 =#
#= # Skipping MacroDefinition: PetscLogEventBarrierEnd ( e , o1 , o2 , o3 , o4 , cm ) 0 =#
#= # Skipping MacroDefinition: PetscLogObjectParents ( p , n , c ) 0 =#
#= # Skipping MacroDefinition: PetscLogObjectCreate ( h ) 0 =#
#= # Skipping MacroDefinition: PetscLogObjectDestroy ( h ) 0 =#
#= # Skipping MacroDefinition: PetscLogDestroy ( ) 0 =#
#= # Skipping MacroDefinition: PetscLogStagePush ( a ) 0 =#
#= # Skipping MacroDefinition: PetscLogStagePop ( ) 0 =#
#= # Skipping MacroDefinition: PetscLogStageRegister ( a , b ) 0 =#
#= # Skipping MacroDefinition: PetscLogStagePrint ( a , flg ) 0 =#
#= # Skipping MacroDefinition: PetscLogView ( viewer ) 0 =#
#= # Skipping MacroDefinition: PetscLogViewFromOptions ( ) 0 =#
#= # Skipping MacroDefinition: PetscLogBegin ( ) 0 =#
#= # Skipping MacroDefinition: PetscLogTraceBegin ( file ) 0 =#
#= # Skipping MacroDefinition: PetscLogSet ( lb , le ) 0 =#
#= # Skipping MacroDefinition: PetscLogAllBegin ( ) 0 =#
#= # Skipping MacroDefinition: PetscLogDump ( c ) 0 =#
#= # Skipping MacroDefinition: PetscLogEventRegister ( a , b , c ) 0 =#
#= # Skipping MacroDefinition: PetscLogObjects ( a ) 0 =#
#= # Skipping MacroDefinition: PetscLogActions ( a ) 0 =#
#= # Skipping MacroDefinition: MPI_Startall_irecv ( count , number , requests ) MPI_Startall ( number , requests ) =#
#= # Skipping MacroDefinition: MPI_Startall_isend ( count , number , requests ) MPI_Startall ( number , requests ) =#
#= # Skipping MacroDefinition: MPI_Start_isend ( count , requests ) MPI_Start ( requests ) =#
#= # Skipping MacroDefinition: PetscLogStageGetId ( a , b ) ( * ( b ) = 0 , 0 ) =#
#= # Skipping MacroDefinition: PetscLogStageSetActive ( a , b ) 0 =#
#= # Skipping MacroDefinition: PetscLogStageGetActive ( a , b ) 0 =#
#= # Skipping MacroDefinition: PetscLogStageGetVisible ( a , b ) 0 =#
#= # Skipping MacroDefinition: PetscLogStageSetVisible ( a , b ) 0 =#
#= # Skipping MacroDefinition: PetscPreLoadBegin ( flag , name ) do { PetscBool PetscPreLoading = flag ; int PetscPreLoadMax , PetscPreLoadIt ; PetscLogStage _stageNum ; PetscErrorCode _3_ierr ; _3_ierr = PetscOptionsGetBool ( NULL , "-preload" , & PetscPreLoading , NULL ) ; CHKERRQ ( _3_ierr ) ; PetscPreLoadMax = ( int ) ( PetscPreLoading ) ; PetscPreLoadingUsed = PetscPreLoading ? PETSC_TRUE : PetscPreLoadingUsed ; for ( PetscPreLoadIt = 0 ; PetscPreLoadIt <= PetscPreLoadMax ; PetscPreLoadIt ++ ) { PetscPreLoadingOn = PetscPreLoading ; _3_ierr = PetscBarrier ( NULL ) ; CHKERRQ ( _3_ierr ) ; if ( PetscPreLoadIt > 0 ) { _3_ierr = PetscLogStageGetId ( name , & _stageNum ) ; CHKERRQ ( _3_ierr ) ; } else { _3_ierr = PetscLogStageRegister ( name , & _stageNum ) ; CHKERRQ ( _3_ierr ) ; } _3_ierr = PetscLogStageSetActive ( _stageNum , ( PetscBool ) ( ! PetscPreLoadMax || PetscPreLoadIt ) ) ; _3_ierr = PetscLogStagePush ( _stageNum ) ; CHKERRQ ( _3_ierr ) ; =#
#= # Skipping MacroDefinition: PetscPreLoadEnd ( ) _3_ierr = PetscLogStagePop ( ) ; CHKERRQ ( _3_ierr ) ; PetscPreLoading = PETSC_FALSE ; } \
#} while ( 0 ) =#
#= # Skipping MacroDefinition: PetscPreLoadStage ( name ) do { _3_ierr = PetscLogStagePop ( ) ; CHKERRQ ( _3_ierr ) ; if ( PetscPreLoadIt > 0 ) { _3_ierr = PetscLogStageGetId ( name , & _stageNum ) ; CHKERRQ ( _3_ierr ) ; } else { _3_ierr = PetscLogStageRegister ( name , & _stageNum ) ; CHKERRQ ( _3_ierr ) ; } _3_ierr = PetscLogStageSetActive ( _stageNum , ( PetscBool ) ( ! PetscPreLoadMax || PetscPreLoadIt ) ) ; _3_ierr = PetscLogStagePush ( _stageNum ) ; CHKERRQ ( _3_ierr ) ; } while ( 0 ) =#
#= # Skipping MacroDefinition: PetscPrefetchBlock ( a , n , rw , t ) do { const char * _p = ( const char * ) ( a ) , * _end = ( const char * ) ( ( a ) + ( n ) ) ; for ( ; _p < _end ; _p += PETSC_LEVEL1_DCACHE_LINESIZE ) PETSC_Prefetch ( _p , ( rw ) , ( t ) ) ; } while ( 0 ) =#
const PETSC_MPI_INT_MAX = 2147483647
const PETSC_MPI_INT_MIN = -2147483647
const PETSC_BLAS_INT_MAX = 2147483647
const PETSC_BLAS_INT_MIN = -2147483647
const PETSC_MAX_PATH_LEN = 4096
const PETSCRAND = symbol("rand")
const PETSCRAND48 = symbol("rand48")
const PETSCSPRNG = symbol("sprng")
const PETSC_BINARY_INT_SIZE = 32 / 8
const PETSC_BINARY_FLOAT_SIZE = 32 / 8
const PETSC_BINARY_CHAR_SIZE = 8 / 8
const PETSC_BINARY_SHORT_SIZE = 16 / 8
const PETSC_BINARY_DOUBLE_SIZE = 64 / 8

#= # Skipping MacroDefinition: PETSC_BINARY_SCALAR_SIZE sizeof ( PetscScalar ) =#
const PETSC_BAG_FILE_CLASSID = 1211219
const PETSCVIEWERSOCKET = symbol("socket")
const PETSCVIEWERASCII = symbol("ascii")
const PETSCVIEWERBINARY = symbol("binary")
const PETSCVIEWERSTRING = symbol("string")
const PETSCVIEWERDRAW = symbol("draw")
const PETSCVIEWERVU = symbol("vu")
const PETSCVIEWERMATHEMATICA = symbol("mathematica")
const PETSCVIEWERNETCDF = symbol("netcdf")
const PETSCVIEWERHDF5 = symbol("hdf5")
const PETSCVIEWERVTK = symbol("vtk")
const PETSCVIEWERMATLAB = symbol("matlab")
const PETSCVIEWERSAWS = symbol("saws")
const PETSC_DRAW_X = symbol("x")
const PETSC_DRAW_GLUT = symbol("glut")
const PETSC_DRAW_OPENGLES = symbol("opengles")
const PETSC_DRAW_NULL = symbol("null")
const PETSC_DRAW_WIN32 = symbol("win32")
const PETSC_DRAW_TIKZ = symbol("tikz")

#= # Skipping MacroDefinition: PetscOptionsViewer ( a , b , c , d , e , f ) PetscOptionsViewer_Private ( PetscOptionsObject , a , b , c , d , e , f ) ; =#
#= # Skipping MacroDefinition: PETSC_VIEWER_STDERR_SELF PETSC_VIEWER_STDERR_ ( PETSC_COMM_SELF ) =#
#= # Skipping MacroDefinition: PETSC_VIEWER_STDERR_WORLD PETSC_VIEWER_STDERR_ ( PETSC_COMM_WORLD ) =#
#= # Skipping MacroDefinition: PETSC_VIEWER_STDOUT_WORLD PETSC_VIEWER_STDOUT_ ( PETSC_COMM_WORLD ) =#
#= # Skipping MacroDefinition: PETSC_VIEWER_STDOUT_SELF PETSC_VIEWER_STDOUT_ ( PETSC_COMM_SELF ) =#
#= # Skipping MacroDefinition: PETSC_VIEWER_DRAW_WORLD PETSC_VIEWER_DRAW_ ( PETSC_COMM_WORLD ) =#
#= # Skipping MacroDefinition: PETSC_VIEWER_DRAW_SELF PETSC_VIEWER_DRAW_ ( PETSC_COMM_SELF ) =#
#= # Skipping MacroDefinition: PETSC_VIEWER_SOCKET_WORLD PETSC_VIEWER_SOCKET_ ( PETSC_COMM_WORLD ) =#
#= # Skipping MacroDefinition: PETSC_VIEWER_SOCKET_SELF PETSC_VIEWER_SOCKET_ ( PETSC_COMM_SELF ) =#
#= # Skipping MacroDefinition: PETSC_VIEWER_BINARY_WORLD PETSC_VIEWER_BINARY_ ( PETSC_COMM_WORLD ) =#
#= # Skipping MacroDefinition: PETSC_VIEWER_BINARY_SELF PETSC_VIEWER_BINARY_ ( PETSC_COMM_SELF ) =#
#= # Skipping MacroDefinition: PETSC_VIEWER_MATLAB_WORLD PETSC_VIEWER_MATLAB_ ( PETSC_COMM_WORLD ) =#
#= # Skipping MacroDefinition: PETSC_VIEWER_MATLAB_SELF PETSC_VIEWER_MATLAB_ ( PETSC_COMM_SELF ) =#
#= # Skipping MacroDefinition: PETSC_VIEWER_MATHEMATICA_WORLD ( PetscViewerInitializeMathematicaWorld_Private ( ) , PETSC_VIEWER_MATHEMATICA_WORLD_PRIVATE ) =#
const PETSC_HASH_FACT = 79943

#= # Skipping MacroDefinition: PETSC_MATLAB_ENGINE_WORLD PETSC_MATLAB_ENGINE_ ( PETSC_COMM_WORLD ) =#
#= # Skipping MacroDefinition: PETSC_MATLAB_ENGINE_SELF PETSC_MATLAB_ENGINE_ ( PETSC_COMM_SELF ) =#
const PETSC_DRAW_BASIC_COLORS = 33
const PETSC_DRAW_ROTATE = -1
const PETSC_DRAW_WHITE = 0
const PETSC_DRAW_BLACK = 1
const PETSC_DRAW_RED = 2
const PETSC_DRAW_GREEN = 3
const PETSC_DRAW_CYAN = 4
const PETSC_DRAW_BLUE = 5
const PETSC_DRAW_MAGENTA = 6
const PETSC_DRAW_AQUAMARINE = 7
const PETSC_DRAW_FORESTGREEN = 8
const PETSC_DRAW_ORANGE = 9
const PETSC_DRAW_VIOLET = 10
const PETSC_DRAW_BROWN = 11
const PETSC_DRAW_PINK = 12
const PETSC_DRAW_CORAL = 13
const PETSC_DRAW_GRAY = 14
const PETSC_DRAW_YELLOW = 15
const PETSC_DRAW_GOLD = 16
const PETSC_DRAW_LIGHTPINK = 17
const PETSC_DRAW_MEDIUMTURQUOISE = 18
const PETSC_DRAW_KHAKI = 19
const PETSC_DRAW_DIMGRAY = 20
const PETSC_DRAW_YELLOWGREEN = 21
const PETSC_DRAW_SKYBLUE = 22
const PETSC_DRAW_DARKGREEN = 23
const PETSC_DRAW_NAVYBLUE = 24
const PETSC_DRAW_SANDYBROWN = 25
const PETSC_DRAW_CADETBLUE = 26
const PETSC_DRAW_POWDERBLUE = 27
const PETSC_DRAW_DEEPPINK = 28
const PETSC_DRAW_THISTLE = 29
const PETSC_DRAW_LIMEGREEN = 30
const PETSC_DRAW_LAVENDERBLUSH = 31
const PETSC_DRAW_PLUM = 32
const PETSC_DRAW_FULL_SIZE = -3
const PETSC_DRAW_HALF_SIZE = -4
const PETSC_DRAW_THIRD_SIZE = -5
const PETSC_DRAW_QUARTER_SIZE = -6
const IS_FILE_CLASSID = 1211218
const ISGENERAL = symbol("general")
const ISSTRIDE = symbol("stride")
const ISBLOCK = symbol("block")
const VECSEQ = symbol("seq")
const VECMPI = symbol("mpi")
const VECSTANDARD = symbol("standard")
const VECSHARED = symbol("shared")
const VECSEQCUSP = symbol("seqcusp")
const VECMPICUSP = symbol("mpicusp")
const VECCUSP = symbol("cusp")
const VECSEQVIENNACL = symbol("seqviennacl")
const VECMPIVIENNACL = symbol("mpiviennacl")
const VECVIENNACL = symbol("viennacl")
const VECNEST = symbol("nest")
const VECSEQPTHREAD = symbol("seqpthread")
const VECMPIPTHREAD = symbol("mpipthread")
const VECPTHREAD = symbol("pthread")
const VEC_FILE_CLASSID = 1211214

#= # begin enum NormType =#
typealias NormType UInt32

const NORM_1 = (UInt32)(0)
const NORM_2 = (UInt32)(1)
const NORM_FROBENIUS = (UInt32)(2)
const NORM_INFINITY = (UInt32)(3)
const NORM_1_AND_2 = (UInt32)(4)

#= # end enum NormType =#
const NORM_MAX = NORM_INFINITY

#= # Skipping MacroDefinition: VecLockGet ( x , arg ) * ( arg ) = 0 =#
#= # Skipping MacroDefinition: VecLockPush ( x ) 0 =#
#= # Skipping MacroDefinition: VecLockPop ( x ) 0 =#
#= # Skipping MacroDefinition: VecLocked ( x , arg ) =#
const MATSAME = symbol("same")
const MATMAIJ = symbol("maij")
const MATSEQMAIJ = symbol("seqmaij")
const MATMPIMAIJ = symbol("mpimaij")
const MATIS = symbol("is")
const MATAIJ = symbol("aij")
const MATSEQAIJ = symbol("seqaij")
const MATSEQAIJPTHREAD = symbol("seqaijpthread")
const MATAIJPTHREAD = symbol("aijpthread")
const MATMPIAIJ = symbol("mpiaij")
const MATAIJCRL = symbol("aijcrl")
const MATSEQAIJCRL = symbol("seqaijcrl")
const MATMPIAIJCRL = symbol("mpiaijcrl")
const MATAIJCUSP = symbol("aijcusp")
const MATSEQAIJCUSP = symbol("seqaijcusp")
const MATMPIAIJCUSP = symbol("mpiaijcusp")
const MATAIJCUSPARSE = symbol("aijcusparse")
const MATSEQAIJCUSPARSE = symbol("seqaijcusparse")
const MATMPIAIJCUSPARSE = symbol("mpiaijcusparse")
const MATAIJVIENNACL = symbol("aijviennacl")
const MATSEQAIJVIENNACL = symbol("seqaijviennacl")
const MATMPIAIJVIENNACL = symbol("mpiaijviennacl")
const MATAIJPERM = symbol("aijperm")
const MATSEQAIJPERM = symbol("seqaijperm")
const MATMPIAIJPERM = symbol("mpiaijperm")
const MATSHELL = symbol("shell")
const MATDENSE = symbol("dense")
const MATSEQDENSE = symbol("seqdense")
const MATMPIDENSE = symbol("mpidense")
const MATELEMENTAL = symbol("elemental")
const MATBAIJ = symbol("baij")
const MATSEQBAIJ = symbol("seqbaij")
const MATMPIBAIJ = symbol("mpibaij")
const MATMPIADJ = symbol("mpiadj")
const MATSBAIJ = symbol("sbaij")
const MATSEQSBAIJ = symbol("seqsbaij")
const MATMPISBAIJ = symbol("mpisbaij")
const MATSEQBSTRM = symbol("seqbstrm")
const MATMPIBSTRM = symbol("mpibstrm")
const MATBSTRM = symbol("bstrm")
const MATSEQSBSTRM = symbol("seqsbstrm")
const MATMPISBSTRM = symbol("mpisbstrm")
const MATSBSTRM = symbol("sbstrm")
const MATDAAD = symbol("daad")
const MATMFFD = symbol("mffd")
const MATNORMAL = symbol("normal")
const MATLRC = symbol("lrc")
const MATSCATTER = symbol("scatter")
const MATBLOCKMAT = symbol("blockmat")
const MATCOMPOSITE = symbol("composite")
const MATFFT = symbol("fft")
const MATFFTW = symbol("fftw")
const MATSEQCUFFT = symbol("seqcufft")
const MATTRANSPOSEMAT = symbol("transpose")
const MATSCHURCOMPLEMENT = symbol("schurcomplement")
const MATPYTHON = symbol("python")
const MATHYPRESTRUCT = symbol("hyprestruct")
const MATHYPRESSTRUCT = symbol("hypresstruct")
const MATSUBMATRIX = symbol("submatrix")
const MATLOCALREF = symbol("localref")
const MATNEST = symbol("nest")

#= # Skipping MacroDefinition: MatSolverPackage char * =#
const MATSOLVERSUPERLU = symbol("superlu")
const MATSOLVERSUPERLU_DIST = symbol("superlu_dist")
const MATSOLVERUMFPACK = symbol("umfpack")
const MATSOLVERCHOLMOD = symbol("cholmod")
const MATSOLVERESSL = symbol("essl")
const MATSOLVERLUSOL = symbol("lusol")
const MATSOLVERMUMPS = symbol("mumps")
const MATSOLVERMKL_PARDISO = symbol("mkl_pardiso")
const MATSOLVERMKL_CPARDISO = symbol("mkl_cpardiso")
const MATSOLVERPASTIX = symbol("pastix")
const MATSOLVERMATLAB = symbol("matlab")
const MATSOLVERPETSC = symbol("petsc")
const MATSOLVERBAS = symbol("bas")
const MATSOLVERCUSPARSE = symbol("cusparse")
const MATSOLVERBSTRM = symbol("bstrm")
const MATSOLVERSBSTRM = symbol("sbstrm")
const MATSOLVERELEMENTAL = symbol("elemental")
const MATSOLVERCLIQUE = symbol("clique")
const MATSOLVERKLU = symbol("klu")
const MAT_FILE_CLASSID = 1211216

#= # Skipping MacroDefinition: MatPreallocateInitialize ( comm , nrows , ncols , dnz , onz ) 0 ; \
#{ PetscErrorCode _4_ierr ; PetscInt __nrows = ( nrows ) , __ctmp = ( ncols ) , __rstart , __start , __end ; _4_ierr = PetscCalloc2 ( __nrows , & dnz , __nrows , & onz ) ; CHKERRQ ( _4_ierr ) ; __start = 0 ; __end = __start ; _4_ierr = MPI_Scan ( & __ctmp , & __end , 1 , MPIU_INT , MPI_SUM , comm ) ; CHKERRQ ( _4_ierr ) ; __start = __end - __ctmp ; _4_ierr = MPI_Scan ( & __nrows , & __rstart , 1 , MPIU_INT , MPI_SUM , comm ) ; CHKERRQ ( _4_ierr ) ; __rstart = __rstart - __nrows ; =#
#= # Skipping MacroDefinition: MatPreallocateSetLocal ( rmap , nrows , rows , cmap , ncols , cols , dnz , onz ) 0 ; \
#{ PetscInt __l ; _4_ierr = ISLocalToGlobalMappingApply ( rmap , nrows , rows , rows ) ; CHKERRQ ( _4_ierr ) ; _4_ierr = ISLocalToGlobalMappingApply ( cmap , ncols , cols , cols ) ; CHKERRQ ( _4_ierr ) ; for ( __l = 0 ; __l < nrows ; __l ++ ) { _4_ierr = MatPreallocateSet ( ( rows ) [ __l ] , ncols , cols , dnz , onz ) ; CHKERRQ ( _4_ierr ) ; } \
#} =#
#= # Skipping MacroDefinition: MatPreallocateSetLocalBlock ( rmap , nrows , rows , cmap , ncols , cols , dnz , onz ) 0 ; \
#{ PetscInt __l ; _4_ierr = ISLocalToGlobalMappingApplyBlock ( rmap , nrows , rows , rows ) ; CHKERRQ ( _4_ierr ) ; _4_ierr = ISLocalToGlobalMappingApplyBlock ( cmap , ncols , cols , cols ) ; CHKERRQ ( _4_ierr ) ; for ( __l = 0 ; __l < nrows ; __l ++ ) { _4_ierr = MatPreallocateSet ( ( rows ) [ __l ] , ncols , cols , dnz , onz ) ; CHKERRQ ( _4_ierr ) ; } \
#} =#
#= # Skipping MacroDefinition: MatPreallocateSymmetricSetLocalBlock ( map , nrows , rows , ncols , cols , dnz , onz ) 0 ; \
#{ PetscInt __l ; _4_ierr = ISLocalToGlobalMappingApplyBlock ( map , nrows , rows , rows ) ; CHKERRQ ( _4_ierr ) ; _4_ierr = ISLocalToGlobalMappingApplyBlock ( map , ncols , cols , cols ) ; CHKERRQ ( _4_ierr ) ; for ( __l = 0 ; __l < nrows ; __l ++ ) { _4_ierr = MatPreallocateSymmetricSetBlock ( ( rows ) [ __l ] , ncols , cols , dnz , onz ) ; CHKERRQ ( _4_ierr ) ; } \
#} =#
#= # Skipping MacroDefinition: MatPreallocateSet ( row , nc , cols , dnz , onz ) 0 ; \
#{ PetscInt __i ; if ( row < __rstart ) SETERRQ2 ( PETSC_COMM_SELF , PETSC_ERR_ARG_OUTOFRANGE , "Trying to set preallocation for row %D less than first local row %D" , row , __rstart ) ; if ( row >= __rstart + __nrows ) SETERRQ2 ( PETSC_COMM_SELF , PETSC_ERR_ARG_OUTOFRANGE , "Trying to set preallocation for row %D greater than last local row %D" , row , __rstart + __nrows - 1 ) ; for ( __i = 0 ; __i < nc ; __i ++ ) { if ( ( cols ) [ __i ] < __start || ( cols ) [ __i ] >= __end ) onz [ row - __rstart ] ++ ; else dnz [ row - __rstart ] ++ ; } \
#} =#
#= # Skipping MacroDefinition: MatPreallocateSymmetricSetBlock ( row , nc , cols , dnz , onz ) 0 ; \
#{ PetscInt __i ; for ( __i = 0 ; __i < nc ; __i ++ ) { if ( cols [ __i ] >= __end ) onz [ row - __rstart ] ++ ; else if ( cols [ __i ] >= row ) dnz [ row - __rstart ] ++ ; } \
#} =#
#= # Skipping MacroDefinition: MatPreallocateLocation ( A , row , ncols , cols , dnz , onz ) 0 ; if ( A ) { ierr = MatSetValues ( A , 1 , & row , ncols , cols , NULL , INSERT_VALUES ) ; CHKERRQ ( ierr ) ; } else { ierr = MatPreallocateSet ( row , ncols , cols , dnz , onz ) ; CHKERRQ ( ierr ) ; } =#
#= # Skipping MacroDefinition: MatPreallocateFinalize ( dnz , onz ) 0 ; _4_ierr = PetscFree2 ( dnz , onz ) ; CHKERRQ ( _4_ierr ) ; } =#
const MAT_SKIP_ALLOCATION = -4
const MATORDERINGNATURAL = symbol("natural")
const MATORDERINGND = symbol("nd")
const MATORDERING1WD = symbol("1wd")
const MATORDERINGRCM = symbol("rcm")
const MATORDERINGQMD = symbol("qmd")
const MATORDERINGROWLENGTH = symbol("rowlength")
const MATORDERINGWBM = symbol("wbm")
const MATORDERINGSPECTRAL = symbol("spectral")
const MATORDERINGAMD = symbol("amd")
const MATCOLORINGJP = symbol("jp")
const MATCOLORINGPOWER = symbol("power")
const MATCOLORINGNATURAL = symbol("natural")
const MATCOLORINGSL = symbol("sl")
const MATCOLORINGLF = symbol("lf")
const MATCOLORINGID = symbol("id")
const MATCOLORINGGREEDY = symbol("greedy")
const MATPARTITIONINGCURRENT = symbol("current")
const MATPARTITIONINGSQUARE = symbol("square")
const MATPARTITIONINGPARMETIS = symbol("parmetis")
const MATPARTITIONINGCHACO = symbol("chaco")
const MATPARTITIONINGPARTY = symbol("party")
const MATPARTITIONINGPTSCOTCH = symbol("ptscotch")
const MP_PARTY_OPT = symbol("opt")
const MP_PARTY_LIN = symbol("lin")
const MP_PARTY_SCA = symbol("sca")
const MP_PARTY_RAN = symbol("ran")
const MP_PARTY_GBF = symbol("gbf")
const MP_PARTY_GCF = symbol("gcf")
const MP_PARTY_BUB = symbol("bub")
const MP_PARTY_DEF = symbol("def")
const MP_PARTY_HELPFUL_SETS = symbol("hs")
const MP_PARTY_KERNIGHAN_LIN = symbol("kl")
const MP_PARTY_NONE = symbol("no")
const MATCOARSENMIS = symbol("mis")
const MATCOARSENHEM = symbol("hem")
const MATRIX_BINARY_FORMAT_DENSE = -1
const MATMFFD_DS = symbol("ds")
const MATMFFD_WP = symbol("wp")
const DMDA = symbol("da")
const DMCOMPOSITE = symbol("composite")
const DMSLICED = symbol("sliced")
const DMSHELL = symbol("shell")
const DMPLEX = symbol("plex")
const DMCARTESIAN = symbol("cartesian")
const DMREDUNDANT = symbol("redundant")
const DMPATCH = symbol("patch")
const DMMOAB = symbol("moab")
const DMNETWORK = symbol("network")
const DM_FILE_CLASSID = 1211221
const PFCONSTANT = symbol("constant")
const PFMAT = symbol("mat")
const PFSTRING = symbol("string")
const PFQUICK = symbol("quick")
const PFIDENTITY = symbol("identity")
const PFMATLAB = symbol("matlab")

#= # Skipping MacroDefinition: PFSetOptionsPrefix ( a , s ) PetscObjectSetOptionsPrefix ( ( PetscObject ) ( a ) , s ) =#
const AOBASIC = symbol("basic")
const AOADVANCED = symbol("advanced")
const AOMAPPING = symbol("mapping")
const AOMEMORYSCALABLE = symbol("memoryscalable")
const PETSCSPACEPOLYNOMIAL = symbol("poly")
const PETSCSPACEDG = symbol("dg")
const PETSCDUALSPACELAGRANGE = symbol("lagrange")
const PETSCDUALSPACESIMPLE = symbol("simple")
const PETSCFEBASIC = symbol("basic")
const PETSCFENONAFFINE = symbol("nonaffine")
const PETSCFEOPENCL = symbol("opencl")
const PETSCFECOMPOSITE = symbol("composite")
const MATSEQUSFFT = symbol("sequsfft")
const PETSCLIMITERSIN = symbol("sin")
const PETSCLIMITERZERO = symbol("zero")
const PETSCLIMITERNONE = symbol("none")
const PETSCLIMITERMINMOD = symbol("minmod")
const PETSCLIMITERVANLEER = symbol("vanleer")
const PETSCLIMITERVANALBADA = symbol("vanalbada")
const PETSCLIMITERSUPERBEE = symbol("superbee")
const PETSCLIMITERMC = symbol("mc")
const PETSCFVUPWIND = symbol("upwind")
const PETSCFVLEASTSQUARES = symbol("leastsquares")
const PETSCPARTITIONERCHACO = symbol("chaco")
const PETSCPARTITIONERPARMETIS = symbol("parmetis")
const PETSCPARTITIONERSHELL = symbol("shell")
const PETSCPARTITIONERSIMPLE = symbol("simple")
const PETSCDSBASIC = symbol("basic")
const CHARACTERISTICDA = symbol("da")
const PCNONE = symbol("none")
const PCJACOBI = symbol("jacobi")
const PCSOR = symbol("sor")
const PCLU = symbol("lu")
const PCSHELL = symbol("shell")
const PCBJACOBI = symbol("bjacobi")
const PCMG = symbol("mg")
const PCEISENSTAT = symbol("eisenstat")
const PCILU = symbol("ilu")
const PCICC = symbol("icc")
const PCASM = symbol("asm")
const PCGASM = symbol("gasm")
const PCKSP = symbol("ksp")
const PCCOMPOSITE = symbol("composite")
const PCREDUNDANT = symbol("redundant")
const PCSPAI = symbol("spai")
const PCNN = symbol("nn")
const PCCHOLESKY = symbol("cholesky")
const PCPBJACOBI = symbol("pbjacobi")
const PCMAT = symbol("mat")
const PCHYPRE = symbol("hypre")
const PCPARMS = symbol("parms")
const PCFIELDSPLIT = symbol("fieldsplit")
const PCTFS = symbol("tfs")
const PCML = symbol("ml")
const PCGALERKIN = symbol("galerkin")
const PCEXOTIC = symbol("exotic")
const PCCP = symbol("cp")
const PCBFBT = symbol("bfbt")
const PCLSC = symbol("lsc")
const PCPYTHON = symbol("python")
const PCPFMG = symbol("pfmg")
const PCSYSPFMG = symbol("syspfmg")
const PCREDISTRIBUTE = symbol("redistribute")
const PCSVD = symbol("svd")
const PCGAMG = symbol("gamg")
const PCSACUSP = symbol("sacusp")
const PCSACUSPPOLY = symbol("sacusppoly")
const PCBICGSTABCUSP = symbol("bicgstabcusp")
const PCAINVCUSP = symbol("ainvcusp")
const PCBDDC = symbol("bddc")
const PCKACZMARZ = symbol("kaczmarz")

#= # begin enum PCSide =#
typealias PCSide Cint

const PC_SIDE_DEFAULT = (Int32)(-1)
const PC_LEFT = (Int32)(0)
const PC_RIGHT = (Int32)(1)
const PC_SYMMETRIC = (Int32)(2)

#= # end enum PCSide =#
const PC_SIDE_MAX = PC_SYMMETRIC + 1
const PCGAMGAGG = symbol("agg")
const PCGAMGGEO = symbol("geo")
const PCGAMGCLASSICAL = symbol("classical")
const PCGAMGCLASSICALDIRECT = symbol("direct")
const PCGAMGCLASSICALSTANDARD = symbol("standard")

#= # begin enum PCMGType =#
typealias PCMGType UInt32

const PC_MG_MULTIPLICATIVE = (UInt32)(0)
const PC_MG_ADDITIVE = (UInt32)(1)
const PC_MG_FULL = (UInt32)(2)
const PC_MG_KASKADE = (UInt32)(3)

#= # end enum PCMGType =#
const PC_MG_CASCADE = PC_MG_KASKADE
const PC_FILE_CLASSID = 1211222
const KSPRICHARDSON = symbol("richardson")
const KSPCHEBYSHEV = symbol("chebyshev")
const KSPCG = symbol("cg")
const KSPGROPPCG = symbol("groppcg")
const KSPPIPECG = symbol("pipecg")
const KSPCGNE = symbol("cgne")
const KSPNASH = symbol("nash")
const KSPSTCG = symbol("stcg")
const KSPGLTR = symbol("gltr")
const KSPFCG = symbol("fcg")
const KSPGMRES = symbol("gmres")
const KSPFGMRES = symbol("fgmres")
const KSPLGMRES = symbol("lgmres")
const KSPDGMRES = symbol("dgmres")
const KSPPGMRES = symbol("pgmres")
const KSPTCQMR = symbol("tcqmr")
const KSPBCGS = symbol("bcgs")
const KSPIBCGS = symbol("ibcgs")
const KSPFBCGS = symbol("fbcgs")
const KSPFBCGSR = symbol("fbcgsr")
const KSPBCGSL = symbol("bcgsl")
const KSPCGS = symbol("cgs")
const KSPTFQMR = symbol("tfqmr")
const KSPCR = symbol("cr")
const KSPPIPECR = symbol("pipecr")
const KSPLSQR = symbol("lsqr")
const KSPPREONLY = symbol("preonly")
const KSPQCG = symbol("qcg")
const KSPBICG = symbol("bicg")
const KSPMINRES = symbol("minres")
const KSPSYMMLQ = symbol("symmlq")
const KSPLCD = symbol("lcd")
const KSPPYTHON = symbol("python")
const KSPGCR = symbol("gcr")
const KSP_FILE_CLASSID = 1211223

#= # begin enum KSPNormType =#
typealias KSPNormType Cint

const KSP_NORM_DEFAULT = (Int32)(-1)
const KSP_NORM_NONE = (Int32)(0)
const KSP_NORM_PRECONDITIONED = (Int32)(1)
const KSP_NORM_UNPRECONDITIONED = (Int32)(2)
const KSP_NORM_NATURAL = (Int32)(3)

#= # end enum KSPNormType =#
const KSP_NORM_MAX = KSP_NORM_NATURAL + 1

#= # Skipping MacroDefinition: KSPDefaultConverged ( KSPDefaultConverged , KSPConvergedDefault ) =#
#= # Skipping MacroDefinition: KSPDefaultConvergedDestroy ( KSPDefaultConvergedDestroy , KSPConvergedDefaultDestroy ) =#
#= # Skipping MacroDefinition: KSPDefaultConvergedCreate ( KSPDefaultConvergedCreate , KSPConvergedDefaultCreate ) =#
#= # Skipping MacroDefinition: KSPDefaultConvergedSetUIRNorm ( KSPDefaultConvergedSetUIRNorm , KSPConvergedDefaultSetUIRNorm ) =#
#= # Skipping MacroDefinition: KSPDefaultConvergedSetUMIRNorm ( KSPDefaultConvergedSetUMIRNorm , KSPConvergedDefaultSetUMIRNorm ) =#
#= # Skipping MacroDefinition: KSPSkipConverged ( KSPSkipConverged , KSPConvergedSkip ) =#
const SNESNEWTONLS = symbol("newtonls")
const SNESNEWTONTR = symbol("newtontr")
const SNESPYTHON = symbol("python")
const SNESTEST = symbol("test")
const SNESNRICHARDSON = symbol("nrichardson")
const SNESKSPONLY = symbol("ksponly")
const SNESVINEWTONRSLS = symbol("vinewtonrsls")
const SNESVINEWTONSSLS = symbol("vinewtonssls")
const SNESNGMRES = symbol("ngmres")
const SNESQN = symbol("qn")
const SNESSHELL = symbol("shell")
const SNESNGS = symbol("ngs")
const SNESNCG = symbol("ncg")
const SNESFAS = symbol("fas")
const SNESMS = symbol("ms")
const SNESNASM = symbol("nasm")
const SNESANDERSON = symbol("anderson")
const SNESASPIN = symbol("aspin")
const SNESCOMPOSITE = symbol("composite")
const SNES_FILE_CLASSID = 1211224

#= # Skipping MacroDefinition: SNESSkipConverged ( SNESSkipConverged , SNESConvergedSkip ) =#
const SNESLINESEARCHBT = symbol("bt")
const SNESLINESEARCHNLEQERR = symbol("nleqerr")
const SNESLINESEARCHBASIC = symbol("basic")
const SNESLINESEARCHL2 = symbol("l2")
const SNESLINESEARCHCP = symbol("cp")
const SNESLINESEARCHSHELL = symbol("shell")
const SNES_LINESEARCH_ORDER_LINEAR = 1
const SNES_LINESEARCH_ORDER_QUADRATIC = 2
const SNES_LINESEARCH_ORDER_CUBIC = 3
const SNESMSM62 = symbol("m62")
const SNESMSEULER = symbol("euler")
const SNESMSJAMESON83 = symbol("jameson83")
const SNESMSVLTP21 = symbol("vltp21")
const SNESMSVLTP31 = symbol("vltp31")
const SNESMSVLTP41 = symbol("vltp41")
const SNESMSVLTP51 = symbol("vltp51")
const SNESMSVLTP61 = symbol("vltp61")
const TSEULER = symbol("euler")
const TSBEULER = symbol("beuler")
const TSPSEUDO = symbol("pseudo")
const TSCN = symbol("cn")
const TSSUNDIALS = symbol("sundials")
const TSRK = symbol("rk")
const TSPYTHON = symbol("python")
const TSTHETA = symbol("theta")
const TSALPHA = symbol("alpha")
const TSGL = symbol("gl")
const TSSSP = symbol("ssp")
const TSARKIMEX = symbol("arkimex")
const TSROSW = symbol("rosw")
const TSEIMEX = symbol("eimex")
const TSMIMEX = symbol("mimex")
const TSTRAJECTORYBASIC = symbol("basic")
const TSTRAJECTORYSINGLEFILE = symbol("singlefile")
const TS_FILE_CLASSID = 1211225
const TSSSPRKS2 = symbol("rks2")
const TSSSPRKS3 = symbol("rks3")
const TSSSPRK104 = symbol("rk104")
const TSADAPTBASIC = symbol("basic")
const TSADAPTNONE = symbol("none")
const TSADAPTCFL = symbol("cfl")
const TSGLADAPT_NONE = symbol("none")
const TSGLADAPT_SIZE = symbol("size")
const TSGLADAPT_BOTH = symbol("both")
const TSGLACCEPT_ALWAYS = symbol("always")
const TSGL_IRKS = symbol("irks")

#= # Skipping MacroDefinition: TSEIMEXType char * =#
const TSRK1FE = symbol("1fe")
const TSRK2A = symbol("2a")
const TSRK3 = symbol("3")
const TSRK3BS = symbol("3bs")
const TSRK4 = symbol("4")
const TSRK5F = symbol("5f")
const TSRK5DP = symbol("5dp")
const TSARKIMEX1BEE = symbol("1bee")
const TSARKIMEXA2 = symbol("a2")
const TSARKIMEXL2 = symbol("l2")
const TSARKIMEXARS122 = symbol("ars122")
const TSARKIMEX2C = symbol("2c")
const TSARKIMEX2D = symbol("2d")
const TSARKIMEX2E = symbol("2e")
const TSARKIMEXPRSSP2 = symbol("prssp2")
const TSARKIMEX3 = symbol("3")
const TSARKIMEXBPR3 = symbol("bpr3")
const TSARKIMEXARS443 = symbol("ars443")
const TSARKIMEX4 = symbol("4")
const TSARKIMEX5 = symbol("5")
const TSROSW2M = symbol("2m")
const TSROSW2P = symbol("2p")
const TSROSWRA3PW = symbol("ra3pw")
const TSROSWRA34PW2 = symbol("ra34pw2")
const TSROSWRODAS3 = symbol("rodas3")
const TSROSWSANDU3 = symbol("sandu3")
const TSROSWASSP3P3S1C = symbol("assp3p3s1c")
const TSROSWLASSP3P4S2C = symbol("lassp3p4s2c")
const TSROSWLLSSP3P4S2C = symbol("llssp3p4s2c")
const TSROSWARK3 = symbol("ark3")
const TSROSWTHETA1 = symbol("theta1")
const TSROSWTHETA2 = symbol("theta2")
const TSROSWGRK4T = symbol("grk4t")
const TSROSWSHAMP4 = symbol("shamp4")
const TSROSWVELDD4 = symbol("veldd4")
const TSROSW4L = symbol("4l")

#= # Skipping MacroDefinition: TaoType char * =#
const TAOLMVM = symbol("lmvm")
const TAONLS = symbol("nls")
const TAONTR = symbol("ntr")
const TAONTL = symbol("ntl")
const TAOCG = symbol("cg")
const TAOTRON = symbol("tron")
const TAOOWLQN = symbol("owlqn")
const TAOBMRM = symbol("bmrm")
const TAOBLMVM = symbol("blmvm")
const TAOBQPIP = symbol("bqpip")
const TAOGPCG = symbol("gpcg")
const TAONM = symbol("nm")
const TAOPOUNDERS = symbol("pounders")
const TAOLCL = symbol("lcl")
const TAOSSILS = symbol("ssils")
const TAOSSFLS = symbol("ssfls")
const TAOASILS = symbol("asils")
const TAOASFLS = symbol("asfls")
const TAOIPM = symbol("ipm")
const TAOTEST = symbol("test")

#= # Skipping MacroDefinition: TaoLineSearchType char * =#
const TAOLINESEARCHUNIT = symbol("unit")
const TAOLINESEARCHMT = symbol("more-thuente")
const TAOLINESEARCHGPCG = symbol("gpcg")
const TAOLINESEARCHARMIJO = symbol("armijo")
const TAOLINESEARCHOWARMIJO = symbol("owarmijo")
const TAOLINESEARCHIPM = symbol("ipm")

typealias PetscErrorCode Cint
typealias PetscClassId Cint
typealias PetscMPIInt Cint

#= # begin enum ANONYMOUS_1 =#
typealias ANONYMOUS_1 UInt32

const ENUM_DUMMY = (UInt32)(0)

#= # end enum ANONYMOUS_1 =#
#= # begin enum PetscEnum =#
typealias PetscEnum UInt32

const ENUM_DUMMY = (UInt32)(0)

#= # end enum PetscEnum =#
typealias Petsc64bitInt Cint

# excluding lhs of  typealias PetscInt Cint
typealias PetscBLASInt Cint

#= # begin enum ANONYMOUS_2 =#
typealias ANONYMOUS_2 UInt32

const PETSC_PRECISION_SINGLE = (UInt32)(4)
const PETSC_PRECISION_DOUBLE = (UInt32)(8)

#= # end enum ANONYMOUS_2 =#
#= # begin enum PetscPrecision =#
typealias PetscPrecision UInt32

const PETSC_PRECISION_SINGLE = (UInt32)(4)
const PETSC_PRECISION_DOUBLE = (UInt32)(8)

#= # end enum PetscPrecision =#
#= # begin enum ANONYMOUS_3 =#
typealias ANONYMOUS_3 UInt32

const PETSC_FALSE = (UInt32)(0)
const PETSC_TRUE = (UInt32)(1)

#= # end enum ANONYMOUS_3 =#
#= # begin enum PetscBool =#
typealias PetscBool UInt32

const PETSC_FALSE = (UInt32)(0)
const PETSC_TRUE = (UInt32)(1)

#= # end enum PetscBool =#
typealias MatReal Float32

#= skipping type declaration with undefined symbols:
immutable petsc_mpiu_2scalar
    a::PetscScalar
    b::PetscScalar
end 
=#
typealias PetscLogDouble Cdouble

#= # begin enum ANONYMOUS_4 =#
typealias ANONYMOUS_4 UInt32

const PETSC_INT = (UInt32)(0)
const PETSC_DOUBLE = (UInt32)(1)
const PETSC_COMPLEX = (UInt32)(2)
const PETSC_LONG = (UInt32)(3)
const PETSC_SHORT = (UInt32)(4)
const PETSC_FLOAT = (UInt32)(5)
const PETSC_CHAR = (UInt32)(6)
const PETSC_BIT_LOGICAL = (UInt32)(7)
const PETSC_ENUM = (UInt32)(8)
const PETSC_BOOL = (UInt32)(9)
const PETSC___FLOAT128 = (UInt32)(10)
const PETSC_OBJECT = (UInt32)(11)
const PETSC_FUNCTION = (UInt32)(12)
const PETSC_STRING = (UInt32)(12)

#= # end enum ANONYMOUS_4 =#
# skipping undefined typealias typealias PetscToken Ptr{_p_PetscToken}
# skipping undefined typealias typealias PetscObject Ptr{_p_PetscObject}
typealias PetscObjectId Petsc64bitInt
typealias PetscObjectState Petsc64bitInt

# skipping undefined typealias typealias PetscFunctionList Ptr{_n_PetscFunctionList}
#= # begin enum ANONYMOUS_5 =#
typealias ANONYMOUS_5 UInt32

const FILE_MODE_READ = (UInt32)(0)
const FILE_MODE_WRITE = (UInt32)(1)
const FILE_MODE_APPEND = (UInt32)(2)
const FILE_MODE_UPDATE = (UInt32)(3)
const FILE_MODE_APPEND_UPDATE = (UInt32)(4)

#= # end enum ANONYMOUS_5 =#
#= # begin enum PetscFileMode =#
typealias PetscFileMode UInt32

const FILE_MODE_READ = (UInt32)(0)
const FILE_MODE_WRITE = (UInt32)(1)
const FILE_MODE_APPEND = (UInt32)(2)
const FILE_MODE_UPDATE = (UInt32)(3)
const FILE_MODE_APPEND_UPDATE = (UInt32)(4)

#= # end enum PetscFileMode =#
#= # begin enum ANONYMOUS_6 =#
typealias ANONYMOUS_6 UInt32

const PETSC_ERROR_INITIAL = (UInt32)(0)
const PETSC_ERROR_REPEAT = (UInt32)(1)
const PETSC_ERROR_IN_CXX = (UInt32)(2)

#= # end enum ANONYMOUS_6 =#
#= # begin enum PetscErrorType =#
typealias PetscErrorType UInt32

const PETSC_ERROR_INITIAL = (UInt32)(0)
const PETSC_ERROR_REPEAT = (UInt32)(1)
const PETSC_ERROR_IN_CXX = (UInt32)(2)

#= # end enum PetscErrorType =#
#= # begin enum ANONYMOUS_7 =#
typealias ANONYMOUS_7 UInt32

const PETSC_FP_TRAP_OFF = (UInt32)(0)
const PETSC_FP_TRAP_ON = (UInt32)(1)

#= # end enum ANONYMOUS_7 =#
#= # begin enum PetscFPTrap =#
typealias PetscFPTrap UInt32

const PETSC_FP_TRAP_OFF = (UInt32)(0)
const PETSC_FP_TRAP_ON = (UInt32)(1)

#= # end enum PetscFPTrap =#
immutable Array_64_Ptr
    d1::Ptr{UInt8}
    d2::Ptr{UInt8}
    d3::Ptr{UInt8}
    d4::Ptr{UInt8}
    d5::Ptr{UInt8}
    d6::Ptr{UInt8}
    d7::Ptr{UInt8}
    d8::Ptr{UInt8}
    d9::Ptr{UInt8}
    d10::Ptr{UInt8}
    d11::Ptr{UInt8}
    d12::Ptr{UInt8}
    d13::Ptr{UInt8}
    d14::Ptr{UInt8}
    d15::Ptr{UInt8}
    d16::Ptr{UInt8}
    d17::Ptr{UInt8}
    d18::Ptr{UInt8}
    d19::Ptr{UInt8}
    d20::Ptr{UInt8}
    d21::Ptr{UInt8}
    d22::Ptr{UInt8}
    d23::Ptr{UInt8}
    d24::Ptr{UInt8}
    d25::Ptr{UInt8}
    d26::Ptr{UInt8}
    d27::Ptr{UInt8}
    d28::Ptr{UInt8}
    d29::Ptr{UInt8}
    d30::Ptr{UInt8}
    d31::Ptr{UInt8}
    d32::Ptr{UInt8}
    d33::Ptr{UInt8}
    d34::Ptr{UInt8}
    d35::Ptr{UInt8}
    d36::Ptr{UInt8}
    d37::Ptr{UInt8}
    d38::Ptr{UInt8}
    d39::Ptr{UInt8}
    d40::Ptr{UInt8}
    d41::Ptr{UInt8}
    d42::Ptr{UInt8}
    d43::Ptr{UInt8}
    d44::Ptr{UInt8}
    d45::Ptr{UInt8}
    d46::Ptr{UInt8}
    d47::Ptr{UInt8}
    d48::Ptr{UInt8}
    d49::Ptr{UInt8}
    d50::Ptr{UInt8}
    d51::Ptr{UInt8}
    d52::Ptr{UInt8}
    d53::Ptr{UInt8}
    d54::Ptr{UInt8}
    d55::Ptr{UInt8}
    d56::Ptr{UInt8}
    d57::Ptr{UInt8}
    d58::Ptr{UInt8}
    d59::Ptr{UInt8}
    d60::Ptr{UInt8}
    d61::Ptr{UInt8}
    d62::Ptr{UInt8}
    d63::Ptr{UInt8}
    d64::Ptr{UInt8}
end

zero(::Type{Array_64_Ptr}) = begin  # /home/jared/.julia/v0.4/Clang_orig/src/wrap_c.jl, line 270:
        Array_64_Ptr(fill(zero(Ptr{UInt8}),64)...)
    end

immutable Array_64_Cint
    d1::Cint
    d2::Cint
    d3::Cint
    d4::Cint
    d5::Cint
    d6::Cint
    d7::Cint
    d8::Cint
    d9::Cint
    d10::Cint
    d11::Cint
    d12::Cint
    d13::Cint
    d14::Cint
    d15::Cint
    d16::Cint
    d17::Cint
    d18::Cint
    d19::Cint
    d20::Cint
    d21::Cint
    d22::Cint
    d23::Cint
    d24::Cint
    d25::Cint
    d26::Cint
    d27::Cint
    d28::Cint
    d29::Cint
    d30::Cint
    d31::Cint
    d32::Cint
    d33::Cint
    d34::Cint
    d35::Cint
    d36::Cint
    d37::Cint
    d38::Cint
    d39::Cint
    d40::Cint
    d41::Cint
    d42::Cint
    d43::Cint
    d44::Cint
    d45::Cint
    d46::Cint
    d47::Cint
    d48::Cint
    d49::Cint
    d50::Cint
    d51::Cint
    d52::Cint
    d53::Cint
    d54::Cint
    d55::Cint
    d56::Cint
    d57::Cint
    d58::Cint
    d59::Cint
    d60::Cint
    d61::Cint
    d62::Cint
    d63::Cint
    d64::Cint
end

zero(::Type{Array_64_Cint}) = begin  # /home/jared/.julia/v0.4/Clang_orig/src/wrap_c.jl, line 270:
        Array_64_Cint(fill(zero(Cint),64)...)
    end

immutable Array_64_PetscBool
    d1::PetscBool
    d2::PetscBool
    d3::PetscBool
    d4::PetscBool
    d5::PetscBool
    d6::PetscBool
    d7::PetscBool
    d8::PetscBool
    d9::PetscBool
    d10::PetscBool
    d11::PetscBool
    d12::PetscBool
    d13::PetscBool
    d14::PetscBool
    d15::PetscBool
    d16::PetscBool
    d17::PetscBool
    d18::PetscBool
    d19::PetscBool
    d20::PetscBool
    d21::PetscBool
    d22::PetscBool
    d23::PetscBool
    d24::PetscBool
    d25::PetscBool
    d26::PetscBool
    d27::PetscBool
    d28::PetscBool
    d29::PetscBool
    d30::PetscBool
    d31::PetscBool
    d32::PetscBool
    d33::PetscBool
    d34::PetscBool
    d35::PetscBool
    d36::PetscBool
    d37::PetscBool
    d38::PetscBool
    d39::PetscBool
    d40::PetscBool
    d41::PetscBool
    d42::PetscBool
    d43::PetscBool
    d44::PetscBool
    d45::PetscBool
    d46::PetscBool
    d47::PetscBool
    d48::PetscBool
    d49::PetscBool
    d50::PetscBool
    d51::PetscBool
    d52::PetscBool
    d53::PetscBool
    d54::PetscBool
    d55::PetscBool
    d56::PetscBool
    d57::PetscBool
    d58::PetscBool
    d59::PetscBool
    d60::PetscBool
    d61::PetscBool
    d62::PetscBool
    d63::PetscBool
    d64::PetscBool
end

zero(::Type{Array_64_PetscBool}) = begin  # /home/jared/.julia/v0.4/Clang_orig/src/wrap_c.jl, line 270:
        Array_64_PetscBool(fill(zero(PetscBool),64)...)
    end

immutable PetscStack
    _function::Array_64_Ptr
    file::Array_64_Ptr
    line::Array_64_Cint
    petscroutine::Array_64_PetscBool
    currentsize::Cint
    hotdepth::Cint
end

typealias PetscVoidStarFunction Ptr{Ptr{Void}}
typealias PetscVoidFunction Ptr{Void}
typealias PetscErrorCodeFunction Ptr{Void}

immutable PetscViewer{T}
    pobj::Ptr{Void}
end

#= # begin enum ANONYMOUS_8 =#
typealias ANONYMOUS_8 UInt32

const OPTION_INT = (UInt32)(0)
const OPTION_BOOL = (UInt32)(1)
const OPTION_REAL = (UInt32)(2)
const OPTION_FLIST = (UInt32)(3)
const OPTION_STRING = (UInt32)(4)
const OPTION_REAL_ARRAY = (UInt32)(5)
const OPTION_SCALAR_ARRAY = (UInt32)(6)
const OPTION_HEAD = (UInt32)(7)
const OPTION_INT_ARRAY = (UInt32)(8)
const OPTION_ELIST = (UInt32)(9)
const OPTION_BOOL_ARRAY = (UInt32)(10)
const OPTION_STRING_ARRAY = (UInt32)(11)

#= # end enum ANONYMOUS_8 =#
#= # begin enum PetscOptionType =#
typealias PetscOptionType UInt32

const OPTION_INT = (UInt32)(0)
const OPTION_BOOL = (UInt32)(1)
const OPTION_REAL = (UInt32)(2)
const OPTION_FLIST = (UInt32)(3)
const OPTION_STRING = (UInt32)(4)
const OPTION_REAL_ARRAY = (UInt32)(5)
const OPTION_SCALAR_ARRAY = (UInt32)(6)
const OPTION_HEAD = (UInt32)(7)
const OPTION_INT_ARRAY = (UInt32)(8)
const OPTION_ELIST = (UInt32)(9)
const OPTION_BOOL_ARRAY = (UInt32)(10)
const OPTION_STRING_ARRAY = (UInt32)(11)

#= # end enum PetscOptionType =#
immutable PetscOption{T}
    pobj::Ptr{Void}
end

#= skipping type declaration with undefined symbols:
immutable PetscOptions
    count::PetscInt
    next::PetscOption
    prefix::Ptr{UInt8}
    pprefix::Ptr{UInt8}
    title::Ptr{UInt8}
    comm::MPI_Comm
    printhelp::PetscBool
    changedmethod::PetscBool
    alreadyprinted::PetscBool
    object::PetscObject
end 
=#
typealias PetscDLHandle Ptr{Void}

#= # begin enum ANONYMOUS_9 =#
typealias ANONYMOUS_9 UInt32

const PETSC_DL_DECIDE = (UInt32)(0)
const PETSC_DL_NOW = (UInt32)(1)
const PETSC_DL_LOCAL = (UInt32)(2)

#= # end enum ANONYMOUS_9 =#
#= # begin enum PetscDLMode =#
typealias PetscDLMode UInt32

const PETSC_DL_DECIDE = (UInt32)(0)
const PETSC_DL_NOW = (UInt32)(1)
const PETSC_DL_LOCAL = (UInt32)(2)

#= # end enum PetscDLMode =#
# skipping undefined typealias typealias PetscObjectList Ptr{_n_PetscObjectList}
# skipping undefined typealias typealias PetscDLLibrary Ptr{_n_PetscDLLibrary}
typealias PetscLogEvent Cint
typealias PetscLogStage Cint

# skipping undefined typealias typealias PetscIntStack Ptr{_n_PetscIntStack}
immutable PetscClassRegInfo
    name::Ptr{UInt8}
    classid::PetscClassId
end

immutable PetscClassPerfInfo
    id::PetscClassId
    creations::Cint
    destructions::Cint
    mem::PetscLogDouble
    descMem::PetscLogDouble
end

# skipping undefined typealias typealias PetscClassRegLog Ptr{_n_PetscClassRegLog}
# skipping undefined typealias typealias PetscClassPerfLog Ptr{_n_PetscClassPerfLog}
immutable PetscEventRegInfo
    name::Ptr{UInt8}
    classid::PetscClassId
end

immutable PetscEventPerfInfo
    id::Cint
    active::PetscBool
    visible::PetscBool
    depth::Cint
    count::Cint
    flops::PetscLogDouble
    flops2::PetscLogDouble
    flopsTmp::PetscLogDouble
    time::PetscLogDouble
    time2::PetscLogDouble
    timeTmp::PetscLogDouble
    numMessages::PetscLogDouble
    messageLength::PetscLogDouble
    numReductions::PetscLogDouble
end

# skipping undefined typealias typealias PetscEventRegLog Ptr{_n_PetscEventRegLog}
# skipping undefined typealias typealias PetscEventPerfLog Ptr{_n_PetscEventPerfLog}
#= skipping type declaration with undefined symbols:
immutable PetscStageInfo
    name::Ptr{UInt8}
    used::PetscBool
    perfInfo::PetscEventPerfInfo
    eventLog::PetscEventPerfLog
    classLog::PetscClassPerfLog
end 
=#
# skipping undefined typealias typealias PetscStageLog Ptr{_n_PetscStageLog}
# skipping undefined typealias typealias PetscContainer Ptr{_p_PetscContainer}
typealias PetscRandomType Symbol

# skipping undefined typealias typealias PetscRandom Ptr{_p_PetscRandom}
#= # begin enum ANONYMOUS_10 =#
typealias ANONYMOUS_10 UInt32

const PETSC_BINARY_SEEK_SET = (UInt32)(0)
const PETSC_BINARY_SEEK_CUR = (UInt32)(1)
const PETSC_BINARY_SEEK_END = (UInt32)(2)

#= # end enum ANONYMOUS_10 =#
#= # begin enum PetscBinarySeekType =#
typealias PetscBinarySeekType UInt32

const PETSC_BINARY_SEEK_SET = (UInt32)(0)
const PETSC_BINARY_SEEK_CUR = (UInt32)(1)
const PETSC_BINARY_SEEK_END = (UInt32)(2)

#= # end enum PetscBinarySeekType =#
#= # begin enum ANONYMOUS_11 =#
typealias ANONYMOUS_11 Cint

const PETSC_BUILDTWOSIDED_NOTSET = (Int32)(-1)
const PETSC_BUILDTWOSIDED_ALLREDUCE = (Int32)(0)
const PETSC_BUILDTWOSIDED_IBARRIER = (Int32)(1)

#= # end enum ANONYMOUS_11 =#
#= # begin enum PetscBuildTwoSidedType =#
typealias PetscBuildTwoSidedType Cint

const PETSC_BUILDTWOSIDED_NOTSET = (Int32)(-1)
const PETSC_BUILDTWOSIDED_ALLREDUCE = (Int32)(0)
const PETSC_BUILDTWOSIDED_IBARRIER = (Int32)(1)

#= # end enum PetscBuildTwoSidedType =#
#= # begin enum ANONYMOUS_12 =#
typealias ANONYMOUS_12 UInt32

const NOT_SET_VALUES = (UInt32)(0)
const INSERT_VALUES = (UInt32)(1)
const ADD_VALUES = (UInt32)(2)
const MAX_VALUES = (UInt32)(3)
const INSERT_ALL_VALUES = (UInt32)(4)
const ADD_ALL_VALUES = (UInt32)(5)
const INSERT_BC_VALUES = (UInt32)(6)
const ADD_BC_VALUES = (UInt32)(7)

#= # end enum ANONYMOUS_12 =#
#= # begin enum InsertMode =#
typealias InsertMode UInt32

const NOT_SET_VALUES = (UInt32)(0)
const INSERT_VALUES = (UInt32)(1)
const ADD_VALUES = (UInt32)(2)
const MAX_VALUES = (UInt32)(3)
const INSERT_ALL_VALUES = (UInt32)(4)
const ADD_ALL_VALUES = (UInt32)(5)
const INSERT_BC_VALUES = (UInt32)(6)
const ADD_BC_VALUES = (UInt32)(7)

#= # end enum InsertMode =#
#= # begin enum ANONYMOUS_13 =#
typealias ANONYMOUS_13 UInt32

const PETSC_SUBCOMM_GENERAL = (UInt32)(0)
const PETSC_SUBCOMM_CONTIGUOUS = (UInt32)(1)
const PETSC_SUBCOMM_INTERLACED = (UInt32)(2)

#= # end enum ANONYMOUS_13 =#
#= # begin enum PetscSubcommType =#
typealias PetscSubcommType UInt32

const PETSC_SUBCOMM_GENERAL = (UInt32)(0)
const PETSC_SUBCOMM_CONTIGUOUS = (UInt32)(1)
const PETSC_SUBCOMM_INTERLACED = (UInt32)(2)

#= # end enum PetscSubcommType =#
# skipping undefined typealias typealias PetscSubcomm Ptr{_n_PetscSubcomm}
# skipping undefined typealias typealias PetscSegBuffer Ptr{_n_PetscSegBuffer}
# skipping undefined typealias typealias PetscBag Ptr{_n_PetscBag}
# skipping undefined typealias typealias PetscBagItem Ptr{_n_PetscBagItem}
typealias PetscViewerType Symbol
typealias PetscDrawType Symbol

# skipping undefined typealias typealias PetscDraw Ptr{_p_PetscDraw}
# skipping undefined typealias typealias PetscDrawAxis Ptr{_p_PetscDrawAxis}
# skipping undefined typealias typealias PetscDrawLG Ptr{_p_PetscDrawLG}
# skipping undefined typealias typealias PetscDrawSP Ptr{_p_PetscDrawSP}
# skipping undefined typealias typealias PetscDrawHG Ptr{_p_PetscDrawHG}
# skipping undefined typealias typealias PetscDrawBar Ptr{_p_PetscDrawBar}
#= # begin enum ANONYMOUS_14 =#
typealias ANONYMOUS_14 UInt32

const PETSC_VIEWER_DEFAULT = (UInt32)(0)
const PETSC_VIEWER_ASCII_MATLAB = (UInt32)(1)
const PETSC_VIEWER_ASCII_MATHEMATICA = (UInt32)(2)
const PETSC_VIEWER_ASCII_IMPL = (UInt32)(3)
const PETSC_VIEWER_ASCII_INFO = (UInt32)(4)
const PETSC_VIEWER_ASCII_INFO_DETAIL = (UInt32)(5)
const PETSC_VIEWER_ASCII_COMMON = (UInt32)(6)
const PETSC_VIEWER_ASCII_SYMMODU = (UInt32)(7)
const PETSC_VIEWER_ASCII_INDEX = (UInt32)(8)
const PETSC_VIEWER_ASCII_DENSE = (UInt32)(9)
const PETSC_VIEWER_ASCII_MATRIXMARKET = (UInt32)(10)
const PETSC_VIEWER_ASCII_VTK = (UInt32)(11)
const PETSC_VIEWER_ASCII_VTK_CELL = (UInt32)(12)
const PETSC_VIEWER_ASCII_VTK_COORDS = (UInt32)(13)
const PETSC_VIEWER_ASCII_PCICE = (UInt32)(14)
const PETSC_VIEWER_ASCII_PYTHON = (UInt32)(15)
const PETSC_VIEWER_ASCII_FACTOR_INFO = (UInt32)(16)
const PETSC_VIEWER_ASCII_LATEX = (UInt32)(17)
const PETSC_VIEWER_DRAW_BASIC = (UInt32)(18)
const PETSC_VIEWER_DRAW_LG = (UInt32)(19)
const PETSC_VIEWER_DRAW_CONTOUR = (UInt32)(20)
const PETSC_VIEWER_DRAW_PORTS = (UInt32)(21)
const PETSC_VIEWER_VTK_VTS = (UInt32)(22)
const PETSC_VIEWER_VTK_VTR = (UInt32)(23)
const PETSC_VIEWER_VTK_VTU = (UInt32)(24)
const PETSC_VIEWER_BINARY_MATLAB = (UInt32)(25)
const PETSC_VIEWER_NATIVE = (UInt32)(26)
const PETSC_VIEWER_HDF5_VIZ = (UInt32)(27)
const PETSC_VIEWER_NOFORMAT = (UInt32)(28)

#= # end enum ANONYMOUS_14 =#
#= # begin enum PetscViewerFormat =#
typealias PetscViewerFormat UInt32

const PETSC_VIEWER_DEFAULT = (UInt32)(0)
const PETSC_VIEWER_ASCII_MATLAB = (UInt32)(1)
const PETSC_VIEWER_ASCII_MATHEMATICA = (UInt32)(2)
const PETSC_VIEWER_ASCII_IMPL = (UInt32)(3)
const PETSC_VIEWER_ASCII_INFO = (UInt32)(4)
const PETSC_VIEWER_ASCII_INFO_DETAIL = (UInt32)(5)
const PETSC_VIEWER_ASCII_COMMON = (UInt32)(6)
const PETSC_VIEWER_ASCII_SYMMODU = (UInt32)(7)
const PETSC_VIEWER_ASCII_INDEX = (UInt32)(8)
const PETSC_VIEWER_ASCII_DENSE = (UInt32)(9)
const PETSC_VIEWER_ASCII_MATRIXMARKET = (UInt32)(10)
const PETSC_VIEWER_ASCII_VTK = (UInt32)(11)
const PETSC_VIEWER_ASCII_VTK_CELL = (UInt32)(12)
const PETSC_VIEWER_ASCII_VTK_COORDS = (UInt32)(13)
const PETSC_VIEWER_ASCII_PCICE = (UInt32)(14)
const PETSC_VIEWER_ASCII_PYTHON = (UInt32)(15)
const PETSC_VIEWER_ASCII_FACTOR_INFO = (UInt32)(16)
const PETSC_VIEWER_ASCII_LATEX = (UInt32)(17)
const PETSC_VIEWER_DRAW_BASIC = (UInt32)(18)
const PETSC_VIEWER_DRAW_LG = (UInt32)(19)
const PETSC_VIEWER_DRAW_CONTOUR = (UInt32)(20)
const PETSC_VIEWER_DRAW_PORTS = (UInt32)(21)
const PETSC_VIEWER_VTK_VTS = (UInt32)(22)
const PETSC_VIEWER_VTK_VTR = (UInt32)(23)
const PETSC_VIEWER_VTK_VTU = (UInt32)(24)
const PETSC_VIEWER_BINARY_MATLAB = (UInt32)(25)
const PETSC_VIEWER_NATIVE = (UInt32)(26)
const PETSC_VIEWER_HDF5_VIZ = (UInt32)(27)
const PETSC_VIEWER_NOFORMAT = (UInt32)(28)

#= # end enum PetscViewerFormat =#
#= # begin enum ANONYMOUS_15 =#
typealias ANONYMOUS_15 UInt32

const PETSC_VTK_POINT_FIELD = (UInt32)(0)
const PETSC_VTK_POINT_VECTOR_FIELD = (UInt32)(1)
const PETSC_VTK_CELL_FIELD = (UInt32)(2)
const PETSC_VTK_CELL_VECTOR_FIELD = (UInt32)(3)

#= # end enum ANONYMOUS_15 =#
#= # begin enum PetscViewerVTKFieldType =#
typealias PetscViewerVTKFieldType UInt32

const PETSC_VTK_POINT_FIELD = (UInt32)(0)
const PETSC_VTK_POINT_VECTOR_FIELD = (UInt32)(1)
const PETSC_VTK_CELL_FIELD = (UInt32)(2)
const PETSC_VTK_CELL_VECTOR_FIELD = (UInt32)(3)

#= # end enum PetscViewerVTKFieldType =#
# skipping undefined typealias typealias PetscViewers Ptr{_n_PetscViewers}
typealias PetscBT Symbol

# skipping undefined typealias typealias PetscTable Ptr{_n_PetscTable}
typealias PetscTablePosition Ptr{Int64}

# skipping undefined typealias typealias PetscMatlabEngine Ptr{_p_PetscMatlabEngine}
#= # begin enum ANONYMOUS_16 =#
typealias ANONYMOUS_16 UInt32

const PETSC_DRAW_MARKER_CROSS = (UInt32)(0)
const PETSC_DRAW_MARKER_POINT = (UInt32)(1)
const PETSC_DRAW_MARKER_PLUS = (UInt32)(2)
const PETSC_DRAW_MARKER_CIRCLE = (UInt32)(3)

#= # end enum ANONYMOUS_16 =#
#= # begin enum PetscDrawMarkerType =#
typealias PetscDrawMarkerType UInt32

const PETSC_DRAW_MARKER_CROSS = (UInt32)(0)
const PETSC_DRAW_MARKER_POINT = (UInt32)(1)
const PETSC_DRAW_MARKER_PLUS = (UInt32)(2)
const PETSC_DRAW_MARKER_CIRCLE = (UInt32)(3)

#= # end enum PetscDrawMarkerType =#
#= # begin enum ANONYMOUS_17 =#
typealias ANONYMOUS_17 UInt32

const PETSC_BUTTON_NONE = (UInt32)(0)
const PETSC_BUTTON_LEFT = (UInt32)(1)
const PETSC_BUTTON_CENTER = (UInt32)(2)
const PETSC_BUTTON_RIGHT = (UInt32)(3)
const PETSC_BUTTON_LEFT_SHIFT = (UInt32)(4)
const PETSC_BUTTON_CENTER_SHIFT = (UInt32)(5)
const PETSC_BUTTON_RIGHT_SHIFT = (UInt32)(6)

#= # end enum ANONYMOUS_17 =#
#= # begin enum PetscDrawButton =#
typealias PetscDrawButton UInt32

const PETSC_BUTTON_NONE = (UInt32)(0)
const PETSC_BUTTON_LEFT = (UInt32)(1)
const PETSC_BUTTON_CENTER = (UInt32)(2)
const PETSC_BUTTON_RIGHT = (UInt32)(3)
const PETSC_BUTTON_LEFT_SHIFT = (UInt32)(4)
const PETSC_BUTTON_CENTER_SHIFT = (UInt32)(5)
const PETSC_BUTTON_RIGHT_SHIFT = (UInt32)(6)

#= # end enum PetscDrawButton =#
#= skipping type declaration with undefined symbols:
immutable PetscDrawViewPorts
    nports::PetscInt
    xl::Ptr{PetscReal}
    xr::Ptr{PetscReal}
    yl::Ptr{PetscReal}
    yr::Ptr{PetscReal}
    draw::PetscDraw
    port_xl::PetscReal
    port_yl::PetscReal
    port_xr::PetscReal
    port_yr::PetscReal
end 
=#
# skipping undefined typealias typealias PetscSF Ptr{_p_PetscSF}
#= skipping type declaration with undefined symbols:
immutable PetscSFNode
    rank::PetscInt
    index::PetscInt
end 
=#
immutable IS{T}
    pobj::Ptr{Void}
end

immutable ISLocalToGlobalMapping{T}
    pobj::Ptr{Void}
end

immutable ISColoring{T}
    pobj::Ptr{Void}
end

immutable PetscLayout{T}
    pobj::Ptr{Void}
end

# skipping undefined typealias typealias PetscSection Ptr{_p_PetscSection}
typealias ISType Symbol

#= # begin enum ANONYMOUS_18 =#
typealias ANONYMOUS_18 UInt32

const IS_GTOLM_MASK = (UInt32)(0)
const IS_GTOLM_DROP = (UInt32)(1)

#= # end enum ANONYMOUS_18 =#
#= # begin enum ISGlobalToLocalMappingType =#
typealias ISGlobalToLocalMappingType UInt32

const IS_GTOLM_MASK = (UInt32)(0)
const IS_GTOLM_DROP = (UInt32)(1)

#= # end enum ISGlobalToLocalMappingType =#
#= # begin enum ANONYMOUS_19 =#
typealias ANONYMOUS_19 UInt32

const IS_COLORING_GLOBAL = (UInt32)(0)
const IS_COLORING_GHOSTED = (UInt32)(1)

#= # end enum ANONYMOUS_19 =#
#= # begin enum ISColoringType =#
typealias ISColoringType UInt32

const IS_COLORING_GLOBAL = (UInt32)(0)
const IS_COLORING_GHOSTED = (UInt32)(1)

#= # end enum ISColoringType =#
typealias PETSC_IS_COLOR_VALUE_TYPE UInt32

immutable Vec{T}
    pobj::Ptr{Void}
end

immutable VecScatter{T}
    pobj::Ptr{Void}
end

#= # begin enum ANONYMOUS_20 =#
typealias ANONYMOUS_20 UInt32

const SCATTER_FORWARD = (UInt32)(0)
const SCATTER_REVERSE = (UInt32)(1)
const SCATTER_FORWARD_LOCAL = (UInt32)(2)
const SCATTER_REVERSE_LOCAL = (UInt32)(3)
const SCATTER_LOCAL = (UInt32)(2)

#= # end enum ANONYMOUS_20 =#
#= # begin enum ScatterMode =#
typealias ScatterMode UInt32

const SCATTER_FORWARD = (UInt32)(0)
const SCATTER_REVERSE = (UInt32)(1)
const SCATTER_FORWARD_LOCAL = (UInt32)(2)
const SCATTER_REVERSE_LOCAL = (UInt32)(3)
const SCATTER_LOCAL = (UInt32)(2)

#= # end enum ScatterMode =#
typealias VecType Symbol

#= # begin enum ANONYMOUS_21 =#
typealias ANONYMOUS_21 UInt32

const NORM_1 = (UInt32)(0)
const NORM_2 = (UInt32)(1)
const NORM_FROBENIUS = (UInt32)(2)
const NORM_INFINITY = (UInt32)(3)
const NORM_1_AND_2 = (UInt32)(4)

#= # end enum ANONYMOUS_21 =#
#= # begin enum ANONYMOUS_22 =#
typealias ANONYMOUS_22 UInt32

const VEC_IGNORE_OFF_PROC_ENTRIES = (UInt32)(0)
const VEC_IGNORE_NEGATIVE_INDICES = (UInt32)(1)

#= # end enum ANONYMOUS_22 =#
#= # begin enum VecOption =#
typealias VecOption UInt32

const VEC_IGNORE_OFF_PROC_ENTRIES = (UInt32)(0)
const VEC_IGNORE_NEGATIVE_INDICES = (UInt32)(1)

#= # end enum VecOption =#
#= # begin enum ANONYMOUS_23 =#
typealias ANONYMOUS_23 UInt32

const VECOP_VIEW = (UInt32)(33)
const VECOP_LOAD = (UInt32)(41)
const VECOP_DUPLICATE = (UInt32)(0)

#= # end enum ANONYMOUS_23 =#
#= # begin enum VecOperation =#
typealias VecOperation UInt32

const VECOP_VIEW = (UInt32)(33)
const VECOP_LOAD = (UInt32)(41)
const VECOP_DUPLICATE = (UInt32)(0)

#= # end enum VecOperation =#
# skipping undefined typealias typealias Vecs Ptr{_n_Vecs}
immutable Mat{T}
    pobj::Ptr{Void}
end

typealias MatType Symbol

#= # begin enum ANONYMOUS_24 =#
typealias ANONYMOUS_24 UInt32

const MAT_FACTOR_NONE = (UInt32)(0)
const MAT_FACTOR_LU = (UInt32)(1)
const MAT_FACTOR_CHOLESKY = (UInt32)(2)
const MAT_FACTOR_ILU = (UInt32)(3)
const MAT_FACTOR_ICC = (UInt32)(4)
const MAT_FACTOR_ILUDT = (UInt32)(5)

#= # end enum ANONYMOUS_24 =#
#= # begin enum MatFactorType =#
typealias MatFactorType UInt32

const MAT_FACTOR_NONE = (UInt32)(0)
const MAT_FACTOR_LU = (UInt32)(1)
const MAT_FACTOR_CHOLESKY = (UInt32)(2)
const MAT_FACTOR_ILU = (UInt32)(3)
const MAT_FACTOR_ICC = (UInt32)(4)
const MAT_FACTOR_ILUDT = (UInt32)(5)

#= # end enum MatFactorType =#
#= # begin enum ANONYMOUS_25 =#
typealias ANONYMOUS_25 UInt32

const MAT_INITIAL_MATRIX = (UInt32)(0)
const MAT_REUSE_MATRIX = (UInt32)(1)
const MAT_IGNORE_MATRIX = (UInt32)(2)

#= # end enum ANONYMOUS_25 =#
#= # begin enum MatReuse =#
typealias MatReuse UInt32

const MAT_INITIAL_MATRIX = (UInt32)(0)
const MAT_REUSE_MATRIX = (UInt32)(1)
const MAT_IGNORE_MATRIX = (UInt32)(2)

#= # end enum MatReuse =#
#= # begin enum ANONYMOUS_26 =#
typealias ANONYMOUS_26 UInt32

const MAT_DO_NOT_GET_VALUES = (UInt32)(0)
const MAT_GET_VALUES = (UInt32)(1)

#= # end enum ANONYMOUS_26 =#
#= # begin enum MatGetSubMatrixOption =#
typealias MatGetSubMatrixOption UInt32

const MAT_DO_NOT_GET_VALUES = (UInt32)(0)
const MAT_GET_VALUES = (UInt32)(1)

#= # end enum MatGetSubMatrixOption =#
#= # begin enum ANONYMOUS_27 =#
typealias ANONYMOUS_27 UInt32

const DIFFERENT_NONZERO_PATTERN = (UInt32)(0)
const SUBSET_NONZERO_PATTERN = (UInt32)(1)
const SAME_NONZERO_PATTERN = (UInt32)(2)

#= # end enum ANONYMOUS_27 =#
#= # begin enum MatStructure =#
typealias MatStructure UInt32

const DIFFERENT_NONZERO_PATTERN = (UInt32)(0)
const SUBSET_NONZERO_PATTERN = (UInt32)(1)
const SAME_NONZERO_PATTERN = (UInt32)(2)

#= # end enum MatStructure =#
#= # begin enum ANONYMOUS_28 =#
typealias ANONYMOUS_28 UInt32

const MAT_COMPOSITE_ADDITIVE = (UInt32)(0)
const MAT_COMPOSITE_MULTIPLICATIVE = (UInt32)(1)

#= # end enum ANONYMOUS_28 =#
#= # begin enum MatCompositeType =#
typealias MatCompositeType UInt32

const MAT_COMPOSITE_ADDITIVE = (UInt32)(0)
const MAT_COMPOSITE_MULTIPLICATIVE = (UInt32)(1)

#= # end enum MatCompositeType =#
#= skipping type declaration with undefined symbols:
immutable MatStencil
    k::PetscInt
    j::PetscInt
    i::PetscInt
    c::PetscInt
end 
=#
#= # begin enum ANONYMOUS_29 =#
typealias ANONYMOUS_29 UInt32

const MAT_FLUSH_ASSEMBLY = (UInt32)(1)
const MAT_FINAL_ASSEMBLY = (UInt32)(0)

#= # end enum ANONYMOUS_29 =#
#= # begin enum MatAssemblyType =#
typealias MatAssemblyType UInt32

const MAT_FLUSH_ASSEMBLY = (UInt32)(1)
const MAT_FINAL_ASSEMBLY = (UInt32)(0)

#= # end enum MatAssemblyType =#
#= # begin enum ANONYMOUS_30 =#
typealias ANONYMOUS_30 Cint

const MAT_OPTION_MIN = (Int32)(-5)
const MAT_NEW_NONZERO_LOCATION_ERR = (Int32)(-4)
const MAT_UNUSED_NONZERO_LOCATION_ERR = (Int32)(-3)
const MAT_NEW_NONZERO_ALLOCATION_ERR = (Int32)(-2)
const MAT_ROW_ORIENTED = (Int32)(-1)
const MAT_SYMMETRIC = (Int32)(1)
const MAT_STRUCTURALLY_SYMMETRIC = (Int32)(2)
const MAT_NEW_DIAGONALS = (Int32)(3)
const MAT_IGNORE_OFF_PROC_ENTRIES = (Int32)(4)
const MAT_USE_HASH_TABLE = (Int32)(5)
const MAT_KEEP_NONZERO_PATTERN = (Int32)(6)
const MAT_IGNORE_ZERO_ENTRIES = (Int32)(7)
const MAT_USE_INODES = (Int32)(8)
const MAT_HERMITIAN = (Int32)(9)
const MAT_SYMMETRY_ETERNAL = (Int32)(10)
const MAT_DUMMY = (Int32)(11)
const MAT_IGNORE_LOWER_TRIANGULAR = (Int32)(12)
const MAT_ERROR_LOWER_TRIANGULAR = (Int32)(13)
const MAT_GETROW_UPPERTRIANGULAR = (Int32)(14)
const MAT_SPD = (Int32)(15)
const MAT_NO_OFF_PROC_ZERO_ROWS = (Int32)(16)
const MAT_NO_OFF_PROC_ENTRIES = (Int32)(17)
const MAT_NEW_NONZERO_LOCATIONS = (Int32)(18)
const MAT_OPTION_MAX = (Int32)(19)

#= # end enum ANONYMOUS_30 =#
#= # begin enum MatOption =#
typealias MatOption Cint

const MAT_OPTION_MIN = (Int32)(-5)
const MAT_NEW_NONZERO_LOCATION_ERR = (Int32)(-4)
const MAT_UNUSED_NONZERO_LOCATION_ERR = (Int32)(-3)
const MAT_NEW_NONZERO_ALLOCATION_ERR = (Int32)(-2)
const MAT_ROW_ORIENTED = (Int32)(-1)
const MAT_SYMMETRIC = (Int32)(1)
const MAT_STRUCTURALLY_SYMMETRIC = (Int32)(2)
const MAT_NEW_DIAGONALS = (Int32)(3)
const MAT_IGNORE_OFF_PROC_ENTRIES = (Int32)(4)
const MAT_USE_HASH_TABLE = (Int32)(5)
const MAT_KEEP_NONZERO_PATTERN = (Int32)(6)
const MAT_IGNORE_ZERO_ENTRIES = (Int32)(7)
const MAT_USE_INODES = (Int32)(8)
const MAT_HERMITIAN = (Int32)(9)
const MAT_SYMMETRY_ETERNAL = (Int32)(10)
const MAT_DUMMY = (Int32)(11)
const MAT_IGNORE_LOWER_TRIANGULAR = (Int32)(12)
const MAT_ERROR_LOWER_TRIANGULAR = (Int32)(13)
const MAT_GETROW_UPPERTRIANGULAR = (Int32)(14)
const MAT_SPD = (Int32)(15)
const MAT_NO_OFF_PROC_ZERO_ROWS = (Int32)(16)
const MAT_NO_OFF_PROC_ENTRIES = (Int32)(17)
const MAT_NEW_NONZERO_LOCATIONS = (Int32)(18)
const MAT_OPTION_MAX = (Int32)(19)

#= # end enum MatOption =#
#= # begin enum ANONYMOUS_31 =#
typealias ANONYMOUS_31 UInt32

const MAT_DO_NOT_COPY_VALUES = (UInt32)(0)
const MAT_COPY_VALUES = (UInt32)(1)
const MAT_SHARE_NONZERO_PATTERN = (UInt32)(2)

#= # end enum ANONYMOUS_31 =#
#= # begin enum MatDuplicateOption =#
typealias MatDuplicateOption UInt32

const MAT_DO_NOT_COPY_VALUES = (UInt32)(0)
const MAT_COPY_VALUES = (UInt32)(1)
const MAT_SHARE_NONZERO_PATTERN = (UInt32)(2)

#= # end enum MatDuplicateOption =#
immutable MatInfo
    block_size::PetscLogDouble
    nz_allocated::PetscLogDouble
    nz_used::PetscLogDouble
    nz_unneeded::PetscLogDouble
    memory::PetscLogDouble
    assemblies::PetscLogDouble
    mallocs::PetscLogDouble
    fill_ratio_given::PetscLogDouble
    fill_ratio_needed::PetscLogDouble
    factor_mallocs::PetscLogDouble
end

#= # begin enum ANONYMOUS_32 =#
typealias ANONYMOUS_32 UInt32

const MAT_LOCAL = (UInt32)(1)
const MAT_GLOBAL_MAX = (UInt32)(2)
const MAT_GLOBAL_SUM = (UInt32)(3)

#= # end enum ANONYMOUS_32 =#
#= # begin enum MatInfoType =#
typealias MatInfoType UInt32

const MAT_LOCAL = (UInt32)(1)
const MAT_GLOBAL_MAX = (UInt32)(2)
const MAT_GLOBAL_SUM = (UInt32)(3)

#= # end enum MatInfoType =#
typealias MatOrderingType Symbol

#= # begin enum ANONYMOUS_33 =#
typealias ANONYMOUS_33 UInt32

const MAT_SHIFT_NONE = (UInt32)(0)
const MAT_SHIFT_NONZERO = (UInt32)(1)
const MAT_SHIFT_POSITIVE_DEFINITE = (UInt32)(2)
const MAT_SHIFT_INBLOCKS = (UInt32)(3)

#= # end enum ANONYMOUS_33 =#
#= # begin enum MatFactorShiftType =#
typealias MatFactorShiftType UInt32

const MAT_SHIFT_NONE = (UInt32)(0)
const MAT_SHIFT_NONZERO = (UInt32)(1)
const MAT_SHIFT_POSITIVE_DEFINITE = (UInt32)(2)
const MAT_SHIFT_INBLOCKS = (UInt32)(3)

#= # end enum MatFactorShiftType =#
#= skipping type declaration with undefined symbols:
immutable MatFactorInfo
    diagonal_fill::PetscReal
    usedt::PetscReal
    dt::PetscReal
    dtcol::PetscReal
    dtcount::PetscReal
    fill::PetscReal
    levels::PetscReal
    pivotinblocks::PetscReal
    zeropivot::PetscReal
    shifttype::PetscReal
    shiftamount::PetscReal
end 
=#
#= # begin enum ANONYMOUS_34 =#
typealias ANONYMOUS_34 UInt32

const SOR_FORWARD_SWEEP = (UInt32)(1)
const SOR_BACKWARD_SWEEP = (UInt32)(2)
const SOR_SYMMETRIC_SWEEP = (UInt32)(3)
const SOR_LOCAL_FORWARD_SWEEP = (UInt32)(4)
const SOR_LOCAL_BACKWARD_SWEEP = (UInt32)(8)
const SOR_LOCAL_SYMMETRIC_SWEEP = (UInt32)(12)
const SOR_ZERO_INITIAL_GUESS = (UInt32)(16)
const SOR_EISENSTAT = (UInt32)(32)
const SOR_APPLY_UPPER = (UInt32)(64)
const SOR_APPLY_LOWER = (UInt32)(128)

#= # end enum ANONYMOUS_34 =#
#= # begin enum MatSORType =#
typealias MatSORType UInt32

const SOR_FORWARD_SWEEP = (UInt32)(1)
const SOR_BACKWARD_SWEEP = (UInt32)(2)
const SOR_SYMMETRIC_SWEEP = (UInt32)(3)
const SOR_LOCAL_FORWARD_SWEEP = (UInt32)(4)
const SOR_LOCAL_BACKWARD_SWEEP = (UInt32)(8)
const SOR_LOCAL_SYMMETRIC_SWEEP = (UInt32)(12)
const SOR_ZERO_INITIAL_GUESS = (UInt32)(16)
const SOR_EISENSTAT = (UInt32)(32)
const SOR_APPLY_UPPER = (UInt32)(64)
const SOR_APPLY_LOWER = (UInt32)(128)

#= # end enum MatSORType =#
# skipping undefined typealias typealias MatColoring Ptr{_p_MatColoring}
typealias MatColoringType Symbol

#= # begin enum ANONYMOUS_35 =#
typealias ANONYMOUS_35 UInt32

const MAT_COLORING_WEIGHT_RANDOM = (UInt32)(0)
const MAT_COLORING_WEIGHT_LEXICAL = (UInt32)(1)
const MAT_COLORING_WEIGHT_LF = (UInt32)(2)
const MAT_COLORING_WEIGHT_SL = (UInt32)(3)

#= # end enum ANONYMOUS_35 =#
#= # begin enum MatColoringWeightType =#
typealias MatColoringWeightType UInt32

const MAT_COLORING_WEIGHT_RANDOM = (UInt32)(0)
const MAT_COLORING_WEIGHT_LEXICAL = (UInt32)(1)
const MAT_COLORING_WEIGHT_LF = (UInt32)(2)
const MAT_COLORING_WEIGHT_SL = (UInt32)(3)

#= # end enum MatColoringWeightType =#
# skipping undefined typealias typealias MatFDColoring Ptr{_p_MatFDColoring}
# skipping undefined typealias typealias MatTransposeColoring Ptr{_p_MatTransposeColoring}
# skipping undefined typealias typealias MatPartitioning Ptr{_p_MatPartitioning}
typealias MatPartitioningType Symbol

#= # begin enum ANONYMOUS_36 =#
typealias ANONYMOUS_36 UInt32

const MP_CHACO_MULTILEVEL = (UInt32)(1)
const MP_CHACO_SPECTRAL = (UInt32)(2)
const MP_CHACO_LINEAR = (UInt32)(4)
const MP_CHACO_RANDOM = (UInt32)(5)
const MP_CHACO_SCATTERED = (UInt32)(6)

#= # end enum ANONYMOUS_36 =#
#= # begin enum MPChacoGlobalType =#
typealias MPChacoGlobalType UInt32

const MP_CHACO_MULTILEVEL = (UInt32)(1)
const MP_CHACO_SPECTRAL = (UInt32)(2)
const MP_CHACO_LINEAR = (UInt32)(4)
const MP_CHACO_RANDOM = (UInt32)(5)
const MP_CHACO_SCATTERED = (UInt32)(6)

#= # end enum MPChacoGlobalType =#
#= # begin enum ANONYMOUS_37 =#
typealias ANONYMOUS_37 UInt32

const MP_CHACO_KERNIGHAN = (UInt32)(1)
const MP_CHACO_NONE = (UInt32)(2)

#= # end enum ANONYMOUS_37 =#
#= # begin enum MPChacoLocalType =#
typealias MPChacoLocalType UInt32

const MP_CHACO_KERNIGHAN = (UInt32)(1)
const MP_CHACO_NONE = (UInt32)(2)

#= # end enum MPChacoLocalType =#
#= # begin enum ANONYMOUS_38 =#
typealias ANONYMOUS_38 UInt32

const MP_CHACO_LANCZOS = (UInt32)(0)
const MP_CHACO_RQI = (UInt32)(1)

#= # end enum ANONYMOUS_38 =#
#= # begin enum MPChacoEigenType =#
typealias MPChacoEigenType UInt32

const MP_CHACO_LANCZOS = (UInt32)(0)
const MP_CHACO_RQI = (UInt32)(1)

#= # end enum MPChacoEigenType =#
#= # begin enum ANONYMOUS_39 =#
typealias ANONYMOUS_39 UInt32

const MP_PTSCOTCH_QUALITY = (UInt32)(0)
const MP_PTSCOTCH_SPEED = (UInt32)(1)
const MP_PTSCOTCH_BALANCE = (UInt32)(2)
const MP_PTSCOTCH_SAFETY = (UInt32)(3)
const MP_PTSCOTCH_SCALABILITY = (UInt32)(4)

#= # end enum ANONYMOUS_39 =#
#= # begin enum MPPTScotchStrategyType =#
typealias MPPTScotchStrategyType UInt32

const MP_PTSCOTCH_QUALITY = (UInt32)(0)
const MP_PTSCOTCH_SPEED = (UInt32)(1)
const MP_PTSCOTCH_BALANCE = (UInt32)(2)
const MP_PTSCOTCH_SAFETY = (UInt32)(3)
const MP_PTSCOTCH_SCALABILITY = (UInt32)(4)

#= # end enum MPPTScotchStrategyType =#
# skipping undefined typealias typealias MatCoarsen Ptr{_p_MatCoarsen}
typealias MatCoarsenType Symbol

#= skipping type declaration with undefined symbols:
immutable PetscCDIntNd
    next::Ptr{_PetscCDIntNd}
    gid::PetscInt
end 
=#
#= skipping type declaration with undefined symbols:
immutable PetscCDArrNd
    next::Ptr{_PetscCDArrNd}
    array::Ptr{_PetscCDIntNd}
end 
=#
#= skipping type declaration with undefined symbols:
immutable PetscCoarsenData
    pool_list::PetscCDArrNd
    new_node::Ptr{PetscCDIntNd}
    new_left::PetscInt
    chk_sz::PetscInt
    extra_nodes::Ptr{PetscCDIntNd}
    array::Ptr{Ptr{PetscCDIntNd}}
    size::PetscInt
    mat::Mat
end 
=#
#= # begin enum ANONYMOUS_40 =#
typealias ANONYMOUS_40 UInt32

const MATOP_SET_VALUES = (UInt32)(0)
const MATOP_GET_ROW = (UInt32)(1)
const MATOP_RESTORE_ROW = (UInt32)(2)
const MATOP_MULT = (UInt32)(3)
const MATOP_MULT_ADD = (UInt32)(4)
const MATOP_MULT_TRANSPOSE = (UInt32)(5)
const MATOP_MULT_TRANSPOSE_ADD = (UInt32)(6)
const MATOP_SOLVE = (UInt32)(7)
const MATOP_SOLVE_ADD = (UInt32)(8)
const MATOP_SOLVE_TRANSPOSE = (UInt32)(9)
const MATOP_SOLVE_TRANSPOSE_ADD = (UInt32)(10)
const MATOP_LUFACTOR = (UInt32)(11)
const MATOP_CHOLESKYFACTOR = (UInt32)(12)
const MATOP_SOR = (UInt32)(13)
const MATOP_TRANSPOSE = (UInt32)(14)
const MATOP_GETINFO = (UInt32)(15)
const MATOP_EQUAL = (UInt32)(16)
const MATOP_GET_DIAGONAL = (UInt32)(17)
const MATOP_DIAGONAL_SCALE = (UInt32)(18)
const MATOP_NORM = (UInt32)(19)
const MATOP_ASSEMBLY_BEGIN = (UInt32)(20)
const MATOP_ASSEMBLY_END = (UInt32)(21)
const MATOP_SET_OPTION = (UInt32)(22)
const MATOP_ZERO_ENTRIES = (UInt32)(23)
const MATOP_ZERO_ROWS = (UInt32)(24)
const MATOP_LUFACTOR_SYMBOLIC = (UInt32)(25)
const MATOP_LUFACTOR_NUMERIC = (UInt32)(26)
const MATOP_CHOLESKY_FACTOR_SYMBOLIC = (UInt32)(27)
const MATOP_CHOLESKY_FACTOR_NUMERIC = (UInt32)(28)
const MATOP_SETUP_PREALLOCATION = (UInt32)(29)
const MATOP_ILUFACTOR_SYMBOLIC = (UInt32)(30)
const MATOP_ICCFACTOR_SYMBOLIC = (UInt32)(31)
const MATOP_DUPLICATE = (UInt32)(34)
const MATOP_FORWARD_SOLVE = (UInt32)(35)
const MATOP_BACKWARD_SOLVE = (UInt32)(36)
const MATOP_ILUFACTOR = (UInt32)(37)
const MATOP_ICCFACTOR = (UInt32)(38)
const MATOP_AXPY = (UInt32)(39)
const MATOP_GET_SUBMATRICES = (UInt32)(40)
const MATOP_INCREASE_OVERLAP = (UInt32)(41)
const MATOP_GET_VALUES = (UInt32)(42)
const MATOP_COPY = (UInt32)(43)
const MATOP_GET_ROW_MAX = (UInt32)(44)
const MATOP_SCALE = (UInt32)(45)
const MATOP_SHIFT = (UInt32)(46)
const MATOP_DIAGONAL_SET = (UInt32)(47)
const MATOP_ZERO_ROWS_COLUMNS = (UInt32)(48)
const MATOP_SET_RANDOM = (UInt32)(49)
const MATOP_GET_ROW_IJ = (UInt32)(50)
const MATOP_RESTORE_ROW_IJ = (UInt32)(51)
const MATOP_GET_COLUMN_IJ = (UInt32)(52)
const MATOP_RESTORE_COLUMN_IJ = (UInt32)(53)
const MATOP_FDCOLORING_CREATE = (UInt32)(54)
const MATOP_COLORING_PATCH = (UInt32)(55)
const MATOP_SET_UNFACTORED = (UInt32)(56)
const MATOP_PERMUTE = (UInt32)(57)
const MATOP_SET_VALUES_BLOCKED = (UInt32)(58)
const MATOP_GET_SUBMATRIX = (UInt32)(59)
const MATOP_DESTROY = (UInt32)(60)
const MATOP_VIEW = (UInt32)(61)
const MATOP_CONVERT_FROM = (UInt32)(62)
const MATOP_MATMAT_MULT = (UInt32)(63)
const MATOP_MATMAT_MULT_SYMBOLIC = (UInt32)(64)
const MATOP_MATMAT_MULT_NUMERIC = (UInt32)(65)
const MATOP_SET_LOCAL_TO_GLOBAL_MAP = (UInt32)(66)
const MATOP_SET_VALUES_LOCAL = (UInt32)(67)
const MATOP_ZERO_ROWS_LOCAL = (UInt32)(68)
const MATOP_GET_ROW_MAX_ABS = (UInt32)(69)
const MATOP_GET_ROW_MIN_ABS = (UInt32)(70)
const MATOP_CONVERT = (UInt32)(71)
const MATOP_SET_COLORING = (UInt32)(72)
const MATOP_SET_VALUES_ADIFOR = (UInt32)(74)
const MATOP_FD_COLORING_APPLY = (UInt32)(75)
const MATOP_SET_FROM_OPTIONS = (UInt32)(76)
const MATOP_MULT_CONSTRAINED = (UInt32)(77)
const MATOP_MULT_TRANSPOSE_CONSTRAIN = (UInt32)(78)
const MATOP_FIND_ZERO_DIAGONALS = (UInt32)(79)
const MATOP_MULT_MULTIPLE = (UInt32)(80)
const MATOP_SOLVE_MULTIPLE = (UInt32)(81)
const MATOP_GET_INERTIA = (UInt32)(82)
const MATOP_LOAD = (UInt32)(83)
const MATOP_IS_SYMMETRIC = (UInt32)(84)
const MATOP_IS_HERMITIAN = (UInt32)(85)
const MATOP_IS_STRUCTURALLY_SYMMETRIC = (UInt32)(86)
const MATOP_SET_VALUES_BLOCKEDLOCAL = (UInt32)(87)
const MATOP_GET_VECS = (UInt32)(88)
const MATOP_MAT_MULT = (UInt32)(89)
const MATOP_MAT_MULT_SYMBOLIC = (UInt32)(90)
const MATOP_MAT_MULT_NUMERIC = (UInt32)(91)
const MATOP_PTAP = (UInt32)(92)
const MATOP_PTAP_SYMBOLIC = (UInt32)(93)
const MATOP_PTAP_NUMERIC = (UInt32)(94)
const MATOP_MAT_TRANSPOSE_MULT = (UInt32)(95)
const MATOP_MAT_TRANSPOSE_MULT_SYMBO = (UInt32)(96)
const MATOP_MAT_TRANSPOSE_MULT_NUMER = (UInt32)(97)
const MATOP_CONJUGATE = (UInt32)(102)
const MATOP_SET_VALUES_ROW = (UInt32)(104)
const MATOP_REAL_PART = (UInt32)(105)
const MATOP_IMAGINARY_PART = (UInt32)(106)
const MATOP_GET_ROW_UPPER_TRIANGULAR = (UInt32)(107)
const MATOP_RESTORE_ROW_UPPER_TRIANG = (UInt32)(108)
const MATOP_MAT_SOLVE = (UInt32)(109)
const MATOP_GET_REDUNDANT_MATRIX = (UInt32)(110)
const MATOP_GET_ROW_MIN = (UInt32)(111)
const MATOP_GET_COLUMN_VECTOR = (UInt32)(112)
const MATOP_MISSING_DIAGONAL = (UInt32)(113)
const MATOP_GET_SEQ_NONZERO_STRUCTUR = (UInt32)(114)
const MATOP_CREATE = (UInt32)(115)
const MATOP_GET_GHOSTS = (UInt32)(116)
const MATOP_GET_LOCAL_SUB_MATRIX = (UInt32)(117)
const MATOP_RESTORE_LOCALSUB_MATRIX = (UInt32)(118)
const MATOP_MULT_DIAGONAL_BLOCK = (UInt32)(119)
const MATOP_HERMITIAN_TRANSPOSE = (UInt32)(120)
const MATOP_MULT_HERMITIAN_TRANSPOSE = (UInt32)(121)
const MATOP_MULT_HERMITIAN_TRANS_ADD = (UInt32)(122)
const MATOP_GET_MULTI_PROC_BLOCK = (UInt32)(123)
const MATOP_FIND_NONZERO_ROWS = (UInt32)(124)
const MATOP_GET_COLUMN_NORMS = (UInt32)(125)
const MATOP_INVERT_BLOCK_DIAGONAL = (UInt32)(126)
const MATOP_GET_SUB_MATRICES_PARALLE = (UInt32)(128)
const MATOP_SET_VALUES_BATCH = (UInt32)(129)
const MATOP_TRANSPOSE_MAT_MULT = (UInt32)(130)
const MATOP_TRANSPOSE_MAT_MULT_SYMBO = (UInt32)(131)
const MATOP_TRANSPOSE_MAT_MULT_NUMER = (UInt32)(132)
const MATOP_TRANSPOSE_COLORING_CREAT = (UInt32)(133)
const MATOP_TRANS_COLORING_APPLY_SPT = (UInt32)(134)
const MATOP_TRANS_COLORING_APPLY_DEN = (UInt32)(135)
const MATOP_RART = (UInt32)(136)
const MATOP_RART_SYMBOLIC = (UInt32)(137)
const MATOP_RART_NUMERIC = (UInt32)(138)
const MATOP_SET_BLOCK_SIZES = (UInt32)(139)
const MATOP_AYPX = (UInt32)(140)
const MATOP_RESIDUAL = (UInt32)(141)
const MATOP_FDCOLORING_SETUP = (UInt32)(142)
const MATOP_MPICONCATENATESEQ = (UInt32)(144)

#= # end enum ANONYMOUS_40 =#
#= # begin enum MatOperation =#
typealias MatOperation UInt32

const MATOP_SET_VALUES = (UInt32)(0)
const MATOP_GET_ROW = (UInt32)(1)
const MATOP_RESTORE_ROW = (UInt32)(2)
const MATOP_MULT = (UInt32)(3)
const MATOP_MULT_ADD = (UInt32)(4)
const MATOP_MULT_TRANSPOSE = (UInt32)(5)
const MATOP_MULT_TRANSPOSE_ADD = (UInt32)(6)
const MATOP_SOLVE = (UInt32)(7)
const MATOP_SOLVE_ADD = (UInt32)(8)
const MATOP_SOLVE_TRANSPOSE = (UInt32)(9)
const MATOP_SOLVE_TRANSPOSE_ADD = (UInt32)(10)
const MATOP_LUFACTOR = (UInt32)(11)
const MATOP_CHOLESKYFACTOR = (UInt32)(12)
const MATOP_SOR = (UInt32)(13)
const MATOP_TRANSPOSE = (UInt32)(14)
const MATOP_GETINFO = (UInt32)(15)
const MATOP_EQUAL = (UInt32)(16)
const MATOP_GET_DIAGONAL = (UInt32)(17)
const MATOP_DIAGONAL_SCALE = (UInt32)(18)
const MATOP_NORM = (UInt32)(19)
const MATOP_ASSEMBLY_BEGIN = (UInt32)(20)
const MATOP_ASSEMBLY_END = (UInt32)(21)
const MATOP_SET_OPTION = (UInt32)(22)
const MATOP_ZERO_ENTRIES = (UInt32)(23)
const MATOP_ZERO_ROWS = (UInt32)(24)
const MATOP_LUFACTOR_SYMBOLIC = (UInt32)(25)
const MATOP_LUFACTOR_NUMERIC = (UInt32)(26)
const MATOP_CHOLESKY_FACTOR_SYMBOLIC = (UInt32)(27)
const MATOP_CHOLESKY_FACTOR_NUMERIC = (UInt32)(28)
const MATOP_SETUP_PREALLOCATION = (UInt32)(29)
const MATOP_ILUFACTOR_SYMBOLIC = (UInt32)(30)
const MATOP_ICCFACTOR_SYMBOLIC = (UInt32)(31)
const MATOP_DUPLICATE = (UInt32)(34)
const MATOP_FORWARD_SOLVE = (UInt32)(35)
const MATOP_BACKWARD_SOLVE = (UInt32)(36)
const MATOP_ILUFACTOR = (UInt32)(37)
const MATOP_ICCFACTOR = (UInt32)(38)
const MATOP_AXPY = (UInt32)(39)
const MATOP_GET_SUBMATRICES = (UInt32)(40)
const MATOP_INCREASE_OVERLAP = (UInt32)(41)
const MATOP_GET_VALUES = (UInt32)(42)
const MATOP_COPY = (UInt32)(43)
const MATOP_GET_ROW_MAX = (UInt32)(44)
const MATOP_SCALE = (UInt32)(45)
const MATOP_SHIFT = (UInt32)(46)
const MATOP_DIAGONAL_SET = (UInt32)(47)
const MATOP_ZERO_ROWS_COLUMNS = (UInt32)(48)
const MATOP_SET_RANDOM = (UInt32)(49)
const MATOP_GET_ROW_IJ = (UInt32)(50)
const MATOP_RESTORE_ROW_IJ = (UInt32)(51)
const MATOP_GET_COLUMN_IJ = (UInt32)(52)
const MATOP_RESTORE_COLUMN_IJ = (UInt32)(53)
const MATOP_FDCOLORING_CREATE = (UInt32)(54)
const MATOP_COLORING_PATCH = (UInt32)(55)
const MATOP_SET_UNFACTORED = (UInt32)(56)
const MATOP_PERMUTE = (UInt32)(57)
const MATOP_SET_VALUES_BLOCKED = (UInt32)(58)
const MATOP_GET_SUBMATRIX = (UInt32)(59)
const MATOP_DESTROY = (UInt32)(60)
const MATOP_VIEW = (UInt32)(61)
const MATOP_CONVERT_FROM = (UInt32)(62)
const MATOP_MATMAT_MULT = (UInt32)(63)
const MATOP_MATMAT_MULT_SYMBOLIC = (UInt32)(64)
const MATOP_MATMAT_MULT_NUMERIC = (UInt32)(65)
const MATOP_SET_LOCAL_TO_GLOBAL_MAP = (UInt32)(66)
const MATOP_SET_VALUES_LOCAL = (UInt32)(67)
const MATOP_ZERO_ROWS_LOCAL = (UInt32)(68)
const MATOP_GET_ROW_MAX_ABS = (UInt32)(69)
const MATOP_GET_ROW_MIN_ABS = (UInt32)(70)
const MATOP_CONVERT = (UInt32)(71)
const MATOP_SET_COLORING = (UInt32)(72)
const MATOP_SET_VALUES_ADIFOR = (UInt32)(74)
const MATOP_FD_COLORING_APPLY = (UInt32)(75)
const MATOP_SET_FROM_OPTIONS = (UInt32)(76)
const MATOP_MULT_CONSTRAINED = (UInt32)(77)
const MATOP_MULT_TRANSPOSE_CONSTRAIN = (UInt32)(78)
const MATOP_FIND_ZERO_DIAGONALS = (UInt32)(79)
const MATOP_MULT_MULTIPLE = (UInt32)(80)
const MATOP_SOLVE_MULTIPLE = (UInt32)(81)
const MATOP_GET_INERTIA = (UInt32)(82)
const MATOP_LOAD = (UInt32)(83)
const MATOP_IS_SYMMETRIC = (UInt32)(84)
const MATOP_IS_HERMITIAN = (UInt32)(85)
const MATOP_IS_STRUCTURALLY_SYMMETRIC = (UInt32)(86)
const MATOP_SET_VALUES_BLOCKEDLOCAL = (UInt32)(87)
const MATOP_GET_VECS = (UInt32)(88)
const MATOP_MAT_MULT = (UInt32)(89)
const MATOP_MAT_MULT_SYMBOLIC = (UInt32)(90)
const MATOP_MAT_MULT_NUMERIC = (UInt32)(91)
const MATOP_PTAP = (UInt32)(92)
const MATOP_PTAP_SYMBOLIC = (UInt32)(93)
const MATOP_PTAP_NUMERIC = (UInt32)(94)
const MATOP_MAT_TRANSPOSE_MULT = (UInt32)(95)
const MATOP_MAT_TRANSPOSE_MULT_SYMBO = (UInt32)(96)
const MATOP_MAT_TRANSPOSE_MULT_NUMER = (UInt32)(97)
const MATOP_CONJUGATE = (UInt32)(102)
const MATOP_SET_VALUES_ROW = (UInt32)(104)
const MATOP_REAL_PART = (UInt32)(105)
const MATOP_IMAGINARY_PART = (UInt32)(106)
const MATOP_GET_ROW_UPPER_TRIANGULAR = (UInt32)(107)
const MATOP_RESTORE_ROW_UPPER_TRIANG = (UInt32)(108)
const MATOP_MAT_SOLVE = (UInt32)(109)
const MATOP_GET_REDUNDANT_MATRIX = (UInt32)(110)
const MATOP_GET_ROW_MIN = (UInt32)(111)
const MATOP_GET_COLUMN_VECTOR = (UInt32)(112)
const MATOP_MISSING_DIAGONAL = (UInt32)(113)
const MATOP_GET_SEQ_NONZERO_STRUCTUR = (UInt32)(114)
const MATOP_CREATE = (UInt32)(115)
const MATOP_GET_GHOSTS = (UInt32)(116)
const MATOP_GET_LOCAL_SUB_MATRIX = (UInt32)(117)
const MATOP_RESTORE_LOCALSUB_MATRIX = (UInt32)(118)
const MATOP_MULT_DIAGONAL_BLOCK = (UInt32)(119)
const MATOP_HERMITIAN_TRANSPOSE = (UInt32)(120)
const MATOP_MULT_HERMITIAN_TRANSPOSE = (UInt32)(121)
const MATOP_MULT_HERMITIAN_TRANS_ADD = (UInt32)(122)
const MATOP_GET_MULTI_PROC_BLOCK = (UInt32)(123)
const MATOP_FIND_NONZERO_ROWS = (UInt32)(124)
const MATOP_GET_COLUMN_NORMS = (UInt32)(125)
const MATOP_INVERT_BLOCK_DIAGONAL = (UInt32)(126)
const MATOP_GET_SUB_MATRICES_PARALLE = (UInt32)(128)
const MATOP_SET_VALUES_BATCH = (UInt32)(129)
const MATOP_TRANSPOSE_MAT_MULT = (UInt32)(130)
const MATOP_TRANSPOSE_MAT_MULT_SYMBO = (UInt32)(131)
const MATOP_TRANSPOSE_MAT_MULT_NUMER = (UInt32)(132)
const MATOP_TRANSPOSE_COLORING_CREAT = (UInt32)(133)
const MATOP_TRANS_COLORING_APPLY_SPT = (UInt32)(134)
const MATOP_TRANS_COLORING_APPLY_DEN = (UInt32)(135)
const MATOP_RART = (UInt32)(136)
const MATOP_RART_SYMBOLIC = (UInt32)(137)
const MATOP_RART_NUMERIC = (UInt32)(138)
const MATOP_SET_BLOCK_SIZES = (UInt32)(139)
const MATOP_AYPX = (UInt32)(140)
const MATOP_RESIDUAL = (UInt32)(141)
const MATOP_FDCOLORING_SETUP = (UInt32)(142)
const MATOP_MPICONCATENATESEQ = (UInt32)(144)

#= # end enum MatOperation =#
# skipping undefined typealias typealias MatNullSpace Ptr{_p_MatNullSpace}
# skipping undefined typealias typealias MatMFFD Ptr{_p_MatMFFD}
typealias MatMFFDType Symbol

# skipping undefined typealias typealias DM Ptr{_p_DM}
#= # begin enum ANONYMOUS_41 =#
typealias ANONYMOUS_41 UInt32

const DM_BOUNDARY_NONE = (UInt32)(0)
const DM_BOUNDARY_GHOSTED = (UInt32)(1)
const DM_BOUNDARY_MIRROR = (UInt32)(2)
const DM_BOUNDARY_PERIODIC = (UInt32)(3)
const DM_BOUNDARY_TWIST = (UInt32)(4)

#= # end enum ANONYMOUS_41 =#
#= # begin enum DMBoundaryType =#
typealias DMBoundaryType UInt32

const DM_BOUNDARY_NONE = (UInt32)(0)
const DM_BOUNDARY_GHOSTED = (UInt32)(1)
const DM_BOUNDARY_MIRROR = (UInt32)(2)
const DM_BOUNDARY_PERIODIC = (UInt32)(3)
const DM_BOUNDARY_TWIST = (UInt32)(4)

#= # end enum DMBoundaryType =#
# skipping undefined typealias typealias PetscPartitioner Ptr{_p_PetscPartitioner}
# skipping undefined typealias typealias PetscSpace Ptr{_p_PetscSpace}
# skipping undefined typealias typealias PetscDualSpace Ptr{_p_PetscDualSpace}
# skipping undefined typealias typealias PetscFE Ptr{_p_PetscFE}
# skipping undefined typealias typealias PetscDS Ptr{_p_PetscDS}
typealias DMType Symbol

immutable NLF_DAAD
end

typealias NLF Ptr{NLF_DAAD}

#= # begin enum ANONYMOUS_42 =#
typealias ANONYMOUS_42 UInt32

const PETSC_UNIT_LENGTH = (UInt32)(0)
const PETSC_UNIT_MASS = (UInt32)(1)
const PETSC_UNIT_TIME = (UInt32)(2)
const PETSC_UNIT_CURRENT = (UInt32)(3)
const PETSC_UNIT_TEMPERATURE = (UInt32)(4)
const PETSC_UNIT_AMOUNT = (UInt32)(5)
const PETSC_UNIT_LUMINOSITY = (UInt32)(6)
const NUM_PETSC_UNITS = (UInt32)(7)

#= # end enum ANONYMOUS_42 =#
#= # begin enum PetscUnit =#
typealias PetscUnit UInt32

const PETSC_UNIT_LENGTH = (UInt32)(0)
const PETSC_UNIT_MASS = (UInt32)(1)
const PETSC_UNIT_TIME = (UInt32)(2)
const PETSC_UNIT_CURRENT = (UInt32)(3)
const PETSC_UNIT_TEMPERATURE = (UInt32)(4)
const PETSC_UNIT_AMOUNT = (UInt32)(5)
const PETSC_UNIT_LUMINOSITY = (UInt32)(6)
const NUM_PETSC_UNITS = (UInt32)(7)

#= # end enum PetscUnit =#
# skipping undefined typealias typealias DMInterpolationInfo Ptr{_DMInterpolationInfo}
#= # begin enum ANONYMOUS_43 =#
typealias ANONYMOUS_43 UInt32

const DMDA_STENCIL_STAR = (UInt32)(0)
const DMDA_STENCIL_BOX = (UInt32)(1)

#= # end enum ANONYMOUS_43 =#
#= # begin enum DMDAStencilType =#
typealias DMDAStencilType UInt32

const DMDA_STENCIL_STAR = (UInt32)(0)
const DMDA_STENCIL_BOX = (UInt32)(1)

#= # end enum DMDAStencilType =#
#= # begin enum ANONYMOUS_44 =#
typealias ANONYMOUS_44 UInt32

const DMDA_Q0 = (UInt32)(0)
const DMDA_Q1 = (UInt32)(1)

#= # end enum ANONYMOUS_44 =#
#= # begin enum DMDAInterpolationType =#
typealias DMDAInterpolationType UInt32

const DMDA_Q0 = (UInt32)(0)
const DMDA_Q1 = (UInt32)(1)

#= # end enum DMDAInterpolationType =#
#= # begin enum ANONYMOUS_45 =#
typealias ANONYMOUS_45 UInt32

const DMDA_ELEMENT_P1 = (UInt32)(0)
const DMDA_ELEMENT_Q1 = (UInt32)(1)

#= # end enum ANONYMOUS_45 =#
#= # begin enum DMDAElementType =#
typealias DMDAElementType UInt32

const DMDA_ELEMENT_P1 = (UInt32)(0)
const DMDA_ELEMENT_Q1 = (UInt32)(1)

#= # end enum DMDAElementType =#
#= skipping type declaration with undefined symbols:
immutable DMDALocalInfo
    dim::PetscInt
    dof::PetscInt
    sw::PetscInt
    mx::PetscInt
    my::PetscInt
    mz::PetscInt
    xs::PetscInt
    ys::PetscInt
    zs::PetscInt
    xm::PetscInt
    ym::PetscInt
    zm::PetscInt
    gxs::PetscInt
    gys::PetscInt
    gzs::PetscInt
    gxm::PetscInt
    gym::PetscInt
    gzm::PetscInt
    bx::DMBoundaryType
    by::DMBoundaryType
    bz::DMBoundaryType
    st::DMDAStencilType
    da::DM
end 
=#
typealias PFType Symbol

# skipping undefined typealias typealias PF Ptr{_p_PF}
typealias AOType Symbol

# skipping undefined typealias typealias PetscQuadrature Ptr{_p_PetscQuadrature}
#= skipping type declaration with undefined symbols:
immutable Array_3_PetscReal
    d1::PetscReal
    d2::PetscReal
    d3::PetscReal
end 
=#
#= skipping undefined expression zero(::Type{Array_3_PetscReal}) = begin  # /home/jared/.julia/v0.4/Clang_orig/src/wrap_c.jl, line 270:
        Array_3_PetscReal(fill(zero(Float32),3)...)
    end =#
#= skipping type declaration with undefined symbols:
immutable Array_9_PetscReal
    d1::PetscReal
    d2::PetscReal
    d3::PetscReal
    d4::PetscReal
    d5::PetscReal
    d6::PetscReal
    d7::PetscReal
    d8::PetscReal
    d9::PetscReal
end 
=#
#= skipping undefined expression zero(::Type{Array_9_PetscReal}) = begin  # /home/jared/.julia/v0.4/Clang_orig/src/wrap_c.jl, line 270:
        Array_9_PetscReal(fill(zero(Float32),9)...)
    end =#
#= skipping type declaration with undefined symbols:
immutable PetscFECellGeom
    v0::Array_3_PetscReal
    J::Array_9_PetscReal
    invJ::Array_9_PetscReal
    detJ::PetscReal
    n::Array_3_PetscReal
    dim::PetscInt
    dimEmbed::PetscInt
end 
=#
typealias PetscSpaceType Symbol
typealias PetscDualSpaceType Symbol
typealias PetscFEType Symbol

#= # begin enum ANONYMOUS_46 =#
typealias ANONYMOUS_46 UInt32

const DMDA_X = (UInt32)(0)
const DMDA_Y = (UInt32)(1)
const DMDA_Z = (UInt32)(2)

#= # end enum ANONYMOUS_46 =#
#= # begin enum DMDADirection =#
typealias DMDADirection UInt32

const DMDA_X = (UInt32)(0)
const DMDA_Y = (UInt32)(1)
const DMDA_Z = (UInt32)(2)

#= # end enum DMDADirection =#
#= skipping type declaration with undefined symbols:
immutable DMDACoor2d
    x::PetscScalar
    y::PetscScalar
end 
=#
#= skipping type declaration with undefined symbols:
immutable DMDACoor3d
    x::PetscScalar
    y::PetscScalar
    z::PetscScalar
end 
=#
# skipping undefined typealias typealias PetscLimiter Ptr{_p_PetscLimiter}
# skipping undefined typealias typealias PetscFV Ptr{_p_PetscFV}
#= skipping type declaration with undefined symbols:
immutable Array_3_PetscScalar
    d1::PetscScalar
    d2::PetscScalar
    d3::PetscScalar
end 
=#
#= skipping undefined expression zero(::Type{Array_3_PetscScalar}) = begin  # /home/jared/.julia/v0.4/Clang_orig/src/wrap_c.jl, line 270:
        Array_3_PetscScalar(fill(zero(Float32),3)...)
    end =#
#= skipping type declaration with undefined symbols:
immutable Array_2_Array_3_PetscScalar
    d1::Array_3_PetscScalar
    d2::Array_3_PetscScalar
end 
=#
#= skipping undefined expression zero(::Type{Array_2_Array_3_PetscScalar}) = begin  # /home/jared/.julia/v0.4/Clang_orig/src/wrap_c.jl, line 270:
        Array_2_Array_3_PetscScalar(fill(zero(Array_3_PetscScalar),2)...)
    end =#
#= skipping type declaration with undefined symbols:
immutable PetscFVFaceGeom
    normal::Array_3_PetscReal
    centroid::Array_3_PetscReal
    grad::Array_2_Array_3_PetscScalar
end 
=#
#= skipping type declaration with undefined symbols:
immutable PetscFVCellGeom
    centroid::Array_3_PetscReal
    volume::PetscReal
end 
=#
typealias PetscLimiterType Symbol
typealias PetscFVType Symbol
typealias PetscPartitionerType Symbol

# skipping undefined typealias typealias DMLabel Ptr{_n_DMLabel}
# skipping undefined typealias typealias DMBoundary Ptr{_n_Boundary}
#= skipping type declaration with undefined symbols:
immutable JacActionCtx
    dm::DM
    u::Vec
    J::Mat
    user::Ptr{Void}
end 
=#
typealias PetscDSType Symbol
typealias PetscPointFunc Ptr{Void}
typealias PetscPointJac Ptr{Void}
typealias PetscBdPointFunc Ptr{Void}
typealias PetscBdPointJac Ptr{Void}
typealias PetscRiemannFunc Ptr{Void}

# skipping undefined typealias typealias Characteristic Ptr{_p_Characteristic}
typealias CharacteristicType Symbol

immutable PC{T}
    pobj::Ptr{Void}
end

typealias PCType Symbol

#= # begin enum ANONYMOUS_47 =#
typealias ANONYMOUS_47 Cint

const PC_SIDE_DEFAULT = (Int32)(-1)
const PC_LEFT = (Int32)(0)
const PC_RIGHT = (Int32)(1)
const PC_SYMMETRIC = (Int32)(2)

#= # end enum ANONYMOUS_47 =#
#= # begin enum ANONYMOUS_48 =#
typealias ANONYMOUS_48 Cint

const PCRICHARDSON_CONVERGED_RTOL = (Int32)(2)
const PCRICHARDSON_CONVERGED_ATOL = (Int32)(3)
const PCRICHARDSON_CONVERGED_ITS = (Int32)(4)
const PCRICHARDSON_DIVERGED_DTOL = (Int32)(-4)

#= # end enum ANONYMOUS_48 =#
#= # begin enum PCRichardsonConvergedReason =#
typealias PCRichardsonConvergedReason Cint

const PCRICHARDSON_CONVERGED_RTOL = (Int32)(2)
const PCRICHARDSON_CONVERGED_ATOL = (Int32)(3)
const PCRICHARDSON_CONVERGED_ITS = (Int32)(4)
const PCRICHARDSON_DIVERGED_DTOL = (Int32)(-4)

#= # end enum PCRichardsonConvergedReason =#
#= # begin enum ANONYMOUS_49 =#
typealias ANONYMOUS_49 UInt32

const PC_JACOBI_DIAGONAL = (UInt32)(0)
const PC_JACOBI_ROWMAX = (UInt32)(1)
const PC_JACOBI_ROWSUM = (UInt32)(2)

#= # end enum ANONYMOUS_49 =#
#= # begin enum PCJacobiType =#
typealias PCJacobiType UInt32

const PC_JACOBI_DIAGONAL = (UInt32)(0)
const PC_JACOBI_ROWMAX = (UInt32)(1)
const PC_JACOBI_ROWSUM = (UInt32)(2)

#= # end enum PCJacobiType =#
#= # begin enum ANONYMOUS_50 =#
typealias ANONYMOUS_50 UInt32

const PC_ASM_BASIC = (UInt32)(3)
const PC_ASM_RESTRICT = (UInt32)(1)
const PC_ASM_INTERPOLATE = (UInt32)(2)
const PC_ASM_NONE = (UInt32)(0)

#= # end enum ANONYMOUS_50 =#
#= # begin enum PCASMType =#
typealias PCASMType UInt32

const PC_ASM_BASIC = (UInt32)(3)
const PC_ASM_RESTRICT = (UInt32)(1)
const PC_ASM_INTERPOLATE = (UInt32)(2)
const PC_ASM_NONE = (UInt32)(0)

#= # end enum PCASMType =#
#= # begin enum ANONYMOUS_51 =#
typealias ANONYMOUS_51 UInt32

const PC_GASM_BASIC = (UInt32)(3)
const PC_GASM_RESTRICT = (UInt32)(1)
const PC_GASM_INTERPOLATE = (UInt32)(2)
const PC_GASM_NONE = (UInt32)(0)

#= # end enum ANONYMOUS_51 =#
#= # begin enum PCGASMType =#
typealias PCGASMType UInt32

const PC_GASM_BASIC = (UInt32)(3)
const PC_GASM_RESTRICT = (UInt32)(1)
const PC_GASM_INTERPOLATE = (UInt32)(2)
const PC_GASM_NONE = (UInt32)(0)

#= # end enum PCGASMType =#
#= # begin enum ANONYMOUS_52 =#
typealias ANONYMOUS_52 UInt32

const PC_COMPOSITE_ADDITIVE = (UInt32)(0)
const PC_COMPOSITE_MULTIPLICATIVE = (UInt32)(1)
const PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE = (UInt32)(2)
const PC_COMPOSITE_SPECIAL = (UInt32)(3)
const PC_COMPOSITE_SCHUR = (UInt32)(4)

#= # end enum ANONYMOUS_52 =#
#= # begin enum PCCompositeType =#
typealias PCCompositeType UInt32

const PC_COMPOSITE_ADDITIVE = (UInt32)(0)
const PC_COMPOSITE_MULTIPLICATIVE = (UInt32)(1)
const PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE = (UInt32)(2)
const PC_COMPOSITE_SPECIAL = (UInt32)(3)
const PC_COMPOSITE_SCHUR = (UInt32)(4)

#= # end enum PCCompositeType =#
#= # begin enum ANONYMOUS_53 =#
typealias ANONYMOUS_53 UInt32

const PC_FIELDSPLIT_SCHUR_PRE_SELF = (UInt32)(0)
const PC_FIELDSPLIT_SCHUR_PRE_SELFP = (UInt32)(1)
const PC_FIELDSPLIT_SCHUR_PRE_A11 = (UInt32)(2)
const PC_FIELDSPLIT_SCHUR_PRE_USER = (UInt32)(3)
const PC_FIELDSPLIT_SCHUR_PRE_FULL = (UInt32)(4)

#= # end enum ANONYMOUS_53 =#
#= # begin enum PCFieldSplitSchurPreType =#
typealias PCFieldSplitSchurPreType UInt32

const PC_FIELDSPLIT_SCHUR_PRE_SELF = (UInt32)(0)
const PC_FIELDSPLIT_SCHUR_PRE_SELFP = (UInt32)(1)
const PC_FIELDSPLIT_SCHUR_PRE_A11 = (UInt32)(2)
const PC_FIELDSPLIT_SCHUR_PRE_USER = (UInt32)(3)
const PC_FIELDSPLIT_SCHUR_PRE_FULL = (UInt32)(4)

#= # end enum PCFieldSplitSchurPreType =#
#= # begin enum ANONYMOUS_54 =#
typealias ANONYMOUS_54 UInt32

const PC_FIELDSPLIT_SCHUR_FACT_DIAG = (UInt32)(0)
const PC_FIELDSPLIT_SCHUR_FACT_LOWER = (UInt32)(1)
const PC_FIELDSPLIT_SCHUR_FACT_UPPER = (UInt32)(2)
const PC_FIELDSPLIT_SCHUR_FACT_FULL = (UInt32)(3)

#= # end enum ANONYMOUS_54 =#
#= # begin enum PCFieldSplitSchurFactType =#
typealias PCFieldSplitSchurFactType UInt32

const PC_FIELDSPLIT_SCHUR_FACT_DIAG = (UInt32)(0)
const PC_FIELDSPLIT_SCHUR_FACT_LOWER = (UInt32)(1)
const PC_FIELDSPLIT_SCHUR_FACT_UPPER = (UInt32)(2)
const PC_FIELDSPLIT_SCHUR_FACT_FULL = (UInt32)(3)

#= # end enum PCFieldSplitSchurFactType =#
#= # begin enum ANONYMOUS_55 =#
typealias ANONYMOUS_55 UInt32

const PC_PARMS_GLOBAL_RAS = (UInt32)(0)
const PC_PARMS_GLOBAL_SCHUR = (UInt32)(1)
const PC_PARMS_GLOBAL_BJ = (UInt32)(2)

#= # end enum ANONYMOUS_55 =#
#= # begin enum PCPARMSGlobalType =#
typealias PCPARMSGlobalType UInt32

const PC_PARMS_GLOBAL_RAS = (UInt32)(0)
const PC_PARMS_GLOBAL_SCHUR = (UInt32)(1)
const PC_PARMS_GLOBAL_BJ = (UInt32)(2)

#= # end enum PCPARMSGlobalType =#
#= # begin enum ANONYMOUS_56 =#
typealias ANONYMOUS_56 UInt32

const PC_PARMS_LOCAL_ILU0 = (UInt32)(0)
const PC_PARMS_LOCAL_ILUK = (UInt32)(1)
const PC_PARMS_LOCAL_ILUT = (UInt32)(2)
const PC_PARMS_LOCAL_ARMS = (UInt32)(3)

#= # end enum ANONYMOUS_56 =#
#= # begin enum PCPARMSLocalType =#
typealias PCPARMSLocalType UInt32

const PC_PARMS_LOCAL_ILU0 = (UInt32)(0)
const PC_PARMS_LOCAL_ILUK = (UInt32)(1)
const PC_PARMS_LOCAL_ILUT = (UInt32)(2)
const PC_PARMS_LOCAL_ARMS = (UInt32)(3)

#= # end enum PCPARMSLocalType =#
typealias PCGAMGType Symbol
typealias PCGAMGClassicalType Symbol

#= # begin enum ANONYMOUS_57 =#
typealias ANONYMOUS_57 UInt32

const PC_MG_MULTIPLICATIVE = (UInt32)(0)
const PC_MG_ADDITIVE = (UInt32)(1)
const PC_MG_FULL = (UInt32)(2)
const PC_MG_KASKADE = (UInt32)(3)

#= # end enum ANONYMOUS_57 =#
#= # begin enum ANONYMOUS_58 =#
typealias ANONYMOUS_58 UInt32

const PC_MG_CYCLE_V = (UInt32)(1)
const PC_MG_CYCLE_W = (UInt32)(2)

#= # end enum ANONYMOUS_58 =#
#= # begin enum PCMGCycleType =#
typealias PCMGCycleType UInt32

const PC_MG_CYCLE_V = (UInt32)(1)
const PC_MG_CYCLE_W = (UInt32)(2)

#= # end enum PCMGCycleType =#
#= # begin enum ANONYMOUS_59 =#
typealias ANONYMOUS_59 UInt32

const PC_EXOTIC_FACE = (UInt32)(0)
const PC_EXOTIC_WIREBASKET = (UInt32)(1)

#= # end enum ANONYMOUS_59 =#
#= # begin enum PCExoticType =#
typealias PCExoticType UInt32

const PC_EXOTIC_FACE = (UInt32)(0)
const PC_EXOTIC_WIREBASKET = (UInt32)(1)

#= # end enum PCExoticType =#
immutable KSP{T}
    pobj::Ptr{Void}
end

typealias KSPType Symbol

#= # begin enum ANONYMOUS_60 =#
typealias ANONYMOUS_60 UInt32

const KSP_FCG_TRUNC_TYPE_STANDARD = (UInt32)(0)
const KSP_FCG_TRUNC_TYPE_NOTAY = (UInt32)(1)

#= # end enum ANONYMOUS_60 =#
#= # begin enum KSPFCGTruncationType =#
typealias KSPFCGTruncationType UInt32

const KSP_FCG_TRUNC_TYPE_STANDARD = (UInt32)(0)
const KSP_FCG_TRUNC_TYPE_NOTAY = (UInt32)(1)

#= # end enum KSPFCGTruncationType =#
#= # begin enum ANONYMOUS_61 =#
typealias ANONYMOUS_61 UInt32

const KSP_GMRES_CGS_REFINE_NEVER = (UInt32)(0)
const KSP_GMRES_CGS_REFINE_IFNEEDED = (UInt32)(1)
const KSP_GMRES_CGS_REFINE_ALWAYS = (UInt32)(2)

#= # end enum ANONYMOUS_61 =#
#= # begin enum KSPGMRESCGSRefinementType =#
typealias KSPGMRESCGSRefinementType UInt32

const KSP_GMRES_CGS_REFINE_NEVER = (UInt32)(0)
const KSP_GMRES_CGS_REFINE_IFNEEDED = (UInt32)(1)
const KSP_GMRES_CGS_REFINE_ALWAYS = (UInt32)(2)

#= # end enum KSPGMRESCGSRefinementType =#
#= # begin enum ANONYMOUS_62 =#
typealias ANONYMOUS_62 Cint

const KSP_NORM_DEFAULT = (Int32)(-1)
const KSP_NORM_NONE = (Int32)(0)
const KSP_NORM_PRECONDITIONED = (Int32)(1)
const KSP_NORM_UNPRECONDITIONED = (Int32)(2)
const KSP_NORM_NATURAL = (Int32)(3)

#= # end enum ANONYMOUS_62 =#
#= # begin enum ANONYMOUS_63 =#
typealias ANONYMOUS_63 Cint

const KSP_CONVERGED_RTOL_NORMAL = (Int32)(1)
const KSP_CONVERGED_ATOL_NORMAL = (Int32)(9)
const KSP_CONVERGED_RTOL = (Int32)(2)
const KSP_CONVERGED_ATOL = (Int32)(3)
const KSP_CONVERGED_ITS = (Int32)(4)
const KSP_CONVERGED_CG_NEG_CURVE = (Int32)(5)
const KSP_CONVERGED_CG_CONSTRAINED = (Int32)(6)
const KSP_CONVERGED_STEP_LENGTH = (Int32)(7)
const KSP_CONVERGED_HAPPY_BREAKDOWN = (Int32)(8)
const KSP_DIVERGED_NULL = (Int32)(-2)
const KSP_DIVERGED_ITS = (Int32)(-3)
const KSP_DIVERGED_DTOL = (Int32)(-4)
const KSP_DIVERGED_BREAKDOWN = (Int32)(-5)
const KSP_DIVERGED_BREAKDOWN_BICG = (Int32)(-6)
const KSP_DIVERGED_NONSYMMETRIC = (Int32)(-7)
const KSP_DIVERGED_INDEFINITE_PC = (Int32)(-8)
const KSP_DIVERGED_NANORINF = (Int32)(-9)
const KSP_DIVERGED_INDEFINITE_MAT = (Int32)(-10)
const KSP_DIVERGED_PCSETUP_FAILED = (Int32)(-11)
const KSP_CONVERGED_ITERATING = (Int32)(0)

#= # end enum ANONYMOUS_63 =#
#= # begin enum KSPConvergedReason =#
typealias KSPConvergedReason Cint

const KSP_CONVERGED_RTOL_NORMAL = (Int32)(1)
const KSP_CONVERGED_ATOL_NORMAL = (Int32)(9)
const KSP_CONVERGED_RTOL = (Int32)(2)
const KSP_CONVERGED_ATOL = (Int32)(3)
const KSP_CONVERGED_ITS = (Int32)(4)
const KSP_CONVERGED_CG_NEG_CURVE = (Int32)(5)
const KSP_CONVERGED_CG_CONSTRAINED = (Int32)(6)
const KSP_CONVERGED_STEP_LENGTH = (Int32)(7)
const KSP_CONVERGED_HAPPY_BREAKDOWN = (Int32)(8)
const KSP_DIVERGED_NULL = (Int32)(-2)
const KSP_DIVERGED_ITS = (Int32)(-3)
const KSP_DIVERGED_DTOL = (Int32)(-4)
const KSP_DIVERGED_BREAKDOWN = (Int32)(-5)
const KSP_DIVERGED_BREAKDOWN_BICG = (Int32)(-6)
const KSP_DIVERGED_NONSYMMETRIC = (Int32)(-7)
const KSP_DIVERGED_INDEFINITE_PC = (Int32)(-8)
const KSP_DIVERGED_NANORINF = (Int32)(-9)
const KSP_DIVERGED_INDEFINITE_MAT = (Int32)(-10)
const KSP_DIVERGED_PCSETUP_FAILED = (Int32)(-11)
const KSP_CONVERGED_ITERATING = (Int32)(0)

#= # end enum KSPConvergedReason =#
#= # begin enum ANONYMOUS_64 =#
typealias ANONYMOUS_64 UInt32

const KSP_CG_SYMMETRIC = (UInt32)(0)
const KSP_CG_HERMITIAN = (UInt32)(1)

#= # end enum ANONYMOUS_64 =#
#= # begin enum KSPCGType =#
typealias KSPCGType UInt32

const KSP_CG_SYMMETRIC = (UInt32)(0)
const KSP_CG_HERMITIAN = (UInt32)(1)

#= # end enum KSPCGType =#
# skipping undefined typealias typealias KSPFischerGuess Ptr{_p_KSPFischerGuess}
#= # begin enum ANONYMOUS_65 =#
typealias ANONYMOUS_65 UInt32

const MAT_SCHUR_COMPLEMENT_AINV_DIAG = (UInt32)(0)
const MAT_SCHUR_COMPLEMENT_AINV_LUMP = (UInt32)(1)

#= # end enum ANONYMOUS_65 =#
#= # begin enum MatSchurComplementAinvType =#
typealias MatSchurComplementAinvType UInt32

const MAT_SCHUR_COMPLEMENT_AINV_DIAG = (UInt32)(0)
const MAT_SCHUR_COMPLEMENT_AINV_LUMP = (UInt32)(1)

#= # end enum MatSchurComplementAinvType =#
# skipping undefined typealias typealias SNES Ptr{_p_SNES}
typealias SNESType Symbol

#= # begin enum ANONYMOUS_66 =#
typealias ANONYMOUS_66 Cint

const SNES_CONVERGED_FNORM_ABS = (Int32)(2)
const SNES_CONVERGED_FNORM_RELATIVE = (Int32)(3)
const SNES_CONVERGED_SNORM_RELATIVE = (Int32)(4)
const SNES_CONVERGED_ITS = (Int32)(5)
const SNES_CONVERGED_TR_DELTA = (Int32)(7)
const SNES_DIVERGED_FUNCTION_DOMAIN = (Int32)(-1)
const SNES_DIVERGED_FUNCTION_COUNT = (Int32)(-2)
const SNES_DIVERGED_LINEAR_SOLVE = (Int32)(-3)
const SNES_DIVERGED_FNORM_NAN = (Int32)(-4)
const SNES_DIVERGED_MAX_IT = (Int32)(-5)
const SNES_DIVERGED_LINE_SEARCH = (Int32)(-6)
const SNES_DIVERGED_INNER = (Int32)(-7)
const SNES_DIVERGED_LOCAL_MIN = (Int32)(-8)
const SNES_CONVERGED_ITERATING = (Int32)(0)

#= # end enum ANONYMOUS_66 =#
#= # begin enum SNESConvergedReason =#
typealias SNESConvergedReason Cint

const SNES_CONVERGED_FNORM_ABS = (Int32)(2)
const SNES_CONVERGED_FNORM_RELATIVE = (Int32)(3)
const SNES_CONVERGED_SNORM_RELATIVE = (Int32)(4)
const SNES_CONVERGED_ITS = (Int32)(5)
const SNES_CONVERGED_TR_DELTA = (Int32)(7)
const SNES_DIVERGED_FUNCTION_DOMAIN = (Int32)(-1)
const SNES_DIVERGED_FUNCTION_COUNT = (Int32)(-2)
const SNES_DIVERGED_LINEAR_SOLVE = (Int32)(-3)
const SNES_DIVERGED_FNORM_NAN = (Int32)(-4)
const SNES_DIVERGED_MAX_IT = (Int32)(-5)
const SNES_DIVERGED_LINE_SEARCH = (Int32)(-6)
const SNES_DIVERGED_INNER = (Int32)(-7)
const SNES_DIVERGED_LOCAL_MIN = (Int32)(-8)
const SNES_CONVERGED_ITERATING = (Int32)(0)

#= # end enum SNESConvergedReason =#
#= # begin enum ANONYMOUS_67 =#
typealias ANONYMOUS_67 Cint

const SNES_NORM_DEFAULT = (Int32)(-1)
const SNES_NORM_NONE = (Int32)(0)
const SNES_NORM_ALWAYS = (Int32)(1)
const SNES_NORM_INITIAL_ONLY = (Int32)(2)
const SNES_NORM_FINAL_ONLY = (Int32)(3)
const SNES_NORM_INITIAL_FINAL_ONLY = (Int32)(4)

#= # end enum ANONYMOUS_67 =#
#= # begin enum SNESNormSchedule =#
typealias SNESNormSchedule Cint

const SNES_NORM_DEFAULT = (Int32)(-1)
const SNES_NORM_NONE = (Int32)(0)
const SNES_NORM_ALWAYS = (Int32)(1)
const SNES_NORM_INITIAL_ONLY = (Int32)(2)
const SNES_NORM_FINAL_ONLY = (Int32)(3)
const SNES_NORM_INITIAL_FINAL_ONLY = (Int32)(4)

#= # end enum SNESNormSchedule =#
#= # begin enum ANONYMOUS_68 =#
typealias ANONYMOUS_68 Cint

const SNES_FUNCTION_DEFAULT = (Int32)(-1)
const SNES_FUNCTION_UNPRECONDITIONED = (Int32)(0)
const SNES_FUNCTION_PRECONDITIONED = (Int32)(1)

#= # end enum ANONYMOUS_68 =#
#= # begin enum SNESFunctionType =#
typealias SNESFunctionType Cint

const SNES_FUNCTION_DEFAULT = (Int32)(-1)
const SNES_FUNCTION_UNPRECONDITIONED = (Int32)(0)
const SNES_FUNCTION_PRECONDITIONED = (Int32)(1)

#= # end enum SNESFunctionType =#
# skipping undefined typealias typealias SNESLineSearch Ptr{_p_LineSearch}
typealias SNESLineSearchType Symbol
typealias SNESLineSearchVIProjectFunc Ptr{Void}
typealias SNESLineSearchVINormFunc Ptr{Void}
typealias SNESLineSearchApplyFunc Ptr{Void}
typealias SNESLineSearchUserFunc Ptr{Void}

#= # begin enum ANONYMOUS_69 =#
typealias ANONYMOUS_69 UInt32

const SNES_LINESEARCH_SUCCEEDED = (UInt32)(0)
const SNES_LINESEARCH_FAILED_NANORINF = (UInt32)(1)
const SNES_LINESEARCH_FAILED_DOMAIN = (UInt32)(2)
const SNES_LINESEARCH_FAILED_REDUCT = (UInt32)(3)
const SNES_LINESEARCH_FAILED_USER = (UInt32)(4)
const SNES_LINESEARCH_FAILED_FUNCTION = (UInt32)(5)

#= # end enum ANONYMOUS_69 =#
#= # begin enum SNESLineSearchReason =#
typealias SNESLineSearchReason UInt32

const SNES_LINESEARCH_SUCCEEDED = (UInt32)(0)
const SNES_LINESEARCH_FAILED_NANORINF = (UInt32)(1)
const SNES_LINESEARCH_FAILED_DOMAIN = (UInt32)(2)
const SNES_LINESEARCH_FAILED_REDUCT = (UInt32)(3)
const SNES_LINESEARCH_FAILED_USER = (UInt32)(4)
const SNES_LINESEARCH_FAILED_FUNCTION = (UInt32)(5)

#= # end enum SNESLineSearchReason =#
typealias DMDASNESFunction Ptr{Void}
typealias DMDASNESJacobian Ptr{Void}
typealias DMDASNESObjective Ptr{Void}
typealias SNESMSType Symbol

#= # begin enum ANONYMOUS_70 =#
typealias ANONYMOUS_70 UInt32

const SNES_NGMRES_RESTART_NONE = (UInt32)(0)
const SNES_NGMRES_RESTART_PERIODIC = (UInt32)(1)
const SNES_NGMRES_RESTART_DIFFERENCE = (UInt32)(2)

#= # end enum ANONYMOUS_70 =#
#= # begin enum SNESNGMRESRestartType =#
typealias SNESNGMRESRestartType UInt32

const SNES_NGMRES_RESTART_NONE = (UInt32)(0)
const SNES_NGMRES_RESTART_PERIODIC = (UInt32)(1)
const SNES_NGMRES_RESTART_DIFFERENCE = (UInt32)(2)

#= # end enum SNESNGMRESRestartType =#
#= # begin enum ANONYMOUS_71 =#
typealias ANONYMOUS_71 UInt32

const SNES_NGMRES_SELECT_NONE = (UInt32)(0)
const SNES_NGMRES_SELECT_DIFFERENCE = (UInt32)(1)
const SNES_NGMRES_SELECT_LINESEARCH = (UInt32)(2)

#= # end enum ANONYMOUS_71 =#
#= # begin enum SNESNGMRESSelectType =#
typealias SNESNGMRESSelectType UInt32

const SNES_NGMRES_SELECT_NONE = (UInt32)(0)
const SNES_NGMRES_SELECT_DIFFERENCE = (UInt32)(1)
const SNES_NGMRES_SELECT_LINESEARCH = (UInt32)(2)

#= # end enum SNESNGMRESSelectType =#
#= # begin enum ANONYMOUS_72 =#
typealias ANONYMOUS_72 UInt32

const SNES_NCG_FR = (UInt32)(0)
const SNES_NCG_PRP = (UInt32)(1)
const SNES_NCG_HS = (UInt32)(2)
const SNES_NCG_DY = (UInt32)(3)
const SNES_NCG_CD = (UInt32)(4)

#= # end enum ANONYMOUS_72 =#
#= # begin enum SNESNCGType =#
typealias SNESNCGType UInt32

const SNES_NCG_FR = (UInt32)(0)
const SNES_NCG_PRP = (UInt32)(1)
const SNES_NCG_HS = (UInt32)(2)
const SNES_NCG_DY = (UInt32)(3)
const SNES_NCG_CD = (UInt32)(4)

#= # end enum SNESNCGType =#
#= # begin enum ANONYMOUS_73 =#
typealias ANONYMOUS_73 UInt32

const SNES_QN_SCALE_DEFAULT = (UInt32)(0)
const SNES_QN_SCALE_NONE = (UInt32)(1)
const SNES_QN_SCALE_SHANNO = (UInt32)(2)
const SNES_QN_SCALE_LINESEARCH = (UInt32)(3)
const SNES_QN_SCALE_JACOBIAN = (UInt32)(4)

#= # end enum ANONYMOUS_73 =#
#= # begin enum SNESQNScaleType =#
typealias SNESQNScaleType UInt32

const SNES_QN_SCALE_DEFAULT = (UInt32)(0)
const SNES_QN_SCALE_NONE = (UInt32)(1)
const SNES_QN_SCALE_SHANNO = (UInt32)(2)
const SNES_QN_SCALE_LINESEARCH = (UInt32)(3)
const SNES_QN_SCALE_JACOBIAN = (UInt32)(4)

#= # end enum SNESQNScaleType =#
#= # begin enum ANONYMOUS_74 =#
typealias ANONYMOUS_74 UInt32

const SNES_QN_RESTART_DEFAULT = (UInt32)(0)
const SNES_QN_RESTART_NONE = (UInt32)(1)
const SNES_QN_RESTART_POWELL = (UInt32)(2)
const SNES_QN_RESTART_PERIODIC = (UInt32)(3)

#= # end enum ANONYMOUS_74 =#
#= # begin enum SNESQNRestartType =#
typealias SNESQNRestartType UInt32

const SNES_QN_RESTART_DEFAULT = (UInt32)(0)
const SNES_QN_RESTART_NONE = (UInt32)(1)
const SNES_QN_RESTART_POWELL = (UInt32)(2)
const SNES_QN_RESTART_PERIODIC = (UInt32)(3)

#= # end enum SNESQNRestartType =#
#= # begin enum ANONYMOUS_75 =#
typealias ANONYMOUS_75 UInt32

const SNES_QN_LBFGS = (UInt32)(0)
const SNES_QN_BROYDEN = (UInt32)(1)
const SNES_QN_BADBROYDEN = (UInt32)(2)

#= # end enum ANONYMOUS_75 =#
#= # begin enum SNESQNType =#
typealias SNESQNType UInt32

const SNES_QN_LBFGS = (UInt32)(0)
const SNES_QN_BROYDEN = (UInt32)(1)
const SNES_QN_BADBROYDEN = (UInt32)(2)

#= # end enum SNESQNType =#
#= # begin enum ANONYMOUS_76 =#
typealias ANONYMOUS_76 UInt32

const SNES_COMPOSITE_ADDITIVE = (UInt32)(0)
const SNES_COMPOSITE_MULTIPLICATIVE = (UInt32)(1)
const SNES_COMPOSITE_ADDITIVEOPTIMAL = (UInt32)(2)

#= # end enum ANONYMOUS_76 =#
#= # begin enum SNESCompositeType =#
typealias SNESCompositeType UInt32

const SNES_COMPOSITE_ADDITIVE = (UInt32)(0)
const SNES_COMPOSITE_MULTIPLICATIVE = (UInt32)(1)
const SNES_COMPOSITE_ADDITIVEOPTIMAL = (UInt32)(2)

#= # end enum SNESCompositeType =#
#= # begin enum ANONYMOUS_77 =#
typealias ANONYMOUS_77 UInt32

const SNES_FAS_MULTIPLICATIVE = (UInt32)(0)
const SNES_FAS_ADDITIVE = (UInt32)(1)
const SNES_FAS_FULL = (UInt32)(2)
const SNES_FAS_KASKADE = (UInt32)(3)

#= # end enum ANONYMOUS_77 =#
#= # begin enum SNESFASType =#
typealias SNESFASType UInt32

const SNES_FAS_MULTIPLICATIVE = (UInt32)(0)
const SNES_FAS_ADDITIVE = (UInt32)(1)
const SNES_FAS_FULL = (UInt32)(2)
const SNES_FAS_KASKADE = (UInt32)(3)

#= # end enum SNESFASType =#
# skipping undefined typealias typealias TS Ptr{_p_TS}
typealias TSType Symbol

#= # begin enum ANONYMOUS_78 =#
typealias ANONYMOUS_78 UInt32

const TS_LINEAR = (UInt32)(0)
const TS_NONLINEAR = (UInt32)(1)

#= # end enum ANONYMOUS_78 =#
#= # begin enum TSProblemType =#
typealias TSProblemType UInt32

const TS_LINEAR = (UInt32)(0)
const TS_NONLINEAR = (UInt32)(1)

#= # end enum TSProblemType =#
#= # begin enum ANONYMOUS_79 =#
typealias ANONYMOUS_79 Cint

const TS_EQ_UNSPECIFIED = (Int32)(-1)
const TS_EQ_EXPLICIT = (Int32)(0)
const TS_EQ_ODE_EXPLICIT = (Int32)(1)
const TS_EQ_DAE_SEMI_EXPLICIT_INDEX1 = (Int32)(100)
const TS_EQ_DAE_SEMI_EXPLICIT_INDEX2 = (Int32)(200)
const TS_EQ_DAE_SEMI_EXPLICIT_INDEX3 = (Int32)(300)
const TS_EQ_DAE_SEMI_EXPLICIT_INDEXHI = (Int32)(500)
const TS_EQ_IMPLICIT = (Int32)(1000)
const TS_EQ_ODE_IMPLICIT = (Int32)(1001)
const TS_EQ_DAE_IMPLICIT_INDEX1 = (Int32)(1100)
const TS_EQ_DAE_IMPLICIT_INDEX2 = (Int32)(1200)
const TS_EQ_DAE_IMPLICIT_INDEX3 = (Int32)(1300)
const TS_EQ_DAE_IMPLICIT_INDEXHI = (Int32)(1500)

#= # end enum ANONYMOUS_79 =#
#= # begin enum TSEquationType =#
typealias TSEquationType Cint

const TS_EQ_UNSPECIFIED = (Int32)(-1)
const TS_EQ_EXPLICIT = (Int32)(0)
const TS_EQ_ODE_EXPLICIT = (Int32)(1)
const TS_EQ_DAE_SEMI_EXPLICIT_INDEX1 = (Int32)(100)
const TS_EQ_DAE_SEMI_EXPLICIT_INDEX2 = (Int32)(200)
const TS_EQ_DAE_SEMI_EXPLICIT_INDEX3 = (Int32)(300)
const TS_EQ_DAE_SEMI_EXPLICIT_INDEXHI = (Int32)(500)
const TS_EQ_IMPLICIT = (Int32)(1000)
const TS_EQ_ODE_IMPLICIT = (Int32)(1001)
const TS_EQ_DAE_IMPLICIT_INDEX1 = (Int32)(1100)
const TS_EQ_DAE_IMPLICIT_INDEX2 = (Int32)(1200)
const TS_EQ_DAE_IMPLICIT_INDEX3 = (Int32)(1300)
const TS_EQ_DAE_IMPLICIT_INDEXHI = (Int32)(1500)

#= # end enum TSEquationType =#
#= # begin enum ANONYMOUS_80 =#
typealias ANONYMOUS_80 Cint

const TS_CONVERGED_ITERATING = (Int32)(0)
const TS_CONVERGED_TIME = (Int32)(1)
const TS_CONVERGED_ITS = (Int32)(2)
const TS_CONVERGED_USER = (Int32)(3)
const TS_CONVERGED_EVENT = (Int32)(4)
const TS_DIVERGED_NONLINEAR_SOLVE = (Int32)(-1)
const TS_DIVERGED_STEP_REJECTED = (Int32)(-2)

#= # end enum ANONYMOUS_80 =#
#= # begin enum TSConvergedReason =#
typealias TSConvergedReason Cint

const TS_CONVERGED_ITERATING = (Int32)(0)
const TS_CONVERGED_TIME = (Int32)(1)
const TS_CONVERGED_ITS = (Int32)(2)
const TS_CONVERGED_USER = (Int32)(3)
const TS_CONVERGED_EVENT = (Int32)(4)
const TS_DIVERGED_NONLINEAR_SOLVE = (Int32)(-1)
const TS_DIVERGED_STEP_REJECTED = (Int32)(-2)

#= # end enum TSConvergedReason =#
#= # begin enum ANONYMOUS_81 =#
typealias ANONYMOUS_81 UInt32

const TS_EXACTFINALTIME_STEPOVER = (UInt32)(0)
const TS_EXACTFINALTIME_INTERPOLATE = (UInt32)(1)
const TS_EXACTFINALTIME_MATCHSTEP = (UInt32)(2)

#= # end enum ANONYMOUS_81 =#
#= # begin enum TSExactFinalTimeOption =#
typealias TSExactFinalTimeOption UInt32

const TS_EXACTFINALTIME_STEPOVER = (UInt32)(0)
const TS_EXACTFINALTIME_INTERPOLATE = (UInt32)(1)
const TS_EXACTFINALTIME_MATCHSTEP = (UInt32)(2)

#= # end enum TSExactFinalTimeOption =#
# skipping undefined typealias typealias TSTrajectory Ptr{_p_TSTrajectory}
typealias TSTrajectoryType Symbol

# skipping undefined typealias typealias TSMonitorDrawCtx Ptr{_n_TSMonitorDrawCtx}
typealias TSRHSFunction Ptr{Void}
typealias TSRHSJacobian Ptr{Void}
typealias TSSolutionFunction Ptr{Void}
typealias TSIFunction Ptr{Void}
typealias TSIJacobian Ptr{Void}
typealias DMDATSRHSFunctionLocal Ptr{Void}
typealias DMDATSRHSJacobianLocal Ptr{Void}
typealias DMDATSIFunctionLocal Ptr{Void}
typealias DMDATSIJacobianLocal Ptr{Void}

# skipping undefined typealias typealias TSMonitorLGCtx Ptr{_n_TSMonitorLGCtx}
#= skipping type declaration with undefined symbols:
immutable TSMonitorDMDARayCtx
    ray::Vec
    scatter::VecScatter
    viewer::PetscViewer
    lgctx::TSMonitorLGCtx
end 
=#
# skipping undefined typealias typealias TSMonitorEnvelopeCtx Ptr{_n_TSMonitorEnvelopeCtx}
# skipping undefined typealias typealias TSMonitorSPEigCtx Ptr{_n_TSMonitorSPEigCtx}
typealias TSSSPType Symbol

# skipping undefined typealias typealias TSAdapt Ptr{_p_TSAdapt}
typealias TSAdaptType Symbol

# skipping undefined typealias typealias TSGLAdapt Ptr{_p_TSGLAdapt}
typealias TSGLAdaptType Symbol
typealias TSGLAcceptType Symbol
typealias TSGLAcceptFunction Ptr{Void}
typealias TSGLType Symbol
typealias TSRKType Symbol
typealias TSARKIMEXType Symbol
typealias TSRosWType Symbol

#= # begin enum ANONYMOUS_82 =#
typealias ANONYMOUS_82 UInt32

const TAO_SUBSET_SUBVEC = (UInt32)(0)
const TAO_SUBSET_MASK = (UInt32)(1)
const TAO_SUBSET_MATRIXFREE = (UInt32)(2)

#= # end enum ANONYMOUS_82 =#
#= # begin enum TaoSubsetType =#
typealias TaoSubsetType UInt32

const TAO_SUBSET_SUBVEC = (UInt32)(0)
const TAO_SUBSET_MASK = (UInt32)(1)
const TAO_SUBSET_MATRIXFREE = (UInt32)(2)

#= # end enum TaoSubsetType =#
# skipping undefined typealias typealias Tao Ptr{_p_Tao}
#= # begin enum ANONYMOUS_83 =#
typealias ANONYMOUS_83 Cint

const TAO_CONVERGED_FATOL = (Int32)(1)
const TAO_CONVERGED_FRTOL = (Int32)(2)
const TAO_CONVERGED_GATOL = (Int32)(3)
const TAO_CONVERGED_GRTOL = (Int32)(4)
const TAO_CONVERGED_GTTOL = (Int32)(5)
const TAO_CONVERGED_STEPTOL = (Int32)(6)
const TAO_CONVERGED_MINF = (Int32)(7)
const TAO_CONVERGED_USER = (Int32)(8)
const TAO_DIVERGED_MAXITS = (Int32)(-2)
const TAO_DIVERGED_NAN = (Int32)(-4)
const TAO_DIVERGED_MAXFCN = (Int32)(-5)
const TAO_DIVERGED_LS_FAILURE = (Int32)(-6)
const TAO_DIVERGED_TR_REDUCTION = (Int32)(-7)
const TAO_DIVERGED_USER = (Int32)(-8)
const TAO_CONTINUE_ITERATING = (Int32)(0)

#= # end enum ANONYMOUS_83 =#
#= # begin enum TaoConvergedReason =#
typealias TaoConvergedReason Cint

const TAO_CONVERGED_FATOL = (Int32)(1)
const TAO_CONVERGED_FRTOL = (Int32)(2)
const TAO_CONVERGED_GATOL = (Int32)(3)
const TAO_CONVERGED_GRTOL = (Int32)(4)
const TAO_CONVERGED_GTTOL = (Int32)(5)
const TAO_CONVERGED_STEPTOL = (Int32)(6)
const TAO_CONVERGED_MINF = (Int32)(7)
const TAO_CONVERGED_USER = (Int32)(8)
const TAO_DIVERGED_MAXITS = (Int32)(-2)
const TAO_DIVERGED_NAN = (Int32)(-4)
const TAO_DIVERGED_MAXFCN = (Int32)(-5)
const TAO_DIVERGED_LS_FAILURE = (Int32)(-6)
const TAO_DIVERGED_TR_REDUCTION = (Int32)(-7)
const TAO_DIVERGED_USER = (Int32)(-8)
const TAO_CONTINUE_ITERATING = (Int32)(0)

#= # end enum TaoConvergedReason =#
# skipping undefined typealias typealias TaoLineSearch Ptr{_p_TaoLineSearch}
#= # begin enum ANONYMOUS_84 =#
typealias ANONYMOUS_84 Cint

const TAOLINESEARCH_FAILED_INFORNAN = (Int32)(-1)
const TAOLINESEARCH_FAILED_BADPARAMETER = (Int32)(-2)
const TAOLINESEARCH_FAILED_ASCENT = (Int32)(-3)
const TAOLINESEARCH_CONTINUE_ITERATING = (Int32)(0)
const TAOLINESEARCH_SUCCESS = (Int32)(1)
const TAOLINESEARCH_SUCCESS_USER = (Int32)(2)
const TAOLINESEARCH_HALTED_OTHER = (Int32)(3)
const TAOLINESEARCH_HALTED_MAXFCN = (Int32)(4)
const TAOLINESEARCH_HALTED_UPPERBOUND = (Int32)(5)
const TAOLINESEARCH_HALTED_LOWERBOUND = (Int32)(6)
const TAOLINESEARCH_HALTED_RTOL = (Int32)(7)
const TAOLINESEARCH_HALTED_USER = (Int32)(8)

#= # end enum ANONYMOUS_84 =#
#= # begin enum TaoLineSearchConvergedReason =#
typealias TaoLineSearchConvergedReason Cint

const TAOLINESEARCH_FAILED_INFORNAN = (Int32)(-1)
const TAOLINESEARCH_FAILED_BADPARAMETER = (Int32)(-2)
const TAOLINESEARCH_FAILED_ASCENT = (Int32)(-3)
const TAOLINESEARCH_CONTINUE_ITERATING = (Int32)(0)
const TAOLINESEARCH_SUCCESS = (Int32)(1)
const TAOLINESEARCH_SUCCESS_USER = (Int32)(2)
const TAOLINESEARCH_HALTED_OTHER = (Int32)(3)
const TAOLINESEARCH_HALTED_MAXFCN = (Int32)(4)
const TAOLINESEARCH_HALTED_UPPERBOUND = (Int32)(5)
const TAOLINESEARCH_HALTED_LOWERBOUND = (Int32)(6)
const TAOLINESEARCH_HALTED_RTOL = (Int32)(7)
const TAOLINESEARCH_HALTED_USER = (Int32)(8)

#= # end enum TaoLineSearchConvergedReason =#
