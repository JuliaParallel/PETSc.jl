# Automatically generated using Clang.jl wrap_c, version 0.0.0

using Compat

const global PETSC_FUNCTION_NAME = PETSC_FUNCTION_NAME_C
const global PETSC_RESTRICT = PETSC_C_RESTRICT
const global PETSC_STATIC_INLINE = PETSC_C_STATIC_INLINE
const global PETSC_VISIBILITY_PUBLIC = PETSC_DLLIMPORT

# Skipping MacroDefinition: PETSC_EXTERN extern PETSC_VISIBILITY_PUBLIC
# Skipping MacroDefinition: PETSC_INTERN extern PETSC_VISIBILITY_INTERNAL

const global PETSC_VERSION_RELEASE = 1
const global PETSC_VERSION_MAJOR = 3
const global PETSC_VERSION_MINOR = 6
const global PETSC_VERSION_SUBMINOR = 0
const global PETSC_VERSION_PATCH = 0
const global PETSC_RELEASE_DATE = "Jun, 9, 2015"
const global PETSC_VERSION_DATE = "Jun, 09, 2015"
const global PETSC_VERSION_GIT = "v3.6"
const global PETSC_VERSION_DATE_GIT = "2015-06-09 16:15:46 -0500"

# Skipping MacroDefinition: PETSC_VERSION_ ( MAJOR , MINOR , SUBMINOR ) ( ( PETSC_VERSION_MAJOR == ( MAJOR ) ) && ( PETSC_VERSION_MINOR == ( MINOR ) ) && ( PETSC_VERSION_SUBMINOR == ( SUBMINOR ) ) && ( PETSC_VERSION_RELEASE == 1 ) )
# Skipping MacroDefinition: PETSC_VERSION_LT ( MAJOR , MINOR , SUBMINOR ) ( PETSC_VERSION_RELEASE == 1 && ( PETSC_VERSION_MAJOR < ( MAJOR ) || ( PETSC_VERSION_MAJOR == ( MAJOR ) && ( PETSC_VERSION_MINOR < ( MINOR ) || ( PETSC_VERSION_MINOR == ( MINOR ) && ( PETSC_VERSION_SUBMINOR < ( SUBMINOR ) ) ) ) ) ) )
# Skipping MacroDefinition: PETSC_VERSION_LE ( MAJOR , MINOR , SUBMINOR ) ( PETSC_VERSION_LT ( MAJOR , MINOR , SUBMINOR ) || PETSC_VERSION_ ( MAJOR , MINOR , SUBMINOR ) )
# Skipping MacroDefinition: PETSC_VERSION_GT ( MAJOR , MINOR , SUBMINOR ) ( 0 == PETSC_VERSION_LE ( MAJOR , MINOR , SUBMINOR ) )
# Skipping MacroDefinition: PETSC_VERSION_GE ( MAJOR , MINOR , SUBMINOR ) ( 0 == PETSC_VERSION_LT ( MAJOR , MINOR , SUBMINOR ) )

const global PETSC_AUTHOR_INFO = "       The PETSc Team\n    petsc-maint@mcs.anl.gov\n http://www.mcs.anl.gov/petsc/\n"
const global MPICH_SKIP_MPICXX = 1
const global OMPI_SKIP_MPICXX = 1

# Skipping MacroDefinition: PetscAttrMPIPointerWithType ( bufno , typeno ) __attribute__ ( ( pointer_with_type_tag ( MPI , bufno , typeno ) ) )
# Skipping MacroDefinition: PetscAttrMPITypeTag ( type ) __attribute__ ( ( type_tag_for_datatype ( MPI , type ) ) )
# Skipping MacroDefinition: PetscAttrMPITypeTagLayoutCompatible ( type ) __attribute__ ( ( type_tag_for_datatype ( MPI , type , layout_compatible ) ) )

const global MPIU_INT = MPI_INT
const global MPIU_INT64 = MPI_LONG_LONG_INT
const global MPIU_SIZE_T = MPI_UNSIGNED

# Skipping MacroDefinition: PetscUnlikely ( cond ) ( cond )
# Skipping MacroDefinition: PetscLikely ( cond ) ( cond )
# Skipping MacroDefinition: PetscExpPassiveScalar ( a ) PetscExpScalar ( )

const global MPIU_SCALAR = MPIU_REAL

# Skipping MacroDefinition: PetscRealPart ( a ) ( a )
# Skipping MacroDefinition: PetscImaginaryPart ( a ) ( ( PetscReal ) 0. )
# Skipping MacroDefinition: PetscConj ( a ) ( a )
# Skipping MacroDefinition: PetscSqrtScalar ( a ) sqrt ( a )
# Skipping MacroDefinition: PetscPowScalar ( a , b ) pow ( a , b )
# Skipping MacroDefinition: PetscExpScalar ( a ) exp ( a )
# Skipping MacroDefinition: PetscLogScalar ( a ) log ( a )
# Skipping MacroDefinition: PetscSinScalar ( a ) sin ( a )
# Skipping MacroDefinition: PetscCosScalar ( a ) cos ( a )
# Skipping MacroDefinition: PetscAsinScalar ( a ) asin ( a )
# Skipping MacroDefinition: PetscAcosScalar ( a ) acos ( a )
# Skipping MacroDefinition: PetscTanScalar ( a ) tan ( a )
# Skipping MacroDefinition: PetscSinhScalar ( a ) sinh ( a )
# Skipping MacroDefinition: PetscCoshScalar ( a ) cosh ( a )
# Skipping MacroDefinition: PetscTanhScalar ( a ) tanh ( a )
# Skipping MacroDefinition: PetscSign ( a ) ( ( ( a ) >= 0 ) ? ( ( a ) == 0 ? 0 : 1 ) : - 1 )
# Skipping MacroDefinition: PetscAbs ( a ) ( ( ( a ) >= 0 ) ? ( a ) : - ( a ) )
# Skipping MacroDefinition: PetscMin ( a , b ) ( ( ( a ) < ( b ) ) ? ( a ) : ( b ) )
# Skipping MacroDefinition: PetscMax ( a , b ) ( ( ( a ) < ( b ) ) ? ( b ) : ( a ) )
# Skipping MacroDefinition: PetscClipInterval ( x , a , b ) ( PetscMax ( ( a ) , PetscMin ( ( x ) , ( b ) ) ) )
# Skipping MacroDefinition: PetscAbsInt ( a ) ( ( ( a ) < 0 ) ? - ( a ) : ( a ) )
# Skipping MacroDefinition: PetscAbsReal ( a ) ( ( ( a ) < 0 ) ? - ( a ) : ( a ) )
# Skipping MacroDefinition: PetscSqr ( a ) ( ( a ) * ( a ) )

const global PETSC_PI = M_PI
const global PETSC_MAX_INT = 2147483647
const global PETSC_MIN_INT = -PETSC_MAX_INT - 1
const global PETSC_INFINITY = PETSC_MAX_REAL / 4.0
const global PETSC_NINFINITY = -PETSC_INFINITY
const global PassiveReal = PetscReal

typealias PetscScalar Cint

const global PassiveScalar = PetscScalar
const global MPIU_MATSCALAR = MPIU_SCALAR
const global MPIU_2INT = MPI_2INT
const global PETSC_NULL = NULL
const global PETSC_IGNORE = NULL
const global PETSC_DECIDE = -1
const global PETSC_DETERMINE = PETSC_DECIDE
const global PETSC_DEFAULT = -2
const global PETSC_COMM_SELF = MPI_COMM_SELF

# Skipping MacroDefinition: PetscMalloc ( a , b ) ( ( * PetscTrMalloc ) ( ( a ) , __LINE__ , PETSC_FUNCTION_NAME , __FILE__ , ( void * * ) ( b ) ) )
# Skipping MacroDefinition: PetscAddrAlign ( a ) ( void * ) ( ( ( ( PETSC_UINTPTR_T ) ( a ) ) + ( PETSC_MEMALIGN - 1 ) ) & ~ ( PETSC_MEMALIGN - 1 ) )
# Skipping MacroDefinition: PetscMalloc1 ( m1 , r1 ) PetscMalloc ( ( m1 ) * sizeof ( * * ( r1 ) ) , r1 )
# Skipping MacroDefinition: PetscCalloc1 ( m1 , r1 ) ( PetscMalloc1 ( ( m1 ) , r1 ) || PetscMemzero ( * ( r1 ) , ( m1 ) * sizeof ( * * ( r1 ) ) ) )
# Skipping MacroDefinition: PetscMalloc2 ( m1 , r1 , m2 , r2 ) ( PetscMalloc1 ( ( m1 ) , ( r1 ) ) || PetscMalloc1 ( ( m2 ) , ( r2 ) ) )
# Skipping MacroDefinition: PetscCalloc2 ( m1 , r1 , m2 , r2 ) ( PetscMalloc2 ( ( m1 ) , ( r1 ) , ( m2 ) , ( r2 ) ) || PetscMemzero ( * ( r1 ) , ( m1 ) * sizeof ( * * ( r1 ) ) ) || PetscMemzero ( * ( r2 ) , ( m2 ) * sizeof ( * * ( r2 ) ) ) )
# Skipping MacroDefinition: PetscMalloc3 ( m1 , r1 , m2 , r2 , m3 , r3 ) ( PetscMalloc1 ( ( m1 ) , ( r1 ) ) || PetscMalloc1 ( ( m2 ) , ( r2 ) ) || PetscMalloc1 ( ( m3 ) , ( r3 ) ) )
# Skipping MacroDefinition: PetscCalloc3 ( m1 , r1 , m2 , r2 , m3 , r3 ) ( PetscMalloc3 ( ( m1 ) , ( r1 ) , ( m2 ) , ( r2 ) , ( m3 ) , ( r3 ) ) || PetscMemzero ( * ( r1 ) , ( m1 ) * sizeof ( * * ( r1 ) ) ) || PetscMemzero ( * ( r2 ) , ( m2 ) * sizeof ( * * ( r2 ) ) ) || PetscMemzero ( * ( r3 ) , ( m3 ) * sizeof ( * * ( r3 ) ) ) )
# Skipping MacroDefinition: PetscMalloc4 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 ) ( PetscMalloc1 ( ( m1 ) , ( r1 ) ) || PetscMalloc1 ( ( m2 ) , ( r2 ) ) || PetscMalloc1 ( ( m3 ) , ( r3 ) ) || PetscMalloc1 ( ( m4 ) , ( r4 ) ) )
# Skipping MacroDefinition: PetscCalloc4 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 ) ( PetscMalloc4 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 ) || PetscMemzero ( * ( r1 ) , ( m1 ) * sizeof ( * * ( r1 ) ) ) || PetscMemzero ( * ( r2 ) , ( m2 ) * sizeof ( * * ( r2 ) ) ) || PetscMemzero ( * ( r3 ) , ( m3 ) * sizeof ( * * ( r3 ) ) ) || PetscMemzero ( * ( r4 ) , ( m4 ) * sizeof ( * * ( r4 ) ) ) )
# Skipping MacroDefinition: PetscMalloc5 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 , m5 , r5 ) ( PetscMalloc1 ( ( m1 ) , ( r1 ) ) || PetscMalloc1 ( ( m2 ) , ( r2 ) ) || PetscMalloc1 ( ( m3 ) , ( r3 ) ) || PetscMalloc1 ( ( m4 ) , ( r4 ) ) || PetscMalloc1 ( ( m5 ) , ( r5 ) ) )
# Skipping MacroDefinition: PetscCalloc5 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 , m5 , r5 ) ( PetscMalloc5 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 , m5 , r5 ) || PetscMemzero ( * ( r1 ) , ( m1 ) * sizeof ( * * ( r1 ) ) ) || PetscMemzero ( * ( r2 ) , ( m2 ) * sizeof ( * * ( r2 ) ) ) || PetscMemzero ( * ( r3 ) , ( m3 ) * sizeof ( * * ( r3 ) ) ) || PetscMemzero ( * ( r4 ) , ( m4 ) * sizeof ( * * ( r4 ) ) ) || PetscMemzero ( * ( r5 ) , ( m5 ) * sizeof ( * * ( r5 ) ) ) )
# Skipping MacroDefinition: PetscMalloc6 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 , m5 , r5 , m6 , r6 ) ( PetscMalloc1 ( ( m1 ) , ( r1 ) ) || PetscMalloc1 ( ( m2 ) , ( r2 ) ) || PetscMalloc1 ( ( m3 ) , ( r3 ) ) || PetscMalloc1 ( ( m4 ) , ( r4 ) ) || PetscMalloc1 ( ( m5 ) , ( r5 ) ) || PetscMalloc1 ( ( m6 ) , ( r6 ) ) )
# Skipping MacroDefinition: PetscCalloc6 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 , m5 , r5 , m6 , r6 ) ( PetscMalloc6 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 , m5 , r5 , m6 , r6 ) || PetscMemzero ( * ( r1 ) , ( m1 ) * sizeof ( * * ( r1 ) ) ) || PetscMemzero ( * ( r2 ) , ( m2 ) * sizeof ( * * ( r2 ) ) ) || PetscMemzero ( * ( r3 ) , ( m3 ) * sizeof ( * * ( r3 ) ) ) || PetscMemzero ( * ( r4 ) , ( m4 ) * sizeof ( * * ( r4 ) ) ) || PetscMemzero ( * ( r5 ) , ( m5 ) * sizeof ( * * ( r5 ) ) ) || PetscMemzero ( * ( r6 ) , ( m6 ) * sizeof ( * * ( r6 ) ) ) )
# Skipping MacroDefinition: PetscMalloc7 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 , m5 , r5 , m6 , r6 , m7 , r7 ) ( PetscMalloc1 ( ( m1 ) , ( r1 ) ) || PetscMalloc1 ( ( m2 ) , ( r2 ) ) || PetscMalloc1 ( ( m3 ) , ( r3 ) ) || PetscMalloc1 ( ( m4 ) , ( r4 ) ) || PetscMalloc1 ( ( m5 ) , ( r5 ) ) || PetscMalloc1 ( ( m6 ) , ( r6 ) ) || PetscMalloc1 ( ( m7 ) , ( r7 ) ) )
# Skipping MacroDefinition: PetscCalloc7 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 , m5 , r5 , m6 , r6 , m7 , r7 ) ( PetscMalloc7 ( m1 , r1 , m2 , r2 , m3 , r3 , m4 , r4 , m5 , r5 , m6 , r6 , m7 , r7 ) || PetscMemzero ( * ( r1 ) , ( m1 ) * sizeof ( * * ( r1 ) ) ) || PetscMemzero ( * ( r2 ) , ( m2 ) * sizeof ( * * ( r2 ) ) ) || PetscMemzero ( * ( r3 ) , ( m3 ) * sizeof ( * * ( r3 ) ) ) || PetscMemzero ( * ( r4 ) , ( m4 ) * sizeof ( * * ( r4 ) ) ) || PetscMemzero ( * ( r5 ) , ( m5 ) * sizeof ( * * ( r5 ) ) ) || PetscMemzero ( * ( r6 ) , ( m6 ) * sizeof ( * * ( r6 ) ) ) || PetscMemzero ( * ( r7 ) , ( m7 ) * sizeof ( * * ( r7 ) ) ) )
# Skipping MacroDefinition: PetscNew ( b ) PetscCalloc1 ( 1 , ( b ) )
# Skipping MacroDefinition: PetscNewLog ( o , b ) ( PetscNew ( ( b ) ) || PetscLogObjectMemory ( ( PetscObject ) o , sizeof ( * * ( b ) ) ) )
# Skipping MacroDefinition: PetscFree ( a ) ( ( * PetscTrFree ) ( ( void * ) ( a ) , __LINE__ , PETSC_FUNCTION_NAME , __FILE__ ) || ( ( a ) = 0 , 0 ) )
# Skipping MacroDefinition: PetscFreeVoid ( a ) ( ( * PetscTrFree ) ( ( a ) , __LINE__ , PETSC_FUNCTION_NAME , __FILE__ ) , ( a ) = 0 )
# Skipping MacroDefinition: PetscFree2 ( m1 , m2 ) ( PetscFree ( m2 ) || PetscFree ( m1 ) )
# Skipping MacroDefinition: PetscFree3 ( m1 , m2 , m3 ) ( PetscFree ( m3 ) || PetscFree ( m2 ) || PetscFree ( m1 ) )
# Skipping MacroDefinition: PetscFree4 ( m1 , m2 , m3 , m4 ) ( PetscFree ( m4 ) || PetscFree ( m3 ) || PetscFree ( m2 ) || PetscFree ( m1 ) )
# Skipping MacroDefinition: PetscFree5 ( m1 , m2 , m3 , m4 , m5 ) ( PetscFree ( m5 ) || PetscFree ( m4 ) || PetscFree ( m3 ) || PetscFree ( m2 ) || PetscFree ( m1 ) )
# Skipping MacroDefinition: PetscFree6 ( m1 , m2 , m3 , m4 , m5 , m6 ) ( PetscFree ( m6 ) || PetscFree ( m5 ) || PetscFree ( m4 ) || PetscFree ( m3 ) || PetscFree ( m2 ) || PetscFree ( m1 ) )
# Skipping MacroDefinition: PetscFree7 ( m1 , m2 , m3 , m4 , m5 , m6 , m7 ) ( PetscFree ( m7 ) || PetscFree ( m6 ) || PetscFree ( m5 ) || PetscFree ( m4 ) || PetscFree ( m3 ) || PetscFree ( m2 ) || PetscFree ( m1 ) )

const global MPIU_PETSCLOGDOUBLE = MPI_DOUBLE

# begin enum PetscDataType
typealias PetscDataType Uint32
const global PETSC_INT = (UInt32)(0)
const global PETSC_DOUBLE = (UInt32)(1)
const global PETSC_COMPLEX = (UInt32)(2)
const global PETSC_LONG = (UInt32)(3)
const global PETSC_SHORT = (UInt32)(4)
const global PETSC_FLOAT = (UInt32)(5)
const global PETSC_CHAR = (UInt32)(6)
const global PETSC_BIT_LOGICAL = (UInt32)(7)
const global PETSC_ENUM = (UInt32)(8)
const global PETSC_BOOL = (UInt32)(9)
const global PETSC___FLOAT128 = (UInt32)(10)
const global PETSC_OBJECT = (UInt32)(11)
const global PETSC_FUNCTION = (UInt32)(12)
const global PETSC_STRING = (UInt32)(12)
# end enum PetscDataType

const global PETSC_SCALAR = PETSC_DOUBLE
const global PETSC_REAL = PETSC_DOUBLE
const global PETSC_FORTRANADDR = PETSC_LONG
const global MPIU_SUM = MPI_SUM
const global MPIU_MAX = MPI_MAX
const global MPIU_MIN = MPI_MIN
const global PETSC_ERR_MIN_VALUE = 54
const global PETSC_ERR_MEM = 55
const global PETSC_ERR_SUP = 56
const global PETSC_ERR_SUP_SYS = 57
const global PETSC_ERR_ORDER = 58
const global PETSC_ERR_SIG = 59
const global PETSC_ERR_FP = 72
const global PETSC_ERR_COR = 74
const global PETSC_ERR_LIB = 76
const global PETSC_ERR_PLIB = 77
const global PETSC_ERR_MEMC = 78
const global PETSC_ERR_CONV_FAILED = 82
const global PETSC_ERR_USER = 83
const global PETSC_ERR_SYS = 88
const global PETSC_ERR_POINTER = 70
const global PETSC_ERR_ARG_SIZ = 60
const global PETSC_ERR_ARG_IDN = 61
const global PETSC_ERR_ARG_WRONG = 62
const global PETSC_ERR_ARG_CORRUPT = 64
const global PETSC_ERR_ARG_OUTOFRANGE = 63
const global PETSC_ERR_ARG_BADPTR = 68
const global PETSC_ERR_ARG_NOTSAMETYPE = 69
const global PETSC_ERR_ARG_NOTSAMECOMM = 80
const global PETSC_ERR_ARG_WRONGSTATE = 73
const global PETSC_ERR_ARG_TYPENOTSET = 89
const global PETSC_ERR_ARG_INCOMP = 75
const global PETSC_ERR_ARG_NULL = 85
const global PETSC_ERR_ARG_UNKNOWN_TYPE = 86
const global PETSC_ERR_FILE_OPEN = 65
const global PETSC_ERR_FILE_READ = 66
const global PETSC_ERR_FILE_WRITE = 67
const global PETSC_ERR_FILE_UNEXPECTED = 79
const global PETSC_ERR_MAT_LU_ZRPVT = 71
const global PETSC_ERR_MAT_CH_ZRPVT = 81
const global PETSC_ERR_INT_OVERFLOW = 84
const global PETSC_ERR_FLOP_COUNT = 90
const global PETSC_ERR_NOT_CONVERGED = 91
const global PETSC_ERR_MISSING_FACTOR = 92
const global PETSC_ERR_MAX_VALUE = 93

# Skipping MacroDefinition: PetscStringizeArg ( a ) # a
# Skipping MacroDefinition: PetscStringize ( a ) PetscStringizeArg ( a )
# Skipping MacroDefinition: SETERRQ ( c , n , s )
# Skipping MacroDefinition: SETERRQ1 ( c , n , s , a1 )
# Skipping MacroDefinition: SETERRQ2 ( c , n , s , a1 , a2 )
# Skipping MacroDefinition: SETERRQ3 ( c , n , s , a1 , a2 , a3 )
# Skipping MacroDefinition: SETERRQ4 ( c , n , s , a1 , a2 , a3 , a4 )
# Skipping MacroDefinition: SETERRQ5 ( c , n , s , a1 , a2 , a3 , a4 , a5 )
# Skipping MacroDefinition: SETERRQ6 ( c , n , s , a1 , a2 , a3 , a4 , a5 , a6 )
# Skipping MacroDefinition: SETERRQ7 ( c , n , s , a1 , a2 , a3 , a4 , a5 , a6 , a7 )
# Skipping MacroDefinition: SETERRQ8 ( c , n , s , a1 , a2 , a3 , a4 , a5 , a6 , a7 , a8 )
# Skipping MacroDefinition: SETERRABORT ( comm , n , s )
# Skipping MacroDefinition: CHKERRQ ( n ) ;
# Skipping MacroDefinition: CHKERRABORT ( comm , n ) ;
# Skipping MacroDefinition: CHKERRCONTINUE ( n ) ;

const global PETSCSTACKSIZE = 64

# Skipping MacroDefinition: PetscStackPushNoCheck ( funct , petsc_routine , hot ) do { } while ( 0 )
# Skipping MacroDefinition: PetscStackPopNoCheck do { } while ( 0 )
# Skipping MacroDefinition: PetscFunctionReturn ( a ) return ( a )
# Skipping MacroDefinition: PetscFunctionReturnVoid ( ) return

const global PetscStackPop = CHKMEMQ

# Skipping MacroDefinition: PetscStackPush ( f ) CHKMEMQ
# Skipping MacroDefinition: PetscStackCall ( name , routine ) do { PetscStackPush ( name ) ; routine ; PetscStackPop ; } while ( 0 )
# Skipping MacroDefinition: PetscStackCallStandard ( func , args ) do { PetscStackPush ( # func ) ; ierr = func args ; PetscStackPop ; if ( ierr ) SETERRQ1 ( PETSC_COMM_SELF , PETSC_ERR_LIB , "Error in %s()" , # func ) ; } while ( 0 )

const global PETSC_SMALLEST_CLASSID = 1211211

# Skipping MacroDefinition: PetscObjectComposeFunction ( a , b , d ) PetscObjectComposeFunction_Private ( a , b , ( PetscVoidFunction ) ( d ) )
# Skipping MacroDefinition: PetscOptionsBegin ( comm , prefix , mess , sec ) 0 ; do { PetscOptions PetscOptionsObjectBase ; PetscOptions * PetscOptionsObject = & PetscOptionsObjectBase ; PetscMemzero ( PetscOptionsObject , sizeof ( PetscOptions ) ) ; for ( PetscOptionsObject -> count = ( PetscOptionsPublish ? - 1 : 1 ) ; PetscOptionsObject -> count < 2 ; PetscOptionsObject -> count ++ ) { PetscErrorCode _5_ierr = PetscOptionsBegin_Private ( PetscOptionsObject , comm , prefix , mess , sec ) ; CHKERRQ ( _5_ierr ) ;
# Skipping MacroDefinition: PetscObjectOptionsBegin ( obj ) 0 ; do { PetscOptions PetscOptionsObjectBase ; PetscOptions * PetscOptionsObject = & PetscOptionsObjectBase ; for ( PetscOptionsObject -> count = ( PetscOptionsPublish ? - 1 : 1 ) ; PetscOptionsObject -> count < 2 ; PetscOptionsObject -> count ++ ) { PetscErrorCode _5_ierr = PetscObjectOptionsBegin_Private ( PetscOptionsObject , obj ) ; CHKERRQ ( _5_ierr ) ;
# Skipping MacroDefinition: PetscOptionsEnd ( ) _5_ierr = PetscOptionsEnd_Private ( PetscOptionsObject ) ; CHKERRQ ( _5_ierr ) ; } } while ( 0 )
# Skipping MacroDefinition: PetscOptionsTail ( ) 0 ; { if ( PetscOptionsObject -> count != 1 ) PetscFunctionReturn ( 0 ) ; }
# Skipping MacroDefinition: PetscOptionsEnum ( a , b , c , d , e , f , g ) PetscOptionsEnum_Private ( PetscOptionsObject , a , b , c , d , e , f , g )
# Skipping MacroDefinition: PetscOptionsInt ( a , b , c , d , e , f ) PetscOptionsInt_Private ( PetscOptionsObject , a , b , c , d , e , f )
# Skipping MacroDefinition: PetscOptionsReal ( a , b , c , d , e , f ) PetscOptionsReal_Private ( PetscOptionsObject , a , b , c , d , e , f )
# Skipping MacroDefinition: PetscOptionsScalar ( a , b , c , d , e , f ) PetscOptionsScalar_Private ( PetscOptionsObject , a , b , c , d , e , f )
# Skipping MacroDefinition: PetscOptionsName ( a , b , c , d ) PetscOptionsName_Private ( PetscOptionsObject , a , b , c , d )
# Skipping MacroDefinition: PetscOptionsString ( a , b , c , d , e , f , g ) PetscOptionsString_Private ( PetscOptionsObject , a , b , c , d , e , f , g )
# Skipping MacroDefinition: PetscOptionsBool ( a , b , c , d , e , f ) PetscOptionsBool_Private ( PetscOptionsObject , a , b , c , d , e , f )
# Skipping MacroDefinition: PetscOptionsBoolGroupBegin ( a , b , c , d ) PetscOptionsBoolGroupBegin_Private ( PetscOptionsObject , a , b , c , d )
# Skipping MacroDefinition: PetscOptionsBoolGroup ( a , b , c , d ) PetscOptionsBoolGroup_Private ( PetscOptionsObject , a , b , c , d )
# Skipping MacroDefinition: PetscOptionsBoolGroupEnd ( a , b , c , d ) PetscOptionsBoolGroupEnd_Private ( PetscOptionsObject , a , b , c , d )
# Skipping MacroDefinition: PetscOptionsFList ( a , b , c , d , e , f , g , h ) PetscOptionsFList_Private ( PetscOptionsObject , a , b , c , d , e , f , g , h )
# Skipping MacroDefinition: PetscOptionsEList ( a , b , c , d , e , f , g , h ) PetscOptionsEList_Private ( PetscOptionsObject , a , b , c , d , e , f , g , h )
# Skipping MacroDefinition: PetscOptionsRealArray ( a , b , c , d , e , f ) PetscOptionsRealArray_Private ( PetscOptionsObject , a , b , c , d , e , f )
# Skipping MacroDefinition: PetscOptionsScalarArray ( a , b , c , d , e , f ) PetscOptionsScalarArray_Private ( PetscOptionsObject , a , b , c , d , e , f )
# Skipping MacroDefinition: PetscOptionsIntArray ( a , b , c , d , e , f ) PetscOptionsIntArray_Private ( PetscOptionsObject , a , b , c , d , e , f )
# Skipping MacroDefinition: PetscOptionsStringArray ( a , b , c , d , e , f ) PetscOptionsStringArray_Private ( PetscOptionsObject , a , b , c , d , e , f )
# Skipping MacroDefinition: PetscOptionsBoolArray ( a , b , c , d , e , f ) PetscOptionsBoolArray_Private ( PetscOptionsObject , a , b , c , d , e , f )
# Skipping MacroDefinition: PetscOptionsEnumArray ( a , b , c , d , e , f , g ) PetscOptionsEnumArray_Private ( PetscOptionsObject , a , b , c , d , e , f , g )
# Skipping MacroDefinition: PetscObjectQueryFunction ( obj , name , fptr ) PetscObjectQueryFunction_Private ( ( obj ) , ( name ) , ( PetscVoidFunction * ) ( fptr ) )
# Skipping MacroDefinition: PetscSAWsBlock ( ) 0
# Skipping MacroDefinition: PetscObjectSAWsViewOff ( obj ) 0
# Skipping MacroDefinition: PetscObjectSAWsSetBlock ( obj , flg ) 0
# Skipping MacroDefinition: PetscObjectSAWsBlock ( obj ) 0
# Skipping MacroDefinition: PetscObjectSAWsGrantAccess ( obj ) 0
# Skipping MacroDefinition: PetscObjectSAWsTakeAccess ( obj ) 0
# Skipping MacroDefinition: PetscStackViewSAWs ( ) 0
# Skipping MacroDefinition: PetscStackSAWsViewOff ( ) 0
# Skipping MacroDefinition: PetscFunctionListAdd ( list , name , fptr ) PetscFunctionListAdd_Private ( ( list ) , ( name ) , ( PetscVoidFunction ) ( fptr ) )
# Skipping MacroDefinition: PetscFunctionListFind ( list , name , fptr ) PetscFunctionListFind_Private ( ( list ) , ( name ) , ( PetscVoidFunction * ) ( fptr ) )
# Skipping MacroDefinition: PetscNot ( a ) ( ( a ) ? PETSC_FALSE : PETSC_TRUE )

const global PETSC_EVENT = 1311311

# Skipping MacroDefinition: PetscInfo ( A , S ) 0
# Skipping MacroDefinition: PetscInfo1 ( A , S , a1 ) 0
# Skipping MacroDefinition: PetscInfo2 ( A , S , a1 , a2 ) 0
# Skipping MacroDefinition: PetscInfo3 ( A , S , a1 , a2 , a3 ) 0
# Skipping MacroDefinition: PetscInfo4 ( A , S , a1 , a2 , a3 , a4 ) 0
# Skipping MacroDefinition: PetscInfo5 ( A , S , a1 , a2 , a3 , a4 , a5 ) 0
# Skipping MacroDefinition: PetscInfo6 ( A , S , a1 , a2 , a3 , a4 , a5 , a6 ) 0
# Skipping MacroDefinition: PetscInfo7 ( A , S , a1 , a2 , a3 , a4 , a5 , a6 , a7 ) 0
# Skipping MacroDefinition: PetscLogFlops ( n ) 0
# Skipping MacroDefinition: PetscLogEventActivate ( a ) 0
# Skipping MacroDefinition: PetscLogEventDeactivate ( a ) 0
# Skipping MacroDefinition: PetscLogEventActivateClass ( a ) 0
# Skipping MacroDefinition: PetscLogEventDeactivateClass ( a ) 0
# Skipping MacroDefinition: PetscLogEventSetActiveAll ( a , b ) 0

const global PetscLogPLB = 0
const global PetscLogPLE = 0
const global PetscLogPHC = 0
const global PetscLogPHD = 0

# Skipping MacroDefinition: PetscGetFlops ( a ) ( * ( a ) = 0.0 , 0 )
# Skipping MacroDefinition: PetscLogEventBegin ( e , o1 , o2 , o3 , o4 ) 0
# Skipping MacroDefinition: PetscLogEventEnd ( e , o1 , o2 , o3 , o4 ) 0
# Skipping MacroDefinition: PetscLogEventBarrierBegin ( e , o1 , o2 , o3 , o4 , cm ) 0
# Skipping MacroDefinition: PetscLogEventBarrierEnd ( e , o1 , o2 , o3 , o4 , cm ) 0
# Skipping MacroDefinition: PetscLogObjectParents ( p , n , c ) 0
# Skipping MacroDefinition: PetscLogObjectCreate ( h ) 0
# Skipping MacroDefinition: PetscLogObjectDestroy ( h ) 0
# Skipping MacroDefinition: PetscLogDestroy ( ) 0
# Skipping MacroDefinition: PetscLogStagePush ( a ) 0
# Skipping MacroDefinition: PetscLogStagePop ( ) 0
# Skipping MacroDefinition: PetscLogStageRegister ( a , b ) 0
# Skipping MacroDefinition: PetscLogStagePrint ( a , flg ) 0
# Skipping MacroDefinition: PetscLogView ( viewer ) 0
# Skipping MacroDefinition: PetscLogViewFromOptions ( ) 0
# Skipping MacroDefinition: PetscLogBegin ( ) 0
# Skipping MacroDefinition: PetscLogTraceBegin ( file ) 0
# Skipping MacroDefinition: PetscLogSet ( lb , le ) 0
# Skipping MacroDefinition: PetscLogAllBegin ( ) 0
# Skipping MacroDefinition: PetscLogDump ( c ) 0
# Skipping MacroDefinition: PetscLogEventRegister ( a , b , c ) 0
# Skipping MacroDefinition: PetscLogObjects ( a ) 0
# Skipping MacroDefinition: PetscLogActions ( a ) 0
# Skipping MacroDefinition: MPI_Startall_irecv ( count , number , requests ) MPI_Startall ( number , requests )
# Skipping MacroDefinition: MPI_Startall_isend ( count , number , requests ) MPI_Startall ( number , requests )
# Skipping MacroDefinition: MPI_Start_isend ( count , requests ) MPI_Start ( requests )
# Skipping MacroDefinition: PetscLogStageGetId ( a , b ) ( * ( b ) = 0 , 0 )
# Skipping MacroDefinition: PetscLogStageSetActive ( a , b ) 0
# Skipping MacroDefinition: PetscLogStageGetActive ( a , b ) 0
# Skipping MacroDefinition: PetscLogStageGetVisible ( a , b ) 0
# Skipping MacroDefinition: PetscLogStageSetVisible ( a , b ) 0
# Skipping MacroDefinition: PetscPreLoadBegin ( flag , name ) do { PetscBool PetscPreLoading = flag ; int PetscPreLoadMax , PetscPreLoadIt ; PetscLogStage _stageNum ; PetscErrorCode _3_ierr ; _3_ierr = PetscOptionsGetBool ( NULL , "-preload" , & PetscPreLoading , NULL ) ; CHKERRQ ( _3_ierr ) ; PetscPreLoadMax = ( int ) ( PetscPreLoading ) ; PetscPreLoadingUsed = PetscPreLoading ? PETSC_TRUE : PetscPreLoadingUsed ; for ( PetscPreLoadIt = 0 ; PetscPreLoadIt <= PetscPreLoadMax ; PetscPreLoadIt ++ ) { PetscPreLoadingOn = PetscPreLoading ; _3_ierr = PetscBarrier ( NULL ) ; CHKERRQ ( _3_ierr ) ; if ( PetscPreLoadIt > 0 ) { _3_ierr = PetscLogStageGetId ( name , & _stageNum ) ; CHKERRQ ( _3_ierr ) ; } else { _3_ierr = PetscLogStageRegister ( name , & _stageNum ) ; CHKERRQ ( _3_ierr ) ; } _3_ierr = PetscLogStageSetActive ( _stageNum , ( PetscBool ) ( ! PetscPreLoadMax || PetscPreLoadIt ) ) ; _3_ierr = PetscLogStagePush ( _stageNum ) ; CHKERRQ ( _3_ierr ) ;
# Skipping MacroDefinition: PetscPreLoadEnd ( ) _3_ierr = PetscLogStagePop ( ) ; CHKERRQ ( _3_ierr ) ; PetscPreLoading = PETSC_FALSE ; } \
} while ( 0 )
# Skipping MacroDefinition: PetscPreLoadStage ( name ) do { _3_ierr = PetscLogStagePop ( ) ; CHKERRQ ( _3_ierr ) ; if ( PetscPreLoadIt > 0 ) { _3_ierr = PetscLogStageGetId ( name , & _stageNum ) ; CHKERRQ ( _3_ierr ) ; } else { _3_ierr = PetscLogStageRegister ( name , & _stageNum ) ; CHKERRQ ( _3_ierr ) ; } _3_ierr = PetscLogStageSetActive ( _stageNum , ( PetscBool ) ( ! PetscPreLoadMax || PetscPreLoadIt ) ) ; _3_ierr = PetscLogStagePush ( _stageNum ) ; CHKERRQ ( _3_ierr ) ; } while ( 0 )
# Skipping MacroDefinition: PetscPrefetchBlock ( a , n , rw , t ) do { const char * _p = ( const char * ) ( a ) , * _end = ( const char * ) ( ( a ) + ( n ) ) ; for ( ; _p < _end ; _p += PETSC_LEVEL1_DCACHE_LINESIZE ) PETSC_Prefetch ( _p , ( rw ) , ( t ) ) ; } while ( 0 )

const global PETSC_MPI_INT_MAX = 2147483647
const global PETSC_MPI_INT_MIN = -2147483647
const global PETSC_BLAS_INT_MAX = 2147483647
const global PETSC_BLAS_INT_MIN = -2147483647
const global PETSC_MAX_PATH_LEN = 4096
const global PETSCRAND = "rand"
const global PETSCRAND48 = "rand48"
const global PETSCSPRNG = "sprng"
const global PETSC_BINARY_INT_SIZE = 32 / 8
const global PETSC_BINARY_FLOAT_SIZE = 32 / 8
const global PETSC_BINARY_CHAR_SIZE = 8 / 8
const global PETSC_BINARY_SHORT_SIZE = 16 / 8
const global PETSC_BINARY_DOUBLE_SIZE = 64 / 8

# Skipping MacroDefinition: PETSC_BINARY_SCALAR_SIZE sizeof ( PetscScalar )

const global IS_FILE_CLASSID = 1211218
const global ISGENERAL = "general"
const global ISSTRIDE = "stride"
const global ISBLOCK = "block"
const global PETSCVIEWERSOCKET = "socket"
const global PETSCVIEWERASCII = "ascii"
const global PETSCVIEWERBINARY = "binary"
const global PETSCVIEWERSTRING = "string"
const global PETSCVIEWERDRAW = "draw"
const global PETSCVIEWERVU = "vu"
const global PETSCVIEWERMATHEMATICA = "mathematica"
const global PETSCVIEWERNETCDF = "netcdf"
const global PETSCVIEWERHDF5 = "hdf5"
const global PETSCVIEWERVTK = "vtk"
const global PETSCVIEWERMATLAB = "matlab"
const global PETSCVIEWERSAWS = "saws"
const global PETSC_DRAW_X = "x"
const global PETSC_DRAW_GLUT = "glut"
const global PETSC_DRAW_OPENGLES = "opengles"
const global PETSC_DRAW_NULL = "null"
const global PETSC_DRAW_WIN32 = "win32"
const global PETSC_DRAW_TIKZ = "tikz"

# Skipping MacroDefinition: PetscOptionsViewer ( a , b , c , d , e , f ) PetscOptionsViewer_Private ( PetscOptionsObject , a , b , c , d , e , f ) ;
# Skipping MacroDefinition: PETSC_VIEWER_STDERR_SELF PETSC_VIEWER_STDERR_ ( PETSC_COMM_SELF )
# Skipping MacroDefinition: PETSC_VIEWER_STDERR_WORLD PETSC_VIEWER_STDERR_ ( PETSC_COMM_WORLD )
# Skipping MacroDefinition: PETSC_VIEWER_STDOUT_WORLD PETSC_VIEWER_STDOUT_ ( PETSC_COMM_WORLD )
# Skipping MacroDefinition: PETSC_VIEWER_STDOUT_SELF PETSC_VIEWER_STDOUT_ ( PETSC_COMM_SELF )
# Skipping MacroDefinition: PETSC_VIEWER_DRAW_WORLD PETSC_VIEWER_DRAW_ ( PETSC_COMM_WORLD )
# Skipping MacroDefinition: PETSC_VIEWER_DRAW_SELF PETSC_VIEWER_DRAW_ ( PETSC_COMM_SELF )
# Skipping MacroDefinition: PETSC_VIEWER_SOCKET_WORLD PETSC_VIEWER_SOCKET_ ( PETSC_COMM_WORLD )
# Skipping MacroDefinition: PETSC_VIEWER_SOCKET_SELF PETSC_VIEWER_SOCKET_ ( PETSC_COMM_SELF )
# Skipping MacroDefinition: PETSC_VIEWER_BINARY_WORLD PETSC_VIEWER_BINARY_ ( PETSC_COMM_WORLD )
# Skipping MacroDefinition: PETSC_VIEWER_BINARY_SELF PETSC_VIEWER_BINARY_ ( PETSC_COMM_SELF )
# Skipping MacroDefinition: PETSC_VIEWER_MATLAB_WORLD PETSC_VIEWER_MATLAB_ ( PETSC_COMM_WORLD )
# Skipping MacroDefinition: PETSC_VIEWER_MATLAB_SELF PETSC_VIEWER_MATLAB_ ( PETSC_COMM_SELF )
# Skipping MacroDefinition: PETSC_VIEWER_MATHEMATICA_WORLD ( PetscViewerInitializeMathematicaWorld_Private ( ) , PETSC_VIEWER_MATHEMATICA_WORLD_PRIVATE )

const global VECSEQ = "seq"
const global VECMPI = "mpi"
const global VECSTANDARD = "standard"
const global VECSHARED = "shared"
const global VECSEQCUSP = "seqcusp"
const global VECMPICUSP = "mpicusp"
const global VECCUSP = "cusp"
const global VECSEQVIENNACL = "seqviennacl"
const global VECMPIVIENNACL = "mpiviennacl"
const global VECVIENNACL = "viennacl"
const global VECNEST = "nest"
const global VECSEQPTHREAD = "seqpthread"
const global VECMPIPTHREAD = "mpipthread"
const global VECPTHREAD = "pthread"
const global VEC_FILE_CLASSID = 1211214

# begin enum NormType
typealias NormType Uint32
const global NORM_1 = (UInt32)(0)
const global NORM_2 = (UInt32)(1)
const global NORM_FROBENIUS = (UInt32)(2)
const global NORM_INFINITY = (UInt32)(3)
const global NORM_1_AND_2 = (UInt32)(4)
# end enum NormType

const global NORM_MAX = NORM_INFINITY

# Skipping MacroDefinition: VecLockGet ( x , arg ) * ( arg ) = 0
# Skipping MacroDefinition: VecLockPush ( x ) 0
# Skipping MacroDefinition: VecLockPop ( x ) 0
# Skipping MacroDefinition: VecLocked ( x , arg )

const global MATSAME = "same"
const global MATMAIJ = "maij"
const global MATSEQMAIJ = "seqmaij"
const global MATMPIMAIJ = "mpimaij"
const global MATIS = "is"
const global MATAIJ = "aij"
const global MATSEQAIJ = "seqaij"
const global MATSEQAIJPTHREAD = "seqaijpthread"
const global MATAIJPTHREAD = "aijpthread"
const global MATMPIAIJ = "mpiaij"
const global MATAIJCRL = "aijcrl"
const global MATSEQAIJCRL = "seqaijcrl"
const global MATMPIAIJCRL = "mpiaijcrl"
const global MATAIJCUSP = "aijcusp"
const global MATSEQAIJCUSP = "seqaijcusp"
const global MATMPIAIJCUSP = "mpiaijcusp"
const global MATAIJCUSPARSE = "aijcusparse"
const global MATSEQAIJCUSPARSE = "seqaijcusparse"
const global MATMPIAIJCUSPARSE = "mpiaijcusparse"
const global MATAIJVIENNACL = "aijviennacl"
const global MATSEQAIJVIENNACL = "seqaijviennacl"
const global MATMPIAIJVIENNACL = "mpiaijviennacl"
const global MATAIJPERM = "aijperm"
const global MATSEQAIJPERM = "seqaijperm"
const global MATMPIAIJPERM = "mpiaijperm"
const global MATSHELL = "shell"
const global MATDENSE = "dense"
const global MATSEQDENSE = "seqdense"
const global MATMPIDENSE = "mpidense"
const global MATELEMENTAL = "elemental"
const global MATBAIJ = "baij"
const global MATSEQBAIJ = "seqbaij"
const global MATMPIBAIJ = "mpibaij"
const global MATMPIADJ = "mpiadj"
const global MATSBAIJ = "sbaij"
const global MATSEQSBAIJ = "seqsbaij"
const global MATMPISBAIJ = "mpisbaij"
const global MATSEQBSTRM = "seqbstrm"
const global MATMPIBSTRM = "mpibstrm"
const global MATBSTRM = "bstrm"
const global MATSEQSBSTRM = "seqsbstrm"
const global MATMPISBSTRM = "mpisbstrm"
const global MATSBSTRM = "sbstrm"
const global MATDAAD = "daad"
const global MATMFFD = "mffd"
const global MATNORMAL = "normal"
const global MATLRC = "lrc"
const global MATSCATTER = "scatter"
const global MATBLOCKMAT = "blockmat"
const global MATCOMPOSITE = "composite"
const global MATFFT = "fft"
const global MATFFTW = "fftw"
const global MATSEQCUFFT = "seqcufft"
const global MATTRANSPOSEMAT = "transpose"
const global MATSCHURCOMPLEMENT = "schurcomplement"
const global MATPYTHON = "python"
const global MATHYPRESTRUCT = "hyprestruct"
const global MATHYPRESSTRUCT = "hypresstruct"
const global MATSUBMATRIX = "submatrix"
const global MATLOCALREF = "localref"
const global MATNEST = "nest"

# Skipping MacroDefinition: MatSolverPackage char *

const global MATSOLVERSUPERLU = "superlu"
const global MATSOLVERSUPERLU_DIST = "superlu_dist"
const global MATSOLVERUMFPACK = "umfpack"
const global MATSOLVERCHOLMOD = "cholmod"
const global MATSOLVERESSL = "essl"
const global MATSOLVERLUSOL = "lusol"
const global MATSOLVERMUMPS = "mumps"
const global MATSOLVERMKL_PARDISO = "mkl_pardiso"
const global MATSOLVERMKL_CPARDISO = "mkl_cpardiso"
const global MATSOLVERPASTIX = "pastix"
const global MATSOLVERMATLAB = "matlab"
const global MATSOLVERPETSC = "petsc"
const global MATSOLVERBAS = "bas"
const global MATSOLVERCUSPARSE = "cusparse"
const global MATSOLVERBSTRM = "bstrm"
const global MATSOLVERSBSTRM = "sbstrm"
const global MATSOLVERELEMENTAL = "elemental"
const global MATSOLVERCLIQUE = "clique"
const global MATSOLVERKLU = "klu"
const global MAT_FILE_CLASSID = 1211216

# Skipping MacroDefinition: MatPreallocateInitialize ( comm , nrows , ncols , dnz , onz ) 0 ; \
{ PetscErrorCode _4_ierr ; PetscInt __nrows = ( nrows ) , __ctmp = ( ncols ) , __rstart , __start , __end ; _4_ierr = PetscCalloc2 ( __nrows , & dnz , __nrows , & onz ) ; CHKERRQ ( _4_ierr ) ; __start = 0 ; __end = __start ; _4_ierr = MPI_Scan ( & __ctmp , & __end , 1 , MPIU_INT , MPI_SUM , comm ) ; CHKERRQ ( _4_ierr ) ; __start = __end - __ctmp ; _4_ierr = MPI_Scan ( & __nrows , & __rstart , 1 , MPIU_INT , MPI_SUM , comm ) ; CHKERRQ ( _4_ierr ) ; __rstart = __rstart - __nrows ;
# Skipping MacroDefinition: MatPreallocateSetLocal ( rmap , nrows , rows , cmap , ncols , cols , dnz , onz ) 0 ; \
{ PetscInt __l ; _4_ierr = ISLocalToGlobalMappingApply ( rmap , nrows , rows , rows ) ; CHKERRQ ( _4_ierr ) ; _4_ierr = ISLocalToGlobalMappingApply ( cmap , ncols , cols , cols ) ; CHKERRQ ( _4_ierr ) ; for ( __l = 0 ; __l < nrows ; __l ++ ) { _4_ierr = MatPreallocateSet ( ( rows ) [ __l ] , ncols , cols , dnz , onz ) ; CHKERRQ ( _4_ierr ) ; } \
}
# Skipping MacroDefinition: MatPreallocateSetLocalBlock ( rmap , nrows , rows , cmap , ncols , cols , dnz , onz ) 0 ; \
{ PetscInt __l ; _4_ierr = ISLocalToGlobalMappingApplyBlock ( rmap , nrows , rows , rows ) ; CHKERRQ ( _4_ierr ) ; _4_ierr = ISLocalToGlobalMappingApplyBlock ( cmap , ncols , cols , cols ) ; CHKERRQ ( _4_ierr ) ; for ( __l = 0 ; __l < nrows ; __l ++ ) { _4_ierr = MatPreallocateSet ( ( rows ) [ __l ] , ncols , cols , dnz , onz ) ; CHKERRQ ( _4_ierr ) ; } \
}
# Skipping MacroDefinition: MatPreallocateSymmetricSetLocalBlock ( map , nrows , rows , ncols , cols , dnz , onz ) 0 ; \
{ PetscInt __l ; _4_ierr = ISLocalToGlobalMappingApplyBlock ( map , nrows , rows , rows ) ; CHKERRQ ( _4_ierr ) ; _4_ierr = ISLocalToGlobalMappingApplyBlock ( map , ncols , cols , cols ) ; CHKERRQ ( _4_ierr ) ; for ( __l = 0 ; __l < nrows ; __l ++ ) { _4_ierr = MatPreallocateSymmetricSetBlock ( ( rows ) [ __l ] , ncols , cols , dnz , onz ) ; CHKERRQ ( _4_ierr ) ; } \
}
# Skipping MacroDefinition: MatPreallocateSet ( row , nc , cols , dnz , onz ) 0 ; \
{ PetscInt __i ; if ( row < __rstart ) SETERRQ2 ( PETSC_COMM_SELF , PETSC_ERR_ARG_OUTOFRANGE , "Trying to set preallocation for row %D less than first local row %D" , row , __rstart ) ; if ( row >= __rstart + __nrows ) SETERRQ2 ( PETSC_COMM_SELF , PETSC_ERR_ARG_OUTOFRANGE , "Trying to set preallocation for row %D greater than last local row %D" , row , __rstart + __nrows - 1 ) ; for ( __i = 0 ; __i < nc ; __i ++ ) { if ( ( cols ) [ __i ] < __start || ( cols ) [ __i ] >= __end ) onz [ row - __rstart ] ++ ; else dnz [ row - __rstart ] ++ ; } \
}
# Skipping MacroDefinition: MatPreallocateSymmetricSetBlock ( row , nc , cols , dnz , onz ) 0 ; \
{ PetscInt __i ; for ( __i = 0 ; __i < nc ; __i ++ ) { if ( cols [ __i ] >= __end ) onz [ row - __rstart ] ++ ; else if ( cols [ __i ] >= row ) dnz [ row - __rstart ] ++ ; } \
}
# Skipping MacroDefinition: MatPreallocateLocation ( A , row , ncols , cols , dnz , onz ) 0 ; if ( A ) { ierr = MatSetValues ( A , 1 , & row , ncols , cols , NULL , INSERT_VALUES ) ; CHKERRQ ( ierr ) ; } else { ierr = MatPreallocateSet ( row , ncols , cols , dnz , onz ) ; CHKERRQ ( ierr ) ; }
# Skipping MacroDefinition: MatPreallocateFinalize ( dnz , onz ) 0 ; _4_ierr = PetscFree2 ( dnz , onz ) ; CHKERRQ ( _4_ierr ) ; }

const global MAT_SKIP_ALLOCATION = -4
const global MATORDERINGNATURAL = "natural"
const global MATORDERINGND = "nd"
const global MATORDERING1WD = "1wd"
const global MATORDERINGRCM = "rcm"
const global MATORDERINGQMD = "qmd"
const global MATORDERINGROWLENGTH = "rowlength"
const global MATORDERINGWBM = "wbm"
const global MATORDERINGSPECTRAL = "spectral"
const global MATORDERINGAMD = "amd"
const global MATCOLORINGJP = "jp"
const global MATCOLORINGPOWER = "power"
const global MATCOLORINGNATURAL = "natural"
const global MATCOLORINGSL = "sl"
const global MATCOLORINGLF = "lf"
const global MATCOLORINGID = "id"
const global MATCOLORINGGREEDY = "greedy"
const global MATPARTITIONINGCURRENT = "current"
const global MATPARTITIONINGSQUARE = "square"
const global MATPARTITIONINGPARMETIS = "parmetis"
const global MATPARTITIONINGCHACO = "chaco"
const global MATPARTITIONINGPARTY = "party"
const global MATPARTITIONINGPTSCOTCH = "ptscotch"
const global MP_PARTY_OPT = "opt"
const global MP_PARTY_LIN = "lin"
const global MP_PARTY_SCA = "sca"
const global MP_PARTY_RAN = "ran"
const global MP_PARTY_GBF = "gbf"
const global MP_PARTY_GCF = "gcf"
const global MP_PARTY_BUB = "bub"
const global MP_PARTY_DEF = "def"
const global MP_PARTY_HELPFUL_SETS = "hs"
const global MP_PARTY_KERNIGHAN_LIN = "kl"
const global MP_PARTY_NONE = "no"
const global MATCOARSENMIS = "mis"
const global MATCOARSENHEM = "hem"
const global MATRIX_BINARY_FORMAT_DENSE = -1
const global MATMFFD_DS = "ds"
const global MATMFFD_WP = "wp"
const global PCNONE = "none"
const global PCJACOBI = "jacobi"
const global PCSOR = "sor"
const global PCLU = "lu"
const global PCSHELL = "shell"
const global PCBJACOBI = "bjacobi"
const global PCMG = "mg"
const global PCEISENSTAT = "eisenstat"
const global PCILU = "ilu"
const global PCICC = "icc"
const global PCASM = "asm"
const global PCGASM = "gasm"
const global PCKSP = "ksp"
const global PCCOMPOSITE = "composite"
const global PCREDUNDANT = "redundant"
const global PCSPAI = "spai"
const global PCNN = "nn"
const global PCCHOLESKY = "cholesky"
const global PCPBJACOBI = "pbjacobi"
const global PCMAT = "mat"
const global PCHYPRE = "hypre"
const global PCPARMS = "parms"
const global PCFIELDSPLIT = "fieldsplit"
const global PCTFS = "tfs"
const global PCML = "ml"
const global PCGALERKIN = "galerkin"
const global PCEXOTIC = "exotic"
const global PCCP = "cp"
const global PCBFBT = "bfbt"
const global PCLSC = "lsc"
const global PCPYTHON = "python"
const global PCPFMG = "pfmg"
const global PCSYSPFMG = "syspfmg"
const global PCREDISTRIBUTE = "redistribute"
const global PCSVD = "svd"
const global PCGAMG = "gamg"
const global PCSACUSP = "sacusp"
const global PCSACUSPPOLY = "sacusppoly"
const global PCBICGSTABCUSP = "bicgstabcusp"
const global PCAINVCUSP = "ainvcusp"
const global PCBDDC = "bddc"
const global PCKACZMARZ = "kaczmarz"

# begin enum PCSide
typealias PCSide Cint
const global PC_SIDE_DEFAULT = (Int32)(-1)
const global PC_LEFT = (Int32)(0)
const global PC_RIGHT = (Int32)(1)
const global PC_SYMMETRIC = (Int32)(2)
# end enum PCSide

const global PC_SIDE_MAX = PC_SYMMETRIC + 1
const global PCGAMGAGG = "agg"
const global PCGAMGGEO = "geo"
const global PCGAMGCLASSICAL = "classical"
const global PCGAMGCLASSICALDIRECT = "direct"
const global PCGAMGCLASSICALSTANDARD = "standard"

# begin enum PCMGType
typealias PCMGType Uint32
const global PC_MG_MULTIPLICATIVE = (UInt32)(0)
const global PC_MG_ADDITIVE = (UInt32)(1)
const global PC_MG_FULL = (UInt32)(2)
const global PC_MG_KASKADE = (UInt32)(3)
# end enum PCMGType

const global PC_MG_CASCADE = PC_MG_KASKADE
const global PC_FILE_CLASSID = 1211222
const global KSPRICHARDSON = "richardson"
const global KSPCHEBYSHEV = "chebyshev"
const global KSPCG = "cg"
const global KSPGROPPCG = "groppcg"
const global KSPPIPECG = "pipecg"
const global KSPCGNE = "cgne"
const global KSPNASH = "nash"
const global KSPSTCG = "stcg"
const global KSPGLTR = "gltr"
const global KSPFCG = "fcg"
const global KSPGMRES = "gmres"
const global KSPFGMRES = "fgmres"
const global KSPLGMRES = "lgmres"
const global KSPDGMRES = "dgmres"
const global KSPPGMRES = "pgmres"
const global KSPTCQMR = "tcqmr"
const global KSPBCGS = "bcgs"
const global KSPIBCGS = "ibcgs"
const global KSPFBCGS = "fbcgs"
const global KSPFBCGSR = "fbcgsr"
const global KSPBCGSL = "bcgsl"
const global KSPCGS = "cgs"
const global KSPTFQMR = "tfqmr"
const global KSPCR = "cr"
const global KSPPIPECR = "pipecr"
const global KSPLSQR = "lsqr"
const global KSPPREONLY = "preonly"
const global KSPQCG = "qcg"
const global KSPBICG = "bicg"
const global KSPMINRES = "minres"
const global KSPSYMMLQ = "symmlq"
const global KSPLCD = "lcd"
const global KSPPYTHON = "python"
const global KSPGCR = "gcr"
const global KSP_FILE_CLASSID = 1211223

# begin enum KSPNormType
typealias KSPNormType Cint
const global KSP_NORM_DEFAULT = (Int32)(-1)
const global KSP_NORM_NONE = (Int32)(0)
const global KSP_NORM_PRECONDITIONED = (Int32)(1)
const global KSP_NORM_UNPRECONDITIONED = (Int32)(2)
const global KSP_NORM_NATURAL = (Int32)(3)
# end enum KSPNormType

const global KSP_NORM_MAX = KSP_NORM_NATURAL + 1

# Skipping MacroDefinition: KSPDefaultConverged ( KSPDefaultConverged , KSPConvergedDefault )
# Skipping MacroDefinition: KSPDefaultConvergedDestroy ( KSPDefaultConvergedDestroy , KSPConvergedDefaultDestroy )
# Skipping MacroDefinition: KSPDefaultConvergedCreate ( KSPDefaultConvergedCreate , KSPConvergedDefaultCreate )
# Skipping MacroDefinition: KSPDefaultConvergedSetUIRNorm ( KSPDefaultConvergedSetUIRNorm , KSPConvergedDefaultSetUIRNorm )
# Skipping MacroDefinition: KSPDefaultConvergedSetUMIRNorm ( KSPDefaultConvergedSetUMIRNorm , KSPConvergedDefaultSetUMIRNorm )
# Skipping MacroDefinition: KSPSkipConverged ( KSPSkipConverged , KSPConvergedSkip )

typealias PetscErrorCode Cint
typealias PetscClassId Cint
typealias PetscMPIInt Cint

# begin enum ANONYMOUS_1
typealias ANONYMOUS_1 Uint32
const global ENUM_DUMMY = (UInt32)(0)
# end enum ANONYMOUS_1

# begin enum PetscEnum
typealias PetscEnum Uint32
const global ENUM_DUMMY = (UInt32)(0)
# end enum PetscEnum

typealias Petsc64bitInt Cint
typealias PetscInt Cint
typealias PetscBLASInt Cint

# begin enum ANONYMOUS_2
typealias ANONYMOUS_2 Uint32
const global PETSC_PRECISION_SINGLE = (UInt32)(4)
const global PETSC_PRECISION_DOUBLE = (UInt32)(8)
# end enum ANONYMOUS_2

# begin enum PetscPrecision
typealias PetscPrecision Uint32
const global PETSC_PRECISION_SINGLE = (UInt32)(4)
const global PETSC_PRECISION_DOUBLE = (UInt32)(8)
# end enum PetscPrecision

# begin enum ANONYMOUS_3
typealias ANONYMOUS_3 Uint32
const global PETSC_FALSE = (UInt32)(0)
const global PETSC_TRUE = (UInt32)(1)
# end enum ANONYMOUS_3

# begin enum PetscBool
typealias PetscBool Uint32
const global PETSC_FALSE = (UInt32)(0)
const global PETSC_TRUE = (UInt32)(1)
# end enum PetscBool

typealias MatReal Cint

type petsc_mpiu_2scalar
    a::PetscScalar
    b::PetscScalar
end

typealias PetscLogDouble Cdouble

# begin enum ANONYMOUS_4
typealias ANONYMOUS_4 Uint32
const global PETSC_INT = (UInt32)(0)
const global PETSC_DOUBLE = (UInt32)(1)
const global PETSC_COMPLEX = (UInt32)(2)
const global PETSC_LONG = (UInt32)(3)
const global PETSC_SHORT = (UInt32)(4)
const global PETSC_FLOAT = (UInt32)(5)
const global PETSC_CHAR = (UInt32)(6)
const global PETSC_BIT_LOGICAL = (UInt32)(7)
const global PETSC_ENUM = (UInt32)(8)
const global PETSC_BOOL = (UInt32)(9)
const global PETSC___FLOAT128 = (UInt32)(10)
const global PETSC_OBJECT = (UInt32)(11)
const global PETSC_FUNCTION = (UInt32)(12)
const global PETSC_STRING = (UInt32)(12)
# end enum ANONYMOUS_4

typealias PetscToken Ptr{Void}
typealias PetscObject Ptr{Void}
typealias PetscObjectId Petsc64bitInt
typealias PetscObjectState Petsc64bitInt
typealias PetscFunctionList Ptr{Void}

# begin enum ANONYMOUS_5
typealias ANONYMOUS_5 Uint32
const global FILE_MODE_READ = (UInt32)(0)
const global FILE_MODE_WRITE = (UInt32)(1)
const global FILE_MODE_APPEND = (UInt32)(2)
const global FILE_MODE_UPDATE = (UInt32)(3)
const global FILE_MODE_APPEND_UPDATE = (UInt32)(4)
# end enum ANONYMOUS_5

# begin enum PetscFileMode
typealias PetscFileMode Uint32
const global FILE_MODE_READ = (UInt32)(0)
const global FILE_MODE_WRITE = (UInt32)(1)
const global FILE_MODE_APPEND = (UInt32)(2)
const global FILE_MODE_UPDATE = (UInt32)(3)
const global FILE_MODE_APPEND_UPDATE = (UInt32)(4)
# end enum PetscFileMode

# begin enum ANONYMOUS_6
typealias ANONYMOUS_6 Uint32
const global PETSC_ERROR_INITIAL = (UInt32)(0)
const global PETSC_ERROR_REPEAT = (UInt32)(1)
const global PETSC_ERROR_IN_CXX = (UInt32)(2)
# end enum ANONYMOUS_6

# begin enum PetscErrorType
typealias PetscErrorType Uint32
const global PETSC_ERROR_INITIAL = (UInt32)(0)
const global PETSC_ERROR_REPEAT = (UInt32)(1)
const global PETSC_ERROR_IN_CXX = (UInt32)(2)
# end enum PetscErrorType

# begin enum ANONYMOUS_7
typealias ANONYMOUS_7 Uint32
const global PETSC_FP_TRAP_OFF = (UInt32)(0)
const global PETSC_FP_TRAP_ON = (UInt32)(1)
# end enum ANONYMOUS_7

# begin enum PetscFPTrap
typealias PetscFPTrap Uint32
const global PETSC_FP_TRAP_OFF = (UInt32)(0)
const global PETSC_FP_TRAP_ON = (UInt32)(1)
# end enum PetscFPTrap

immutable Array_64_Ptr
    d1::Ptr{Uint8}
    d2::Ptr{Uint8}
    d3::Ptr{Uint8}
    d4::Ptr{Uint8}
    d5::Ptr{Uint8}
    d6::Ptr{Uint8}
    d7::Ptr{Uint8}
    d8::Ptr{Uint8}
    d9::Ptr{Uint8}
    d10::Ptr{Uint8}
    d11::Ptr{Uint8}
    d12::Ptr{Uint8}
    d13::Ptr{Uint8}
    d14::Ptr{Uint8}
    d15::Ptr{Uint8}
    d16::Ptr{Uint8}
    d17::Ptr{Uint8}
    d18::Ptr{Uint8}
    d19::Ptr{Uint8}
    d20::Ptr{Uint8}
    d21::Ptr{Uint8}
    d22::Ptr{Uint8}
    d23::Ptr{Uint8}
    d24::Ptr{Uint8}
    d25::Ptr{Uint8}
    d26::Ptr{Uint8}
    d27::Ptr{Uint8}
    d28::Ptr{Uint8}
    d29::Ptr{Uint8}
    d30::Ptr{Uint8}
    d31::Ptr{Uint8}
    d32::Ptr{Uint8}
    d33::Ptr{Uint8}
    d34::Ptr{Uint8}
    d35::Ptr{Uint8}
    d36::Ptr{Uint8}
    d37::Ptr{Uint8}
    d38::Ptr{Uint8}
    d39::Ptr{Uint8}
    d40::Ptr{Uint8}
    d41::Ptr{Uint8}
    d42::Ptr{Uint8}
    d43::Ptr{Uint8}
    d44::Ptr{Uint8}
    d45::Ptr{Uint8}
    d46::Ptr{Uint8}
    d47::Ptr{Uint8}
    d48::Ptr{Uint8}
    d49::Ptr{Uint8}
    d50::Ptr{Uint8}
    d51::Ptr{Uint8}
    d52::Ptr{Uint8}
    d53::Ptr{Uint8}
    d54::Ptr{Uint8}
    d55::Ptr{Uint8}
    d56::Ptr{Uint8}
    d57::Ptr{Uint8}
    d58::Ptr{Uint8}
    d59::Ptr{Uint8}
    d60::Ptr{Uint8}
    d61::Ptr{Uint8}
    d62::Ptr{Uint8}
    d63::Ptr{Uint8}
    d64::Ptr{Uint8}
end

zero(::Type{Array_64_Ptr}) = begin  # /home/jared/.julia/v0.4/Clang/src/wrap_c.jl, line 264:
        Array_64_Ptr(fill(zero(Ptr{Uint8}),64)...)
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

zero(::Type{Array_64_Cint}) = begin  # /home/jared/.julia/v0.4/Clang/src/wrap_c.jl, line 264:
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

zero(::Type{Array_64_PetscBool}) = begin  # /home/jared/.julia/v0.4/Clang/src/wrap_c.jl, line 264:
        Array_64_PetscBool(fill(zero(PetscBool),64)...)
    end

type PetscStack
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
typealias PetscViewer Ptr{Void}

# begin enum ANONYMOUS_8
typealias ANONYMOUS_8 Uint32
const global OPTION_INT = (UInt32)(0)
const global OPTION_BOOL = (UInt32)(1)
const global OPTION_REAL = (UInt32)(2)
const global OPTION_FLIST = (UInt32)(3)
const global OPTION_STRING = (UInt32)(4)
const global OPTION_REAL_ARRAY = (UInt32)(5)
const global OPTION_SCALAR_ARRAY = (UInt32)(6)
const global OPTION_HEAD = (UInt32)(7)
const global OPTION_INT_ARRAY = (UInt32)(8)
const global OPTION_ELIST = (UInt32)(9)
const global OPTION_BOOL_ARRAY = (UInt32)(10)
const global OPTION_STRING_ARRAY = (UInt32)(11)
# end enum ANONYMOUS_8

# begin enum PetscOptionType
typealias PetscOptionType Uint32
const global OPTION_INT = (UInt32)(0)
const global OPTION_BOOL = (UInt32)(1)
const global OPTION_REAL = (UInt32)(2)
const global OPTION_FLIST = (UInt32)(3)
const global OPTION_STRING = (UInt32)(4)
const global OPTION_REAL_ARRAY = (UInt32)(5)
const global OPTION_SCALAR_ARRAY = (UInt32)(6)
const global OPTION_HEAD = (UInt32)(7)
const global OPTION_INT_ARRAY = (UInt32)(8)
const global OPTION_ELIST = (UInt32)(9)
const global OPTION_BOOL_ARRAY = (UInt32)(10)
const global OPTION_STRING_ARRAY = (UInt32)(11)
# end enum PetscOptionType

typealias PetscOption Ptr{Void}

type PetscOptions
    count::PetscInt
    next::PetscOption
    prefix::Ptr{Uint8}
    pprefix::Ptr{Uint8}
    title::Ptr{Uint8}
    comm::MPI_Comm
    printhelp::PetscBool
    changedmethod::PetscBool
    alreadyprinted::PetscBool
    object::PetscObject
end

typealias PetscDLHandle Ptr{Void}

# begin enum ANONYMOUS_9
typealias ANONYMOUS_9 Uint32
const global PETSC_DL_DECIDE = (UInt32)(0)
const global PETSC_DL_NOW = (UInt32)(1)
const global PETSC_DL_LOCAL = (UInt32)(2)
# end enum ANONYMOUS_9

# begin enum PetscDLMode
typealias PetscDLMode Uint32
const global PETSC_DL_DECIDE = (UInt32)(0)
const global PETSC_DL_NOW = (UInt32)(1)
const global PETSC_DL_LOCAL = (UInt32)(2)
# end enum PetscDLMode

typealias PetscObjectList Ptr{Void}
typealias PetscDLLibrary Ptr{Void}
typealias PetscLogEvent Cint
typealias PetscLogStage Cint
typealias PetscIntStack Ptr{Void}

type PetscClassRegInfo
    name::Ptr{Uint8}
    classid::PetscClassId
end

type PetscClassPerfInfo
    id::PetscClassId
    creations::Cint
    destructions::Cint
    mem::PetscLogDouble
    descMem::PetscLogDouble
end

typealias PetscClassRegLog Ptr{Void}
typealias PetscClassPerfLog Ptr{Void}

type PetscEventRegInfo
    name::Ptr{Uint8}
    classid::PetscClassId
end

type PetscEventPerfInfo
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

typealias PetscEventRegLog Ptr{Void}
typealias PetscEventPerfLog Ptr{Void}

type PetscStageInfo
    name::Ptr{Uint8}
    used::PetscBool
    perfInfo::PetscEventPerfInfo
    eventLog::PetscEventPerfLog
    classLog::PetscClassPerfLog
end

typealias PetscStageLog Ptr{Void}
typealias PetscContainer Ptr{Void}
typealias PetscRandomType Ptr{Uint8}
typealias PetscRandom Ptr{Void}

# begin enum ANONYMOUS_10
typealias ANONYMOUS_10 Uint32
const global PETSC_BINARY_SEEK_SET = (UInt32)(0)
const global PETSC_BINARY_SEEK_CUR = (UInt32)(1)
const global PETSC_BINARY_SEEK_END = (UInt32)(2)
# end enum ANONYMOUS_10

# begin enum PetscBinarySeekType
typealias PetscBinarySeekType Uint32
const global PETSC_BINARY_SEEK_SET = (UInt32)(0)
const global PETSC_BINARY_SEEK_CUR = (UInt32)(1)
const global PETSC_BINARY_SEEK_END = (UInt32)(2)
# end enum PetscBinarySeekType

# begin enum ANONYMOUS_11
typealias ANONYMOUS_11 Cint
const global PETSC_BUILDTWOSIDED_NOTSET = (Int32)(-1)
const global PETSC_BUILDTWOSIDED_ALLREDUCE = (Int32)(0)
const global PETSC_BUILDTWOSIDED_IBARRIER = (Int32)(1)
# end enum ANONYMOUS_11

# begin enum PetscBuildTwoSidedType
typealias PetscBuildTwoSidedType Cint
const global PETSC_BUILDTWOSIDED_NOTSET = (Int32)(-1)
const global PETSC_BUILDTWOSIDED_ALLREDUCE = (Int32)(0)
const global PETSC_BUILDTWOSIDED_IBARRIER = (Int32)(1)
# end enum PetscBuildTwoSidedType

# begin enum ANONYMOUS_12
typealias ANONYMOUS_12 Uint32
const global NOT_SET_VALUES = (UInt32)(0)
const global INSERT_VALUES = (UInt32)(1)
const global ADD_VALUES = (UInt32)(2)
const global MAX_VALUES = (UInt32)(3)
const global INSERT_ALL_VALUES = (UInt32)(4)
const global ADD_ALL_VALUES = (UInt32)(5)
const global INSERT_BC_VALUES = (UInt32)(6)
const global ADD_BC_VALUES = (UInt32)(7)
# end enum ANONYMOUS_12

# begin enum InsertMode
typealias InsertMode Uint32
const global NOT_SET_VALUES = (UInt32)(0)
const global INSERT_VALUES = (UInt32)(1)
const global ADD_VALUES = (UInt32)(2)
const global MAX_VALUES = (UInt32)(3)
const global INSERT_ALL_VALUES = (UInt32)(4)
const global ADD_ALL_VALUES = (UInt32)(5)
const global INSERT_BC_VALUES = (UInt32)(6)
const global ADD_BC_VALUES = (UInt32)(7)
# end enum InsertMode

# begin enum ANONYMOUS_13
typealias ANONYMOUS_13 Uint32
const global PETSC_SUBCOMM_GENERAL = (UInt32)(0)
const global PETSC_SUBCOMM_CONTIGUOUS = (UInt32)(1)
const global PETSC_SUBCOMM_INTERLACED = (UInt32)(2)
# end enum ANONYMOUS_13

# begin enum PetscSubcommType
typealias PetscSubcommType Uint32
const global PETSC_SUBCOMM_GENERAL = (UInt32)(0)
const global PETSC_SUBCOMM_CONTIGUOUS = (UInt32)(1)
const global PETSC_SUBCOMM_INTERLACED = (UInt32)(2)
# end enum PetscSubcommType

typealias PetscSubcomm Ptr{Void}
typealias PetscSegBuffer Ptr{Void}
typealias PetscSF Ptr{Void}

type PetscSFNode
    rank::PetscInt
    index::PetscInt
end

typealias IS Ptr{Void}
typealias ISLocalToGlobalMapping Ptr{Void}
typealias ISColoring Ptr{Void}
typealias PetscLayout Ptr{Void}
typealias PetscSection Ptr{Void}
typealias ISType Ptr{Uint8}

# begin enum ANONYMOUS_14
typealias ANONYMOUS_14 Uint32
const global IS_GTOLM_MASK = (UInt32)(0)
const global IS_GTOLM_DROP = (UInt32)(1)
# end enum ANONYMOUS_14

# begin enum ISGlobalToLocalMappingType
typealias ISGlobalToLocalMappingType Uint32
const global IS_GTOLM_MASK = (UInt32)(0)
const global IS_GTOLM_DROP = (UInt32)(1)
# end enum ISGlobalToLocalMappingType

# begin enum ANONYMOUS_15
typealias ANONYMOUS_15 Uint32
const global IS_COLORING_GLOBAL = (UInt32)(0)
const global IS_COLORING_GHOSTED = (UInt32)(1)
# end enum ANONYMOUS_15

# begin enum ISColoringType
typealias ISColoringType Uint32
const global IS_COLORING_GLOBAL = (UInt32)(0)
const global IS_COLORING_GHOSTED = (UInt32)(1)
# end enum ISColoringType

typealias PETSC_IS_COLOR_VALUE_TYPE Uint32
typealias PetscViewerType Ptr{Uint8}
typealias PetscDrawType Ptr{Uint8}
typealias PetscDraw Ptr{Void}
typealias PetscDrawAxis Ptr{Void}
typealias PetscDrawLG Ptr{Void}
typealias PetscDrawSP Ptr{Void}
typealias PetscDrawHG Ptr{Void}
typealias PetscDrawBar Ptr{Void}

# begin enum ANONYMOUS_16
typealias ANONYMOUS_16 Uint32
const global PETSC_VIEWER_DEFAULT = (UInt32)(0)
const global PETSC_VIEWER_ASCII_MATLAB = (UInt32)(1)
const global PETSC_VIEWER_ASCII_MATHEMATICA = (UInt32)(2)
const global PETSC_VIEWER_ASCII_IMPL = (UInt32)(3)
const global PETSC_VIEWER_ASCII_INFO = (UInt32)(4)
const global PETSC_VIEWER_ASCII_INFO_DETAIL = (UInt32)(5)
const global PETSC_VIEWER_ASCII_COMMON = (UInt32)(6)
const global PETSC_VIEWER_ASCII_SYMMODU = (UInt32)(7)
const global PETSC_VIEWER_ASCII_INDEX = (UInt32)(8)
const global PETSC_VIEWER_ASCII_DENSE = (UInt32)(9)
const global PETSC_VIEWER_ASCII_MATRIXMARKET = (UInt32)(10)
const global PETSC_VIEWER_ASCII_VTK = (UInt32)(11)
const global PETSC_VIEWER_ASCII_VTK_CELL = (UInt32)(12)
const global PETSC_VIEWER_ASCII_VTK_COORDS = (UInt32)(13)
const global PETSC_VIEWER_ASCII_PCICE = (UInt32)(14)
const global PETSC_VIEWER_ASCII_PYTHON = (UInt32)(15)
const global PETSC_VIEWER_ASCII_FACTOR_INFO = (UInt32)(16)
const global PETSC_VIEWER_ASCII_LATEX = (UInt32)(17)
const global PETSC_VIEWER_DRAW_BASIC = (UInt32)(18)
const global PETSC_VIEWER_DRAW_LG = (UInt32)(19)
const global PETSC_VIEWER_DRAW_CONTOUR = (UInt32)(20)
const global PETSC_VIEWER_DRAW_PORTS = (UInt32)(21)
const global PETSC_VIEWER_VTK_VTS = (UInt32)(22)
const global PETSC_VIEWER_VTK_VTR = (UInt32)(23)
const global PETSC_VIEWER_VTK_VTU = (UInt32)(24)
const global PETSC_VIEWER_BINARY_MATLAB = (UInt32)(25)
const global PETSC_VIEWER_NATIVE = (UInt32)(26)
const global PETSC_VIEWER_HDF5_VIZ = (UInt32)(27)
const global PETSC_VIEWER_NOFORMAT = (UInt32)(28)
# end enum ANONYMOUS_16

# begin enum PetscViewerFormat
typealias PetscViewerFormat Uint32
const global PETSC_VIEWER_DEFAULT = (UInt32)(0)
const global PETSC_VIEWER_ASCII_MATLAB = (UInt32)(1)
const global PETSC_VIEWER_ASCII_MATHEMATICA = (UInt32)(2)
const global PETSC_VIEWER_ASCII_IMPL = (UInt32)(3)
const global PETSC_VIEWER_ASCII_INFO = (UInt32)(4)
const global PETSC_VIEWER_ASCII_INFO_DETAIL = (UInt32)(5)
const global PETSC_VIEWER_ASCII_COMMON = (UInt32)(6)
const global PETSC_VIEWER_ASCII_SYMMODU = (UInt32)(7)
const global PETSC_VIEWER_ASCII_INDEX = (UInt32)(8)
const global PETSC_VIEWER_ASCII_DENSE = (UInt32)(9)
const global PETSC_VIEWER_ASCII_MATRIXMARKET = (UInt32)(10)
const global PETSC_VIEWER_ASCII_VTK = (UInt32)(11)
const global PETSC_VIEWER_ASCII_VTK_CELL = (UInt32)(12)
const global PETSC_VIEWER_ASCII_VTK_COORDS = (UInt32)(13)
const global PETSC_VIEWER_ASCII_PCICE = (UInt32)(14)
const global PETSC_VIEWER_ASCII_PYTHON = (UInt32)(15)
const global PETSC_VIEWER_ASCII_FACTOR_INFO = (UInt32)(16)
const global PETSC_VIEWER_ASCII_LATEX = (UInt32)(17)
const global PETSC_VIEWER_DRAW_BASIC = (UInt32)(18)
const global PETSC_VIEWER_DRAW_LG = (UInt32)(19)
const global PETSC_VIEWER_DRAW_CONTOUR = (UInt32)(20)
const global PETSC_VIEWER_DRAW_PORTS = (UInt32)(21)
const global PETSC_VIEWER_VTK_VTS = (UInt32)(22)
const global PETSC_VIEWER_VTK_VTR = (UInt32)(23)
const global PETSC_VIEWER_VTK_VTU = (UInt32)(24)
const global PETSC_VIEWER_BINARY_MATLAB = (UInt32)(25)
const global PETSC_VIEWER_NATIVE = (UInt32)(26)
const global PETSC_VIEWER_HDF5_VIZ = (UInt32)(27)
const global PETSC_VIEWER_NOFORMAT = (UInt32)(28)
# end enum PetscViewerFormat

# begin enum ANONYMOUS_17
typealias ANONYMOUS_17 Uint32
const global PETSC_VTK_POINT_FIELD = (UInt32)(0)
const global PETSC_VTK_POINT_VECTOR_FIELD = (UInt32)(1)
const global PETSC_VTK_CELL_FIELD = (UInt32)(2)
const global PETSC_VTK_CELL_VECTOR_FIELD = (UInt32)(3)
# end enum ANONYMOUS_17

# begin enum PetscViewerVTKFieldType
typealias PetscViewerVTKFieldType Uint32
const global PETSC_VTK_POINT_FIELD = (UInt32)(0)
const global PETSC_VTK_POINT_VECTOR_FIELD = (UInt32)(1)
const global PETSC_VTK_CELL_FIELD = (UInt32)(2)
const global PETSC_VTK_CELL_VECTOR_FIELD = (UInt32)(3)
# end enum PetscViewerVTKFieldType

typealias PetscViewers Ptr{Void}
typealias Vec Ptr{Void}
typealias VecScatter Ptr{Void}

# begin enum ANONYMOUS_18
typealias ANONYMOUS_18 Uint32
const global SCATTER_FORWARD = (UInt32)(0)
const global SCATTER_REVERSE = (UInt32)(1)
const global SCATTER_FORWARD_LOCAL = (UInt32)(2)
const global SCATTER_REVERSE_LOCAL = (UInt32)(3)
const global SCATTER_LOCAL = (UInt32)(2)
# end enum ANONYMOUS_18

# begin enum ScatterMode
typealias ScatterMode Uint32
const global SCATTER_FORWARD = (UInt32)(0)
const global SCATTER_REVERSE = (UInt32)(1)
const global SCATTER_FORWARD_LOCAL = (UInt32)(2)
const global SCATTER_REVERSE_LOCAL = (UInt32)(3)
const global SCATTER_LOCAL = (UInt32)(2)
# end enum ScatterMode

typealias VecType Ptr{Uint8}

# begin enum ANONYMOUS_19
typealias ANONYMOUS_19 Uint32
const global NORM_1 = (UInt32)(0)
const global NORM_2 = (UInt32)(1)
const global NORM_FROBENIUS = (UInt32)(2)
const global NORM_INFINITY = (UInt32)(3)
const global NORM_1_AND_2 = (UInt32)(4)
# end enum ANONYMOUS_19

# begin enum ANONYMOUS_20
typealias ANONYMOUS_20 Uint32
const global VEC_IGNORE_OFF_PROC_ENTRIES = (UInt32)(0)
const global VEC_IGNORE_NEGATIVE_INDICES = (UInt32)(1)
# end enum ANONYMOUS_20

# begin enum VecOption
typealias VecOption Uint32
const global VEC_IGNORE_OFF_PROC_ENTRIES = (UInt32)(0)
const global VEC_IGNORE_NEGATIVE_INDICES = (UInt32)(1)
# end enum VecOption

# begin enum ANONYMOUS_21
typealias ANONYMOUS_21 Uint32
const global VECOP_VIEW = (UInt32)(33)
const global VECOP_LOAD = (UInt32)(41)
const global VECOP_DUPLICATE = (UInt32)(0)
# end enum ANONYMOUS_21

# begin enum VecOperation
typealias VecOperation Uint32
const global VECOP_VIEW = (UInt32)(33)
const global VECOP_LOAD = (UInt32)(41)
const global VECOP_DUPLICATE = (UInt32)(0)
# end enum VecOperation

typealias Vecs Ptr{Void}
typealias Mat Ptr{Void}
typealias MatType Ptr{Uint8}

# begin enum ANONYMOUS_22
typealias ANONYMOUS_22 Uint32
const global MAT_FACTOR_NONE = (UInt32)(0)
const global MAT_FACTOR_LU = (UInt32)(1)
const global MAT_FACTOR_CHOLESKY = (UInt32)(2)
const global MAT_FACTOR_ILU = (UInt32)(3)
const global MAT_FACTOR_ICC = (UInt32)(4)
const global MAT_FACTOR_ILUDT = (UInt32)(5)
# end enum ANONYMOUS_22

# begin enum MatFactorType
typealias MatFactorType Uint32
const global MAT_FACTOR_NONE = (UInt32)(0)
const global MAT_FACTOR_LU = (UInt32)(1)
const global MAT_FACTOR_CHOLESKY = (UInt32)(2)
const global MAT_FACTOR_ILU = (UInt32)(3)
const global MAT_FACTOR_ICC = (UInt32)(4)
const global MAT_FACTOR_ILUDT = (UInt32)(5)
# end enum MatFactorType

# begin enum ANONYMOUS_23
typealias ANONYMOUS_23 Uint32
const global MAT_INITIAL_MATRIX = (UInt32)(0)
const global MAT_REUSE_MATRIX = (UInt32)(1)
const global MAT_IGNORE_MATRIX = (UInt32)(2)
# end enum ANONYMOUS_23

# begin enum MatReuse
typealias MatReuse Uint32
const global MAT_INITIAL_MATRIX = (UInt32)(0)
const global MAT_REUSE_MATRIX = (UInt32)(1)
const global MAT_IGNORE_MATRIX = (UInt32)(2)
# end enum MatReuse

# begin enum ANONYMOUS_24
typealias ANONYMOUS_24 Uint32
const global MAT_DO_NOT_GET_VALUES = (UInt32)(0)
const global MAT_GET_VALUES = (UInt32)(1)
# end enum ANONYMOUS_24

# begin enum MatGetSubMatrixOption
typealias MatGetSubMatrixOption Uint32
const global MAT_DO_NOT_GET_VALUES = (UInt32)(0)
const global MAT_GET_VALUES = (UInt32)(1)
# end enum MatGetSubMatrixOption

# begin enum ANONYMOUS_25
typealias ANONYMOUS_25 Uint32
const global DIFFERENT_NONZERO_PATTERN = (UInt32)(0)
const global SUBSET_NONZERO_PATTERN = (UInt32)(1)
const global SAME_NONZERO_PATTERN = (UInt32)(2)
# end enum ANONYMOUS_25

# begin enum MatStructure
typealias MatStructure Uint32
const global DIFFERENT_NONZERO_PATTERN = (UInt32)(0)
const global SUBSET_NONZERO_PATTERN = (UInt32)(1)
const global SAME_NONZERO_PATTERN = (UInt32)(2)
# end enum MatStructure

# begin enum ANONYMOUS_26
typealias ANONYMOUS_26 Uint32
const global MAT_COMPOSITE_ADDITIVE = (UInt32)(0)
const global MAT_COMPOSITE_MULTIPLICATIVE = (UInt32)(1)
# end enum ANONYMOUS_26

# begin enum MatCompositeType
typealias MatCompositeType Uint32
const global MAT_COMPOSITE_ADDITIVE = (UInt32)(0)
const global MAT_COMPOSITE_MULTIPLICATIVE = (UInt32)(1)
# end enum MatCompositeType

type MatStencil
    k::PetscInt
    j::PetscInt
    i::PetscInt
    c::PetscInt
end

# begin enum ANONYMOUS_27
typealias ANONYMOUS_27 Uint32
const global MAT_FLUSH_ASSEMBLY = (UInt32)(1)
const global MAT_FINAL_ASSEMBLY = (UInt32)(0)
# end enum ANONYMOUS_27

# begin enum MatAssemblyType
typealias MatAssemblyType Uint32
const global MAT_FLUSH_ASSEMBLY = (UInt32)(1)
const global MAT_FINAL_ASSEMBLY = (UInt32)(0)
# end enum MatAssemblyType

# begin enum ANONYMOUS_28
typealias ANONYMOUS_28 Cint
const global MAT_OPTION_MIN = (Int32)(-5)
const global MAT_NEW_NONZERO_LOCATION_ERR = (Int32)(-4)
const global MAT_UNUSED_NONZERO_LOCATION_ERR = (Int32)(-3)
const global MAT_NEW_NONZERO_ALLOCATION_ERR = (Int32)(-2)
const global MAT_ROW_ORIENTED = (Int32)(-1)
const global MAT_SYMMETRIC = (Int32)(1)
const global MAT_STRUCTURALLY_SYMMETRIC = (Int32)(2)
const global MAT_NEW_DIAGONALS = (Int32)(3)
const global MAT_IGNORE_OFF_PROC_ENTRIES = (Int32)(4)
const global MAT_USE_HASH_TABLE = (Int32)(5)
const global MAT_KEEP_NONZERO_PATTERN = (Int32)(6)
const global MAT_IGNORE_ZERO_ENTRIES = (Int32)(7)
const global MAT_USE_INODES = (Int32)(8)
const global MAT_HERMITIAN = (Int32)(9)
const global MAT_SYMMETRY_ETERNAL = (Int32)(10)
const global MAT_DUMMY = (Int32)(11)
const global MAT_IGNORE_LOWER_TRIANGULAR = (Int32)(12)
const global MAT_ERROR_LOWER_TRIANGULAR = (Int32)(13)
const global MAT_GETROW_UPPERTRIANGULAR = (Int32)(14)
const global MAT_SPD = (Int32)(15)
const global MAT_NO_OFF_PROC_ZERO_ROWS = (Int32)(16)
const global MAT_NO_OFF_PROC_ENTRIES = (Int32)(17)
const global MAT_NEW_NONZERO_LOCATIONS = (Int32)(18)
const global MAT_OPTION_MAX = (Int32)(19)
# end enum ANONYMOUS_28

# begin enum MatOption
typealias MatOption Cint
const global MAT_OPTION_MIN = (Int32)(-5)
const global MAT_NEW_NONZERO_LOCATION_ERR = (Int32)(-4)
const global MAT_UNUSED_NONZERO_LOCATION_ERR = (Int32)(-3)
const global MAT_NEW_NONZERO_ALLOCATION_ERR = (Int32)(-2)
const global MAT_ROW_ORIENTED = (Int32)(-1)
const global MAT_SYMMETRIC = (Int32)(1)
const global MAT_STRUCTURALLY_SYMMETRIC = (Int32)(2)
const global MAT_NEW_DIAGONALS = (Int32)(3)
const global MAT_IGNORE_OFF_PROC_ENTRIES = (Int32)(4)
const global MAT_USE_HASH_TABLE = (Int32)(5)
const global MAT_KEEP_NONZERO_PATTERN = (Int32)(6)
const global MAT_IGNORE_ZERO_ENTRIES = (Int32)(7)
const global MAT_USE_INODES = (Int32)(8)
const global MAT_HERMITIAN = (Int32)(9)
const global MAT_SYMMETRY_ETERNAL = (Int32)(10)
const global MAT_DUMMY = (Int32)(11)
const global MAT_IGNORE_LOWER_TRIANGULAR = (Int32)(12)
const global MAT_ERROR_LOWER_TRIANGULAR = (Int32)(13)
const global MAT_GETROW_UPPERTRIANGULAR = (Int32)(14)
const global MAT_SPD = (Int32)(15)
const global MAT_NO_OFF_PROC_ZERO_ROWS = (Int32)(16)
const global MAT_NO_OFF_PROC_ENTRIES = (Int32)(17)
const global MAT_NEW_NONZERO_LOCATIONS = (Int32)(18)
const global MAT_OPTION_MAX = (Int32)(19)
# end enum MatOption

# begin enum ANONYMOUS_29
typealias ANONYMOUS_29 Uint32
const global MAT_DO_NOT_COPY_VALUES = (UInt32)(0)
const global MAT_COPY_VALUES = (UInt32)(1)
const global MAT_SHARE_NONZERO_PATTERN = (UInt32)(2)
# end enum ANONYMOUS_29

# begin enum MatDuplicateOption
typealias MatDuplicateOption Uint32
const global MAT_DO_NOT_COPY_VALUES = (UInt32)(0)
const global MAT_COPY_VALUES = (UInt32)(1)
const global MAT_SHARE_NONZERO_PATTERN = (UInt32)(2)
# end enum MatDuplicateOption

type MatInfo
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

# begin enum ANONYMOUS_30
typealias ANONYMOUS_30 Uint32
const global MAT_LOCAL = (UInt32)(1)
const global MAT_GLOBAL_MAX = (UInt32)(2)
const global MAT_GLOBAL_SUM = (UInt32)(3)
# end enum ANONYMOUS_30

# begin enum MatInfoType
typealias MatInfoType Uint32
const global MAT_LOCAL = (UInt32)(1)
const global MAT_GLOBAL_MAX = (UInt32)(2)
const global MAT_GLOBAL_SUM = (UInt32)(3)
# end enum MatInfoType

typealias MatOrderingType Ptr{Uint8}

# begin enum ANONYMOUS_31
typealias ANONYMOUS_31 Uint32
const global MAT_SHIFT_NONE = (UInt32)(0)
const global MAT_SHIFT_NONZERO = (UInt32)(1)
const global MAT_SHIFT_POSITIVE_DEFINITE = (UInt32)(2)
const global MAT_SHIFT_INBLOCKS = (UInt32)(3)
# end enum ANONYMOUS_31

# begin enum MatFactorShiftType
typealias MatFactorShiftType Uint32
const global MAT_SHIFT_NONE = (UInt32)(0)
const global MAT_SHIFT_NONZERO = (UInt32)(1)
const global MAT_SHIFT_POSITIVE_DEFINITE = (UInt32)(2)
const global MAT_SHIFT_INBLOCKS = (UInt32)(3)
# end enum MatFactorShiftType

type MatFactorInfo
    diagonal_fill::Cint
    usedt::Cint
    dt::Cint
    dtcol::Cint
    dtcount::Cint
    fill::Cint
    levels::Cint
    pivotinblocks::Cint
    zeropivot::Cint
    shifttype::Cint
    shiftamount::Cint
end

# begin enum ANONYMOUS_32
typealias ANONYMOUS_32 Uint32
const global SOR_FORWARD_SWEEP = (UInt32)(1)
const global SOR_BACKWARD_SWEEP = (UInt32)(2)
const global SOR_SYMMETRIC_SWEEP = (UInt32)(3)
const global SOR_LOCAL_FORWARD_SWEEP = (UInt32)(4)
const global SOR_LOCAL_BACKWARD_SWEEP = (UInt32)(8)
const global SOR_LOCAL_SYMMETRIC_SWEEP = (UInt32)(12)
const global SOR_ZERO_INITIAL_GUESS = (UInt32)(16)
const global SOR_EISENSTAT = (UInt32)(32)
const global SOR_APPLY_UPPER = (UInt32)(64)
const global SOR_APPLY_LOWER = (UInt32)(128)
# end enum ANONYMOUS_32

# begin enum MatSORType
typealias MatSORType Uint32
const global SOR_FORWARD_SWEEP = (UInt32)(1)
const global SOR_BACKWARD_SWEEP = (UInt32)(2)
const global SOR_SYMMETRIC_SWEEP = (UInt32)(3)
const global SOR_LOCAL_FORWARD_SWEEP = (UInt32)(4)
const global SOR_LOCAL_BACKWARD_SWEEP = (UInt32)(8)
const global SOR_LOCAL_SYMMETRIC_SWEEP = (UInt32)(12)
const global SOR_ZERO_INITIAL_GUESS = (UInt32)(16)
const global SOR_EISENSTAT = (UInt32)(32)
const global SOR_APPLY_UPPER = (UInt32)(64)
const global SOR_APPLY_LOWER = (UInt32)(128)
# end enum MatSORType

typealias MatColoring Ptr{Void}
typealias MatColoringType Ptr{Uint8}

# begin enum ANONYMOUS_33
typealias ANONYMOUS_33 Uint32
const global MAT_COLORING_WEIGHT_RANDOM = (UInt32)(0)
const global MAT_COLORING_WEIGHT_LEXICAL = (UInt32)(1)
const global MAT_COLORING_WEIGHT_LF = (UInt32)(2)
const global MAT_COLORING_WEIGHT_SL = (UInt32)(3)
# end enum ANONYMOUS_33

# begin enum MatColoringWeightType
typealias MatColoringWeightType Uint32
const global MAT_COLORING_WEIGHT_RANDOM = (UInt32)(0)
const global MAT_COLORING_WEIGHT_LEXICAL = (UInt32)(1)
const global MAT_COLORING_WEIGHT_LF = (UInt32)(2)
const global MAT_COLORING_WEIGHT_SL = (UInt32)(3)
# end enum MatColoringWeightType

typealias MatFDColoring Ptr{Void}
typealias MatTransposeColoring Ptr{Void}
typealias MatPartitioning Ptr{Void}
typealias MatPartitioningType Ptr{Uint8}

# begin enum ANONYMOUS_34
typealias ANONYMOUS_34 Uint32
const global MP_CHACO_MULTILEVEL = (UInt32)(1)
const global MP_CHACO_SPECTRAL = (UInt32)(2)
const global MP_CHACO_LINEAR = (UInt32)(4)
const global MP_CHACO_RANDOM = (UInt32)(5)
const global MP_CHACO_SCATTERED = (UInt32)(6)
# end enum ANONYMOUS_34

# begin enum MPChacoGlobalType
typealias MPChacoGlobalType Uint32
const global MP_CHACO_MULTILEVEL = (UInt32)(1)
const global MP_CHACO_SPECTRAL = (UInt32)(2)
const global MP_CHACO_LINEAR = (UInt32)(4)
const global MP_CHACO_RANDOM = (UInt32)(5)
const global MP_CHACO_SCATTERED = (UInt32)(6)
# end enum MPChacoGlobalType

# begin enum ANONYMOUS_35
typealias ANONYMOUS_35 Uint32
const global MP_CHACO_KERNIGHAN = (UInt32)(1)
const global MP_CHACO_NONE = (UInt32)(2)
# end enum ANONYMOUS_35

# begin enum MPChacoLocalType
typealias MPChacoLocalType Uint32
const global MP_CHACO_KERNIGHAN = (UInt32)(1)
const global MP_CHACO_NONE = (UInt32)(2)
# end enum MPChacoLocalType

# begin enum ANONYMOUS_36
typealias ANONYMOUS_36 Uint32
const global MP_CHACO_LANCZOS = (UInt32)(0)
const global MP_CHACO_RQI = (UInt32)(1)
# end enum ANONYMOUS_36

# begin enum MPChacoEigenType
typealias MPChacoEigenType Uint32
const global MP_CHACO_LANCZOS = (UInt32)(0)
const global MP_CHACO_RQI = (UInt32)(1)
# end enum MPChacoEigenType

# begin enum ANONYMOUS_37
typealias ANONYMOUS_37 Uint32
const global MP_PTSCOTCH_QUALITY = (UInt32)(0)
const global MP_PTSCOTCH_SPEED = (UInt32)(1)
const global MP_PTSCOTCH_BALANCE = (UInt32)(2)
const global MP_PTSCOTCH_SAFETY = (UInt32)(3)
const global MP_PTSCOTCH_SCALABILITY = (UInt32)(4)
# end enum ANONYMOUS_37

# begin enum MPPTScotchStrategyType
typealias MPPTScotchStrategyType Uint32
const global MP_PTSCOTCH_QUALITY = (UInt32)(0)
const global MP_PTSCOTCH_SPEED = (UInt32)(1)
const global MP_PTSCOTCH_BALANCE = (UInt32)(2)
const global MP_PTSCOTCH_SAFETY = (UInt32)(3)
const global MP_PTSCOTCH_SCALABILITY = (UInt32)(4)
# end enum MPPTScotchStrategyType

typealias MatCoarsen Ptr{Void}
typealias MatCoarsenType Ptr{Uint8}

type PetscCDIntNd
    next::Ptr{_PetscCDIntNd}
    gid::PetscInt
end

type PetscCDArrNd
    next::Ptr{_PetscCDArrNd}
    array::Ptr{_PetscCDIntNd}
end

type PetscCoarsenData
    pool_list::PetscCDArrNd
    new_node::Ptr{PetscCDIntNd}
    new_left::PetscInt
    chk_sz::PetscInt
    extra_nodes::Ptr{PetscCDIntNd}
    array::Ptr{Ptr{PetscCDIntNd}}
    size::PetscInt
    mat::Mat
end

# begin enum ANONYMOUS_38
typealias ANONYMOUS_38 Uint32
const global MATOP_SET_VALUES = (UInt32)(0)
const global MATOP_GET_ROW = (UInt32)(1)
const global MATOP_RESTORE_ROW = (UInt32)(2)
const global MATOP_MULT = (UInt32)(3)
const global MATOP_MULT_ADD = (UInt32)(4)
const global MATOP_MULT_TRANSPOSE = (UInt32)(5)
const global MATOP_MULT_TRANSPOSE_ADD = (UInt32)(6)
const global MATOP_SOLVE = (UInt32)(7)
const global MATOP_SOLVE_ADD = (UInt32)(8)
const global MATOP_SOLVE_TRANSPOSE = (UInt32)(9)
const global MATOP_SOLVE_TRANSPOSE_ADD = (UInt32)(10)
const global MATOP_LUFACTOR = (UInt32)(11)
const global MATOP_CHOLESKYFACTOR = (UInt32)(12)
const global MATOP_SOR = (UInt32)(13)
const global MATOP_TRANSPOSE = (UInt32)(14)
const global MATOP_GETINFO = (UInt32)(15)
const global MATOP_EQUAL = (UInt32)(16)
const global MATOP_GET_DIAGONAL = (UInt32)(17)
const global MATOP_DIAGONAL_SCALE = (UInt32)(18)
const global MATOP_NORM = (UInt32)(19)
const global MATOP_ASSEMBLY_BEGIN = (UInt32)(20)
const global MATOP_ASSEMBLY_END = (UInt32)(21)
const global MATOP_SET_OPTION = (UInt32)(22)
const global MATOP_ZERO_ENTRIES = (UInt32)(23)
const global MATOP_ZERO_ROWS = (UInt32)(24)
const global MATOP_LUFACTOR_SYMBOLIC = (UInt32)(25)
const global MATOP_LUFACTOR_NUMERIC = (UInt32)(26)
const global MATOP_CHOLESKY_FACTOR_SYMBOLIC = (UInt32)(27)
const global MATOP_CHOLESKY_FACTOR_NUMERIC = (UInt32)(28)
const global MATOP_SETUP_PREALLOCATION = (UInt32)(29)
const global MATOP_ILUFACTOR_SYMBOLIC = (UInt32)(30)
const global MATOP_ICCFACTOR_SYMBOLIC = (UInt32)(31)
const global MATOP_DUPLICATE = (UInt32)(34)
const global MATOP_FORWARD_SOLVE = (UInt32)(35)
const global MATOP_BACKWARD_SOLVE = (UInt32)(36)
const global MATOP_ILUFACTOR = (UInt32)(37)
const global MATOP_ICCFACTOR = (UInt32)(38)
const global MATOP_AXPY = (UInt32)(39)
const global MATOP_GET_SUBMATRICES = (UInt32)(40)
const global MATOP_INCREASE_OVERLAP = (UInt32)(41)
const global MATOP_GET_VALUES = (UInt32)(42)
const global MATOP_COPY = (UInt32)(43)
const global MATOP_GET_ROW_MAX = (UInt32)(44)
const global MATOP_SCALE = (UInt32)(45)
const global MATOP_SHIFT = (UInt32)(46)
const global MATOP_DIAGONAL_SET = (UInt32)(47)
const global MATOP_ZERO_ROWS_COLUMNS = (UInt32)(48)
const global MATOP_SET_RANDOM = (UInt32)(49)
const global MATOP_GET_ROW_IJ = (UInt32)(50)
const global MATOP_RESTORE_ROW_IJ = (UInt32)(51)
const global MATOP_GET_COLUMN_IJ = (UInt32)(52)
const global MATOP_RESTORE_COLUMN_IJ = (UInt32)(53)
const global MATOP_FDCOLORING_CREATE = (UInt32)(54)
const global MATOP_COLORING_PATCH = (UInt32)(55)
const global MATOP_SET_UNFACTORED = (UInt32)(56)
const global MATOP_PERMUTE = (UInt32)(57)
const global MATOP_SET_VALUES_BLOCKED = (UInt32)(58)
const global MATOP_GET_SUBMATRIX = (UInt32)(59)
const global MATOP_DESTROY = (UInt32)(60)
const global MATOP_VIEW = (UInt32)(61)
const global MATOP_CONVERT_FROM = (UInt32)(62)
const global MATOP_MATMAT_MULT = (UInt32)(63)
const global MATOP_MATMAT_MULT_SYMBOLIC = (UInt32)(64)
const global MATOP_MATMAT_MULT_NUMERIC = (UInt32)(65)
const global MATOP_SET_LOCAL_TO_GLOBAL_MAP = (UInt32)(66)
const global MATOP_SET_VALUES_LOCAL = (UInt32)(67)
const global MATOP_ZERO_ROWS_LOCAL = (UInt32)(68)
const global MATOP_GET_ROW_MAX_ABS = (UInt32)(69)
const global MATOP_GET_ROW_MIN_ABS = (UInt32)(70)
const global MATOP_CONVERT = (UInt32)(71)
const global MATOP_SET_COLORING = (UInt32)(72)
const global MATOP_SET_VALUES_ADIFOR = (UInt32)(74)
const global MATOP_FD_COLORING_APPLY = (UInt32)(75)
const global MATOP_SET_FROM_OPTIONS = (UInt32)(76)
const global MATOP_MULT_CONSTRAINED = (UInt32)(77)
const global MATOP_MULT_TRANSPOSE_CONSTRAIN = (UInt32)(78)
const global MATOP_FIND_ZERO_DIAGONALS = (UInt32)(79)
const global MATOP_MULT_MULTIPLE = (UInt32)(80)
const global MATOP_SOLVE_MULTIPLE = (UInt32)(81)
const global MATOP_GET_INERTIA = (UInt32)(82)
const global MATOP_LOAD = (UInt32)(83)
const global MATOP_IS_SYMMETRIC = (UInt32)(84)
const global MATOP_IS_HERMITIAN = (UInt32)(85)
const global MATOP_IS_STRUCTURALLY_SYMMETRIC = (UInt32)(86)
const global MATOP_SET_VALUES_BLOCKEDLOCAL = (UInt32)(87)
const global MATOP_GET_VECS = (UInt32)(88)
const global MATOP_MAT_MULT = (UInt32)(89)
const global MATOP_MAT_MULT_SYMBOLIC = (UInt32)(90)
const global MATOP_MAT_MULT_NUMERIC = (UInt32)(91)
const global MATOP_PTAP = (UInt32)(92)
const global MATOP_PTAP_SYMBOLIC = (UInt32)(93)
const global MATOP_PTAP_NUMERIC = (UInt32)(94)
const global MATOP_MAT_TRANSPOSE_MULT = (UInt32)(95)
const global MATOP_MAT_TRANSPOSE_MULT_SYMBO = (UInt32)(96)
const global MATOP_MAT_TRANSPOSE_MULT_NUMER = (UInt32)(97)
const global MATOP_CONJUGATE = (UInt32)(102)
const global MATOP_SET_VALUES_ROW = (UInt32)(104)
const global MATOP_REAL_PART = (UInt32)(105)
const global MATOP_IMAGINARY_PART = (UInt32)(106)
const global MATOP_GET_ROW_UPPER_TRIANGULAR = (UInt32)(107)
const global MATOP_RESTORE_ROW_UPPER_TRIANG = (UInt32)(108)
const global MATOP_MAT_SOLVE = (UInt32)(109)
const global MATOP_GET_REDUNDANT_MATRIX = (UInt32)(110)
const global MATOP_GET_ROW_MIN = (UInt32)(111)
const global MATOP_GET_COLUMN_VECTOR = (UInt32)(112)
const global MATOP_MISSING_DIAGONAL = (UInt32)(113)
const global MATOP_GET_SEQ_NONZERO_STRUCTUR = (UInt32)(114)
const global MATOP_CREATE = (UInt32)(115)
const global MATOP_GET_GHOSTS = (UInt32)(116)
const global MATOP_GET_LOCAL_SUB_MATRIX = (UInt32)(117)
const global MATOP_RESTORE_LOCALSUB_MATRIX = (UInt32)(118)
const global MATOP_MULT_DIAGONAL_BLOCK = (UInt32)(119)
const global MATOP_HERMITIAN_TRANSPOSE = (UInt32)(120)
const global MATOP_MULT_HERMITIAN_TRANSPOSE = (UInt32)(121)
const global MATOP_MULT_HERMITIAN_TRANS_ADD = (UInt32)(122)
const global MATOP_GET_MULTI_PROC_BLOCK = (UInt32)(123)
const global MATOP_FIND_NONZERO_ROWS = (UInt32)(124)
const global MATOP_GET_COLUMN_NORMS = (UInt32)(125)
const global MATOP_INVERT_BLOCK_DIAGONAL = (UInt32)(126)
const global MATOP_GET_SUB_MATRICES_PARALLE = (UInt32)(128)
const global MATOP_SET_VALUES_BATCH = (UInt32)(129)
const global MATOP_TRANSPOSE_MAT_MULT = (UInt32)(130)
const global MATOP_TRANSPOSE_MAT_MULT_SYMBO = (UInt32)(131)
const global MATOP_TRANSPOSE_MAT_MULT_NUMER = (UInt32)(132)
const global MATOP_TRANSPOSE_COLORING_CREAT = (UInt32)(133)
const global MATOP_TRANS_COLORING_APPLY_SPT = (UInt32)(134)
const global MATOP_TRANS_COLORING_APPLY_DEN = (UInt32)(135)
const global MATOP_RART = (UInt32)(136)
const global MATOP_RART_SYMBOLIC = (UInt32)(137)
const global MATOP_RART_NUMERIC = (UInt32)(138)
const global MATOP_SET_BLOCK_SIZES = (UInt32)(139)
const global MATOP_AYPX = (UInt32)(140)
const global MATOP_RESIDUAL = (UInt32)(141)
const global MATOP_FDCOLORING_SETUP = (UInt32)(142)
const global MATOP_MPICONCATENATESEQ = (UInt32)(144)
# end enum ANONYMOUS_38

# begin enum MatOperation
typealias MatOperation Uint32
const global MATOP_SET_VALUES = (UInt32)(0)
const global MATOP_GET_ROW = (UInt32)(1)
const global MATOP_RESTORE_ROW = (UInt32)(2)
const global MATOP_MULT = (UInt32)(3)
const global MATOP_MULT_ADD = (UInt32)(4)
const global MATOP_MULT_TRANSPOSE = (UInt32)(5)
const global MATOP_MULT_TRANSPOSE_ADD = (UInt32)(6)
const global MATOP_SOLVE = (UInt32)(7)
const global MATOP_SOLVE_ADD = (UInt32)(8)
const global MATOP_SOLVE_TRANSPOSE = (UInt32)(9)
const global MATOP_SOLVE_TRANSPOSE_ADD = (UInt32)(10)
const global MATOP_LUFACTOR = (UInt32)(11)
const global MATOP_CHOLESKYFACTOR = (UInt32)(12)
const global MATOP_SOR = (UInt32)(13)
const global MATOP_TRANSPOSE = (UInt32)(14)
const global MATOP_GETINFO = (UInt32)(15)
const global MATOP_EQUAL = (UInt32)(16)
const global MATOP_GET_DIAGONAL = (UInt32)(17)
const global MATOP_DIAGONAL_SCALE = (UInt32)(18)
const global MATOP_NORM = (UInt32)(19)
const global MATOP_ASSEMBLY_BEGIN = (UInt32)(20)
const global MATOP_ASSEMBLY_END = (UInt32)(21)
const global MATOP_SET_OPTION = (UInt32)(22)
const global MATOP_ZERO_ENTRIES = (UInt32)(23)
const global MATOP_ZERO_ROWS = (UInt32)(24)
const global MATOP_LUFACTOR_SYMBOLIC = (UInt32)(25)
const global MATOP_LUFACTOR_NUMERIC = (UInt32)(26)
const global MATOP_CHOLESKY_FACTOR_SYMBOLIC = (UInt32)(27)
const global MATOP_CHOLESKY_FACTOR_NUMERIC = (UInt32)(28)
const global MATOP_SETUP_PREALLOCATION = (UInt32)(29)
const global MATOP_ILUFACTOR_SYMBOLIC = (UInt32)(30)
const global MATOP_ICCFACTOR_SYMBOLIC = (UInt32)(31)
const global MATOP_DUPLICATE = (UInt32)(34)
const global MATOP_FORWARD_SOLVE = (UInt32)(35)
const global MATOP_BACKWARD_SOLVE = (UInt32)(36)
const global MATOP_ILUFACTOR = (UInt32)(37)
const global MATOP_ICCFACTOR = (UInt32)(38)
const global MATOP_AXPY = (UInt32)(39)
const global MATOP_GET_SUBMATRICES = (UInt32)(40)
const global MATOP_INCREASE_OVERLAP = (UInt32)(41)
const global MATOP_GET_VALUES = (UInt32)(42)
const global MATOP_COPY = (UInt32)(43)
const global MATOP_GET_ROW_MAX = (UInt32)(44)
const global MATOP_SCALE = (UInt32)(45)
const global MATOP_SHIFT = (UInt32)(46)
const global MATOP_DIAGONAL_SET = (UInt32)(47)
const global MATOP_ZERO_ROWS_COLUMNS = (UInt32)(48)
const global MATOP_SET_RANDOM = (UInt32)(49)
const global MATOP_GET_ROW_IJ = (UInt32)(50)
const global MATOP_RESTORE_ROW_IJ = (UInt32)(51)
const global MATOP_GET_COLUMN_IJ = (UInt32)(52)
const global MATOP_RESTORE_COLUMN_IJ = (UInt32)(53)
const global MATOP_FDCOLORING_CREATE = (UInt32)(54)
const global MATOP_COLORING_PATCH = (UInt32)(55)
const global MATOP_SET_UNFACTORED = (UInt32)(56)
const global MATOP_PERMUTE = (UInt32)(57)
const global MATOP_SET_VALUES_BLOCKED = (UInt32)(58)
const global MATOP_GET_SUBMATRIX = (UInt32)(59)
const global MATOP_DESTROY = (UInt32)(60)
const global MATOP_VIEW = (UInt32)(61)
const global MATOP_CONVERT_FROM = (UInt32)(62)
const global MATOP_MATMAT_MULT = (UInt32)(63)
const global MATOP_MATMAT_MULT_SYMBOLIC = (UInt32)(64)
const global MATOP_MATMAT_MULT_NUMERIC = (UInt32)(65)
const global MATOP_SET_LOCAL_TO_GLOBAL_MAP = (UInt32)(66)
const global MATOP_SET_VALUES_LOCAL = (UInt32)(67)
const global MATOP_ZERO_ROWS_LOCAL = (UInt32)(68)
const global MATOP_GET_ROW_MAX_ABS = (UInt32)(69)
const global MATOP_GET_ROW_MIN_ABS = (UInt32)(70)
const global MATOP_CONVERT = (UInt32)(71)
const global MATOP_SET_COLORING = (UInt32)(72)
const global MATOP_SET_VALUES_ADIFOR = (UInt32)(74)
const global MATOP_FD_COLORING_APPLY = (UInt32)(75)
const global MATOP_SET_FROM_OPTIONS = (UInt32)(76)
const global MATOP_MULT_CONSTRAINED = (UInt32)(77)
const global MATOP_MULT_TRANSPOSE_CONSTRAIN = (UInt32)(78)
const global MATOP_FIND_ZERO_DIAGONALS = (UInt32)(79)
const global MATOP_MULT_MULTIPLE = (UInt32)(80)
const global MATOP_SOLVE_MULTIPLE = (UInt32)(81)
const global MATOP_GET_INERTIA = (UInt32)(82)
const global MATOP_LOAD = (UInt32)(83)
const global MATOP_IS_SYMMETRIC = (UInt32)(84)
const global MATOP_IS_HERMITIAN = (UInt32)(85)
const global MATOP_IS_STRUCTURALLY_SYMMETRIC = (UInt32)(86)
const global MATOP_SET_VALUES_BLOCKEDLOCAL = (UInt32)(87)
const global MATOP_GET_VECS = (UInt32)(88)
const global MATOP_MAT_MULT = (UInt32)(89)
const global MATOP_MAT_MULT_SYMBOLIC = (UInt32)(90)
const global MATOP_MAT_MULT_NUMERIC = (UInt32)(91)
const global MATOP_PTAP = (UInt32)(92)
const global MATOP_PTAP_SYMBOLIC = (UInt32)(93)
const global MATOP_PTAP_NUMERIC = (UInt32)(94)
const global MATOP_MAT_TRANSPOSE_MULT = (UInt32)(95)
const global MATOP_MAT_TRANSPOSE_MULT_SYMBO = (UInt32)(96)
const global MATOP_MAT_TRANSPOSE_MULT_NUMER = (UInt32)(97)
const global MATOP_CONJUGATE = (UInt32)(102)
const global MATOP_SET_VALUES_ROW = (UInt32)(104)
const global MATOP_REAL_PART = (UInt32)(105)
const global MATOP_IMAGINARY_PART = (UInt32)(106)
const global MATOP_GET_ROW_UPPER_TRIANGULAR = (UInt32)(107)
const global MATOP_RESTORE_ROW_UPPER_TRIANG = (UInt32)(108)
const global MATOP_MAT_SOLVE = (UInt32)(109)
const global MATOP_GET_REDUNDANT_MATRIX = (UInt32)(110)
const global MATOP_GET_ROW_MIN = (UInt32)(111)
const global MATOP_GET_COLUMN_VECTOR = (UInt32)(112)
const global MATOP_MISSING_DIAGONAL = (UInt32)(113)
const global MATOP_GET_SEQ_NONZERO_STRUCTUR = (UInt32)(114)
const global MATOP_CREATE = (UInt32)(115)
const global MATOP_GET_GHOSTS = (UInt32)(116)
const global MATOP_GET_LOCAL_SUB_MATRIX = (UInt32)(117)
const global MATOP_RESTORE_LOCALSUB_MATRIX = (UInt32)(118)
const global MATOP_MULT_DIAGONAL_BLOCK = (UInt32)(119)
const global MATOP_HERMITIAN_TRANSPOSE = (UInt32)(120)
const global MATOP_MULT_HERMITIAN_TRANSPOSE = (UInt32)(121)
const global MATOP_MULT_HERMITIAN_TRANS_ADD = (UInt32)(122)
const global MATOP_GET_MULTI_PROC_BLOCK = (UInt32)(123)
const global MATOP_FIND_NONZERO_ROWS = (UInt32)(124)
const global MATOP_GET_COLUMN_NORMS = (UInt32)(125)
const global MATOP_INVERT_BLOCK_DIAGONAL = (UInt32)(126)
const global MATOP_GET_SUB_MATRICES_PARALLE = (UInt32)(128)
const global MATOP_SET_VALUES_BATCH = (UInt32)(129)
const global MATOP_TRANSPOSE_MAT_MULT = (UInt32)(130)
const global MATOP_TRANSPOSE_MAT_MULT_SYMBO = (UInt32)(131)
const global MATOP_TRANSPOSE_MAT_MULT_NUMER = (UInt32)(132)
const global MATOP_TRANSPOSE_COLORING_CREAT = (UInt32)(133)
const global MATOP_TRANS_COLORING_APPLY_SPT = (UInt32)(134)
const global MATOP_TRANS_COLORING_APPLY_DEN = (UInt32)(135)
const global MATOP_RART = (UInt32)(136)
const global MATOP_RART_SYMBOLIC = (UInt32)(137)
const global MATOP_RART_NUMERIC = (UInt32)(138)
const global MATOP_SET_BLOCK_SIZES = (UInt32)(139)
const global MATOP_AYPX = (UInt32)(140)
const global MATOP_RESIDUAL = (UInt32)(141)
const global MATOP_FDCOLORING_SETUP = (UInt32)(142)
const global MATOP_MPICONCATENATESEQ = (UInt32)(144)
# end enum MatOperation

typealias MatNullSpace Ptr{Void}
typealias MatMFFD Ptr{Void}
typealias MatMFFDType Ptr{Uint8}
typealias DM Ptr{Void}

# begin enum ANONYMOUS_39
typealias ANONYMOUS_39 Uint32
const global DM_BOUNDARY_NONE = (UInt32)(0)
const global DM_BOUNDARY_GHOSTED = (UInt32)(1)
const global DM_BOUNDARY_MIRROR = (UInt32)(2)
const global DM_BOUNDARY_PERIODIC = (UInt32)(3)
const global DM_BOUNDARY_TWIST = (UInt32)(4)
# end enum ANONYMOUS_39

# begin enum DMBoundaryType
typealias DMBoundaryType Uint32
const global DM_BOUNDARY_NONE = (UInt32)(0)
const global DM_BOUNDARY_GHOSTED = (UInt32)(1)
const global DM_BOUNDARY_MIRROR = (UInt32)(2)
const global DM_BOUNDARY_PERIODIC = (UInt32)(3)
const global DM_BOUNDARY_TWIST = (UInt32)(4)
# end enum DMBoundaryType

typealias PetscPartitioner Ptr{Void}
typealias PC Ptr{Void}
typealias PCType Ptr{Uint8}

# begin enum ANONYMOUS_40
typealias ANONYMOUS_40 Cint
const global PC_SIDE_DEFAULT = (Int32)(-1)
const global PC_LEFT = (Int32)(0)
const global PC_RIGHT = (Int32)(1)
const global PC_SYMMETRIC = (Int32)(2)
# end enum ANONYMOUS_40

# begin enum ANONYMOUS_41
typealias ANONYMOUS_41 Cint
const global PCRICHARDSON_CONVERGED_RTOL = (Int32)(2)
const global PCRICHARDSON_CONVERGED_ATOL = (Int32)(3)
const global PCRICHARDSON_CONVERGED_ITS = (Int32)(4)
const global PCRICHARDSON_DIVERGED_DTOL = (Int32)(-4)
# end enum ANONYMOUS_41

# begin enum PCRichardsonConvergedReason
typealias PCRichardsonConvergedReason Cint
const global PCRICHARDSON_CONVERGED_RTOL = (Int32)(2)
const global PCRICHARDSON_CONVERGED_ATOL = (Int32)(3)
const global PCRICHARDSON_CONVERGED_ITS = (Int32)(4)
const global PCRICHARDSON_DIVERGED_DTOL = (Int32)(-4)
# end enum PCRichardsonConvergedReason

# begin enum ANONYMOUS_42
typealias ANONYMOUS_42 Uint32
const global PC_JACOBI_DIAGONAL = (UInt32)(0)
const global PC_JACOBI_ROWMAX = (UInt32)(1)
const global PC_JACOBI_ROWSUM = (UInt32)(2)
# end enum ANONYMOUS_42

# begin enum PCJacobiType
typealias PCJacobiType Uint32
const global PC_JACOBI_DIAGONAL = (UInt32)(0)
const global PC_JACOBI_ROWMAX = (UInt32)(1)
const global PC_JACOBI_ROWSUM = (UInt32)(2)
# end enum PCJacobiType

# begin enum ANONYMOUS_43
typealias ANONYMOUS_43 Uint32
const global PC_ASM_BASIC = (UInt32)(3)
const global PC_ASM_RESTRICT = (UInt32)(1)
const global PC_ASM_INTERPOLATE = (UInt32)(2)
const global PC_ASM_NONE = (UInt32)(0)
# end enum ANONYMOUS_43

# begin enum PCASMType
typealias PCASMType Uint32
const global PC_ASM_BASIC = (UInt32)(3)
const global PC_ASM_RESTRICT = (UInt32)(1)
const global PC_ASM_INTERPOLATE = (UInt32)(2)
const global PC_ASM_NONE = (UInt32)(0)
# end enum PCASMType

# begin enum ANONYMOUS_44
typealias ANONYMOUS_44 Uint32
const global PC_GASM_BASIC = (UInt32)(3)
const global PC_GASM_RESTRICT = (UInt32)(1)
const global PC_GASM_INTERPOLATE = (UInt32)(2)
const global PC_GASM_NONE = (UInt32)(0)
# end enum ANONYMOUS_44

# begin enum PCGASMType
typealias PCGASMType Uint32
const global PC_GASM_BASIC = (UInt32)(3)
const global PC_GASM_RESTRICT = (UInt32)(1)
const global PC_GASM_INTERPOLATE = (UInt32)(2)
const global PC_GASM_NONE = (UInt32)(0)
# end enum PCGASMType

# begin enum ANONYMOUS_45
typealias ANONYMOUS_45 Uint32
const global PC_COMPOSITE_ADDITIVE = (UInt32)(0)
const global PC_COMPOSITE_MULTIPLICATIVE = (UInt32)(1)
const global PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE = (UInt32)(2)
const global PC_COMPOSITE_SPECIAL = (UInt32)(3)
const global PC_COMPOSITE_SCHUR = (UInt32)(4)
# end enum ANONYMOUS_45

# begin enum PCCompositeType
typealias PCCompositeType Uint32
const global PC_COMPOSITE_ADDITIVE = (UInt32)(0)
const global PC_COMPOSITE_MULTIPLICATIVE = (UInt32)(1)
const global PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE = (UInt32)(2)
const global PC_COMPOSITE_SPECIAL = (UInt32)(3)
const global PC_COMPOSITE_SCHUR = (UInt32)(4)
# end enum PCCompositeType

# begin enum ANONYMOUS_46
typealias ANONYMOUS_46 Uint32
const global PC_FIELDSPLIT_SCHUR_PRE_SELF = (UInt32)(0)
const global PC_FIELDSPLIT_SCHUR_PRE_SELFP = (UInt32)(1)
const global PC_FIELDSPLIT_SCHUR_PRE_A11 = (UInt32)(2)
const global PC_FIELDSPLIT_SCHUR_PRE_USER = (UInt32)(3)
const global PC_FIELDSPLIT_SCHUR_PRE_FULL = (UInt32)(4)
# end enum ANONYMOUS_46

# begin enum PCFieldSplitSchurPreType
typealias PCFieldSplitSchurPreType Uint32
const global PC_FIELDSPLIT_SCHUR_PRE_SELF = (UInt32)(0)
const global PC_FIELDSPLIT_SCHUR_PRE_SELFP = (UInt32)(1)
const global PC_FIELDSPLIT_SCHUR_PRE_A11 = (UInt32)(2)
const global PC_FIELDSPLIT_SCHUR_PRE_USER = (UInt32)(3)
const global PC_FIELDSPLIT_SCHUR_PRE_FULL = (UInt32)(4)
# end enum PCFieldSplitSchurPreType

# begin enum ANONYMOUS_47
typealias ANONYMOUS_47 Uint32
const global PC_FIELDSPLIT_SCHUR_FACT_DIAG = (UInt32)(0)
const global PC_FIELDSPLIT_SCHUR_FACT_LOWER = (UInt32)(1)
const global PC_FIELDSPLIT_SCHUR_FACT_UPPER = (UInt32)(2)
const global PC_FIELDSPLIT_SCHUR_FACT_FULL = (UInt32)(3)
# end enum ANONYMOUS_47

# begin enum PCFieldSplitSchurFactType
typealias PCFieldSplitSchurFactType Uint32
const global PC_FIELDSPLIT_SCHUR_FACT_DIAG = (UInt32)(0)
const global PC_FIELDSPLIT_SCHUR_FACT_LOWER = (UInt32)(1)
const global PC_FIELDSPLIT_SCHUR_FACT_UPPER = (UInt32)(2)
const global PC_FIELDSPLIT_SCHUR_FACT_FULL = (UInt32)(3)
# end enum PCFieldSplitSchurFactType

# begin enum ANONYMOUS_48
typealias ANONYMOUS_48 Uint32
const global PC_PARMS_GLOBAL_RAS = (UInt32)(0)
const global PC_PARMS_GLOBAL_SCHUR = (UInt32)(1)
const global PC_PARMS_GLOBAL_BJ = (UInt32)(2)
# end enum ANONYMOUS_48

# begin enum PCPARMSGlobalType
typealias PCPARMSGlobalType Uint32
const global PC_PARMS_GLOBAL_RAS = (UInt32)(0)
const global PC_PARMS_GLOBAL_SCHUR = (UInt32)(1)
const global PC_PARMS_GLOBAL_BJ = (UInt32)(2)
# end enum PCPARMSGlobalType

# begin enum ANONYMOUS_49
typealias ANONYMOUS_49 Uint32
const global PC_PARMS_LOCAL_ILU0 = (UInt32)(0)
const global PC_PARMS_LOCAL_ILUK = (UInt32)(1)
const global PC_PARMS_LOCAL_ILUT = (UInt32)(2)
const global PC_PARMS_LOCAL_ARMS = (UInt32)(3)
# end enum ANONYMOUS_49

# begin enum PCPARMSLocalType
typealias PCPARMSLocalType Uint32
const global PC_PARMS_LOCAL_ILU0 = (UInt32)(0)
const global PC_PARMS_LOCAL_ILUK = (UInt32)(1)
const global PC_PARMS_LOCAL_ILUT = (UInt32)(2)
const global PC_PARMS_LOCAL_ARMS = (UInt32)(3)
# end enum PCPARMSLocalType

typealias PCGAMGType Ptr{Uint8}
typealias PCGAMGClassicalType Ptr{Uint8}

# begin enum ANONYMOUS_50
typealias ANONYMOUS_50 Uint32
const global PC_MG_MULTIPLICATIVE = (UInt32)(0)
const global PC_MG_ADDITIVE = (UInt32)(1)
const global PC_MG_FULL = (UInt32)(2)
const global PC_MG_KASKADE = (UInt32)(3)
# end enum ANONYMOUS_50

# begin enum ANONYMOUS_51
typealias ANONYMOUS_51 Uint32
const global PC_MG_CYCLE_V = (UInt32)(1)
const global PC_MG_CYCLE_W = (UInt32)(2)
# end enum ANONYMOUS_51

# begin enum PCMGCycleType
typealias PCMGCycleType Uint32
const global PC_MG_CYCLE_V = (UInt32)(1)
const global PC_MG_CYCLE_W = (UInt32)(2)
# end enum PCMGCycleType

# begin enum ANONYMOUS_52
typealias ANONYMOUS_52 Uint32
const global PC_EXOTIC_FACE = (UInt32)(0)
const global PC_EXOTIC_WIREBASKET = (UInt32)(1)
# end enum ANONYMOUS_52

# begin enum PCExoticType
typealias PCExoticType Uint32
const global PC_EXOTIC_FACE = (UInt32)(0)
const global PC_EXOTIC_WIREBASKET = (UInt32)(1)
# end enum PCExoticType

typealias KSP Ptr{Void}
typealias KSPType Ptr{Uint8}

# begin enum ANONYMOUS_53
typealias ANONYMOUS_53 Uint32
const global KSP_FCG_TRUNC_TYPE_STANDARD = (UInt32)(0)
const global KSP_FCG_TRUNC_TYPE_NOTAY = (UInt32)(1)
# end enum ANONYMOUS_53

# begin enum KSPFCGTruncationType
typealias KSPFCGTruncationType Uint32
const global KSP_FCG_TRUNC_TYPE_STANDARD = (UInt32)(0)
const global KSP_FCG_TRUNC_TYPE_NOTAY = (UInt32)(1)
# end enum KSPFCGTruncationType

# begin enum ANONYMOUS_54
typealias ANONYMOUS_54 Uint32
const global KSP_GMRES_CGS_REFINE_NEVER = (UInt32)(0)
const global KSP_GMRES_CGS_REFINE_IFNEEDED = (UInt32)(1)
const global KSP_GMRES_CGS_REFINE_ALWAYS = (UInt32)(2)
# end enum ANONYMOUS_54

# begin enum KSPGMRESCGSRefinementType
typealias KSPGMRESCGSRefinementType Uint32
const global KSP_GMRES_CGS_REFINE_NEVER = (UInt32)(0)
const global KSP_GMRES_CGS_REFINE_IFNEEDED = (UInt32)(1)
const global KSP_GMRES_CGS_REFINE_ALWAYS = (UInt32)(2)
# end enum KSPGMRESCGSRefinementType

# begin enum ANONYMOUS_55
typealias ANONYMOUS_55 Cint
const global KSP_NORM_DEFAULT = (Int32)(-1)
const global KSP_NORM_NONE = (Int32)(0)
const global KSP_NORM_PRECONDITIONED = (Int32)(1)
const global KSP_NORM_UNPRECONDITIONED = (Int32)(2)
const global KSP_NORM_NATURAL = (Int32)(3)
# end enum ANONYMOUS_55

# begin enum ANONYMOUS_56
typealias ANONYMOUS_56 Cint
const global KSP_CONVERGED_RTOL_NORMAL = (Int32)(1)
const global KSP_CONVERGED_ATOL_NORMAL = (Int32)(9)
const global KSP_CONVERGED_RTOL = (Int32)(2)
const global KSP_CONVERGED_ATOL = (Int32)(3)
const global KSP_CONVERGED_ITS = (Int32)(4)
const global KSP_CONVERGED_CG_NEG_CURVE = (Int32)(5)
const global KSP_CONVERGED_CG_CONSTRAINED = (Int32)(6)
const global KSP_CONVERGED_STEP_LENGTH = (Int32)(7)
const global KSP_CONVERGED_HAPPY_BREAKDOWN = (Int32)(8)
const global KSP_DIVERGED_NULL = (Int32)(-2)
const global KSP_DIVERGED_ITS = (Int32)(-3)
const global KSP_DIVERGED_DTOL = (Int32)(-4)
const global KSP_DIVERGED_BREAKDOWN = (Int32)(-5)
const global KSP_DIVERGED_BREAKDOWN_BICG = (Int32)(-6)
const global KSP_DIVERGED_NONSYMMETRIC = (Int32)(-7)
const global KSP_DIVERGED_INDEFINITE_PC = (Int32)(-8)
const global KSP_DIVERGED_NANORINF = (Int32)(-9)
const global KSP_DIVERGED_INDEFINITE_MAT = (Int32)(-10)
const global KSP_DIVERGED_PCSETUP_FAILED = (Int32)(-11)
const global KSP_CONVERGED_ITERATING = (Int32)(0)
# end enum ANONYMOUS_56

# begin enum KSPConvergedReason
typealias KSPConvergedReason Cint
const global KSP_CONVERGED_RTOL_NORMAL = (Int32)(1)
const global KSP_CONVERGED_ATOL_NORMAL = (Int32)(9)
const global KSP_CONVERGED_RTOL = (Int32)(2)
const global KSP_CONVERGED_ATOL = (Int32)(3)
const global KSP_CONVERGED_ITS = (Int32)(4)
const global KSP_CONVERGED_CG_NEG_CURVE = (Int32)(5)
const global KSP_CONVERGED_CG_CONSTRAINED = (Int32)(6)
const global KSP_CONVERGED_STEP_LENGTH = (Int32)(7)
const global KSP_CONVERGED_HAPPY_BREAKDOWN = (Int32)(8)
const global KSP_DIVERGED_NULL = (Int32)(-2)
const global KSP_DIVERGED_ITS = (Int32)(-3)
const global KSP_DIVERGED_DTOL = (Int32)(-4)
const global KSP_DIVERGED_BREAKDOWN = (Int32)(-5)
const global KSP_DIVERGED_BREAKDOWN_BICG = (Int32)(-6)
const global KSP_DIVERGED_NONSYMMETRIC = (Int32)(-7)
const global KSP_DIVERGED_INDEFINITE_PC = (Int32)(-8)
const global KSP_DIVERGED_NANORINF = (Int32)(-9)
const global KSP_DIVERGED_INDEFINITE_MAT = (Int32)(-10)
const global KSP_DIVERGED_PCSETUP_FAILED = (Int32)(-11)
const global KSP_CONVERGED_ITERATING = (Int32)(0)
# end enum KSPConvergedReason

# begin enum ANONYMOUS_57
typealias ANONYMOUS_57 Uint32
const global KSP_CG_SYMMETRIC = (UInt32)(0)
const global KSP_CG_HERMITIAN = (UInt32)(1)
# end enum ANONYMOUS_57

# begin enum KSPCGType
typealias KSPCGType Uint32
const global KSP_CG_SYMMETRIC = (UInt32)(0)
const global KSP_CG_HERMITIAN = (UInt32)(1)
# end enum KSPCGType

typealias KSPFischerGuess Ptr{Void}

# begin enum ANONYMOUS_58
typealias ANONYMOUS_58 Uint32
const global MAT_SCHUR_COMPLEMENT_AINV_DIAG = (UInt32)(0)
const global MAT_SCHUR_COMPLEMENT_AINV_LUMP = (UInt32)(1)
# end enum ANONYMOUS_58

# begin enum MatSchurComplementAinvType
typealias MatSchurComplementAinvType Uint32
const global MAT_SCHUR_COMPLEMENT_AINV_DIAG = (UInt32)(0)
const global MAT_SCHUR_COMPLEMENT_AINV_LUMP = (UInt32)(1)
# end enum MatSchurComplementAinvType
