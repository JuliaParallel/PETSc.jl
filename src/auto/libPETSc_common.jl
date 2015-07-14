# Automatically generated using Clang.jl wrap_c, version 0.0.0

using Compat
using MPI
#=
const PETSC_FUNCTION_NAME = PETSC_FUNCTION_NAME_C
const PETSC_RESTRICT = PETSC_C_RESTRICT
const PETSC_STATIC_INLINE = PETSC_C_STATIC_INLINE
const PETSC_VISIBILITY_PUBLIC = PETSC_DLLIMPORT
=#
# Skipping MacroDefinition: PETSC_EXTERN extern PETSC_VISIBILITY_PUBLIC
# Skipping MacroDefinition: PETSC_INTERN extern PETSC_VISIBILITY_INTERNAL

const PETSC_VERSION_RELEASE = 1
const PETSC_VERSION_MAJOR = 3
const PETSC_VERSION_MINOR = 6
const PETSC_VERSION_SUBMINOR = 0
const PETSC_VERSION_PATCH = 0
const PETSC_RELEASE_DATE = "Jun, 9, 2015"
const PETSC_VERSION_DATE = "Jun, 09, 2015"
const PETSC_VERSION_GIT = "v3.6"
const PETSC_VERSION_DATE_GIT = "2015-06-09 16:15:46 -0500"

# Skipping MacroDefinition: PETSC_VERSION_ ( MAJOR , MINOR , SUBMINOR ) ( ( PETSC_VERSION_MAJOR == ( MAJOR ) ) && ( PETSC_VERSION_MINOR == ( MINOR ) ) && ( PETSC_VERSION_SUBMINOR == ( SUBMINOR ) ) && ( PETSC_VERSION_RELEASE == 1 ) )
# Skipping MacroDefinition: PETSC_VERSION_LT ( MAJOR , MINOR , SUBMINOR ) ( PETSC_VERSION_RELEASE == 1 && ( PETSC_VERSION_MAJOR < ( MAJOR ) || ( PETSC_VERSION_MAJOR == ( MAJOR ) && ( PETSC_VERSION_MINOR < ( MINOR ) || ( PETSC_VERSION_MINOR == ( MINOR ) && ( PETSC_VERSION_SUBMINOR < ( SUBMINOR ) ) ) ) ) ) )
# Skipping MacroDefinition: PETSC_VERSION_LE ( MAJOR , MINOR , SUBMINOR ) ( PETSC_VERSION_LT ( MAJOR , MINOR , SUBMINOR ) || PETSC_VERSION_ ( MAJOR , MINOR , SUBMINOR ) )
# Skipping MacroDefinition: PETSC_VERSION_GT ( MAJOR , MINOR , SUBMINOR ) ( 0 == PETSC_VERSION_LE ( MAJOR , MINOR , SUBMINOR ) )
# Skipping MacroDefinition: PETSC_VERSION_GE ( MAJOR , MINOR , SUBMINOR ) ( 0 == PETSC_VERSION_LT ( MAJOR , MINOR , SUBMINOR ) )

const PETSC_AUTHOR_INFO = "       The PETSc Team\n    petsc-maint@mcs.anl.gov\n http://www.mcs.anl.gov/petsc/\n"
const MPICH_SKIP_MPICXX = 1
const OMPI_SKIP_MPICXX = 1

# Skipping MacroDefinition: PetscAttrMPIPointerWithType ( bufno , typeno ) __attribute__ ( ( pointer_with_type_tag ( MPI , bufno , typeno ) ) )
# Skipping MacroDefinition: PetscAttrMPITypeTag ( type ) __attribute__ ( ( type_tag_for_datatype ( MPI , type ) ) )
# Skipping MacroDefinition: PetscAttrMPITypeTagLayoutCompatible ( type ) __attribute__ ( ( type_tag_for_datatype ( MPI , type , layout_compatible ) ) )

const MPIU_INT = Int32
const MPIU_INT64 = Int64
const MPIU_SIZE_T = UInt64

# Skipping MacroDefinition: PetscUnlikely ( cond ) ( cond )
# Skipping MacroDefinition: PetscLikely ( cond ) ( cond )
# Skipping MacroDefinition: PetscExpPassiveScalar ( a ) PetscExpScalar ( )

const MPIU_SCALAR = Int32

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

const PETSC_PI = pi
const PETSC_MAX_INT = 2147483647
const PETSC_MIN_INT = -PETSC_MAX_INT - 1
const PETSC_MAX_REAL = PETSC_MAX_INT  # made up
const PETSC_INFINITY = PETSC_MAX_REAL / 4.0
const PETSC_NINFINITY = -PETSC_INFINITY
const PassiveReal = Cint

typealias PetscScalar Cint

const PassiveScalar = PetscScalar
const MPIU_MATSCALAR = MPIU_SCALAR
const MPIU_2INT = Cint
const PETSC_NULL = C_NULL
const PETSC_IGNORE = C_NULL
const PETSC_DECIDE = -1
const PETSC_DETERMINE = PETSC_DECIDE
const PETSC_DEFAULT = -2
const PETSC_COMM_SELF = MPI.COMM_SELF

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

const MPIU_PETSCLOGDOUBLE = Cdouble

# begin enum PetscDataType
typealias PetscDataType Uint32
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
# end enum PetscDataType

const PETSC_SCALAR = PETSC_DOUBLE
const PETSC_REAL = PETSC_DOUBLE
const PETSC_FORTRANADDR = PETSC_LONG
const MPIU_SUM = MPI.SUM
const MPIU_MAX = MPI.MAX
const MPIU_MIN = MPI.MIN
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

const PETSCSTACKSIZE = 64

# Skipping MacroDefinition: PetscStackPushNoCheck ( funct , petsc_routine , hot ) do { } while ( 0 )
# Skipping MacroDefinition: PetscStackPopNoCheck do { } while ( 0 )
# Skipping MacroDefinition: PetscFunctionReturn ( a ) return ( a )
# Skipping MacroDefinition: PetscFunctionReturnVoid ( ) return

# what is this?
#const PetscStackPop = CHKMEMQ

# Skipping MacroDefinition: PetscStackPush ( f ) CHKMEMQ
# Skipping MacroDefinition: PetscStackCall ( name , routine ) do { PetscStackPush ( name ) ; routine ; PetscStackPop ; } while ( 0 )
# Skipping MacroDefinition: PetscStackCallStandard ( func , args ) do { PetscStackPush ( # func ) ; ierr = func args ; PetscStackPop ; if ( ierr ) SETERRQ1 ( PETSC_COMM_SELF , PETSC_ERR_LIB , "Error in %s()" , # func ) ; } while ( 0 )

const PETSC_SMALLEST_CLASSID = 1211211

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

const PETSC_EVENT = 1311311

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

const PetscLogPLB = 0
const PetscLogPLE = 0
const PetscLogPHC = 0
const PetscLogPHD = 0

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
#} while ( 0 )
# Skipping MacroDefinition: PetscPreLoadStage ( name ) do { _3_ierr = PetscLogStagePop ( ) ; CHKERRQ ( _3_ierr ) ; if ( PetscPreLoadIt > 0 ) { _3_ierr = PetscLogStageGetId ( name , & _stageNum ) ; CHKERRQ ( _3_ierr ) ; } else { _3_ierr = PetscLogStageRegister ( name , & _stageNum ) ; CHKERRQ ( _3_ierr ) ; } _3_ierr = PetscLogStageSetActive ( _stageNum , ( PetscBool ) ( ! PetscPreLoadMax || PetscPreLoadIt ) ) ; _3_ierr = PetscLogStagePush ( _stageNum ) ; CHKERRQ ( _3_ierr ) ; } while ( 0 )
# Skipping MacroDefinition: PetscPrefetchBlock ( a , n , rw , t ) do { const char * _p = ( const char * ) ( a ) , * _end = ( const char * ) ( ( a ) + ( n ) ) ; for ( ; _p < _end ; _p += PETSC_LEVEL1_DCACHE_LINESIZE ) PETSC_Prefetch ( _p , ( rw ) , ( t ) ) ; } while ( 0 )

const PETSC_MPI_INT_MAX = 2147483647
const PETSC_MPI_INT_MIN = -2147483647
const PETSC_BLAS_INT_MAX = 2147483647
const PETSC_BLAS_INT_MIN = -2147483647
const PETSC_MAX_PATH_LEN = 4096
const PETSCRAND = "rand"
const PETSCRAND48 = "rand48"
const PETSCSPRNG = "sprng"
const PETSC_BINARY_INT_SIZE = 32 / 8
const PETSC_BINARY_FLOAT_SIZE = 32 / 8
const PETSC_BINARY_CHAR_SIZE = 8 / 8
const PETSC_BINARY_SHORT_SIZE = 16 / 8
const PETSC_BINARY_DOUBLE_SIZE = 64 / 8

# Skipping MacroDefinition: PETSC_BINARY_SCALAR_SIZE sizeof ( PetscScalar )

const PETSC_BAG_FILE_CLASSID = 1211219
const PETSCVIEWERSOCKET = "socket"
const PETSCVIEWERASCII = "ascii"
const PETSCVIEWERBINARY = "binary"
const PETSCVIEWERSTRING = "string"
const PETSCVIEWERDRAW = "draw"
const PETSCVIEWERVU = "vu"
const PETSCVIEWERMATHEMATICA = "mathematica"
const PETSCVIEWERNETCDF = "netcdf"
const PETSCVIEWERHDF5 = "hdf5"
const PETSCVIEWERVTK = "vtk"
const PETSCVIEWERMATLAB = "matlab"
const PETSCVIEWERSAWS = "saws"
const PETSC_DRAW_X = "x"
const PETSC_DRAW_GLUT = "glut"
const PETSC_DRAW_OPENGLES = "opengles"
const PETSC_DRAW_NULL = "null"
const PETSC_DRAW_WIN32 = "win32"
const PETSC_DRAW_TIKZ = "tikz"

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

const PETSC_HASH_FACT = 79943

# Skipping MacroDefinition: PETSC_MATLAB_ENGINE_WORLD PETSC_MATLAB_ENGINE_ ( PETSC_COMM_WORLD )
# Skipping MacroDefinition: PETSC_MATLAB_ENGINE_SELF PETSC_MATLAB_ENGINE_ ( PETSC_COMM_SELF )

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
const ISGENERAL = "general"
const ISSTRIDE = "stride"
const ISBLOCK = "block"

# vector formats
const VECSEQ = "seq"
const VECMPI = "mpi"
const VECSTANDARD = "standard"
const VECSHARED = "shared"
const VECSEQCUSP = "seqcusp"
const VECMPICUSP = "mpicusp"
const VECCUSP = "cusp"
const VECSEQVIENNACL = "seqviennacl"
const VECMPIVIENNACL = "mpiviennacl"
const VECVIENNACL = "viennacl"
const VECNEST = "nest"
const VECSEQPTHREAD = "seqpthread"
const VECMPIPTHREAD = "mpipthread"
const VECPTHREAD = "pthread"


const VEC_FILE_CLASSID = 1211214

# begin enum NormType
typealias NormType Uint32
const NORM_1 = (UInt32)(0)
const NORM_2 = (UInt32)(1)
const NORM_FROBENIUS = (UInt32)(2)
const NORM_INFINITY = (UInt32)(3)
const NORM_1_AND_2 = (UInt32)(4)
# end enum NormType

const NORM_MAX = NORM_INFINITY

# Skipping MacroDefinition: VecLockGet ( x , arg ) * ( arg ) = 0
# Skipping MacroDefinition: VecLockPush ( x ) 0
# Skipping MacroDefinition: VecLockPop ( x ) 0
# Skipping MacroDefinition: VecLocked ( x , arg )

# matrix formats
const MATSAME = "same"
const MATMAIJ = "maij"
const MATSEQMAIJ = "seqmaij"
const MATMPIMAIJ = "mpimaij"
const MATIS = "is"
const MATAIJ = "aij"
const MATSEQAIJ = "seqaij"
const MATSEQAIJPTHREAD = "seqaijpthread"
const MATAIJPTHREAD = "aijpthread"
const MATMPIAIJ = "mpiaij"
const MATAIJCRL = "aijcrl"
const MATSEQAIJCRL = "seqaijcrl"
const MATMPIAIJCRL = "mpiaijcrl"
const MATAIJCUSP = "aijcusp"
const MATSEQAIJCUSP = "seqaijcusp"
const MATMPIAIJCUSP = "mpiaijcusp"
const MATAIJCUSPARSE = "aijcusparse"
const MATSEQAIJCUSPARSE = "seqaijcusparse"
const MATMPIAIJCUSPARSE = "mpiaijcusparse"
const MATAIJVIENNACL = "aijviennacl"
const MATSEQAIJVIENNACL = "seqaijviennacl"
const MATMPIAIJVIENNACL = "mpiaijviennacl"
const MATAIJPERM = "aijperm"
const MATSEQAIJPERM = "seqaijperm"
const MATMPIAIJPERM = "mpiaijperm"
const MATSHELL = "shell"
const MATDENSE = "dense"
const MATSEQDENSE = "seqdense"
const MATMPIDENSE = "mpidense"
const MATELEMENTAL = "elemental"
const MATBAIJ = "baij"
const MATSEQBAIJ = "seqbaij"
const MATMPIBAIJ = "mpibaij"
const MATMPIADJ = "mpiadj"
const MATSBAIJ = "sbaij"
const MATSEQSBAIJ = "seqsbaij"
const MATMPISBAIJ = "mpisbaij"
const MATSEQBSTRM = "seqbstrm"
const MATMPIBSTRM = "mpibstrm"
const MATBSTRM = "bstrm"
const MATSEQSBSTRM = "seqsbstrm"
const MATMPISBSTRM = "mpisbstrm"
const MATSBSTRM = "sbstrm"
const MATDAAD = "daad"
const MATMFFD = "mffd"
const MATNORMAL = "normal"
const MATLRC = "lrc"
const MATSCATTER = "scatter"
const MATBLOCKMAT = "blockmat"
const MATCOMPOSITE = "composite"
const MATFFT = "fft"
const MATFFTW = "fftw"
const MATSEQCUFFT = "seqcufft"
const MATTRANSPOSEMAT = "transpose"
const MATSCHURCOMPLEMENT = "schurcomplement"
const MATPYTHON = "python"
const MATHYPRESTRUCT = "hyprestruct"
const MATHYPRESSTRUCT = "hypresstruct"
const MATSUBMATRIX = "submatrix"
const MATLOCALREF = "localref"
const MATNEST = "nest"

# Skipping MacroDefinition: MatSolverPackage char *

const MATSOLVERSUPERLU = "superlu"
const MATSOLVERSUPERLU_DIST = "superlu_dist"
const MATSOLVERUMFPACK = "umfpack"
const MATSOLVERCHOLMOD = "cholmod"
const MATSOLVERESSL = "essl"
const MATSOLVERLUSOL = "lusol"
const MATSOLVERMUMPS = "mumps"
const MATSOLVERMKL_PARDISO = "mkl_pardiso"
const MATSOLVERMKL_CPARDISO = "mkl_cpardiso"
const MATSOLVERPASTIX = "pastix"
const MATSOLVERMATLAB = "matlab"
const MATSOLVERPETSC = "petsc"
const MATSOLVERBAS = "bas"
const MATSOLVERCUSPARSE = "cusparse"
const MATSOLVERBSTRM = "bstrm"
const MATSOLVERSBSTRM = "sbstrm"
const MATSOLVERELEMENTAL = "elemental"
const MATSOLVERCLIQUE = "clique"
const MATSOLVERKLU = "klu"
const MAT_FILE_CLASSID = 1211216


#=
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
=#
const MAT_SKIP_ALLOCATION = -4
const MATORDERINGNATURAL = "natural"
const MATORDERINGND = "nd"
const MATORDERING1WD = "1wd"
const MATORDERINGRCM = "rcm"
const MATORDERINGQMD = "qmd"
const MATORDERINGROWLENGTH = "rowlength"
const MATORDERINGWBM = "wbm"
const MATORDERINGSPECTRAL = "spectral"
const MATORDERINGAMD = "amd"
const MATCOLORINGJP = "jp"
const MATCOLORINGPOWER = "power"
const MATCOLORINGNATURAL = "natural"
const MATCOLORINGSL = "sl"
const MATCOLORINGLF = "lf"
const MATCOLORINGID = "id"
const MATCOLORINGGREEDY = "greedy"
const MATPARTITIONINGCURRENT = "current"
const MATPARTITIONINGSQUARE = "square"
const MATPARTITIONINGPARMETIS = "parmetis"
const MATPARTITIONINGCHACO = "chaco"
const MATPARTITIONINGPARTY = "party"
const MATPARTITIONINGPTSCOTCH = "ptscotch"
const MP_PARTY_OPT = "opt"
const MP_PARTY_LIN = "lin"
const MP_PARTY_SCA = "sca"
const MP_PARTY_RAN = "ran"
const MP_PARTY_GBF = "gbf"
const MP_PARTY_GCF = "gcf"
const MP_PARTY_BUB = "bub"
const MP_PARTY_DEF = "def"
const MP_PARTY_HELPFUL_SETS = "hs"
const MP_PARTY_KERNIGHAN_LIN = "kl"
const MP_PARTY_NONE = "no"
const MATCOARSENMIS = "mis"
const MATCOARSENHEM = "hem"
const MATRIX_BINARY_FORMAT_DENSE = -1
const MATMFFD_DS = "ds"
const MATMFFD_WP = "wp"
const DMDA = "da"
const DMCOMPOSITE = "composite"
const DMSLICED = "sliced"
const DMSHELL = "shell"
const DMPLEX = "plex"
const DMCARTESIAN = "cartesian"
const DMREDUNDANT = "redundant"
const DMPATCH = "patch"
const DMMOAB = "moab"
const DMNETWORK = "network"
const DM_FILE_CLASSID = 1211221
const PFCONSTANT = "constant"
const PFMAT = "mat"
const PFSTRING = "string"
const PFQUICK = "quick"
const PFIDENTITY = "identity"
const PFMATLAB = "matlab"

# Skipping MacroDefinition: PFSetOptionsPrefix ( a , s ) PetscObjectSetOptionsPrefix ( ( PetscObject ) ( a ) , s )

const AOBASIC = "basic"
const AOADVANCED = "advanced"
const AOMAPPING = "mapping"
const AOMEMORYSCALABLE = "memoryscalable"
const PETSCSPACEPOLYNOMIAL = "poly"
const PETSCSPACEDG = "dg"
const PETSCDUALSPACELAGRANGE = "lagrange"
const PETSCDUALSPACESIMPLE = "simple"
const PETSCFEBASIC = "basic"
const PETSCFENONAFFINE = "nonaffine"
const PETSCFEOPENCL = "opencl"
const PETSCFECOMPOSITE = "composite"
const MATSEQUSFFT = "sequsfft"
const PETSCLIMITERSIN = "sin"
const PETSCLIMITERZERO = "zero"
const PETSCLIMITERNONE = "none"
const PETSCLIMITERMINMOD = "minmod"
const PETSCLIMITERVANLEER = "vanleer"
const PETSCLIMITERVANALBADA = "vanalbada"
const PETSCLIMITERSUPERBEE = "superbee"
const PETSCLIMITERMC = "mc"
const PETSCFVUPWIND = "upwind"
const PETSCFVLEASTSQUARES = "leastsquares"
const PETSCPARTITIONERCHACO = "chaco"
const PETSCPARTITIONERPARMETIS = "parmetis"
const PETSCPARTITIONERSHELL = "shell"
const PETSCPARTITIONERSIMPLE = "simple"
const PETSCDSBASIC = "basic"
const CHARACTERISTICDA = "da"
const PCNONE = "none"
const PCJACOBI = "jacobi"
const PCSOR = "sor"
const PCLU = "lu"
const PCSHELL = "shell"
const PCBJACOBI = "bjacobi"
const PCMG = "mg"
const PCEISENSTAT = "eisenstat"
const PCILU = "ilu"
const PCICC = "icc"
const PCASM = "asm"
const PCGASM = "gasm"
const PCKSP = "ksp"
const PCCOMPOSITE = "composite"
const PCREDUNDANT = "redundant"
const PCSPAI = "spai"
const PCNN = "nn"
const PCCHOLESKY = "cholesky"
const PCPBJACOBI = "pbjacobi"
const PCMAT = "mat"
const PCHYPRE = "hypre"
const PCPARMS = "parms"
const PCFIELDSPLIT = "fieldsplit"
const PCTFS = "tfs"
const PCML = "ml"
const PCGALERKIN = "galerkin"
const PCEXOTIC = "exotic"
const PCCP = "cp"
const PCBFBT = "bfbt"
const PCLSC = "lsc"
const PCPYTHON = "python"
const PCPFMG = "pfmg"
const PCSYSPFMG = "syspfmg"
const PCREDISTRIBUTE = "redistribute"
const PCSVD = "svd"
const PCGAMG = "gamg"
const PCSACUSP = "sacusp"
const PCSACUSPPOLY = "sacusppoly"
const PCBICGSTABCUSP = "bicgstabcusp"
const PCAINVCUSP = "ainvcusp"
const PCBDDC = "bddc"
const PCKACZMARZ = "kaczmarz"

# begin enum PCSide
typealias PCSide Cint
const PC_SIDE_DEFAULT = (Int32)(-1)
const PC_LEFT = (Int32)(0)
const PC_RIGHT = (Int32)(1)
const PC_SYMMETRIC = (Int32)(2)
# end enum PCSide

const PC_SIDE_MAX = PC_SYMMETRIC + 1
const PCGAMGAGG = "agg"
const PCGAMGGEO = "geo"
const PCGAMGCLASSICAL = "classical"
const PCGAMGCLASSICALDIRECT = "direct"
const PCGAMGCLASSICALSTANDARD = "standard"

# begin enum PCMGType
typealias PCMGType Uint32
const PC_MG_MULTIPLICATIVE = (UInt32)(0)
const PC_MG_ADDITIVE = (UInt32)(1)
const PC_MG_FULL = (UInt32)(2)
const PC_MG_KASKADE = (UInt32)(3)
# end enum PCMGType

const PC_MG_CASCADE = PC_MG_KASKADE
const PC_FILE_CLASSID = 1211222
const KSPRICHARDSON = "richardson"
const KSPCHEBYSHEV = "chebyshev"
const KSPCG = "cg"
const KSPGROPPCG = "groppcg"
const KSPPIPECG = "pipecg"
const KSPCGNE = "cgne"
const KSPNASH = "nash"
const KSPSTCG = "stcg"
const KSPGLTR = "gltr"
const KSPFCG = "fcg"
const KSPGMRES = "gmres"
const KSPFGMRES = "fgmres"
const KSPLGMRES = "lgmres"
const KSPDGMRES = "dgmres"
const KSPPGMRES = "pgmres"
const KSPTCQMR = "tcqmr"
const KSPBCGS = "bcgs"
const KSPIBCGS = "ibcgs"
const KSPFBCGS = "fbcgs"
const KSPFBCGSR = "fbcgsr"
const KSPBCGSL = "bcgsl"
const KSPCGS = "cgs"
const KSPTFQMR = "tfqmr"
const KSPCR = "cr"
const KSPPIPECR = "pipecr"
const KSPLSQR = "lsqr"
const KSPPREONLY = "preonly"
const KSPQCG = "qcg"
const KSPBICG = "bicg"
const KSPMINRES = "minres"
const KSPSYMMLQ = "symmlq"
const KSPLCD = "lcd"
const KSPPYTHON = "python"
const KSPGCR = "gcr"
const KSP_FILE_CLASSID = 1211223

# begin enum KSPNormType
typealias KSPNormType Cint
const KSP_NORM_DEFAULT = (Int32)(-1)
const KSP_NORM_NONE = (Int32)(0)
const KSP_NORM_PRECONDITIONED = (Int32)(1)
const KSP_NORM_UNPRECONDITIONED = (Int32)(2)
const KSP_NORM_NATURAL = (Int32)(3)
# end enum KSPNormType

const KSP_NORM_MAX = KSP_NORM_NATURAL + 1

# Skipping MacroDefinition: KSPDefaultConverged ( KSPDefaultConverged , KSPConvergedDefault )
# Skipping MacroDefinition: KSPDefaultConvergedDestroy ( KSPDefaultConvergedDestroy , KSPConvergedDefaultDestroy )
# Skipping MacroDefinition: KSPDefaultConvergedCreate ( KSPDefaultConvergedCreate , KSPConvergedDefaultCreate )
# Skipping MacroDefinition: KSPDefaultConvergedSetUIRNorm ( KSPDefaultConvergedSetUIRNorm , KSPConvergedDefaultSetUIRNorm )
# Skipping MacroDefinition: KSPDefaultConvergedSetUMIRNorm ( KSPDefaultConvergedSetUMIRNorm , KSPConvergedDefaultSetUMIRNorm )
# Skipping MacroDefinition: KSPSkipConverged ( KSPSkipConverged , KSPConvergedSkip )

const SNESNEWTONLS = "newtonls"
const SNESNEWTONTR = "newtontr"
const SNESPYTHON = "python"
const SNESTEST = "test"
const SNESNRICHARDSON = "nrichardson"
const SNESKSPONLY = "ksponly"
const SNESVINEWTONRSLS = "vinewtonrsls"
const SNESVINEWTONSSLS = "vinewtonssls"
const SNESNGMRES = "ngmres"
const SNESQN = "qn"
const SNESSHELL = "shell"
const SNESNGS = "ngs"
const SNESNCG = "ncg"
const SNESFAS = "fas"
const SNESMS = "ms"
const SNESNASM = "nasm"
const SNESANDERSON = "anderson"
const SNESASPIN = "aspin"
const SNESCOMPOSITE = "composite"
const SNES_FILE_CLASSID = 1211224

# Skipping MacroDefinition: SNESSkipConverged ( SNESSkipConverged , SNESConvergedSkip )

const SNESLINESEARCHBT = "bt"
const SNESLINESEARCHNLEQERR = "nleqerr"
const SNESLINESEARCHBASIC = "basic"
const SNESLINESEARCHL2 = "l2"
const SNESLINESEARCHCP = "cp"
const SNESLINESEARCHSHELL = "shell"
const SNES_LINESEARCH_ORDER_LINEAR = 1
const SNES_LINESEARCH_ORDER_QUADRATIC = 2
const SNES_LINESEARCH_ORDER_CUBIC = 3
const SNESMSM62 = "m62"
const SNESMSEULER = "euler"
const SNESMSJAMESON83 = "jameson83"
const SNESMSVLTP21 = "vltp21"
const SNESMSVLTP31 = "vltp31"
const SNESMSVLTP41 = "vltp41"
const SNESMSVLTP51 = "vltp51"
const SNESMSVLTP61 = "vltp61"
const TSEULER = "euler"
const TSBEULER = "beuler"
const TSPSEUDO = "pseudo"
const TSCN = "cn"
const TSSUNDIALS = "sundials"
const TSRK = "rk"
const TSPYTHON = "python"
const TSTHETA = "theta"
const TSALPHA = "alpha"
const TSGL = "gl"
const TSSSP = "ssp"
const TSARKIMEX = "arkimex"
const TSROSW = "rosw"
const TSEIMEX = "eimex"
const TSMIMEX = "mimex"
const TSTRAJECTORYBASIC = "basic"
const TSTRAJECTORYSINGLEFILE = "singlefile"
const TS_FILE_CLASSID = 1211225
const TSSSPRKS2 = "rks2"
const TSSSPRKS3 = "rks3"
const TSSSPRK104 = "rk104"
const TSADAPTBASIC = "basic"
const TSADAPTNONE = "none"
const TSADAPTCFL = "cfl"
const TSGLADAPT_NONE = "none"
const TSGLADAPT_SIZE = "size"
const TSGLADAPT_BOTH = "both"
const TSGLACCEPT_ALWAYS = "always"
const TSGL_IRKS = "irks"

# Skipping MacroDefinition: TSEIMEXType char *

const TSRK1FE = "1fe"
const TSRK2A = "2a"
const TSRK3 = "3"
const TSRK3BS = "3bs"
const TSRK4 = "4"
const TSRK5F = "5f"
const TSRK5DP = "5dp"
const TSARKIMEX1BEE = "1bee"
const TSARKIMEXA2 = "a2"
const TSARKIMEXL2 = "l2"
const TSARKIMEXARS122 = "ars122"
const TSARKIMEX2C = "2c"
const TSARKIMEX2D = "2d"
const TSARKIMEX2E = "2e"
const TSARKIMEXPRSSP2 = "prssp2"
const TSARKIMEX3 = "3"
const TSARKIMEXBPR3 = "bpr3"
const TSARKIMEXARS443 = "ars443"
const TSARKIMEX4 = "4"
const TSARKIMEX5 = "5"
const TSROSW2M = "2m"
const TSROSW2P = "2p"
const TSROSWRA3PW = "ra3pw"
const TSROSWRA34PW2 = "ra34pw2"
const TSROSWRODAS3 = "rodas3"
const TSROSWSANDU3 = "sandu3"
const TSROSWASSP3P3S1C = "assp3p3s1c"
const TSROSWLASSP3P4S2C = "lassp3p4s2c"
const TSROSWLLSSP3P4S2C = "llssp3p4s2c"
const TSROSWARK3 = "ark3"
const TSROSWTHETA1 = "theta1"
const TSROSWTHETA2 = "theta2"
const TSROSWGRK4T = "grk4t"
const TSROSWSHAMP4 = "shamp4"
const TSROSWVELDD4 = "veldd4"
const TSROSW4L = "4l"

# Skipping MacroDefinition: TaoType char *

const TAOLMVM = "lmvm"
const TAONLS = "nls"
const TAONTR = "ntr"
const TAONTL = "ntl"
const TAOCG = "cg"
const TAOTRON = "tron"
const TAOOWLQN = "owlqn"
const TAOBMRM = "bmrm"
const TAOBLMVM = "blmvm"
const TAOBQPIP = "bqpip"
const TAOGPCG = "gpcg"
const TAONM = "nm"
const TAOPOUNDERS = "pounders"
const TAOLCL = "lcl"
const TAOSSILS = "ssils"
const TAOSSFLS = "ssfls"
const TAOASILS = "asils"
const TAOASFLS = "asfls"
const TAOIPM = "ipm"
const TAOTEST = "test"

# Skipping MacroDefinition: TaoLineSearchType char *

const TAOLINESEARCHUNIT = "unit"
const TAOLINESEARCHMT = "more-thuente"
const TAOLINESEARCHGPCG = "gpcg"
const TAOLINESEARCHARMIJO = "armijo"
const TAOLINESEARCHOWARMIJO = "owarmijo"
const TAOLINESEARCHIPM = "ipm"

typealias PetscErrorCode Cint
typealias PetscClassId Cint
typealias PetscMPIInt Cint

# begin enum ANONYMOUS_1
typealias ANONYMOUS_1 Uint32
const ENUM_DUMMY = (UInt32)(0)
# end enum ANONYMOUS_1

# begin enum PetscEnum
typealias PetscEnum Uint32
const ENUM_DUMMY = (UInt32)(0)
# end enum PetscEnum

typealias Petsc64bitInt Cint
typealias PetscInt Cint
typealias PetscBLASInt Cint

# begin enum ANONYMOUS_2
typealias ANONYMOUS_2 Uint32
const PETSC_PRECISION_SINGLE = (UInt32)(4)
const PETSC_PRECISION_DOUBLE = (UInt32)(8)
# end enum ANONYMOUS_2

# begin enum PetscPrecision
typealias PetscPrecision Uint32
const PETSC_PRECISION_SINGLE = (UInt32)(4)
const PETSC_PRECISION_DOUBLE = (UInt32)(8)
# end enum PetscPrecision

# begin enum ANONYMOUS_3
typealias ANONYMOUS_3 Uint32
const PETSC_FALSE = (UInt32)(0)
const PETSC_TRUE = (UInt32)(1)
# end enum ANONYMOUS_3

# begin enum PetscBool
typealias PetscBool Uint32
const PETSC_FALSE = (UInt32)(0)
const PETSC_TRUE = (UInt32)(1)
# end enum PetscBool

typealias MatReal Cint

type petsc_mpiu_2scalar
    a::PetscScalar
    b::PetscScalar
end

typealias PetscLogDouble Cdouble

# begin enum ANONYMOUS_4
typealias ANONYMOUS_4 Uint32
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
# end enum ANONYMOUS_4

type _p_PetscToken
end

typealias PetscToken Ptr{_p_PetscToken}

type _p_PetscObject
end

typealias PetscObject Ptr{_p_PetscObject}
typealias PetscObjectId Petsc64bitInt
typealias PetscObjectState Petsc64bitInt

type _n_PetscFunctionList
end

typealias PetscFunctionList Ptr{_n_PetscFunctionList}

# begin enum ANONYMOUS_5
typealias ANONYMOUS_5 Uint32
const FILE_MODE_READ = (UInt32)(0)
const FILE_MODE_WRITE = (UInt32)(1)
const FILE_MODE_APPEND = (UInt32)(2)
const FILE_MODE_UPDATE = (UInt32)(3)
const FILE_MODE_APPEND_UPDATE = (UInt32)(4)
# end enum ANONYMOUS_5

# begin enum PetscFileMode
typealias PetscFileMode Uint32
const FILE_MODE_READ = (UInt32)(0)
const FILE_MODE_WRITE = (UInt32)(1)
const FILE_MODE_APPEND = (UInt32)(2)
const FILE_MODE_UPDATE = (UInt32)(3)
const FILE_MODE_APPEND_UPDATE = (UInt32)(4)
# end enum PetscFileMode

# begin enum ANONYMOUS_6
typealias ANONYMOUS_6 Uint32
const PETSC_ERROR_INITIAL = (UInt32)(0)
const PETSC_ERROR_REPEAT = (UInt32)(1)
const PETSC_ERROR_IN_CXX = (UInt32)(2)
# end enum ANONYMOUS_6

# begin enum PetscErrorType
typealias PetscErrorType Uint32
const PETSC_ERROR_INITIAL = (UInt32)(0)
const PETSC_ERROR_REPEAT = (UInt32)(1)
const PETSC_ERROR_IN_CXX = (UInt32)(2)
# end enum PetscErrorType

# begin enum ANONYMOUS_7
typealias ANONYMOUS_7 Uint32
const PETSC_FP_TRAP_OFF = (UInt32)(0)
const PETSC_FP_TRAP_ON = (UInt32)(1)
# end enum ANONYMOUS_7

# begin enum PetscFPTrap
typealias PetscFPTrap Uint32
const PETSC_FP_TRAP_OFF = (UInt32)(0)
const PETSC_FP_TRAP_ON = (UInt32)(1)
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

type _p_PetscViewer
end

typealias PetscViewer Ptr{_p_PetscViewer}

# begin enum ANONYMOUS_8
typealias ANONYMOUS_8 Uint32
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
# end enum ANONYMOUS_8

# begin enum PetscOptionType
typealias PetscOptionType Uint32
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
# end enum PetscOptionType

#typealias PetscOption Ptr{_n_PetscOption}

type _n_PetscOption
    option::Ptr{Uint8}
    text::Ptr{Uint8}
    data::Ptr{Void}
    flist::PetscFunctionList
    list::Ptr{Ptr{Uint8}}
    nlist::Uint8
    man::Ptr{Uint8}
    arraylength::Cint
    set::PetscBool
    _type::PetscOptionType
    next::Ptr{Void}  # was PetscOption
    pman::Ptr{Uint8}
    edata::Ptr{Void}
end


typealias PetscOption Ptr{_n_PetscOption}

typealias MPI_Comm MPI.Comm

type _p_PetscOptions
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
const PETSC_DL_DECIDE = (UInt32)(0)
const PETSC_DL_NOW = (UInt32)(1)
const PETSC_DL_LOCAL = (UInt32)(2)
# end enum ANONYMOUS_9

# begin enum PetscDLMode
typealias PetscDLMode Uint32
const PETSC_DL_DECIDE = (UInt32)(0)
const PETSC_DL_NOW = (UInt32)(1)
const PETSC_DL_LOCAL = (UInt32)(2)
# end enum PetscDLMode

type _n_PetscObjectList
end

typealias PetscObjectList Ptr{_n_PetscObjectList}

type _n_PetscDLLibrary
end

typealias PetscDLLibrary Ptr{_n_PetscDLLibrary}
typealias PetscLogEvent Cint
typealias PetscLogStage Cint

type _n_PetscIntStack
end

typealias PetscIntStack Ptr{_n_PetscIntStack}

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

type _n_PetscClassRegLog
    numClasses::Cint
    maxClasses::Cint
    classInfo::Ptr{PetscClassRegInfo}
end

typealias PetscClassRegLog Ptr{_n_PetscClassRegLog}

type _n_PetscClassPerfLog
    numClasses::Cint
    maxClasses::Cint
    classInfo::Ptr{PetscClassPerfInfo}
end

typealias PetscClassPerfLog Ptr{_n_PetscClassPerfLog}

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

type _n_PetscEventRegLog
    numEvents::Cint
    maxEvents::Cint
    eventInfo::Ptr{PetscEventRegInfo}
end

typealias PetscEventRegLog Ptr{_n_PetscEventRegLog}

type _n_PetscEventPerfLog
    numEvents::Cint
    maxEvents::Cint
    eventInfo::Ptr{PetscEventPerfInfo}
end

typealias PetscEventPerfLog Ptr{_n_PetscEventPerfLog}

type _PetscStageInfo
    name::Ptr{Uint8}
    used::PetscBool
    perfInfo::PetscEventPerfInfo
    eventLog::PetscEventPerfLog
    classLog::PetscClassPerfLog
end

type PetscStageInfo
    name::Ptr{Uint8}
    used::PetscBool
    perfInfo::PetscEventPerfInfo
    eventLog::PetscEventPerfLog
    classLog::PetscClassPerfLog
end

type _n_PetscStageLog
    numStages::Cint
    maxStages::Cint
    stack::PetscIntStack
    curStage::Cint
    stageInfo::Ptr{PetscStageInfo}
    eventLog::PetscEventRegLog
    classLog::PetscClassRegLog
end

typealias PetscStageLog Ptr{_n_PetscStageLog}

type _p_PetscContainer
end

typealias PetscContainer Ptr{_p_PetscContainer}
typealias PetscRandomType Ptr{Uint8}

type _p_PetscRandom
end

typealias PetscRandom Ptr{_p_PetscRandom}

# begin enum ANONYMOUS_10
typealias ANONYMOUS_10 Uint32
const PETSC_BINARY_SEEK_SET = (UInt32)(0)
const PETSC_BINARY_SEEK_CUR = (UInt32)(1)
const PETSC_BINARY_SEEK_END = (UInt32)(2)
# end enum ANONYMOUS_10

# begin enum PetscBinarySeekType
typealias PetscBinarySeekType Uint32
const PETSC_BINARY_SEEK_SET = (UInt32)(0)
const PETSC_BINARY_SEEK_CUR = (UInt32)(1)
const PETSC_BINARY_SEEK_END = (UInt32)(2)
# end enum PetscBinarySeekType

# begin enum ANONYMOUS_11
typealias ANONYMOUS_11 Cint
const PETSC_BUILDTWOSIDED_NOTSET = (Int32)(-1)
const PETSC_BUILDTWOSIDED_ALLREDUCE = (Int32)(0)
const PETSC_BUILDTWOSIDED_IBARRIER = (Int32)(1)
# end enum ANONYMOUS_11

# begin enum PetscBuildTwoSidedType
typealias PetscBuildTwoSidedType Cint
const PETSC_BUILDTWOSIDED_NOTSET = (Int32)(-1)
const PETSC_BUILDTWOSIDED_ALLREDUCE = (Int32)(0)
const PETSC_BUILDTWOSIDED_IBARRIER = (Int32)(1)
# end enum PetscBuildTwoSidedType

# begin enum ANONYMOUS_12
typealias ANONYMOUS_12 Uint32
const NOT_SET_VALUES = (UInt32)(0)
const INSERT_VALUES = (UInt32)(1)
const ADD_VALUES = (UInt32)(2)
const MAX_VALUES = (UInt32)(3)
const INSERT_ALL_VALUES = (UInt32)(4)
const ADD_ALL_VALUES = (UInt32)(5)
const INSERT_BC_VALUES = (UInt32)(6)
const ADD_BC_VALUES = (UInt32)(7)
# end enum ANONYMOUS_12

# begin enum InsertMode
typealias InsertMode Uint32
const NOT_SET_VALUES = (UInt32)(0)
const INSERT_VALUES = (UInt32)(1)
const ADD_VALUES = (UInt32)(2)
const MAX_VALUES = (UInt32)(3)
const INSERT_ALL_VALUES = (UInt32)(4)
const ADD_ALL_VALUES = (UInt32)(5)
const INSERT_BC_VALUES = (UInt32)(6)
const ADD_BC_VALUES = (UInt32)(7)
# end enum InsertMode

# begin enum ANONYMOUS_13
typealias ANONYMOUS_13 Uint32
const PETSC_SUBCOMM_GENERAL = (UInt32)(0)
const PETSC_SUBCOMM_CONTIGUOUS = (UInt32)(1)
const PETSC_SUBCOMM_INTERLACED = (UInt32)(2)
# end enum ANONYMOUS_13

# begin enum PetscSubcommType
typealias PetscSubcommType Uint32
const PETSC_SUBCOMM_GENERAL = (UInt32)(0)
const PETSC_SUBCOMM_CONTIGUOUS = (UInt32)(1)
const PETSC_SUBCOMM_INTERLACED = (UInt32)(2)
# end enum PetscSubcommType

type _n_PetscSubcomm
    parent::MPI_Comm
    dupparent::MPI_Comm
    child::MPI_Comm
    n::PetscMPIInt
    color::PetscMPIInt
    subsize::Ptr{PetscMPIInt}
    _type::PetscSubcommType
end

typealias PetscSubcomm Ptr{_n_PetscSubcomm}

type _n_PetscSegBuffer
end

typealias PetscSegBuffer Ptr{_n_PetscSegBuffer}

type _n_PetscBag
end

typealias PetscBag Ptr{_n_PetscBag}

type _n_PetscBagItem
end

typealias PetscBagItem Ptr{_n_PetscBagItem}
typealias PetscViewerType Ptr{Uint8}
typealias PetscDrawType Ptr{Uint8}

type _p_PetscDraw
end

typealias PetscDraw Ptr{_p_PetscDraw}

type _p_PetscDrawAxis
end

typealias PetscDrawAxis Ptr{_p_PetscDrawAxis}

type _p_PetscDrawLG
end

typealias PetscDrawLG Ptr{_p_PetscDrawLG}

type _p_PetscDrawSP
end

typealias PetscDrawSP Ptr{_p_PetscDrawSP}

type _p_PetscDrawHG
end

typealias PetscDrawHG Ptr{_p_PetscDrawHG}

type _p_PetscDrawBar
end

typealias PetscDrawBar Ptr{_p_PetscDrawBar}

# begin enum ANONYMOUS_14
typealias ANONYMOUS_14 Uint32
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
# end enum ANONYMOUS_14

# begin enum PetscViewerFormat
typealias PetscViewerFormat Uint32
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
# end enum PetscViewerFormat

# begin enum ANONYMOUS_15
typealias ANONYMOUS_15 Uint32
const PETSC_VTK_POINT_FIELD = (UInt32)(0)
const PETSC_VTK_POINT_VECTOR_FIELD = (UInt32)(1)
const PETSC_VTK_CELL_FIELD = (UInt32)(2)
const PETSC_VTK_CELL_VECTOR_FIELD = (UInt32)(3)
# end enum ANONYMOUS_15

# begin enum PetscViewerVTKFieldType
typealias PetscViewerVTKFieldType Uint32
const PETSC_VTK_POINT_FIELD = (UInt32)(0)
const PETSC_VTK_POINT_VECTOR_FIELD = (UInt32)(1)
const PETSC_VTK_CELL_FIELD = (UInt32)(2)
const PETSC_VTK_CELL_VECTOR_FIELD = (UInt32)(3)
# end enum PetscViewerVTKFieldType

type _n_PetscViewers
end

typealias PetscViewers Ptr{_n_PetscViewers}
typealias PetscBT Ptr{Uint8}

type _n_PetscTable
end

typealias PetscTable Ptr{_n_PetscTable}
typealias PetscTablePosition Ptr{PetscInt}

type _p_PetscMatlabEngine
end

typealias PetscMatlabEngine Ptr{_p_PetscMatlabEngine}

# begin enum ANONYMOUS_16
typealias ANONYMOUS_16 Uint32
const PETSC_DRAW_MARKER_CROSS = (UInt32)(0)
const PETSC_DRAW_MARKER_POINT = (UInt32)(1)
const PETSC_DRAW_MARKER_PLUS = (UInt32)(2)
const PETSC_DRAW_MARKER_CIRCLE = (UInt32)(3)
# end enum ANONYMOUS_16

# begin enum PetscDrawMarkerType
typealias PetscDrawMarkerType Uint32
const PETSC_DRAW_MARKER_CROSS = (UInt32)(0)
const PETSC_DRAW_MARKER_POINT = (UInt32)(1)
const PETSC_DRAW_MARKER_PLUS = (UInt32)(2)
const PETSC_DRAW_MARKER_CIRCLE = (UInt32)(3)
# end enum PetscDrawMarkerType

# begin enum ANONYMOUS_17
typealias ANONYMOUS_17 Uint32
const PETSC_BUTTON_NONE = (UInt32)(0)
const PETSC_BUTTON_LEFT = (UInt32)(1)
const PETSC_BUTTON_CENTER = (UInt32)(2)
const PETSC_BUTTON_RIGHT = (UInt32)(3)
const PETSC_BUTTON_LEFT_SHIFT = (UInt32)(4)
const PETSC_BUTTON_CENTER_SHIFT = (UInt32)(5)
const PETSC_BUTTON_RIGHT_SHIFT = (UInt32)(6)
# end enum ANONYMOUS_17

# begin enum PetscDrawButton
typealias PetscDrawButton Uint32
const PETSC_BUTTON_NONE = (UInt32)(0)
const PETSC_BUTTON_LEFT = (UInt32)(1)
const PETSC_BUTTON_CENTER = (UInt32)(2)
const PETSC_BUTTON_RIGHT = (UInt32)(3)
const PETSC_BUTTON_LEFT_SHIFT = (UInt32)(4)
const PETSC_BUTTON_CENTER_SHIFT = (UInt32)(5)
const PETSC_BUTTON_RIGHT_SHIFT = (UInt32)(6)
# end enum PetscDrawButton

type PetscDrawViewPorts
    nports::PetscInt
    xl::Ptr{Cint}
    xr::Ptr{Cint}
    yl::Ptr{Cint}
    yr::Ptr{Cint}
    draw::PetscDraw
    port_xl::Cint
    port_yl::Cint
    port_xr::Cint
    port_yr::Cint
end

type _p_PetscSF
end

typealias PetscSF Ptr{_p_PetscSF}

type PetscSFNode
    rank::PetscInt
    index::PetscInt
end

type _p_IS
end

typealias IS Ptr{_p_IS}

type _p_ISLocalToGlobalMapping
end

typealias ISLocalToGlobalMapping Ptr{_p_ISLocalToGlobalMapping}

type _n_ISColoring
end

typealias ISColoring Ptr{_n_ISColoring}

type _n_PetscLayout
    comm::MPI_Comm
    n::PetscInt
    N::PetscInt
    rstart::PetscInt
    rend::PetscInt
    range::Ptr{PetscInt}
    bs::PetscInt
    refcnt::PetscInt
    mapping::ISLocalToGlobalMapping
    trstarts::Ptr{PetscInt}
end

typealias PetscLayout Ptr{_n_PetscLayout}

type _p_PetscSection
end

typealias PetscSection Ptr{_p_PetscSection}
typealias ISType Ptr{Uint8}

# begin enum ANONYMOUS_18
typealias ANONYMOUS_18 Uint32
const IS_GTOLM_MASK = (UInt32)(0)
const IS_GTOLM_DROP = (UInt32)(1)
# end enum ANONYMOUS_18

# begin enum ISGlobalToLocalMappingType
typealias ISGlobalToLocalMappingType Uint32
const IS_GTOLM_MASK = (UInt32)(0)
const IS_GTOLM_DROP = (UInt32)(1)
# end enum ISGlobalToLocalMappingType

# begin enum ANONYMOUS_19
typealias ANONYMOUS_19 Uint32
const IS_COLORING_GLOBAL = (UInt32)(0)
const IS_COLORING_GHOSTED = (UInt32)(1)
# end enum ANONYMOUS_19

# begin enum ISColoringType
typealias ISColoringType Uint32
const IS_COLORING_GLOBAL = (UInt32)(0)
const IS_COLORING_GHOSTED = (UInt32)(1)
# end enum ISColoringType

typealias PETSC_IS_COLOR_VALUE_TYPE Uint32

type _p_Vec
end

typealias Vec Ptr{_p_Vec}

type _p_VecScatter
end

typealias VecScatter Ptr{_p_VecScatter}

# begin enum ANONYMOUS_20
typealias ANONYMOUS_20 Uint32
const SCATTER_FORWARD = (UInt32)(0)
const SCATTER_REVERSE = (UInt32)(1)
const SCATTER_FORWARD_LOCAL = (UInt32)(2)
const SCATTER_REVERSE_LOCAL = (UInt32)(3)
const SCATTER_LOCAL = (UInt32)(2)
# end enum ANONYMOUS_20

# begin enum ScatterMode
typealias ScatterMode Uint32
const SCATTER_FORWARD = (UInt32)(0)
const SCATTER_REVERSE = (UInt32)(1)
const SCATTER_FORWARD_LOCAL = (UInt32)(2)
const SCATTER_REVERSE_LOCAL = (UInt32)(3)
const SCATTER_LOCAL = (UInt32)(2)
# end enum ScatterMode

typealias VecType Ptr{Uint8}

# begin enum ANONYMOUS_21
typealias ANONYMOUS_21 Uint32
const NORM_1 = (UInt32)(0)
const NORM_2 = (UInt32)(1)
const NORM_FROBENIUS = (UInt32)(2)
const NORM_INFINITY = (UInt32)(3)
const NORM_1_AND_2 = (UInt32)(4)
# end enum ANONYMOUS_21

# begin enum ANONYMOUS_22
typealias ANONYMOUS_22 Uint32
const VEC_IGNORE_OFF_PROC_ENTRIES = (UInt32)(0)
const VEC_IGNORE_NEGATIVE_INDICES = (UInt32)(1)
# end enum ANONYMOUS_22

# begin enum VecOption
typealias VecOption Uint32
const VEC_IGNORE_OFF_PROC_ENTRIES = (UInt32)(0)
const VEC_IGNORE_NEGATIVE_INDICES = (UInt32)(1)
# end enum VecOption

# begin enum ANONYMOUS_23
typealias ANONYMOUS_23 Uint32
const VECOP_VIEW = (UInt32)(33)
const VECOP_LOAD = (UInt32)(41)
const VECOP_DUPLICATE = (UInt32)(0)
# end enum ANONYMOUS_23

# begin enum VecOperation
typealias VecOperation Uint32
const VECOP_VIEW = (UInt32)(33)
const VECOP_LOAD = (UInt32)(41)
const VECOP_DUPLICATE = (UInt32)(0)
# end enum VecOperation

type _n_Vecs
    n::PetscInt
    v::Vec
end

typealias Vecs Ptr{_n_Vecs}

type _p_Mat
end

typealias Mat Ptr{_p_Mat}
typealias MatType Ptr{Uint8}

# begin enum ANONYMOUS_24
typealias ANONYMOUS_24 Uint32
const MAT_FACTOR_NONE = (UInt32)(0)
const MAT_FACTOR_LU = (UInt32)(1)
const MAT_FACTOR_CHOLESKY = (UInt32)(2)
const MAT_FACTOR_ILU = (UInt32)(3)
const MAT_FACTOR_ICC = (UInt32)(4)
const MAT_FACTOR_ILUDT = (UInt32)(5)
# end enum ANONYMOUS_24

# begin enum MatFactorType
typealias MatFactorType Uint32
const MAT_FACTOR_NONE = (UInt32)(0)
const MAT_FACTOR_LU = (UInt32)(1)
const MAT_FACTOR_CHOLESKY = (UInt32)(2)
const MAT_FACTOR_ILU = (UInt32)(3)
const MAT_FACTOR_ICC = (UInt32)(4)
const MAT_FACTOR_ILUDT = (UInt32)(5)
# end enum MatFactorType

# begin enum ANONYMOUS_25
typealias ANONYMOUS_25 Uint32
const MAT_INITIAL_MATRIX = (UInt32)(0)
const MAT_REUSE_MATRIX = (UInt32)(1)
const MAT_IGNORE_MATRIX = (UInt32)(2)
# end enum ANONYMOUS_25

# begin enum MatReuse
typealias MatReuse Uint32
const MAT_INITIAL_MATRIX = (UInt32)(0)
const MAT_REUSE_MATRIX = (UInt32)(1)
const MAT_IGNORE_MATRIX = (UInt32)(2)
# end enum MatReuse

# begin enum ANONYMOUS_26
typealias ANONYMOUS_26 Uint32
const MAT_DO_NOT_GET_VALUES = (UInt32)(0)
const MAT_GET_VALUES = (UInt32)(1)
# end enum ANONYMOUS_26

# begin enum MatGetSubMatrixOption
typealias MatGetSubMatrixOption Uint32
const MAT_DO_NOT_GET_VALUES = (UInt32)(0)
const MAT_GET_VALUES = (UInt32)(1)
# end enum MatGetSubMatrixOption

# begin enum ANONYMOUS_27
typealias ANONYMOUS_27 Uint32
const DIFFERENT_NONZERO_PATTERN = (UInt32)(0)
const SUBSET_NONZERO_PATTERN = (UInt32)(1)
const SAME_NONZERO_PATTERN = (UInt32)(2)
# end enum ANONYMOUS_27

# begin enum MatStructure
typealias MatStructure Uint32
const DIFFERENT_NONZERO_PATTERN = (UInt32)(0)
const SUBSET_NONZERO_PATTERN = (UInt32)(1)
const SAME_NONZERO_PATTERN = (UInt32)(2)
# end enum MatStructure

# begin enum ANONYMOUS_28
typealias ANONYMOUS_28 Uint32
const MAT_COMPOSITE_ADDITIVE = (UInt32)(0)
const MAT_COMPOSITE_MULTIPLICATIVE = (UInt32)(1)
# end enum ANONYMOUS_28

# begin enum MatCompositeType
typealias MatCompositeType Uint32
const MAT_COMPOSITE_ADDITIVE = (UInt32)(0)
const MAT_COMPOSITE_MULTIPLICATIVE = (UInt32)(1)
# end enum MatCompositeType

type MatStencil
    k::PetscInt
    j::PetscInt
    i::PetscInt
    c::PetscInt
end

# begin enum ANONYMOUS_29
typealias ANONYMOUS_29 Uint32
const MAT_FLUSH_ASSEMBLY = (UInt32)(1)
const MAT_FINAL_ASSEMBLY = (UInt32)(0)
# end enum ANONYMOUS_29

# begin enum MatAssemblyType
typealias MatAssemblyType Uint32
const MAT_FLUSH_ASSEMBLY = (UInt32)(1)
const MAT_FINAL_ASSEMBLY = (UInt32)(0)
# end enum MatAssemblyType

# begin enum ANONYMOUS_30
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
# end enum ANONYMOUS_30

# begin enum MatOption
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
# end enum MatOption

# begin enum ANONYMOUS_31
typealias ANONYMOUS_31 Uint32
const MAT_DO_NOT_COPY_VALUES = (UInt32)(0)
const MAT_COPY_VALUES = (UInt32)(1)
const MAT_SHARE_NONZERO_PATTERN = (UInt32)(2)
# end enum ANONYMOUS_31

# begin enum MatDuplicateOption
typealias MatDuplicateOption Uint32
const MAT_DO_NOT_COPY_VALUES = (UInt32)(0)
const MAT_COPY_VALUES = (UInt32)(1)
const MAT_SHARE_NONZERO_PATTERN = (UInt32)(2)
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

# begin enum ANONYMOUS_32
typealias ANONYMOUS_32 Uint32
const MAT_LOCAL = (UInt32)(1)
const MAT_GLOBAL_MAX = (UInt32)(2)
const MAT_GLOBAL_SUM = (UInt32)(3)
# end enum ANONYMOUS_32

# begin enum MatInfoType
typealias MatInfoType Uint32
const MAT_LOCAL = (UInt32)(1)
const MAT_GLOBAL_MAX = (UInt32)(2)
const MAT_GLOBAL_SUM = (UInt32)(3)
# end enum MatInfoType

typealias MatOrderingType Ptr{Uint8}

# begin enum ANONYMOUS_33
typealias ANONYMOUS_33 Uint32
const MAT_SHIFT_NONE = (UInt32)(0)
const MAT_SHIFT_NONZERO = (UInt32)(1)
const MAT_SHIFT_POSITIVE_DEFINITE = (UInt32)(2)
const MAT_SHIFT_INBLOCKS = (UInt32)(3)
# end enum ANONYMOUS_33

# begin enum MatFactorShiftType
typealias MatFactorShiftType Uint32
const MAT_SHIFT_NONE = (UInt32)(0)
const MAT_SHIFT_NONZERO = (UInt32)(1)
const MAT_SHIFT_POSITIVE_DEFINITE = (UInt32)(2)
const MAT_SHIFT_INBLOCKS = (UInt32)(3)
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

# begin enum ANONYMOUS_34
typealias ANONYMOUS_34 Uint32
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
# end enum ANONYMOUS_34

# begin enum MatSORType
typealias MatSORType Uint32
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
# end enum MatSORType

type _p_MatColoring
end

typealias MatColoring Ptr{_p_MatColoring}
typealias MatColoringType Ptr{Uint8}

# begin enum ANONYMOUS_35
typealias ANONYMOUS_35 Uint32
const MAT_COLORING_WEIGHT_RANDOM = (UInt32)(0)
const MAT_COLORING_WEIGHT_LEXICAL = (UInt32)(1)
const MAT_COLORING_WEIGHT_LF = (UInt32)(2)
const MAT_COLORING_WEIGHT_SL = (UInt32)(3)
# end enum ANONYMOUS_35

# begin enum MatColoringWeightType
typealias MatColoringWeightType Uint32
const MAT_COLORING_WEIGHT_RANDOM = (UInt32)(0)
const MAT_COLORING_WEIGHT_LEXICAL = (UInt32)(1)
const MAT_COLORING_WEIGHT_LF = (UInt32)(2)
const MAT_COLORING_WEIGHT_SL = (UInt32)(3)
# end enum MatColoringWeightType

type _p_MatFDColoring
end

typealias MatFDColoring Ptr{_p_MatFDColoring}

type _p_MatTransposeColoring
end

typealias MatTransposeColoring Ptr{_p_MatTransposeColoring}

type _p_MatPartitioning
end

typealias MatPartitioning Ptr{_p_MatPartitioning}
typealias MatPartitioningType Ptr{Uint8}

# begin enum ANONYMOUS_36
typealias ANONYMOUS_36 Uint32
const MP_CHACO_MULTILEVEL = (UInt32)(1)
const MP_CHACO_SPECTRAL = (UInt32)(2)
const MP_CHACO_LINEAR = (UInt32)(4)
const MP_CHACO_RANDOM = (UInt32)(5)
const MP_CHACO_SCATTERED = (UInt32)(6)
# end enum ANONYMOUS_36

# begin enum MPChacoGlobalType
typealias MPChacoGlobalType Uint32
const MP_CHACO_MULTILEVEL = (UInt32)(1)
const MP_CHACO_SPECTRAL = (UInt32)(2)
const MP_CHACO_LINEAR = (UInt32)(4)
const MP_CHACO_RANDOM = (UInt32)(5)
const MP_CHACO_SCATTERED = (UInt32)(6)
# end enum MPChacoGlobalType

# begin enum ANONYMOUS_37
typealias ANONYMOUS_37 Uint32
const MP_CHACO_KERNIGHAN = (UInt32)(1)
const MP_CHACO_NONE = (UInt32)(2)
# end enum ANONYMOUS_37

# begin enum MPChacoLocalType
typealias MPChacoLocalType Uint32
const MP_CHACO_KERNIGHAN = (UInt32)(1)
const MP_CHACO_NONE = (UInt32)(2)
# end enum MPChacoLocalType

# begin enum ANONYMOUS_38
typealias ANONYMOUS_38 Uint32
const MP_CHACO_LANCZOS = (UInt32)(0)
const MP_CHACO_RQI = (UInt32)(1)
# end enum ANONYMOUS_38

# begin enum MPChacoEigenType
typealias MPChacoEigenType Uint32
const MP_CHACO_LANCZOS = (UInt32)(0)
const MP_CHACO_RQI = (UInt32)(1)
# end enum MPChacoEigenType

# begin enum ANONYMOUS_39
typealias ANONYMOUS_39 Uint32
const MP_PTSCOTCH_QUALITY = (UInt32)(0)
const MP_PTSCOTCH_SPEED = (UInt32)(1)
const MP_PTSCOTCH_BALANCE = (UInt32)(2)
const MP_PTSCOTCH_SAFETY = (UInt32)(3)
const MP_PTSCOTCH_SCALABILITY = (UInt32)(4)
# end enum ANONYMOUS_39

# begin enum MPPTScotchStrategyType
typealias MPPTScotchStrategyType Uint32
const MP_PTSCOTCH_QUALITY = (UInt32)(0)
const MP_PTSCOTCH_SPEED = (UInt32)(1)
const MP_PTSCOTCH_BALANCE = (UInt32)(2)
const MP_PTSCOTCH_SAFETY = (UInt32)(3)
const MP_PTSCOTCH_SCALABILITY = (UInt32)(4)
# end enum MPPTScotchStrategyType

type _p_MatCoarsen
end

typealias MatCoarsen Ptr{_p_MatCoarsen}
typealias MatCoarsenType Ptr{Uint8}

type _PetscCDIntNd
    next::Ptr{_PetscCDIntNd}
    gid::PetscInt
end

type PetscCDIntNd
    next::Ptr{_PetscCDIntNd}
    gid::PetscInt
end

type _PetscCDArrNd
    next::Ptr{_PetscCDArrNd}
    array::Ptr{_PetscCDIntNd}
end

type PetscCDArrNd
    next::Ptr{_PetscCDArrNd}
    array::Ptr{_PetscCDIntNd}
end

type _PetscCoarsenData
    pool_list::PetscCDArrNd
    new_node::Ptr{PetscCDIntNd}
    new_left::PetscInt
    chk_sz::PetscInt
    extra_nodes::Ptr{PetscCDIntNd}
    array::Ptr{Ptr{PetscCDIntNd}}
    size::PetscInt
    mat::Mat
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

# begin enum ANONYMOUS_40
typealias ANONYMOUS_40 Uint32
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
# end enum ANONYMOUS_40

# begin enum MatOperation
typealias MatOperation Uint32
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
# end enum MatOperation

type _p_MatNullSpace
end

typealias MatNullSpace Ptr{_p_MatNullSpace}

type _p_MatMFFD
end

typealias MatMFFD Ptr{_p_MatMFFD}
typealias MatMFFDType Ptr{Uint8}

type _p_DM
end

typealias DM Ptr{_p_DM}

# begin enum ANONYMOUS_41
typealias ANONYMOUS_41 Uint32
const DM_BOUNDARY_NONE = (UInt32)(0)
const DM_BOUNDARY_GHOSTED = (UInt32)(1)
const DM_BOUNDARY_MIRROR = (UInt32)(2)
const DM_BOUNDARY_PERIODIC = (UInt32)(3)
const DM_BOUNDARY_TWIST = (UInt32)(4)
# end enum ANONYMOUS_41

# begin enum DMBoundaryType
typealias DMBoundaryType Uint32
const DM_BOUNDARY_NONE = (UInt32)(0)
const DM_BOUNDARY_GHOSTED = (UInt32)(1)
const DM_BOUNDARY_MIRROR = (UInt32)(2)
const DM_BOUNDARY_PERIODIC = (UInt32)(3)
const DM_BOUNDARY_TWIST = (UInt32)(4)
# end enum DMBoundaryType

type _p_PetscPartitioner
end

typealias PetscPartitioner Ptr{_p_PetscPartitioner}

type _p_PetscSpace
end

typealias PetscSpace Ptr{_p_PetscSpace}

type _p_PetscDualSpace
end

typealias PetscDualSpace Ptr{_p_PetscDualSpace}

type _p_PetscFE
end

typealias PetscFE Ptr{_p_PetscFE}

type _p_PetscDS
end

typealias PetscDS Ptr{_p_PetscDS}
typealias DMType Ptr{Uint8}

type NLF_DAAD
end

typealias NLF Ptr{NLF_DAAD}

# begin enum ANONYMOUS_42
typealias ANONYMOUS_42 Uint32
const PETSC_UNIT_LENGTH = (UInt32)(0)
const PETSC_UNIT_MASS = (UInt32)(1)
const PETSC_UNIT_TIME = (UInt32)(2)
const PETSC_UNIT_CURRENT = (UInt32)(3)
const PETSC_UNIT_TEMPERATURE = (UInt32)(4)
const PETSC_UNIT_AMOUNT = (UInt32)(5)
const PETSC_UNIT_LUMINOSITY = (UInt32)(6)
const NUM_PETSC_UNITS = (UInt32)(7)
# end enum ANONYMOUS_42

# begin enum PetscUnit
typealias PetscUnit Uint32
const PETSC_UNIT_LENGTH = (UInt32)(0)
const PETSC_UNIT_MASS = (UInt32)(1)
const PETSC_UNIT_TIME = (UInt32)(2)
const PETSC_UNIT_CURRENT = (UInt32)(3)
const PETSC_UNIT_TEMPERATURE = (UInt32)(4)
const PETSC_UNIT_AMOUNT = (UInt32)(5)
const PETSC_UNIT_LUMINOSITY = (UInt32)(6)
const NUM_PETSC_UNITS = (UInt32)(7)
# end enum PetscUnit

type _DMInterpolationInfo
    comm::MPI_Comm
    dim::PetscInt
    nInput::PetscInt
    points::Ptr{Cint}
    cells::Ptr{PetscInt}
    n::PetscInt
    coords::Vec
    dof::PetscInt
end

typealias DMInterpolationInfo Ptr{_DMInterpolationInfo}

# begin enum ANONYMOUS_43
typealias ANONYMOUS_43 Uint32
const DMDA_STENCIL_STAR = (UInt32)(0)
const DMDA_STENCIL_BOX = (UInt32)(1)
# end enum ANONYMOUS_43

# begin enum DMDAStencilType
typealias DMDAStencilType Uint32
const DMDA_STENCIL_STAR = (UInt32)(0)
const DMDA_STENCIL_BOX = (UInt32)(1)
# end enum DMDAStencilType

# begin enum ANONYMOUS_44
typealias ANONYMOUS_44 Uint32
const DMDA_Q0 = (UInt32)(0)
const DMDA_Q1 = (UInt32)(1)
# end enum ANONYMOUS_44

# begin enum DMDAInterpolationType
typealias DMDAInterpolationType Uint32
const DMDA_Q0 = (UInt32)(0)
const DMDA_Q1 = (UInt32)(1)
# end enum DMDAInterpolationType

# begin enum ANONYMOUS_45
typealias ANONYMOUS_45 Uint32
const DMDA_ELEMENT_P1 = (UInt32)(0)
const DMDA_ELEMENT_Q1 = (UInt32)(1)
# end enum ANONYMOUS_45

# begin enum DMDAElementType
typealias DMDAElementType Uint32
const DMDA_ELEMENT_P1 = (UInt32)(0)
const DMDA_ELEMENT_Q1 = (UInt32)(1)
# end enum DMDAElementType

type DMDALocalInfo
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

typealias PFType Ptr{Uint8}

type _p_PF
end

typealias PF Ptr{_p_PF}
typealias AOType Ptr{Uint8}

type _p_PetscQuadrature
end

typealias PetscQuadrature Ptr{_p_PetscQuadrature}

immutable Array_3_Cint
    d1::Cint
    d2::Cint
    d3::Cint
end

zero(::Type{Array_3_Cint}) = begin  # /home/jared/.julia/v0.4/Clang/src/wrap_c.jl, line 264:
        Array_3_Cint(fill(zero(Cint),3)...)
    end

immutable Array_9_Cint
    d1::Cint
    d2::Cint
    d3::Cint
    d4::Cint
    d5::Cint
    d6::Cint
    d7::Cint
    d8::Cint
    d9::Cint
end

zero(::Type{Array_9_Cint}) = begin  # /home/jared/.julia/v0.4/Clang/src/wrap_c.jl, line 264:
        Array_9_Cint(fill(zero(Cint),9)...)
    end

type PetscFECellGeom
    v0::Array_3_Cint
    J::Array_9_Cint
    invJ::Array_9_Cint
    detJ::Cint
    n::Array_3_Cint
    dim::PetscInt
    dimEmbed::PetscInt
end

typealias PetscSpaceType Ptr{Uint8}
typealias PetscDualSpaceType Ptr{Uint8}
typealias PetscFEType Ptr{Uint8}

# begin enum ANONYMOUS_46
typealias ANONYMOUS_46 Uint32
const DMDA_X = (UInt32)(0)
const DMDA_Y = (UInt32)(1)
const DMDA_Z = (UInt32)(2)
# end enum ANONYMOUS_46

# begin enum DMDADirection
typealias DMDADirection Uint32
const DMDA_X = (UInt32)(0)
const DMDA_Y = (UInt32)(1)
const DMDA_Z = (UInt32)(2)
# end enum DMDADirection

type DMDACoor2d
    x::PetscScalar
    y::PetscScalar
end

type DMDACoor3d
    x::PetscScalar
    y::PetscScalar
    z::PetscScalar
end

type _p_PetscLimiter
end

typealias PetscLimiter Ptr{_p_PetscLimiter}

type _p_PetscFV
end

typealias PetscFV Ptr{_p_PetscFV}

immutable Array_3_PetscScalar
    d1::PetscScalar
    d2::PetscScalar
    d3::PetscScalar
end

zero(::Type{Array_3_PetscScalar}) = begin  # /home/jared/.julia/v0.4/Clang/src/wrap_c.jl, line 264:
        Array_3_PetscScalar(fill(zero(PetscScalar),3)...)
    end

immutable Array_2_Array_3_PetscScalar
    d1::Array_3_PetscScalar
    d2::Array_3_PetscScalar
end

zero(::Type{Array_2_Array_3_PetscScalar}) = begin  # /home/jared/.julia/v0.4/Clang/src/wrap_c.jl, line 264:
        Array_2_Array_3_PetscScalar(fill(zero(Array_3_PetscScalar),2)...)
    end

type PetscFVFaceGeom
    normal::Array_3_Cint
    centroid::Array_3_Cint
    grad::Array_2_Array_3_PetscScalar
end

type PetscFVCellGeom
    centroid::Array_3_Cint
    volume::Cint
end

typealias PetscLimiterType Ptr{Uint8}
typealias PetscFVType Ptr{Uint8}
typealias PetscPartitionerType Ptr{Uint8}

type _n_DMLabel
end

typealias DMLabel Ptr{_n_DMLabel}

type _n_Boundary
end

typealias DMBoundary Ptr{_n_Boundary}

type JacActionCtx
    dm::DM
    u::Vec
    J::Mat
    user::Ptr{Void}
end

typealias PetscDSType Ptr{Uint8}
typealias PetscPointFunc Ptr{Void}
typealias PetscPointJac Ptr{Void}
typealias PetscBdPointFunc Ptr{Void}
typealias PetscBdPointJac Ptr{Void}
typealias PetscRiemannFunc Ptr{Void}

type _p_Characteristic
end

typealias Characteristic Ptr{_p_Characteristic}
typealias CharacteristicType Ptr{Uint8}

type _p_PC
end

typealias PC Ptr{_p_PC}
typealias PCType Ptr{Uint8}

# begin enum ANONYMOUS_47
typealias ANONYMOUS_47 Cint
const PC_SIDE_DEFAULT = (Int32)(-1)
const PC_LEFT = (Int32)(0)
const PC_RIGHT = (Int32)(1)
const PC_SYMMETRIC = (Int32)(2)
# end enum ANONYMOUS_47

# begin enum ANONYMOUS_48
typealias ANONYMOUS_48 Cint
const PCRICHARDSON_CONVERGED_RTOL = (Int32)(2)
const PCRICHARDSON_CONVERGED_ATOL = (Int32)(3)
const PCRICHARDSON_CONVERGED_ITS = (Int32)(4)
const PCRICHARDSON_DIVERGED_DTOL = (Int32)(-4)
# end enum ANONYMOUS_48

# begin enum PCRichardsonConvergedReason
typealias PCRichardsonConvergedReason Cint
const PCRICHARDSON_CONVERGED_RTOL = (Int32)(2)
const PCRICHARDSON_CONVERGED_ATOL = (Int32)(3)
const PCRICHARDSON_CONVERGED_ITS = (Int32)(4)
const PCRICHARDSON_DIVERGED_DTOL = (Int32)(-4)
# end enum PCRichardsonConvergedReason

# begin enum ANONYMOUS_49
typealias ANONYMOUS_49 Uint32
const PC_JACOBI_DIAGONAL = (UInt32)(0)
const PC_JACOBI_ROWMAX = (UInt32)(1)
const PC_JACOBI_ROWSUM = (UInt32)(2)
# end enum ANONYMOUS_49

# begin enum PCJacobiType
typealias PCJacobiType Uint32
const PC_JACOBI_DIAGONAL = (UInt32)(0)
const PC_JACOBI_ROWMAX = (UInt32)(1)
const PC_JACOBI_ROWSUM = (UInt32)(2)
# end enum PCJacobiType

# begin enum ANONYMOUS_50
typealias ANONYMOUS_50 Uint32
const PC_ASM_BASIC = (UInt32)(3)
const PC_ASM_RESTRICT = (UInt32)(1)
const PC_ASM_INTERPOLATE = (UInt32)(2)
const PC_ASM_NONE = (UInt32)(0)
# end enum ANONYMOUS_50

# begin enum PCASMType
typealias PCASMType Uint32
const PC_ASM_BASIC = (UInt32)(3)
const PC_ASM_RESTRICT = (UInt32)(1)
const PC_ASM_INTERPOLATE = (UInt32)(2)
const PC_ASM_NONE = (UInt32)(0)
# end enum PCASMType

# begin enum ANONYMOUS_51
typealias ANONYMOUS_51 Uint32
const PC_GASM_BASIC = (UInt32)(3)
const PC_GASM_RESTRICT = (UInt32)(1)
const PC_GASM_INTERPOLATE = (UInt32)(2)
const PC_GASM_NONE = (UInt32)(0)
# end enum ANONYMOUS_51

# begin enum PCGASMType
typealias PCGASMType Uint32
const PC_GASM_BASIC = (UInt32)(3)
const PC_GASM_RESTRICT = (UInt32)(1)
const PC_GASM_INTERPOLATE = (UInt32)(2)
const PC_GASM_NONE = (UInt32)(0)
# end enum PCGASMType

# begin enum ANONYMOUS_52
typealias ANONYMOUS_52 Uint32
const PC_COMPOSITE_ADDITIVE = (UInt32)(0)
const PC_COMPOSITE_MULTIPLICATIVE = (UInt32)(1)
const PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE = (UInt32)(2)
const PC_COMPOSITE_SPECIAL = (UInt32)(3)
const PC_COMPOSITE_SCHUR = (UInt32)(4)
# end enum ANONYMOUS_52

# begin enum PCCompositeType
typealias PCCompositeType Uint32
const PC_COMPOSITE_ADDITIVE = (UInt32)(0)
const PC_COMPOSITE_MULTIPLICATIVE = (UInt32)(1)
const PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE = (UInt32)(2)
const PC_COMPOSITE_SPECIAL = (UInt32)(3)
const PC_COMPOSITE_SCHUR = (UInt32)(4)
# end enum PCCompositeType

# begin enum ANONYMOUS_53
typealias ANONYMOUS_53 Uint32
const PC_FIELDSPLIT_SCHUR_PRE_SELF = (UInt32)(0)
const PC_FIELDSPLIT_SCHUR_PRE_SELFP = (UInt32)(1)
const PC_FIELDSPLIT_SCHUR_PRE_A11 = (UInt32)(2)
const PC_FIELDSPLIT_SCHUR_PRE_USER = (UInt32)(3)
const PC_FIELDSPLIT_SCHUR_PRE_FULL = (UInt32)(4)
# end enum ANONYMOUS_53

# begin enum PCFieldSplitSchurPreType
typealias PCFieldSplitSchurPreType Uint32
const PC_FIELDSPLIT_SCHUR_PRE_SELF = (UInt32)(0)
const PC_FIELDSPLIT_SCHUR_PRE_SELFP = (UInt32)(1)
const PC_FIELDSPLIT_SCHUR_PRE_A11 = (UInt32)(2)
const PC_FIELDSPLIT_SCHUR_PRE_USER = (UInt32)(3)
const PC_FIELDSPLIT_SCHUR_PRE_FULL = (UInt32)(4)
# end enum PCFieldSplitSchurPreType

# begin enum ANONYMOUS_54
typealias ANONYMOUS_54 Uint32
const PC_FIELDSPLIT_SCHUR_FACT_DIAG = (UInt32)(0)
const PC_FIELDSPLIT_SCHUR_FACT_LOWER = (UInt32)(1)
const PC_FIELDSPLIT_SCHUR_FACT_UPPER = (UInt32)(2)
const PC_FIELDSPLIT_SCHUR_FACT_FULL = (UInt32)(3)
# end enum ANONYMOUS_54

# begin enum PCFieldSplitSchurFactType
typealias PCFieldSplitSchurFactType Uint32
const PC_FIELDSPLIT_SCHUR_FACT_DIAG = (UInt32)(0)
const PC_FIELDSPLIT_SCHUR_FACT_LOWER = (UInt32)(1)
const PC_FIELDSPLIT_SCHUR_FACT_UPPER = (UInt32)(2)
const PC_FIELDSPLIT_SCHUR_FACT_FULL = (UInt32)(3)
# end enum PCFieldSplitSchurFactType

# begin enum ANONYMOUS_55
typealias ANONYMOUS_55 Uint32
const PC_PARMS_GLOBAL_RAS = (UInt32)(0)
const PC_PARMS_GLOBAL_SCHUR = (UInt32)(1)
const PC_PARMS_GLOBAL_BJ = (UInt32)(2)
# end enum ANONYMOUS_55

# begin enum PCPARMSGlobalType
typealias PCPARMSGlobalType Uint32
const PC_PARMS_GLOBAL_RAS = (UInt32)(0)
const PC_PARMS_GLOBAL_SCHUR = (UInt32)(1)
const PC_PARMS_GLOBAL_BJ = (UInt32)(2)
# end enum PCPARMSGlobalType

# begin enum ANONYMOUS_56
typealias ANONYMOUS_56 Uint32
const PC_PARMS_LOCAL_ILU0 = (UInt32)(0)
const PC_PARMS_LOCAL_ILUK = (UInt32)(1)
const PC_PARMS_LOCAL_ILUT = (UInt32)(2)
const PC_PARMS_LOCAL_ARMS = (UInt32)(3)
# end enum ANONYMOUS_56

# begin enum PCPARMSLocalType
typealias PCPARMSLocalType Uint32
const PC_PARMS_LOCAL_ILU0 = (UInt32)(0)
const PC_PARMS_LOCAL_ILUK = (UInt32)(1)
const PC_PARMS_LOCAL_ILUT = (UInt32)(2)
const PC_PARMS_LOCAL_ARMS = (UInt32)(3)
# end enum PCPARMSLocalType

typealias PCGAMGType Ptr{Uint8}
typealias PCGAMGClassicalType Ptr{Uint8}

# begin enum ANONYMOUS_57
typealias ANONYMOUS_57 Uint32
const PC_MG_MULTIPLICATIVE = (UInt32)(0)
const PC_MG_ADDITIVE = (UInt32)(1)
const PC_MG_FULL = (UInt32)(2)
const PC_MG_KASKADE = (UInt32)(3)
# end enum ANONYMOUS_57

# begin enum ANONYMOUS_58
typealias ANONYMOUS_58 Uint32
const PC_MG_CYCLE_V = (UInt32)(1)
const PC_MG_CYCLE_W = (UInt32)(2)
# end enum ANONYMOUS_58

# begin enum PCMGCycleType
typealias PCMGCycleType Uint32
const PC_MG_CYCLE_V = (UInt32)(1)
const PC_MG_CYCLE_W = (UInt32)(2)
# end enum PCMGCycleType

# begin enum ANONYMOUS_59
typealias ANONYMOUS_59 Uint32
const PC_EXOTIC_FACE = (UInt32)(0)
const PC_EXOTIC_WIREBASKET = (UInt32)(1)
# end enum ANONYMOUS_59

# begin enum PCExoticType
typealias PCExoticType Uint32
const PC_EXOTIC_FACE = (UInt32)(0)
const PC_EXOTIC_WIREBASKET = (UInt32)(1)
# end enum PCExoticType

type _p_KSP
end

typealias KSP Ptr{_p_KSP}
typealias KSPType Ptr{Uint8}

# begin enum ANONYMOUS_60
typealias ANONYMOUS_60 Uint32
const KSP_FCG_TRUNC_TYPE_STANDARD = (UInt32)(0)
const KSP_FCG_TRUNC_TYPE_NOTAY = (UInt32)(1)
# end enum ANONYMOUS_60

# begin enum KSPFCGTruncationType
typealias KSPFCGTruncationType Uint32
const KSP_FCG_TRUNC_TYPE_STANDARD = (UInt32)(0)
const KSP_FCG_TRUNC_TYPE_NOTAY = (UInt32)(1)
# end enum KSPFCGTruncationType

# begin enum ANONYMOUS_61
typealias ANONYMOUS_61 Uint32
const KSP_GMRES_CGS_REFINE_NEVER = (UInt32)(0)
const KSP_GMRES_CGS_REFINE_IFNEEDED = (UInt32)(1)
const KSP_GMRES_CGS_REFINE_ALWAYS = (UInt32)(2)
# end enum ANONYMOUS_61

# begin enum KSPGMRESCGSRefinementType
typealias KSPGMRESCGSRefinementType Uint32
const KSP_GMRES_CGS_REFINE_NEVER = (UInt32)(0)
const KSP_GMRES_CGS_REFINE_IFNEEDED = (UInt32)(1)
const KSP_GMRES_CGS_REFINE_ALWAYS = (UInt32)(2)
# end enum KSPGMRESCGSRefinementType

# begin enum ANONYMOUS_62
typealias ANONYMOUS_62 Cint
const KSP_NORM_DEFAULT = (Int32)(-1)
const KSP_NORM_NONE = (Int32)(0)
const KSP_NORM_PRECONDITIONED = (Int32)(1)
const KSP_NORM_UNPRECONDITIONED = (Int32)(2)
const KSP_NORM_NATURAL = (Int32)(3)
# end enum ANONYMOUS_62

# begin enum ANONYMOUS_63
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
# end enum ANONYMOUS_63

# begin enum KSPConvergedReason
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
# end enum KSPConvergedReason

# begin enum ANONYMOUS_64
typealias ANONYMOUS_64 Uint32
const KSP_CG_SYMMETRIC = (UInt32)(0)
const KSP_CG_HERMITIAN = (UInt32)(1)
# end enum ANONYMOUS_64

# begin enum KSPCGType
typealias KSPCGType Uint32
const KSP_CG_SYMMETRIC = (UInt32)(0)
const KSP_CG_HERMITIAN = (UInt32)(1)
# end enum KSPCGType

type _p_KSPFischerGuess
    method::PetscInt
    curl::PetscInt
    maxl::PetscInt
    refcnt::PetscInt
    monitor::PetscBool
    mat::Mat
    ksp::KSP
end

typealias KSPFischerGuess Ptr{_p_KSPFischerGuess}

# begin enum ANONYMOUS_65
typealias ANONYMOUS_65 Uint32
const MAT_SCHUR_COMPLEMENT_AINV_DIAG = (UInt32)(0)
const MAT_SCHUR_COMPLEMENT_AINV_LUMP = (UInt32)(1)
# end enum ANONYMOUS_65

# begin enum MatSchurComplementAinvType
typealias MatSchurComplementAinvType Uint32
const MAT_SCHUR_COMPLEMENT_AINV_DIAG = (UInt32)(0)
const MAT_SCHUR_COMPLEMENT_AINV_LUMP = (UInt32)(1)
# end enum MatSchurComplementAinvType

type _p_SNES
end

typealias SNES Ptr{_p_SNES}
typealias SNESType Ptr{Uint8}

# begin enum ANONYMOUS_66
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
# end enum ANONYMOUS_66

# begin enum SNESConvergedReason
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
# end enum SNESConvergedReason

# begin enum ANONYMOUS_67
typealias ANONYMOUS_67 Cint
const SNES_NORM_DEFAULT = (Int32)(-1)
const SNES_NORM_NONE = (Int32)(0)
const SNES_NORM_ALWAYS = (Int32)(1)
const SNES_NORM_INITIAL_ONLY = (Int32)(2)
const SNES_NORM_FINAL_ONLY = (Int32)(3)
const SNES_NORM_INITIAL_FINAL_ONLY = (Int32)(4)
# end enum ANONYMOUS_67

# begin enum SNESNormSchedule
typealias SNESNormSchedule Cint
const SNES_NORM_DEFAULT = (Int32)(-1)
const SNES_NORM_NONE = (Int32)(0)
const SNES_NORM_ALWAYS = (Int32)(1)
const SNES_NORM_INITIAL_ONLY = (Int32)(2)
const SNES_NORM_FINAL_ONLY = (Int32)(3)
const SNES_NORM_INITIAL_FINAL_ONLY = (Int32)(4)
# end enum SNESNormSchedule

# begin enum ANONYMOUS_68
typealias ANONYMOUS_68 Cint
const SNES_FUNCTION_DEFAULT = (Int32)(-1)
const SNES_FUNCTION_UNPRECONDITIONED = (Int32)(0)
const SNES_FUNCTION_PRECONDITIONED = (Int32)(1)
# end enum ANONYMOUS_68

# begin enum SNESFunctionType
typealias SNESFunctionType Cint
const SNES_FUNCTION_DEFAULT = (Int32)(-1)
const SNES_FUNCTION_UNPRECONDITIONED = (Int32)(0)
const SNES_FUNCTION_PRECONDITIONED = (Int32)(1)
# end enum SNESFunctionType

type _p_LineSearch
end

typealias SNESLineSearch Ptr{_p_LineSearch}
typealias SNESLineSearchType Ptr{Uint8}
typealias SNESLineSearchVIProjectFunc Ptr{Void}
typealias SNESLineSearchVINormFunc Ptr{Void}
typealias SNESLineSearchApplyFunc Ptr{Void}
typealias SNESLineSearchUserFunc Ptr{Void}

# begin enum ANONYMOUS_69
typealias ANONYMOUS_69 Uint32
const SNES_LINESEARCH_SUCCEEDED = (UInt32)(0)
const SNES_LINESEARCH_FAILED_NANORINF = (UInt32)(1)
const SNES_LINESEARCH_FAILED_DOMAIN = (UInt32)(2)
const SNES_LINESEARCH_FAILED_REDUCT = (UInt32)(3)
const SNES_LINESEARCH_FAILED_USER = (UInt32)(4)
const SNES_LINESEARCH_FAILED_FUNCTION = (UInt32)(5)
# end enum ANONYMOUS_69

# begin enum SNESLineSearchReason
typealias SNESLineSearchReason Uint32
const SNES_LINESEARCH_SUCCEEDED = (UInt32)(0)
const SNES_LINESEARCH_FAILED_NANORINF = (UInt32)(1)
const SNES_LINESEARCH_FAILED_DOMAIN = (UInt32)(2)
const SNES_LINESEARCH_FAILED_REDUCT = (UInt32)(3)
const SNES_LINESEARCH_FAILED_USER = (UInt32)(4)
const SNES_LINESEARCH_FAILED_FUNCTION = (UInt32)(5)
# end enum SNESLineSearchReason

typealias DMDASNESFunction Ptr{Void}
typealias DMDASNESJacobian Ptr{Void}
typealias DMDASNESObjective Ptr{Void}
typealias SNESMSType Ptr{Uint8}

# begin enum ANONYMOUS_70
typealias ANONYMOUS_70 Uint32
const SNES_NGMRES_RESTART_NONE = (UInt32)(0)
const SNES_NGMRES_RESTART_PERIODIC = (UInt32)(1)
const SNES_NGMRES_RESTART_DIFFERENCE = (UInt32)(2)
# end enum ANONYMOUS_70

# begin enum SNESNGMRESRestartType
typealias SNESNGMRESRestartType Uint32
const SNES_NGMRES_RESTART_NONE = (UInt32)(0)
const SNES_NGMRES_RESTART_PERIODIC = (UInt32)(1)
const SNES_NGMRES_RESTART_DIFFERENCE = (UInt32)(2)
# end enum SNESNGMRESRestartType

# begin enum ANONYMOUS_71
typealias ANONYMOUS_71 Uint32
const SNES_NGMRES_SELECT_NONE = (UInt32)(0)
const SNES_NGMRES_SELECT_DIFFERENCE = (UInt32)(1)
const SNES_NGMRES_SELECT_LINESEARCH = (UInt32)(2)
# end enum ANONYMOUS_71

# begin enum SNESNGMRESSelectType
typealias SNESNGMRESSelectType Uint32
const SNES_NGMRES_SELECT_NONE = (UInt32)(0)
const SNES_NGMRES_SELECT_DIFFERENCE = (UInt32)(1)
const SNES_NGMRES_SELECT_LINESEARCH = (UInt32)(2)
# end enum SNESNGMRESSelectType

# begin enum ANONYMOUS_72
typealias ANONYMOUS_72 Uint32
const SNES_NCG_FR = (UInt32)(0)
const SNES_NCG_PRP = (UInt32)(1)
const SNES_NCG_HS = (UInt32)(2)
const SNES_NCG_DY = (UInt32)(3)
const SNES_NCG_CD = (UInt32)(4)
# end enum ANONYMOUS_72

# begin enum SNESNCGType
typealias SNESNCGType Uint32
const SNES_NCG_FR = (UInt32)(0)
const SNES_NCG_PRP = (UInt32)(1)
const SNES_NCG_HS = (UInt32)(2)
const SNES_NCG_DY = (UInt32)(3)
const SNES_NCG_CD = (UInt32)(4)
# end enum SNESNCGType

# begin enum ANONYMOUS_73
typealias ANONYMOUS_73 Uint32
const SNES_QN_SCALE_DEFAULT = (UInt32)(0)
const SNES_QN_SCALE_NONE = (UInt32)(1)
const SNES_QN_SCALE_SHANNO = (UInt32)(2)
const SNES_QN_SCALE_LINESEARCH = (UInt32)(3)
const SNES_QN_SCALE_JACOBIAN = (UInt32)(4)
# end enum ANONYMOUS_73

# begin enum SNESQNScaleType
typealias SNESQNScaleType Uint32
const SNES_QN_SCALE_DEFAULT = (UInt32)(0)
const SNES_QN_SCALE_NONE = (UInt32)(1)
const SNES_QN_SCALE_SHANNO = (UInt32)(2)
const SNES_QN_SCALE_LINESEARCH = (UInt32)(3)
const SNES_QN_SCALE_JACOBIAN = (UInt32)(4)
# end enum SNESQNScaleType

# begin enum ANONYMOUS_74
typealias ANONYMOUS_74 Uint32
const SNES_QN_RESTART_DEFAULT = (UInt32)(0)
const SNES_QN_RESTART_NONE = (UInt32)(1)
const SNES_QN_RESTART_POWELL = (UInt32)(2)
const SNES_QN_RESTART_PERIODIC = (UInt32)(3)
# end enum ANONYMOUS_74

# begin enum SNESQNRestartType
typealias SNESQNRestartType Uint32
const SNES_QN_RESTART_DEFAULT = (UInt32)(0)
const SNES_QN_RESTART_NONE = (UInt32)(1)
const SNES_QN_RESTART_POWELL = (UInt32)(2)
const SNES_QN_RESTART_PERIODIC = (UInt32)(3)
# end enum SNESQNRestartType

# begin enum ANONYMOUS_75
typealias ANONYMOUS_75 Uint32
const SNES_QN_LBFGS = (UInt32)(0)
const SNES_QN_BROYDEN = (UInt32)(1)
const SNES_QN_BADBROYDEN = (UInt32)(2)
# end enum ANONYMOUS_75

# begin enum SNESQNType
typealias SNESQNType Uint32
const SNES_QN_LBFGS = (UInt32)(0)
const SNES_QN_BROYDEN = (UInt32)(1)
const SNES_QN_BADBROYDEN = (UInt32)(2)
# end enum SNESQNType

# begin enum ANONYMOUS_76
typealias ANONYMOUS_76 Uint32
const SNES_COMPOSITE_ADDITIVE = (UInt32)(0)
const SNES_COMPOSITE_MULTIPLICATIVE = (UInt32)(1)
const SNES_COMPOSITE_ADDITIVEOPTIMAL = (UInt32)(2)
# end enum ANONYMOUS_76

# begin enum SNESCompositeType
typealias SNESCompositeType Uint32
const SNES_COMPOSITE_ADDITIVE = (UInt32)(0)
const SNES_COMPOSITE_MULTIPLICATIVE = (UInt32)(1)
const SNES_COMPOSITE_ADDITIVEOPTIMAL = (UInt32)(2)
# end enum SNESCompositeType

# begin enum ANONYMOUS_77
typealias ANONYMOUS_77 Uint32
const SNES_FAS_MULTIPLICATIVE = (UInt32)(0)
const SNES_FAS_ADDITIVE = (UInt32)(1)
const SNES_FAS_FULL = (UInt32)(2)
const SNES_FAS_KASKADE = (UInt32)(3)
# end enum ANONYMOUS_77

# begin enum SNESFASType
typealias SNESFASType Uint32
const SNES_FAS_MULTIPLICATIVE = (UInt32)(0)
const SNES_FAS_ADDITIVE = (UInt32)(1)
const SNES_FAS_FULL = (UInt32)(2)
const SNES_FAS_KASKADE = (UInt32)(3)
# end enum SNESFASType

type _p_TS
end

typealias TS Ptr{_p_TS}
typealias TSType Ptr{Uint8}

# begin enum ANONYMOUS_78
typealias ANONYMOUS_78 Uint32
const TS_LINEAR = (UInt32)(0)
const TS_NONLINEAR = (UInt32)(1)
# end enum ANONYMOUS_78

# begin enum TSProblemType
typealias TSProblemType Uint32
const TS_LINEAR = (UInt32)(0)
const TS_NONLINEAR = (UInt32)(1)
# end enum TSProblemType

# begin enum ANONYMOUS_79
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
# end enum ANONYMOUS_79

# begin enum TSEquationType
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
# end enum TSEquationType

# begin enum ANONYMOUS_80
typealias ANONYMOUS_80 Cint
const TS_CONVERGED_ITERATING = (Int32)(0)
const TS_CONVERGED_TIME = (Int32)(1)
const TS_CONVERGED_ITS = (Int32)(2)
const TS_CONVERGED_USER = (Int32)(3)
const TS_CONVERGED_EVENT = (Int32)(4)
const TS_DIVERGED_NONLINEAR_SOLVE = (Int32)(-1)
const TS_DIVERGED_STEP_REJECTED = (Int32)(-2)
# end enum ANONYMOUS_80

# begin enum TSConvergedReason
typealias TSConvergedReason Cint
const TS_CONVERGED_ITERATING = (Int32)(0)
const TS_CONVERGED_TIME = (Int32)(1)
const TS_CONVERGED_ITS = (Int32)(2)
const TS_CONVERGED_USER = (Int32)(3)
const TS_CONVERGED_EVENT = (Int32)(4)
const TS_DIVERGED_NONLINEAR_SOLVE = (Int32)(-1)
const TS_DIVERGED_STEP_REJECTED = (Int32)(-2)
# end enum TSConvergedReason

# begin enum ANONYMOUS_81
typealias ANONYMOUS_81 Uint32
const TS_EXACTFINALTIME_STEPOVER = (UInt32)(0)
const TS_EXACTFINALTIME_INTERPOLATE = (UInt32)(1)
const TS_EXACTFINALTIME_MATCHSTEP = (UInt32)(2)
# end enum ANONYMOUS_81

# begin enum TSExactFinalTimeOption
typealias TSExactFinalTimeOption Uint32
const TS_EXACTFINALTIME_STEPOVER = (UInt32)(0)
const TS_EXACTFINALTIME_INTERPOLATE = (UInt32)(1)
const TS_EXACTFINALTIME_MATCHSTEP = (UInt32)(2)
# end enum TSExactFinalTimeOption

type _p_TSTrajectory
end

typealias TSTrajectory Ptr{_p_TSTrajectory}
typealias TSTrajectoryType Ptr{Uint8}

type _n_TSMonitorDrawCtx
end

typealias TSMonitorDrawCtx Ptr{_n_TSMonitorDrawCtx}
typealias TSRHSFunction Ptr{Void}
typealias TSRHSJacobian Ptr{Void}
typealias TSSolutionFunction Ptr{Void}
typealias TSIFunction Ptr{Void}
typealias TSIJacobian Ptr{Void}
typealias DMDATSRHSFunctionLocal Ptr{Void}
typealias DMDATSRHSJacobianLocal Ptr{Void}
typealias DMDATSIFunctionLocal Ptr{Void}
typealias DMDATSIJacobianLocal Ptr{Void}

type _n_TSMonitorLGCtx
end

typealias TSMonitorLGCtx Ptr{_n_TSMonitorLGCtx}

type TSMonitorDMDARayCtx
    ray::Vec
    scatter::VecScatter
    viewer::PetscViewer
    lgctx::TSMonitorLGCtx
end

type _n_TSMonitorEnvelopeCtx
end

typealias TSMonitorEnvelopeCtx Ptr{_n_TSMonitorEnvelopeCtx}

type _n_TSMonitorSPEigCtx
end

typealias TSMonitorSPEigCtx Ptr{_n_TSMonitorSPEigCtx}
typealias TSSSPType Ptr{Uint8}

type _p_TSAdapt
end

typealias TSAdapt Ptr{_p_TSAdapt}
typealias TSAdaptType Ptr{Uint8}

type _p_TSGLAdapt
end

typealias TSGLAdapt Ptr{_p_TSGLAdapt}
typealias TSGLAdaptType Ptr{Uint8}
typealias TSGLAcceptType Ptr{Uint8}
typealias TSGLAcceptFunction Ptr{Void}
typealias TSGLType Ptr{Uint8}
typealias TSRKType Ptr{Uint8}
typealias TSARKIMEXType Ptr{Uint8}
typealias TSRosWType Ptr{Uint8}

# begin enum ANONYMOUS_82
typealias ANONYMOUS_82 Uint32
const TAO_SUBSET_SUBVEC = (UInt32)(0)
const TAO_SUBSET_MASK = (UInt32)(1)
const TAO_SUBSET_MATRIXFREE = (UInt32)(2)
# end enum ANONYMOUS_82

# begin enum TaoSubsetType
typealias TaoSubsetType Uint32
const TAO_SUBSET_SUBVEC = (UInt32)(0)
const TAO_SUBSET_MASK = (UInt32)(1)
const TAO_SUBSET_MATRIXFREE = (UInt32)(2)
# end enum TaoSubsetType

type _p_Tao
end

typealias Tao Ptr{_p_Tao}

# begin enum ANONYMOUS_83
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
# end enum ANONYMOUS_83

# begin enum TaoConvergedReason
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
# end enum TaoConvergedReason

type _p_TaoLineSearch
end

typealias TaoLineSearch Ptr{_p_TaoLineSearch}

# begin enum ANONYMOUS_84
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
# end enum ANONYMOUS_84

# begin enum TaoLineSearchConvergedReason
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
# end enum TaoLineSearchConvergedReason
