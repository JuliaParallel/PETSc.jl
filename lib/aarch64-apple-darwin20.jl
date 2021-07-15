const __darwin_off_t = Int64

const off_t = __darwin_off_t

const PetscErrorCode = Cint

const PetscClassId = Cint

const PetscMPIInt = Cint

@enum PetscEnum::UInt32 begin
    ENUM_DUMMY = 0
end

const PetscShort = Cshort

const PetscChar = Cchar

const PetscFloat = Cfloat

const PetscInt64 = Int64

const PetscInt = PetscInt64

const PetscBLASInt = Cint

@enum PetscBool::UInt32 begin
    PETSC_FALSE = 0
    PETSC_TRUE = 1
end

const PetscReal = Cfloat

const PetscComplex = ComplexF32

const PetscScalar = PetscComplex

@enum PetscCopyMode::UInt32 begin
    PETSC_COPY_VALUES = 0
    PETSC_OWN_POINTER = 1
    PETSC_USE_POINTER = 2
end

const PetscLogDouble = Cdouble

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

mutable struct _p_PetscToken end

const PetscToken = Ptr{_p_PetscToken}

mutable struct _p_PetscObject end

const PetscObject = Ptr{_p_PetscObject}

const PetscObjectId = PetscInt64

const PetscObjectState = PetscInt64

mutable struct _n_PetscFunctionList end

const PetscFunctionList = Ptr{_n_PetscFunctionList}

@enum PetscFileMode::Int32 begin
    FILE_MODE_UNDEFINED = -1
    FILE_MODE_READ = 0
    FILE_MODE_WRITE = 1
    FILE_MODE_APPEND = 2
    FILE_MODE_UPDATE = 3
    FILE_MODE_APPEND_UPDATE = 4
end

const PetscDLHandle = Ptr{Cvoid}

@enum PetscDLMode::UInt32 begin
    PETSC_DL_DECIDE = 0
    PETSC_DL_NOW = 1
    PETSC_DL_LOCAL = 2
end

mutable struct _n_PetscObjectList end

const PetscObjectList = Ptr{_n_PetscObjectList}

mutable struct _n_PetscDLLibrary end

const PetscDLLibrary = Ptr{_n_PetscDLLibrary}

mutable struct _p_PetscContainer end

const PetscContainer = Ptr{_p_PetscContainer}

mutable struct _p_PetscRandom end

const PetscRandom = Ptr{_p_PetscRandom}

@enum PetscBinarySeekType::UInt32 begin
    PETSC_BINARY_SEEK_SET = 0
    PETSC_BINARY_SEEK_CUR = 1
    PETSC_BINARY_SEEK_END = 2
end

@enum PetscBuildTwoSidedType::Int32 begin
    PETSC_BUILDTWOSIDED_NOTSET = -1
    PETSC_BUILDTWOSIDED_ALLREDUCE = 0
    PETSC_BUILDTWOSIDED_IBARRIER = 1
    PETSC_BUILDTWOSIDED_REDSCATTER = 2
end

@enum InsertMode::UInt32 begin
    NOT_SET_VALUES = 0
    INSERT_VALUES = 1
    ADD_VALUES = 2
    MAX_VALUES = 3
    MIN_VALUES = 4
    INSERT_ALL_VALUES = 5
    ADD_ALL_VALUES = 6
    INSERT_BC_VALUES = 7
    ADD_BC_VALUES = 8
end

@enum PetscSubcommType::UInt32 begin
    PETSC_SUBCOMM_GENERAL = 0
    PETSC_SUBCOMM_CONTIGUOUS = 1
    PETSC_SUBCOMM_INTERLACED = 2
end

mutable struct _n_PetscSubcomm
    parent::MPI_Comm
    dupparent::MPI_Comm
    child::MPI_Comm
    n::PetscMPIInt
    color::PetscMPIInt
    subsize::Ptr{PetscMPIInt}
    type::PetscSubcommType
    subcommprefix::Ptr{Cchar}
    _n_PetscSubcomm() = new()
end

const PetscSubcomm = Ptr{_n_PetscSubcomm}

mutable struct _PetscHeap end

const PetscHeap = Ptr{_PetscHeap}

mutable struct _n_PetscShmComm end

const PetscShmComm = Ptr{_n_PetscShmComm}

mutable struct _n_PetscOmpCtrl end

const PetscOmpCtrl = Ptr{_n_PetscOmpCtrl}

mutable struct _n_PetscSegBuffer end

const PetscSegBuffer = Ptr{_n_PetscSegBuffer}

mutable struct _n_PetscOptionsHelpPrinted end

const PetscOptionsHelpPrinted = Ptr{_n_PetscOptionsHelpPrinted}

@enum PetscMemType::UInt32 begin
    PETSC_MEMTYPE_HOST = 0
    PETSC_MEMTYPE_DEVICE = 1
    # PETSC_MEMTYPE_CUDA = 1
    PETSC_MEMTYPE_NVSHMEM = 17
    PETSC_MEMTYPE_HIP = 3
end

function PetscSignReal(a)
    @ccall libpetsc.PetscSignReal(a::PetscReal)::PetscReal
end

function PetscCMPLX(x, y)
    @ccall libpetsc.PetscCMPLX(x::PetscReal, y::PetscReal)::PetscComplex
end

@enum PetscScalarPrecision::UInt32 begin
    PETSC_SCALAR_DOUBLE = 0
    PETSC_SCALAR_SINGLE = 1
    PETSC_SCALAR_LONG_DOUBLE = 2
    PETSC_SCALAR_HALF = 3
end

function PetscIsInfReal(arg1)
    @ccall libpetsc.PetscIsInfReal(arg1::PetscReal)::PetscBool
end

function PetscIsNanReal(arg1)
    @ccall libpetsc.PetscIsNanReal(arg1::PetscReal)::PetscBool
end

function PetscIsNormalReal(arg1)
    @ccall libpetsc.PetscIsNormalReal(arg1::PetscReal)::PetscBool
end

function PetscIsInfOrNanReal(v)
    @ccall libpetsc.PetscIsInfOrNanReal(v::PetscReal)::PetscBool
end

function PetscIsInfScalar(v)
    @ccall libpetsc.PetscIsInfScalar(v::PetscScalar)::PetscBool
end

function PetscIsNanScalar(v)
    @ccall libpetsc.PetscIsNanScalar(v::PetscScalar)::PetscBool
end

function PetscIsInfOrNanScalar(v)
    @ccall libpetsc.PetscIsInfOrNanScalar(v::PetscScalar)::PetscBool
end

function PetscIsNormalScalar(v)
    @ccall libpetsc.PetscIsNormalScalar(v::PetscScalar)::PetscBool
end

function PetscIsCloseAtTol(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscIsCloseAtTol(arg1::PetscReal, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal)::PetscBool
end

function PetscEqualReal(arg1, arg2)
    @ccall libpetsc.PetscEqualReal(arg1::PetscReal, arg2::PetscReal)::PetscBool
end

function PetscEqualScalar(arg1, arg2)
    @ccall libpetsc.PetscEqualScalar(arg1::PetscScalar, arg2::PetscScalar)::PetscBool
end

const MatScalar = PetscScalar

const MatReal = PetscReal

mutable struct petsc_mpiu_2scalar
    a::PetscScalar
    b::PetscScalar
    petsc_mpiu_2scalar() = new()
end

mutable struct petsc_mpiu_2int
    a::PetscInt
    b::PetscInt
    petsc_mpiu_2int() = new()
end

function PetscPowInt(base, power)
    @ccall libpetsc.PetscPowInt(base::PetscInt, power::PetscInt)::PetscInt
end

function PetscPowInt64(base, power)
    @ccall libpetsc.PetscPowInt64(base::PetscInt, power::PetscInt)::PetscInt64
end

function PetscPowRealInt(base, power)
    @ccall libpetsc.PetscPowRealInt(base::PetscReal, power::PetscInt)::PetscReal
end

function PetscPowScalarInt(base, power)
    @ccall libpetsc.PetscPowScalarInt(base::PetscScalar, power::PetscInt)::PetscScalar
end

function PetscPowScalarReal(base, power)
    @ccall libpetsc.PetscPowScalarReal(base::PetscScalar, power::PetscReal)::PetscScalar
end

function PetscLinearRegression(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscLinearRegression(arg1::PetscInt, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal})::PetscErrorCode
end

function PetscSetHelpVersionFunctions(arg1, arg2)
    @ccall libpetsc.PetscSetHelpVersionFunctions(arg1::Ptr{Cvoid}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscCommDuplicate(arg1, arg2, arg3)
    @ccall libpetsc.PetscCommDuplicate(arg1::MPI_Comm, arg2::Ptr{MPI_Comm}, arg3::Ptr{Cint})::PetscErrorCode
end

function PetscCommDestroy(arg1)
    @ccall libpetsc.PetscCommDestroy(arg1::Ptr{MPI_Comm})::PetscErrorCode
end

function PetscMallocSetCoalesce(arg1)
    @ccall libpetsc.PetscMallocSetCoalesce(arg1::PetscBool)::PetscErrorCode
end

function PetscMallocSet(arg1, arg2, arg3)
    @ccall libpetsc.PetscMallocSet(arg1::Ptr{Cvoid}, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function PetscMallocClear()
    @ccall libpetsc.PetscMallocClear()::PetscErrorCode
end

function PetscMallocSetDRAM()
    @ccall libpetsc.PetscMallocSetDRAM()::PetscErrorCode
end

function PetscMallocResetDRAM()
    @ccall libpetsc.PetscMallocResetDRAM()::PetscErrorCode
end

function PetscMallocDump(arg1)
    @ccall libpetsc.PetscMallocDump(arg1::Ptr{Libc.FILE})::PetscErrorCode
end

function PetscMallocView(arg1)
    @ccall libpetsc.PetscMallocView(arg1::Ptr{Libc.FILE})::PetscErrorCode
end

function PetscMallocGetCurrentUsage(arg1)
    @ccall libpetsc.PetscMallocGetCurrentUsage(arg1::Ptr{PetscLogDouble})::PetscErrorCode
end

function PetscMallocGetMaximumUsage(arg1)
    @ccall libpetsc.PetscMallocGetMaximumUsage(arg1::Ptr{PetscLogDouble})::PetscErrorCode
end

function PetscMallocPushMaximumUsage(arg1)
    @ccall libpetsc.PetscMallocPushMaximumUsage(arg1::Cint)::PetscErrorCode
end

function PetscMallocPopMaximumUsage(arg1, arg2)
    @ccall libpetsc.PetscMallocPopMaximumUsage(arg1::Cint, arg2::Ptr{PetscLogDouble})::PetscErrorCode
end

function PetscMallocSetDebug(arg1, arg2)
    @ccall libpetsc.PetscMallocSetDebug(arg1::PetscBool, arg2::PetscBool)::PetscErrorCode
end

function PetscMallocGetDebug(arg1, arg2, arg3)
    @ccall libpetsc.PetscMallocGetDebug(arg1::Ptr{PetscBool}, arg2::Ptr{PetscBool}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function PetscMallocValidate(arg1, arg2, arg3)
    @ccall libpetsc.PetscMallocValidate(arg1::Cint, arg2::Ptr{Cchar}, arg3::Ptr{Cchar})::PetscErrorCode
end

function PetscMallocViewSet(arg1)
    @ccall libpetsc.PetscMallocViewSet(arg1::PetscLogDouble)::PetscErrorCode
end

function PetscMallocViewGet(arg1)
    @ccall libpetsc.PetscMallocViewGet(arg1::Ptr{PetscBool})::PetscErrorCode
end

function PetscMallocLogRequestedSizeSet(arg1)
    @ccall libpetsc.PetscMallocLogRequestedSizeSet(arg1::PetscBool)::PetscErrorCode
end

function PetscMallocLogRequestedSizeGet(arg1)
    @ccall libpetsc.PetscMallocLogRequestedSizeGet(arg1::Ptr{PetscBool})::PetscErrorCode
end

function PetscDataTypeToMPIDataType(arg1, arg2)
    @ccall libpetsc.PetscDataTypeToMPIDataType(arg1::PetscDataType, arg2::Ptr{MPI_Datatype})::PetscErrorCode
end

function PetscMPIDataTypeToPetscDataType(arg1, arg2)
    @ccall libpetsc.PetscMPIDataTypeToPetscDataType(arg1::MPI_Datatype, arg2::Ptr{PetscDataType})::PetscErrorCode
end

function PetscDataTypeGetSize(arg1, arg2)
    @ccall libpetsc.PetscDataTypeGetSize(arg1::PetscDataType, arg2::Ptr{Csize_t})::PetscErrorCode
end

function PetscDataTypeFromString(arg1, arg2, arg3)
    @ccall libpetsc.PetscDataTypeFromString(arg1::Ptr{Cchar}, arg2::Ptr{PetscDataType}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function PetscMemcmp(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscMemcmp(arg1::Ptr{Cvoid}, arg2::Ptr{Cvoid}, arg3::Csize_t, arg4::Ptr{PetscBool})::PetscErrorCode
end

function PetscStrlen(arg1, arg2)
    @ccall libpetsc.PetscStrlen(arg1::Ptr{Cchar}, arg2::Ptr{Csize_t})::PetscErrorCode
end

function PetscStrToArray(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscStrToArray(arg1::Ptr{Cchar}, arg2::Cchar, arg3::Ptr{Cint}, arg4::Ptr{Ptr{Ptr{Cchar}}})::PetscErrorCode
end

function PetscStrToArrayDestroy(arg1, arg2)
    @ccall libpetsc.PetscStrToArrayDestroy(arg1::Cint, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscStrcmp(arg1, arg2, arg3)
    @ccall libpetsc.PetscStrcmp(arg1::Ptr{Cchar}, arg2::Ptr{Cchar}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function PetscStrgrt(arg1, arg2, arg3)
    @ccall libpetsc.PetscStrgrt(arg1::Ptr{Cchar}, arg2::Ptr{Cchar}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function PetscStrcasecmp(arg1, arg2, arg3)
    @ccall libpetsc.PetscStrcasecmp(arg1::Ptr{Cchar}, arg2::Ptr{Cchar}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function PetscStrncmp(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscStrncmp(arg1::Ptr{Cchar}, arg2::Ptr{Cchar}, arg3::Csize_t, arg4::Ptr{PetscBool})::PetscErrorCode
end

function PetscStrcpy(arg1, arg2)
    @ccall libpetsc.PetscStrcpy(arg1::Ptr{Cchar}, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscStrcat(arg1, arg2)
    @ccall libpetsc.PetscStrcat(arg1::Ptr{Cchar}, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscStrlcat(arg1, arg2, arg3)
    @ccall libpetsc.PetscStrlcat(arg1::Ptr{Cchar}, arg2::Ptr{Cchar}, arg3::Csize_t)::PetscErrorCode
end

function PetscStrncpy(arg1, arg2, arg3)
    @ccall libpetsc.PetscStrncpy(arg1::Ptr{Cchar}, arg2::Ptr{Cchar}, arg3::Csize_t)::PetscErrorCode
end

function PetscStrchr(arg1, arg2, arg3)
    @ccall libpetsc.PetscStrchr(arg1::Ptr{Cchar}, arg2::Cchar, arg3::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscStrtolower(arg1)
    @ccall libpetsc.PetscStrtolower(arg1::Ptr{Cchar})::PetscErrorCode
end

function PetscStrtoupper(arg1)
    @ccall libpetsc.PetscStrtoupper(arg1::Ptr{Cchar})::PetscErrorCode
end

function PetscStrrchr(arg1, arg2, arg3)
    @ccall libpetsc.PetscStrrchr(arg1::Ptr{Cchar}, arg2::Cchar, arg3::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscStrstr(arg1, arg2, arg3)
    @ccall libpetsc.PetscStrstr(arg1::Ptr{Cchar}, arg2::Ptr{Cchar}, arg3::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscStrrstr(arg1, arg2, arg3)
    @ccall libpetsc.PetscStrrstr(arg1::Ptr{Cchar}, arg2::Ptr{Cchar}, arg3::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscStrendswith(arg1, arg2, arg3)
    @ccall libpetsc.PetscStrendswith(arg1::Ptr{Cchar}, arg2::Ptr{Cchar}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function PetscStrbeginswith(arg1, arg2, arg3)
    @ccall libpetsc.PetscStrbeginswith(arg1::Ptr{Cchar}, arg2::Ptr{Cchar}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function PetscStrendswithwhich(arg1, arg2, arg3)
    @ccall libpetsc.PetscStrendswithwhich(arg1::Ptr{Cchar}, arg2::Ptr{Ptr{Cchar}}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscStrallocpy(arg1, arg2)
    @ccall libpetsc.PetscStrallocpy(arg1::Ptr{Cchar}, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscStrArrayallocpy(arg1, arg2)
    @ccall libpetsc.PetscStrArrayallocpy(arg1::Ptr{Ptr{Cchar}}, arg2::Ptr{Ptr{Ptr{Cchar}}})::PetscErrorCode
end

function PetscStrArrayDestroy(arg1)
    @ccall libpetsc.PetscStrArrayDestroy(arg1::Ptr{Ptr{Ptr{Cchar}}})::PetscErrorCode
end

function PetscStrNArrayallocpy(arg1, arg2, arg3)
    @ccall libpetsc.PetscStrNArrayallocpy(arg1::PetscInt, arg2::Ptr{Ptr{Cchar}}, arg3::Ptr{Ptr{Ptr{Cchar}}})::PetscErrorCode
end

function PetscStrNArrayDestroy(arg1, arg2)
    @ccall libpetsc.PetscStrNArrayDestroy(arg1::PetscInt, arg2::Ptr{Ptr{Ptr{Cchar}}})::PetscErrorCode
end

function PetscStrreplace(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscStrreplace(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Csize_t)::PetscErrorCode
end

function PetscStrcmpNoError(arg1, arg2, arg3)
    @ccall libpetsc.PetscStrcmpNoError(arg1::Ptr{Cchar}, arg2::Ptr{Cchar}, arg3::Ptr{PetscBool})::Cvoid
end

function PetscTokenCreate(arg1, arg2, arg3)
    @ccall libpetsc.PetscTokenCreate(arg1::Ptr{Cchar}, arg2::Cchar, arg3::Ptr{PetscToken})::PetscErrorCode
end

function PetscTokenFind(arg1, arg2)
    @ccall libpetsc.PetscTokenFind(arg1::PetscToken, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscTokenDestroy(arg1)
    @ccall libpetsc.PetscTokenDestroy(arg1::Ptr{PetscToken})::PetscErrorCode
end

function PetscStrInList(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscStrInList(arg1::Ptr{Cchar}, arg2::Ptr{Cchar}, arg3::Cchar, arg4::Ptr{PetscBool})::PetscErrorCode
end

function PetscEListFind(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscEListFind(arg1::PetscInt, arg2::Ptr{Ptr{Cchar}}, arg3::Ptr{Cchar}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscBool})::PetscErrorCode
end

function PetscEnumFind(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscEnumFind(arg1::Ptr{Ptr{Cchar}}, arg2::Ptr{Cchar}, arg3::Ptr{PetscEnum}, arg4::Ptr{PetscBool})::PetscErrorCode
end

function PetscMaxSum(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscMaxSum(arg1::MPI_Comm, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function MPIULong_Send(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MPIULong_Send(arg1::Ptr{Cvoid}, arg2::PetscInt, arg3::MPI_Datatype, arg4::PetscMPIInt, arg5::PetscMPIInt, arg6::MPI_Comm)::PetscErrorCode
end

function MPIULong_Recv(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MPIULong_Recv(arg1::Ptr{Cvoid}, arg2::PetscInt, arg3::MPI_Datatype, arg4::PetscMPIInt, arg5::PetscMPIInt, arg6::MPI_Comm)::PetscErrorCode
end

function PetscAbortFindSourceFile_Private(arg1, arg2)
    @ccall libpetsc.PetscAbortFindSourceFile_Private(arg1::Ptr{Cchar}, arg2::Ptr{PetscInt})::PetscErrorCode
end

@enum PetscErrorType::UInt32 begin
    PETSC_ERROR_INITIAL = 0
    PETSC_ERROR_REPEAT = 1
    PETSC_ERROR_IN_CXX = 2
end

function PetscErrorPrintfInitialize()
    @ccall libpetsc.PetscErrorPrintfInitialize()::PetscErrorCode
end

function PetscErrorMessage(arg1, arg2, arg3)
    @ccall libpetsc.PetscErrorMessage(arg1::Cint, arg2::Ptr{Ptr{Cchar}}, arg3::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscTraceBackErrorHandler(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscTraceBackErrorHandler(arg1::MPI_Comm, arg2::Cint, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::PetscErrorCode, arg6::PetscErrorType, arg7::Ptr{Cchar}, arg8::Ptr{Cvoid})::PetscErrorCode
end

function PetscIgnoreErrorHandler(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscIgnoreErrorHandler(arg1::MPI_Comm, arg2::Cint, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::PetscErrorCode, arg6::PetscErrorType, arg7::Ptr{Cchar}, arg8::Ptr{Cvoid})::PetscErrorCode
end

function PetscEmacsClientErrorHandler(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscEmacsClientErrorHandler(arg1::MPI_Comm, arg2::Cint, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::PetscErrorCode, arg6::PetscErrorType, arg7::Ptr{Cchar}, arg8::Ptr{Cvoid})::PetscErrorCode
end

function PetscMPIAbortErrorHandler(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscMPIAbortErrorHandler(arg1::MPI_Comm, arg2::Cint, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::PetscErrorCode, arg6::PetscErrorType, arg7::Ptr{Cchar}, arg8::Ptr{Cvoid})::PetscErrorCode
end

function PetscAbortErrorHandler(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscAbortErrorHandler(arg1::MPI_Comm, arg2::Cint, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::PetscErrorCode, arg6::PetscErrorType, arg7::Ptr{Cchar}, arg8::Ptr{Cvoid})::PetscErrorCode
end

function PetscAttachDebuggerErrorHandler(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscAttachDebuggerErrorHandler(arg1::MPI_Comm, arg2::Cint, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::PetscErrorCode, arg6::PetscErrorType, arg7::Ptr{Cchar}, arg8::Ptr{Cvoid})::PetscErrorCode
end

function PetscReturnErrorHandler(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscReturnErrorHandler(arg1::MPI_Comm, arg2::Cint, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::PetscErrorCode, arg6::PetscErrorType, arg7::Ptr{Cchar}, arg8::Ptr{Cvoid})::PetscErrorCode
end

function PetscPushErrorHandler(handler, arg2)
    @ccall libpetsc.PetscPushErrorHandler(handler::Ptr{Cvoid}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscPopErrorHandler()
    @ccall libpetsc.PetscPopErrorHandler()::PetscErrorCode
end

function PetscSignalHandlerDefault(arg1, arg2)
    @ccall libpetsc.PetscSignalHandlerDefault(arg1::Cint, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscPushSignalHandler(arg1, arg2)
    @ccall libpetsc.PetscPushSignalHandler(arg1::Ptr{Cvoid}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscPopSignalHandler()
    @ccall libpetsc.PetscPopSignalHandler()::PetscErrorCode
end

function PetscCheckPointerSetIntensity(arg1)
    @ccall libpetsc.PetscCheckPointerSetIntensity(arg1::PetscInt)::PetscErrorCode
end

function PetscSignalSegvCheckPointerOrMpi()
    @ccall libpetsc.PetscSignalSegvCheckPointerOrMpi()::Cvoid
end

function PetscSignalSegvCheckPointer()
    @ccall libpetsc.PetscSignalSegvCheckPointer()::Cvoid
end

@enum PetscFPTrap::UInt32 begin
    PETSC_FP_TRAP_OFF = 0
    PETSC_FP_TRAP_ON = 1
end

function PetscSetFPTrap(arg1)
    @ccall libpetsc.PetscSetFPTrap(arg1::PetscFPTrap)::PetscErrorCode
end

function PetscFPTrapPush(arg1)
    @ccall libpetsc.PetscFPTrapPush(arg1::PetscFPTrap)::PetscErrorCode
end

function PetscFPTrapPop()
    @ccall libpetsc.PetscFPTrapPop()::PetscErrorCode
end

function PetscDetermineInitialFPTrap()
    @ccall libpetsc.PetscDetermineInitialFPTrap()::PetscErrorCode
end

mutable struct PetscStack
    _function::NTuple{64, Ptr{Cchar}}
    file::NTuple{64, Ptr{Cchar}}
    line::NTuple{64, Cint}
    petscroutine::NTuple{64, PetscBool}
    currentsize::Cint
    hotdepth::Cint
    PetscStack() = new()
end

function PetscStackCopy(arg1, arg2)
    @ccall libpetsc.PetscStackCopy(arg1::Ptr{PetscStack}, arg2::Ptr{PetscStack})::PetscErrorCode
end

function PetscStackPrint(arg1, arg2)
    @ccall libpetsc.PetscStackPrint(arg1::Ptr{PetscStack}, arg2::Ptr{Libc.FILE})::PetscErrorCode
end

function PetscStackActive()
    @ccall libpetsc.PetscStackActive()::PetscBool
end

function PetscStackCreate()
    @ccall libpetsc.PetscStackCreate()::PetscErrorCode
end

function PetscStackView(arg1)
    @ccall libpetsc.PetscStackView(arg1::Ptr{Libc.FILE})::PetscErrorCode
end

function PetscStackDestroy()
    @ccall libpetsc.PetscStackDestroy()::PetscErrorCode
end

function PetscClassIdRegister(arg1, arg2)
    @ccall libpetsc.PetscClassIdRegister(arg1::Ptr{Cchar}, arg2::Ptr{PetscClassId})::PetscErrorCode
end

function PetscObjectGetId(arg1, arg2)
    @ccall libpetsc.PetscObjectGetId(arg1::PetscObject, arg2::Ptr{PetscObjectId})::PetscErrorCode
end

function PetscObjectCompareId(arg1, arg2, arg3)
    @ccall libpetsc.PetscObjectCompareId(arg1::PetscObject, arg2::PetscObjectId, arg3::Ptr{PetscBool})::PetscErrorCode
end

function PetscMemoryGetCurrentUsage(arg1)
    @ccall libpetsc.PetscMemoryGetCurrentUsage(arg1::Ptr{PetscLogDouble})::PetscErrorCode
end

function PetscMemoryGetMaximumUsage(arg1)
    @ccall libpetsc.PetscMemoryGetMaximumUsage(arg1::Ptr{PetscLogDouble})::PetscErrorCode
end

function PetscMemorySetGetMaximumUsage()
    @ccall libpetsc.PetscMemorySetGetMaximumUsage()::PetscErrorCode
end

function PetscMemoryTrace(arg1)
    @ccall libpetsc.PetscMemoryTrace(arg1::Ptr{Cchar})::PetscErrorCode
end

function PetscSleep(arg1)
    @ccall libpetsc.PetscSleep(arg1::PetscReal)::PetscErrorCode
end

function PetscInitialize(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscInitialize(arg1::Ptr{Cint}, arg2::Ptr{Ptr{Ptr{Cchar}}}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar})::PetscErrorCode
end

function PetscInitializeNoPointers(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscInitializeNoPointers(arg1::Cint, arg2::Ptr{Ptr{Cchar}}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar})::PetscErrorCode
end

function PetscInitializeNoArguments()
    @ccall libpetsc.PetscInitializeNoArguments()::PetscErrorCode
end

function PetscInitialized(arg1)
    @ccall libpetsc.PetscInitialized(arg1::Ptr{PetscBool})::PetscErrorCode
end

function PetscFinalized(arg1)
    @ccall libpetsc.PetscFinalized(arg1::Ptr{PetscBool})::PetscErrorCode
end

function PetscFinalize()
    @ccall libpetsc.PetscFinalize()::PetscErrorCode
end

function PetscInitializeFortran()
    @ccall libpetsc.PetscInitializeFortran()::PetscErrorCode
end

function PetscGetArgs(arg1, arg2)
    @ccall libpetsc.PetscGetArgs(arg1::Ptr{Cint}, arg2::Ptr{Ptr{Ptr{Cchar}}})::PetscErrorCode
end

function PetscGetArguments(arg1)
    @ccall libpetsc.PetscGetArguments(arg1::Ptr{Ptr{Ptr{Cchar}}})::PetscErrorCode
end

function PetscFreeArguments(arg1)
    @ccall libpetsc.PetscFreeArguments(arg1::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscEnd()
    @ccall libpetsc.PetscEnd()::PetscErrorCode
end

function PetscSysInitializePackage()
    @ccall libpetsc.PetscSysInitializePackage()::PetscErrorCode
end

function PetscPythonInitialize(arg1, arg2)
    @ccall libpetsc.PetscPythonInitialize(arg1::Ptr{Cchar}, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscPythonFinalize()
    @ccall libpetsc.PetscPythonFinalize()::PetscErrorCode
end

function PetscPythonPrintError()
    @ccall libpetsc.PetscPythonPrintError()::PetscErrorCode
end

function PetscPythonMonitorSet(arg1, arg2)
    @ccall libpetsc.PetscPythonMonitorSet(arg1::PetscObject, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscMonitorCompare(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscMonitorCompare(arg1::Ptr{Cvoid}, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid}, arg6::Ptr{Cvoid}, arg7::Ptr{PetscBool})::PetscErrorCode
end

# typedef void ( * * PetscVoidStarFunction ) ( void )
const PetscVoidStarFunction = Ptr{Ptr{Cvoid}}

# typedef void ( * PetscVoidFunction ) ( void )
const PetscVoidFunction = Ptr{Cvoid}

# typedef PetscErrorCode ( * PetscErrorCodeFunction ) ( void )
const PetscErrorCodeFunction = Ptr{Cvoid}

function PetscObjectDestroy(arg1)
    @ccall libpetsc.PetscObjectDestroy(arg1::Ptr{PetscObject})::PetscErrorCode
end

function PetscObjectGetComm(arg1, arg2)
    @ccall libpetsc.PetscObjectGetComm(arg1::PetscObject, arg2::Ptr{MPI_Comm})::PetscErrorCode
end

function PetscObjectGetClassId(arg1, arg2)
    @ccall libpetsc.PetscObjectGetClassId(arg1::PetscObject, arg2::Ptr{PetscClassId})::PetscErrorCode
end

function PetscObjectGetClassName(arg1, arg2)
    @ccall libpetsc.PetscObjectGetClassName(arg1::PetscObject, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscObjectSetType(arg1, arg2)
    @ccall libpetsc.PetscObjectSetType(arg1::PetscObject, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscObjectGetType(arg1, arg2)
    @ccall libpetsc.PetscObjectGetType(arg1::PetscObject, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscObjectSetName(arg1, arg2)
    @ccall libpetsc.PetscObjectSetName(arg1::PetscObject, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscObjectGetName(arg1, arg2)
    @ccall libpetsc.PetscObjectGetName(arg1::PetscObject, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscObjectSetTabLevel(arg1, arg2)
    @ccall libpetsc.PetscObjectSetTabLevel(arg1::PetscObject, arg2::PetscInt)::PetscErrorCode
end

function PetscObjectGetTabLevel(arg1, arg2)
    @ccall libpetsc.PetscObjectGetTabLevel(arg1::PetscObject, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscObjectIncrementTabLevel(arg1, arg2, arg3)
    @ccall libpetsc.PetscObjectIncrementTabLevel(arg1::PetscObject, arg2::PetscObject, arg3::PetscInt)::PetscErrorCode
end

function PetscObjectReference(arg1)
    @ccall libpetsc.PetscObjectReference(arg1::PetscObject)::PetscErrorCode
end

function PetscObjectGetReference(arg1, arg2)
    @ccall libpetsc.PetscObjectGetReference(arg1::PetscObject, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscObjectDereference(arg1)
    @ccall libpetsc.PetscObjectDereference(arg1::PetscObject)::PetscErrorCode
end

function PetscObjectGetNewTag(arg1, arg2)
    @ccall libpetsc.PetscObjectGetNewTag(arg1::PetscObject, arg2::Ptr{PetscMPIInt})::PetscErrorCode
end

function PetscObjectCompose(arg1, arg2, arg3)
    @ccall libpetsc.PetscObjectCompose(arg1::PetscObject, arg2::Ptr{Cchar}, arg3::PetscObject)::PetscErrorCode
end

function PetscObjectRemoveReference(arg1, arg2)
    @ccall libpetsc.PetscObjectRemoveReference(arg1::PetscObject, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscObjectQuery(arg1, arg2, arg3)
    @ccall libpetsc.PetscObjectQuery(arg1::PetscObject, arg2::Ptr{Cchar}, arg3::Ptr{PetscObject})::PetscErrorCode
end

function PetscObjectComposeFunction_Private(arg1, arg2, arg3)
    @ccall libpetsc.PetscObjectComposeFunction_Private(arg1::PetscObject, arg2::Ptr{Cchar}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function PetscObjectSetFromOptions(arg1)
    @ccall libpetsc.PetscObjectSetFromOptions(arg1::PetscObject)::PetscErrorCode
end

function PetscObjectSetUp(arg1)
    @ccall libpetsc.PetscObjectSetUp(arg1::PetscObject)::PetscErrorCode
end

function PetscObjectSetPrintedOptions(arg1)
    @ccall libpetsc.PetscObjectSetPrintedOptions(arg1::PetscObject)::PetscErrorCode
end

function PetscObjectInheritPrintedOptions(arg1, arg2)
    @ccall libpetsc.PetscObjectInheritPrintedOptions(arg1::PetscObject, arg2::PetscObject)::PetscErrorCode
end

function PetscCommGetNewTag(arg1, arg2)
    @ccall libpetsc.PetscCommGetNewTag(arg1::MPI_Comm, arg2::Ptr{PetscMPIInt})::PetscErrorCode
end

mutable struct _p_PetscViewer end

const PetscViewer = Ptr{_p_PetscViewer}

mutable struct _n_PetscOptions end

const PetscOptions = Ptr{_n_PetscOptions}

function PetscOptionsCreate(arg1)
    @ccall libpetsc.PetscOptionsCreate(arg1::Ptr{PetscOptions})::PetscErrorCode
end

function PetscOptionsPush(arg1)
    @ccall libpetsc.PetscOptionsPush(arg1::PetscOptions)::PetscErrorCode
end

function PetscOptionsPop()
    @ccall libpetsc.PetscOptionsPop()::PetscErrorCode
end

function PetscOptionsDestroy(arg1)
    @ccall libpetsc.PetscOptionsDestroy(arg1::Ptr{PetscOptions})::PetscErrorCode
end

function PetscOptionsCreateDefault()
    @ccall libpetsc.PetscOptionsCreateDefault()::PetscErrorCode
end

function PetscOptionsDestroyDefault()
    @ccall libpetsc.PetscOptionsDestroyDefault()::PetscErrorCode
end

function PetscOptionsHasHelp(arg1, arg2)
    @ccall libpetsc.PetscOptionsHasHelp(arg1::PetscOptions, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsHasName(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscOptionsHasName(arg1::PetscOptions, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsGetBool(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscOptionsGetBool(arg1::PetscOptions, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{PetscBool}, arg5::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsGetInt(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscOptionsGetInt(arg1::PetscOptions, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsGetEnum(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscOptionsGetEnum(arg1::PetscOptions, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Ptr{Cchar}}, arg5::Ptr{PetscEnum}, arg6::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsGetEList(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscOptionsGetEList(arg1::PetscOptions, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Ptr{Cchar}}, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsGetReal(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscOptionsGetReal(arg1::PetscOptions, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsGetScalar(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscOptionsGetScalar(arg1::PetscOptions, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{PetscScalar}, arg5::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsGetString(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscOptionsGetString(arg1::PetscOptions, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Csize_t, arg6::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsGetBoolArray(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscOptionsGetBoolArray(arg1::PetscOptions, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{PetscBool}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsGetEnumArray(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscOptionsGetEnumArray(arg1::PetscOptions, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Ptr{Cchar}}, arg5::Ptr{PetscEnum}, arg6::Ptr{PetscInt}, arg7::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsGetIntArray(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscOptionsGetIntArray(arg1::PetscOptions, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsGetRealArray(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscOptionsGetRealArray(arg1::PetscOptions, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsGetScalarArray(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscOptionsGetScalarArray(arg1::PetscOptions, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{PetscScalar}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsGetStringArray(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscOptionsGetStringArray(arg1::PetscOptions, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Ptr{Cchar}}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsValidKey(arg1, arg2)
    @ccall libpetsc.PetscOptionsValidKey(arg1::Ptr{Cchar}, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsSetAlias(arg1, arg2, arg3)
    @ccall libpetsc.PetscOptionsSetAlias(arg1::PetscOptions, arg2::Ptr{Cchar}, arg3::Ptr{Cchar})::PetscErrorCode
end

function PetscOptionsSetValue(arg1, arg2, arg3)
    @ccall libpetsc.PetscOptionsSetValue(arg1::PetscOptions, arg2::Ptr{Cchar}, arg3::Ptr{Cchar})::PetscErrorCode
end

function PetscOptionsClearValue(arg1, arg2)
    @ccall libpetsc.PetscOptionsClearValue(arg1::PetscOptions, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscOptionsFindPair(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscOptionsFindPair(arg1::PetscOptions, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Ptr{Cchar}}, arg5::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsGetAll(arg1, arg2)
    @ccall libpetsc.PetscOptionsGetAll(arg1::PetscOptions, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscOptionsAllUsed(arg1, arg2)
    @ccall libpetsc.PetscOptionsAllUsed(arg1::PetscOptions, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscOptionsUsed(arg1, arg2, arg3)
    @ccall libpetsc.PetscOptionsUsed(arg1::PetscOptions, arg2::Ptr{Cchar}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsLeft(arg1)
    @ccall libpetsc.PetscOptionsLeft(arg1::PetscOptions)::PetscErrorCode
end

function PetscOptionsLeftGet(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscOptionsLeftGet(arg1::PetscOptions, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{Ptr{Cchar}}}, arg4::Ptr{Ptr{Ptr{Cchar}}})::PetscErrorCode
end

function PetscOptionsLeftRestore(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscOptionsLeftRestore(arg1::PetscOptions, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{Ptr{Cchar}}}, arg4::Ptr{Ptr{Ptr{Cchar}}})::PetscErrorCode
end

function PetscOptionsView(arg1, arg2)
    @ccall libpetsc.PetscOptionsView(arg1::PetscOptions, arg2::PetscViewer)::PetscErrorCode
end

function PetscOptionsReject(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscOptionsReject(arg1::PetscOptions, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar})::PetscErrorCode
end

function PetscOptionsInsert(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscOptionsInsert(arg1::PetscOptions, arg2::Ptr{Cint}, arg3::Ptr{Ptr{Ptr{Cchar}}}, arg4::Ptr{Cchar})::PetscErrorCode
end

function PetscOptionsInsertFile(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscOptionsInsertFile(arg1::MPI_Comm, arg2::PetscOptions, arg3::Ptr{Cchar}, arg4::PetscBool)::PetscErrorCode
end

function PetscOptionsInsertFileYAML(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscOptionsInsertFileYAML(arg1::MPI_Comm, arg2::PetscOptions, arg3::Ptr{Cchar}, arg4::PetscBool)::PetscErrorCode
end

function PetscOptionsInsertString(arg1, arg2)
    @ccall libpetsc.PetscOptionsInsertString(arg1::PetscOptions, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscOptionsInsertStringYAML(arg1, arg2)
    @ccall libpetsc.PetscOptionsInsertStringYAML(arg1::PetscOptions, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscOptionsInsertArgs(arg1, arg2, arg3)
    @ccall libpetsc.PetscOptionsInsertArgs(arg1::PetscOptions, arg2::Cint, arg3::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscOptionsClear(arg1)
    @ccall libpetsc.PetscOptionsClear(arg1::PetscOptions)::PetscErrorCode
end

function PetscOptionsPrefixPush(arg1, arg2)
    @ccall libpetsc.PetscOptionsPrefixPush(arg1::PetscOptions, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscOptionsPrefixPop(arg1)
    @ccall libpetsc.PetscOptionsPrefixPop(arg1::PetscOptions)::PetscErrorCode
end

function PetscOptionsGetenv(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscOptionsGetenv(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Csize_t, arg5::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsStringToBool(arg1, arg2)
    @ccall libpetsc.PetscOptionsStringToBool(arg1::Ptr{Cchar}, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsStringToInt(arg1, arg2)
    @ccall libpetsc.PetscOptionsStringToInt(arg1::Ptr{Cchar}, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscOptionsStringToReal(arg1, arg2)
    @ccall libpetsc.PetscOptionsStringToReal(arg1::Ptr{Cchar}, arg2::Ptr{PetscReal})::PetscErrorCode
end

function PetscOptionsStringToScalar(arg1, arg2)
    @ccall libpetsc.PetscOptionsStringToScalar(arg1::Ptr{Cchar}, arg2::Ptr{PetscScalar})::PetscErrorCode
end

function PetscOptionsMonitorSet(arg1, arg2, arg3)
    @ccall libpetsc.PetscOptionsMonitorSet(arg1::Ptr{Cvoid}, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function PetscOptionsMonitorDefault(arg1, arg2, arg3)
    @ccall libpetsc.PetscOptionsMonitorDefault(arg1::Ptr{Cchar}, arg2::Ptr{Cchar}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function PetscObjectSetOptions(arg1, arg2)
    @ccall libpetsc.PetscObjectSetOptions(arg1::PetscObject, arg2::PetscOptions)::PetscErrorCode
end

function PetscObjectGetOptions(arg1, arg2)
    @ccall libpetsc.PetscObjectGetOptions(arg1::PetscObject, arg2::Ptr{PetscOptions})::PetscErrorCode
end

@enum PetscOptionType::UInt32 begin
    OPTION_INT = 0
    OPTION_BOOL = 1
    OPTION_REAL = 2
    OPTION_FLIST = 3
    OPTION_STRING = 4
    OPTION_REAL_ARRAY = 5
    OPTION_SCALAR_ARRAY = 6
    OPTION_HEAD = 7
    OPTION_INT_ARRAY = 8
    OPTION_ELIST = 9
    OPTION_BOOL_ARRAY = 10
    OPTION_STRING_ARRAY = 11
end

mutable struct __JL__n_PetscOptionItem
end

function Base.unsafe_load(x::Ptr{__JL__n_PetscOptionItem})
    unsafe_load(Ptr{_n_PetscOptionItem}(x))
end

function Base.getproperty(x::Ptr{__JL__n_PetscOptionItem}, f::Symbol)
    getproperty(Ptr{_n_PetscOptionItem}(x), f)
end

function Base.setproperty!(x::Ptr{__JL__n_PetscOptionItem}, f::Symbol, v)
    setproperty!(Ptr{_n_PetscOptionItem}(x), f, v)
end

const PetscOptionItem = Ptr{__JL__n_PetscOptionItem}

mutable struct _n_PetscOptionItem
    option::Ptr{Cchar}
    text::Ptr{Cchar}
    data::Ptr{Cvoid}
    flist::PetscFunctionList
    list::Ptr{Ptr{Cchar}}
    nlist::Cchar
    man::Ptr{Cchar}
    arraylength::Csize_t
    set::PetscBool
    type::PetscOptionType
    next::PetscOptionItem
    pman::Ptr{Cchar}
    edata::Ptr{Cvoid}
    _n_PetscOptionItem() = new()
end

mutable struct _p_PetscOptionItems
    count::PetscInt
    next::PetscOptionItem
    prefix::Ptr{Cchar}
    pprefix::Ptr{Cchar}
    title::Ptr{Cchar}
    comm::MPI_Comm
    printhelp::PetscBool
    changedmethod::PetscBool
    alreadyprinted::PetscBool
    object::PetscObject
    options::PetscOptions
    _p_PetscOptionItems() = new()
end

const PetscOptionItems = _p_PetscOptionItems

function PetscOptionsBegin_Private(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscOptionsBegin_Private(arg1::Ptr{PetscOptionItems}, arg2::MPI_Comm, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{Cchar})::PetscErrorCode
end

function PetscObjectOptionsBegin_Private(arg1, arg2)
    @ccall libpetsc.PetscObjectOptionsBegin_Private(arg1::Ptr{PetscOptionItems}, arg2::PetscObject)::PetscErrorCode
end

function PetscOptionsEnd_Private(arg1)
    @ccall libpetsc.PetscOptionsEnd_Private(arg1::Ptr{PetscOptionItems})::PetscErrorCode
end

function PetscOptionsHead(arg1, arg2)
    @ccall libpetsc.PetscOptionsHead(arg1::Ptr{PetscOptionItems}, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscOptionsEnum_Private(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscOptionsEnum_Private(arg1::Ptr{PetscOptionItems}, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{Ptr{Cchar}}, arg6::PetscEnum, arg7::Ptr{PetscEnum}, arg8::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsInt_Private(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.PetscOptionsInt_Private(arg1::Ptr{PetscOptionItems}, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::Ptr{PetscBool}, arg8::PetscInt, arg9::PetscInt)::PetscErrorCode
end

function PetscOptionsReal_Private(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscOptionsReal_Private(arg1::Ptr{PetscOptionItems}, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::PetscReal, arg6::Ptr{PetscReal}, arg7::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsScalar_Private(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscOptionsScalar_Private(arg1::Ptr{PetscOptionItems}, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::PetscScalar, arg6::Ptr{PetscScalar}, arg7::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsName_Private(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscOptionsName_Private(arg1::Ptr{PetscOptionItems}, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsString_Private(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscOptionsString_Private(arg1::Ptr{PetscOptionItems}, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{Cchar}, arg6::Ptr{Cchar}, arg7::Csize_t, arg8::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsBool_Private(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscOptionsBool_Private(arg1::Ptr{PetscOptionItems}, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::PetscBool, arg6::Ptr{PetscBool}, arg7::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsBoolGroupBegin_Private(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscOptionsBoolGroupBegin_Private(arg1::Ptr{PetscOptionItems}, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsBoolGroup_Private(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscOptionsBoolGroup_Private(arg1::Ptr{PetscOptionItems}, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsBoolGroupEnd_Private(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscOptionsBoolGroupEnd_Private(arg1::Ptr{PetscOptionItems}, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsFList_Private(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.PetscOptionsFList_Private(arg1::Ptr{PetscOptionItems}, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::PetscFunctionList, arg6::Ptr{Cchar}, arg7::Ptr{Cchar}, arg8::Csize_t, arg9::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsEList_Private(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.PetscOptionsEList_Private(arg1::Ptr{PetscOptionItems}, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{Ptr{Cchar}}, arg6::PetscInt, arg7::Ptr{Cchar}, arg8::Ptr{PetscInt}, arg9::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsRealArray_Private(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscOptionsRealArray_Private(arg1::Ptr{PetscOptionItems}, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{PetscReal}, arg6::Ptr{PetscInt}, arg7::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsScalarArray_Private(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscOptionsScalarArray_Private(arg1::Ptr{PetscOptionItems}, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{PetscScalar}, arg6::Ptr{PetscInt}, arg7::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsIntArray_Private(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscOptionsIntArray_Private(arg1::Ptr{PetscOptionItems}, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscInt}, arg7::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsStringArray_Private(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscOptionsStringArray_Private(arg1::Ptr{PetscOptionItems}, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{Ptr{Cchar}}, arg6::Ptr{PetscInt}, arg7::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsBoolArray_Private(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscOptionsBoolArray_Private(arg1::Ptr{PetscOptionItems}, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{PetscBool}, arg6::Ptr{PetscInt}, arg7::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsEnumArray_Private(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscOptionsEnumArray_Private(arg1::Ptr{PetscOptionItems}, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{Ptr{Cchar}}, arg6::Ptr{PetscEnum}, arg7::Ptr{PetscInt}, arg8::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsDeprecated_Private(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscOptionsDeprecated_Private(arg1::Ptr{PetscOptionItems}, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{Cchar})::PetscErrorCode
end

function PetscOptionsSAWsDestroy()
    @ccall libpetsc.PetscOptionsSAWsDestroy()::PetscErrorCode
end

function PetscObjectAddOptionsHandler(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscObjectAddOptionsHandler(arg1::PetscObject, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function PetscObjectProcessOptionsHandlers(arg1, arg2)
    @ccall libpetsc.PetscObjectProcessOptionsHandlers(arg1::Ptr{PetscOptionItems}, arg2::PetscObject)::PetscErrorCode
end

function PetscObjectDestroyOptionsHandlers(arg1)
    @ccall libpetsc.PetscObjectDestroyOptionsHandlers(arg1::PetscObject)::PetscErrorCode
end

function PetscMallocTraceSet(arg1, arg2, arg3)
    @ccall libpetsc.PetscMallocTraceSet(arg1::PetscViewer, arg2::PetscBool, arg3::PetscLogDouble)::PetscErrorCode
end

function PetscMallocTraceGet(arg1)
    @ccall libpetsc.PetscMallocTraceGet(arg1::Ptr{PetscBool})::PetscErrorCode
end

function PetscObjectsListGetGlobalNumbering(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscObjectsListGetGlobalNumbering(arg1::MPI_Comm, arg2::PetscInt, arg3::Ptr{PetscObject}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt})::PetscErrorCode
end

function PetscMemoryView(arg1, arg2)
    @ccall libpetsc.PetscMemoryView(arg1::PetscViewer, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscObjectPrintClassNamePrefixType(arg1, arg2)
    @ccall libpetsc.PetscObjectPrintClassNamePrefixType(arg1::PetscObject, arg2::PetscViewer)::PetscErrorCode
end

function PetscObjectView(arg1, arg2)
    @ccall libpetsc.PetscObjectView(arg1::PetscObject, arg2::PetscViewer)::PetscErrorCode
end

function PetscObjectQueryFunction_Private(arg1, arg2, arg3)
    @ccall libpetsc.PetscObjectQueryFunction_Private(arg1::PetscObject, arg2::Ptr{Cchar}, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscObjectSetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.PetscObjectSetOptionsPrefix(arg1::PetscObject, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscObjectAppendOptionsPrefix(arg1, arg2)
    @ccall libpetsc.PetscObjectAppendOptionsPrefix(arg1::PetscObject, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscObjectPrependOptionsPrefix(arg1, arg2)
    @ccall libpetsc.PetscObjectPrependOptionsPrefix(arg1::PetscObject, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscObjectGetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.PetscObjectGetOptionsPrefix(arg1::PetscObject, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscObjectChangeTypeName(arg1, arg2)
    @ccall libpetsc.PetscObjectChangeTypeName(arg1::PetscObject, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscObjectRegisterDestroy(arg1)
    @ccall libpetsc.PetscObjectRegisterDestroy(arg1::PetscObject)::PetscErrorCode
end

function PetscObjectRegisterDestroyAll()
    @ccall libpetsc.PetscObjectRegisterDestroyAll()::PetscErrorCode
end

function PetscObjectViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.PetscObjectViewFromOptions(arg1::PetscObject, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function PetscObjectName(arg1)
    @ccall libpetsc.PetscObjectName(arg1::PetscObject)::PetscErrorCode
end

function PetscObjectTypeCompare(arg1, arg2, arg3)
    @ccall libpetsc.PetscObjectTypeCompare(arg1::PetscObject, arg2::Ptr{Cchar}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function PetscObjectBaseTypeCompare(arg1, arg2, arg3)
    @ccall libpetsc.PetscObjectBaseTypeCompare(arg1::PetscObject, arg2::Ptr{Cchar}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function PetscRegisterFinalize(arg1)
    @ccall libpetsc.PetscRegisterFinalize(arg1::Ptr{Cvoid})::PetscErrorCode
end

function PetscRegisterFinalizeAll()
    @ccall libpetsc.PetscRegisterFinalizeAll()::PetscErrorCode
end

function PetscDLOpen(arg1, arg2, arg3)
    @ccall libpetsc.PetscDLOpen(arg1::Ptr{Cchar}, arg2::PetscDLMode, arg3::Ptr{PetscDLHandle})::PetscErrorCode
end

function PetscDLClose(arg1)
    @ccall libpetsc.PetscDLClose(arg1::Ptr{PetscDLHandle})::PetscErrorCode
end

function PetscDLSym(arg1, arg2, arg3)
    @ccall libpetsc.PetscDLSym(arg1::PetscDLHandle, arg2::Ptr{Cchar}, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscDLAddr(arg1, arg2)
    @ccall libpetsc.PetscDLAddr(arg1::Ptr{Cvoid}, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscObjectsDump(arg1, arg2)
    @ccall libpetsc.PetscObjectsDump(arg1::Ptr{Libc.FILE}, arg2::PetscBool)::PetscErrorCode
end

function PetscObjectListDestroy(arg1)
    @ccall libpetsc.PetscObjectListDestroy(arg1::Ptr{PetscObjectList})::PetscErrorCode
end

function PetscObjectListFind(arg1, arg2, arg3)
    @ccall libpetsc.PetscObjectListFind(arg1::PetscObjectList, arg2::Ptr{Cchar}, arg3::Ptr{PetscObject})::PetscErrorCode
end

function PetscObjectListReverseFind(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscObjectListReverseFind(arg1::PetscObjectList, arg2::PetscObject, arg3::Ptr{Ptr{Cchar}}, arg4::Ptr{PetscBool})::PetscErrorCode
end

function PetscObjectListAdd(arg1, arg2, arg3)
    @ccall libpetsc.PetscObjectListAdd(arg1::Ptr{PetscObjectList}, arg2::Ptr{Cchar}, arg3::PetscObject)::PetscErrorCode
end

function PetscObjectListRemoveReference(arg1, arg2)
    @ccall libpetsc.PetscObjectListRemoveReference(arg1::Ptr{PetscObjectList}, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscObjectListDuplicate(arg1, arg2)
    @ccall libpetsc.PetscObjectListDuplicate(arg1::PetscObjectList, arg2::Ptr{PetscObjectList})::PetscErrorCode
end

function PetscFunctionListAdd_Private(arg1, arg2, arg3)
    @ccall libpetsc.PetscFunctionListAdd_Private(arg1::Ptr{PetscFunctionList}, arg2::Ptr{Cchar}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function PetscFunctionListDestroy(arg1)
    @ccall libpetsc.PetscFunctionListDestroy(arg1::Ptr{PetscFunctionList})::PetscErrorCode
end

function PetscFunctionListFind_Private(arg1, arg2, arg3)
    @ccall libpetsc.PetscFunctionListFind_Private(arg1::PetscFunctionList, arg2::Ptr{Cchar}, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscFunctionListPrintTypes(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.PetscFunctionListPrintTypes(arg1::MPI_Comm, arg2::Ptr{Libc.FILE}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{Cchar}, arg6::Ptr{Cchar}, arg7::PetscFunctionList, arg8::Ptr{Cchar}, arg9::Ptr{Cchar})::PetscErrorCode
end

function PetscFunctionListDuplicate(arg1, arg2)
    @ccall libpetsc.PetscFunctionListDuplicate(arg1::PetscFunctionList, arg2::Ptr{PetscFunctionList})::PetscErrorCode
end

function PetscFunctionListView(arg1, arg2)
    @ccall libpetsc.PetscFunctionListView(arg1::PetscFunctionList, arg2::PetscViewer)::PetscErrorCode
end

function PetscFunctionListGet(arg1, arg2, arg3)
    @ccall libpetsc.PetscFunctionListGet(arg1::PetscFunctionList, arg2::Ptr{Ptr{Ptr{Cchar}}}, arg3::Ptr{Cint})::PetscErrorCode
end

function PetscDLLibraryAppend(arg1, arg2, arg3)
    @ccall libpetsc.PetscDLLibraryAppend(arg1::MPI_Comm, arg2::Ptr{PetscDLLibrary}, arg3::Ptr{Cchar})::PetscErrorCode
end

function PetscDLLibraryPrepend(arg1, arg2, arg3)
    @ccall libpetsc.PetscDLLibraryPrepend(arg1::MPI_Comm, arg2::Ptr{PetscDLLibrary}, arg3::Ptr{Cchar})::PetscErrorCode
end

function PetscDLLibrarySym(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDLLibrarySym(arg1::MPI_Comm, arg2::Ptr{PetscDLLibrary}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscDLLibraryPrintPath(arg1)
    @ccall libpetsc.PetscDLLibraryPrintPath(arg1::PetscDLLibrary)::PetscErrorCode
end

function PetscDLLibraryRetrieve(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDLLibraryRetrieve(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Csize_t, arg5::Ptr{PetscBool})::PetscErrorCode
end

function PetscDLLibraryOpen(arg1, arg2, arg3)
    @ccall libpetsc.PetscDLLibraryOpen(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{PetscDLLibrary})::PetscErrorCode
end

function PetscDLLibraryClose(arg1)
    @ccall libpetsc.PetscDLLibraryClose(arg1::PetscDLLibrary)::PetscErrorCode
end

function PetscSplitOwnership(arg1, arg2, arg3)
    @ccall libpetsc.PetscSplitOwnership(arg1::MPI_Comm, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscSplitOwnershipBlock(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSplitOwnershipBlock(arg1::MPI_Comm, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function PetscSplitOwnershipEqual(arg1, arg2, arg3)
    @ccall libpetsc.PetscSplitOwnershipEqual(arg1::MPI_Comm, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscSequentialPhaseBegin(arg1, arg2)
    @ccall libpetsc.PetscSequentialPhaseBegin(arg1::MPI_Comm, arg2::PetscMPIInt)::PetscErrorCode
end

function PetscSequentialPhaseEnd(arg1, arg2)
    @ccall libpetsc.PetscSequentialPhaseEnd(arg1::MPI_Comm, arg2::PetscMPIInt)::PetscErrorCode
end

function PetscBarrier(arg1)
    @ccall libpetsc.PetscBarrier(arg1::PetscObject)::PetscErrorCode
end

function PetscMPIDump(arg1)
    @ccall libpetsc.PetscMPIDump(arg1::Ptr{Libc.FILE})::PetscErrorCode
end

function PetscGlobalMinMaxInt(arg1, arg2, arg3)
    @ccall libpetsc.PetscGlobalMinMaxInt(arg1::MPI_Comm, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscGlobalMinMaxReal(arg1, arg2, arg3)
    @ccall libpetsc.PetscGlobalMinMaxReal(arg1::MPI_Comm, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal})::PetscErrorCode
end

function PetscGetCPUTime(arg1)
    @ccall libpetsc.PetscGetCPUTime(arg1::Ptr{PetscLogDouble})::PetscErrorCode
end

function PetscTime(v)
    @ccall libpetsc.PetscTime(v::Ptr{PetscLogDouble})::PetscErrorCode
end

function PetscTimeSubtract(v)
    @ccall libpetsc.PetscTimeSubtract(v::Ptr{PetscLogDouble})::PetscErrorCode
end

function PetscTimeAdd(v)
    @ccall libpetsc.PetscTimeAdd(v::Ptr{PetscLogDouble})::PetscErrorCode
end

@enum PetscInfoCommFlag::Int32 begin
    PETSC_INFO_COMM_ALL = -1
    PETSC_INFO_COMM_NO_SELF = 0
    PETSC_INFO_COMM_ONLY_SELF = 1
end

function PetscInfoDeactivateClass(arg1)
    @ccall libpetsc.PetscInfoDeactivateClass(arg1::PetscClassId)::PetscErrorCode
end

function PetscInfoActivateClass(arg1)
    @ccall libpetsc.PetscInfoActivateClass(arg1::PetscClassId)::PetscErrorCode
end

function PetscInfoEnabled(arg1, arg2)
    @ccall libpetsc.PetscInfoEnabled(arg1::PetscClassId, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscInfoAllow(arg1)
    @ccall libpetsc.PetscInfoAllow(arg1::PetscBool)::PetscErrorCode
end

function PetscInfoSetFile(arg1, arg2)
    @ccall libpetsc.PetscInfoSetFile(arg1::Ptr{Cchar}, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscInfoGetFile(arg1, arg2)
    @ccall libpetsc.PetscInfoGetFile(arg1::Ptr{Ptr{Cchar}}, arg2::Ptr{Ptr{Libc.FILE}})::PetscErrorCode
end

function PetscInfoSetClasses(arg1, arg2, arg3)
    @ccall libpetsc.PetscInfoSetClasses(arg1::PetscBool, arg2::PetscInt, arg3::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscInfoGetClass(arg1, arg2)
    @ccall libpetsc.PetscInfoGetClass(arg1::Ptr{Cchar}, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscInfoGetInfo(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscInfoGetInfo(arg1::Ptr{PetscBool}, arg2::Ptr{PetscBool}, arg3::Ptr{PetscBool}, arg4::Ptr{PetscBool}, arg5::Ptr{PetscInfoCommFlag})::PetscErrorCode
end

function PetscInfoProcessClass(arg1, arg2, arg3)
    @ccall libpetsc.PetscInfoProcessClass(arg1::Ptr{Cchar}, arg2::PetscInt, arg3::Ptr{PetscClassId})::PetscErrorCode
end

function PetscInfoSetFilterCommSelf(arg1)
    @ccall libpetsc.PetscInfoSetFilterCommSelf(arg1::PetscInfoCommFlag)::PetscErrorCode
end

function PetscInfoSetFromOptions(arg1)
    @ccall libpetsc.PetscInfoSetFromOptions(arg1::PetscOptions)::PetscErrorCode
end

function PetscInfoDestroy()
    @ccall libpetsc.PetscInfoDestroy()::PetscErrorCode
end

const PetscLogEvent = Cint

const PetscLogStage = Cint

mutable struct _n_PetscIntStack end

const PetscIntStack = Ptr{_n_PetscIntStack}

mutable struct PetscClassRegInfo
    name::Ptr{Cchar}
    classid::PetscClassId
    PetscClassRegInfo() = new()
end

mutable struct PetscClassPerfInfo
    id::PetscClassId
    creations::Cint
    destructions::Cint
    mem::PetscLogDouble
    descMem::PetscLogDouble
    PetscClassPerfInfo() = new()
end

struct _n_PetscClassRegLog
    numClasses::Cint
    maxClasses::Cint
    classInfo::Ptr{PetscClassRegInfo}
end

const PetscClassRegLog = Ptr{_n_PetscClassRegLog}

struct _n_PetscClassPerfLog
    numClasses::Cint
    maxClasses::Cint
    classInfo::Ptr{PetscClassPerfInfo}
end

const PetscClassPerfLog = Ptr{_n_PetscClassPerfLog}

mutable struct PetscEventRegInfo
    name::Ptr{Cchar}
    classid::PetscClassId
    collective::PetscBool
    PetscEventRegInfo() = new()
end

mutable struct PetscEventPerfInfo
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
    syncTime::PetscLogDouble
    dof::NTuple{8, PetscLogDouble}
    errors::NTuple{8, PetscLogDouble}
    numMessages::PetscLogDouble
    messageLength::PetscLogDouble
    numReductions::PetscLogDouble
    memIncrease::PetscLogDouble
    mallocIncrease::PetscLogDouble
    mallocSpace::PetscLogDouble
    mallocIncreaseEvent::PetscLogDouble
    PetscEventPerfInfo() = new()
end

struct _n_PetscEventRegLog
    numEvents::Cint
    maxEvents::Cint
    eventInfo::Ptr{PetscEventRegInfo}
end

const PetscEventRegLog = Ptr{_n_PetscEventRegLog}

struct _n_PetscEventPerfLog
    numEvents::Cint
    maxEvents::Cint
    eventInfo::Ptr{PetscEventPerfInfo}
end

const PetscEventPerfLog = Ptr{_n_PetscEventPerfLog}

struct _PetscStageInfo
    name::Ptr{Cchar}
    used::PetscBool
    perfInfo::PetscEventPerfInfo
    eventLog::PetscEventPerfLog
    classLog::PetscClassPerfLog
end

const PetscStageInfo = _PetscStageInfo

mutable struct _n_PetscStageLog
    numStages::Cint
    maxStages::Cint
    stack::PetscIntStack
    curStage::Cint
    stageInfo::Ptr{PetscStageInfo}
    eventLog::PetscEventRegLog
    classLog::PetscClassRegLog
    _n_PetscStageLog() = new()
end

const PetscStageLog = Ptr{_n_PetscStageLog}

function PetscLogObjectParent(arg1, arg2)
    @ccall libpetsc.PetscLogObjectParent(arg1::PetscObject, arg2::PetscObject)::PetscErrorCode
end

function PetscLogObjectMemory(arg1, arg2)
    @ccall libpetsc.PetscLogObjectMemory(arg1::PetscObject, arg2::PetscLogDouble)::PetscErrorCode
end

function PetscLogGetStageLog(arg1)
    @ccall libpetsc.PetscLogGetStageLog(arg1::Ptr{PetscStageLog})::PetscErrorCode
end

function PetscStageLogGetCurrent(arg1, arg2)
    @ccall libpetsc.PetscStageLogGetCurrent(arg1::PetscStageLog, arg2::Ptr{Cint})::PetscErrorCode
end

function PetscStageLogGetEventPerfLog(arg1, arg2, arg3)
    @ccall libpetsc.PetscStageLogGetEventPerfLog(arg1::PetscStageLog, arg2::Cint, arg3::Ptr{PetscEventPerfLog})::PetscErrorCode
end

function PetscLogFlops(n)
    @ccall libpetsc.PetscLogFlops(n::PetscLogDouble)::PetscErrorCode
end

function PetscGetFlops(arg1)
    @ccall libpetsc.PetscGetFlops(arg1::Ptr{PetscLogDouble})::PetscErrorCode
end

function PetscLogDefaultBegin()
    @ccall libpetsc.PetscLogDefaultBegin()::PetscErrorCode
end

function PetscLogAllBegin()
    @ccall libpetsc.PetscLogAllBegin()::PetscErrorCode
end

function PetscLogNestedBegin()
    @ccall libpetsc.PetscLogNestedBegin()::PetscErrorCode
end

function PetscLogTraceBegin(arg1)
    @ccall libpetsc.PetscLogTraceBegin(arg1::Ptr{Libc.FILE})::PetscErrorCode
end

function PetscLogActions(arg1)
    @ccall libpetsc.PetscLogActions(arg1::PetscBool)::PetscErrorCode
end

function PetscLogObjects(arg1)
    @ccall libpetsc.PetscLogObjects(arg1::PetscBool)::PetscErrorCode
end

function PetscLogSetThreshold(arg1, arg2)
    @ccall libpetsc.PetscLogSetThreshold(arg1::PetscLogDouble, arg2::Ptr{PetscLogDouble})::PetscErrorCode
end

function PetscLogSet(arg1, arg2)
    @ccall libpetsc.PetscLogSet(arg1::Ptr{Cvoid}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscLogView(arg1)
    @ccall libpetsc.PetscLogView(arg1::PetscViewer)::PetscErrorCode
end

function PetscLogViewFromOptions()
    @ccall libpetsc.PetscLogViewFromOptions()::PetscErrorCode
end

function PetscLogDump(arg1)
    @ccall libpetsc.PetscLogDump(arg1::Ptr{Cchar})::PetscErrorCode
end

function PetscLogStageRegister(arg1, arg2)
    @ccall libpetsc.PetscLogStageRegister(arg1::Ptr{Cchar}, arg2::Ptr{PetscLogStage})::PetscErrorCode
end

function PetscLogStagePush(arg1)
    @ccall libpetsc.PetscLogStagePush(arg1::PetscLogStage)::PetscErrorCode
end

function PetscLogStagePop()
    @ccall libpetsc.PetscLogStagePop()::PetscErrorCode
end

function PetscLogStageSetActive(arg1, arg2)
    @ccall libpetsc.PetscLogStageSetActive(arg1::PetscLogStage, arg2::PetscBool)::PetscErrorCode
end

function PetscLogStageGetActive(arg1, arg2)
    @ccall libpetsc.PetscLogStageGetActive(arg1::PetscLogStage, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscLogStageSetVisible(arg1, arg2)
    @ccall libpetsc.PetscLogStageSetVisible(arg1::PetscLogStage, arg2::PetscBool)::PetscErrorCode
end

function PetscLogStageGetVisible(arg1, arg2)
    @ccall libpetsc.PetscLogStageGetVisible(arg1::PetscLogStage, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscLogStageGetId(arg1, arg2)
    @ccall libpetsc.PetscLogStageGetId(arg1::Ptr{Cchar}, arg2::Ptr{PetscLogStage})::PetscErrorCode
end

function PetscLogEventRegister(arg1, arg2, arg3)
    @ccall libpetsc.PetscLogEventRegister(arg1::Ptr{Cchar}, arg2::PetscClassId, arg3::Ptr{PetscLogEvent})::PetscErrorCode
end

function PetscLogEventSetCollective(arg1, arg2)
    @ccall libpetsc.PetscLogEventSetCollective(arg1::PetscLogEvent, arg2::PetscBool)::PetscErrorCode
end

function PetscLogEventIncludeClass(arg1)
    @ccall libpetsc.PetscLogEventIncludeClass(arg1::PetscClassId)::PetscErrorCode
end

function PetscLogEventExcludeClass(arg1)
    @ccall libpetsc.PetscLogEventExcludeClass(arg1::PetscClassId)::PetscErrorCode
end

function PetscLogEventActivate(arg1)
    @ccall libpetsc.PetscLogEventActivate(arg1::PetscLogEvent)::PetscErrorCode
end

function PetscLogEventDeactivate(arg1)
    @ccall libpetsc.PetscLogEventDeactivate(arg1::PetscLogEvent)::PetscErrorCode
end

function PetscLogEventDeactivatePush(arg1)
    @ccall libpetsc.PetscLogEventDeactivatePush(arg1::PetscLogEvent)::PetscErrorCode
end

function PetscLogEventDeactivatePop(arg1)
    @ccall libpetsc.PetscLogEventDeactivatePop(arg1::PetscLogEvent)::PetscErrorCode
end

function PetscLogEventSetActiveAll(arg1, arg2)
    @ccall libpetsc.PetscLogEventSetActiveAll(arg1::PetscLogEvent, arg2::PetscBool)::PetscErrorCode
end

function PetscLogEventActivateClass(arg1)
    @ccall libpetsc.PetscLogEventActivateClass(arg1::PetscClassId)::PetscErrorCode
end

function PetscLogEventDeactivateClass(arg1)
    @ccall libpetsc.PetscLogEventDeactivateClass(arg1::PetscClassId)::PetscErrorCode
end

function PetscLogEventGetId(arg1, arg2)
    @ccall libpetsc.PetscLogEventGetId(arg1::Ptr{Cchar}, arg2::Ptr{PetscLogEvent})::PetscErrorCode
end

function PetscLogEventGetPerfInfo(arg1, arg2, arg3)
    @ccall libpetsc.PetscLogEventGetPerfInfo(arg1::Cint, arg2::PetscLogEvent, arg3::Ptr{PetscEventPerfInfo})::PetscErrorCode
end

function PetscLogEventSetDof(arg1, arg2, arg3)
    @ccall libpetsc.PetscLogEventSetDof(arg1::PetscLogEvent, arg2::PetscInt, arg3::PetscLogDouble)::PetscErrorCode
end

function PetscLogEventSetError(arg1, arg2, arg3)
    @ccall libpetsc.PetscLogEventSetError(arg1::PetscLogEvent, arg2::PetscInt, arg3::PetscLogDouble)::PetscErrorCode
end

function PetscLogEventSynchronize(arg1, arg2)
    @ccall libpetsc.PetscLogEventSynchronize(arg1::PetscLogEvent, arg2::MPI_Comm)::PetscErrorCode
end

function PetscLogEventGetFlops(arg1, arg2)
    @ccall libpetsc.PetscLogEventGetFlops(arg1::PetscLogEvent, arg2::Ptr{PetscLogDouble})::PetscErrorCode
end

function PetscLogEventZeroFlops(arg1)
    @ccall libpetsc.PetscLogEventZeroFlops(arg1::PetscLogEvent)::PetscErrorCode
end

function PetscMPITypeSize(count, type, length)
    @ccall libpetsc.PetscMPITypeSize(count::PetscInt, type::MPI_Datatype, length::Ptr{PetscLogDouble})::PetscErrorCode
end

function PetscMPITypeSizeComm(comm, counts, type, length)
    @ccall libpetsc.PetscMPITypeSizeComm(comm::MPI_Comm, counts::Ptr{PetscMPIInt}, type::MPI_Datatype, length::Ptr{PetscLogDouble})::PetscErrorCode
end

function PetscMPITypeSizeCount(n, counts, type, length)
    @ccall libpetsc.PetscMPITypeSizeCount(n::PetscInt, counts::Ptr{PetscMPIInt}, type::MPI_Datatype, length::Ptr{PetscLogDouble})::PetscErrorCode
end

function PetscMPIParallelComm(comm)
    @ccall libpetsc.PetscMPIParallelComm(comm::MPI_Comm)::Cint
end

function PetscFixFilename(arg1, arg2)
    @ccall libpetsc.PetscFixFilename(arg1::Ptr{Cchar}, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscFOpen(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscFOpen(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Ptr{Libc.FILE}})::PetscErrorCode
end

function PetscFClose(arg1, arg2)
    @ccall libpetsc.PetscFClose(arg1::MPI_Comm, arg2::Ptr{Libc.FILE})::PetscErrorCode
end

function PetscFormatRealArray(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscFormatRealArray(arg1::Ptr{Cchar}, arg2::Csize_t, arg3::Ptr{Cchar}, arg4::PetscInt, arg5::Ptr{PetscReal})::PetscErrorCode
end

function PetscFormatConvertGetSize(arg1, arg2)
    @ccall libpetsc.PetscFormatConvertGetSize(arg1::Ptr{Cchar}, arg2::Ptr{Csize_t})::PetscErrorCode
end

function PetscFormatConvert(arg1, arg2)
    @ccall libpetsc.PetscFormatConvert(arg1::Ptr{Cchar}, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscPOpen(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscPOpen(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{Ptr{Libc.FILE}})::PetscErrorCode
end

function PetscPClose(arg1, arg2)
    @ccall libpetsc.PetscPClose(arg1::MPI_Comm, arg2::Ptr{Libc.FILE})::PetscErrorCode
end

function PetscPOpenSetMachine(arg1)
    @ccall libpetsc.PetscPOpenSetMachine(arg1::Ptr{Cchar})::PetscErrorCode
end

function PetscSynchronizedFlush(arg1, arg2)
    @ccall libpetsc.PetscSynchronizedFlush(arg1::MPI_Comm, arg2::Ptr{Libc.FILE})::PetscErrorCode
end

function PetscSynchronizedFGets(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSynchronizedFGets(arg1::MPI_Comm, arg2::Ptr{Libc.FILE}, arg3::Csize_t, arg4::Ptr{Cchar})::PetscErrorCode
end

function PetscStartMatlab(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscStartMatlab(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Ptr{Libc.FILE}})::PetscErrorCode
end

function PetscStartJava(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscStartJava(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Ptr{Libc.FILE}})::PetscErrorCode
end

function PetscGetPetscDir(arg1)
    @ccall libpetsc.PetscGetPetscDir(arg1::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscContainerGetPointer(arg1, arg2)
    @ccall libpetsc.PetscContainerGetPointer(arg1::PetscContainer, arg2::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscContainerSetPointer(arg1, arg2)
    @ccall libpetsc.PetscContainerSetPointer(arg1::PetscContainer, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscContainerDestroy(arg1)
    @ccall libpetsc.PetscContainerDestroy(arg1::Ptr{PetscContainer})::PetscErrorCode
end

function PetscContainerCreate(arg1, arg2)
    @ccall libpetsc.PetscContainerCreate(arg1::MPI_Comm, arg2::Ptr{PetscContainer})::PetscErrorCode
end

function PetscContainerSetUserDestroy(arg1, arg2)
    @ccall libpetsc.PetscContainerSetUserDestroy(arg1::PetscContainer, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscContainerUserDestroyDefault(arg1)
    @ccall libpetsc.PetscContainerUserDestroyDefault(arg1::Ptr{Cvoid})::PetscErrorCode
end

function PetscIntView(arg1, arg2, arg3)
    @ccall libpetsc.PetscIntView(arg1::PetscInt, arg2::Ptr{PetscInt}, arg3::PetscViewer)::PetscErrorCode
end

function PetscRealView(arg1, arg2, arg3)
    @ccall libpetsc.PetscRealView(arg1::PetscInt, arg2::Ptr{PetscReal}, arg3::PetscViewer)::PetscErrorCode
end

function PetscScalarView(arg1, arg2, arg3)
    @ccall libpetsc.PetscScalarView(arg1::PetscInt, arg2::Ptr{PetscScalar}, arg3::PetscViewer)::PetscErrorCode
end

function PetscMemmove(a, b, n)
    @ccall libpetsc.PetscMemmove(a::Ptr{Cvoid}, b::Ptr{Cvoid}, n::Csize_t)::PetscErrorCode
end

function PetscMemcpy(a, b, n)
    @ccall libpetsc.PetscMemcpy(a::Ptr{Cvoid}, b::Ptr{Cvoid}, n::Csize_t)::PetscErrorCode
end

function PetscMemzero(a, n)
    @ccall libpetsc.PetscMemzero(a::Ptr{Cvoid}, n::Csize_t)::PetscErrorCode
end

function MPIU_File_write_all(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MPIU_File_write_all(arg1::MPI_File, arg2::Ptr{Cvoid}, arg3::PetscMPIInt, arg4::MPI_Datatype, arg5::Ptr{MPI_Status})::PetscErrorCode
end

function MPIU_File_read_all(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MPIU_File_read_all(arg1::MPI_File, arg2::Ptr{Cvoid}, arg3::PetscMPIInt, arg4::MPI_Datatype, arg5::Ptr{MPI_Status})::PetscErrorCode
end

function MPIU_File_write_at(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MPIU_File_write_at(arg1::MPI_File, arg2::MPI_Offset, arg3::Ptr{Cvoid}, arg4::PetscMPIInt, arg5::MPI_Datatype, arg6::Ptr{MPI_Status})::PetscErrorCode
end

function MPIU_File_read_at(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MPIU_File_read_at(arg1::MPI_File, arg2::MPI_Offset, arg3::Ptr{Cvoid}, arg4::PetscMPIInt, arg5::MPI_Datatype, arg6::Ptr{MPI_Status})::PetscErrorCode
end

function MPIU_File_write_at_all(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MPIU_File_write_at_all(arg1::MPI_File, arg2::MPI_Offset, arg3::Ptr{Cvoid}, arg4::PetscMPIInt, arg5::MPI_Datatype, arg6::Ptr{MPI_Status})::PetscErrorCode
end

function MPIU_File_read_at_all(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MPIU_File_read_at_all(arg1::MPI_File, arg2::MPI_Offset, arg3::Ptr{Cvoid}, arg4::PetscMPIInt, arg5::MPI_Datatype, arg6::Ptr{MPI_Status})::PetscErrorCode
end

function PetscIntCast(a, b)
    @ccall libpetsc.PetscIntCast(a::PetscInt64, b::Ptr{PetscInt})::PetscErrorCode
end

function PetscBLASIntCast(a, b)
    @ccall libpetsc.PetscBLASIntCast(a::PetscInt, b::Ptr{PetscBLASInt})::PetscErrorCode
end

function PetscMPIIntCast(a, b)
    @ccall libpetsc.PetscMPIIntCast(a::PetscInt, b::Ptr{PetscMPIInt})::PetscErrorCode
end

function PetscRealIntMultTruncate(a, b)
    @ccall libpetsc.PetscRealIntMultTruncate(a::PetscReal, b::PetscInt)::PetscInt
end

function PetscIntMultTruncate(a, b)
    @ccall libpetsc.PetscIntMultTruncate(a::PetscInt, b::PetscInt)::PetscInt
end

function PetscIntSumTruncate(a, b)
    @ccall libpetsc.PetscIntSumTruncate(a::PetscInt, b::PetscInt)::PetscInt
end

function PetscIntMultError(a, b, result)
    @ccall libpetsc.PetscIntMultError(a::PetscInt, b::PetscInt, result::Ptr{PetscInt})::PetscErrorCode
end

function PetscIntSumError(a, b, result)
    @ccall libpetsc.PetscIntSumError(a::PetscInt, b::PetscInt, result::Ptr{PetscInt})::PetscErrorCode
end

function PetscGetArchType(arg1, arg2)
    @ccall libpetsc.PetscGetArchType(arg1::Ptr{Cchar}, arg2::Csize_t)::PetscErrorCode
end

function PetscGetHostName(arg1, arg2)
    @ccall libpetsc.PetscGetHostName(arg1::Ptr{Cchar}, arg2::Csize_t)::PetscErrorCode
end

function PetscGetUserName(arg1, arg2)
    @ccall libpetsc.PetscGetUserName(arg1::Ptr{Cchar}, arg2::Csize_t)::PetscErrorCode
end

function PetscGetProgramName(arg1, arg2)
    @ccall libpetsc.PetscGetProgramName(arg1::Ptr{Cchar}, arg2::Csize_t)::PetscErrorCode
end

function PetscSetProgramName(arg1)
    @ccall libpetsc.PetscSetProgramName(arg1::Ptr{Cchar})::PetscErrorCode
end

function PetscGetDate(arg1, arg2)
    @ccall libpetsc.PetscGetDate(arg1::Ptr{Cchar}, arg2::Csize_t)::PetscErrorCode
end

function PetscGetVersion(arg1, arg2)
    @ccall libpetsc.PetscGetVersion(arg1::Ptr{Cchar}, arg2::Csize_t)::PetscErrorCode
end

function PetscGetVersionNumber(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscGetVersionNumber(arg1::Ptr{PetscInt}, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function PetscSortedInt(arg1, arg2, arg3)
    @ccall libpetsc.PetscSortedInt(arg1::PetscInt, arg2::Ptr{PetscInt}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function PetscSortedMPIInt(arg1, arg2, arg3)
    @ccall libpetsc.PetscSortedMPIInt(arg1::PetscInt, arg2::Ptr{PetscMPIInt}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function PetscSortedReal(arg1, arg2, arg3)
    @ccall libpetsc.PetscSortedReal(arg1::PetscInt, arg2::Ptr{PetscReal}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function PetscSortInt(arg1, arg2)
    @ccall libpetsc.PetscSortInt(arg1::PetscInt, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscSortReverseInt(arg1, arg2)
    @ccall libpetsc.PetscSortReverseInt(arg1::PetscInt, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscSortedRemoveDupsInt(arg1, arg2)
    @ccall libpetsc.PetscSortedRemoveDupsInt(arg1::Ptr{PetscInt}, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscSortRemoveDupsInt(arg1, arg2)
    @ccall libpetsc.PetscSortRemoveDupsInt(arg1::Ptr{PetscInt}, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscCheckDupsInt(arg1, arg2, arg3)
    @ccall libpetsc.PetscCheckDupsInt(arg1::PetscInt, arg2::Ptr{PetscInt}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function PetscFindInt(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscFindInt(arg1::PetscInt, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function PetscFindMPIInt(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscFindMPIInt(arg1::PetscMPIInt, arg2::PetscInt, arg3::Ptr{PetscMPIInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function PetscSortIntWithPermutation(arg1, arg2, arg3)
    @ccall libpetsc.PetscSortIntWithPermutation(arg1::PetscInt, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscSortStrWithPermutation(arg1, arg2, arg3)
    @ccall libpetsc.PetscSortStrWithPermutation(arg1::PetscInt, arg2::Ptr{Ptr{Cchar}}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscSortIntWithArray(arg1, arg2, arg3)
    @ccall libpetsc.PetscSortIntWithArray(arg1::PetscInt, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscSortIntWithArrayPair(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSortIntWithArrayPair(arg1::PetscInt, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function PetscSortMPIInt(arg1, arg2)
    @ccall libpetsc.PetscSortMPIInt(arg1::PetscInt, arg2::Ptr{PetscMPIInt})::PetscErrorCode
end

function PetscSortRemoveDupsMPIInt(arg1, arg2)
    @ccall libpetsc.PetscSortRemoveDupsMPIInt(arg1::Ptr{PetscInt}, arg2::Ptr{PetscMPIInt})::PetscErrorCode
end

function PetscSortMPIIntWithArray(arg1, arg2, arg3)
    @ccall libpetsc.PetscSortMPIIntWithArray(arg1::PetscMPIInt, arg2::Ptr{PetscMPIInt}, arg3::Ptr{PetscMPIInt})::PetscErrorCode
end

function PetscSortMPIIntWithIntArray(arg1, arg2, arg3)
    @ccall libpetsc.PetscSortMPIIntWithIntArray(arg1::PetscMPIInt, arg2::Ptr{PetscMPIInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscSortIntWithScalarArray(arg1, arg2, arg3)
    @ccall libpetsc.PetscSortIntWithScalarArray(arg1::PetscInt, arg2::Ptr{PetscInt}, arg3::Ptr{PetscScalar})::PetscErrorCode
end

function PetscSortIntWithDataArray(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscSortIntWithDataArray(arg1::PetscInt, arg2::Ptr{PetscInt}, arg3::Ptr{Cvoid}, arg4::Csize_t, arg5::Ptr{Cvoid})::PetscErrorCode
end

function PetscSortReal(arg1, arg2)
    @ccall libpetsc.PetscSortReal(arg1::PetscInt, arg2::Ptr{PetscReal})::PetscErrorCode
end

function PetscSortRealWithArrayInt(arg1, arg2, arg3)
    @ccall libpetsc.PetscSortRealWithArrayInt(arg1::PetscInt, arg2::Ptr{PetscReal}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscSortRealWithPermutation(arg1, arg2, arg3)
    @ccall libpetsc.PetscSortRealWithPermutation(arg1::PetscInt, arg2::Ptr{PetscReal}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscSortRemoveDupsReal(arg1, arg2)
    @ccall libpetsc.PetscSortRemoveDupsReal(arg1::Ptr{PetscInt}, arg2::Ptr{PetscReal})::PetscErrorCode
end

function PetscFindReal(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscFindReal(arg1::PetscReal, arg2::PetscInt, arg3::Ptr{PetscReal}, arg4::PetscReal, arg5::Ptr{PetscInt})::PetscErrorCode
end

function PetscSortSplit(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSortSplit(arg1::PetscInt, arg2::PetscInt, arg3::Ptr{PetscScalar}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function PetscSortSplitReal(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSortSplitReal(arg1::PetscInt, arg2::PetscInt, arg3::Ptr{PetscReal}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function PetscProcessTree(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscProcessTree(arg1::PetscInt, arg2::Ptr{PetscBool}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{Ptr{PetscInt}}, arg6::Ptr{Ptr{PetscInt}}, arg7::Ptr{Ptr{PetscInt}}, arg8::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function PetscMergeIntArrayPair(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.PetscMergeIntArrayPair(arg1::PetscInt, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{PetscInt}, arg7::Ptr{PetscInt}, arg8::Ptr{Ptr{PetscInt}}, arg9::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function PetscMergeIntArray(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscMergeIntArray(arg1::PetscInt, arg2::Ptr{PetscInt}, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function PetscMergeMPIIntArray(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscMergeMPIIntArray(arg1::PetscInt, arg2::Ptr{PetscMPIInt}, arg3::PetscInt, arg4::Ptr{PetscMPIInt}, arg5::Ptr{PetscInt}, arg6::Ptr{Ptr{PetscMPIInt}})::PetscErrorCode
end

function PetscParallelSortedInt(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscParallelSortedInt(arg1::MPI_Comm, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscBool})::PetscErrorCode
end

function PetscTimSort(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscTimSort(arg1::PetscInt, arg2::Ptr{Cvoid}, arg3::Csize_t, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid})::PetscErrorCode
end

function PetscIntSortSemiOrdered(arg1, arg2)
    @ccall libpetsc.PetscIntSortSemiOrdered(arg1::PetscInt, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscMPIIntSortSemiOrdered(arg1, arg2)
    @ccall libpetsc.PetscMPIIntSortSemiOrdered(arg1::PetscInt, arg2::Ptr{PetscMPIInt})::PetscErrorCode
end

function PetscRealSortSemiOrdered(arg1, arg2)
    @ccall libpetsc.PetscRealSortSemiOrdered(arg1::PetscInt, arg2::Ptr{PetscReal})::PetscErrorCode
end

function PetscTimSortWithArray(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscTimSortWithArray(arg1::PetscInt, arg2::Ptr{Cvoid}, arg3::Csize_t, arg4::Ptr{Cvoid}, arg5::Csize_t, arg6::Ptr{Cvoid}, arg7::Ptr{Cvoid})::PetscErrorCode
end

function PetscIntSortSemiOrderedWithArray(arg1, arg2, arg3)
    @ccall libpetsc.PetscIntSortSemiOrderedWithArray(arg1::PetscInt, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscMPIIntSortSemiOrderedWithArray(arg1, arg2, arg3)
    @ccall libpetsc.PetscMPIIntSortSemiOrderedWithArray(arg1::PetscInt, arg2::Ptr{PetscMPIInt}, arg3::Ptr{PetscMPIInt})::PetscErrorCode
end

function PetscRealSortSemiOrderedWithArrayInt(arg1, arg2, arg3)
    @ccall libpetsc.PetscRealSortSemiOrderedWithArrayInt(arg1::PetscInt, arg2::Ptr{PetscReal}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscSetDisplay()
    @ccall libpetsc.PetscSetDisplay()::PetscErrorCode
end

function PetscGetDisplay(arg1, arg2)
    @ccall libpetsc.PetscGetDisplay(arg1::Ptr{Cchar}, arg2::Csize_t)::PetscErrorCode
end

const PetscRandomType = Ptr{Cchar}

function PetscRandomInitializePackage()
    @ccall libpetsc.PetscRandomInitializePackage()::PetscErrorCode
end

function PetscRandomRegister(arg1, arg2)
    @ccall libpetsc.PetscRandomRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscRandomSetType(arg1, arg2)
    @ccall libpetsc.PetscRandomSetType(arg1::PetscRandom, arg2::PetscRandomType)::PetscErrorCode
end

function PetscRandomSetFromOptions(arg1)
    @ccall libpetsc.PetscRandomSetFromOptions(arg1::PetscRandom)::PetscErrorCode
end

function PetscRandomGetType(arg1, arg2)
    @ccall libpetsc.PetscRandomGetType(arg1::PetscRandom, arg2::Ptr{PetscRandomType})::PetscErrorCode
end

function PetscRandomViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.PetscRandomViewFromOptions(arg1::PetscRandom, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function PetscRandomView(arg1, arg2)
    @ccall libpetsc.PetscRandomView(arg1::PetscRandom, arg2::PetscViewer)::PetscErrorCode
end

function PetscRandomCreate(arg1, arg2)
    @ccall libpetsc.PetscRandomCreate(arg1::MPI_Comm, arg2::Ptr{PetscRandom})::PetscErrorCode
end

function PetscRandomGetValue(arg1, arg2)
    @ccall libpetsc.PetscRandomGetValue(arg1::PetscRandom, arg2::Ptr{PetscScalar})::PetscErrorCode
end

function PetscRandomGetValueReal(arg1, arg2)
    @ccall libpetsc.PetscRandomGetValueReal(arg1::PetscRandom, arg2::Ptr{PetscReal})::PetscErrorCode
end

function PetscRandomGetValues(arg1, arg2, arg3)
    @ccall libpetsc.PetscRandomGetValues(arg1::PetscRandom, arg2::PetscInt, arg3::Ptr{PetscScalar})::PetscErrorCode
end

function PetscRandomGetValuesReal(arg1, arg2, arg3)
    @ccall libpetsc.PetscRandomGetValuesReal(arg1::PetscRandom, arg2::PetscInt, arg3::Ptr{PetscReal})::PetscErrorCode
end

function PetscRandomGetInterval(arg1, arg2, arg3)
    @ccall libpetsc.PetscRandomGetInterval(arg1::PetscRandom, arg2::Ptr{PetscScalar}, arg3::Ptr{PetscScalar})::PetscErrorCode
end

function PetscRandomSetInterval(arg1, arg2, arg3)
    @ccall libpetsc.PetscRandomSetInterval(arg1::PetscRandom, arg2::PetscScalar, arg3::PetscScalar)::PetscErrorCode
end

function PetscRandomSetSeed(arg1, arg2)
    @ccall libpetsc.PetscRandomSetSeed(arg1::PetscRandom, arg2::Culong)::PetscErrorCode
end

function PetscRandomGetSeed(arg1, arg2)
    @ccall libpetsc.PetscRandomGetSeed(arg1::PetscRandom, arg2::Ptr{Culong})::PetscErrorCode
end

function PetscRandomSeed(arg1)
    @ccall libpetsc.PetscRandomSeed(arg1::PetscRandom)::PetscErrorCode
end

function PetscRandomDestroy(arg1)
    @ccall libpetsc.PetscRandomDestroy(arg1::Ptr{PetscRandom})::PetscErrorCode
end

function PetscGetFullPath(arg1, arg2, arg3)
    @ccall libpetsc.PetscGetFullPath(arg1::Ptr{Cchar}, arg2::Ptr{Cchar}, arg3::Csize_t)::PetscErrorCode
end

function PetscGetRelativePath(arg1, arg2, arg3)
    @ccall libpetsc.PetscGetRelativePath(arg1::Ptr{Cchar}, arg2::Ptr{Cchar}, arg3::Csize_t)::PetscErrorCode
end

function PetscGetWorkingDirectory(arg1, arg2)
    @ccall libpetsc.PetscGetWorkingDirectory(arg1::Ptr{Cchar}, arg2::Csize_t)::PetscErrorCode
end

function PetscGetRealPath(arg1, arg2)
    @ccall libpetsc.PetscGetRealPath(arg1::Ptr{Cchar}, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscGetHomeDirectory(arg1, arg2)
    @ccall libpetsc.PetscGetHomeDirectory(arg1::Ptr{Cchar}, arg2::Csize_t)::PetscErrorCode
end

function PetscTestFile(arg1, arg2, arg3)
    @ccall libpetsc.PetscTestFile(arg1::Ptr{Cchar}, arg2::Cchar, arg3::Ptr{PetscBool})::PetscErrorCode
end

function PetscTestDirectory(arg1, arg2, arg3)
    @ccall libpetsc.PetscTestDirectory(arg1::Ptr{Cchar}, arg2::Cchar, arg3::Ptr{PetscBool})::PetscErrorCode
end

function PetscMkdir(arg1)
    @ccall libpetsc.PetscMkdir(arg1::Ptr{Cchar})::PetscErrorCode
end

function PetscMkdtemp(arg1)
    @ccall libpetsc.PetscMkdtemp(arg1::Ptr{Cchar})::PetscErrorCode
end

function PetscRMTree(arg1)
    @ccall libpetsc.PetscRMTree(arg1::Ptr{Cchar})::PetscErrorCode
end

function PetscBinaryBigEndian()
    @ccall libpetsc.PetscBinaryBigEndian()::PetscBool
end

function PetscBinaryRead(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscBinaryRead(arg1::Cint, arg2::Ptr{Cvoid}, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::PetscDataType)::PetscErrorCode
end

function PetscBinarySynchronizedRead(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscBinarySynchronizedRead(arg1::MPI_Comm, arg2::Cint, arg3::Ptr{Cvoid}, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::PetscDataType)::PetscErrorCode
end

function PetscBinaryWrite(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscBinaryWrite(arg1::Cint, arg2::Ptr{Cvoid}, arg3::PetscInt, arg4::PetscDataType)::PetscErrorCode
end

function PetscBinarySynchronizedWrite(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscBinarySynchronizedWrite(arg1::MPI_Comm, arg2::Cint, arg3::Ptr{Cvoid}, arg4::PetscInt, arg5::PetscDataType)::PetscErrorCode
end

function PetscBinaryOpen(arg1, arg2, arg3)
    @ccall libpetsc.PetscBinaryOpen(arg1::Ptr{Cchar}, arg2::PetscFileMode, arg3::Ptr{Cint})::PetscErrorCode
end

function PetscBinaryClose(arg1)
    @ccall libpetsc.PetscBinaryClose(arg1::Cint)::PetscErrorCode
end

function PetscSharedTmp(arg1, arg2)
    @ccall libpetsc.PetscSharedTmp(arg1::MPI_Comm, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscSharedWorkingDirectory(arg1, arg2)
    @ccall libpetsc.PetscSharedWorkingDirectory(arg1::MPI_Comm, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscGetTmp(arg1, arg2, arg3)
    @ccall libpetsc.PetscGetTmp(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Csize_t)::PetscErrorCode
end

function PetscFileRetrieve(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscFileRetrieve(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Csize_t, arg5::Ptr{PetscBool})::PetscErrorCode
end

function PetscLs(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscLs(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Csize_t, arg5::Ptr{PetscBool})::PetscErrorCode
end

function PetscOpenSocket(arg1, arg2, arg3)
    @ccall libpetsc.PetscOpenSocket(arg1::Ptr{Cchar}, arg2::Cint, arg3::Ptr{Cint})::PetscErrorCode
end

function PetscBinarySeek(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscBinarySeek(arg1::Cint, arg2::off_t, arg3::PetscBinarySeekType, arg4::Ptr{off_t})::PetscErrorCode
end

function PetscBinarySynchronizedSeek(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscBinarySynchronizedSeek(arg1::MPI_Comm, arg2::Cint, arg3::off_t, arg4::PetscBinarySeekType, arg5::Ptr{off_t})::PetscErrorCode
end

function PetscByteSwap(arg1, arg2, arg3)
    @ccall libpetsc.PetscByteSwap(arg1::Ptr{Cvoid}, arg2::PetscDataType, arg3::PetscInt)::PetscErrorCode
end

function PetscSetDebugTerminal(arg1)
    @ccall libpetsc.PetscSetDebugTerminal(arg1::Ptr{Cchar})::PetscErrorCode
end

function PetscSetDebugger(arg1, arg2)
    @ccall libpetsc.PetscSetDebugger(arg1::Ptr{Cchar}, arg2::PetscBool)::PetscErrorCode
end

function PetscSetDefaultDebugger()
    @ccall libpetsc.PetscSetDefaultDebugger()::PetscErrorCode
end

function PetscSetDebuggerFromString(arg1)
    @ccall libpetsc.PetscSetDebuggerFromString(arg1::Ptr{Cchar})::PetscErrorCode
end

function PetscAttachDebugger()
    @ccall libpetsc.PetscAttachDebugger()::PetscErrorCode
end

function PetscStopForDebugger()
    @ccall libpetsc.PetscStopForDebugger()::PetscErrorCode
end

function PetscWaitOnError()
    @ccall libpetsc.PetscWaitOnError()::PetscErrorCode
end

function PetscGatherNumberOfMessages(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscGatherNumberOfMessages(arg1::MPI_Comm, arg2::Ptr{PetscMPIInt}, arg3::Ptr{PetscMPIInt}, arg4::Ptr{PetscMPIInt})::PetscErrorCode
end

function PetscGatherMessageLengths(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscGatherMessageLengths(arg1::MPI_Comm, arg2::PetscMPIInt, arg3::PetscMPIInt, arg4::Ptr{PetscMPIInt}, arg5::Ptr{Ptr{PetscMPIInt}}, arg6::Ptr{Ptr{PetscMPIInt}})::PetscErrorCode
end

function PetscGatherMessageLengths2(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscGatherMessageLengths2(arg1::MPI_Comm, arg2::PetscMPIInt, arg3::PetscMPIInt, arg4::Ptr{PetscMPIInt}, arg5::Ptr{PetscMPIInt}, arg6::Ptr{Ptr{PetscMPIInt}}, arg7::Ptr{Ptr{PetscMPIInt}}, arg8::Ptr{Ptr{PetscMPIInt}})::PetscErrorCode
end

function PetscPostIrecvInt(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscPostIrecvInt(arg1::MPI_Comm, arg2::PetscMPIInt, arg3::PetscMPIInt, arg4::Ptr{PetscMPIInt}, arg5::Ptr{PetscMPIInt}, arg6::Ptr{Ptr{Ptr{PetscInt}}}, arg7::Ptr{Ptr{MPI_Request}})::PetscErrorCode
end

function PetscPostIrecvScalar(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscPostIrecvScalar(arg1::MPI_Comm, arg2::PetscMPIInt, arg3::PetscMPIInt, arg4::Ptr{PetscMPIInt}, arg5::Ptr{PetscMPIInt}, arg6::Ptr{Ptr{Ptr{PetscScalar}}}, arg7::Ptr{Ptr{MPI_Request}})::PetscErrorCode
end

function PetscCommBuildTwoSided(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.PetscCommBuildTwoSided(arg1::MPI_Comm, arg2::PetscMPIInt, arg3::MPI_Datatype, arg4::PetscMPIInt, arg5::Ptr{PetscMPIInt}, arg6::Ptr{Cvoid}, arg7::Ptr{PetscMPIInt}, arg8::Ptr{Ptr{PetscMPIInt}}, arg9::Ptr{Cvoid})::PetscErrorCode
end

function PetscCommBuildTwoSidedF(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, send, recv, ctx)
    @ccall libpetsc.PetscCommBuildTwoSidedF(arg1::MPI_Comm, arg2::PetscMPIInt, arg3::MPI_Datatype, arg4::PetscMPIInt, arg5::Ptr{PetscMPIInt}, arg6::Ptr{Cvoid}, arg7::Ptr{PetscMPIInt}, arg8::Ptr{Ptr{PetscMPIInt}}, arg9::Ptr{Cvoid}, arg10::PetscMPIInt, send::Ptr{Cvoid}, recv::Ptr{Cvoid}, ctx::Ptr{Cvoid})::PetscErrorCode
end

function PetscCommBuildTwoSidedFReq(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, send, recv, ctx)
    @ccall libpetsc.PetscCommBuildTwoSidedFReq(arg1::MPI_Comm, arg2::PetscMPIInt, arg3::MPI_Datatype, arg4::PetscMPIInt, arg5::Ptr{PetscMPIInt}, arg6::Ptr{Cvoid}, arg7::Ptr{PetscMPIInt}, arg8::Ptr{Ptr{PetscMPIInt}}, arg9::Ptr{Cvoid}, arg10::PetscMPIInt, arg11::Ptr{Ptr{MPI_Request}}, arg12::Ptr{Ptr{MPI_Request}}, send::Ptr{Cvoid}, recv::Ptr{Cvoid}, ctx::Ptr{Cvoid})::PetscErrorCode
end

function PetscCommBuildTwoSidedSetType(arg1, arg2)
    @ccall libpetsc.PetscCommBuildTwoSidedSetType(arg1::MPI_Comm, arg2::PetscBuildTwoSidedType)::PetscErrorCode
end

function PetscCommBuildTwoSidedGetType(arg1, arg2)
    @ccall libpetsc.PetscCommBuildTwoSidedGetType(arg1::MPI_Comm, arg2::Ptr{PetscBuildTwoSidedType})::PetscErrorCode
end

function PetscSSEIsEnabled(arg1, arg2, arg3)
    @ccall libpetsc.PetscSSEIsEnabled(arg1::MPI_Comm, arg2::Ptr{PetscBool}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function PetscObjectComm(arg1)
    @ccall libpetsc.PetscObjectComm(arg1::PetscObject)::MPI_Comm
end

function PetscSubcommParent(scomm)
    @ccall libpetsc.PetscSubcommParent(scomm::PetscSubcomm)::MPI_Comm
end

function PetscSubcommChild(scomm)
    @ccall libpetsc.PetscSubcommChild(scomm::PetscSubcomm)::MPI_Comm
end

function PetscSubcommContiguousParent(scomm)
    @ccall libpetsc.PetscSubcommContiguousParent(scomm::PetscSubcomm)::MPI_Comm
end

function PetscSubcommCreate(arg1, arg2)
    @ccall libpetsc.PetscSubcommCreate(arg1::MPI_Comm, arg2::Ptr{PetscSubcomm})::PetscErrorCode
end

function PetscSubcommDestroy(arg1)
    @ccall libpetsc.PetscSubcommDestroy(arg1::Ptr{PetscSubcomm})::PetscErrorCode
end

function PetscSubcommSetNumber(arg1, arg2)
    @ccall libpetsc.PetscSubcommSetNumber(arg1::PetscSubcomm, arg2::PetscInt)::PetscErrorCode
end

function PetscSubcommSetType(arg1, arg2)
    @ccall libpetsc.PetscSubcommSetType(arg1::PetscSubcomm, arg2::PetscSubcommType)::PetscErrorCode
end

function PetscSubcommSetTypeGeneral(arg1, arg2, arg3)
    @ccall libpetsc.PetscSubcommSetTypeGeneral(arg1::PetscSubcomm, arg2::PetscMPIInt, arg3::PetscMPIInt)::PetscErrorCode
end

function PetscSubcommView(arg1, arg2)
    @ccall libpetsc.PetscSubcommView(arg1::PetscSubcomm, arg2::PetscViewer)::PetscErrorCode
end

function PetscSubcommSetFromOptions(arg1)
    @ccall libpetsc.PetscSubcommSetFromOptions(arg1::PetscSubcomm)::PetscErrorCode
end

function PetscSubcommSetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.PetscSubcommSetOptionsPrefix(arg1::PetscSubcomm, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscSubcommGetParent(arg1, arg2)
    @ccall libpetsc.PetscSubcommGetParent(arg1::PetscSubcomm, arg2::Ptr{MPI_Comm})::PetscErrorCode
end

function PetscSubcommGetContiguousParent(arg1, arg2)
    @ccall libpetsc.PetscSubcommGetContiguousParent(arg1::PetscSubcomm, arg2::Ptr{MPI_Comm})::PetscErrorCode
end

function PetscSubcommGetChild(arg1, arg2)
    @ccall libpetsc.PetscSubcommGetChild(arg1::PetscSubcomm, arg2::Ptr{MPI_Comm})::PetscErrorCode
end

function PetscHeapCreate(arg1, arg2)
    @ccall libpetsc.PetscHeapCreate(arg1::PetscInt, arg2::Ptr{PetscHeap})::PetscErrorCode
end

function PetscHeapAdd(arg1, arg2, arg3)
    @ccall libpetsc.PetscHeapAdd(arg1::PetscHeap, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function PetscHeapPop(arg1, arg2, arg3)
    @ccall libpetsc.PetscHeapPop(arg1::PetscHeap, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscHeapPeek(arg1, arg2, arg3)
    @ccall libpetsc.PetscHeapPeek(arg1::PetscHeap, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscHeapStash(arg1, arg2, arg3)
    @ccall libpetsc.PetscHeapStash(arg1::PetscHeap, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function PetscHeapUnstash(arg1)
    @ccall libpetsc.PetscHeapUnstash(arg1::PetscHeap)::PetscErrorCode
end

function PetscHeapDestroy(arg1)
    @ccall libpetsc.PetscHeapDestroy(arg1::Ptr{PetscHeap})::PetscErrorCode
end

function PetscHeapView(arg1, arg2)
    @ccall libpetsc.PetscHeapView(arg1::PetscHeap, arg2::PetscViewer)::PetscErrorCode
end

function PetscProcessPlacementView(arg1)
    @ccall libpetsc.PetscProcessPlacementView(arg1::PetscViewer)::PetscErrorCode
end

function PetscShmCommGet(arg1, arg2)
    @ccall libpetsc.PetscShmCommGet(arg1::MPI_Comm, arg2::Ptr{PetscShmComm})::PetscErrorCode
end

function PetscShmCommGlobalToLocal(arg1, arg2, arg3)
    @ccall libpetsc.PetscShmCommGlobalToLocal(arg1::PetscShmComm, arg2::PetscMPIInt, arg3::Ptr{PetscMPIInt})::PetscErrorCode
end

function PetscShmCommLocalToGlobal(arg1, arg2, arg3)
    @ccall libpetsc.PetscShmCommLocalToGlobal(arg1::PetscShmComm, arg2::PetscMPIInt, arg3::Ptr{PetscMPIInt})::PetscErrorCode
end

function PetscShmCommGetMpiShmComm(arg1, arg2)
    @ccall libpetsc.PetscShmCommGetMpiShmComm(arg1::PetscShmComm, arg2::Ptr{MPI_Comm})::PetscErrorCode
end

function PetscOmpCtrlCreate(arg1, arg2, arg3)
    @ccall libpetsc.PetscOmpCtrlCreate(arg1::MPI_Comm, arg2::PetscInt, arg3::Ptr{PetscOmpCtrl})::PetscErrorCode
end

function PetscOmpCtrlGetOmpComms(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscOmpCtrlGetOmpComms(arg1::PetscOmpCtrl, arg2::Ptr{MPI_Comm}, arg3::Ptr{MPI_Comm}, arg4::Ptr{PetscBool})::PetscErrorCode
end

function PetscOmpCtrlDestroy(arg1)
    @ccall libpetsc.PetscOmpCtrlDestroy(arg1::Ptr{PetscOmpCtrl})::PetscErrorCode
end

function PetscOmpCtrlBarrier(arg1)
    @ccall libpetsc.PetscOmpCtrlBarrier(arg1::PetscOmpCtrl)::PetscErrorCode
end

function PetscOmpCtrlOmpRegionOnMasterBegin(arg1)
    @ccall libpetsc.PetscOmpCtrlOmpRegionOnMasterBegin(arg1::PetscOmpCtrl)::PetscErrorCode
end

function PetscOmpCtrlOmpRegionOnMasterEnd(arg1)
    @ccall libpetsc.PetscOmpCtrlOmpRegionOnMasterEnd(arg1::PetscOmpCtrl)::PetscErrorCode
end

function PetscSegBufferCreate(arg1, arg2, arg3)
    @ccall libpetsc.PetscSegBufferCreate(arg1::Csize_t, arg2::Csize_t, arg3::Ptr{PetscSegBuffer})::PetscErrorCode
end

function PetscSegBufferDestroy(arg1)
    @ccall libpetsc.PetscSegBufferDestroy(arg1::Ptr{PetscSegBuffer})::PetscErrorCode
end

function PetscSegBufferGet(arg1, arg2, arg3)
    @ccall libpetsc.PetscSegBufferGet(arg1::PetscSegBuffer, arg2::Csize_t, arg3::Ptr{Cvoid})::PetscErrorCode
end

function PetscSegBufferExtractAlloc(arg1, arg2)
    @ccall libpetsc.PetscSegBufferExtractAlloc(arg1::PetscSegBuffer, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscSegBufferExtractTo(arg1, arg2)
    @ccall libpetsc.PetscSegBufferExtractTo(arg1::PetscSegBuffer, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscSegBufferExtractInPlace(arg1, arg2)
    @ccall libpetsc.PetscSegBufferExtractInPlace(arg1::PetscSegBuffer, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscSegBufferGetSize(arg1, arg2)
    @ccall libpetsc.PetscSegBufferGetSize(arg1::PetscSegBuffer, arg2::Ptr{Csize_t})::PetscErrorCode
end

function PetscSegBufferUnuse(arg1, arg2)
    @ccall libpetsc.PetscSegBufferUnuse(arg1::PetscSegBuffer, arg2::Csize_t)::PetscErrorCode
end

function PetscSegBufferGetInts(seg, count, slot)
    @ccall libpetsc.PetscSegBufferGetInts(seg::PetscSegBuffer, count::Csize_t, slot::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function PetscOptionsHelpPrintedDestroy(arg1)
    @ccall libpetsc.PetscOptionsHelpPrintedDestroy(arg1::Ptr{PetscOptionsHelpPrinted})::PetscErrorCode
end

function PetscOptionsHelpPrintedCreate(arg1)
    @ccall libpetsc.PetscOptionsHelpPrintedCreate(arg1::Ptr{PetscOptionsHelpPrinted})::PetscErrorCode
end

function PetscOptionsHelpPrintedCheck(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscOptionsHelpPrintedCheck(arg1::PetscOptionsHelpPrinted, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{PetscBool})::PetscErrorCode
end

function PetscCitationsRegister(cit, set)
    @ccall libpetsc.PetscCitationsRegister(cit::Ptr{Cchar}, set::Ptr{PetscBool})::PetscErrorCode
end

function PetscURLShorten(arg1, arg2, arg3)
    @ccall libpetsc.PetscURLShorten(arg1::Ptr{Cchar}, arg2::Ptr{Cchar}, arg3::Csize_t)::PetscErrorCode
end

function PetscGoogleDriveAuthorize(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscGoogleDriveAuthorize(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Csize_t)::PetscErrorCode
end

function PetscGoogleDriveRefresh(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscGoogleDriveRefresh(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Csize_t)::PetscErrorCode
end

function PetscGoogleDriveUpload(arg1, arg2, arg3)
    @ccall libpetsc.PetscGoogleDriveUpload(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar})::PetscErrorCode
end

function PetscBoxAuthorize(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscBoxAuthorize(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Csize_t)::PetscErrorCode
end

function PetscBoxRefresh(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscBoxRefresh(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Csize_t)::PetscErrorCode
end

function PetscGlobusGetTransfers(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscGlobusGetTransfers(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Csize_t)::PetscErrorCode
end

function PetscTextBelt(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscTextBelt(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{PetscBool})::PetscErrorCode
end

function PetscTellMyCell(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscTellMyCell(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{PetscBool})::PetscErrorCode
end

function PetscPullJSONValue(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscPullJSONValue(arg1::Ptr{Cchar}, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Csize_t, arg5::Ptr{PetscBool})::PetscErrorCode
end

function PetscPushJSONValue(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscPushJSONValue(arg1::Ptr{Cchar}, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Csize_t)::PetscErrorCode
end

function MPIU_Win_allocate_shared(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MPIU_Win_allocate_shared(arg1::MPI_Aint, arg2::PetscMPIInt, arg3::MPI_Info, arg4::MPI_Comm, arg5::Ptr{Cvoid}, arg6::Ptr{MPI_Win})::PetscErrorCode
end

function MPIU_Win_shared_query(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MPIU_Win_shared_query(arg1::MPI_Win, arg2::PetscMPIInt, arg3::Ptr{MPI_Aint}, arg4::Ptr{PetscMPIInt}, arg5::Ptr{Cvoid})::PetscErrorCode
end

function PetscHasExternalPackage(arg1, arg2)
    @ccall libpetsc.PetscHasExternalPackage(arg1::Ptr{Cchar}, arg2::Ptr{PetscBool})::PetscErrorCode
end

mutable struct _n_PetscBag end

const PetscBag = Ptr{_n_PetscBag}

mutable struct _n_PetscBagItem end

const PetscBagItem = Ptr{_n_PetscBagItem}

function PetscBagCreate(arg1, arg2, arg3)
    @ccall libpetsc.PetscBagCreate(arg1::MPI_Comm, arg2::Csize_t, arg3::Ptr{PetscBag})::PetscErrorCode
end

function PetscBagDestroy(arg1)
    @ccall libpetsc.PetscBagDestroy(arg1::Ptr{PetscBag})::PetscErrorCode
end

function PetscBagGetData(arg1, arg2)
    @ccall libpetsc.PetscBagGetData(arg1::PetscBag, arg2::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscBagRegisterReal(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscBagRegisterReal(arg1::PetscBag, arg2::Ptr{Cvoid}, arg3::PetscReal, arg4::Ptr{Cchar}, arg5::Ptr{Cchar})::PetscErrorCode
end

function PetscBagRegisterRealArray(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscBagRegisterRealArray(arg1::PetscBag, arg2::Ptr{Cvoid}, arg3::PetscInt, arg4::Ptr{Cchar}, arg5::Ptr{Cchar})::PetscErrorCode
end

function PetscBagRegisterString(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscBagRegisterString(arg1::PetscBag, arg2::Ptr{Cvoid}, arg3::PetscInt, arg4::Ptr{Cchar}, arg5::Ptr{Cchar}, arg6::Ptr{Cchar})::PetscErrorCode
end

function PetscBagRegisterScalar(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscBagRegisterScalar(arg1::PetscBag, arg2::Ptr{Cvoid}, arg3::PetscScalar, arg4::Ptr{Cchar}, arg5::Ptr{Cchar})::PetscErrorCode
end

function PetscBagRegisterInt(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscBagRegisterInt(arg1::PetscBag, arg2::Ptr{Cvoid}, arg3::PetscInt, arg4::Ptr{Cchar}, arg5::Ptr{Cchar})::PetscErrorCode
end

function PetscBagRegisterInt64(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscBagRegisterInt64(arg1::PetscBag, arg2::Ptr{Cvoid}, arg3::PetscInt64, arg4::Ptr{Cchar}, arg5::Ptr{Cchar})::PetscErrorCode
end

function PetscBagRegisterIntArray(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscBagRegisterIntArray(arg1::PetscBag, arg2::Ptr{Cvoid}, arg3::PetscInt, arg4::Ptr{Cchar}, arg5::Ptr{Cchar})::PetscErrorCode
end

function PetscBagRegisterEnum(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscBagRegisterEnum(arg1::PetscBag, arg2::Ptr{Cvoid}, arg3::Ptr{Ptr{Cchar}}, arg4::PetscEnum, arg5::Ptr{Cchar}, arg6::Ptr{Cchar})::PetscErrorCode
end

function PetscBagRegisterBool(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscBagRegisterBool(arg1::PetscBag, arg2::Ptr{Cvoid}, arg3::PetscBool, arg4::Ptr{Cchar}, arg5::Ptr{Cchar})::PetscErrorCode
end

function PetscBagRegisterBoolArray(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscBagRegisterBoolArray(arg1::PetscBag, arg2::Ptr{Cvoid}, arg3::PetscInt, arg4::Ptr{Cchar}, arg5::Ptr{Cchar})::PetscErrorCode
end

function PetscBagGetNames(arg1, arg2)
    @ccall libpetsc.PetscBagGetNames(arg1::PetscBag, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscBagSetFromOptions(arg1)
    @ccall libpetsc.PetscBagSetFromOptions(arg1::PetscBag)::PetscErrorCode
end

function PetscBagGetName(arg1, arg2)
    @ccall libpetsc.PetscBagGetName(arg1::PetscBag, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscBagSetName(arg1, arg2, arg3)
    @ccall libpetsc.PetscBagSetName(arg1::PetscBag, arg2::Ptr{Cchar}, arg3::Ptr{Cchar})::PetscErrorCode
end

function PetscBagSetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.PetscBagSetOptionsPrefix(arg1::PetscBag, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscBagView(arg1, arg2)
    @ccall libpetsc.PetscBagView(arg1::PetscBag, arg2::PetscViewer)::PetscErrorCode
end

function PetscBagLoad(arg1, arg2)
    @ccall libpetsc.PetscBagLoad(arg1::PetscViewer, arg2::PetscBag)::PetscErrorCode
end

function PetscBagViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.PetscBagViewFromOptions(arg1::PetscBag, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function PetscBagSetViewer(arg1, arg2)
    @ccall libpetsc.PetscBagSetViewer(arg1::PetscBag, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscBagSetLoader(arg1, arg2)
    @ccall libpetsc.PetscBagSetLoader(arg1::PetscBag, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscBagSetDestroy(arg1, arg2)
    @ccall libpetsc.PetscBagSetDestroy(arg1::PetscBag, arg2::Ptr{Cvoid})::PetscErrorCode
end

const PetscDrawType = Ptr{Cchar}

mutable struct _p_PetscDraw end

const PetscDraw = Ptr{_p_PetscDraw}

mutable struct _p_PetscDrawAxis end

const PetscDrawAxis = Ptr{_p_PetscDrawAxis}

mutable struct _p_PetscDrawLG end

const PetscDrawLG = Ptr{_p_PetscDrawLG}

mutable struct _p_PetscDrawSP end

const PetscDrawSP = Ptr{_p_PetscDrawSP}

mutable struct _p_PetscDrawHG end

const PetscDrawHG = Ptr{_p_PetscDrawHG}

mutable struct _p_PetscDrawBar end

const PetscDrawBar = Ptr{_p_PetscDrawBar}

const PetscViewerType = Ptr{Cchar}

function PetscViewerInitializePackage()
    @ccall libpetsc.PetscViewerInitializePackage()::PetscErrorCode
end

function PetscViewerRegister(arg1, arg2)
    @ccall libpetsc.PetscViewerRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscViewerCreate(arg1, arg2)
    @ccall libpetsc.PetscViewerCreate(arg1::MPI_Comm, arg2::Ptr{PetscViewer})::PetscErrorCode
end

function PetscViewerSetFromOptions(arg1)
    @ccall libpetsc.PetscViewerSetFromOptions(arg1::PetscViewer)::PetscErrorCode
end

function PetscViewerASCIIOpenWithFILE(arg1, arg2, arg3)
    @ccall libpetsc.PetscViewerASCIIOpenWithFILE(arg1::MPI_Comm, arg2::Ptr{Libc.FILE}, arg3::Ptr{PetscViewer})::PetscErrorCode
end

function PetscViewerASCIIOpen(arg1, arg2, arg3)
    @ccall libpetsc.PetscViewerASCIIOpen(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{PetscViewer})::PetscErrorCode
end

function PetscViewerASCIISetFILE(arg1, arg2)
    @ccall libpetsc.PetscViewerASCIISetFILE(arg1::PetscViewer, arg2::Ptr{Libc.FILE})::PetscErrorCode
end

function PetscViewerBinaryOpen(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscViewerBinaryOpen(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::PetscFileMode, arg4::Ptr{PetscViewer})::PetscErrorCode
end

function PetscViewerADIOSOpen(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscViewerADIOSOpen(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::PetscFileMode, arg4::Ptr{PetscViewer})::PetscErrorCode
end

function PetscViewerADIOS2Open(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscViewerADIOS2Open(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::PetscFileMode, arg4::Ptr{PetscViewer})::PetscErrorCode
end

function PetscViewerBinaryGetFlowControl(arg1, arg2)
    @ccall libpetsc.PetscViewerBinaryGetFlowControl(arg1::PetscViewer, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscViewerBinarySetFlowControl(arg1, arg2)
    @ccall libpetsc.PetscViewerBinarySetFlowControl(arg1::PetscViewer, arg2::PetscInt)::PetscErrorCode
end

function PetscViewerBinarySetUseMPIIO(arg1, arg2)
    @ccall libpetsc.PetscViewerBinarySetUseMPIIO(arg1::PetscViewer, arg2::PetscBool)::PetscErrorCode
end

function PetscViewerBinaryGetUseMPIIO(arg1, arg2)
    @ccall libpetsc.PetscViewerBinaryGetUseMPIIO(arg1::PetscViewer, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscViewerBinaryGetMPIIODescriptor(arg1, arg2)
    @ccall libpetsc.PetscViewerBinaryGetMPIIODescriptor(arg1::PetscViewer, arg2::Ptr{MPI_File})::PetscErrorCode
end

function PetscViewerBinaryGetMPIIOOffset(arg1, arg2)
    @ccall libpetsc.PetscViewerBinaryGetMPIIOOffset(arg1::PetscViewer, arg2::Ptr{MPI_Offset})::PetscErrorCode
end

function PetscViewerBinaryAddMPIIOOffset(arg1, arg2)
    @ccall libpetsc.PetscViewerBinaryAddMPIIOOffset(arg1::PetscViewer, arg2::MPI_Offset)::PetscErrorCode
end

function PetscViewerSocketOpen(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscViewerSocketOpen(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Cint, arg4::Ptr{PetscViewer})::PetscErrorCode
end

function PetscViewerStringOpen(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscViewerStringOpen(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Csize_t, arg4::Ptr{PetscViewer})::PetscErrorCode
end

function PetscViewerDrawOpen(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscViewerDrawOpen(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Cint, arg5::Cint, arg6::Cint, arg7::Cint, arg8::Ptr{PetscViewer})::PetscErrorCode
end

function PetscViewerDrawSetDrawType(arg1, arg2)
    @ccall libpetsc.PetscViewerDrawSetDrawType(arg1::PetscViewer, arg2::PetscDrawType)::PetscErrorCode
end

function PetscViewerDrawGetDrawType(arg1, arg2)
    @ccall libpetsc.PetscViewerDrawGetDrawType(arg1::PetscViewer, arg2::Ptr{PetscDrawType})::PetscErrorCode
end

function PetscViewerDrawSetTitle(arg1, arg2)
    @ccall libpetsc.PetscViewerDrawSetTitle(arg1::PetscViewer, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscViewerDrawGetTitle(arg1, arg2)
    @ccall libpetsc.PetscViewerDrawGetTitle(arg1::PetscViewer, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscViewerDrawGetDraw(arg1, arg2, arg3)
    @ccall libpetsc.PetscViewerDrawGetDraw(arg1::PetscViewer, arg2::PetscInt, arg3::Ptr{PetscDraw})::PetscErrorCode
end

function PetscViewerDrawBaseAdd(arg1, arg2)
    @ccall libpetsc.PetscViewerDrawBaseAdd(arg1::PetscViewer, arg2::PetscInt)::PetscErrorCode
end

function PetscViewerDrawBaseSet(arg1, arg2)
    @ccall libpetsc.PetscViewerDrawBaseSet(arg1::PetscViewer, arg2::PetscInt)::PetscErrorCode
end

function PetscViewerDrawGetDrawLG(arg1, arg2, arg3)
    @ccall libpetsc.PetscViewerDrawGetDrawLG(arg1::PetscViewer, arg2::PetscInt, arg3::Ptr{PetscDrawLG})::PetscErrorCode
end

function PetscViewerDrawGetDrawAxis(arg1, arg2, arg3)
    @ccall libpetsc.PetscViewerDrawGetDrawAxis(arg1::PetscViewer, arg2::PetscInt, arg3::Ptr{PetscDrawAxis})::PetscErrorCode
end

function PetscViewerMathematicaOpen(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscViewerMathematicaOpen(arg1::MPI_Comm, arg2::Cint, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{PetscViewer})::PetscErrorCode
end

function PetscViewerSiloOpen(arg1, arg2, arg3)
    @ccall libpetsc.PetscViewerSiloOpen(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{PetscViewer})::PetscErrorCode
end

function PetscViewerMatlabOpen(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscViewerMatlabOpen(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::PetscFileMode, arg4::Ptr{PetscViewer})::PetscErrorCode
end

@enum PetscViewerGLVisType::UInt32 begin
    PETSC_VIEWER_GLVIS_DUMP = 0
    PETSC_VIEWER_GLVIS_SOCKET = 1
end

function PetscViewerGLVisOpen(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscViewerGLVisOpen(arg1::MPI_Comm, arg2::PetscViewerGLVisType, arg3::Ptr{Cchar}, arg4::PetscInt, arg5::Ptr{PetscViewer})::PetscErrorCode
end

function PetscViewerGLVisSetPrecision(arg1, arg2)
    @ccall libpetsc.PetscViewerGLVisSetPrecision(arg1::PetscViewer, arg2::PetscInt)::PetscErrorCode
end

function PetscViewerGLVisSetSnapId(arg1, arg2)
    @ccall libpetsc.PetscViewerGLVisSetSnapId(arg1::PetscViewer, arg2::PetscInt)::PetscErrorCode
end

function PetscViewerGLVisSetFields(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscViewerGLVisSetFields(arg1::PetscViewer, arg2::PetscInt, arg3::Ptr{Ptr{Cchar}}, arg4::Ptr{PetscInt}, arg5::Ptr{Cvoid}, arg6::Ptr{PetscObject}, arg7::Ptr{Cvoid}, arg8::Ptr{Cvoid})::PetscErrorCode
end

function PetscViewerGetType(arg1, arg2)
    @ccall libpetsc.PetscViewerGetType(arg1::PetscViewer, arg2::Ptr{PetscViewerType})::PetscErrorCode
end

function PetscViewerSetType(arg1, arg2)
    @ccall libpetsc.PetscViewerSetType(arg1::PetscViewer, arg2::PetscViewerType)::PetscErrorCode
end

function PetscViewerDestroy(arg1)
    @ccall libpetsc.PetscViewerDestroy(arg1::Ptr{PetscViewer})::PetscErrorCode
end

function PetscViewerGetSubViewer(arg1, arg2, arg3)
    @ccall libpetsc.PetscViewerGetSubViewer(arg1::PetscViewer, arg2::MPI_Comm, arg3::Ptr{PetscViewer})::PetscErrorCode
end

function PetscViewerRestoreSubViewer(arg1, arg2, arg3)
    @ccall libpetsc.PetscViewerRestoreSubViewer(arg1::PetscViewer, arg2::MPI_Comm, arg3::Ptr{PetscViewer})::PetscErrorCode
end

function PetscViewerSetUp(arg1)
    @ccall libpetsc.PetscViewerSetUp(arg1::PetscViewer)::PetscErrorCode
end

function PetscViewerView(arg1, arg2)
    @ccall libpetsc.PetscViewerView(arg1::PetscViewer, arg2::PetscViewer)::PetscErrorCode
end

function PetscViewerViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.PetscViewerViewFromOptions(arg1::PetscViewer, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function PetscViewerSetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.PetscViewerSetOptionsPrefix(arg1::PetscViewer, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscViewerAppendOptionsPrefix(arg1, arg2)
    @ccall libpetsc.PetscViewerAppendOptionsPrefix(arg1::PetscViewer, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscViewerGetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.PetscViewerGetOptionsPrefix(arg1::PetscViewer, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscViewerReadable(arg1, arg2)
    @ccall libpetsc.PetscViewerReadable(arg1::PetscViewer, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscViewerWritable(arg1, arg2)
    @ccall libpetsc.PetscViewerWritable(arg1::PetscViewer, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscViewerCheckReadable(arg1)
    @ccall libpetsc.PetscViewerCheckReadable(arg1::PetscViewer)::PetscErrorCode
end

function PetscViewerCheckWritable(arg1)
    @ccall libpetsc.PetscViewerCheckWritable(arg1::PetscViewer)::PetscErrorCode
end

@enum PetscViewerFormat::UInt32 begin
    PETSC_VIEWER_DEFAULT = 0
    PETSC_VIEWER_ASCII_MATLAB = 1
    PETSC_VIEWER_ASCII_MATHEMATICA = 2
    PETSC_VIEWER_ASCII_IMPL = 3
    PETSC_VIEWER_ASCII_INFO = 4
    PETSC_VIEWER_ASCII_INFO_DETAIL = 5
    PETSC_VIEWER_ASCII_COMMON = 6
    PETSC_VIEWER_ASCII_SYMMODU = 7
    PETSC_VIEWER_ASCII_INDEX = 8
    PETSC_VIEWER_ASCII_DENSE = 9
    PETSC_VIEWER_ASCII_MATRIXMARKET = 10
    PETSC_VIEWER_ASCII_VTK_DEPRECATED = 11
    # PETSC_VIEWER_ASCII_VTK = 11
    PETSC_VIEWER_ASCII_VTK_CELL_DEPRECATED = 12
    # PETSC_VIEWER_ASCII_VTK_CELL = 12
    PETSC_VIEWER_ASCII_VTK_COORDS_DEPRECATED = 13
    # PETSC_VIEWER_ASCII_VTK_COORDS = 13
    PETSC_VIEWER_ASCII_PCICE = 14
    PETSC_VIEWER_ASCII_PYTHON = 15
    PETSC_VIEWER_ASCII_FACTOR_INFO = 16
    PETSC_VIEWER_ASCII_LATEX = 17
    PETSC_VIEWER_ASCII_XML = 18
    PETSC_VIEWER_ASCII_FLAMEGRAPH = 19
    PETSC_VIEWER_ASCII_GLVIS = 20
    PETSC_VIEWER_ASCII_CSV = 21
    PETSC_VIEWER_DRAW_BASIC = 22
    PETSC_VIEWER_DRAW_LG = 23
    PETSC_VIEWER_DRAW_LG_XRANGE = 24
    PETSC_VIEWER_DRAW_CONTOUR = 25
    PETSC_VIEWER_DRAW_PORTS = 26
    PETSC_VIEWER_VTK_VTS = 27
    PETSC_VIEWER_VTK_VTR = 28
    PETSC_VIEWER_VTK_VTU = 29
    PETSC_VIEWER_BINARY_MATLAB = 30
    PETSC_VIEWER_NATIVE = 31
    PETSC_VIEWER_HDF5_PETSC = 32
    PETSC_VIEWER_HDF5_VIZ = 33
    PETSC_VIEWER_HDF5_XDMF = 34
    PETSC_VIEWER_HDF5_MAT = 35
    PETSC_VIEWER_NOFORMAT = 36
    PETSC_VIEWER_LOAD_BALANCE = 37
    PETSC_VIEWER_FAILED = 38
end

function PetscViewerSetFormat(arg1, arg2)
    @ccall libpetsc.PetscViewerSetFormat(arg1::PetscViewer, arg2::PetscViewerFormat)::PetscErrorCode
end

function PetscViewerPushFormat(arg1, arg2)
    @ccall libpetsc.PetscViewerPushFormat(arg1::PetscViewer, arg2::PetscViewerFormat)::PetscErrorCode
end

function PetscViewerPopFormat(arg1)
    @ccall libpetsc.PetscViewerPopFormat(arg1::PetscViewer)::PetscErrorCode
end

function PetscViewerGetFormat(arg1, arg2)
    @ccall libpetsc.PetscViewerGetFormat(arg1::PetscViewer, arg2::Ptr{PetscViewerFormat})::PetscErrorCode
end

function PetscViewerFlush(arg1)
    @ccall libpetsc.PetscViewerFlush(arg1::PetscViewer)::PetscErrorCode
end

function PetscOptionsPushGetViewerOff(arg1)
    @ccall libpetsc.PetscOptionsPushGetViewerOff(arg1::PetscBool)::PetscErrorCode
end

function PetscOptionsPopGetViewerOff()
    @ccall libpetsc.PetscOptionsPopGetViewerOff()::PetscErrorCode
end

function PetscOptionsGetViewerOff(arg1)
    @ccall libpetsc.PetscOptionsGetViewerOff(arg1::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsGetViewer(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscOptionsGetViewer(arg1::MPI_Comm, arg2::PetscOptions, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{PetscViewer}, arg6::Ptr{PetscViewerFormat}, arg7::Ptr{PetscBool})::PetscErrorCode
end

function PetscOptionsViewer_Private(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscOptionsViewer_Private(arg1::Ptr{PetscOptionItems}, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{PetscViewer}, arg6::Ptr{PetscViewerFormat}, arg7::Ptr{PetscBool})::PetscErrorCode
end

mutable struct PetscViewerAndFormat
    viewer::PetscViewer
    format::PetscViewerFormat
    lg::PetscDrawLG
    data::Ptr{Cvoid}
    PetscViewerAndFormat() = new()
end

function PetscViewerAndFormatCreate(arg1, arg2, arg3)
    @ccall libpetsc.PetscViewerAndFormatCreate(arg1::PetscViewer, arg2::PetscViewerFormat, arg3::Ptr{Ptr{PetscViewerAndFormat}})::PetscErrorCode
end

function PetscViewerAndFormatDestroy(arg1)
    @ccall libpetsc.PetscViewerAndFormatDestroy(arg1::Ptr{Ptr{PetscViewerAndFormat}})::PetscErrorCode
end

function PetscViewerASCIIGetPointer(arg1, arg2)
    @ccall libpetsc.PetscViewerASCIIGetPointer(arg1::PetscViewer, arg2::Ptr{Ptr{Libc.FILE}})::PetscErrorCode
end

function PetscViewerFileGetMode(arg1, arg2)
    @ccall libpetsc.PetscViewerFileGetMode(arg1::PetscViewer, arg2::Ptr{PetscFileMode})::PetscErrorCode
end

function PetscViewerFileSetMode(arg1, arg2)
    @ccall libpetsc.PetscViewerFileSetMode(arg1::PetscViewer, arg2::PetscFileMode)::PetscErrorCode
end

function PetscViewerRead(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscViewerRead(arg1::PetscViewer, arg2::Ptr{Cvoid}, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::PetscDataType)::PetscErrorCode
end

function PetscViewerASCIIPushSynchronized(arg1)
    @ccall libpetsc.PetscViewerASCIIPushSynchronized(arg1::PetscViewer)::PetscErrorCode
end

function PetscViewerASCIIPopSynchronized(arg1)
    @ccall libpetsc.PetscViewerASCIIPopSynchronized(arg1::PetscViewer)::PetscErrorCode
end

function PetscViewerASCIIPushTab(arg1)
    @ccall libpetsc.PetscViewerASCIIPushTab(arg1::PetscViewer)::PetscErrorCode
end

function PetscViewerASCIIPopTab(arg1)
    @ccall libpetsc.PetscViewerASCIIPopTab(arg1::PetscViewer)::PetscErrorCode
end

function PetscViewerASCIIUseTabs(arg1, arg2)
    @ccall libpetsc.PetscViewerASCIIUseTabs(arg1::PetscViewer, arg2::PetscBool)::PetscErrorCode
end

function PetscViewerASCIISetTab(arg1, arg2)
    @ccall libpetsc.PetscViewerASCIISetTab(arg1::PetscViewer, arg2::PetscInt)::PetscErrorCode
end

function PetscViewerASCIIGetTab(arg1, arg2)
    @ccall libpetsc.PetscViewerASCIIGetTab(arg1::PetscViewer, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscViewerASCIIAddTab(arg1, arg2)
    @ccall libpetsc.PetscViewerASCIIAddTab(arg1::PetscViewer, arg2::PetscInt)::PetscErrorCode
end

function PetscViewerASCIISubtractTab(arg1, arg2)
    @ccall libpetsc.PetscViewerASCIISubtractTab(arg1::PetscViewer, arg2::PetscInt)::PetscErrorCode
end

function PetscViewerASCIIRead(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscViewerASCIIRead(arg1::PetscViewer, arg2::Ptr{Cvoid}, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::PetscDataType)::PetscErrorCode
end

function PetscViewerBinaryGetDescriptor(arg1, arg2)
    @ccall libpetsc.PetscViewerBinaryGetDescriptor(arg1::PetscViewer, arg2::Ptr{Cint})::PetscErrorCode
end

function PetscViewerBinaryGetInfoPointer(arg1, arg2)
    @ccall libpetsc.PetscViewerBinaryGetInfoPointer(arg1::PetscViewer, arg2::Ptr{Ptr{Libc.FILE}})::PetscErrorCode
end

function PetscViewerBinaryRead(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscViewerBinaryRead(arg1::PetscViewer, arg2::Ptr{Cvoid}, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::PetscDataType)::PetscErrorCode
end

function PetscViewerBinaryWrite(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscViewerBinaryWrite(arg1::PetscViewer, arg2::Ptr{Cvoid}, arg3::PetscInt, arg4::PetscDataType)::PetscErrorCode
end

function PetscViewerBinaryReadAll(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscViewerBinaryReadAll(arg1::PetscViewer, arg2::Ptr{Cvoid}, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscDataType)::PetscErrorCode
end

function PetscViewerBinaryWriteAll(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscViewerBinaryWriteAll(arg1::PetscViewer, arg2::Ptr{Cvoid}, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscDataType)::PetscErrorCode
end

function PetscViewerStringSetString(arg1, arg2, arg3)
    @ccall libpetsc.PetscViewerStringSetString(arg1::PetscViewer, arg2::Ptr{Cchar}, arg3::Csize_t)::PetscErrorCode
end

function PetscViewerStringGetStringRead(arg1, arg2, arg3)
    @ccall libpetsc.PetscViewerStringGetStringRead(arg1::PetscViewer, arg2::Ptr{Ptr{Cchar}}, arg3::Ptr{Csize_t})::PetscErrorCode
end

function PetscViewerStringSetOwnString(arg1)
    @ccall libpetsc.PetscViewerStringSetOwnString(arg1::PetscViewer)::PetscErrorCode
end

function PetscViewerDrawClear(arg1)
    @ccall libpetsc.PetscViewerDrawClear(arg1::PetscViewer)::PetscErrorCode
end

function PetscViewerDrawSetHold(arg1, arg2)
    @ccall libpetsc.PetscViewerDrawSetHold(arg1::PetscViewer, arg2::PetscBool)::PetscErrorCode
end

function PetscViewerDrawGetHold(arg1, arg2)
    @ccall libpetsc.PetscViewerDrawGetHold(arg1::PetscViewer, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscViewerDrawSetPause(arg1, arg2)
    @ccall libpetsc.PetscViewerDrawSetPause(arg1::PetscViewer, arg2::PetscReal)::PetscErrorCode
end

function PetscViewerDrawGetPause(arg1, arg2)
    @ccall libpetsc.PetscViewerDrawGetPause(arg1::PetscViewer, arg2::Ptr{PetscReal})::PetscErrorCode
end

function PetscViewerDrawSetInfo(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscViewerDrawSetInfo(arg1::PetscViewer, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Cint, arg5::Cint, arg6::Cint, arg7::Cint)::PetscErrorCode
end

function PetscViewerDrawResize(arg1, arg2, arg3)
    @ccall libpetsc.PetscViewerDrawResize(arg1::PetscViewer, arg2::Cint, arg3::Cint)::PetscErrorCode
end

function PetscViewerDrawSetBounds(arg1, arg2, arg3)
    @ccall libpetsc.PetscViewerDrawSetBounds(arg1::PetscViewer, arg2::PetscInt, arg3::Ptr{PetscReal})::PetscErrorCode
end

function PetscViewerDrawGetBounds(arg1, arg2, arg3)
    @ccall libpetsc.PetscViewerDrawGetBounds(arg1::PetscViewer, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{PetscReal}})::PetscErrorCode
end

function PetscViewerSocketSetConnection(arg1, arg2, arg3)
    @ccall libpetsc.PetscViewerSocketSetConnection(arg1::PetscViewer, arg2::Ptr{Cchar}, arg3::Cint)::PetscErrorCode
end

function PetscViewerBinarySkipInfo(arg1)
    @ccall libpetsc.PetscViewerBinarySkipInfo(arg1::PetscViewer)::PetscErrorCode
end

function PetscViewerBinarySetSkipInfo(arg1, arg2)
    @ccall libpetsc.PetscViewerBinarySetSkipInfo(arg1::PetscViewer, arg2::PetscBool)::PetscErrorCode
end

function PetscViewerBinaryGetSkipInfo(arg1, arg2)
    @ccall libpetsc.PetscViewerBinaryGetSkipInfo(arg1::PetscViewer, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscViewerBinarySetSkipOptions(arg1, arg2)
    @ccall libpetsc.PetscViewerBinarySetSkipOptions(arg1::PetscViewer, arg2::PetscBool)::PetscErrorCode
end

function PetscViewerBinaryGetSkipOptions(arg1, arg2)
    @ccall libpetsc.PetscViewerBinaryGetSkipOptions(arg1::PetscViewer, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscViewerBinarySetSkipHeader(arg1, arg2)
    @ccall libpetsc.PetscViewerBinarySetSkipHeader(arg1::PetscViewer, arg2::PetscBool)::PetscErrorCode
end

function PetscViewerBinaryGetSkipHeader(arg1, arg2)
    @ccall libpetsc.PetscViewerBinaryGetSkipHeader(arg1::PetscViewer, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscViewerBinaryReadStringArray(arg1, arg2)
    @ccall libpetsc.PetscViewerBinaryReadStringArray(arg1::PetscViewer, arg2::Ptr{Ptr{Ptr{Cchar}}})::PetscErrorCode
end

function PetscViewerBinaryWriteStringArray(arg1, arg2)
    @ccall libpetsc.PetscViewerBinaryWriteStringArray(arg1::PetscViewer, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscViewerFileSetName(arg1, arg2)
    @ccall libpetsc.PetscViewerFileSetName(arg1::PetscViewer, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscViewerFileGetName(arg1, arg2)
    @ccall libpetsc.PetscViewerFileGetName(arg1::PetscViewer, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscViewerVUGetPointer(arg1, arg2)
    @ccall libpetsc.PetscViewerVUGetPointer(arg1::PetscViewer, arg2::Ptr{Ptr{Libc.FILE}})::PetscErrorCode
end

function PetscViewerVUSetVecSeen(arg1, arg2)
    @ccall libpetsc.PetscViewerVUSetVecSeen(arg1::PetscViewer, arg2::PetscBool)::PetscErrorCode
end

function PetscViewerVUGetVecSeen(arg1, arg2)
    @ccall libpetsc.PetscViewerVUGetVecSeen(arg1::PetscViewer, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscViewerVUFlushDeferred(arg1)
    @ccall libpetsc.PetscViewerVUFlushDeferred(arg1::PetscViewer)::PetscErrorCode
end

function PetscViewerVUSetMode(viewer, mode)
    @ccall libpetsc.PetscViewerVUSetMode(viewer::PetscViewer, mode::PetscFileMode)::PetscErrorCode
end

function PetscViewerMathematicaInitializePackage()
    @ccall libpetsc.PetscViewerMathematicaInitializePackage()::PetscErrorCode
end

function PetscViewerMathematicaFinalizePackage()
    @ccall libpetsc.PetscViewerMathematicaFinalizePackage()::PetscErrorCode
end

function PetscViewerMathematicaGetName(arg1, arg2)
    @ccall libpetsc.PetscViewerMathematicaGetName(arg1::PetscViewer, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscViewerMathematicaSetName(arg1, arg2)
    @ccall libpetsc.PetscViewerMathematicaSetName(arg1::PetscViewer, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscViewerMathematicaClearName(arg1)
    @ccall libpetsc.PetscViewerMathematicaClearName(arg1::PetscViewer)::PetscErrorCode
end

function PetscViewerMathematicaSkipPackets(arg1, arg2)
    @ccall libpetsc.PetscViewerMathematicaSkipPackets(arg1::PetscViewer, arg2::Cint)::PetscErrorCode
end

function PetscViewerSiloGetName(arg1, arg2)
    @ccall libpetsc.PetscViewerSiloGetName(arg1::PetscViewer, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscViewerSiloSetName(arg1, arg2)
    @ccall libpetsc.PetscViewerSiloSetName(arg1::PetscViewer, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscViewerSiloClearName(arg1)
    @ccall libpetsc.PetscViewerSiloClearName(arg1::PetscViewer)::PetscErrorCode
end

function PetscViewerSiloGetMeshName(arg1, arg2)
    @ccall libpetsc.PetscViewerSiloGetMeshName(arg1::PetscViewer, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscViewerSiloSetMeshName(arg1, arg2)
    @ccall libpetsc.PetscViewerSiloSetMeshName(arg1::PetscViewer, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscViewerSiloClearMeshName(arg1)
    @ccall libpetsc.PetscViewerSiloClearMeshName(arg1::PetscViewer)::PetscErrorCode
end

@enum PetscViewerVTKFieldType::UInt32 begin
    PETSC_VTK_INVALID = 0
    PETSC_VTK_POINT_FIELD = 1
    PETSC_VTK_POINT_VECTOR_FIELD = 2
    PETSC_VTK_CELL_FIELD = 3
    PETSC_VTK_CELL_VECTOR_FIELD = 4
end

function PetscViewerVTKAddField(arg1, arg2, PetscViewerVTKWriteFunction, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscViewerVTKAddField(arg1::PetscViewer, arg2::PetscObject, PetscViewerVTKWriteFunction::Ptr{Cvoid}, arg4::PetscInt, arg5::PetscViewerVTKFieldType, arg6::PetscBool, arg7::PetscObject)::PetscErrorCode
end

function PetscViewerVTKGetDM(arg1, arg2)
    @ccall libpetsc.PetscViewerVTKGetDM(arg1::PetscViewer, arg2::Ptr{PetscObject})::PetscErrorCode
end

function PetscViewerVTKOpen(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscViewerVTKOpen(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::PetscFileMode, arg4::Ptr{PetscViewer})::PetscErrorCode
end

function PETSC_VIEWER_STDOUT_(arg1)
    @ccall libpetsc.PETSC_VIEWER_STDOUT_(arg1::MPI_Comm)::PetscViewer
end

function PetscViewerASCIIGetStdout(arg1, arg2)
    @ccall libpetsc.PetscViewerASCIIGetStdout(arg1::MPI_Comm, arg2::Ptr{PetscViewer})::PetscErrorCode
end

function PETSC_VIEWER_STDERR_(arg1)
    @ccall libpetsc.PETSC_VIEWER_STDERR_(arg1::MPI_Comm)::PetscViewer
end

function PetscViewerASCIIGetStderr(arg1, arg2)
    @ccall libpetsc.PetscViewerASCIIGetStderr(arg1::MPI_Comm, arg2::Ptr{PetscViewer})::PetscErrorCode
end

function PETSC_VIEWER_DRAW_(arg1)
    @ccall libpetsc.PETSC_VIEWER_DRAW_(arg1::MPI_Comm)::PetscViewer
end

function PETSC_VIEWER_SOCKET_(arg1)
    @ccall libpetsc.PETSC_VIEWER_SOCKET_(arg1::MPI_Comm)::PetscViewer
end

function PETSC_VIEWER_BINARY_(arg1)
    @ccall libpetsc.PETSC_VIEWER_BINARY_(arg1::MPI_Comm)::PetscViewer
end

function PETSC_VIEWER_MATLAB_(arg1)
    @ccall libpetsc.PETSC_VIEWER_MATLAB_(arg1::MPI_Comm)::PetscViewer
end

function PETSC_VIEWER_HDF5_(arg1)
    @ccall libpetsc.PETSC_VIEWER_HDF5_(arg1::MPI_Comm)::PetscViewer
end

function PETSC_VIEWER_GLVIS_(arg1)
    @ccall libpetsc.PETSC_VIEWER_GLVIS_(arg1::MPI_Comm)::PetscViewer
end

function PETSC_VIEWER_EXODUSII_(arg1)
    @ccall libpetsc.PETSC_VIEWER_EXODUSII_(arg1::MPI_Comm)::PetscViewer
end

function PetscViewerFlowControlStart(arg1, arg2, arg3)
    @ccall libpetsc.PetscViewerFlowControlStart(arg1::PetscViewer, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscViewerFlowControlStepMain(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscViewerFlowControlStepMain(arg1::PetscViewer, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::PetscInt)::PetscErrorCode
end

function PetscViewerFlowControlEndMain(arg1, arg2)
    @ccall libpetsc.PetscViewerFlowControlEndMain(arg1::PetscViewer, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscViewerFlowControlStepWorker(arg1, arg2, arg3)
    @ccall libpetsc.PetscViewerFlowControlStepWorker(arg1::PetscViewer, arg2::PetscMPIInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscViewerFlowControlEndWorker(arg1, arg2)
    @ccall libpetsc.PetscViewerFlowControlEndWorker(arg1::PetscViewer, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscViewerFlowControlStepMaster(viewer, i, mcnt, cnt)
    @ccall libpetsc.PetscViewerFlowControlStepMaster(viewer::PetscViewer, i::PetscInt, mcnt::Ptr{PetscInt}, cnt::PetscInt)::PetscErrorCode
end

function PetscViewerFlowControlEndMaster(viewer, mcnt)
    @ccall libpetsc.PetscViewerFlowControlEndMaster(viewer::PetscViewer, mcnt::Ptr{PetscInt})::PetscErrorCode
end

function PetscViewerMatlabPutArray(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscViewerMatlabPutArray(arg1::PetscViewer, arg2::Cint, arg3::Cint, arg4::Ptr{PetscScalar}, arg5::Ptr{Cchar})::PetscErrorCode
end

function PetscViewerMatlabGetArray(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscViewerMatlabGetArray(arg1::PetscViewer, arg2::Cint, arg3::Cint, arg4::Ptr{PetscScalar}, arg5::Ptr{Cchar})::PetscErrorCode
end

function PetscViewerMatlabPutVariable(arg1, arg2, arg3)
    @ccall libpetsc.PetscViewerMatlabPutVariable(arg1::PetscViewer, arg2::Ptr{Cchar}, arg3::Ptr{Cvoid})::PetscErrorCode
end

mutable struct _n_PetscViewers end

const PetscViewers = Ptr{_n_PetscViewers}

function PetscViewersCreate(arg1, arg2)
    @ccall libpetsc.PetscViewersCreate(arg1::MPI_Comm, arg2::Ptr{PetscViewers})::PetscErrorCode
end

function PetscViewersDestroy(arg1)
    @ccall libpetsc.PetscViewersDestroy(arg1::Ptr{PetscViewers})::PetscErrorCode
end

function PetscViewersGetViewer(arg1, arg2, arg3)
    @ccall libpetsc.PetscViewersGetViewer(arg1::PetscViewers, arg2::PetscInt, arg3::Ptr{PetscViewer})::PetscErrorCode
end

const PetscBT = Ptr{Cchar}

function PetscBTLength(m)
    @ccall libpetsc.PetscBTLength(m::PetscInt)::PetscInt
end

function PetscBTMemzero(m, array)
    @ccall libpetsc.PetscBTMemzero(m::PetscInt, array::PetscBT)::PetscErrorCode
end

function PetscBTDestroy(array)
    @ccall libpetsc.PetscBTDestroy(array::Ptr{PetscBT})::PetscErrorCode
end

function PetscBTLookup(array, index)
    @ccall libpetsc.PetscBTLookup(array::PetscBT, index::PetscInt)::Cchar
end

function PetscBTView(m, bt, viewer)
    @ccall libpetsc.PetscBTView(m::PetscInt, bt::PetscBT, viewer::PetscViewer)::PetscErrorCode
end

function PetscBTCreate(m, array)
    @ccall libpetsc.PetscBTCreate(m::PetscInt, array::Ptr{PetscBT})::PetscErrorCode
end

function PetscBTLookupSet(array, index)
    @ccall libpetsc.PetscBTLookupSet(array::PetscBT, index::PetscInt)::Cchar
end

function PetscBTSet(array, index)
    @ccall libpetsc.PetscBTSet(array::PetscBT, index::PetscInt)::PetscErrorCode
end

function PetscBTNegate(array, index)
    @ccall libpetsc.PetscBTNegate(array::PetscBT, index::PetscInt)::PetscErrorCode
end

function PetscBTLookupClear(array, index)
    @ccall libpetsc.PetscBTLookupClear(array::PetscBT, index::PetscInt)::Cchar
end

function PetscBTClear(array, index)
    @ccall libpetsc.PetscBTClear(array::PetscBT, index::PetscInt)::PetscErrorCode
end

mutable struct _n_PetscTable
    keytable::Ptr{PetscInt}
    table::Ptr{PetscInt}
    count::PetscInt
    tablesize::PetscInt
    head::PetscInt
    maxkey::PetscInt
    _n_PetscTable() = new()
end

const PetscTable = Ptr{_n_PetscTable}

const PetscTablePosition = Ptr{PetscInt}

function PetscHash(ta, x)
    @ccall libpetsc.PetscHash(ta::PetscTable, x::Culong)::Culong
end

function PetscHashStep(ta, x)
    @ccall libpetsc.PetscHashStep(ta::PetscTable, x::Culong)::Culong
end

function PetscTableCreate(arg1, arg2, arg3)
    @ccall libpetsc.PetscTableCreate(arg1::PetscInt, arg2::PetscInt, arg3::Ptr{PetscTable})::PetscErrorCode
end

function PetscTableCreateCopy(arg1, arg2)
    @ccall libpetsc.PetscTableCreateCopy(arg1::PetscTable, arg2::Ptr{PetscTable})::PetscErrorCode
end

function PetscTableDestroy(arg1)
    @ccall libpetsc.PetscTableDestroy(arg1::Ptr{PetscTable})::PetscErrorCode
end

function PetscTableGetCount(arg1, arg2)
    @ccall libpetsc.PetscTableGetCount(arg1::PetscTable, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscTableIsEmpty(arg1, arg2)
    @ccall libpetsc.PetscTableIsEmpty(arg1::PetscTable, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscTableAddExpand(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscTableAddExpand(arg1::PetscTable, arg2::PetscInt, arg3::PetscInt, arg4::InsertMode)::PetscErrorCode
end

function PetscTableAddCountExpand(arg1, arg2)
    @ccall libpetsc.PetscTableAddCountExpand(arg1::PetscTable, arg2::PetscInt)::PetscErrorCode
end

function PetscTableGetHeadPosition(arg1, arg2)
    @ccall libpetsc.PetscTableGetHeadPosition(arg1::PetscTable, arg2::Ptr{PetscTablePosition})::PetscErrorCode
end

function PetscTableGetNext(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscTableGetNext(arg1::PetscTable, arg2::Ptr{PetscTablePosition}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function PetscTableRemoveAll(arg1)
    @ccall libpetsc.PetscTableRemoveAll(arg1::PetscTable)::PetscErrorCode
end

function PetscTableAdd(ta, key, data, imode)
    @ccall libpetsc.PetscTableAdd(ta::PetscTable, key::PetscInt, data::PetscInt, imode::InsertMode)::PetscErrorCode
end

function PetscTableAddCount(ta, key)
    @ccall libpetsc.PetscTableAddCount(ta::PetscTable, key::PetscInt)::PetscErrorCode
end

function PetscTableFind(ta, key, data)
    @ccall libpetsc.PetscTableFind(ta::PetscTable, key::PetscInt, data::Ptr{PetscInt})::PetscErrorCode
end

mutable struct _p_PetscMatlabEngine end

const PetscMatlabEngine = Ptr{_p_PetscMatlabEngine}

function PetscMatlabEngineCreate(arg1, arg2, arg3)
    @ccall libpetsc.PetscMatlabEngineCreate(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{PetscMatlabEngine})::PetscErrorCode
end

function PetscMatlabEngineDestroy(arg1)
    @ccall libpetsc.PetscMatlabEngineDestroy(arg1::Ptr{PetscMatlabEngine})::PetscErrorCode
end

function PetscMatlabEngineGetOutput(arg1, arg2)
    @ccall libpetsc.PetscMatlabEngineGetOutput(arg1::PetscMatlabEngine, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscMatlabEnginePrintOutput(arg1, arg2)
    @ccall libpetsc.PetscMatlabEnginePrintOutput(arg1::PetscMatlabEngine, arg2::Ptr{Libc.FILE})::PetscErrorCode
end

function PetscMatlabEnginePut(arg1, arg2)
    @ccall libpetsc.PetscMatlabEnginePut(arg1::PetscMatlabEngine, arg2::PetscObject)::PetscErrorCode
end

function PetscMatlabEngineGet(arg1, arg2)
    @ccall libpetsc.PetscMatlabEngineGet(arg1::PetscMatlabEngine, arg2::PetscObject)::PetscErrorCode
end

function PetscMatlabEnginePutArray(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscMatlabEnginePutArray(arg1::PetscMatlabEngine, arg2::Cint, arg3::Cint, arg4::Ptr{PetscScalar}, arg5::Ptr{Cchar})::PetscErrorCode
end

function PetscMatlabEngineGetArray(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscMatlabEngineGetArray(arg1::PetscMatlabEngine, arg2::Cint, arg3::Cint, arg4::Ptr{PetscScalar}, arg5::Ptr{Cchar})::PetscErrorCode
end

function PETSC_MATLAB_ENGINE_(arg1)
    @ccall libpetsc.PETSC_MATLAB_ENGINE_(arg1::MPI_Comm)::PetscMatlabEngine
end

function PetscDrawInitializePackage()
    @ccall libpetsc.PetscDrawInitializePackage()::PetscErrorCode
end

function PetscDrawFinalizePackage()
    @ccall libpetsc.PetscDrawFinalizePackage()::PetscErrorCode
end

function PetscDrawRegister(arg1, arg2)
    @ccall libpetsc.PetscDrawRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscDrawGetType(arg1, arg2)
    @ccall libpetsc.PetscDrawGetType(arg1::PetscDraw, arg2::Ptr{PetscDrawType})::PetscErrorCode
end

function PetscDrawSetType(arg1, arg2)
    @ccall libpetsc.PetscDrawSetType(arg1::PetscDraw, arg2::PetscDrawType)::PetscErrorCode
end

function PetscDrawCreate(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscDrawCreate(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Cint, arg5::Cint, arg6::Cint, arg7::Cint, arg8::Ptr{PetscDraw})::PetscErrorCode
end

function PetscDrawSetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.PetscDrawSetOptionsPrefix(arg1::PetscDraw, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscDrawSetFromOptions(arg1)
    @ccall libpetsc.PetscDrawSetFromOptions(arg1::PetscDraw)::PetscErrorCode
end

function PetscDrawSetSave(arg1, arg2)
    @ccall libpetsc.PetscDrawSetSave(arg1::PetscDraw, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscDrawSetSaveMovie(arg1, arg2)
    @ccall libpetsc.PetscDrawSetSaveMovie(arg1::PetscDraw, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscDrawSetSaveFinalImage(arg1, arg2)
    @ccall libpetsc.PetscDrawSetSaveFinalImage(arg1::PetscDraw, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscDrawView(arg1, arg2)
    @ccall libpetsc.PetscDrawView(arg1::PetscDraw, arg2::PetscViewer)::PetscErrorCode
end

function PetscDrawViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.PetscDrawViewFromOptions(arg1::PetscDraw, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function PetscDrawRealToColor(value, min, max)
    @ccall libpetsc.PetscDrawRealToColor(value::PetscReal, min::PetscReal, max::PetscReal)::Cint
end

function PetscDrawOpenX(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscDrawOpenX(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Cint, arg5::Cint, arg6::Cint, arg7::Cint, arg8::Ptr{PetscDraw})::PetscErrorCode
end

function PetscDrawOpenImage(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDrawOpenImage(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Cint, arg4::Cint, arg5::Ptr{PetscDraw})::PetscErrorCode
end

function PetscDrawOpenNull(arg1, arg2)
    @ccall libpetsc.PetscDrawOpenNull(arg1::MPI_Comm, arg2::Ptr{PetscDraw})::PetscErrorCode
end

function PetscDrawDestroy(arg1)
    @ccall libpetsc.PetscDrawDestroy(arg1::Ptr{PetscDraw})::PetscErrorCode
end

function PetscDrawIsNull(arg1, arg2)
    @ccall libpetsc.PetscDrawIsNull(arg1::PetscDraw, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscDrawGetPopup(arg1, arg2)
    @ccall libpetsc.PetscDrawGetPopup(arg1::PetscDraw, arg2::Ptr{PetscDraw})::PetscErrorCode
end

function PetscDrawScalePopup(arg1, arg2, arg3)
    @ccall libpetsc.PetscDrawScalePopup(arg1::PetscDraw, arg2::PetscReal, arg3::PetscReal)::PetscErrorCode
end

function PetscDrawCheckResizedWindow(arg1)
    @ccall libpetsc.PetscDrawCheckResizedWindow(arg1::PetscDraw)::PetscErrorCode
end

function PetscDrawResizeWindow(arg1, arg2, arg3)
    @ccall libpetsc.PetscDrawResizeWindow(arg1::PetscDraw, arg2::Cint, arg3::Cint)::PetscErrorCode
end

function PetscDrawGetWindowSize(arg1, arg2, arg3)
    @ccall libpetsc.PetscDrawGetWindowSize(arg1::PetscDraw, arg2::Ptr{Cint}, arg3::Ptr{Cint})::PetscErrorCode
end

function PetscDrawPixelToCoordinate(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDrawPixelToCoordinate(arg1::PetscDraw, arg2::Cint, arg3::Cint, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal})::PetscErrorCode
end

function PetscDrawCoordinateToPixel(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDrawCoordinateToPixel(arg1::PetscDraw, arg2::PetscReal, arg3::PetscReal, arg4::Ptr{Cint}, arg5::Ptr{Cint})::PetscErrorCode
end

function PetscDrawIndicatorFunction(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscDrawIndicatorFunction(arg1::PetscDraw, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal, arg6::Cint, arg7::Ptr{Cvoid}, arg8::Ptr{Cvoid})::PetscErrorCode
end

function PetscDrawLine(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscDrawLine(arg1::PetscDraw, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal, arg6::Cint)::PetscErrorCode
end

function PetscDrawArrow(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscDrawArrow(arg1::PetscDraw, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal, arg6::Cint)::PetscErrorCode
end

function PetscDrawLineSetWidth(arg1, arg2)
    @ccall libpetsc.PetscDrawLineSetWidth(arg1::PetscDraw, arg2::PetscReal)::PetscErrorCode
end

function PetscDrawLineGetWidth(arg1, arg2)
    @ccall libpetsc.PetscDrawLineGetWidth(arg1::PetscDraw, arg2::Ptr{PetscReal})::PetscErrorCode
end

@enum PetscDrawMarkerType::UInt32 begin
    PETSC_DRAW_MARKER_CROSS = 0
    PETSC_DRAW_MARKER_POINT = 1
    PETSC_DRAW_MARKER_PLUS = 2
    PETSC_DRAW_MARKER_CIRCLE = 3
end

function PetscDrawMarker(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDrawMarker(arg1::PetscDraw, arg2::PetscReal, arg3::PetscReal, arg4::Cint)::PetscErrorCode
end

function PetscDrawSetMarkerType(arg1, arg2)
    @ccall libpetsc.PetscDrawSetMarkerType(arg1::PetscDraw, arg2::PetscDrawMarkerType)::PetscErrorCode
end

function PetscDrawGetMarkerType(arg1, arg2)
    @ccall libpetsc.PetscDrawGetMarkerType(arg1::PetscDraw, arg2::Ptr{PetscDrawMarkerType})::PetscErrorCode
end

function PetscDrawPoint(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDrawPoint(arg1::PetscDraw, arg2::PetscReal, arg3::PetscReal, arg4::Cint)::PetscErrorCode
end

function PetscDrawPointPixel(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDrawPointPixel(arg1::PetscDraw, arg2::Cint, arg3::Cint, arg4::Cint)::PetscErrorCode
end

function PetscDrawPointSetSize(arg1, arg2)
    @ccall libpetsc.PetscDrawPointSetSize(arg1::PetscDraw, arg2::PetscReal)::PetscErrorCode
end

function PetscDrawRectangle(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.PetscDrawRectangle(arg1::PetscDraw, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal, arg6::Cint, arg7::Cint, arg8::Cint, arg9::Cint)::PetscErrorCode
end

function PetscDrawTriangle(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.PetscDrawTriangle(arg1::PetscDraw, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal, arg6::PetscReal, arg7::PetscReal, arg8::Cint, arg9::Cint, arg10::Cint)::PetscErrorCode
end

function PetscDrawEllipse(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscDrawEllipse(arg1::PetscDraw, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal, arg6::Cint)::PetscErrorCode
end

function PetscDrawTensorContourPatch(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscDrawTensorContourPatch(arg1::PetscDraw, arg2::Cint, arg3::Cint, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal}, arg6::PetscReal, arg7::PetscReal, arg8::Ptr{PetscReal})::PetscErrorCode
end

function PetscDrawTensorContour(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscDrawTensorContour(arg1::PetscDraw, arg2::Cint, arg3::Cint, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal})::PetscErrorCode
end

function PetscDrawString(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDrawString(arg1::PetscDraw, arg2::PetscReal, arg3::PetscReal, arg4::Cint, arg5::Ptr{Cchar})::PetscErrorCode
end

function PetscDrawStringCentered(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDrawStringCentered(arg1::PetscDraw, arg2::PetscReal, arg3::PetscReal, arg4::Cint, arg5::Ptr{Cchar})::PetscErrorCode
end

function PetscDrawStringBoxed(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscDrawStringBoxed(arg1::PetscDraw, arg2::PetscReal, arg3::PetscReal, arg4::Cint, arg5::Cint, arg6::Ptr{Cchar}, arg7::Ptr{PetscReal}, arg8::Ptr{PetscReal})::PetscErrorCode
end

function PetscDrawStringVertical(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDrawStringVertical(arg1::PetscDraw, arg2::PetscReal, arg3::PetscReal, arg4::Cint, arg5::Ptr{Cchar})::PetscErrorCode
end

function PetscDrawStringSetSize(arg1, arg2, arg3)
    @ccall libpetsc.PetscDrawStringSetSize(arg1::PetscDraw, arg2::PetscReal, arg3::PetscReal)::PetscErrorCode
end

function PetscDrawStringGetSize(arg1, arg2, arg3)
    @ccall libpetsc.PetscDrawStringGetSize(arg1::PetscDraw, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal})::PetscErrorCode
end

function PetscDrawSetViewPort(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDrawSetViewPort(arg1::PetscDraw, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal)::PetscErrorCode
end

function PetscDrawGetViewPort(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDrawGetViewPort(arg1::PetscDraw, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal})::PetscErrorCode
end

function PetscDrawSplitViewPort(arg1)
    @ccall libpetsc.PetscDrawSplitViewPort(arg1::PetscDraw)::PetscErrorCode
end

function PetscDrawSetCoordinates(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDrawSetCoordinates(arg1::PetscDraw, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal)::PetscErrorCode
end

function PetscDrawGetCoordinates(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDrawGetCoordinates(arg1::PetscDraw, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal})::PetscErrorCode
end

function PetscDrawSetTitle(arg1, arg2)
    @ccall libpetsc.PetscDrawSetTitle(arg1::PetscDraw, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscDrawAppendTitle(arg1, arg2)
    @ccall libpetsc.PetscDrawAppendTitle(arg1::PetscDraw, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscDrawGetTitle(arg1, arg2)
    @ccall libpetsc.PetscDrawGetTitle(arg1::PetscDraw, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscDrawSetPause(arg1, arg2)
    @ccall libpetsc.PetscDrawSetPause(arg1::PetscDraw, arg2::PetscReal)::PetscErrorCode
end

function PetscDrawGetPause(arg1, arg2)
    @ccall libpetsc.PetscDrawGetPause(arg1::PetscDraw, arg2::Ptr{PetscReal})::PetscErrorCode
end

function PetscDrawPause(arg1)
    @ccall libpetsc.PetscDrawPause(arg1::PetscDraw)::PetscErrorCode
end

function PetscDrawSetDoubleBuffer(arg1)
    @ccall libpetsc.PetscDrawSetDoubleBuffer(arg1::PetscDraw)::PetscErrorCode
end

function PetscDrawClear(arg1)
    @ccall libpetsc.PetscDrawClear(arg1::PetscDraw)::PetscErrorCode
end

function PetscDrawFlush(arg1)
    @ccall libpetsc.PetscDrawFlush(arg1::PetscDraw)::PetscErrorCode
end

function PetscDrawSave(arg1)
    @ccall libpetsc.PetscDrawSave(arg1::PetscDraw)::PetscErrorCode
end

function PetscDrawSaveMovie(arg1)
    @ccall libpetsc.PetscDrawSaveMovie(arg1::PetscDraw)::PetscErrorCode
end

function PetscDrawBOP(arg1)
    @ccall libpetsc.PetscDrawBOP(arg1::PetscDraw)::PetscErrorCode
end

function PetscDrawEOP(arg1)
    @ccall libpetsc.PetscDrawEOP(arg1::PetscDraw)::PetscErrorCode
end

function PetscDrawSetDisplay(arg1, arg2)
    @ccall libpetsc.PetscDrawSetDisplay(arg1::PetscDraw, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscDrawGetSingleton(arg1, arg2)
    @ccall libpetsc.PetscDrawGetSingleton(arg1::PetscDraw, arg2::Ptr{PetscDraw})::PetscErrorCode
end

function PetscDrawRestoreSingleton(arg1, arg2)
    @ccall libpetsc.PetscDrawRestoreSingleton(arg1::PetscDraw, arg2::Ptr{PetscDraw})::PetscErrorCode
end

function PetscDrawGetCurrentPoint(arg1, arg2, arg3)
    @ccall libpetsc.PetscDrawGetCurrentPoint(arg1::PetscDraw, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal})::PetscErrorCode
end

function PetscDrawSetCurrentPoint(arg1, arg2, arg3)
    @ccall libpetsc.PetscDrawSetCurrentPoint(arg1::PetscDraw, arg2::PetscReal, arg3::PetscReal)::PetscErrorCode
end

function PetscDrawPushCurrentPoint(arg1, arg2, arg3)
    @ccall libpetsc.PetscDrawPushCurrentPoint(arg1::PetscDraw, arg2::PetscReal, arg3::PetscReal)::PetscErrorCode
end

function PetscDrawPopCurrentPoint(arg1)
    @ccall libpetsc.PetscDrawPopCurrentPoint(arg1::PetscDraw)::PetscErrorCode
end

function PetscDrawGetBoundingBox(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDrawGetBoundingBox(arg1::PetscDraw, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal})::PetscErrorCode
end

@enum PetscDrawButton::UInt32 begin
    PETSC_BUTTON_NONE = 0
    PETSC_BUTTON_LEFT = 1
    PETSC_BUTTON_CENTER = 2
    PETSC_BUTTON_RIGHT = 3
    PETSC_BUTTON_WHEEL_UP = 4
    PETSC_BUTTON_WHEEL_DOWN = 5
    PETSC_BUTTON_LEFT_SHIFT = 6
    PETSC_BUTTON_CENTER_SHIFT = 7
    PETSC_BUTTON_RIGHT_SHIFT = 8
end

function PetscDrawGetMouseButton(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscDrawGetMouseButton(arg1::PetscDraw, arg2::Ptr{PetscDrawButton}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal})::PetscErrorCode
end

function PetscDrawZoom(arg1, arg2, arg3)
    @ccall libpetsc.PetscDrawZoom(arg1::PetscDraw, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

mutable struct PetscDrawViewPorts
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
    PetscDrawViewPorts() = new()
end

function PetscDrawViewPortsCreate(arg1, arg2, arg3)
    @ccall libpetsc.PetscDrawViewPortsCreate(arg1::PetscDraw, arg2::PetscInt, arg3::Ptr{Ptr{PetscDrawViewPorts}})::PetscErrorCode
end

function PetscDrawViewPortsCreateRect(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDrawViewPortsCreateRect(arg1::PetscDraw, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Ptr{PetscDrawViewPorts}})::PetscErrorCode
end

function PetscDrawViewPortsDestroy(arg1)
    @ccall libpetsc.PetscDrawViewPortsDestroy(arg1::Ptr{PetscDrawViewPorts})::PetscErrorCode
end

function PetscDrawViewPortsSet(arg1, arg2)
    @ccall libpetsc.PetscDrawViewPortsSet(arg1::Ptr{PetscDrawViewPorts}, arg2::PetscInt)::PetscErrorCode
end

function PetscDrawAxisCreate(arg1, arg2)
    @ccall libpetsc.PetscDrawAxisCreate(arg1::PetscDraw, arg2::Ptr{PetscDrawAxis})::PetscErrorCode
end

function PetscDrawAxisDestroy(arg1)
    @ccall libpetsc.PetscDrawAxisDestroy(arg1::Ptr{PetscDrawAxis})::PetscErrorCode
end

function PetscDrawAxisDraw(arg1)
    @ccall libpetsc.PetscDrawAxisDraw(arg1::PetscDrawAxis)::PetscErrorCode
end

function PetscDrawAxisSetLimits(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDrawAxisSetLimits(arg1::PetscDrawAxis, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal)::PetscErrorCode
end

function PetscDrawAxisGetLimits(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDrawAxisGetLimits(arg1::PetscDrawAxis, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal})::PetscErrorCode
end

function PetscDrawAxisSetHoldLimits(arg1, arg2)
    @ccall libpetsc.PetscDrawAxisSetHoldLimits(arg1::PetscDrawAxis, arg2::PetscBool)::PetscErrorCode
end

function PetscDrawAxisSetColors(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDrawAxisSetColors(arg1::PetscDrawAxis, arg2::Cint, arg3::Cint, arg4::Cint)::PetscErrorCode
end

function PetscDrawAxisSetLabels(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDrawAxisSetLabels(arg1::PetscDrawAxis, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar})::PetscErrorCode
end

function PetscDrawLGCreate(arg1, arg2, arg3)
    @ccall libpetsc.PetscDrawLGCreate(arg1::PetscDraw, arg2::PetscInt, arg3::Ptr{PetscDrawLG})::PetscErrorCode
end

function PetscDrawLGDestroy(arg1)
    @ccall libpetsc.PetscDrawLGDestroy(arg1::Ptr{PetscDrawLG})::PetscErrorCode
end

function PetscDrawLGAddPoint(arg1, arg2, arg3)
    @ccall libpetsc.PetscDrawLGAddPoint(arg1::PetscDrawLG, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal})::PetscErrorCode
end

function PetscDrawLGAddCommonPoint(arg1, arg2, arg3)
    @ccall libpetsc.PetscDrawLGAddCommonPoint(arg1::PetscDrawLG, arg2::PetscReal, arg3::Ptr{PetscReal})::PetscErrorCode
end

function PetscDrawLGAddPoints(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDrawLGAddPoints(arg1::PetscDrawLG, arg2::PetscInt, arg3::Ptr{Ptr{PetscReal}}, arg4::Ptr{Ptr{PetscReal}})::PetscErrorCode
end

function PetscDrawLGDraw(arg1)
    @ccall libpetsc.PetscDrawLGDraw(arg1::PetscDrawLG)::PetscErrorCode
end

function PetscDrawLGSave(arg1)
    @ccall libpetsc.PetscDrawLGSave(arg1::PetscDrawLG)::PetscErrorCode
end

function PetscDrawLGView(arg1, arg2)
    @ccall libpetsc.PetscDrawLGView(arg1::PetscDrawLG, arg2::PetscViewer)::PetscErrorCode
end

function PetscDrawLGReset(arg1)
    @ccall libpetsc.PetscDrawLGReset(arg1::PetscDrawLG)::PetscErrorCode
end

function PetscDrawLGSetDimension(arg1, arg2)
    @ccall libpetsc.PetscDrawLGSetDimension(arg1::PetscDrawLG, arg2::PetscInt)::PetscErrorCode
end

function PetscDrawLGGetDimension(arg1, arg2)
    @ccall libpetsc.PetscDrawLGGetDimension(arg1::PetscDrawLG, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscDrawLGSetLegend(arg1, arg2)
    @ccall libpetsc.PetscDrawLGSetLegend(arg1::PetscDrawLG, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscDrawLGGetAxis(arg1, arg2)
    @ccall libpetsc.PetscDrawLGGetAxis(arg1::PetscDrawLG, arg2::Ptr{PetscDrawAxis})::PetscErrorCode
end

function PetscDrawLGGetDraw(arg1, arg2)
    @ccall libpetsc.PetscDrawLGGetDraw(arg1::PetscDrawLG, arg2::Ptr{PetscDraw})::PetscErrorCode
end

function PetscDrawLGSetUseMarkers(arg1, arg2)
    @ccall libpetsc.PetscDrawLGSetUseMarkers(arg1::PetscDrawLG, arg2::PetscBool)::PetscErrorCode
end

function PetscDrawLGSetLimits(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDrawLGSetLimits(arg1::PetscDrawLG, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal)::PetscErrorCode
end

function PetscDrawLGSetColors(arg1, arg2)
    @ccall libpetsc.PetscDrawLGSetColors(arg1::PetscDrawLG, arg2::Ptr{Cint})::PetscErrorCode
end

function PetscDrawLGSetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.PetscDrawLGSetOptionsPrefix(arg1::PetscDrawLG, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscDrawLGSetFromOptions(arg1)
    @ccall libpetsc.PetscDrawLGSetFromOptions(arg1::PetscDrawLG)::PetscErrorCode
end

function PetscDrawSPCreate(arg1, arg2, arg3)
    @ccall libpetsc.PetscDrawSPCreate(arg1::PetscDraw, arg2::Cint, arg3::Ptr{PetscDrawSP})::PetscErrorCode
end

function PetscDrawSPDestroy(arg1)
    @ccall libpetsc.PetscDrawSPDestroy(arg1::Ptr{PetscDrawSP})::PetscErrorCode
end

function PetscDrawSPAddPoint(arg1, arg2, arg3)
    @ccall libpetsc.PetscDrawSPAddPoint(arg1::PetscDrawSP, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal})::PetscErrorCode
end

function PetscDrawSPAddPoints(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDrawSPAddPoints(arg1::PetscDrawSP, arg2::Cint, arg3::Ptr{Ptr{PetscReal}}, arg4::Ptr{Ptr{PetscReal}})::PetscErrorCode
end

function PetscDrawSPDraw(arg1, arg2)
    @ccall libpetsc.PetscDrawSPDraw(arg1::PetscDrawSP, arg2::PetscBool)::PetscErrorCode
end

function PetscDrawSPSave(arg1)
    @ccall libpetsc.PetscDrawSPSave(arg1::PetscDrawSP)::PetscErrorCode
end

function PetscDrawSPReset(arg1)
    @ccall libpetsc.PetscDrawSPReset(arg1::PetscDrawSP)::PetscErrorCode
end

function PetscDrawSPSetDimension(arg1, arg2)
    @ccall libpetsc.PetscDrawSPSetDimension(arg1::PetscDrawSP, arg2::Cint)::PetscErrorCode
end

function PetscDrawSPGetAxis(arg1, arg2)
    @ccall libpetsc.PetscDrawSPGetAxis(arg1::PetscDrawSP, arg2::Ptr{PetscDrawAxis})::PetscErrorCode
end

function PetscDrawSPGetDraw(arg1, arg2)
    @ccall libpetsc.PetscDrawSPGetDraw(arg1::PetscDrawSP, arg2::Ptr{PetscDraw})::PetscErrorCode
end

function PetscDrawSPSetLimits(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDrawSPSetLimits(arg1::PetscDrawSP, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal)::PetscErrorCode
end

function PetscDrawLGSPDraw(arg1, arg2)
    @ccall libpetsc.PetscDrawLGSPDraw(arg1::PetscDrawLG, arg2::PetscDrawSP)::PetscErrorCode
end

function PetscDrawHGCreate(arg1, arg2, arg3)
    @ccall libpetsc.PetscDrawHGCreate(arg1::PetscDraw, arg2::Cint, arg3::Ptr{PetscDrawHG})::PetscErrorCode
end

function PetscDrawHGDestroy(arg1)
    @ccall libpetsc.PetscDrawHGDestroy(arg1::Ptr{PetscDrawHG})::PetscErrorCode
end

function PetscDrawHGAddValue(arg1, arg2)
    @ccall libpetsc.PetscDrawHGAddValue(arg1::PetscDrawHG, arg2::PetscReal)::PetscErrorCode
end

function PetscDrawHGDraw(arg1)
    @ccall libpetsc.PetscDrawHGDraw(arg1::PetscDrawHG)::PetscErrorCode
end

function PetscDrawHGSave(arg1)
    @ccall libpetsc.PetscDrawHGSave(arg1::PetscDrawHG)::PetscErrorCode
end

function PetscDrawHGView(arg1, arg2)
    @ccall libpetsc.PetscDrawHGView(arg1::PetscDrawHG, arg2::PetscViewer)::PetscErrorCode
end

function PetscDrawHGReset(arg1)
    @ccall libpetsc.PetscDrawHGReset(arg1::PetscDrawHG)::PetscErrorCode
end

function PetscDrawHGGetAxis(arg1, arg2)
    @ccall libpetsc.PetscDrawHGGetAxis(arg1::PetscDrawHG, arg2::Ptr{PetscDrawAxis})::PetscErrorCode
end

function PetscDrawHGGetDraw(arg1, arg2)
    @ccall libpetsc.PetscDrawHGGetDraw(arg1::PetscDrawHG, arg2::Ptr{PetscDraw})::PetscErrorCode
end

function PetscDrawHGSetLimits(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDrawHGSetLimits(arg1::PetscDrawHG, arg2::PetscReal, arg3::PetscReal, arg4::Cint, arg5::Cint)::PetscErrorCode
end

function PetscDrawHGSetNumberBins(arg1, arg2)
    @ccall libpetsc.PetscDrawHGSetNumberBins(arg1::PetscDrawHG, arg2::Cint)::PetscErrorCode
end

function PetscDrawHGSetColor(arg1, arg2)
    @ccall libpetsc.PetscDrawHGSetColor(arg1::PetscDrawHG, arg2::Cint)::PetscErrorCode
end

function PetscDrawHGCalcStats(arg1, arg2)
    @ccall libpetsc.PetscDrawHGCalcStats(arg1::PetscDrawHG, arg2::PetscBool)::PetscErrorCode
end

function PetscDrawHGIntegerBins(arg1, arg2)
    @ccall libpetsc.PetscDrawHGIntegerBins(arg1::PetscDrawHG, arg2::PetscBool)::PetscErrorCode
end

function PetscDrawBarCreate(arg1, arg2)
    @ccall libpetsc.PetscDrawBarCreate(arg1::PetscDraw, arg2::Ptr{PetscDrawBar})::PetscErrorCode
end

function PetscDrawBarSetData(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDrawBarSetData(arg1::PetscDrawBar, arg2::PetscInt, arg3::Ptr{PetscReal}, arg4::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscDrawBarDestroy(arg1)
    @ccall libpetsc.PetscDrawBarDestroy(arg1::Ptr{PetscDrawBar})::PetscErrorCode
end

function PetscDrawBarDraw(arg1)
    @ccall libpetsc.PetscDrawBarDraw(arg1::PetscDrawBar)::PetscErrorCode
end

function PetscDrawBarSave(arg1)
    @ccall libpetsc.PetscDrawBarSave(arg1::PetscDrawBar)::PetscErrorCode
end

function PetscDrawBarSetColor(arg1, arg2)
    @ccall libpetsc.PetscDrawBarSetColor(arg1::PetscDrawBar, arg2::Cint)::PetscErrorCode
end

function PetscDrawBarSetLimits(arg1, arg2, arg3)
    @ccall libpetsc.PetscDrawBarSetLimits(arg1::PetscDrawBar, arg2::PetscReal, arg3::PetscReal)::PetscErrorCode
end

function PetscDrawBarSort(arg1, arg2, arg3)
    @ccall libpetsc.PetscDrawBarSort(arg1::PetscDrawBar, arg2::PetscBool, arg3::PetscReal)::PetscErrorCode
end

function PetscDrawBarSetFromOptions(arg1)
    @ccall libpetsc.PetscDrawBarSetFromOptions(arg1::PetscDrawBar)::PetscErrorCode
end

function PetscDrawBarGetAxis(arg1, arg2)
    @ccall libpetsc.PetscDrawBarGetAxis(arg1::PetscDrawBar, arg2::Ptr{PetscDrawAxis})::PetscErrorCode
end

function PetscDrawBarGetDraw(arg1, arg2)
    @ccall libpetsc.PetscDrawBarGetDraw(arg1::PetscDrawBar, arg2::Ptr{PetscDraw})::PetscErrorCode
end

function PetscDrawUtilitySetCmap(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDrawUtilitySetCmap(arg1::Ptr{Cchar}, arg2::Cint, arg3::Ptr{Cuchar}, arg4::Ptr{Cuchar}, arg5::Ptr{Cuchar})::PetscErrorCode
end

function PetscDrawUtilitySetGamma(arg1)
    @ccall libpetsc.PetscDrawUtilitySetGamma(arg1::PetscReal)::PetscErrorCode
end

mutable struct _p_PetscSF end

const PetscSF = Ptr{_p_PetscSF}

const PetscSFType = Ptr{Cchar}

mutable struct PetscSFNode
    rank::PetscInt
    index::PetscInt
    PetscSFNode() = new()
end

const VecScatter = PetscSF

const VecScatterType = PetscSFType

mutable struct _p_PetscSection end

const PetscSection = Ptr{_p_PetscSection}

mutable struct _p_PetscSectionSym end

const PetscSectionSym = Ptr{_p_PetscSectionSym}

const PetscSectionSymType = Ptr{Cchar}

mutable struct _p_IS end

const IS = Ptr{_p_IS}

mutable struct _p_ISLocalToGlobalMapping end

const ISLocalToGlobalMapping = Ptr{_p_ISLocalToGlobalMapping}

mutable struct _n_ISColoring end

const ISColoring = Ptr{_n_ISColoring}

mutable struct _n_PetscLayout
    comm::MPI_Comm
    size::PetscMPIInt
    n::PetscInt
    N::PetscInt
    rstart::PetscInt
    rend::PetscInt
    range::Ptr{PetscInt}
    range_alloc::PetscBool
    bs::PetscInt
    refcnt::PetscInt
    mapping::ISLocalToGlobalMapping
    setupcalled::PetscBool
    oldn::PetscInt
    oldN::PetscInt
    oldbs::PetscInt
    _n_PetscLayout() = new()
end

const PetscLayout = Ptr{_n_PetscLayout}

function ISInitializePackage()
    @ccall libpetsc.ISInitializePackage()::PetscErrorCode
end

const ISType = Ptr{Cchar}

function ISSetType(arg1, arg2)
    @ccall libpetsc.ISSetType(arg1::IS, arg2::ISType)::PetscErrorCode
end

function ISGetType(arg1, arg2)
    @ccall libpetsc.ISGetType(arg1::IS, arg2::Ptr{ISType})::PetscErrorCode
end

function ISRegister(arg1, arg2)
    @ccall libpetsc.ISRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function ISCreate(arg1, arg2)
    @ccall libpetsc.ISCreate(arg1::MPI_Comm, arg2::Ptr{IS})::PetscErrorCode
end

function ISCreateGeneral(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.ISCreateGeneral(arg1::MPI_Comm, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::PetscCopyMode, arg5::Ptr{IS})::PetscErrorCode
end

function ISGeneralSetIndices(arg1, arg2, arg3, arg4)
    @ccall libpetsc.ISGeneralSetIndices(arg1::IS, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::PetscCopyMode)::PetscErrorCode
end

function ISCreateBlock(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.ISCreateBlock(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::PetscCopyMode, arg6::Ptr{IS})::PetscErrorCode
end

function ISBlockSetIndices(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.ISBlockSetIndices(arg1::IS, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::PetscCopyMode)::PetscErrorCode
end

function ISCreateStride(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.ISCreateStride(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{IS})::PetscErrorCode
end

function ISStrideSetStride(arg1, arg2, arg3, arg4)
    @ccall libpetsc.ISStrideSetStride(arg1::IS, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt)::PetscErrorCode
end

function ISDestroy(arg1)
    @ccall libpetsc.ISDestroy(arg1::Ptr{IS})::PetscErrorCode
end

function ISSetPermutation(arg1)
    @ccall libpetsc.ISSetPermutation(arg1::IS)::PetscErrorCode
end

function ISPermutation(arg1, arg2)
    @ccall libpetsc.ISPermutation(arg1::IS, arg2::Ptr{PetscBool})::PetscErrorCode
end

function ISSetIdentity(arg1)
    @ccall libpetsc.ISSetIdentity(arg1::IS)::PetscErrorCode
end

function ISIdentity(arg1, arg2)
    @ccall libpetsc.ISIdentity(arg1::IS, arg2::Ptr{PetscBool})::PetscErrorCode
end

function ISContiguousLocal(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.ISContiguousLocal(arg1::IS, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::Ptr{PetscBool})::PetscErrorCode
end

@enum ISInfo::Int32 begin
    IS_INFO_MIN = -1
    IS_SORTED = 0
    IS_UNIQUE = 1
    IS_PERMUTATION = 2
    IS_INTERVAL = 3
    IS_IDENTITY = 4
    IS_INFO_MAX = 5
end

@enum ISInfoType::UInt32 begin
    IS_LOCAL = 0
    IS_GLOBAL = 1
end

function ISSetInfo(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.ISSetInfo(arg1::IS, arg2::ISInfo, arg3::ISInfoType, arg4::PetscBool, arg5::PetscBool)::PetscErrorCode
end

function ISGetInfo(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.ISGetInfo(arg1::IS, arg2::ISInfo, arg3::ISInfoType, arg4::PetscBool, arg5::Ptr{PetscBool})::PetscErrorCode
end

function ISClearInfoCache(arg1, arg2)
    @ccall libpetsc.ISClearInfoCache(arg1::IS, arg2::PetscBool)::PetscErrorCode
end

function ISGetIndices(arg1, arg2)
    @ccall libpetsc.ISGetIndices(arg1::IS, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function ISRestoreIndices(arg1, arg2)
    @ccall libpetsc.ISRestoreIndices(arg1::IS, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function ISGetTotalIndices(arg1, arg2)
    @ccall libpetsc.ISGetTotalIndices(arg1::IS, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function ISRestoreTotalIndices(arg1, arg2)
    @ccall libpetsc.ISRestoreTotalIndices(arg1::IS, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function ISGetNonlocalIndices(arg1, arg2)
    @ccall libpetsc.ISGetNonlocalIndices(arg1::IS, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function ISRestoreNonlocalIndices(arg1, arg2)
    @ccall libpetsc.ISRestoreNonlocalIndices(arg1::IS, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function ISGetNonlocalIS(arg1, is)
    @ccall libpetsc.ISGetNonlocalIS(arg1::IS, is::Ptr{IS})::PetscErrorCode
end

function ISRestoreNonlocalIS(arg1, is)
    @ccall libpetsc.ISRestoreNonlocalIS(arg1::IS, is::Ptr{IS})::PetscErrorCode
end

function ISGetSize(arg1, arg2)
    @ccall libpetsc.ISGetSize(arg1::IS, arg2::Ptr{PetscInt})::PetscErrorCode
end

function ISGetLocalSize(arg1, arg2)
    @ccall libpetsc.ISGetLocalSize(arg1::IS, arg2::Ptr{PetscInt})::PetscErrorCode
end

function ISInvertPermutation(arg1, arg2, arg3)
    @ccall libpetsc.ISInvertPermutation(arg1::IS, arg2::PetscInt, arg3::Ptr{IS})::PetscErrorCode
end

function ISView(arg1, arg2)
    @ccall libpetsc.ISView(arg1::IS, arg2::PetscViewer)::PetscErrorCode
end

function ISViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.ISViewFromOptions(arg1::IS, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function ISLoad(arg1, arg2)
    @ccall libpetsc.ISLoad(arg1::IS, arg2::PetscViewer)::PetscErrorCode
end

function ISEqual(arg1, arg2, arg3)
    @ccall libpetsc.ISEqual(arg1::IS, arg2::IS, arg3::Ptr{PetscBool})::PetscErrorCode
end

function ISEqualUnsorted(arg1, arg2, arg3)
    @ccall libpetsc.ISEqualUnsorted(arg1::IS, arg2::IS, arg3::Ptr{PetscBool})::PetscErrorCode
end

function ISSort(arg1)
    @ccall libpetsc.ISSort(arg1::IS)::PetscErrorCode
end

function ISSortRemoveDups(arg1)
    @ccall libpetsc.ISSortRemoveDups(arg1::IS)::PetscErrorCode
end

function ISSorted(arg1, arg2)
    @ccall libpetsc.ISSorted(arg1::IS, arg2::Ptr{PetscBool})::PetscErrorCode
end

function ISDifference(arg1, arg2, arg3)
    @ccall libpetsc.ISDifference(arg1::IS, arg2::IS, arg3::Ptr{IS})::PetscErrorCode
end

function ISSum(arg1, arg2, arg3)
    @ccall libpetsc.ISSum(arg1::IS, arg2::IS, arg3::Ptr{IS})::PetscErrorCode
end

function ISExpand(arg1, arg2, arg3)
    @ccall libpetsc.ISExpand(arg1::IS, arg2::IS, arg3::Ptr{IS})::PetscErrorCode
end

function ISIntersect(arg1, arg2, arg3)
    @ccall libpetsc.ISIntersect(arg1::IS, arg2::IS, arg3::Ptr{IS})::PetscErrorCode
end

function ISGetMinMax(arg1, arg2, arg3)
    @ccall libpetsc.ISGetMinMax(arg1::IS, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function ISLocate(arg1, arg2, arg3)
    @ccall libpetsc.ISLocate(arg1::IS, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function ISGetPointRange(arg1, arg2, arg3, arg4)
    @ccall libpetsc.ISGetPointRange(arg1::IS, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function ISRestorePointRange(arg1, arg2, arg3, arg4)
    @ccall libpetsc.ISRestorePointRange(arg1::IS, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function ISGetPointSubrange(arg1, arg2, arg3, arg4)
    @ccall libpetsc.ISGetPointSubrange(arg1::IS, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt})::PetscErrorCode
end

function ISBlockGetIndices(arg1, arg2)
    @ccall libpetsc.ISBlockGetIndices(arg1::IS, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function ISBlockRestoreIndices(arg1, arg2)
    @ccall libpetsc.ISBlockRestoreIndices(arg1::IS, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function ISBlockGetLocalSize(arg1, arg2)
    @ccall libpetsc.ISBlockGetLocalSize(arg1::IS, arg2::Ptr{PetscInt})::PetscErrorCode
end

function ISBlockGetSize(arg1, arg2)
    @ccall libpetsc.ISBlockGetSize(arg1::IS, arg2::Ptr{PetscInt})::PetscErrorCode
end

function ISGetBlockSize(arg1, arg2)
    @ccall libpetsc.ISGetBlockSize(arg1::IS, arg2::Ptr{PetscInt})::PetscErrorCode
end

function ISSetBlockSize(arg1, arg2)
    @ccall libpetsc.ISSetBlockSize(arg1::IS, arg2::PetscInt)::PetscErrorCode
end

function ISStrideGetInfo(arg1, arg2, arg3)
    @ccall libpetsc.ISStrideGetInfo(arg1::IS, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function ISToGeneral(arg1)
    @ccall libpetsc.ISToGeneral(arg1::IS)::PetscErrorCode
end

function ISDuplicate(arg1, arg2)
    @ccall libpetsc.ISDuplicate(arg1::IS, arg2::Ptr{IS})::PetscErrorCode
end

function ISCopy(arg1, arg2)
    @ccall libpetsc.ISCopy(arg1::IS, arg2::IS)::PetscErrorCode
end

function ISAllGather(arg1, arg2)
    @ccall libpetsc.ISAllGather(arg1::IS, arg2::Ptr{IS})::PetscErrorCode
end

function ISComplement(arg1, arg2, arg3, arg4)
    @ccall libpetsc.ISComplement(arg1::IS, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{IS})::PetscErrorCode
end

function ISConcatenate(arg1, arg2, arg3, arg4)
    @ccall libpetsc.ISConcatenate(arg1::MPI_Comm, arg2::PetscInt, arg3::Ptr{IS}, arg4::Ptr{IS})::PetscErrorCode
end

function ISListToPair(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.ISListToPair(arg1::MPI_Comm, arg2::PetscInt, arg3::Ptr{IS}, arg4::Ptr{IS}, arg5::Ptr{IS})::PetscErrorCode
end

function ISPairToList(arg1, arg2, arg3, arg4)
    @ccall libpetsc.ISPairToList(arg1::IS, arg2::IS, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{IS}})::PetscErrorCode
end

function ISEmbed(arg1, arg2, arg3, arg4)
    @ccall libpetsc.ISEmbed(arg1::IS, arg2::IS, arg3::PetscBool, arg4::Ptr{IS})::PetscErrorCode
end

function ISSortPermutation(arg1, arg2, arg3)
    @ccall libpetsc.ISSortPermutation(arg1::IS, arg2::PetscBool, arg3::Ptr{IS})::PetscErrorCode
end

function ISOnComm(arg1, arg2, arg3, arg4)
    @ccall libpetsc.ISOnComm(arg1::IS, arg2::MPI_Comm, arg3::PetscCopyMode, arg4::Ptr{IS})::PetscErrorCode
end

function ISRenumber(arg1, arg2, arg3, arg4)
    @ccall libpetsc.ISRenumber(arg1::IS, arg2::IS, arg3::Ptr{PetscInt}, arg4::Ptr{IS})::PetscErrorCode
end

function ISCreateSubIS(arg1, arg2, arg3)
    @ccall libpetsc.ISCreateSubIS(arg1::IS, arg2::IS, arg3::Ptr{IS})::PetscErrorCode
end

function ISGeneralFilter(arg1, arg2, arg3)
    @ccall libpetsc.ISGeneralFilter(arg1::IS, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

@enum ISGlobalToLocalMappingMode::UInt32 begin
    IS_GTOLM_MASK = 0
    IS_GTOLM_DROP = 1
end

const ISLocalToGlobalMappingType = Ptr{Cchar}

function ISLocalToGlobalMappingSetType(arg1, arg2)
    @ccall libpetsc.ISLocalToGlobalMappingSetType(arg1::ISLocalToGlobalMapping, arg2::ISLocalToGlobalMappingType)::PetscErrorCode
end

function ISLocalToGlobalMappingRegister(arg1, arg2)
    @ccall libpetsc.ISLocalToGlobalMappingRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function ISLocalToGlobalMappingRegisterAll()
    @ccall libpetsc.ISLocalToGlobalMappingRegisterAll()::PetscErrorCode
end

function ISLocalToGlobalMappingCreate(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.ISLocalToGlobalMappingCreate(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::PetscCopyMode, arg6::Ptr{ISLocalToGlobalMapping})::PetscErrorCode
end

function ISLocalToGlobalMappingCreateIS(arg1, arg2)
    @ccall libpetsc.ISLocalToGlobalMappingCreateIS(arg1::IS, arg2::Ptr{ISLocalToGlobalMapping})::PetscErrorCode
end

function ISLocalToGlobalMappingCreateSF(arg1, arg2, arg3)
    @ccall libpetsc.ISLocalToGlobalMappingCreateSF(arg1::PetscSF, arg2::PetscInt, arg3::Ptr{ISLocalToGlobalMapping})::PetscErrorCode
end

function ISLocalToGlobalMappingSetFromOptions(arg1)
    @ccall libpetsc.ISLocalToGlobalMappingSetFromOptions(arg1::ISLocalToGlobalMapping)::PetscErrorCode
end

function ISLocalToGlobalMappingSetUp(arg1)
    @ccall libpetsc.ISLocalToGlobalMappingSetUp(arg1::ISLocalToGlobalMapping)::PetscErrorCode
end

function ISLocalToGlobalMappingView(arg1, arg2)
    @ccall libpetsc.ISLocalToGlobalMappingView(arg1::ISLocalToGlobalMapping, arg2::PetscViewer)::PetscErrorCode
end

function ISLocalToGlobalMappingViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.ISLocalToGlobalMappingViewFromOptions(arg1::ISLocalToGlobalMapping, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function ISLocalToGlobalMappingDestroy(arg1)
    @ccall libpetsc.ISLocalToGlobalMappingDestroy(arg1::Ptr{ISLocalToGlobalMapping})::PetscErrorCode
end

function ISLocalToGlobalMappingApply(arg1, arg2, arg3, arg4)
    @ccall libpetsc.ISLocalToGlobalMappingApply(arg1::ISLocalToGlobalMapping, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function ISLocalToGlobalMappingApplyBlock(arg1, arg2, arg3, arg4)
    @ccall libpetsc.ISLocalToGlobalMappingApplyBlock(arg1::ISLocalToGlobalMapping, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function ISLocalToGlobalMappingApplyIS(arg1, arg2, arg3)
    @ccall libpetsc.ISLocalToGlobalMappingApplyIS(arg1::ISLocalToGlobalMapping, arg2::IS, arg3::Ptr{IS})::PetscErrorCode
end

function ISGlobalToLocalMappingApply(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.ISGlobalToLocalMappingApply(arg1::ISLocalToGlobalMapping, arg2::ISGlobalToLocalMappingMode, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscInt})::PetscErrorCode
end

function ISGlobalToLocalMappingApplyBlock(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.ISGlobalToLocalMappingApplyBlock(arg1::ISLocalToGlobalMapping, arg2::ISGlobalToLocalMappingMode, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscInt})::PetscErrorCode
end

function ISGlobalToLocalMappingApplyIS(arg1, arg2, arg3, arg4)
    @ccall libpetsc.ISGlobalToLocalMappingApplyIS(arg1::ISLocalToGlobalMapping, arg2::ISGlobalToLocalMappingMode, arg3::IS, arg4::Ptr{IS})::PetscErrorCode
end

function ISLocalToGlobalMappingGetSize(arg1, arg2)
    @ccall libpetsc.ISLocalToGlobalMappingGetSize(arg1::ISLocalToGlobalMapping, arg2::Ptr{PetscInt})::PetscErrorCode
end

function ISLocalToGlobalMappingGetNodeInfo(arg1, arg2, arg3, arg4)
    @ccall libpetsc.ISLocalToGlobalMappingGetNodeInfo(arg1::ISLocalToGlobalMapping, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{PetscInt}}, arg4::Ptr{Ptr{Ptr{PetscInt}}})::PetscErrorCode
end

function ISLocalToGlobalMappingRestoreNodeInfo(arg1, arg2, arg3, arg4)
    @ccall libpetsc.ISLocalToGlobalMappingRestoreNodeInfo(arg1::ISLocalToGlobalMapping, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{PetscInt}}, arg4::Ptr{Ptr{Ptr{PetscInt}}})::PetscErrorCode
end

function ISLocalToGlobalMappingGetInfo(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.ISLocalToGlobalMappingGetInfo(arg1::ISLocalToGlobalMapping, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{PetscInt}}, arg4::Ptr{Ptr{PetscInt}}, arg5::Ptr{Ptr{Ptr{PetscInt}}})::PetscErrorCode
end

function ISLocalToGlobalMappingRestoreInfo(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.ISLocalToGlobalMappingRestoreInfo(arg1::ISLocalToGlobalMapping, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{PetscInt}}, arg4::Ptr{Ptr{PetscInt}}, arg5::Ptr{Ptr{Ptr{PetscInt}}})::PetscErrorCode
end

function ISLocalToGlobalMappingGetBlockInfo(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.ISLocalToGlobalMappingGetBlockInfo(arg1::ISLocalToGlobalMapping, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{PetscInt}}, arg4::Ptr{Ptr{PetscInt}}, arg5::Ptr{Ptr{Ptr{PetscInt}}})::PetscErrorCode
end

function ISLocalToGlobalMappingRestoreBlockInfo(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.ISLocalToGlobalMappingRestoreBlockInfo(arg1::ISLocalToGlobalMapping, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{PetscInt}}, arg4::Ptr{Ptr{PetscInt}}, arg5::Ptr{Ptr{Ptr{PetscInt}}})::PetscErrorCode
end

function ISLocalToGlobalMappingGetIndices(arg1, arg2)
    @ccall libpetsc.ISLocalToGlobalMappingGetIndices(arg1::ISLocalToGlobalMapping, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function ISLocalToGlobalMappingRestoreIndices(arg1, arg2)
    @ccall libpetsc.ISLocalToGlobalMappingRestoreIndices(arg1::ISLocalToGlobalMapping, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function ISLocalToGlobalMappingGetBlockIndices(arg1, arg2)
    @ccall libpetsc.ISLocalToGlobalMappingGetBlockIndices(arg1::ISLocalToGlobalMapping, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function ISLocalToGlobalMappingRestoreBlockIndices(arg1, arg2)
    @ccall libpetsc.ISLocalToGlobalMappingRestoreBlockIndices(arg1::ISLocalToGlobalMapping, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function ISLocalToGlobalMappingConcatenate(arg1, arg2, arg3, arg4)
    @ccall libpetsc.ISLocalToGlobalMappingConcatenate(arg1::MPI_Comm, arg2::PetscInt, arg3::Ptr{ISLocalToGlobalMapping}, arg4::Ptr{ISLocalToGlobalMapping})::PetscErrorCode
end

function ISLocalToGlobalMappingGetBlockSize(arg1, arg2)
    @ccall libpetsc.ISLocalToGlobalMappingGetBlockSize(arg1::ISLocalToGlobalMapping, arg2::Ptr{PetscInt})::PetscErrorCode
end

function ISLocalToGlobalMappingSetBlockSize(arg1, arg2)
    @ccall libpetsc.ISLocalToGlobalMappingSetBlockSize(arg1::ISLocalToGlobalMapping, arg2::PetscInt)::PetscErrorCode
end

function ISLocalToGlobalMappingDuplicate(arg1, arg2)
    @ccall libpetsc.ISLocalToGlobalMappingDuplicate(arg1::ISLocalToGlobalMapping, arg2::Ptr{ISLocalToGlobalMapping})::PetscErrorCode
end

@enum ISColoringType::UInt32 begin
    IS_COLORING_GLOBAL = 0
    IS_COLORING_LOCAL = 1
end

const ISColoringValue = Cushort

function ISAllGatherColors(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.ISAllGatherColors(arg1::MPI_Comm, arg2::PetscInt, arg3::Ptr{ISColoringValue}, arg4::Ptr{PetscInt}, arg5::Ptr{Ptr{ISColoringValue}})::PetscErrorCode
end

function ISColoringCreate(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.ISColoringCreate(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{ISColoringValue}, arg5::PetscCopyMode, arg6::Ptr{ISColoring})::PetscErrorCode
end

function ISColoringDestroy(arg1)
    @ccall libpetsc.ISColoringDestroy(arg1::Ptr{ISColoring})::PetscErrorCode
end

function ISColoringView(arg1, arg2)
    @ccall libpetsc.ISColoringView(arg1::ISColoring, arg2::PetscViewer)::PetscErrorCode
end

function ISColoringViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.ISColoringViewFromOptions(arg1::ISColoring, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function ISColoringGetIS(arg1, arg2, arg3, arg4)
    @ccall libpetsc.ISColoringGetIS(arg1::ISColoring, arg2::PetscCopyMode, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{IS}})::PetscErrorCode
end

function ISColoringRestoreIS(arg1, arg2, arg3)
    @ccall libpetsc.ISColoringRestoreIS(arg1::ISColoring, arg2::PetscCopyMode, arg3::Ptr{Ptr{IS}})::PetscErrorCode
end

function ISColoringReference(arg1)
    @ccall libpetsc.ISColoringReference(arg1::ISColoring)::PetscErrorCode
end

function ISColoringSetType(arg1, arg2)
    @ccall libpetsc.ISColoringSetType(arg1::ISColoring, arg2::ISColoringType)::PetscErrorCode
end

function ISColoringGetType(arg1, arg2)
    @ccall libpetsc.ISColoringGetType(arg1::ISColoring, arg2::Ptr{ISColoringType})::PetscErrorCode
end

function ISColoringGetColors(arg1, arg2, arg3, arg4)
    @ccall libpetsc.ISColoringGetColors(arg1::ISColoring, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{ISColoringValue}})::PetscErrorCode
end

function ISBuildTwoSided(arg1, arg2, arg3)
    @ccall libpetsc.ISBuildTwoSided(arg1::IS, arg2::IS, arg3::Ptr{IS})::PetscErrorCode
end

function ISPartitioningToNumbering(arg1, arg2)
    @ccall libpetsc.ISPartitioningToNumbering(arg1::IS, arg2::Ptr{IS})::PetscErrorCode
end

function ISPartitioningCount(arg1, arg2, arg3)
    @ccall libpetsc.ISPartitioningCount(arg1::IS, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function ISCompressIndicesGeneral(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.ISCompressIndicesGeneral(arg1::PetscInt, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{IS}, arg6::Ptr{IS})::PetscErrorCode
end

function ISCompressIndicesSorted(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.ISCompressIndicesSorted(arg1::PetscInt, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{IS}, arg5::Ptr{IS})::PetscErrorCode
end

function ISExpandIndicesGeneral(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.ISExpandIndicesGeneral(arg1::PetscInt, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{IS}, arg6::Ptr{IS})::PetscErrorCode
end

function PetscLayoutFindOwner(map, idx, owner)
    @ccall libpetsc.PetscLayoutFindOwner(map::PetscLayout, idx::PetscInt, owner::Ptr{PetscMPIInt})::PetscErrorCode
end

function PetscLayoutFindOwnerIndex(map, idx, owner, lidx)
    @ccall libpetsc.PetscLayoutFindOwnerIndex(map::PetscLayout, idx::PetscInt, owner::Ptr{PetscMPIInt}, lidx::Ptr{PetscInt})::PetscErrorCode
end

function PetscLayoutCreate(arg1, arg2)
    @ccall libpetsc.PetscLayoutCreate(arg1::MPI_Comm, arg2::Ptr{PetscLayout})::PetscErrorCode
end

function PetscLayoutCreateFromSizes(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscLayoutCreateFromSizes(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscLayout})::PetscErrorCode
end

function PetscLayoutCreateFromRanges(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscLayoutCreateFromRanges(arg1::MPI_Comm, arg2::Ptr{PetscInt}, arg3::PetscCopyMode, arg4::PetscInt, arg5::Ptr{PetscLayout})::PetscErrorCode
end

function PetscLayoutSetUp(arg1)
    @ccall libpetsc.PetscLayoutSetUp(arg1::PetscLayout)::PetscErrorCode
end

function PetscLayoutDestroy(arg1)
    @ccall libpetsc.PetscLayoutDestroy(arg1::Ptr{PetscLayout})::PetscErrorCode
end

function PetscLayoutDuplicate(arg1, arg2)
    @ccall libpetsc.PetscLayoutDuplicate(arg1::PetscLayout, arg2::Ptr{PetscLayout})::PetscErrorCode
end

function PetscLayoutReference(arg1, arg2)
    @ccall libpetsc.PetscLayoutReference(arg1::PetscLayout, arg2::Ptr{PetscLayout})::PetscErrorCode
end

function PetscLayoutSetLocalSize(arg1, arg2)
    @ccall libpetsc.PetscLayoutSetLocalSize(arg1::PetscLayout, arg2::PetscInt)::PetscErrorCode
end

function PetscLayoutGetLocalSize(arg1, arg2)
    @ccall libpetsc.PetscLayoutGetLocalSize(arg1::PetscLayout, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscLayoutSetSize(arg1, arg2)
    @ccall libpetsc.PetscLayoutSetSize(arg1::PetscLayout, arg2::PetscInt)::PetscErrorCode
end

function PetscLayoutGetSize(arg1, arg2)
    @ccall libpetsc.PetscLayoutGetSize(arg1::PetscLayout, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscLayoutSetBlockSize(arg1, arg2)
    @ccall libpetsc.PetscLayoutSetBlockSize(arg1::PetscLayout, arg2::PetscInt)::PetscErrorCode
end

function PetscLayoutGetBlockSize(arg1, arg2)
    @ccall libpetsc.PetscLayoutGetBlockSize(arg1::PetscLayout, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscLayoutGetRange(arg1, arg2, arg3)
    @ccall libpetsc.PetscLayoutGetRange(arg1::PetscLayout, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscLayoutGetRanges(arg1, arg2)
    @ccall libpetsc.PetscLayoutGetRanges(arg1::PetscLayout, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function PetscLayoutCompare(arg1, arg2, arg3)
    @ccall libpetsc.PetscLayoutCompare(arg1::PetscLayout, arg2::PetscLayout, arg3::Ptr{PetscBool})::PetscErrorCode
end

function PetscLayoutSetISLocalToGlobalMapping(arg1, arg2)
    @ccall libpetsc.PetscLayoutSetISLocalToGlobalMapping(arg1::PetscLayout, arg2::ISLocalToGlobalMapping)::PetscErrorCode
end

function PetscLayoutMapLocal(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscLayoutMapLocal(arg1::PetscLayout, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{Ptr{PetscInt}}, arg6::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function PetscParallelSortInt(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscParallelSortInt(arg1::PetscLayout, arg2::PetscLayout, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function ISGetLayout(arg1, arg2)
    @ccall libpetsc.ISGetLayout(arg1::IS, arg2::Ptr{PetscLayout})::PetscErrorCode
end

mutable struct _p_Vec end

const Vec = Ptr{_p_Vec}

@enum ScatterMode::UInt32 begin
    SCATTER_FORWARD = 0
    SCATTER_REVERSE = 1
    SCATTER_FORWARD_LOCAL = 2
    SCATTER_REVERSE_LOCAL = 3
    # SCATTER_LOCAL = 2
end

const VecType = Ptr{Cchar}

function VecScatterSetType(arg1, arg2)
    @ccall libpetsc.VecScatterSetType(arg1::VecScatter, arg2::VecScatterType)::PetscErrorCode
end

function VecScatterGetType(arg1, arg2)
    @ccall libpetsc.VecScatterGetType(arg1::VecScatter, arg2::Ptr{VecScatterType})::PetscErrorCode
end

function VecScatterSetFromOptions(arg1)
    @ccall libpetsc.VecScatterSetFromOptions(arg1::VecScatter)::PetscErrorCode
end

function VecScatterRegister(arg1, arg2)
    @ccall libpetsc.VecScatterRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function VecScatterCreate(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.VecScatterCreate(arg1::Vec, arg2::IS, arg3::Vec, arg4::IS, arg5::Ptr{VecScatter})::PetscErrorCode
end

function VecInitializePackage()
    @ccall libpetsc.VecInitializePackage()::PetscErrorCode
end

function VecFinalizePackage()
    @ccall libpetsc.VecFinalizePackage()::PetscErrorCode
end

function VecCreate(arg1, arg2)
    @ccall libpetsc.VecCreate(arg1::MPI_Comm, arg2::Ptr{Vec})::PetscErrorCode
end

function VecCreateSeq(arg1, arg2, arg3)
    @ccall libpetsc.VecCreateSeq(arg1::MPI_Comm, arg2::PetscInt, arg3::Ptr{Vec})::PetscErrorCode
end

function VecCreateMPI(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecCreateMPI(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Vec})::PetscErrorCode
end

function VecCreateSeqWithArray(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.VecCreateSeqWithArray(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscScalar}, arg5::Ptr{Vec})::PetscErrorCode
end

function VecCreateMPIWithArray(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.VecCreateMPIWithArray(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscScalar}, arg6::Ptr{Vec})::PetscErrorCode
end

function VecCreateShared(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecCreateShared(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Vec})::PetscErrorCode
end

function VecCreateNode(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecCreateNode(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Vec})::PetscErrorCode
end

function VecSetFromOptions(arg1)
    @ccall libpetsc.VecSetFromOptions(arg1::Vec)::PetscErrorCode
end

function VecViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.VecViewFromOptions(arg1::Vec, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function VecSetUp(arg1)
    @ccall libpetsc.VecSetUp(arg1::Vec)::PetscErrorCode
end

function VecDestroy(arg1)
    @ccall libpetsc.VecDestroy(arg1::Ptr{Vec})::PetscErrorCode
end

function VecZeroEntries(arg1)
    @ccall libpetsc.VecZeroEntries(arg1::Vec)::PetscErrorCode
end

function VecSetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.VecSetOptionsPrefix(arg1::Vec, arg2::Ptr{Cchar})::PetscErrorCode
end

function VecAppendOptionsPrefix(arg1, arg2)
    @ccall libpetsc.VecAppendOptionsPrefix(arg1::Vec, arg2::Ptr{Cchar})::PetscErrorCode
end

function VecGetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.VecGetOptionsPrefix(arg1::Vec, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function VecSetSizes(arg1, arg2, arg3)
    @ccall libpetsc.VecSetSizes(arg1::Vec, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function VecDotNorm2(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecDotNorm2(arg1::Vec, arg2::Vec, arg3::Ptr{PetscScalar}, arg4::Ptr{PetscReal})::PetscErrorCode
end

function VecDot(arg1, arg2, arg3)
    @ccall libpetsc.VecDot(arg1::Vec, arg2::Vec, arg3::Ptr{PetscScalar})::PetscErrorCode
end

function VecDotRealPart(arg1, arg2, arg3)
    @ccall libpetsc.VecDotRealPart(arg1::Vec, arg2::Vec, arg3::Ptr{PetscReal})::PetscErrorCode
end

function VecTDot(arg1, arg2, arg3)
    @ccall libpetsc.VecTDot(arg1::Vec, arg2::Vec, arg3::Ptr{PetscScalar})::PetscErrorCode
end

function VecMDot(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecMDot(arg1::Vec, arg2::PetscInt, arg3::Ptr{Vec}, arg4::Ptr{PetscScalar})::PetscErrorCode
end

function VecMTDot(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecMTDot(arg1::Vec, arg2::PetscInt, arg3::Ptr{Vec}, arg4::Ptr{PetscScalar})::PetscErrorCode
end

function VecGetSubVector(arg1, arg2, arg3)
    @ccall libpetsc.VecGetSubVector(arg1::Vec, arg2::IS, arg3::Ptr{Vec})::PetscErrorCode
end

function VecRestoreSubVector(arg1, arg2, arg3)
    @ccall libpetsc.VecRestoreSubVector(arg1::Vec, arg2::IS, arg3::Ptr{Vec})::PetscErrorCode
end

function VecConcatenate(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecConcatenate(arg1::PetscInt, arg2::Ptr{Vec}, arg3::Ptr{Vec}, arg4::Ptr{Ptr{IS}})::PetscErrorCode
end

@enum NormType::UInt32 begin
    NORM_1 = 0
    NORM_2 = 1
    NORM_FROBENIUS = 2
    NORM_INFINITY = 3
    NORM_1_AND_2 = 4
end

function VecNorm(arg1, arg2, arg3)
    @ccall libpetsc.VecNorm(arg1::Vec, arg2::NormType, arg3::Ptr{PetscReal})::PetscErrorCode
end

function VecNormAvailable(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecNormAvailable(arg1::Vec, arg2::NormType, arg3::Ptr{PetscBool}, arg4::Ptr{PetscReal})::PetscErrorCode
end

function VecNormalize(arg1, arg2)
    @ccall libpetsc.VecNormalize(arg1::Vec, arg2::Ptr{PetscReal})::PetscErrorCode
end

function VecSum(arg1, arg2)
    @ccall libpetsc.VecSum(arg1::Vec, arg2::Ptr{PetscScalar})::PetscErrorCode
end

function VecMax(arg1, arg2, arg3)
    @ccall libpetsc.VecMax(arg1::Vec, arg2::Ptr{PetscInt}, arg3::Ptr{PetscReal})::PetscErrorCode
end

function VecMin(arg1, arg2, arg3)
    @ccall libpetsc.VecMin(arg1::Vec, arg2::Ptr{PetscInt}, arg3::Ptr{PetscReal})::PetscErrorCode
end

function VecScale(arg1, arg2)
    @ccall libpetsc.VecScale(arg1::Vec, arg2::PetscScalar)::PetscErrorCode
end

function VecCopy(arg1, arg2)
    @ccall libpetsc.VecCopy(arg1::Vec, arg2::Vec)::PetscErrorCode
end

function VecSetRandom(arg1, arg2)
    @ccall libpetsc.VecSetRandom(arg1::Vec, arg2::PetscRandom)::PetscErrorCode
end

function VecSet(arg1, arg2)
    @ccall libpetsc.VecSet(arg1::Vec, arg2::PetscScalar)::PetscErrorCode
end

function VecSetInf(arg1)
    @ccall libpetsc.VecSetInf(arg1::Vec)::PetscErrorCode
end

function VecSwap(arg1, arg2)
    @ccall libpetsc.VecSwap(arg1::Vec, arg2::Vec)::PetscErrorCode
end

function VecAXPY(arg1, arg2, arg3)
    @ccall libpetsc.VecAXPY(arg1::Vec, arg2::PetscScalar, arg3::Vec)::PetscErrorCode
end

function VecAXPBY(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecAXPBY(arg1::Vec, arg2::PetscScalar, arg3::PetscScalar, arg4::Vec)::PetscErrorCode
end

function VecMAXPY(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecMAXPY(arg1::Vec, arg2::PetscInt, arg3::Ptr{PetscScalar}, arg4::Ptr{Vec})::PetscErrorCode
end

function VecAYPX(arg1, arg2, arg3)
    @ccall libpetsc.VecAYPX(arg1::Vec, arg2::PetscScalar, arg3::Vec)::PetscErrorCode
end

function VecWAXPY(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecWAXPY(arg1::Vec, arg2::PetscScalar, arg3::Vec, arg4::Vec)::PetscErrorCode
end

function VecAXPBYPCZ(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.VecAXPBYPCZ(arg1::Vec, arg2::PetscScalar, arg3::PetscScalar, arg4::PetscScalar, arg5::Vec, arg6::Vec)::PetscErrorCode
end

function VecPointwiseMax(arg1, arg2, arg3)
    @ccall libpetsc.VecPointwiseMax(arg1::Vec, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function VecPointwiseMaxAbs(arg1, arg2, arg3)
    @ccall libpetsc.VecPointwiseMaxAbs(arg1::Vec, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function VecPointwiseMin(arg1, arg2, arg3)
    @ccall libpetsc.VecPointwiseMin(arg1::Vec, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function VecPointwiseMult(arg1, arg2, arg3)
    @ccall libpetsc.VecPointwiseMult(arg1::Vec, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function VecPointwiseDivide(arg1, arg2, arg3)
    @ccall libpetsc.VecPointwiseDivide(arg1::Vec, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function VecMaxPointwiseDivide(arg1, arg2, arg3)
    @ccall libpetsc.VecMaxPointwiseDivide(arg1::Vec, arg2::Vec, arg3::Ptr{PetscReal})::PetscErrorCode
end

function VecShift(arg1, arg2)
    @ccall libpetsc.VecShift(arg1::Vec, arg2::PetscScalar)::PetscErrorCode
end

function VecReciprocal(arg1)
    @ccall libpetsc.VecReciprocal(arg1::Vec)::PetscErrorCode
end

function VecPermute(arg1, arg2, arg3)
    @ccall libpetsc.VecPermute(arg1::Vec, arg2::IS, arg3::PetscBool)::PetscErrorCode
end

function VecSqrtAbs(arg1)
    @ccall libpetsc.VecSqrtAbs(arg1::Vec)::PetscErrorCode
end

function VecLog(arg1)
    @ccall libpetsc.VecLog(arg1::Vec)::PetscErrorCode
end

function VecExp(arg1)
    @ccall libpetsc.VecExp(arg1::Vec)::PetscErrorCode
end

function VecAbs(arg1)
    @ccall libpetsc.VecAbs(arg1::Vec)::PetscErrorCode
end

function VecDuplicate(arg1, arg2)
    @ccall libpetsc.VecDuplicate(arg1::Vec, arg2::Ptr{Vec})::PetscErrorCode
end

function VecDuplicateVecs(arg1, arg2, arg3)
    @ccall libpetsc.VecDuplicateVecs(arg1::Vec, arg2::PetscInt, arg3::Ptr{Ptr{Vec}})::PetscErrorCode
end

function VecDestroyVecs(arg1, arg2)
    @ccall libpetsc.VecDestroyVecs(arg1::PetscInt, arg2::Ptr{Ptr{Vec}})::PetscErrorCode
end

function VecStrideNormAll(arg1, arg2, arg3)
    @ccall libpetsc.VecStrideNormAll(arg1::Vec, arg2::NormType, arg3::Ptr{PetscReal})::PetscErrorCode
end

function VecStrideMaxAll(arg1, arg2, arg3)
    @ccall libpetsc.VecStrideMaxAll(arg1::Vec, arg2::Ptr{PetscInt}, arg3::Ptr{PetscReal})::PetscErrorCode
end

function VecStrideMinAll(arg1, arg2, arg3)
    @ccall libpetsc.VecStrideMinAll(arg1::Vec, arg2::Ptr{PetscInt}, arg3::Ptr{PetscReal})::PetscErrorCode
end

function VecStrideScaleAll(arg1, arg2)
    @ccall libpetsc.VecStrideScaleAll(arg1::Vec, arg2::Ptr{PetscScalar})::PetscErrorCode
end

function VecUniqueEntries(arg1, arg2, arg3)
    @ccall libpetsc.VecUniqueEntries(arg1::Vec, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecStrideNorm(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecStrideNorm(arg1::Vec, arg2::PetscInt, arg3::NormType, arg4::Ptr{PetscReal})::PetscErrorCode
end

function VecStrideMax(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecStrideMax(arg1::Vec, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscReal})::PetscErrorCode
end

function VecStrideMin(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecStrideMin(arg1::Vec, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscReal})::PetscErrorCode
end

function VecStrideScale(arg1, arg2, arg3)
    @ccall libpetsc.VecStrideScale(arg1::Vec, arg2::PetscInt, arg3::PetscScalar)::PetscErrorCode
end

function VecStrideSet(arg1, arg2, arg3)
    @ccall libpetsc.VecStrideSet(arg1::Vec, arg2::PetscInt, arg3::PetscScalar)::PetscErrorCode
end

function VecStrideGather(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecStrideGather(arg1::Vec, arg2::PetscInt, arg3::Vec, arg4::InsertMode)::PetscErrorCode
end

function VecStrideScatter(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecStrideScatter(arg1::Vec, arg2::PetscInt, arg3::Vec, arg4::InsertMode)::PetscErrorCode
end

function VecStrideGatherAll(arg1, arg2, arg3)
    @ccall libpetsc.VecStrideGatherAll(arg1::Vec, arg2::Ptr{Vec}, arg3::InsertMode)::PetscErrorCode
end

function VecStrideScatterAll(arg1, arg2, arg3)
    @ccall libpetsc.VecStrideScatterAll(arg1::Ptr{Vec}, arg2::Vec, arg3::InsertMode)::PetscErrorCode
end

function VecStrideSubSetScatter(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.VecStrideSubSetScatter(arg1::Vec, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Vec, arg6::InsertMode)::PetscErrorCode
end

function VecStrideSubSetGather(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.VecStrideSubSetGather(arg1::Vec, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Vec, arg6::InsertMode)::PetscErrorCode
end

function VecSetValues(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.VecSetValues(arg1::Vec, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscScalar}, arg5::InsertMode)::PetscErrorCode
end

function VecGetValues(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecGetValues(arg1::Vec, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscScalar})::PetscErrorCode
end

function VecAssemblyBegin(arg1)
    @ccall libpetsc.VecAssemblyBegin(arg1::Vec)::PetscErrorCode
end

function VecAssemblyEnd(arg1)
    @ccall libpetsc.VecAssemblyEnd(arg1::Vec)::PetscErrorCode
end

function VecStashSetInitialSize(arg1, arg2, arg3)
    @ccall libpetsc.VecStashSetInitialSize(arg1::Vec, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function VecStashView(arg1, arg2)
    @ccall libpetsc.VecStashView(arg1::Vec, arg2::PetscViewer)::PetscErrorCode
end

function VecStashViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.VecStashViewFromOptions(arg1::Vec, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function VecStashGetInfo(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.VecStashGetInfo(arg1::Vec, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt})::PetscErrorCode
end

function VecSetValue(v, i, va, mode)
    @ccall libpetsc.VecSetValue(v::Vec, i::PetscInt, va::PetscScalar, mode::InsertMode)::PetscErrorCode
end

function VecSetBlockSize(arg1, arg2)
    @ccall libpetsc.VecSetBlockSize(arg1::Vec, arg2::PetscInt)::PetscErrorCode
end

function VecGetBlockSize(arg1, arg2)
    @ccall libpetsc.VecGetBlockSize(arg1::Vec, arg2::Ptr{PetscInt})::PetscErrorCode
end

function VecSetValuesBlocked(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.VecSetValuesBlocked(arg1::Vec, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscScalar}, arg5::InsertMode)::PetscErrorCode
end

function VecSetType(arg1, arg2)
    @ccall libpetsc.VecSetType(arg1::Vec, arg2::VecType)::PetscErrorCode
end

function VecGetType(arg1, arg2)
    @ccall libpetsc.VecGetType(arg1::Vec, arg2::Ptr{VecType})::PetscErrorCode
end

function VecRegister(arg1, arg2)
    @ccall libpetsc.VecRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function VecScatterBegin(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.VecScatterBegin(arg1::VecScatter, arg2::Vec, arg3::Vec, arg4::InsertMode, arg5::ScatterMode)::PetscErrorCode
end

function VecScatterEnd(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.VecScatterEnd(arg1::VecScatter, arg2::Vec, arg3::Vec, arg4::InsertMode, arg5::ScatterMode)::PetscErrorCode
end

function VecScatterDestroy(arg1)
    @ccall libpetsc.VecScatterDestroy(arg1::Ptr{VecScatter})::PetscErrorCode
end

function VecScatterSetUp(arg1)
    @ccall libpetsc.VecScatterSetUp(arg1::VecScatter)::PetscErrorCode
end

function VecScatterCopy(arg1, arg2)
    @ccall libpetsc.VecScatterCopy(arg1::VecScatter, arg2::Ptr{VecScatter})::PetscErrorCode
end

function VecScatterView(arg1, arg2)
    @ccall libpetsc.VecScatterView(arg1::VecScatter, arg2::PetscViewer)::PetscErrorCode
end

function VecScatterViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.VecScatterViewFromOptions(arg1::VecScatter, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function VecScatterRemap(arg1, arg2, arg3)
    @ccall libpetsc.VecScatterRemap(arg1::VecScatter, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function VecScatterGetMerged(arg1, arg2)
    @ccall libpetsc.VecScatterGetMerged(arg1::VecScatter, arg2::Ptr{PetscBool})::PetscErrorCode
end

function VecGetArray4d(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.VecGetArray4d(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::PetscInt, arg8::PetscInt, arg9::PetscInt, arg10::Ptr{Ptr{Ptr{Ptr{Ptr{PetscScalar}}}}})::PetscErrorCode
end

function VecRestoreArray4d(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.VecRestoreArray4d(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::PetscInt, arg8::PetscInt, arg9::PetscInt, arg10::Ptr{Ptr{Ptr{Ptr{Ptr{PetscScalar}}}}})::PetscErrorCode
end

function VecGetArray3d(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.VecGetArray3d(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::PetscInt, arg8::Ptr{Ptr{Ptr{Ptr{PetscScalar}}}})::PetscErrorCode
end

function VecRestoreArray3d(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.VecRestoreArray3d(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::PetscInt, arg8::Ptr{Ptr{Ptr{Ptr{PetscScalar}}}})::PetscErrorCode
end

function VecGetArray2d(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.VecGetArray2d(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{Ptr{Ptr{PetscScalar}}})::PetscErrorCode
end

function VecRestoreArray2d(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.VecRestoreArray2d(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{Ptr{Ptr{PetscScalar}}})::PetscErrorCode
end

function VecGetArray1d(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecGetArray1d(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecRestoreArray1d(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecRestoreArray1d(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecGetArray4dWrite(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.VecGetArray4dWrite(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::PetscInt, arg8::PetscInt, arg9::PetscInt, arg10::Ptr{Ptr{Ptr{Ptr{Ptr{PetscScalar}}}}})::PetscErrorCode
end

function VecRestoreArray4dWrite(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.VecRestoreArray4dWrite(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::PetscInt, arg8::PetscInt, arg9::PetscInt, arg10::Ptr{Ptr{Ptr{Ptr{Ptr{PetscScalar}}}}})::PetscErrorCode
end

function VecGetArray3dWrite(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.VecGetArray3dWrite(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::PetscInt, arg8::Ptr{Ptr{Ptr{Ptr{PetscScalar}}}})::PetscErrorCode
end

function VecRestoreArray3dWrite(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.VecRestoreArray3dWrite(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::PetscInt, arg8::Ptr{Ptr{Ptr{Ptr{PetscScalar}}}})::PetscErrorCode
end

function VecGetArray2dWrite(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.VecGetArray2dWrite(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{Ptr{Ptr{PetscScalar}}})::PetscErrorCode
end

function VecRestoreArray2dWrite(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.VecRestoreArray2dWrite(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{Ptr{Ptr{PetscScalar}}})::PetscErrorCode
end

function VecGetArray1dWrite(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecGetArray1dWrite(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecRestoreArray1dWrite(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecRestoreArray1dWrite(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecGetArray4dRead(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.VecGetArray4dRead(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::PetscInt, arg8::PetscInt, arg9::PetscInt, arg10::Ptr{Ptr{Ptr{Ptr{Ptr{PetscScalar}}}}})::PetscErrorCode
end

function VecRestoreArray4dRead(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.VecRestoreArray4dRead(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::PetscInt, arg8::PetscInt, arg9::PetscInt, arg10::Ptr{Ptr{Ptr{Ptr{Ptr{PetscScalar}}}}})::PetscErrorCode
end

function VecGetArray3dRead(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.VecGetArray3dRead(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::PetscInt, arg8::Ptr{Ptr{Ptr{Ptr{PetscScalar}}}})::PetscErrorCode
end

function VecRestoreArray3dRead(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.VecRestoreArray3dRead(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::PetscInt, arg8::Ptr{Ptr{Ptr{Ptr{PetscScalar}}}})::PetscErrorCode
end

function VecGetArray2dRead(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.VecGetArray2dRead(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{Ptr{Ptr{PetscScalar}}})::PetscErrorCode
end

function VecRestoreArray2dRead(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.VecRestoreArray2dRead(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{Ptr{Ptr{PetscScalar}}})::PetscErrorCode
end

function VecGetArray1dRead(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecGetArray1dRead(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecRestoreArray1dRead(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecRestoreArray1dRead(arg1::Vec, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecPlaceArray(arg1, arg2)
    @ccall libpetsc.VecPlaceArray(arg1::Vec, arg2::Ptr{PetscScalar})::PetscErrorCode
end

function VecResetArray(arg1)
    @ccall libpetsc.VecResetArray(arg1::Vec)::PetscErrorCode
end

function VecReplaceArray(arg1, arg2)
    @ccall libpetsc.VecReplaceArray(arg1::Vec, arg2::Ptr{PetscScalar})::PetscErrorCode
end

function VecGetArrays(arg1, arg2, arg3)
    @ccall libpetsc.VecGetArrays(arg1::Ptr{Vec}, arg2::PetscInt, arg3::Ptr{Ptr{Ptr{PetscScalar}}})::PetscErrorCode
end

function VecRestoreArrays(arg1, arg2, arg3)
    @ccall libpetsc.VecRestoreArrays(arg1::Ptr{Vec}, arg2::PetscInt, arg3::Ptr{Ptr{Ptr{PetscScalar}}})::PetscErrorCode
end

function VecView(arg1, arg2)
    @ccall libpetsc.VecView(arg1::Vec, arg2::PetscViewer)::PetscErrorCode
end

function VecEqual(arg1, arg2, arg3)
    @ccall libpetsc.VecEqual(arg1::Vec, arg2::Vec, arg3::Ptr{PetscBool})::PetscErrorCode
end

function VecLoad(arg1, arg2)
    @ccall libpetsc.VecLoad(arg1::Vec, arg2::PetscViewer)::PetscErrorCode
end

function VecGetSize(arg1, arg2)
    @ccall libpetsc.VecGetSize(arg1::Vec, arg2::Ptr{PetscInt})::PetscErrorCode
end

function VecGetLocalSize(arg1, arg2)
    @ccall libpetsc.VecGetLocalSize(arg1::Vec, arg2::Ptr{PetscInt})::PetscErrorCode
end

function VecGetOwnershipRange(arg1, arg2, arg3)
    @ccall libpetsc.VecGetOwnershipRange(arg1::Vec, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function VecGetOwnershipRanges(arg1, arg2)
    @ccall libpetsc.VecGetOwnershipRanges(arg1::Vec, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function VecSetLocalToGlobalMapping(arg1, arg2)
    @ccall libpetsc.VecSetLocalToGlobalMapping(arg1::Vec, arg2::ISLocalToGlobalMapping)::PetscErrorCode
end

function VecSetValuesLocal(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.VecSetValuesLocal(arg1::Vec, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscScalar}, arg5::InsertMode)::PetscErrorCode
end

function VecCUDAGetArray(arg1, arg2)
    @ccall libpetsc.VecCUDAGetArray(arg1::Vec, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecCUDARestoreArray(arg1, arg2)
    @ccall libpetsc.VecCUDARestoreArray(arg1::Vec, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecCUDAGetArrayRead(arg1, arg2)
    @ccall libpetsc.VecCUDAGetArrayRead(arg1::Vec, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecCUDARestoreArrayRead(arg1, arg2)
    @ccall libpetsc.VecCUDARestoreArrayRead(arg1::Vec, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecCUDAGetArrayWrite(arg1, arg2)
    @ccall libpetsc.VecCUDAGetArrayWrite(arg1::Vec, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecCUDARestoreArrayWrite(arg1, arg2)
    @ccall libpetsc.VecCUDARestoreArrayWrite(arg1::Vec, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecCUDAPlaceArray(arg1, arg2)
    @ccall libpetsc.VecCUDAPlaceArray(arg1::Vec, arg2::Ptr{PetscScalar})::PetscErrorCode
end

function VecCUDAReplaceArray(arg1, arg2)
    @ccall libpetsc.VecCUDAReplaceArray(arg1::Vec, arg2::Ptr{PetscScalar})::PetscErrorCode
end

function VecCUDAResetArray(arg1)
    @ccall libpetsc.VecCUDAResetArray(arg1::Vec)::PetscErrorCode
end

function VecHIPGetArray(arg1, arg2)
    @ccall libpetsc.VecHIPGetArray(arg1::Vec, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecHIPRestoreArray(arg1, arg2)
    @ccall libpetsc.VecHIPRestoreArray(arg1::Vec, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecHIPGetArrayRead(arg1, arg2)
    @ccall libpetsc.VecHIPGetArrayRead(arg1::Vec, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecHIPRestoreArrayRead(arg1, arg2)
    @ccall libpetsc.VecHIPRestoreArrayRead(arg1::Vec, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecHIPGetArrayWrite(arg1, arg2)
    @ccall libpetsc.VecHIPGetArrayWrite(arg1::Vec, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecHIPRestoreArrayWrite(arg1, arg2)
    @ccall libpetsc.VecHIPRestoreArrayWrite(arg1::Vec, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecHIPPlaceArray(arg1, arg2)
    @ccall libpetsc.VecHIPPlaceArray(arg1::Vec, arg2::Ptr{PetscScalar})::PetscErrorCode
end

function VecHIPReplaceArray(arg1, arg2)
    @ccall libpetsc.VecHIPReplaceArray(arg1::Vec, arg2::Ptr{PetscScalar})::PetscErrorCode
end

function VecHIPResetArray(arg1)
    @ccall libpetsc.VecHIPResetArray(arg1::Vec)::PetscErrorCode
end

function VecViennaCLGetCLContext(arg1, arg2)
    @ccall libpetsc.VecViennaCLGetCLContext(arg1::Vec, arg2::Ptr{Csize_t})::PetscErrorCode
end

function VecViennaCLGetCLQueue(arg1, arg2)
    @ccall libpetsc.VecViennaCLGetCLQueue(arg1::Vec, arg2::Ptr{Csize_t})::PetscErrorCode
end

function VecViennaCLGetCLMemRead(arg1, arg2)
    @ccall libpetsc.VecViennaCLGetCLMemRead(arg1::Vec, arg2::Ptr{Csize_t})::PetscErrorCode
end

function VecViennaCLGetCLMemWrite(arg1, arg2)
    @ccall libpetsc.VecViennaCLGetCLMemWrite(arg1::Vec, arg2::Ptr{Csize_t})::PetscErrorCode
end

function VecViennaCLRestoreCLMemWrite(arg1)
    @ccall libpetsc.VecViennaCLRestoreCLMemWrite(arg1::Vec)::PetscErrorCode
end

function VecViennaCLGetCLMem(arg1, arg2)
    @ccall libpetsc.VecViennaCLGetCLMem(arg1::Vec, arg2::Ptr{Csize_t})::PetscErrorCode
end

function VecViennaCLRestoreCLMem(arg1)
    @ccall libpetsc.VecViennaCLRestoreCLMem(arg1::Vec)::PetscErrorCode
end

function VecSetValueLocal(v, i, va, mode)
    @ccall libpetsc.VecSetValueLocal(v::Vec, i::PetscInt, va::PetscScalar, mode::InsertMode)::PetscErrorCode
end

function VecSetValuesBlockedLocal(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.VecSetValuesBlockedLocal(arg1::Vec, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscScalar}, arg5::InsertMode)::PetscErrorCode
end

function VecGetLocalToGlobalMapping(arg1, arg2)
    @ccall libpetsc.VecGetLocalToGlobalMapping(arg1::Vec, arg2::Ptr{ISLocalToGlobalMapping})::PetscErrorCode
end

function VecDotBegin(arg1, arg2, arg3)
    @ccall libpetsc.VecDotBegin(arg1::Vec, arg2::Vec, arg3::Ptr{PetscScalar})::PetscErrorCode
end

function VecDotEnd(arg1, arg2, arg3)
    @ccall libpetsc.VecDotEnd(arg1::Vec, arg2::Vec, arg3::Ptr{PetscScalar})::PetscErrorCode
end

function VecTDotBegin(arg1, arg2, arg3)
    @ccall libpetsc.VecTDotBegin(arg1::Vec, arg2::Vec, arg3::Ptr{PetscScalar})::PetscErrorCode
end

function VecTDotEnd(arg1, arg2, arg3)
    @ccall libpetsc.VecTDotEnd(arg1::Vec, arg2::Vec, arg3::Ptr{PetscScalar})::PetscErrorCode
end

function VecNormBegin(arg1, arg2, arg3)
    @ccall libpetsc.VecNormBegin(arg1::Vec, arg2::NormType, arg3::Ptr{PetscReal})::PetscErrorCode
end

function VecNormEnd(arg1, arg2, arg3)
    @ccall libpetsc.VecNormEnd(arg1::Vec, arg2::NormType, arg3::Ptr{PetscReal})::PetscErrorCode
end

function VecMDotBegin(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecMDotBegin(arg1::Vec, arg2::PetscInt, arg3::Ptr{Vec}, arg4::Ptr{PetscScalar})::PetscErrorCode
end

function VecMDotEnd(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecMDotEnd(arg1::Vec, arg2::PetscInt, arg3::Ptr{Vec}, arg4::Ptr{PetscScalar})::PetscErrorCode
end

function VecMTDotBegin(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecMTDotBegin(arg1::Vec, arg2::PetscInt, arg3::Ptr{Vec}, arg4::Ptr{PetscScalar})::PetscErrorCode
end

function VecMTDotEnd(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecMTDotEnd(arg1::Vec, arg2::PetscInt, arg3::Ptr{Vec}, arg4::Ptr{PetscScalar})::PetscErrorCode
end

function PetscCommSplitReductionBegin(arg1)
    @ccall libpetsc.PetscCommSplitReductionBegin(arg1::MPI_Comm)::PetscErrorCode
end

function VecBindToCPU(arg1, arg2)
    @ccall libpetsc.VecBindToCPU(arg1::Vec, arg2::PetscBool)::PetscErrorCode
end

function VecPinToCPU(v, flg)
    @ccall libpetsc.VecPinToCPU(v::Vec, flg::PetscBool)::PetscErrorCode
end

function VecSetPinnedMemoryMin(arg1, arg2)
    @ccall libpetsc.VecSetPinnedMemoryMin(arg1::Vec, arg2::Csize_t)::PetscErrorCode
end

function VecGetPinnedMemoryMin(arg1, arg2)
    @ccall libpetsc.VecGetPinnedMemoryMin(arg1::Vec, arg2::Ptr{Csize_t})::PetscErrorCode
end

@enum PetscOffloadMask::UInt32 begin
    PETSC_OFFLOAD_UNALLOCATED = 0
    PETSC_OFFLOAD_CPU = 1
    PETSC_OFFLOAD_GPU = 2
    PETSC_OFFLOAD_BOTH = 3
    PETSC_OFFLOAD_VECKOKKOS = 256
end

function VecGetOffloadMask(arg1, arg2)
    @ccall libpetsc.VecGetOffloadMask(arg1::Vec, arg2::Ptr{PetscOffloadMask})::PetscErrorCode
end

@enum VecOption::UInt32 begin
    VEC_IGNORE_OFF_PROC_ENTRIES = 0
    VEC_IGNORE_NEGATIVE_INDICES = 1
    VEC_SUBSET_OFF_PROC_ENTRIES = 2
end

function VecSetOption(arg1, arg2, arg3)
    @ccall libpetsc.VecSetOption(arg1::Vec, arg2::VecOption, arg3::PetscBool)::PetscErrorCode
end

function VecGetArray(arg1, arg2)
    @ccall libpetsc.VecGetArray(arg1::Vec, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecGetArrayWrite(arg1, arg2)
    @ccall libpetsc.VecGetArrayWrite(arg1::Vec, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecGetArrayRead(arg1, arg2)
    @ccall libpetsc.VecGetArrayRead(arg1::Vec, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecRestoreArray(arg1, arg2)
    @ccall libpetsc.VecRestoreArray(arg1::Vec, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecRestoreArrayWrite(arg1, arg2)
    @ccall libpetsc.VecRestoreArrayWrite(arg1::Vec, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecRestoreArrayRead(arg1, arg2)
    @ccall libpetsc.VecRestoreArrayRead(arg1::Vec, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecGetLocalVector(arg1, arg2)
    @ccall libpetsc.VecGetLocalVector(arg1::Vec, arg2::Vec)::PetscErrorCode
end

function VecRestoreLocalVector(arg1, arg2)
    @ccall libpetsc.VecRestoreLocalVector(arg1::Vec, arg2::Vec)::PetscErrorCode
end

function VecGetLocalVectorRead(arg1, arg2)
    @ccall libpetsc.VecGetLocalVectorRead(arg1::Vec, arg2::Vec)::PetscErrorCode
end

function VecRestoreLocalVectorRead(arg1, arg2)
    @ccall libpetsc.VecRestoreLocalVectorRead(arg1::Vec, arg2::Vec)::PetscErrorCode
end

function VecGetArrayAndMemType(arg1, arg2, arg3)
    @ccall libpetsc.VecGetArrayAndMemType(arg1::Vec, arg2::Ptr{Ptr{PetscScalar}}, arg3::Ptr{PetscMemType})::PetscErrorCode
end

function VecRestoreArrayAndMemType(arg1, arg2)
    @ccall libpetsc.VecRestoreArrayAndMemType(arg1::Vec, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecGetArrayReadAndMemType(arg1, arg2, arg3)
    @ccall libpetsc.VecGetArrayReadAndMemType(arg1::Vec, arg2::Ptr{Ptr{PetscScalar}}, arg3::Ptr{PetscMemType})::PetscErrorCode
end

function VecRestoreArrayReadAndMemType(arg1, arg2)
    @ccall libpetsc.VecRestoreArrayReadAndMemType(arg1::Vec, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecGetArrayPair(x, y, xv, yv)
    @ccall libpetsc.VecGetArrayPair(x::Vec, y::Vec, xv::Ptr{Ptr{PetscScalar}}, yv::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecRestoreArrayPair(x, y, xv, yv)
    @ccall libpetsc.VecRestoreArrayPair(x::Vec, y::Vec, xv::Ptr{Ptr{PetscScalar}}, yv::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecValidValues(arg1, arg2, arg3)
    @ccall libpetsc.VecValidValues(arg1::Vec, arg2::PetscInt, arg3::PetscBool)::PetscErrorCode
end

@enum VecOperation::UInt32 begin
    VECOP_DUPLICATE = 0
    VECOP_VIEW = 33
    VECOP_LOAD = 41
    VECOP_VIEWNATIVE = 68
    VECOP_LOADNATIVE = 69
end

function VecSetOperation(arg1, arg2, arg3)
    @ccall libpetsc.VecSetOperation(arg1::Vec, arg2::VecOperation, arg3::Ptr{Cvoid})::PetscErrorCode
end

function VecMPISetGhost(arg1, arg2, arg3)
    @ccall libpetsc.VecMPISetGhost(arg1::Vec, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function VecCreateGhost(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.VecCreateGhost(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{Vec})::PetscErrorCode
end

function VecCreateGhostWithArray(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.VecCreateGhostWithArray(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{PetscScalar}, arg7::Ptr{Vec})::PetscErrorCode
end

function VecCreateGhostBlock(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.VecCreateGhostBlock(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::Ptr{Vec})::PetscErrorCode
end

function VecCreateGhostBlockWithArray(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.VecCreateGhostBlockWithArray(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::Ptr{PetscScalar}, arg8::Ptr{Vec})::PetscErrorCode
end

function VecGhostGetLocalForm(arg1, arg2)
    @ccall libpetsc.VecGhostGetLocalForm(arg1::Vec, arg2::Ptr{Vec})::PetscErrorCode
end

function VecGhostRestoreLocalForm(arg1, arg2)
    @ccall libpetsc.VecGhostRestoreLocalForm(arg1::Vec, arg2::Ptr{Vec})::PetscErrorCode
end

function VecGhostIsLocalForm(arg1, arg2, arg3)
    @ccall libpetsc.VecGhostIsLocalForm(arg1::Vec, arg2::Vec, arg3::Ptr{PetscBool})::PetscErrorCode
end

function VecGhostUpdateBegin(arg1, arg2, arg3)
    @ccall libpetsc.VecGhostUpdateBegin(arg1::Vec, arg2::InsertMode, arg3::ScatterMode)::PetscErrorCode
end

function VecGhostUpdateEnd(arg1, arg2, arg3)
    @ccall libpetsc.VecGhostUpdateEnd(arg1::Vec, arg2::InsertMode, arg3::ScatterMode)::PetscErrorCode
end

function VecConjugate(arg1)
    @ccall libpetsc.VecConjugate(arg1::Vec)::PetscErrorCode
end

function VecImaginaryPart(arg1)
    @ccall libpetsc.VecImaginaryPart(arg1::Vec)::PetscErrorCode
end

function VecRealPart(arg1)
    @ccall libpetsc.VecRealPart(arg1::Vec)::PetscErrorCode
end

function VecScatterCreateToAll(arg1, arg2, arg3)
    @ccall libpetsc.VecScatterCreateToAll(arg1::Vec, arg2::Ptr{VecScatter}, arg3::Ptr{Vec})::PetscErrorCode
end

function VecScatterCreateToZero(arg1, arg2, arg3)
    @ccall libpetsc.VecScatterCreateToZero(arg1::Vec, arg2::Ptr{VecScatter}, arg3::Ptr{Vec})::PetscErrorCode
end

function ISComplementVec(arg1, arg2, arg3)
    @ccall libpetsc.ISComplementVec(arg1::IS, arg2::Vec, arg3::Ptr{IS})::PetscErrorCode
end

function VecPow(arg1, arg2)
    @ccall libpetsc.VecPow(arg1::Vec, arg2::PetscScalar)::PetscErrorCode
end

function VecMedian(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecMedian(arg1::Vec, arg2::Vec, arg3::Vec, arg4::Vec)::PetscErrorCode
end

function VecWhichInactive(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.VecWhichInactive(arg1::Vec, arg2::Vec, arg3::Vec, arg4::Vec, arg5::PetscBool, arg6::Ptr{IS})::PetscErrorCode
end

function VecWhichBetween(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecWhichBetween(arg1::Vec, arg2::Vec, arg3::Vec, arg4::Ptr{IS})::PetscErrorCode
end

function VecWhichBetweenOrEqual(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecWhichBetweenOrEqual(arg1::Vec, arg2::Vec, arg3::Vec, arg4::Ptr{IS})::PetscErrorCode
end

function VecWhichGreaterThan(arg1, arg2, arg3)
    @ccall libpetsc.VecWhichGreaterThan(arg1::Vec, arg2::Vec, arg3::Ptr{IS})::PetscErrorCode
end

function VecWhichLessThan(arg1, arg2, arg3)
    @ccall libpetsc.VecWhichLessThan(arg1::Vec, arg2::Vec, arg3::Ptr{IS})::PetscErrorCode
end

function VecWhichEqual(arg1, arg2, arg3)
    @ccall libpetsc.VecWhichEqual(arg1::Vec, arg2::Vec, arg3::Ptr{IS})::PetscErrorCode
end

function VecISAXPY(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecISAXPY(arg1::Vec, arg2::IS, arg3::PetscScalar, arg4::Vec)::PetscErrorCode
end

function VecISCopy(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecISCopy(arg1::Vec, arg2::IS, arg3::ScatterMode, arg4::Vec)::PetscErrorCode
end

function VecISSet(arg1, arg2, arg3)
    @ccall libpetsc.VecISSet(arg1::Vec, arg2::IS, arg3::PetscScalar)::PetscErrorCode
end

function VecBoundGradientProjection(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.VecBoundGradientProjection(arg1::Vec, arg2::Vec, arg3::Vec, arg4::Vec, arg5::Vec)::PetscErrorCode
end

function VecStepBoundInfo(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.VecStepBoundInfo(arg1::Vec, arg2::Vec, arg3::Vec, arg4::Vec, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal}, arg7::Ptr{PetscReal})::PetscErrorCode
end

function VecStepMax(arg1, arg2, arg3)
    @ccall libpetsc.VecStepMax(arg1::Vec, arg2::Vec, arg3::Ptr{PetscReal})::PetscErrorCode
end

function VecStepMaxBounded(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.VecStepMaxBounded(arg1::Vec, arg2::Vec, arg3::Vec, arg4::Vec, arg5::Ptr{PetscReal})::PetscErrorCode
end

function PetscViewerMathematicaGetVector(arg1, arg2)
    @ccall libpetsc.PetscViewerMathematicaGetVector(arg1::PetscViewer, arg2::Vec)::PetscErrorCode
end

function PetscViewerMathematicaPutVector(arg1, arg2)
    @ccall libpetsc.PetscViewerMathematicaPutVector(arg1::PetscViewer, arg2::Vec)::PetscErrorCode
end

mutable struct _n_Vecs
    n::PetscInt
    v::Vec
    _n_Vecs() = new()
end

const Vecs = Ptr{_n_Vecs}

function VecsDestroy(arg1)
    @ccall libpetsc.VecsDestroy(arg1::Vecs)::PetscErrorCode
end

function VecsCreateSeq(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecsCreateSeq(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Vecs})::PetscErrorCode
end

function VecsCreateSeqWithArray(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.VecsCreateSeqWithArray(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscScalar}, arg5::Ptr{Vecs})::PetscErrorCode
end

function VecsDuplicate(arg1, arg2)
    @ccall libpetsc.VecsDuplicate(arg1::Vecs, arg2::Ptr{Vecs})::PetscErrorCode
end

function VecNestGetSubVecs(arg1, arg2, arg3)
    @ccall libpetsc.VecNestGetSubVecs(arg1::Vec, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{Vec}})::PetscErrorCode
end

function VecNestGetSubVec(arg1, arg2, arg3)
    @ccall libpetsc.VecNestGetSubVec(arg1::Vec, arg2::PetscInt, arg3::Ptr{Vec})::PetscErrorCode
end

function VecNestSetSubVecs(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecNestSetSubVecs(arg1::Vec, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{Vec})::PetscErrorCode
end

function VecNestSetSubVec(arg1, arg2, arg3)
    @ccall libpetsc.VecNestSetSubVec(arg1::Vec, arg2::PetscInt, arg3::Vec)::PetscErrorCode
end

function VecCreateNest(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.VecCreateNest(arg1::MPI_Comm, arg2::PetscInt, arg3::Ptr{IS}, arg4::Ptr{Vec}, arg5::Ptr{Vec})::PetscErrorCode
end

function VecNestGetSize(arg1, arg2)
    @ccall libpetsc.VecNestGetSize(arg1::Vec, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscOptionsGetVec(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscOptionsGetVec(arg1::PetscOptions, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Vec, arg5::Ptr{PetscBool})::PetscErrorCode
end

function VecChop(arg1, arg2)
    @ccall libpetsc.VecChop(arg1::Vec, arg2::PetscReal)::PetscErrorCode
end

function VecGetLayout(arg1, arg2)
    @ccall libpetsc.VecGetLayout(arg1::Vec, arg2::Ptr{PetscLayout})::PetscErrorCode
end

function VecSetLayout(arg1, arg2)
    @ccall libpetsc.VecSetLayout(arg1::Vec, arg2::PetscLayout)::PetscErrorCode
end

function PetscSectionVecView(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionVecView(arg1::PetscSection, arg2::Vec, arg3::PetscViewer)::PetscErrorCode
end

function VecGetValuesSection(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecGetValuesSection(arg1::Vec, arg2::PetscSection, arg3::PetscInt, arg4::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function VecSetValuesSection(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.VecSetValuesSection(arg1::Vec, arg2::PetscSection, arg3::PetscInt, arg4::Ptr{PetscScalar}, arg5::InsertMode)::PetscErrorCode
end

function PetscSectionVecNorm(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscSectionVecNorm(arg1::PetscSection, arg2::PetscSection, arg3::Vec, arg4::NormType, arg5::Ptr{PetscReal})::PetscErrorCode
end

mutable struct _p_VecTagger end

const VecTagger = Ptr{_p_VecTagger}

const VecTaggerType = Ptr{Cchar}

function VecTaggerRegister(arg1, arg2)
    @ccall libpetsc.VecTaggerRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function VecTaggerCreate(arg1, arg2)
    @ccall libpetsc.VecTaggerCreate(arg1::MPI_Comm, arg2::Ptr{VecTagger})::PetscErrorCode
end

function VecTaggerSetBlockSize(arg1, arg2)
    @ccall libpetsc.VecTaggerSetBlockSize(arg1::VecTagger, arg2::PetscInt)::PetscErrorCode
end

function VecTaggerGetBlockSize(arg1, arg2)
    @ccall libpetsc.VecTaggerGetBlockSize(arg1::VecTagger, arg2::Ptr{PetscInt})::PetscErrorCode
end

function VecTaggerSetType(arg1, arg2)
    @ccall libpetsc.VecTaggerSetType(arg1::VecTagger, arg2::VecTaggerType)::PetscErrorCode
end

function VecTaggerGetType(arg1, arg2)
    @ccall libpetsc.VecTaggerGetType(arg1::VecTagger, arg2::Ptr{VecTaggerType})::PetscErrorCode
end

function VecTaggerSetInvert(arg1, arg2)
    @ccall libpetsc.VecTaggerSetInvert(arg1::VecTagger, arg2::PetscBool)::PetscErrorCode
end

function VecTaggerGetInvert(arg1, arg2)
    @ccall libpetsc.VecTaggerGetInvert(arg1::VecTagger, arg2::Ptr{PetscBool})::PetscErrorCode
end

function VecTaggerSetFromOptions(arg1)
    @ccall libpetsc.VecTaggerSetFromOptions(arg1::VecTagger)::PetscErrorCode
end

function VecTaggerSetUp(arg1)
    @ccall libpetsc.VecTaggerSetUp(arg1::VecTagger)::PetscErrorCode
end

function VecTaggerView(arg1, arg2)
    @ccall libpetsc.VecTaggerView(arg1::VecTagger, arg2::PetscViewer)::PetscErrorCode
end

function VecTaggerComputeIS(arg1, arg2, arg3)
    @ccall libpetsc.VecTaggerComputeIS(arg1::VecTagger, arg2::Vec, arg3::Ptr{IS})::PetscErrorCode
end

function VecTaggerDestroy(arg1)
    @ccall libpetsc.VecTaggerDestroy(arg1::Ptr{VecTagger})::PetscErrorCode
end

mutable struct VecTaggerBox
    min::PetscScalar
    max::PetscScalar
    VecTaggerBox() = new()
end

function VecTaggerComputeBoxes(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecTaggerComputeBoxes(arg1::VecTagger, arg2::Vec, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{VecTaggerBox}})::PetscErrorCode
end

function VecTaggerAbsoluteSetBox(arg1, arg2)
    @ccall libpetsc.VecTaggerAbsoluteSetBox(arg1::VecTagger, arg2::Ptr{VecTaggerBox})::PetscErrorCode
end

function VecTaggerAbsoluteGetBox(arg1, arg2)
    @ccall libpetsc.VecTaggerAbsoluteGetBox(arg1::VecTagger, arg2::Ptr{Ptr{VecTaggerBox}})::PetscErrorCode
end

function VecTaggerRelativeSetBox(arg1, arg2)
    @ccall libpetsc.VecTaggerRelativeSetBox(arg1::VecTagger, arg2::Ptr{VecTaggerBox})::PetscErrorCode
end

function VecTaggerRelativeGetBox(arg1, arg2)
    @ccall libpetsc.VecTaggerRelativeGetBox(arg1::VecTagger, arg2::Ptr{Ptr{VecTaggerBox}})::PetscErrorCode
end

function VecTaggerCDFSetBox(arg1, arg2)
    @ccall libpetsc.VecTaggerCDFSetBox(arg1::VecTagger, arg2::Ptr{VecTaggerBox})::PetscErrorCode
end

function VecTaggerCDFGetBox(arg1, arg2)
    @ccall libpetsc.VecTaggerCDFGetBox(arg1::VecTagger, arg2::Ptr{Ptr{VecTaggerBox}})::PetscErrorCode
end

@enum VecTaggerCDFMethod::UInt32 begin
    VECTAGGER_CDF_GATHER = 0
    VECTAGGER_CDF_ITERATIVE = 1
    VECTAGGER_CDF_NUM_METHODS = 2
end

function VecTaggerCDFSetMethod(arg1, arg2)
    @ccall libpetsc.VecTaggerCDFSetMethod(arg1::VecTagger, arg2::VecTaggerCDFMethod)::PetscErrorCode
end

function VecTaggerCDFGetMethod(arg1, arg2)
    @ccall libpetsc.VecTaggerCDFGetMethod(arg1::VecTagger, arg2::Ptr{VecTaggerCDFMethod})::PetscErrorCode
end

function VecTaggerCDFIterativeSetTolerances(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecTaggerCDFIterativeSetTolerances(arg1::VecTagger, arg2::PetscInt, arg3::PetscReal, arg4::PetscReal)::PetscErrorCode
end

function VecTaggerCDFIterativeGetTolerances(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecTaggerCDFIterativeGetTolerances(arg1::VecTagger, arg2::Ptr{PetscInt}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal})::PetscErrorCode
end

function VecTaggerOrSetSubs(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecTaggerOrSetSubs(arg1::VecTagger, arg2::PetscInt, arg3::Ptr{VecTagger}, arg4::PetscCopyMode)::PetscErrorCode
end

function VecTaggerOrGetSubs(arg1, arg2, arg3)
    @ccall libpetsc.VecTaggerOrGetSubs(arg1::VecTagger, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{VecTagger}})::PetscErrorCode
end

function VecTaggerAndSetSubs(arg1, arg2, arg3, arg4)
    @ccall libpetsc.VecTaggerAndSetSubs(arg1::VecTagger, arg2::PetscInt, arg3::Ptr{VecTagger}, arg4::PetscCopyMode)::PetscErrorCode
end

function VecTaggerAndGetSubs(arg1, arg2, arg3)
    @ccall libpetsc.VecTaggerAndGetSubs(arg1::VecTagger, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{VecTagger}})::PetscErrorCode
end

function VecTaggerInitializePackage()
    @ccall libpetsc.VecTaggerInitializePackage()::PetscErrorCode
end

function VecTaggerFinalizePackage()
    @ccall libpetsc.VecTaggerFinalizePackage()::PetscErrorCode
end

@enum PetscSFPattern::UInt32 begin
    PETSCSF_PATTERN_GENERAL = 0
    PETSCSF_PATTERN_ALLGATHER = 1
    PETSCSF_PATTERN_GATHER = 2
    PETSCSF_PATTERN_ALLTOALL = 3
end

@enum PetscSFWindowSyncType::UInt32 begin
    PETSCSF_WINDOW_SYNC_FENCE = 0
    PETSCSF_WINDOW_SYNC_LOCK = 1
    PETSCSF_WINDOW_SYNC_ACTIVE = 2
end

@enum PetscSFWindowFlavorType::UInt32 begin
    PETSCSF_WINDOW_FLAVOR_CREATE = 0
    PETSCSF_WINDOW_FLAVOR_DYNAMIC = 1
    PETSCSF_WINDOW_FLAVOR_ALLOCATE = 2
    PETSCSF_WINDOW_FLAVOR_SHARED = 3
end

@enum PetscSFDuplicateOption::UInt32 begin
    PETSCSF_DUPLICATE_CONFONLY = 0
    PETSCSF_DUPLICATE_RANKS = 1
    PETSCSF_DUPLICATE_GRAPH = 2
end

function PetscSFRegister(arg1, arg2)
    @ccall libpetsc.PetscSFRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscSFInitializePackage()
    @ccall libpetsc.PetscSFInitializePackage()::PetscErrorCode
end

function PetscSFFinalizePackage()
    @ccall libpetsc.PetscSFFinalizePackage()::PetscErrorCode
end

function PetscSFCreate(arg1, arg2)
    @ccall libpetsc.PetscSFCreate(arg1::MPI_Comm, arg2::Ptr{PetscSF})::PetscErrorCode
end

function PetscSFDestroy(arg1)
    @ccall libpetsc.PetscSFDestroy(arg1::Ptr{PetscSF})::PetscErrorCode
end

function PetscSFSetType(arg1, arg2)
    @ccall libpetsc.PetscSFSetType(arg1::PetscSF, arg2::PetscSFType)::PetscErrorCode
end

function PetscSFGetType(arg1, arg2)
    @ccall libpetsc.PetscSFGetType(arg1::PetscSF, arg2::Ptr{PetscSFType})::PetscErrorCode
end

function PetscSFView(arg1, arg2)
    @ccall libpetsc.PetscSFView(arg1::PetscSF, arg2::PetscViewer)::PetscErrorCode
end

function PetscSFViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.PetscSFViewFromOptions(arg1::PetscSF, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function PetscSFSetUp(arg1)
    @ccall libpetsc.PetscSFSetUp(arg1::PetscSF)::PetscErrorCode
end

function PetscSFSetFromOptions(arg1)
    @ccall libpetsc.PetscSFSetFromOptions(arg1::PetscSF)::PetscErrorCode
end

function PetscSFDuplicate(arg1, arg2, arg3)
    @ccall libpetsc.PetscSFDuplicate(arg1::PetscSF, arg2::PetscSFDuplicateOption, arg3::Ptr{PetscSF})::PetscErrorCode
end

function PetscSFWindowSetSyncType(arg1, arg2)
    @ccall libpetsc.PetscSFWindowSetSyncType(arg1::PetscSF, arg2::PetscSFWindowSyncType)::PetscErrorCode
end

function PetscSFWindowGetSyncType(arg1, arg2)
    @ccall libpetsc.PetscSFWindowGetSyncType(arg1::PetscSF, arg2::Ptr{PetscSFWindowSyncType})::PetscErrorCode
end

function PetscSFWindowSetFlavorType(arg1, arg2)
    @ccall libpetsc.PetscSFWindowSetFlavorType(arg1::PetscSF, arg2::PetscSFWindowFlavorType)::PetscErrorCode
end

function PetscSFWindowGetFlavorType(arg1, arg2)
    @ccall libpetsc.PetscSFWindowGetFlavorType(arg1::PetscSF, arg2::Ptr{PetscSFWindowFlavorType})::PetscErrorCode
end

function PetscSFWindowSetInfo(arg1, arg2)
    @ccall libpetsc.PetscSFWindowSetInfo(arg1::PetscSF, arg2::MPI_Info)::PetscErrorCode
end

function PetscSFWindowGetInfo(arg1, arg2)
    @ccall libpetsc.PetscSFWindowGetInfo(arg1::PetscSF, arg2::Ptr{MPI_Info})::PetscErrorCode
end

function PetscSFSetRankOrder(arg1, arg2)
    @ccall libpetsc.PetscSFSetRankOrder(arg1::PetscSF, arg2::PetscBool)::PetscErrorCode
end

function PetscSFSetGraph(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscSFSetGraph(arg1::PetscSF, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::PetscCopyMode, arg6::Ptr{PetscSFNode}, arg7::PetscCopyMode)::PetscErrorCode
end

function PetscSFSetGraphWithPattern(arg1, arg2, arg3)
    @ccall libpetsc.PetscSFSetGraphWithPattern(arg1::PetscSF, arg2::PetscLayout, arg3::PetscSFPattern)::PetscErrorCode
end

function PetscSFGetGraph(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscSFGetGraph(arg1::PetscSF, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{PetscInt}}, arg5::Ptr{Ptr{PetscSFNode}})::PetscErrorCode
end

function PetscSFGetLeafRange(arg1, arg2, arg3)
    @ccall libpetsc.PetscSFGetLeafRange(arg1::PetscSF, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscSFCreateEmbeddedRootSF(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSFCreateEmbeddedRootSF(arg1::PetscSF, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscSF})::PetscErrorCode
end

function PetscSFCreateEmbeddedLeafSF(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSFCreateEmbeddedLeafSF(arg1::PetscSF, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscSF})::PetscErrorCode
end

function PetscSFReset(arg1)
    @ccall libpetsc.PetscSFReset(arg1::PetscSF)::PetscErrorCode
end

function PetscSFSetUpRanks(arg1, arg2)
    @ccall libpetsc.PetscSFSetUpRanks(arg1::PetscSF, arg2::MPI_Group)::PetscErrorCode
end

function PetscSFGetRootRanks(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscSFGetRootRanks(arg1::PetscSF, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{PetscMPIInt}}, arg4::Ptr{Ptr{PetscInt}}, arg5::Ptr{Ptr{PetscInt}}, arg6::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function PetscSFGetLeafRanks(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscSFGetLeafRanks(arg1::PetscSF, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{PetscMPIInt}}, arg4::Ptr{Ptr{PetscInt}}, arg5::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function PetscSFGetGroups(arg1, arg2, arg3)
    @ccall libpetsc.PetscSFGetGroups(arg1::PetscSF, arg2::Ptr{MPI_Group}, arg3::Ptr{MPI_Group})::PetscErrorCode
end

function PetscSFGetMultiSF(arg1, arg2)
    @ccall libpetsc.PetscSFGetMultiSF(arg1::PetscSF, arg2::Ptr{PetscSF})::PetscErrorCode
end

function PetscSFCreateInverseSF(arg1, arg2)
    @ccall libpetsc.PetscSFCreateInverseSF(arg1::PetscSF, arg2::Ptr{PetscSF})::PetscErrorCode
end

function PetscSFSetGraphLayout(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscSFSetGraphLayout(arg1::PetscSF, arg2::PetscLayout, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::PetscCopyMode, arg6::Ptr{PetscInt})::PetscErrorCode
end

function PetscSFCreateFromLayouts(arg1, arg2, arg3)
    @ccall libpetsc.PetscSFCreateFromLayouts(arg1::PetscLayout, arg2::PetscLayout, arg3::Ptr{PetscSF})::PetscErrorCode
end

function PetscLayoutsCreateSF(rmap, lmap, sf)
    @ccall libpetsc.PetscLayoutsCreateSF(rmap::PetscLayout, lmap::PetscLayout, sf::Ptr{PetscSF})::PetscErrorCode
end

function PetscSFCreateByMatchingIndices(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11)
    @ccall libpetsc.PetscSFCreateByMatchingIndices(arg1::PetscLayout, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{PetscInt}, arg8::Ptr{PetscInt}, arg9::PetscInt, arg10::Ptr{PetscSF}, arg11::Ptr{PetscSF})::PetscErrorCode
end

function PetscSFSetGraphSection(arg1, arg2, arg3)
    @ccall libpetsc.PetscSFSetGraphSection(arg1::PetscSF, arg2::PetscSection, arg3::PetscSection)::PetscErrorCode
end

function PetscSFCreateRemoteOffsets(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSFCreateRemoteOffsets(arg1::PetscSF, arg2::PetscSection, arg3::PetscSection, arg4::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function PetscSFDistributeSection(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSFDistributeSection(arg1::PetscSF, arg2::PetscSection, arg3::Ptr{Ptr{PetscInt}}, arg4::PetscSection)::PetscErrorCode
end

function PetscSFCreateSectionSF(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscSFCreateSectionSF(arg1::PetscSF, arg2::PetscSection, arg3::Ptr{PetscInt}, arg4::PetscSection, arg5::Ptr{PetscSF})::PetscErrorCode
end

function PetscSFBcastBegin(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscSFBcastBegin(arg1::PetscSF, arg2::MPI_Datatype, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid}, arg5::MPI_Op)::PetscErrorCode
end

function PetscSFBcastEnd(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscSFBcastEnd(arg1::PetscSF, arg2::MPI_Datatype, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid}, arg5::MPI_Op)::PetscErrorCode
end

function PetscSFBcastWithMemTypeBegin(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscSFBcastWithMemTypeBegin(arg1::PetscSF, arg2::MPI_Datatype, arg3::PetscMemType, arg4::Ptr{Cvoid}, arg5::PetscMemType, arg6::Ptr{Cvoid}, arg7::MPI_Op)::PetscErrorCode
end

function PetscSFReduceBegin(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscSFReduceBegin(arg1::PetscSF, arg2::MPI_Datatype, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid}, arg5::MPI_Op)::PetscErrorCode
end

function PetscSFReduceEnd(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscSFReduceEnd(arg1::PetscSF, arg2::MPI_Datatype, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid}, arg5::MPI_Op)::PetscErrorCode
end

function PetscSFReduceWithMemTypeBegin(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscSFReduceWithMemTypeBegin(arg1::PetscSF, arg2::MPI_Datatype, arg3::PetscMemType, arg4::Ptr{Cvoid}, arg5::PetscMemType, arg6::Ptr{Cvoid}, arg7::MPI_Op)::PetscErrorCode
end

function PetscSFFetchAndOpBegin(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscSFFetchAndOpBegin(arg1::PetscSF, arg2::MPI_Datatype, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid}, arg6::MPI_Op)::PetscErrorCode
end

function PetscSFFetchAndOpEnd(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscSFFetchAndOpEnd(arg1::PetscSF, arg2::MPI_Datatype, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid}, arg6::MPI_Op)::PetscErrorCode
end

function PetscSFComputeDegreeBegin(arg1, arg2)
    @ccall libpetsc.PetscSFComputeDegreeBegin(arg1::PetscSF, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function PetscSFComputeDegreeEnd(arg1, arg2)
    @ccall libpetsc.PetscSFComputeDegreeEnd(arg1::PetscSF, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function PetscSFComputeMultiRootOriginalNumbering(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSFComputeMultiRootOriginalNumbering(arg1::PetscSF, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function PetscSFGatherBegin(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSFGatherBegin(arg1::PetscSF, arg2::MPI_Datatype, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function PetscSFGatherEnd(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSFGatherEnd(arg1::PetscSF, arg2::MPI_Datatype, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function PetscSFScatterBegin(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSFScatterBegin(arg1::PetscSF, arg2::MPI_Datatype, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function PetscSFScatterEnd(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSFScatterEnd(arg1::PetscSF, arg2::MPI_Datatype, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function PetscSFCompose(arg1, arg2, arg3)
    @ccall libpetsc.PetscSFCompose(arg1::PetscSF, arg2::PetscSF, arg3::Ptr{PetscSF})::PetscErrorCode
end

function PetscSFComposeInverse(arg1, arg2, arg3)
    @ccall libpetsc.PetscSFComposeInverse(arg1::PetscSF, arg2::PetscSF, arg3::Ptr{PetscSF})::PetscErrorCode
end

function PetscSFGetRanks(sf, nranks, ranks, roffset, rmine, rremote)
    @ccall libpetsc.PetscSFGetRanks(sf::PetscSF, nranks::Ptr{PetscInt}, ranks::Ptr{Ptr{PetscMPIInt}}, roffset::Ptr{Ptr{PetscInt}}, rmine::Ptr{Ptr{PetscInt}}, rremote::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function PetscSFCreateEmbeddedSF(sf, nselected, selected, esf)
    @ccall libpetsc.PetscSFCreateEmbeddedSF(sf::PetscSF, nselected::PetscInt, selected::Ptr{PetscInt}, esf::Ptr{PetscSF})::PetscErrorCode
end

function PetscSFBcastAndOpBegin(sf, unit, rootdata, leafdata, op)
    @ccall libpetsc.PetscSFBcastAndOpBegin(sf::PetscSF, unit::MPI_Datatype, rootdata::Ptr{Cvoid}, leafdata::Ptr{Cvoid}, op::MPI_Op)::PetscErrorCode
end

function PetscSFBcastAndOpEnd(sf, unit, rootdata, leafdata, op)
    @ccall libpetsc.PetscSFBcastAndOpEnd(sf::PetscSF, unit::MPI_Datatype, rootdata::Ptr{Cvoid}, leafdata::Ptr{Cvoid}, op::MPI_Op)::PetscErrorCode
end

function PetscSFBcastAndOpWithMemtypeBegin(sf, unit, rootmtype, rootdata, leafmtype, leafdata, op)
    @ccall libpetsc.PetscSFBcastAndOpWithMemtypeBegin(sf::PetscSF, unit::MPI_Datatype, rootmtype::PetscMemType, rootdata::Ptr{Cvoid}, leafmtype::PetscMemType, leafdata::Ptr{Cvoid}, op::MPI_Op)::PetscErrorCode
end

function PetscSectionCreate(arg1, arg2)
    @ccall libpetsc.PetscSectionCreate(arg1::MPI_Comm, arg2::Ptr{PetscSection})::PetscErrorCode
end

function PetscSectionClone(arg1, arg2)
    @ccall libpetsc.PetscSectionClone(arg1::PetscSection, arg2::Ptr{PetscSection})::PetscErrorCode
end

function PetscSectionSetFromOptions(arg1)
    @ccall libpetsc.PetscSectionSetFromOptions(arg1::PetscSection)::PetscErrorCode
end

function PetscSectionCopy(arg1, arg2)
    @ccall libpetsc.PetscSectionCopy(arg1::PetscSection, arg2::PetscSection)::PetscErrorCode
end

function PetscSectionCompare(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionCompare(arg1::PetscSection, arg2::PetscSection, arg3::Ptr{PetscBool})::PetscErrorCode
end

function PetscSectionGetNumFields(arg1, arg2)
    @ccall libpetsc.PetscSectionGetNumFields(arg1::PetscSection, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscSectionSetNumFields(arg1, arg2)
    @ccall libpetsc.PetscSectionSetNumFields(arg1::PetscSection, arg2::PetscInt)::PetscErrorCode
end

function PetscSectionGetFieldName(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionGetFieldName(arg1::PetscSection, arg2::PetscInt, arg3::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscSectionSetFieldName(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionSetFieldName(arg1::PetscSection, arg2::PetscInt, arg3::Ptr{Cchar})::PetscErrorCode
end

function PetscSectionGetComponentName(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSectionGetComponentName(arg1::PetscSection, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscSectionSetComponentName(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSectionSetComponentName(arg1::PetscSection, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Cchar})::PetscErrorCode
end

function PetscSectionGetFieldComponents(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionGetFieldComponents(arg1::PetscSection, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscSectionSetFieldComponents(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionSetFieldComponents(arg1::PetscSection, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function PetscSectionGetChart(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionGetChart(arg1::PetscSection, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscSectionSetChart(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionSetChart(arg1::PetscSection, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function PetscSectionGetPermutation(arg1, arg2)
    @ccall libpetsc.PetscSectionGetPermutation(arg1::PetscSection, arg2::Ptr{IS})::PetscErrorCode
end

function PetscSectionSetPermutation(arg1, arg2)
    @ccall libpetsc.PetscSectionSetPermutation(arg1::PetscSection, arg2::IS)::PetscErrorCode
end

function PetscSectionGetPointMajor(arg1, arg2)
    @ccall libpetsc.PetscSectionGetPointMajor(arg1::PetscSection, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscSectionSetPointMajor(arg1, arg2)
    @ccall libpetsc.PetscSectionSetPointMajor(arg1::PetscSection, arg2::PetscBool)::PetscErrorCode
end

function PetscSectionGetDof(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionGetDof(arg1::PetscSection, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscSectionSetDof(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionSetDof(arg1::PetscSection, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function PetscSectionAddDof(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionAddDof(arg1::PetscSection, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function PetscSectionGetFieldDof(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSectionGetFieldDof(arg1::PetscSection, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt})::PetscErrorCode
end

function PetscSectionSetFieldDof(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSectionSetFieldDof(arg1::PetscSection, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt)::PetscErrorCode
end

function PetscSectionAddFieldDof(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSectionAddFieldDof(arg1::PetscSection, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt)::PetscErrorCode
end

function PetscSectionHasConstraints(arg1, arg2)
    @ccall libpetsc.PetscSectionHasConstraints(arg1::PetscSection, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscSectionGetConstraintDof(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionGetConstraintDof(arg1::PetscSection, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscSectionSetConstraintDof(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionSetConstraintDof(arg1::PetscSection, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function PetscSectionAddConstraintDof(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionAddConstraintDof(arg1::PetscSection, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function PetscSectionGetFieldConstraintDof(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSectionGetFieldConstraintDof(arg1::PetscSection, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt})::PetscErrorCode
end

function PetscSectionSetFieldConstraintDof(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSectionSetFieldConstraintDof(arg1::PetscSection, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt)::PetscErrorCode
end

function PetscSectionAddFieldConstraintDof(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSectionAddFieldConstraintDof(arg1::PetscSection, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt)::PetscErrorCode
end

function PetscSectionGetConstraintIndices(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionGetConstraintIndices(arg1::PetscSection, arg2::PetscInt, arg3::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function PetscSectionSetConstraintIndices(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionSetConstraintIndices(arg1::PetscSection, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscSectionGetFieldConstraintIndices(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSectionGetFieldConstraintIndices(arg1::PetscSection, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function PetscSectionSetFieldConstraintIndices(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSectionSetFieldConstraintIndices(arg1::PetscSection, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt})::PetscErrorCode
end

function PetscSectionSetUpBC(arg1)
    @ccall libpetsc.PetscSectionSetUpBC(arg1::PetscSection)::PetscErrorCode
end

function PetscSectionSetUp(arg1)
    @ccall libpetsc.PetscSectionSetUp(arg1::PetscSection)::PetscErrorCode
end

function PetscSectionGetMaxDof(arg1, arg2)
    @ccall libpetsc.PetscSectionGetMaxDof(arg1::PetscSection, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscSectionGetStorageSize(arg1, arg2)
    @ccall libpetsc.PetscSectionGetStorageSize(arg1::PetscSection, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscSectionGetConstrainedStorageSize(arg1, arg2)
    @ccall libpetsc.PetscSectionGetConstrainedStorageSize(arg1::PetscSection, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscSectionGetOffset(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionGetOffset(arg1::PetscSection, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscSectionSetOffset(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionSetOffset(arg1::PetscSection, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function PetscSectionGetFieldOffset(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSectionGetFieldOffset(arg1::PetscSection, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt})::PetscErrorCode
end

function PetscSectionSetFieldOffset(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSectionSetFieldOffset(arg1::PetscSection, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt)::PetscErrorCode
end

function PetscSectionGetFieldPointOffset(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSectionGetFieldPointOffset(arg1::PetscSection, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt})::PetscErrorCode
end

function PetscSectionGetOffsetRange(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionGetOffsetRange(arg1::PetscSection, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscSectionView(arg1, arg2)
    @ccall libpetsc.PetscSectionView(arg1::PetscSection, arg2::PetscViewer)::PetscErrorCode
end

function PetscSectionViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionViewFromOptions(arg1::PetscSection, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function PetscSectionReset(arg1)
    @ccall libpetsc.PetscSectionReset(arg1::PetscSection)::PetscErrorCode
end

function PetscSectionDestroy(arg1)
    @ccall libpetsc.PetscSectionDestroy(arg1::Ptr{PetscSection})::PetscErrorCode
end

function PetscSectionCreateGlobalSection(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscSectionCreateGlobalSection(arg1::PetscSection, arg2::PetscSF, arg3::PetscBool, arg4::PetscBool, arg5::Ptr{PetscSection})::PetscErrorCode
end

function PetscSectionCreateGlobalSectionCensored(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscSectionCreateGlobalSectionCensored(arg1::PetscSection, arg2::PetscSF, arg3::PetscBool, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{PetscSection})::PetscErrorCode
end

function PetscSectionCreateSubsection(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSectionCreateSubsection(arg1::PetscSection, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscSection})::PetscErrorCode
end

function PetscSectionCreateSupersection(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionCreateSupersection(arg1::Ptr{PetscSection}, arg2::PetscInt, arg3::Ptr{PetscSection})::PetscErrorCode
end

function PetscSectionCreateSubmeshSection(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionCreateSubmeshSection(arg1::PetscSection, arg2::IS, arg3::Ptr{PetscSection})::PetscErrorCode
end

function PetscSectionGetPointLayout(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionGetPointLayout(arg1::MPI_Comm, arg2::PetscSection, arg3::Ptr{PetscLayout})::PetscErrorCode
end

function PetscSectionGetValueLayout(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionGetValueLayout(arg1::MPI_Comm, arg2::PetscSection, arg3::Ptr{PetscLayout})::PetscErrorCode
end

function PetscSectionPermute(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionPermute(arg1::PetscSection, arg2::IS, arg3::Ptr{PetscSection})::PetscErrorCode
end

function PetscSectionGetField(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionGetField(arg1::PetscSection, arg2::PetscInt, arg3::Ptr{PetscSection})::PetscErrorCode
end

function PetscSectionSetUseFieldOffsets(arg1, arg2)
    @ccall libpetsc.PetscSectionSetUseFieldOffsets(arg1::PetscSection, arg2::PetscBool)::PetscErrorCode
end

function PetscSectionGetUseFieldOffsets(arg1, arg2)
    @ccall libpetsc.PetscSectionGetUseFieldOffsets(arg1::PetscSection, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscSectionExtractDofsFromArray(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscSectionExtractDofsFromArray(arg1::PetscSection, arg2::MPI_Datatype, arg3::Ptr{Cvoid}, arg4::IS, arg5::Ptr{PetscSection}, arg6::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscSectionSetClosureIndex(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSectionSetClosureIndex(arg1::PetscSection, arg2::PetscObject, arg3::PetscSection, arg4::IS)::PetscErrorCode
end

function PetscSectionGetClosureIndex(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSectionGetClosureIndex(arg1::PetscSection, arg2::PetscObject, arg3::Ptr{PetscSection}, arg4::Ptr{IS})::PetscErrorCode
end

function PetscSectionSetClosurePermutation(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscSectionSetClosurePermutation(arg1::PetscSection, arg2::PetscObject, arg3::PetscInt, arg4::IS)::PetscErrorCode
end

function PetscSectionGetClosurePermutation(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscSectionGetClosurePermutation(arg1::PetscSection, arg2::PetscObject, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{IS})::PetscErrorCode
end

function PetscSectionGetClosureInversePermutation(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscSectionGetClosureInversePermutation(arg1::PetscSection, arg2::PetscObject, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{IS})::PetscErrorCode
end

function PetscSectionSymSetType(arg1, arg2)
    @ccall libpetsc.PetscSectionSymSetType(arg1::PetscSectionSym, arg2::PetscSectionSymType)::PetscErrorCode
end

function PetscSectionSymGetType(arg1, arg2)
    @ccall libpetsc.PetscSectionSymGetType(arg1::PetscSectionSym, arg2::Ptr{PetscSectionSymType})::PetscErrorCode
end

function PetscSectionSymRegister(arg1, arg2)
    @ccall libpetsc.PetscSectionSymRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscSectionSymCreate(arg1, arg2)
    @ccall libpetsc.PetscSectionSymCreate(arg1::MPI_Comm, arg2::Ptr{PetscSectionSym})::PetscErrorCode
end

function PetscSectionSymDestroy(arg1)
    @ccall libpetsc.PetscSectionSymDestroy(arg1::Ptr{PetscSectionSym})::PetscErrorCode
end

function PetscSectionSymView(arg1, arg2)
    @ccall libpetsc.PetscSectionSymView(arg1::PetscSectionSym, arg2::PetscViewer)::PetscErrorCode
end

function PetscSectionSetSym(arg1, arg2)
    @ccall libpetsc.PetscSectionSetSym(arg1::PetscSection, arg2::PetscSectionSym)::PetscErrorCode
end

function PetscSectionGetSym(arg1, arg2)
    @ccall libpetsc.PetscSectionGetSym(arg1::PetscSection, arg2::Ptr{PetscSectionSym})::PetscErrorCode
end

function PetscSectionSetFieldSym(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionSetFieldSym(arg1::PetscSection, arg2::PetscInt, arg3::PetscSectionSym)::PetscErrorCode
end

function PetscSectionGetFieldSym(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionGetFieldSym(arg1::PetscSection, arg2::PetscInt, arg3::Ptr{PetscSectionSym})::PetscErrorCode
end

function PetscSectionGetPointSyms(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscSectionGetPointSyms(arg1::PetscSection, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{Ptr{PetscInt}}}, arg5::Ptr{Ptr{Ptr{PetscScalar}}})::PetscErrorCode
end

function PetscSectionRestorePointSyms(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscSectionRestorePointSyms(arg1::PetscSection, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{Ptr{PetscInt}}}, arg5::Ptr{Ptr{Ptr{PetscScalar}}})::PetscErrorCode
end

function PetscSectionGetFieldPointSyms(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscSectionGetFieldPointSyms(arg1::PetscSection, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::Ptr{Ptr{Ptr{PetscInt}}}, arg6::Ptr{Ptr{Ptr{PetscScalar}}})::PetscErrorCode
end

function PetscSectionRestoreFieldPointSyms(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscSectionRestoreFieldPointSyms(arg1::PetscSection, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::Ptr{Ptr{Ptr{PetscInt}}}, arg6::Ptr{Ptr{Ptr{PetscScalar}}})::PetscErrorCode
end

mutable struct _p_Mat end

const Mat = Ptr{_p_Mat}

const MatType = Ptr{Cchar}

const MatSolverType = Ptr{Cchar}

@enum MatFactorType::UInt32 begin
    MAT_FACTOR_NONE = 0
    MAT_FACTOR_LU = 1
    MAT_FACTOR_CHOLESKY = 2
    MAT_FACTOR_ILU = 3
    MAT_FACTOR_ICC = 4
    MAT_FACTOR_ILUDT = 5
    MAT_FACTOR_QR = 6
    MAT_FACTOR_NUM_TYPES = 7
end

function MatGetFactor(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatGetFactor(arg1::Mat, arg2::MatSolverType, arg3::MatFactorType, arg4::Ptr{Mat})::PetscErrorCode
end

function MatGetFactorAvailable(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatGetFactorAvailable(arg1::Mat, arg2::MatSolverType, arg3::MatFactorType, arg4::Ptr{PetscBool})::PetscErrorCode
end

function MatFactorGetUseOrdering(arg1, arg2)
    @ccall libpetsc.MatFactorGetUseOrdering(arg1::Mat, arg2::Ptr{PetscBool})::PetscErrorCode
end

function MatFactorGetSolverType(arg1, arg2)
    @ccall libpetsc.MatFactorGetSolverType(arg1::Mat, arg2::Ptr{MatSolverType})::PetscErrorCode
end

function MatGetFactorType(arg1, arg2)
    @ccall libpetsc.MatGetFactorType(arg1::Mat, arg2::Ptr{MatFactorType})::PetscErrorCode
end

function MatSetFactorType(arg1, arg2)
    @ccall libpetsc.MatSetFactorType(arg1::Mat, arg2::MatFactorType)::PetscErrorCode
end

# typedef PetscErrorCode ( * MatSolverFunction ) ( Mat , MatFactorType , Mat * )
const MatSolverFunction = Ptr{Cvoid}

function MatSolverTypeRegister(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatSolverTypeRegister(arg1::MatSolverType, arg2::MatType, arg3::MatFactorType, arg4::MatSolverFunction)::PetscErrorCode
end

function MatSolverTypeGet(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatSolverTypeGet(arg1::MatSolverType, arg2::MatType, arg3::MatFactorType, arg4::Ptr{PetscBool}, arg5::Ptr{PetscBool}, arg6::Ptr{MatSolverFunction})::PetscErrorCode
end

const MatSolverPackage = MatSolverType

function MatSolverPackageRegister(stype, mtype, ftype, f)
    @ccall libpetsc.MatSolverPackageRegister(stype::MatSolverType, mtype::MatType, ftype::MatFactorType, f::MatSolverFunction)::PetscErrorCode
end

function MatSolverPackageGet(stype, mtype, ftype, foundmtype, foundstype, f)
    @ccall libpetsc.MatSolverPackageGet(stype::MatSolverType, mtype::MatType, ftype::MatFactorType, foundmtype::Ptr{PetscBool}, foundstype::Ptr{PetscBool}, f::Ptr{MatSolverFunction})::PetscErrorCode
end

@enum MatProductType::UInt32 begin
    MATPRODUCT_UNSPECIFIED = 0
    MATPRODUCT_AB = 1
    MATPRODUCT_AtB = 2
    MATPRODUCT_ABt = 3
    MATPRODUCT_PtAP = 4
    MATPRODUCT_RARt = 5
    MATPRODUCT_ABC = 6
end

const MatProductAlgorithm = Ptr{Cchar}

function MatProductCreate(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatProductCreate(arg1::Mat, arg2::Mat, arg3::Mat, arg4::Ptr{Mat})::PetscErrorCode
end

function MatProductCreateWithMat(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatProductCreateWithMat(arg1::Mat, arg2::Mat, arg3::Mat, arg4::Mat)::PetscErrorCode
end

function MatProductSetType(arg1, arg2)
    @ccall libpetsc.MatProductSetType(arg1::Mat, arg2::MatProductType)::PetscErrorCode
end

function MatProductSetAlgorithm(arg1, arg2)
    @ccall libpetsc.MatProductSetAlgorithm(arg1::Mat, arg2::MatProductAlgorithm)::PetscErrorCode
end

function MatProductSetFill(arg1, arg2)
    @ccall libpetsc.MatProductSetFill(arg1::Mat, arg2::PetscReal)::PetscErrorCode
end

function MatProductSetFromOptions(arg1)
    @ccall libpetsc.MatProductSetFromOptions(arg1::Mat)::PetscErrorCode
end

function MatProductSymbolic(arg1)
    @ccall libpetsc.MatProductSymbolic(arg1::Mat)::PetscErrorCode
end

function MatProductNumeric(arg1)
    @ccall libpetsc.MatProductNumeric(arg1::Mat)::PetscErrorCode
end

function MatProductReplaceMats(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatProductReplaceMats(arg1::Mat, arg2::Mat, arg3::Mat, arg4::Mat)::PetscErrorCode
end

function MatProductClear(arg1)
    @ccall libpetsc.MatProductClear(arg1::Mat)::PetscErrorCode
end

function MatProductView(arg1, arg2)
    @ccall libpetsc.MatProductView(arg1::Mat, arg2::PetscViewer)::PetscErrorCode
end

@enum MatReuse::UInt32 begin
    MAT_INITIAL_MATRIX = 0
    MAT_REUSE_MATRIX = 1
    MAT_IGNORE_MATRIX = 2
    MAT_INPLACE_MATRIX = 3
end

@enum MatCreateSubMatrixOption::UInt32 begin
    MAT_DO_NOT_GET_VALUES = 0
    MAT_GET_VALUES = 1
end

function MatInitializePackage()
    @ccall libpetsc.MatInitializePackage()::PetscErrorCode
end

function MatCreate(arg1, arg2)
    @ccall libpetsc.MatCreate(arg1::MPI_Comm, arg2::Ptr{Mat})::PetscErrorCode
end

function MatSetSizes(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatSetSizes(arg1::Mat, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt)::PetscErrorCode
end

function MatSetType(arg1, arg2)
    @ccall libpetsc.MatSetType(arg1::Mat, arg2::MatType)::PetscErrorCode
end

function MatGetVecType(arg1, arg2)
    @ccall libpetsc.MatGetVecType(arg1::Mat, arg2::Ptr{VecType})::PetscErrorCode
end

function MatSetVecType(arg1, arg2)
    @ccall libpetsc.MatSetVecType(arg1::Mat, arg2::VecType)::PetscErrorCode
end

function MatSetFromOptions(arg1)
    @ccall libpetsc.MatSetFromOptions(arg1::Mat)::PetscErrorCode
end

function MatViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.MatViewFromOptions(arg1::Mat, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function MatRegister(arg1, arg2)
    @ccall libpetsc.MatRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function MatRegisterRootName(arg1, arg2, arg3)
    @ccall libpetsc.MatRegisterRootName(arg1::Ptr{Cchar}, arg2::Ptr{Cchar}, arg3::Ptr{Cchar})::PetscErrorCode
end

function MatSetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.MatSetOptionsPrefix(arg1::Mat, arg2::Ptr{Cchar})::PetscErrorCode
end

function MatAppendOptionsPrefix(arg1, arg2)
    @ccall libpetsc.MatAppendOptionsPrefix(arg1::Mat, arg2::Ptr{Cchar})::PetscErrorCode
end

function MatGetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.MatGetOptionsPrefix(arg1::Mat, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function MatSetErrorIfFailure(arg1, arg2)
    @ccall libpetsc.MatSetErrorIfFailure(arg1::Mat, arg2::PetscBool)::PetscErrorCode
end

@enum MatStructure::UInt32 begin
    DIFFERENT_NONZERO_PATTERN = 0
    SUBSET_NONZERO_PATTERN = 1
    SAME_NONZERO_PATTERN = 2
    UNKNOWN_NONZERO_PATTERN = 3
end

function MatCreateSeqSELL(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatCreateSeqSELL(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{Mat})::PetscErrorCode
end

function MatCreateSELL(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.MatCreateSELL(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{PetscInt}, arg8::PetscInt, arg9::Ptr{PetscInt}, arg10::Ptr{Mat})::PetscErrorCode
end

function MatSeqSELLSetPreallocation(arg1, arg2, arg3)
    @ccall libpetsc.MatSeqSELLSetPreallocation(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function MatMPISELLSetPreallocation(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatMPISELLSetPreallocation(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::PetscInt, arg5::Ptr{PetscInt})::PetscErrorCode
end

function MatCreateSeqDense(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatCreateSeqDense(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscScalar}, arg5::Ptr{Mat})::PetscErrorCode
end

function MatCreateDense(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.MatCreateDense(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{PetscScalar}, arg7::Ptr{Mat})::PetscErrorCode
end

function MatCreateSeqAIJ(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatCreateSeqAIJ(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{Mat})::PetscErrorCode
end

function MatCreateAIJ(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.MatCreateAIJ(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{PetscInt}, arg8::PetscInt, arg9::Ptr{PetscInt}, arg10::Ptr{Mat})::PetscErrorCode
end

function MatCreateMPIAIJWithArrays(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.MatCreateMPIAIJWithArrays(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::Ptr{PetscInt}, arg8::Ptr{PetscScalar}, arg9::Ptr{Mat})::PetscErrorCode
end

function MatUpdateMPIAIJWithArrays(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.MatUpdateMPIAIJWithArrays(arg1::Mat, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::Ptr{PetscInt}, arg8::Ptr{PetscScalar})::PetscErrorCode
end

function MatCreateMPIAIJWithSplitArrays(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12)
    @ccall libpetsc.MatCreateMPIAIJWithSplitArrays(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::Ptr{PetscInt}, arg8::Ptr{PetscScalar}, arg9::Ptr{PetscInt}, arg10::Ptr{PetscInt}, arg11::Ptr{PetscScalar}, arg12::Ptr{Mat})::PetscErrorCode
end

function MatCreateMPIAIJWithSeqAIJ(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatCreateMPIAIJWithSeqAIJ(arg1::MPI_Comm, arg2::Mat, arg3::Mat, arg4::Ptr{PetscInt}, arg5::Ptr{Mat})::PetscErrorCode
end

function MatCreateSeqBAIJ(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.MatCreateSeqBAIJ(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::Ptr{Mat})::PetscErrorCode
end

function MatCreateBAIJ(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11)
    @ccall libpetsc.MatCreateBAIJ(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::PetscInt, arg8::Ptr{PetscInt}, arg9::PetscInt, arg10::Ptr{PetscInt}, arg11::Ptr{Mat})::PetscErrorCode
end

function MatCreateMPIBAIJWithArrays(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.MatCreateMPIBAIJWithArrays(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{PetscInt}, arg8::Ptr{PetscInt}, arg9::Ptr{PetscScalar}, arg10::Ptr{Mat})::PetscErrorCode
end

function MatSetPreallocationCOO(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatSetPreallocationCOO(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function MatSetValuesCOO(arg1, arg2, arg3)
    @ccall libpetsc.MatSetValuesCOO(arg1::Mat, arg2::Ptr{PetscScalar}, arg3::InsertMode)::PetscErrorCode
end

function MatCreateMPIAdj(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.MatCreateMPIAdj(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscInt}, arg7::Ptr{Mat})::PetscErrorCode
end

function MatCreateSeqSBAIJ(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.MatCreateSeqSBAIJ(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::Ptr{Mat})::PetscErrorCode
end

function MatCreateSBAIJ(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11)
    @ccall libpetsc.MatCreateSBAIJ(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::PetscInt, arg8::Ptr{PetscInt}, arg9::PetscInt, arg10::Ptr{PetscInt}, arg11::Ptr{Mat})::PetscErrorCode
end

function MatCreateMPISBAIJWithArrays(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.MatCreateMPISBAIJWithArrays(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{PetscInt}, arg8::Ptr{PetscInt}, arg9::Ptr{PetscScalar}, arg10::Ptr{Mat})::PetscErrorCode
end

function MatSeqSBAIJSetPreallocationCSR(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatSeqSBAIJSetPreallocationCSR(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscScalar})::PetscErrorCode
end

function MatMPISBAIJSetPreallocationCSR(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatMPISBAIJSetPreallocationCSR(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscScalar})::PetscErrorCode
end

function MatXAIJSetPreallocation(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatXAIJSetPreallocation(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscInt})::PetscErrorCode
end

function MatCreateShell(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.MatCreateShell(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{Cvoid}, arg7::Ptr{Mat})::PetscErrorCode
end

function MatCreateNormal(arg1, arg2)
    @ccall libpetsc.MatCreateNormal(arg1::Mat, arg2::Ptr{Mat})::PetscErrorCode
end

function MatCreateNormalHermitian(arg1, arg2)
    @ccall libpetsc.MatCreateNormalHermitian(arg1::Mat, arg2::Ptr{Mat})::PetscErrorCode
end

function MatCreateLRC(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatCreateLRC(arg1::Mat, arg2::Mat, arg3::Vec, arg4::Mat, arg5::Ptr{Mat})::PetscErrorCode
end

function MatLRCGetMats(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatLRCGetMats(arg1::Mat, arg2::Ptr{Mat}, arg3::Ptr{Mat}, arg4::Ptr{Vec}, arg5::Ptr{Mat})::PetscErrorCode
end

function MatCreateIS(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.MatCreateIS(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::ISLocalToGlobalMapping, arg8::ISLocalToGlobalMapping, arg9::Ptr{Mat})::PetscErrorCode
end

function MatCreateSeqAIJCRL(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatCreateSeqAIJCRL(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{Mat})::PetscErrorCode
end

function MatCreateMPIAIJCRL(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.MatCreateMPIAIJCRL(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::PetscInt, arg7::Ptr{PetscInt}, arg8::Ptr{Mat})::PetscErrorCode
end

function MatCreateScatter(arg1, arg2, arg3)
    @ccall libpetsc.MatCreateScatter(arg1::MPI_Comm, arg2::VecScatter, arg3::Ptr{Mat})::PetscErrorCode
end

function MatScatterSetVecScatter(arg1, arg2)
    @ccall libpetsc.MatScatterSetVecScatter(arg1::Mat, arg2::VecScatter)::PetscErrorCode
end

function MatScatterGetVecScatter(arg1, arg2)
    @ccall libpetsc.MatScatterGetVecScatter(arg1::Mat, arg2::Ptr{VecScatter})::PetscErrorCode
end

function MatCreateBlockMat(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.MatCreateBlockMat(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::Ptr{Mat})::PetscErrorCode
end

function MatCompositeAddMat(arg1, arg2)
    @ccall libpetsc.MatCompositeAddMat(arg1::Mat, arg2::Mat)::PetscErrorCode
end

function MatCompositeMerge(arg1)
    @ccall libpetsc.MatCompositeMerge(arg1::Mat)::PetscErrorCode
end

@enum MatCompositeMergeType::UInt32 begin
    MAT_COMPOSITE_MERGE_RIGHT = 0
    MAT_COMPOSITE_MERGE_LEFT = 1
end

function MatCompositeSetMergeType(arg1, arg2)
    @ccall libpetsc.MatCompositeSetMergeType(arg1::Mat, arg2::MatCompositeMergeType)::PetscErrorCode
end

function MatCreateComposite(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatCreateComposite(arg1::MPI_Comm, arg2::PetscInt, arg3::Ptr{Mat}, arg4::Ptr{Mat})::PetscErrorCode
end

@enum MatCompositeType::UInt32 begin
    MAT_COMPOSITE_ADDITIVE = 0
    MAT_COMPOSITE_MULTIPLICATIVE = 1
end

function MatCompositeSetType(arg1, arg2)
    @ccall libpetsc.MatCompositeSetType(arg1::Mat, arg2::MatCompositeType)::PetscErrorCode
end

function MatCompositeGetType(arg1, arg2)
    @ccall libpetsc.MatCompositeGetType(arg1::Mat, arg2::Ptr{MatCompositeType})::PetscErrorCode
end

function MatCompositeSetMatStructure(arg1, arg2)
    @ccall libpetsc.MatCompositeSetMatStructure(arg1::Mat, arg2::MatStructure)::PetscErrorCode
end

function MatCompositeGetMatStructure(arg1, arg2)
    @ccall libpetsc.MatCompositeGetMatStructure(arg1::Mat, arg2::Ptr{MatStructure})::PetscErrorCode
end

function MatCompositeGetNumberMat(arg1, arg2)
    @ccall libpetsc.MatCompositeGetNumberMat(arg1::Mat, arg2::Ptr{PetscInt})::PetscErrorCode
end

function MatCompositeGetMat(arg1, arg2, arg3)
    @ccall libpetsc.MatCompositeGetMat(arg1::Mat, arg2::PetscInt, arg3::Ptr{Mat})::PetscErrorCode
end

function MatCompositeSetScalings(arg1, arg2)
    @ccall libpetsc.MatCompositeSetScalings(arg1::Mat, arg2::Ptr{PetscScalar})::PetscErrorCode
end

function MatCreateFFT(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatCreateFFT(arg1::MPI_Comm, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::MatType, arg5::Ptr{Mat})::PetscErrorCode
end

function MatCreateSeqCUFFT(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatCreateSeqCUFFT(arg1::MPI_Comm, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{Mat})::PetscErrorCode
end

function MatCreateTranspose(arg1, arg2)
    @ccall libpetsc.MatCreateTranspose(arg1::Mat, arg2::Ptr{Mat})::PetscErrorCode
end

function MatTransposeGetMat(arg1, arg2)
    @ccall libpetsc.MatTransposeGetMat(arg1::Mat, arg2::Ptr{Mat})::PetscErrorCode
end

function MatCreateHermitianTranspose(arg1, arg2)
    @ccall libpetsc.MatCreateHermitianTranspose(arg1::Mat, arg2::Ptr{Mat})::PetscErrorCode
end

function MatHermitianTransposeGetMat(arg1, arg2)
    @ccall libpetsc.MatHermitianTransposeGetMat(arg1::Mat, arg2::Ptr{Mat})::PetscErrorCode
end

function MatCreateSubMatrixVirtual(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatCreateSubMatrixVirtual(arg1::Mat, arg2::IS, arg3::IS, arg4::Ptr{Mat})::PetscErrorCode
end

function MatSubMatrixVirtualUpdate(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatSubMatrixVirtualUpdate(arg1::Mat, arg2::Mat, arg3::IS, arg4::IS)::PetscErrorCode
end

function MatCreateLocalRef(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatCreateLocalRef(arg1::Mat, arg2::IS, arg3::IS, arg4::Ptr{Mat})::PetscErrorCode
end

function MatCreateConstantDiagonal(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.MatCreateConstantDiagonal(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscScalar, arg7::Ptr{Mat})::PetscErrorCode
end

function MatPythonSetType(arg1, arg2)
    @ccall libpetsc.MatPythonSetType(arg1::Mat, arg2::Ptr{Cchar})::PetscErrorCode
end

function MatResetPreallocation(arg1)
    @ccall libpetsc.MatResetPreallocation(arg1::Mat)::PetscErrorCode
end

function MatSetUp(arg1)
    @ccall libpetsc.MatSetUp(arg1::Mat)::PetscErrorCode
end

function MatDestroy(arg1)
    @ccall libpetsc.MatDestroy(arg1::Ptr{Mat})::PetscErrorCode
end

function MatGetNonzeroState(arg1, arg2)
    @ccall libpetsc.MatGetNonzeroState(arg1::Mat, arg2::Ptr{PetscObjectState})::PetscErrorCode
end

function MatConjugate(arg1)
    @ccall libpetsc.MatConjugate(arg1::Mat)::PetscErrorCode
end

function MatRealPart(arg1)
    @ccall libpetsc.MatRealPart(arg1::Mat)::PetscErrorCode
end

function MatImaginaryPart(arg1)
    @ccall libpetsc.MatImaginaryPart(arg1::Mat)::PetscErrorCode
end

function MatGetDiagonalBlock(arg1, arg2)
    @ccall libpetsc.MatGetDiagonalBlock(arg1::Mat, arg2::Ptr{Mat})::PetscErrorCode
end

function MatGetTrace(arg1, arg2)
    @ccall libpetsc.MatGetTrace(arg1::Mat, arg2::Ptr{PetscScalar})::PetscErrorCode
end

function MatInvertBlockDiagonal(arg1, arg2)
    @ccall libpetsc.MatInvertBlockDiagonal(arg1::Mat, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatInvertVariableBlockDiagonal(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatInvertVariableBlockDiagonal(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscScalar})::PetscErrorCode
end

function MatInvertBlockDiagonalMat(arg1, arg2)
    @ccall libpetsc.MatInvertBlockDiagonalMat(arg1::Mat, arg2::Mat)::PetscErrorCode
end

function MatSetValues(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.MatSetValues(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{PetscScalar}, arg7::InsertMode)::PetscErrorCode
end

function MatSetValuesBlocked(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.MatSetValuesBlocked(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{PetscScalar}, arg7::InsertMode)::PetscErrorCode
end

function MatSetValuesRow(arg1, arg2, arg3)
    @ccall libpetsc.MatSetValuesRow(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscScalar})::PetscErrorCode
end

function MatSetValuesRowLocal(arg1, arg2, arg3)
    @ccall libpetsc.MatSetValuesRowLocal(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscScalar})::PetscErrorCode
end

function MatSetValuesBatch(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatSetValuesBatch(arg1::Mat, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::Ptr{PetscScalar})::PetscErrorCode
end

function MatSetRandom(arg1, arg2)
    @ccall libpetsc.MatSetRandom(arg1::Mat, arg2::PetscRandom)::PetscErrorCode
end

mutable struct MatStencil
    k::PetscInt
    j::PetscInt
    i::PetscInt
    c::PetscInt
    MatStencil() = new()
end

function MatSetValuesStencil(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.MatSetValuesStencil(arg1::Mat, arg2::PetscInt, arg3::Ptr{MatStencil}, arg4::PetscInt, arg5::Ptr{MatStencil}, arg6::Ptr{PetscScalar}, arg7::InsertMode)::PetscErrorCode
end

function MatSetValuesBlockedStencil(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.MatSetValuesBlockedStencil(arg1::Mat, arg2::PetscInt, arg3::Ptr{MatStencil}, arg4::PetscInt, arg5::Ptr{MatStencil}, arg6::Ptr{PetscScalar}, arg7::InsertMode)::PetscErrorCode
end

function MatSetStencil(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatSetStencil(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::PetscInt)::PetscErrorCode
end

@enum MatAssemblyType::UInt32 begin
    MAT_FLUSH_ASSEMBLY = 1
    MAT_FINAL_ASSEMBLY = 0
end

function MatAssemblyBegin(arg1, arg2)
    @ccall libpetsc.MatAssemblyBegin(arg1::Mat, arg2::MatAssemblyType)::PetscErrorCode
end

function MatAssemblyEnd(arg1, arg2)
    @ccall libpetsc.MatAssemblyEnd(arg1::Mat, arg2::MatAssemblyType)::PetscErrorCode
end

function MatAssembled(arg1, arg2)
    @ccall libpetsc.MatAssembled(arg1::Mat, arg2::Ptr{PetscBool})::PetscErrorCode
end

@enum MatOption::Int32 begin
    MAT_OPTION_MIN = -3
    MAT_UNUSED_NONZERO_LOCATION_ERR = -2
    MAT_ROW_ORIENTED = -1
    MAT_SYMMETRIC = 1
    MAT_STRUCTURALLY_SYMMETRIC = 2
    MAT_FORCE_DIAGONAL_ENTRIES = 3
    MAT_IGNORE_OFF_PROC_ENTRIES = 4
    MAT_USE_HASH_TABLE = 5
    MAT_KEEP_NONZERO_PATTERN = 6
    MAT_IGNORE_ZERO_ENTRIES = 7
    MAT_USE_INODES = 8
    MAT_HERMITIAN = 9
    MAT_SYMMETRY_ETERNAL = 10
    MAT_NEW_NONZERO_LOCATION_ERR = 11
    MAT_IGNORE_LOWER_TRIANGULAR = 12
    MAT_ERROR_LOWER_TRIANGULAR = 13
    MAT_GETROW_UPPERTRIANGULAR = 14
    MAT_SPD = 15
    MAT_NO_OFF_PROC_ZERO_ROWS = 16
    MAT_NO_OFF_PROC_ENTRIES = 17
    MAT_NEW_NONZERO_LOCATIONS = 18
    MAT_NEW_NONZERO_ALLOCATION_ERR = 19
    MAT_SUBSET_OFF_PROC_ENTRIES = 20
    MAT_SUBMAT_SINGLEIS = 21
    MAT_STRUCTURE_ONLY = 22
    MAT_SORTED_FULL = 23
    MAT_FORM_EXPLICIT_TRANSPOSE = 24
    MAT_OPTION_MAX = 25
end

function MatSetOption(arg1, arg2, arg3)
    @ccall libpetsc.MatSetOption(arg1::Mat, arg2::MatOption, arg3::PetscBool)::PetscErrorCode
end

function MatGetOption(arg1, arg2, arg3)
    @ccall libpetsc.MatGetOption(arg1::Mat, arg2::MatOption, arg3::Ptr{PetscBool})::PetscErrorCode
end

function MatPropagateSymmetryOptions(arg1, arg2)
    @ccall libpetsc.MatPropagateSymmetryOptions(arg1::Mat, arg2::Mat)::PetscErrorCode
end

function MatGetType(arg1, arg2)
    @ccall libpetsc.MatGetType(arg1::Mat, arg2::Ptr{MatType})::PetscErrorCode
end

function MatGetValues(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatGetValues(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{PetscScalar})::PetscErrorCode
end

function MatGetRow(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatGetRow(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{PetscInt}}, arg5::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatRestoreRow(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatRestoreRow(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{PetscInt}}, arg5::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatGetRowUpperTriangular(arg1)
    @ccall libpetsc.MatGetRowUpperTriangular(arg1::Mat)::PetscErrorCode
end

function MatRestoreRowUpperTriangular(arg1)
    @ccall libpetsc.MatRestoreRowUpperTriangular(arg1::Mat)::PetscErrorCode
end

function MatGetColumnVector(arg1, arg2, arg3)
    @ccall libpetsc.MatGetColumnVector(arg1::Mat, arg2::Vec, arg3::PetscInt)::PetscErrorCode
end

function MatSeqAIJGetArray(arg1, arg2)
    @ccall libpetsc.MatSeqAIJGetArray(arg1::Mat, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatSeqAIJGetArrayRead(arg1, arg2)
    @ccall libpetsc.MatSeqAIJGetArrayRead(arg1::Mat, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatSeqAIJRestoreArray(arg1, arg2)
    @ccall libpetsc.MatSeqAIJRestoreArray(arg1::Mat, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatSeqAIJRestoreArrayRead(arg1, arg2)
    @ccall libpetsc.MatSeqAIJRestoreArrayRead(arg1::Mat, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatSeqAIJGetMaxRowNonzeros(arg1, arg2)
    @ccall libpetsc.MatSeqAIJGetMaxRowNonzeros(arg1::Mat, arg2::Ptr{PetscInt})::PetscErrorCode
end

function MatSeqAIJSetValuesLocalFast(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.MatSeqAIJSetValuesLocalFast(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{PetscScalar}, arg7::InsertMode)::PetscErrorCode
end

function MatSeqAIJSetType(arg1, arg2)
    @ccall libpetsc.MatSeqAIJSetType(arg1::Mat, arg2::MatType)::PetscErrorCode
end

function MatSeqAIJRegister(arg1, arg2)
    @ccall libpetsc.MatSeqAIJRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function MatSeqBAIJGetArray(arg1, arg2)
    @ccall libpetsc.MatSeqBAIJGetArray(arg1::Mat, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatSeqBAIJRestoreArray(arg1, arg2)
    @ccall libpetsc.MatSeqBAIJRestoreArray(arg1::Mat, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatSeqSBAIJGetArray(arg1, arg2)
    @ccall libpetsc.MatSeqSBAIJGetArray(arg1::Mat, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatSeqSBAIJRestoreArray(arg1, arg2)
    @ccall libpetsc.MatSeqSBAIJRestoreArray(arg1::Mat, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatDenseGetArray(arg1, arg2)
    @ccall libpetsc.MatDenseGetArray(arg1::Mat, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatDenseRestoreArray(arg1, arg2)
    @ccall libpetsc.MatDenseRestoreArray(arg1::Mat, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatDensePlaceArray(arg1, arg2)
    @ccall libpetsc.MatDensePlaceArray(arg1::Mat, arg2::Ptr{PetscScalar})::PetscErrorCode
end

function MatDenseReplaceArray(arg1, arg2)
    @ccall libpetsc.MatDenseReplaceArray(arg1::Mat, arg2::Ptr{PetscScalar})::PetscErrorCode
end

function MatDenseResetArray(arg1)
    @ccall libpetsc.MatDenseResetArray(arg1::Mat)::PetscErrorCode
end

function MatDenseGetArrayRead(arg1, arg2)
    @ccall libpetsc.MatDenseGetArrayRead(arg1::Mat, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatDenseRestoreArrayRead(arg1, arg2)
    @ccall libpetsc.MatDenseRestoreArrayRead(arg1::Mat, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatDenseGetArrayWrite(arg1, arg2)
    @ccall libpetsc.MatDenseGetArrayWrite(arg1::Mat, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatDenseRestoreArrayWrite(arg1, arg2)
    @ccall libpetsc.MatDenseRestoreArrayWrite(arg1::Mat, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatGetBlockSize(arg1, arg2)
    @ccall libpetsc.MatGetBlockSize(arg1::Mat, arg2::Ptr{PetscInt})::PetscErrorCode
end

function MatSetBlockSize(arg1, arg2)
    @ccall libpetsc.MatSetBlockSize(arg1::Mat, arg2::PetscInt)::PetscErrorCode
end

function MatGetBlockSizes(arg1, arg2, arg3)
    @ccall libpetsc.MatGetBlockSizes(arg1::Mat, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function MatSetBlockSizes(arg1, arg2, arg3)
    @ccall libpetsc.MatSetBlockSizes(arg1::Mat, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function MatSetBlockSizesFromMats(arg1, arg2, arg3)
    @ccall libpetsc.MatSetBlockSizesFromMats(arg1::Mat, arg2::Mat, arg3::Mat)::PetscErrorCode
end

function MatSetVariableBlockSizes(arg1, arg2, arg3)
    @ccall libpetsc.MatSetVariableBlockSizes(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function MatGetVariableBlockSizes(arg1, arg2, arg3)
    @ccall libpetsc.MatGetVariableBlockSizes(arg1::Mat, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function MatDenseGetColumn(arg1, arg2, arg3)
    @ccall libpetsc.MatDenseGetColumn(arg1::Mat, arg2::PetscInt, arg3::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatDenseRestoreColumn(arg1, arg2)
    @ccall libpetsc.MatDenseRestoreColumn(arg1::Mat, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatDenseGetColumnVec(arg1, arg2, arg3)
    @ccall libpetsc.MatDenseGetColumnVec(arg1::Mat, arg2::PetscInt, arg3::Ptr{Vec})::PetscErrorCode
end

function MatDenseRestoreColumnVec(arg1, arg2, arg3)
    @ccall libpetsc.MatDenseRestoreColumnVec(arg1::Mat, arg2::PetscInt, arg3::Ptr{Vec})::PetscErrorCode
end

function MatDenseGetColumnVecRead(arg1, arg2, arg3)
    @ccall libpetsc.MatDenseGetColumnVecRead(arg1::Mat, arg2::PetscInt, arg3::Ptr{Vec})::PetscErrorCode
end

function MatDenseRestoreColumnVecRead(arg1, arg2, arg3)
    @ccall libpetsc.MatDenseRestoreColumnVecRead(arg1::Mat, arg2::PetscInt, arg3::Ptr{Vec})::PetscErrorCode
end

function MatDenseGetColumnVecWrite(arg1, arg2, arg3)
    @ccall libpetsc.MatDenseGetColumnVecWrite(arg1::Mat, arg2::PetscInt, arg3::Ptr{Vec})::PetscErrorCode
end

function MatDenseRestoreColumnVecWrite(arg1, arg2, arg3)
    @ccall libpetsc.MatDenseRestoreColumnVecWrite(arg1::Mat, arg2::PetscInt, arg3::Ptr{Vec})::PetscErrorCode
end

function MatDenseGetSubMatrix(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatDenseGetSubMatrix(arg1::Mat, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Mat})::PetscErrorCode
end

function MatDenseRestoreSubMatrix(arg1, arg2)
    @ccall libpetsc.MatDenseRestoreSubMatrix(arg1::Mat, arg2::Ptr{Mat})::PetscErrorCode
end

function MatMult(arg1, arg2, arg3)
    @ccall libpetsc.MatMult(arg1::Mat, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function MatMultDiagonalBlock(arg1, arg2, arg3)
    @ccall libpetsc.MatMultDiagonalBlock(arg1::Mat, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function MatMultAdd(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatMultAdd(arg1::Mat, arg2::Vec, arg3::Vec, arg4::Vec)::PetscErrorCode
end

function MatMultTranspose(arg1, arg2, arg3)
    @ccall libpetsc.MatMultTranspose(arg1::Mat, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function MatMultHermitianTranspose(arg1, arg2, arg3)
    @ccall libpetsc.MatMultHermitianTranspose(arg1::Mat, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function MatIsTranspose(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatIsTranspose(arg1::Mat, arg2::Mat, arg3::PetscReal, arg4::Ptr{PetscBool})::PetscErrorCode
end

function MatIsHermitianTranspose(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatIsHermitianTranspose(arg1::Mat, arg2::Mat, arg3::PetscReal, arg4::Ptr{PetscBool})::PetscErrorCode
end

function MatMultTransposeAdd(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatMultTransposeAdd(arg1::Mat, arg2::Vec, arg3::Vec, arg4::Vec)::PetscErrorCode
end

function MatMultHermitianTransposeAdd(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatMultHermitianTransposeAdd(arg1::Mat, arg2::Vec, arg3::Vec, arg4::Vec)::PetscErrorCode
end

function MatMultConstrained(arg1, arg2, arg3)
    @ccall libpetsc.MatMultConstrained(arg1::Mat, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function MatMultTransposeConstrained(arg1, arg2, arg3)
    @ccall libpetsc.MatMultTransposeConstrained(arg1::Mat, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function MatMatSolve(arg1, arg2, arg3)
    @ccall libpetsc.MatMatSolve(arg1::Mat, arg2::Mat, arg3::Mat)::PetscErrorCode
end

function MatMatSolveTranspose(arg1, arg2, arg3)
    @ccall libpetsc.MatMatSolveTranspose(arg1::Mat, arg2::Mat, arg3::Mat)::PetscErrorCode
end

function MatMatTransposeSolve(arg1, arg2, arg3)
    @ccall libpetsc.MatMatTransposeSolve(arg1::Mat, arg2::Mat, arg3::Mat)::PetscErrorCode
end

function MatResidual(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatResidual(arg1::Mat, arg2::Vec, arg3::Vec, arg4::Vec)::PetscErrorCode
end

@enum MatDuplicateOption::UInt32 begin
    MAT_DO_NOT_COPY_VALUES = 0
    MAT_COPY_VALUES = 1
    MAT_SHARE_NONZERO_PATTERN = 2
end

function MatConvert(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatConvert(arg1::Mat, arg2::MatType, arg3::MatReuse, arg4::Ptr{Mat})::PetscErrorCode
end

function MatDuplicate(arg1, arg2, arg3)
    @ccall libpetsc.MatDuplicate(arg1::Mat, arg2::MatDuplicateOption, arg3::Ptr{Mat})::PetscErrorCode
end

function MatCopy(arg1, arg2, arg3)
    @ccall libpetsc.MatCopy(arg1::Mat, arg2::Mat, arg3::MatStructure)::PetscErrorCode
end

function MatView(arg1, arg2)
    @ccall libpetsc.MatView(arg1::Mat, arg2::PetscViewer)::PetscErrorCode
end

function MatIsSymmetric(arg1, arg2, arg3)
    @ccall libpetsc.MatIsSymmetric(arg1::Mat, arg2::PetscReal, arg3::Ptr{PetscBool})::PetscErrorCode
end

function MatIsStructurallySymmetric(arg1, arg2)
    @ccall libpetsc.MatIsStructurallySymmetric(arg1::Mat, arg2::Ptr{PetscBool})::PetscErrorCode
end

function MatIsHermitian(arg1, arg2, arg3)
    @ccall libpetsc.MatIsHermitian(arg1::Mat, arg2::PetscReal, arg3::Ptr{PetscBool})::PetscErrorCode
end

function MatIsSymmetricKnown(arg1, arg2, arg3)
    @ccall libpetsc.MatIsSymmetricKnown(arg1::Mat, arg2::Ptr{PetscBool}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function MatIsHermitianKnown(arg1, arg2, arg3)
    @ccall libpetsc.MatIsHermitianKnown(arg1::Mat, arg2::Ptr{PetscBool}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function MatMissingDiagonal(arg1, arg2, arg3)
    @ccall libpetsc.MatMissingDiagonal(arg1::Mat, arg2::Ptr{PetscBool}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function MatLoad(arg1, arg2)
    @ccall libpetsc.MatLoad(arg1::Mat, arg2::PetscViewer)::PetscErrorCode
end

function MatGetRowIJ(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.MatGetRowIJ(arg1::Mat, arg2::PetscInt, arg3::PetscBool, arg4::PetscBool, arg5::Ptr{PetscInt}, arg6::Ptr{Ptr{PetscInt}}, arg7::Ptr{Ptr{PetscInt}}, arg8::Ptr{PetscBool})::PetscErrorCode
end

function MatRestoreRowIJ(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.MatRestoreRowIJ(arg1::Mat, arg2::PetscInt, arg3::PetscBool, arg4::PetscBool, arg5::Ptr{PetscInt}, arg6::Ptr{Ptr{PetscInt}}, arg7::Ptr{Ptr{PetscInt}}, arg8::Ptr{PetscBool})::PetscErrorCode
end

function MatGetColumnIJ(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.MatGetColumnIJ(arg1::Mat, arg2::PetscInt, arg3::PetscBool, arg4::PetscBool, arg5::Ptr{PetscInt}, arg6::Ptr{Ptr{PetscInt}}, arg7::Ptr{Ptr{PetscInt}}, arg8::Ptr{PetscBool})::PetscErrorCode
end

function MatRestoreColumnIJ(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.MatRestoreColumnIJ(arg1::Mat, arg2::PetscInt, arg3::PetscBool, arg4::PetscBool, arg5::Ptr{PetscInt}, arg6::Ptr{Ptr{PetscInt}}, arg7::Ptr{Ptr{PetscInt}}, arg8::Ptr{PetscBool})::PetscErrorCode
end

mutable struct MatInfo
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
    MatInfo() = new()
end

@enum MatInfoType::UInt32 begin
    MAT_LOCAL = 1
    MAT_GLOBAL_MAX = 2
    MAT_GLOBAL_SUM = 3
end

function MatGetInfo(arg1, arg2, arg3)
    @ccall libpetsc.MatGetInfo(arg1::Mat, arg2::MatInfoType, arg3::Ptr{MatInfo})::PetscErrorCode
end

function MatGetDiagonal(arg1, arg2)
    @ccall libpetsc.MatGetDiagonal(arg1::Mat, arg2::Vec)::PetscErrorCode
end

function MatGetRowMax(arg1, arg2, arg3)
    @ccall libpetsc.MatGetRowMax(arg1::Mat, arg2::Vec, arg3::Ptr{PetscInt})::PetscErrorCode
end

function MatGetRowMin(arg1, arg2, arg3)
    @ccall libpetsc.MatGetRowMin(arg1::Mat, arg2::Vec, arg3::Ptr{PetscInt})::PetscErrorCode
end

function MatGetRowMaxAbs(arg1, arg2, arg3)
    @ccall libpetsc.MatGetRowMaxAbs(arg1::Mat, arg2::Vec, arg3::Ptr{PetscInt})::PetscErrorCode
end

function MatGetRowMinAbs(arg1, arg2, arg3)
    @ccall libpetsc.MatGetRowMinAbs(arg1::Mat, arg2::Vec, arg3::Ptr{PetscInt})::PetscErrorCode
end

function MatGetRowSum(arg1, arg2)
    @ccall libpetsc.MatGetRowSum(arg1::Mat, arg2::Vec)::PetscErrorCode
end

function MatTranspose(arg1, arg2, arg3)
    @ccall libpetsc.MatTranspose(arg1::Mat, arg2::MatReuse, arg3::Ptr{Mat})::PetscErrorCode
end

function MatHermitianTranspose(arg1, arg2, arg3)
    @ccall libpetsc.MatHermitianTranspose(arg1::Mat, arg2::MatReuse, arg3::Ptr{Mat})::PetscErrorCode
end

function MatPermute(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatPermute(arg1::Mat, arg2::IS, arg3::IS, arg4::Ptr{Mat})::PetscErrorCode
end

function MatDiagonalScale(arg1, arg2, arg3)
    @ccall libpetsc.MatDiagonalScale(arg1::Mat, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function MatDiagonalSet(arg1, arg2, arg3)
    @ccall libpetsc.MatDiagonalSet(arg1::Mat, arg2::Vec, arg3::InsertMode)::PetscErrorCode
end

function MatEqual(arg1, arg2, arg3)
    @ccall libpetsc.MatEqual(arg1::Mat, arg2::Mat, arg3::Ptr{PetscBool})::PetscErrorCode
end

function MatMultEqual(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatMultEqual(arg1::Mat, arg2::Mat, arg3::PetscInt, arg4::Ptr{PetscBool})::PetscErrorCode
end

function MatMultAddEqual(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatMultAddEqual(arg1::Mat, arg2::Mat, arg3::PetscInt, arg4::Ptr{PetscBool})::PetscErrorCode
end

function MatMultTransposeEqual(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatMultTransposeEqual(arg1::Mat, arg2::Mat, arg3::PetscInt, arg4::Ptr{PetscBool})::PetscErrorCode
end

function MatMultTransposeAddEqual(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatMultTransposeAddEqual(arg1::Mat, arg2::Mat, arg3::PetscInt, arg4::Ptr{PetscBool})::PetscErrorCode
end

function MatMatMultEqual(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatMatMultEqual(arg1::Mat, arg2::Mat, arg3::Mat, arg4::PetscInt, arg5::Ptr{PetscBool})::PetscErrorCode
end

function MatTransposeMatMultEqual(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatTransposeMatMultEqual(arg1::Mat, arg2::Mat, arg3::Mat, arg4::PetscInt, arg5::Ptr{PetscBool})::PetscErrorCode
end

function MatMatTransposeMultEqual(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatMatTransposeMultEqual(arg1::Mat, arg2::Mat, arg3::Mat, arg4::PetscInt, arg5::Ptr{PetscBool})::PetscErrorCode
end

function MatPtAPMultEqual(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatPtAPMultEqual(arg1::Mat, arg2::Mat, arg3::Mat, arg4::PetscInt, arg5::Ptr{PetscBool})::PetscErrorCode
end

function MatRARtMultEqual(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatRARtMultEqual(arg1::Mat, arg2::Mat, arg3::Mat, arg4::PetscInt, arg5::Ptr{PetscBool})::PetscErrorCode
end

function MatIsLinear(arg1, arg2, arg3)
    @ccall libpetsc.MatIsLinear(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscBool})::PetscErrorCode
end

function MatNorm(arg1, arg2, arg3)
    @ccall libpetsc.MatNorm(arg1::Mat, arg2::NormType, arg3::Ptr{PetscReal})::PetscErrorCode
end

function MatGetColumnNorms(arg1, arg2, arg3)
    @ccall libpetsc.MatGetColumnNorms(arg1::Mat, arg2::NormType, arg3::Ptr{PetscReal})::PetscErrorCode
end

function MatZeroEntries(arg1)
    @ccall libpetsc.MatZeroEntries(arg1::Mat)::PetscErrorCode
end

function MatSetInf(arg1)
    @ccall libpetsc.MatSetInf(arg1::Mat)::PetscErrorCode
end

function MatZeroRows(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatZeroRows(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::PetscScalar, arg5::Vec, arg6::Vec)::PetscErrorCode
end

function MatZeroRowsIS(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatZeroRowsIS(arg1::Mat, arg2::IS, arg3::PetscScalar, arg4::Vec, arg5::Vec)::PetscErrorCode
end

function MatZeroRowsStencil(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatZeroRowsStencil(arg1::Mat, arg2::PetscInt, arg3::Ptr{MatStencil}, arg4::PetscScalar, arg5::Vec, arg6::Vec)::PetscErrorCode
end

function MatZeroRowsColumnsStencil(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatZeroRowsColumnsStencil(arg1::Mat, arg2::PetscInt, arg3::Ptr{MatStencil}, arg4::PetscScalar, arg5::Vec, arg6::Vec)::PetscErrorCode
end

function MatZeroRowsColumns(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatZeroRowsColumns(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::PetscScalar, arg5::Vec, arg6::Vec)::PetscErrorCode
end

function MatZeroRowsColumnsIS(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatZeroRowsColumnsIS(arg1::Mat, arg2::IS, arg3::PetscScalar, arg4::Vec, arg5::Vec)::PetscErrorCode
end

function MatGetSize(arg1, arg2, arg3)
    @ccall libpetsc.MatGetSize(arg1::Mat, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function MatGetLocalSize(arg1, arg2, arg3)
    @ccall libpetsc.MatGetLocalSize(arg1::Mat, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function MatGetOwnershipRange(arg1, arg2, arg3)
    @ccall libpetsc.MatGetOwnershipRange(arg1::Mat, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function MatGetOwnershipRanges(arg1, arg2)
    @ccall libpetsc.MatGetOwnershipRanges(arg1::Mat, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function MatGetOwnershipRangeColumn(arg1, arg2, arg3)
    @ccall libpetsc.MatGetOwnershipRangeColumn(arg1::Mat, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function MatGetOwnershipRangesColumn(arg1, arg2)
    @ccall libpetsc.MatGetOwnershipRangesColumn(arg1::Mat, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function MatGetOwnershipIS(arg1, arg2, arg3)
    @ccall libpetsc.MatGetOwnershipIS(arg1::Mat, arg2::Ptr{IS}, arg3::Ptr{IS})::PetscErrorCode
end

function MatCreateSubMatrices(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatCreateSubMatrices(arg1::Mat, arg2::PetscInt, arg3::Ptr{IS}, arg4::Ptr{IS}, arg5::MatReuse, arg6::Ptr{Ptr{Mat}})::PetscErrorCode
end

function MatGetSubMatrices(mat, n, irow, icol, scall, submat)
    @ccall libpetsc.MatGetSubMatrices(mat::Mat, n::PetscInt, irow::Ptr{IS}, icol::Ptr{IS}, scall::MatReuse, submat::Ptr{Ptr{Mat}})::PetscErrorCode
end

function MatCreateSubMatricesMPI(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatCreateSubMatricesMPI(arg1::Mat, arg2::PetscInt, arg3::Ptr{IS}, arg4::Ptr{IS}, arg5::MatReuse, arg6::Ptr{Ptr{Mat}})::PetscErrorCode
end

function MatGetSubMatricesMPI(mat, n, irow, icol, scall, submat)
    @ccall libpetsc.MatGetSubMatricesMPI(mat::Mat, n::PetscInt, irow::Ptr{IS}, icol::Ptr{IS}, scall::MatReuse, submat::Ptr{Ptr{Mat}})::PetscErrorCode
end

function MatDestroyMatrices(arg1, arg2)
    @ccall libpetsc.MatDestroyMatrices(arg1::PetscInt, arg2::Ptr{Ptr{Mat}})::PetscErrorCode
end

function MatDestroySubMatrices(arg1, arg2)
    @ccall libpetsc.MatDestroySubMatrices(arg1::PetscInt, arg2::Ptr{Ptr{Mat}})::PetscErrorCode
end

function MatCreateSubMatrix(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatCreateSubMatrix(arg1::Mat, arg2::IS, arg3::IS, arg4::MatReuse, arg5::Ptr{Mat})::PetscErrorCode
end

function MatGetSubMatrix(mat, isrow, iscol, cll, newmat)
    @ccall libpetsc.MatGetSubMatrix(mat::Mat, isrow::IS, iscol::IS, cll::MatReuse, newmat::Ptr{Mat})::PetscErrorCode
end

function MatGetLocalSubMatrix(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatGetLocalSubMatrix(arg1::Mat, arg2::IS, arg3::IS, arg4::Ptr{Mat})::PetscErrorCode
end

function MatRestoreLocalSubMatrix(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatRestoreLocalSubMatrix(arg1::Mat, arg2::IS, arg3::IS, arg4::Ptr{Mat})::PetscErrorCode
end

function MatGetSeqNonzeroStructure(arg1, arg2)
    @ccall libpetsc.MatGetSeqNonzeroStructure(arg1::Mat, arg2::Ptr{Mat})::PetscErrorCode
end

function MatDestroySeqNonzeroStructure(arg1)
    @ccall libpetsc.MatDestroySeqNonzeroStructure(arg1::Ptr{Mat})::PetscErrorCode
end

function MatCreateMPIAIJSumSeqAIJ(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatCreateMPIAIJSumSeqAIJ(arg1::MPI_Comm, arg2::Mat, arg3::PetscInt, arg4::PetscInt, arg5::MatReuse, arg6::Ptr{Mat})::PetscErrorCode
end

function MatCreateMPIAIJSumSeqAIJSymbolic(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatCreateMPIAIJSumSeqAIJSymbolic(arg1::MPI_Comm, arg2::Mat, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{Mat})::PetscErrorCode
end

function MatCreateMPIAIJSumSeqAIJNumeric(arg1, arg2)
    @ccall libpetsc.MatCreateMPIAIJSumSeqAIJNumeric(arg1::Mat, arg2::Mat)::PetscErrorCode
end

function MatMPIAIJGetLocalMat(arg1, arg2, arg3)
    @ccall libpetsc.MatMPIAIJGetLocalMat(arg1::Mat, arg2::MatReuse, arg3::Ptr{Mat})::PetscErrorCode
end

function MatMPIAIJGetLocalMatCondensed(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatMPIAIJGetLocalMatCondensed(arg1::Mat, arg2::MatReuse, arg3::Ptr{IS}, arg4::Ptr{IS}, arg5::Ptr{Mat})::PetscErrorCode
end

function MatMPIAIJGetLocalMatMerge(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatMPIAIJGetLocalMatMerge(arg1::Mat, arg2::MatReuse, arg3::Ptr{IS}, arg4::Ptr{Mat})::PetscErrorCode
end

function MatGetBrowsOfAcols(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatGetBrowsOfAcols(arg1::Mat, arg2::Mat, arg3::MatReuse, arg4::Ptr{IS}, arg5::Ptr{IS}, arg6::Ptr{Mat})::PetscErrorCode
end

function MatGetGhosts(arg1, arg2, arg3)
    @ccall libpetsc.MatGetGhosts(arg1::Mat, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function MatIncreaseOverlap(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatIncreaseOverlap(arg1::Mat, arg2::PetscInt, arg3::Ptr{IS}, arg4::PetscInt)::PetscErrorCode
end

function MatIncreaseOverlapSplit(mat, n, is, ov)
    @ccall libpetsc.MatIncreaseOverlapSplit(mat::Mat, n::PetscInt, is::Ptr{IS}, ov::PetscInt)::PetscErrorCode
end

function MatMPIAIJSetUseScalableIncreaseOverlap(arg1, arg2)
    @ccall libpetsc.MatMPIAIJSetUseScalableIncreaseOverlap(arg1::Mat, arg2::PetscBool)::PetscErrorCode
end

function MatMatMult(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatMatMult(arg1::Mat, arg2::Mat, arg3::MatReuse, arg4::PetscReal, arg5::Ptr{Mat})::PetscErrorCode
end

function MatMatMatMult(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatMatMatMult(arg1::Mat, arg2::Mat, arg3::Mat, arg4::MatReuse, arg5::PetscReal, arg6::Ptr{Mat})::PetscErrorCode
end

function MatGalerkin(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatGalerkin(arg1::Mat, arg2::Mat, arg3::Mat, arg4::MatReuse, arg5::PetscReal, arg6::Ptr{Mat})::PetscErrorCode
end

function MatPtAP(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatPtAP(arg1::Mat, arg2::Mat, arg3::MatReuse, arg4::PetscReal, arg5::Ptr{Mat})::PetscErrorCode
end

function MatRARt(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatRARt(arg1::Mat, arg2::Mat, arg3::MatReuse, arg4::PetscReal, arg5::Ptr{Mat})::PetscErrorCode
end

function MatTransposeMatMult(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatTransposeMatMult(arg1::Mat, arg2::Mat, arg3::MatReuse, arg4::PetscReal, arg5::Ptr{Mat})::PetscErrorCode
end

function MatMatTransposeMult(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatMatTransposeMult(arg1::Mat, arg2::Mat, arg3::MatReuse, arg4::PetscReal, arg5::Ptr{Mat})::PetscErrorCode
end

function MatAXPY(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatAXPY(arg1::Mat, arg2::PetscScalar, arg3::Mat, arg4::MatStructure)::PetscErrorCode
end

function MatAYPX(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatAYPX(arg1::Mat, arg2::PetscScalar, arg3::Mat, arg4::MatStructure)::PetscErrorCode
end

function MatScale(arg1, arg2)
    @ccall libpetsc.MatScale(arg1::Mat, arg2::PetscScalar)::PetscErrorCode
end

function MatShift(arg1, arg2)
    @ccall libpetsc.MatShift(arg1::Mat, arg2::PetscScalar)::PetscErrorCode
end

function MatSetLocalToGlobalMapping(arg1, arg2, arg3)
    @ccall libpetsc.MatSetLocalToGlobalMapping(arg1::Mat, arg2::ISLocalToGlobalMapping, arg3::ISLocalToGlobalMapping)::PetscErrorCode
end

function MatGetLocalToGlobalMapping(arg1, arg2, arg3)
    @ccall libpetsc.MatGetLocalToGlobalMapping(arg1::Mat, arg2::Ptr{ISLocalToGlobalMapping}, arg3::Ptr{ISLocalToGlobalMapping})::PetscErrorCode
end

function MatGetLayouts(arg1, arg2, arg3)
    @ccall libpetsc.MatGetLayouts(arg1::Mat, arg2::Ptr{PetscLayout}, arg3::Ptr{PetscLayout})::PetscErrorCode
end

function MatSetLayouts(arg1, arg2, arg3)
    @ccall libpetsc.MatSetLayouts(arg1::Mat, arg2::PetscLayout, arg3::PetscLayout)::PetscErrorCode
end

function MatZeroRowsLocal(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatZeroRowsLocal(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::PetscScalar, arg5::Vec, arg6::Vec)::PetscErrorCode
end

function MatZeroRowsLocalIS(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatZeroRowsLocalIS(arg1::Mat, arg2::IS, arg3::PetscScalar, arg4::Vec, arg5::Vec)::PetscErrorCode
end

function MatZeroRowsColumnsLocal(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatZeroRowsColumnsLocal(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::PetscScalar, arg5::Vec, arg6::Vec)::PetscErrorCode
end

function MatZeroRowsColumnsLocalIS(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatZeroRowsColumnsLocalIS(arg1::Mat, arg2::IS, arg3::PetscScalar, arg4::Vec, arg5::Vec)::PetscErrorCode
end

function MatGetValuesLocal(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatGetValuesLocal(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{PetscScalar})::PetscErrorCode
end

function MatSetValuesLocal(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.MatSetValuesLocal(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{PetscScalar}, arg7::InsertMode)::PetscErrorCode
end

function MatSetValuesBlockedLocal(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.MatSetValuesBlockedLocal(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{PetscScalar}, arg7::InsertMode)::PetscErrorCode
end

function MatStashSetInitialSize(arg1, arg2, arg3)
    @ccall libpetsc.MatStashSetInitialSize(arg1::Mat, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function MatStashGetInfo(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatStashGetInfo(arg1::Mat, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt})::PetscErrorCode
end

function MatInterpolate(arg1, arg2, arg3)
    @ccall libpetsc.MatInterpolate(arg1::Mat, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function MatInterpolateAdd(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatInterpolateAdd(arg1::Mat, arg2::Vec, arg3::Vec, arg4::Vec)::PetscErrorCode
end

function MatRestrict(arg1, arg2, arg3)
    @ccall libpetsc.MatRestrict(arg1::Mat, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function MatMatInterpolate(arg1, arg2, arg3)
    @ccall libpetsc.MatMatInterpolate(arg1::Mat, arg2::Mat, arg3::Ptr{Mat})::PetscErrorCode
end

function MatMatInterpolateAdd(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatMatInterpolateAdd(arg1::Mat, arg2::Mat, arg3::Mat, arg4::Ptr{Mat})::PetscErrorCode
end

function MatMatRestrict(arg1, arg2, arg3)
    @ccall libpetsc.MatMatRestrict(arg1::Mat, arg2::Mat, arg3::Ptr{Mat})::PetscErrorCode
end

function MatCreateVecs(arg1, arg2, arg3)
    @ccall libpetsc.MatCreateVecs(arg1::Mat, arg2::Ptr{Vec}, arg3::Ptr{Vec})::PetscErrorCode
end

function MatGetVecs(mat, x, y)
    @ccall libpetsc.MatGetVecs(mat::Mat, x::Ptr{Vec}, y::Ptr{Vec})::PetscErrorCode
end

function MatCreateRedundantMatrix(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatCreateRedundantMatrix(arg1::Mat, arg2::PetscInt, arg3::MPI_Comm, arg4::MatReuse, arg5::Ptr{Mat})::PetscErrorCode
end

function MatGetMultiProcBlock(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatGetMultiProcBlock(arg1::Mat, arg2::MPI_Comm, arg3::MatReuse, arg4::Ptr{Mat})::PetscErrorCode
end

function MatFindZeroDiagonals(arg1, arg2)
    @ccall libpetsc.MatFindZeroDiagonals(arg1::Mat, arg2::Ptr{IS})::PetscErrorCode
end

function MatFindOffBlockDiagonalEntries(arg1, arg2)
    @ccall libpetsc.MatFindOffBlockDiagonalEntries(arg1::Mat, arg2::Ptr{IS})::PetscErrorCode
end

function MatCreateMPIMatConcatenateSeqMat(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatCreateMPIMatConcatenateSeqMat(arg1::MPI_Comm, arg2::Mat, arg3::PetscInt, arg4::MatReuse, arg5::Ptr{Mat})::PetscErrorCode
end

function MatSetValue(v, i, j, va, mode)
    @ccall libpetsc.MatSetValue(v::Mat, i::PetscInt, j::PetscInt, va::PetscScalar, mode::InsertMode)::PetscErrorCode
end

function MatGetValue(v, i, j, va)
    @ccall libpetsc.MatGetValue(v::Mat, i::PetscInt, j::PetscInt, va::Ptr{PetscScalar})::PetscErrorCode
end

function MatSetValueLocal(v, i, j, va, mode)
    @ccall libpetsc.MatSetValueLocal(v::Mat, i::PetscInt, j::PetscInt, va::PetscScalar, mode::InsertMode)::PetscErrorCode
end

function MatShellGetContext(arg1, arg2)
    @ccall libpetsc.MatShellGetContext(arg1::Mat, arg2::Ptr{Cvoid})::PetscErrorCode
end

function MatInodeAdjustForInodes(arg1, arg2, arg3)
    @ccall libpetsc.MatInodeAdjustForInodes(arg1::Mat, arg2::Ptr{IS}, arg3::Ptr{IS})::PetscErrorCode
end

function MatInodeGetInodeSizes(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatInodeGetInodeSizes(arg1::Mat, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{PetscInt}}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function MatSeqAIJSetColumnIndices(arg1, arg2)
    @ccall libpetsc.MatSeqAIJSetColumnIndices(arg1::Mat, arg2::Ptr{PetscInt})::PetscErrorCode
end

function MatSeqBAIJSetColumnIndices(arg1, arg2)
    @ccall libpetsc.MatSeqBAIJSetColumnIndices(arg1::Mat, arg2::Ptr{PetscInt})::PetscErrorCode
end

function MatCreateSeqAIJWithArrays(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.MatCreateSeqAIJWithArrays(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscScalar}, arg7::Ptr{Mat})::PetscErrorCode
end

function MatCreateSeqBAIJWithArrays(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.MatCreateSeqBAIJWithArrays(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{PetscInt}, arg7::Ptr{PetscScalar}, arg8::Ptr{Mat})::PetscErrorCode
end

function MatCreateSeqSBAIJWithArrays(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.MatCreateSeqSBAIJWithArrays(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{PetscInt}, arg7::Ptr{PetscScalar}, arg8::Ptr{Mat})::PetscErrorCode
end

function MatCreateSeqAIJFromTriple(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.MatCreateSeqAIJFromTriple(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscScalar}, arg7::Ptr{Mat}, arg8::PetscInt, arg9::PetscBool)::PetscErrorCode
end

function MatSeqBAIJSetPreallocation(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatSeqBAIJSetPreallocation(arg1::Mat, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt})::PetscErrorCode
end

function MatSeqSBAIJSetPreallocation(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatSeqSBAIJSetPreallocation(arg1::Mat, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt})::PetscErrorCode
end

function MatSeqAIJSetPreallocation(arg1, arg2, arg3)
    @ccall libpetsc.MatSeqAIJSetPreallocation(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function MatSeqAIJSetTotalPreallocation(arg1, arg2)
    @ccall libpetsc.MatSeqAIJSetTotalPreallocation(arg1::Mat, arg2::PetscInt)::PetscErrorCode
end

function MatMPIBAIJSetPreallocation(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatMPIBAIJSetPreallocation(arg1::Mat, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::PetscInt, arg6::Ptr{PetscInt})::PetscErrorCode
end

function MatMPISBAIJSetPreallocation(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatMPISBAIJSetPreallocation(arg1::Mat, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::PetscInt, arg6::Ptr{PetscInt})::PetscErrorCode
end

function MatMPIAIJSetPreallocation(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatMPIAIJSetPreallocation(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::PetscInt, arg5::Ptr{PetscInt})::PetscErrorCode
end

function MatSeqAIJSetPreallocationCSR(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatSeqAIJSetPreallocationCSR(arg1::Mat, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscScalar})::PetscErrorCode
end

function MatSeqBAIJSetPreallocationCSR(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatSeqBAIJSetPreallocationCSR(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscScalar})::PetscErrorCode
end

function MatMPIAIJSetPreallocationCSR(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatMPIAIJSetPreallocationCSR(arg1::Mat, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscScalar})::PetscErrorCode
end

function MatMPIBAIJSetPreallocationCSR(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatMPIBAIJSetPreallocationCSR(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscScalar})::PetscErrorCode
end

function MatMPIAdjSetPreallocation(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatMPIAdjSetPreallocation(arg1::Mat, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function MatMPIAdjToSeq(arg1, arg2)
    @ccall libpetsc.MatMPIAdjToSeq(arg1::Mat, arg2::Ptr{Mat})::PetscErrorCode
end

function MatMPIDenseSetPreallocation(arg1, arg2)
    @ccall libpetsc.MatMPIDenseSetPreallocation(arg1::Mat, arg2::Ptr{PetscScalar})::PetscErrorCode
end

function MatSeqDenseSetPreallocation(arg1, arg2)
    @ccall libpetsc.MatSeqDenseSetPreallocation(arg1::Mat, arg2::Ptr{PetscScalar})::PetscErrorCode
end

function MatMPIAIJGetSeqAIJ(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatMPIAIJGetSeqAIJ(arg1::Mat, arg2::Ptr{Mat}, arg3::Ptr{Mat}, arg4::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function MatMPIBAIJGetSeqBAIJ(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatMPIBAIJGetSeqBAIJ(arg1::Mat, arg2::Ptr{Mat}, arg3::Ptr{Mat}, arg4::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function MatMPIAdjCreateNonemptySubcommMat(arg1, arg2)
    @ccall libpetsc.MatMPIAdjCreateNonemptySubcommMat(arg1::Mat, arg2::Ptr{Mat})::PetscErrorCode
end

function MatDenseGetLDA(arg1, arg2)
    @ccall libpetsc.MatDenseGetLDA(arg1::Mat, arg2::Ptr{PetscInt})::PetscErrorCode
end

function MatDenseSetLDA(arg1, arg2)
    @ccall libpetsc.MatDenseSetLDA(arg1::Mat, arg2::PetscInt)::PetscErrorCode
end

function MatSeqDenseSetLDA(A, lda)
    @ccall libpetsc.MatSeqDenseSetLDA(A::Mat, lda::PetscInt)::PetscErrorCode
end

function MatDenseGetLocalMatrix(arg1, arg2)
    @ccall libpetsc.MatDenseGetLocalMatrix(arg1::Mat, arg2::Ptr{Mat})::PetscErrorCode
end

function MatBlockMatSetPreallocation(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatBlockMatSetPreallocation(arg1::Mat, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt})::PetscErrorCode
end

function MatStoreValues(arg1)
    @ccall libpetsc.MatStoreValues(arg1::Mat)::PetscErrorCode
end

function MatRetrieveValues(arg1)
    @ccall libpetsc.MatRetrieveValues(arg1::Mat)::PetscErrorCode
end

function MatFindNonzeroRows(arg1, arg2)
    @ccall libpetsc.MatFindNonzeroRows(arg1::Mat, arg2::Ptr{IS})::PetscErrorCode
end

function MatFindZeroRows(arg1, arg2)
    @ccall libpetsc.MatFindZeroRows(arg1::Mat, arg2::Ptr{IS})::PetscErrorCode
end

const MatOrderingType = Ptr{Cchar}

function MatGetOrdering(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatGetOrdering(arg1::Mat, arg2::MatOrderingType, arg3::Ptr{IS}, arg4::Ptr{IS})::PetscErrorCode
end

function MatGetOrderingList(arg1)
    @ccall libpetsc.MatGetOrderingList(arg1::Ptr{PetscFunctionList})::PetscErrorCode
end

function MatOrderingRegister(arg1, arg2)
    @ccall libpetsc.MatOrderingRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function MatReorderForNonzeroDiagonal(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatReorderForNonzeroDiagonal(arg1::Mat, arg2::PetscReal, arg3::IS, arg4::IS)::PetscErrorCode
end

function MatCreateLaplacian(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatCreateLaplacian(arg1::Mat, arg2::PetscReal, arg3::PetscBool, arg4::Ptr{Mat})::PetscErrorCode
end

@enum MatFactorShiftType::UInt32 begin
    MAT_SHIFT_NONE = 0
    MAT_SHIFT_NONZERO = 1
    MAT_SHIFT_POSITIVE_DEFINITE = 2
    MAT_SHIFT_INBLOCKS = 3
end

@enum MatFactorError::UInt32 begin
    MAT_FACTOR_NOERROR = 0
    MAT_FACTOR_STRUCT_ZEROPIVOT = 1
    MAT_FACTOR_NUMERIC_ZEROPIVOT = 2
    MAT_FACTOR_OUTMEMORY = 3
    MAT_FACTOR_OTHER = 4
end

function MatFactorGetError(arg1, arg2)
    @ccall libpetsc.MatFactorGetError(arg1::Mat, arg2::Ptr{MatFactorError})::PetscErrorCode
end

function MatFactorClearError(arg1)
    @ccall libpetsc.MatFactorClearError(arg1::Mat)::PetscErrorCode
end

function MatFactorGetErrorZeroPivot(arg1, arg2, arg3)
    @ccall libpetsc.MatFactorGetErrorZeroPivot(arg1::Mat, arg2::Ptr{PetscReal}, arg3::Ptr{PetscInt})::PetscErrorCode
end

mutable struct MatFactorInfo
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
    MatFactorInfo() = new()
end

function MatFactorInfoInitialize(arg1)
    @ccall libpetsc.MatFactorInfoInitialize(arg1::Ptr{MatFactorInfo})::PetscErrorCode
end

function MatCholeskyFactor(arg1, arg2, arg3)
    @ccall libpetsc.MatCholeskyFactor(arg1::Mat, arg2::IS, arg3::Ptr{MatFactorInfo})::PetscErrorCode
end

function MatCholeskyFactorSymbolic(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatCholeskyFactorSymbolic(arg1::Mat, arg2::Mat, arg3::IS, arg4::Ptr{MatFactorInfo})::PetscErrorCode
end

function MatCholeskyFactorNumeric(arg1, arg2, arg3)
    @ccall libpetsc.MatCholeskyFactorNumeric(arg1::Mat, arg2::Mat, arg3::Ptr{MatFactorInfo})::PetscErrorCode
end

function MatLUFactor(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatLUFactor(arg1::Mat, arg2::IS, arg3::IS, arg4::Ptr{MatFactorInfo})::PetscErrorCode
end

function MatILUFactor(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatILUFactor(arg1::Mat, arg2::IS, arg3::IS, arg4::Ptr{MatFactorInfo})::PetscErrorCode
end

function MatLUFactorSymbolic(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatLUFactorSymbolic(arg1::Mat, arg2::Mat, arg3::IS, arg4::IS, arg5::Ptr{MatFactorInfo})::PetscErrorCode
end

function MatILUFactorSymbolic(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatILUFactorSymbolic(arg1::Mat, arg2::Mat, arg3::IS, arg4::IS, arg5::Ptr{MatFactorInfo})::PetscErrorCode
end

function MatICCFactorSymbolic(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatICCFactorSymbolic(arg1::Mat, arg2::Mat, arg3::IS, arg4::Ptr{MatFactorInfo})::PetscErrorCode
end

function MatICCFactor(arg1, arg2, arg3)
    @ccall libpetsc.MatICCFactor(arg1::Mat, arg2::IS, arg3::Ptr{MatFactorInfo})::PetscErrorCode
end

function MatLUFactorNumeric(arg1, arg2, arg3)
    @ccall libpetsc.MatLUFactorNumeric(arg1::Mat, arg2::Mat, arg3::Ptr{MatFactorInfo})::PetscErrorCode
end

function MatQRFactor(arg1, arg2, arg3)
    @ccall libpetsc.MatQRFactor(arg1::Mat, arg2::IS, arg3::Ptr{MatFactorInfo})::PetscErrorCode
end

function MatQRFactorSymbolic(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatQRFactorSymbolic(arg1::Mat, arg2::Mat, arg3::IS, arg4::Ptr{MatFactorInfo})::PetscErrorCode
end

function MatQRFactorNumeric(arg1, arg2, arg3)
    @ccall libpetsc.MatQRFactorNumeric(arg1::Mat, arg2::Mat, arg3::Ptr{MatFactorInfo})::PetscErrorCode
end

function MatGetInertia(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatGetInertia(arg1::Mat, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function MatSolve(arg1, arg2, arg3)
    @ccall libpetsc.MatSolve(arg1::Mat, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function MatForwardSolve(arg1, arg2, arg3)
    @ccall libpetsc.MatForwardSolve(arg1::Mat, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function MatBackwardSolve(arg1, arg2, arg3)
    @ccall libpetsc.MatBackwardSolve(arg1::Mat, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function MatSolveAdd(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatSolveAdd(arg1::Mat, arg2::Vec, arg3::Vec, arg4::Vec)::PetscErrorCode
end

function MatSolveTranspose(arg1, arg2, arg3)
    @ccall libpetsc.MatSolveTranspose(arg1::Mat, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function MatSolveTransposeAdd(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatSolveTransposeAdd(arg1::Mat, arg2::Vec, arg3::Vec, arg4::Vec)::PetscErrorCode
end

function MatSolves(arg1, arg2, arg3)
    @ccall libpetsc.MatSolves(arg1::Mat, arg2::Vecs, arg3::Vecs)::PetscErrorCode
end

function MatSetUnfactored(arg1)
    @ccall libpetsc.MatSetUnfactored(arg1::Mat)::PetscErrorCode
end

@enum MatFactorSchurStatus::UInt32 begin
    MAT_FACTOR_SCHUR_UNFACTORED = 0
    MAT_FACTOR_SCHUR_FACTORED = 1
    MAT_FACTOR_SCHUR_INVERTED = 2
end

function MatFactorSetSchurIS(arg1, arg2)
    @ccall libpetsc.MatFactorSetSchurIS(arg1::Mat, arg2::IS)::PetscErrorCode
end

function MatFactorGetSchurComplement(arg1, arg2, arg3)
    @ccall libpetsc.MatFactorGetSchurComplement(arg1::Mat, arg2::Ptr{Mat}, arg3::Ptr{MatFactorSchurStatus})::PetscErrorCode
end

function MatFactorRestoreSchurComplement(arg1, arg2, arg3)
    @ccall libpetsc.MatFactorRestoreSchurComplement(arg1::Mat, arg2::Ptr{Mat}, arg3::MatFactorSchurStatus)::PetscErrorCode
end

function MatFactorInvertSchurComplement(arg1)
    @ccall libpetsc.MatFactorInvertSchurComplement(arg1::Mat)::PetscErrorCode
end

function MatFactorCreateSchurComplement(arg1, arg2, arg3)
    @ccall libpetsc.MatFactorCreateSchurComplement(arg1::Mat, arg2::Ptr{Mat}, arg3::Ptr{MatFactorSchurStatus})::PetscErrorCode
end

function MatFactorSolveSchurComplement(arg1, arg2, arg3)
    @ccall libpetsc.MatFactorSolveSchurComplement(arg1::Mat, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function MatFactorSolveSchurComplementTranspose(arg1, arg2, arg3)
    @ccall libpetsc.MatFactorSolveSchurComplementTranspose(arg1::Mat, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function MatFactorFactorizeSchurComplement(arg1)
    @ccall libpetsc.MatFactorFactorizeSchurComplement(arg1::Mat)::PetscErrorCode
end

@enum MatSORType::UInt32 begin
    SOR_FORWARD_SWEEP = 1
    SOR_BACKWARD_SWEEP = 2
    SOR_SYMMETRIC_SWEEP = 3
    SOR_LOCAL_FORWARD_SWEEP = 4
    SOR_LOCAL_BACKWARD_SWEEP = 8
    SOR_LOCAL_SYMMETRIC_SWEEP = 12
    SOR_ZERO_INITIAL_GUESS = 16
    SOR_EISENSTAT = 32
    SOR_APPLY_UPPER = 64
    SOR_APPLY_LOWER = 128
end

function MatSOR(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.MatSOR(arg1::Mat, arg2::Vec, arg3::PetscReal, arg4::MatSORType, arg5::PetscReal, arg6::PetscInt, arg7::PetscInt, arg8::Vec)::PetscErrorCode
end

mutable struct _p_MatColoring end

const MatColoring = Ptr{_p_MatColoring}

const MatColoringType = Ptr{Cchar}

@enum MatColoringWeightType::UInt32 begin
    MAT_COLORING_WEIGHT_RANDOM = 0
    MAT_COLORING_WEIGHT_LEXICAL = 1
    MAT_COLORING_WEIGHT_LF = 2
    MAT_COLORING_WEIGHT_SL = 3
end

function MatColoringCreate(arg1, arg2)
    @ccall libpetsc.MatColoringCreate(arg1::Mat, arg2::Ptr{MatColoring})::PetscErrorCode
end

function MatColoringGetDegrees(arg1, arg2, arg3)
    @ccall libpetsc.MatColoringGetDegrees(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function MatColoringDestroy(arg1)
    @ccall libpetsc.MatColoringDestroy(arg1::Ptr{MatColoring})::PetscErrorCode
end

function MatColoringView(arg1, arg2)
    @ccall libpetsc.MatColoringView(arg1::MatColoring, arg2::PetscViewer)::PetscErrorCode
end

function MatColoringSetType(arg1, arg2)
    @ccall libpetsc.MatColoringSetType(arg1::MatColoring, arg2::MatColoringType)::PetscErrorCode
end

function MatColoringSetFromOptions(arg1)
    @ccall libpetsc.MatColoringSetFromOptions(arg1::MatColoring)::PetscErrorCode
end

function MatColoringSetDistance(arg1, arg2)
    @ccall libpetsc.MatColoringSetDistance(arg1::MatColoring, arg2::PetscInt)::PetscErrorCode
end

function MatColoringGetDistance(arg1, arg2)
    @ccall libpetsc.MatColoringGetDistance(arg1::MatColoring, arg2::Ptr{PetscInt})::PetscErrorCode
end

function MatColoringSetMaxColors(arg1, arg2)
    @ccall libpetsc.MatColoringSetMaxColors(arg1::MatColoring, arg2::PetscInt)::PetscErrorCode
end

function MatColoringGetMaxColors(arg1, arg2)
    @ccall libpetsc.MatColoringGetMaxColors(arg1::MatColoring, arg2::Ptr{PetscInt})::PetscErrorCode
end

function MatColoringApply(arg1, arg2)
    @ccall libpetsc.MatColoringApply(arg1::MatColoring, arg2::Ptr{ISColoring})::PetscErrorCode
end

function MatColoringRegister(arg1, arg2)
    @ccall libpetsc.MatColoringRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function MatColoringPatch(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatColoringPatch(arg1::Mat, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{ISColoringValue}, arg5::Ptr{ISColoring})::PetscErrorCode
end

function MatColoringSetWeightType(arg1, arg2)
    @ccall libpetsc.MatColoringSetWeightType(arg1::MatColoring, arg2::MatColoringWeightType)::PetscErrorCode
end

function MatColoringSetWeights(arg1, arg2, arg3)
    @ccall libpetsc.MatColoringSetWeights(arg1::MatColoring, arg2::Ptr{PetscReal}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function MatColoringCreateWeights(arg1, arg2, lperm)
    @ccall libpetsc.MatColoringCreateWeights(arg1::MatColoring, arg2::Ptr{Ptr{PetscReal}}, lperm::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function MatColoringTest(arg1, arg2)
    @ccall libpetsc.MatColoringTest(arg1::MatColoring, arg2::ISColoring)::PetscErrorCode
end

function MatColoringTestValid(matcoloring, iscoloring)
    @ccall libpetsc.MatColoringTestValid(matcoloring::MatColoring, iscoloring::ISColoring)::PetscErrorCode
end

function MatISColoringTest(arg1, arg2)
    @ccall libpetsc.MatISColoringTest(arg1::Mat, arg2::ISColoring)::PetscErrorCode
end

mutable struct _p_MatFDColoring end

const MatFDColoring = Ptr{_p_MatFDColoring}

function MatFDColoringCreate(arg1, arg2, arg3)
    @ccall libpetsc.MatFDColoringCreate(arg1::Mat, arg2::ISColoring, arg3::Ptr{MatFDColoring})::PetscErrorCode
end

function MatFDColoringDestroy(arg1)
    @ccall libpetsc.MatFDColoringDestroy(arg1::Ptr{MatFDColoring})::PetscErrorCode
end

function MatFDColoringView(arg1, arg2)
    @ccall libpetsc.MatFDColoringView(arg1::MatFDColoring, arg2::PetscViewer)::PetscErrorCode
end

function MatFDColoringSetFunction(arg1, arg2, arg3)
    @ccall libpetsc.MatFDColoringSetFunction(arg1::MatFDColoring, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function MatFDColoringGetFunction(arg1, arg2, arg3)
    @ccall libpetsc.MatFDColoringGetFunction(arg1::MatFDColoring, arg2::Ptr{Ptr{Cvoid}}, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function MatFDColoringSetParameters(arg1, arg2, arg3)
    @ccall libpetsc.MatFDColoringSetParameters(arg1::MatFDColoring, arg2::PetscReal, arg3::PetscReal)::PetscErrorCode
end

function MatFDColoringSetFromOptions(arg1)
    @ccall libpetsc.MatFDColoringSetFromOptions(arg1::MatFDColoring)::PetscErrorCode
end

function MatFDColoringApply(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatFDColoringApply(arg1::Mat, arg2::MatFDColoring, arg3::Vec, arg4::Ptr{Cvoid})::PetscErrorCode
end

function MatFDColoringSetF(arg1, arg2)
    @ccall libpetsc.MatFDColoringSetF(arg1::MatFDColoring, arg2::Vec)::PetscErrorCode
end

function MatFDColoringGetPerturbedColumns(arg1, arg2, arg3)
    @ccall libpetsc.MatFDColoringGetPerturbedColumns(arg1::MatFDColoring, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function MatFDColoringSetUp(arg1, arg2, arg3)
    @ccall libpetsc.MatFDColoringSetUp(arg1::Mat, arg2::ISColoring, arg3::MatFDColoring)::PetscErrorCode
end

function MatFDColoringSetBlockSize(arg1, arg2, arg3)
    @ccall libpetsc.MatFDColoringSetBlockSize(arg1::MatFDColoring, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function MatFDColoringSetValues(arg1, arg2, arg3)
    @ccall libpetsc.MatFDColoringSetValues(arg1::Mat, arg2::MatFDColoring, arg3::Ptr{PetscScalar})::PetscErrorCode
end

mutable struct _p_MatTransposeColoring end

const MatTransposeColoring = Ptr{_p_MatTransposeColoring}

function MatTransposeColoringCreate(arg1, arg2, arg3)
    @ccall libpetsc.MatTransposeColoringCreate(arg1::Mat, arg2::ISColoring, arg3::Ptr{MatTransposeColoring})::PetscErrorCode
end

function MatTransColoringApplySpToDen(arg1, arg2, arg3)
    @ccall libpetsc.MatTransColoringApplySpToDen(arg1::MatTransposeColoring, arg2::Mat, arg3::Mat)::PetscErrorCode
end

function MatTransColoringApplyDenToSp(arg1, arg2, arg3)
    @ccall libpetsc.MatTransColoringApplyDenToSp(arg1::MatTransposeColoring, arg2::Mat, arg3::Mat)::PetscErrorCode
end

function MatTransposeColoringDestroy(arg1)
    @ccall libpetsc.MatTransposeColoringDestroy(arg1::Ptr{MatTransposeColoring})::PetscErrorCode
end

mutable struct _p_MatPartitioning end

const MatPartitioning = Ptr{_p_MatPartitioning}

const MatPartitioningType = Ptr{Cchar}

function MatPartitioningCreate(arg1, arg2)
    @ccall libpetsc.MatPartitioningCreate(arg1::MPI_Comm, arg2::Ptr{MatPartitioning})::PetscErrorCode
end

function MatPartitioningSetType(arg1, arg2)
    @ccall libpetsc.MatPartitioningSetType(arg1::MatPartitioning, arg2::MatPartitioningType)::PetscErrorCode
end

function MatPartitioningSetNParts(arg1, arg2)
    @ccall libpetsc.MatPartitioningSetNParts(arg1::MatPartitioning, arg2::PetscInt)::PetscErrorCode
end

function MatPartitioningSetAdjacency(arg1, arg2)
    @ccall libpetsc.MatPartitioningSetAdjacency(arg1::MatPartitioning, arg2::Mat)::PetscErrorCode
end

function MatPartitioningSetVertexWeights(arg1, arg2)
    @ccall libpetsc.MatPartitioningSetVertexWeights(arg1::MatPartitioning, arg2::Ptr{PetscInt})::PetscErrorCode
end

function MatPartitioningSetPartitionWeights(arg1, arg2)
    @ccall libpetsc.MatPartitioningSetPartitionWeights(arg1::MatPartitioning, arg2::Ptr{PetscReal})::PetscErrorCode
end

function MatPartitioningSetUseEdgeWeights(arg1, arg2)
    @ccall libpetsc.MatPartitioningSetUseEdgeWeights(arg1::MatPartitioning, arg2::PetscBool)::PetscErrorCode
end

function MatPartitioningGetUseEdgeWeights(arg1, arg2)
    @ccall libpetsc.MatPartitioningGetUseEdgeWeights(arg1::MatPartitioning, arg2::Ptr{PetscBool})::PetscErrorCode
end

function MatPartitioningApply(arg1, arg2)
    @ccall libpetsc.MatPartitioningApply(arg1::MatPartitioning, arg2::Ptr{IS})::PetscErrorCode
end

function MatPartitioningImprove(arg1, arg2)
    @ccall libpetsc.MatPartitioningImprove(arg1::MatPartitioning, arg2::Ptr{IS})::PetscErrorCode
end

function MatPartitioningViewImbalance(arg1, arg2)
    @ccall libpetsc.MatPartitioningViewImbalance(arg1::MatPartitioning, arg2::IS)::PetscErrorCode
end

function MatPartitioningApplyND(arg1, arg2)
    @ccall libpetsc.MatPartitioningApplyND(arg1::MatPartitioning, arg2::Ptr{IS})::PetscErrorCode
end

function MatPartitioningDestroy(arg1)
    @ccall libpetsc.MatPartitioningDestroy(arg1::Ptr{MatPartitioning})::PetscErrorCode
end

function MatPartitioningRegister(arg1, arg2)
    @ccall libpetsc.MatPartitioningRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function MatPartitioningView(arg1, arg2)
    @ccall libpetsc.MatPartitioningView(arg1::MatPartitioning, arg2::PetscViewer)::PetscErrorCode
end

function MatPartitioningViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.MatPartitioningViewFromOptions(arg1::MatPartitioning, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function MatPartitioningSetFromOptions(arg1)
    @ccall libpetsc.MatPartitioningSetFromOptions(arg1::MatPartitioning)::PetscErrorCode
end

function MatPartitioningGetType(arg1, arg2)
    @ccall libpetsc.MatPartitioningGetType(arg1::MatPartitioning, arg2::Ptr{MatPartitioningType})::PetscErrorCode
end

function MatPartitioningParmetisSetRepartition(part)
    @ccall libpetsc.MatPartitioningParmetisSetRepartition(part::MatPartitioning)::PetscErrorCode
end

function MatPartitioningParmetisSetCoarseSequential(arg1)
    @ccall libpetsc.MatPartitioningParmetisSetCoarseSequential(arg1::MatPartitioning)::PetscErrorCode
end

function MatPartitioningParmetisGetEdgeCut(arg1, arg2)
    @ccall libpetsc.MatPartitioningParmetisGetEdgeCut(arg1::MatPartitioning, arg2::Ptr{PetscInt})::PetscErrorCode
end

@enum MPChacoGlobalType::UInt32 begin
    MP_CHACO_MULTILEVEL = 1
    MP_CHACO_SPECTRAL = 2
    MP_CHACO_LINEAR = 4
    MP_CHACO_RANDOM = 5
    MP_CHACO_SCATTERED = 6
end

@enum MPChacoLocalType::UInt32 begin
    MP_CHACO_KERNIGHAN = 1
    MP_CHACO_NONE = 2
end

@enum MPChacoEigenType::UInt32 begin
    MP_CHACO_LANCZOS = 0
    MP_CHACO_RQI = 1
end

function MatPartitioningChacoSetGlobal(arg1, arg2)
    @ccall libpetsc.MatPartitioningChacoSetGlobal(arg1::MatPartitioning, arg2::MPChacoGlobalType)::PetscErrorCode
end

function MatPartitioningChacoGetGlobal(arg1, arg2)
    @ccall libpetsc.MatPartitioningChacoGetGlobal(arg1::MatPartitioning, arg2::Ptr{MPChacoGlobalType})::PetscErrorCode
end

function MatPartitioningChacoSetLocal(arg1, arg2)
    @ccall libpetsc.MatPartitioningChacoSetLocal(arg1::MatPartitioning, arg2::MPChacoLocalType)::PetscErrorCode
end

function MatPartitioningChacoGetLocal(arg1, arg2)
    @ccall libpetsc.MatPartitioningChacoGetLocal(arg1::MatPartitioning, arg2::Ptr{MPChacoLocalType})::PetscErrorCode
end

function MatPartitioningChacoSetCoarseLevel(arg1, arg2)
    @ccall libpetsc.MatPartitioningChacoSetCoarseLevel(arg1::MatPartitioning, arg2::PetscReal)::PetscErrorCode
end

function MatPartitioningChacoSetEigenSolver(arg1, arg2)
    @ccall libpetsc.MatPartitioningChacoSetEigenSolver(arg1::MatPartitioning, arg2::MPChacoEigenType)::PetscErrorCode
end

function MatPartitioningChacoGetEigenSolver(arg1, arg2)
    @ccall libpetsc.MatPartitioningChacoGetEigenSolver(arg1::MatPartitioning, arg2::Ptr{MPChacoEigenType})::PetscErrorCode
end

function MatPartitioningChacoSetEigenTol(arg1, arg2)
    @ccall libpetsc.MatPartitioningChacoSetEigenTol(arg1::MatPartitioning, arg2::PetscReal)::PetscErrorCode
end

function MatPartitioningChacoGetEigenTol(arg1, arg2)
    @ccall libpetsc.MatPartitioningChacoGetEigenTol(arg1::MatPartitioning, arg2::Ptr{PetscReal})::PetscErrorCode
end

function MatPartitioningChacoSetEigenNumber(arg1, arg2)
    @ccall libpetsc.MatPartitioningChacoSetEigenNumber(arg1::MatPartitioning, arg2::PetscInt)::PetscErrorCode
end

function MatPartitioningChacoGetEigenNumber(arg1, arg2)
    @ccall libpetsc.MatPartitioningChacoGetEigenNumber(arg1::MatPartitioning, arg2::Ptr{PetscInt})::PetscErrorCode
end

function MatPartitioningPartySetGlobal(arg1, arg2)
    @ccall libpetsc.MatPartitioningPartySetGlobal(arg1::MatPartitioning, arg2::Ptr{Cchar})::PetscErrorCode
end

function MatPartitioningPartySetLocal(arg1, arg2)
    @ccall libpetsc.MatPartitioningPartySetLocal(arg1::MatPartitioning, arg2::Ptr{Cchar})::PetscErrorCode
end

function MatPartitioningPartySetCoarseLevel(arg1, arg2)
    @ccall libpetsc.MatPartitioningPartySetCoarseLevel(arg1::MatPartitioning, arg2::PetscReal)::PetscErrorCode
end

function MatPartitioningPartySetBipart(arg1, arg2)
    @ccall libpetsc.MatPartitioningPartySetBipart(arg1::MatPartitioning, arg2::PetscBool)::PetscErrorCode
end

function MatPartitioningPartySetMatchOptimization(arg1, arg2)
    @ccall libpetsc.MatPartitioningPartySetMatchOptimization(arg1::MatPartitioning, arg2::PetscBool)::PetscErrorCode
end

@enum MPPTScotchStrategyType::UInt32 begin
    MP_PTSCOTCH_DEFAULT = 0
    MP_PTSCOTCH_QUALITY = 1
    MP_PTSCOTCH_SPEED = 2
    MP_PTSCOTCH_BALANCE = 3
    MP_PTSCOTCH_SAFETY = 4
    MP_PTSCOTCH_SCALABILITY = 5
end

function MatPartitioningPTScotchSetImbalance(arg1, arg2)
    @ccall libpetsc.MatPartitioningPTScotchSetImbalance(arg1::MatPartitioning, arg2::PetscReal)::PetscErrorCode
end

function MatPartitioningPTScotchGetImbalance(arg1, arg2)
    @ccall libpetsc.MatPartitioningPTScotchGetImbalance(arg1::MatPartitioning, arg2::Ptr{PetscReal})::PetscErrorCode
end

function MatPartitioningPTScotchSetStrategy(arg1, arg2)
    @ccall libpetsc.MatPartitioningPTScotchSetStrategy(arg1::MatPartitioning, arg2::MPPTScotchStrategyType)::PetscErrorCode
end

function MatPartitioningPTScotchGetStrategy(arg1, arg2)
    @ccall libpetsc.MatPartitioningPTScotchGetStrategy(arg1::MatPartitioning, arg2::Ptr{MPPTScotchStrategyType})::PetscErrorCode
end

function MatPartitioningHierarchicalGetFineparts(arg1, arg2)
    @ccall libpetsc.MatPartitioningHierarchicalGetFineparts(arg1::MatPartitioning, arg2::Ptr{IS})::PetscErrorCode
end

function MatPartitioningHierarchicalGetCoarseparts(arg1, arg2)
    @ccall libpetsc.MatPartitioningHierarchicalGetCoarseparts(arg1::MatPartitioning, arg2::Ptr{IS})::PetscErrorCode
end

function MatPartitioningHierarchicalSetNcoarseparts(arg1, arg2)
    @ccall libpetsc.MatPartitioningHierarchicalSetNcoarseparts(arg1::MatPartitioning, arg2::PetscInt)::PetscErrorCode
end

function MatPartitioningHierarchicalSetNfineparts(arg1, arg2)
    @ccall libpetsc.MatPartitioningHierarchicalSetNfineparts(arg1::MatPartitioning, arg2::PetscInt)::PetscErrorCode
end

function MatMeshToVertexGraph(arg1, arg2, arg3)
    @ccall libpetsc.MatMeshToVertexGraph(arg1::Mat, arg2::PetscInt, arg3::Ptr{Mat})::PetscErrorCode
end

function MatMeshToCellGraph(arg1, arg2, arg3)
    @ccall libpetsc.MatMeshToCellGraph(arg1::Mat, arg2::PetscInt, arg3::Ptr{Mat})::PetscErrorCode
end

@enum MatOperation::UInt32 begin
    MATOP_SET_VALUES = 0
    MATOP_GET_ROW = 1
    MATOP_RESTORE_ROW = 2
    MATOP_MULT = 3
    MATOP_MULT_ADD = 4
    MATOP_MULT_TRANSPOSE = 5
    MATOP_MULT_TRANSPOSE_ADD = 6
    MATOP_SOLVE = 7
    MATOP_SOLVE_ADD = 8
    MATOP_SOLVE_TRANSPOSE = 9
    MATOP_SOLVE_TRANSPOSE_ADD = 10
    MATOP_LUFACTOR = 11
    MATOP_CHOLESKYFACTOR = 12
    MATOP_SOR = 13
    MATOP_TRANSPOSE = 14
    MATOP_GETINFO = 15
    MATOP_EQUAL = 16
    MATOP_GET_DIAGONAL = 17
    MATOP_DIAGONAL_SCALE = 18
    MATOP_NORM = 19
    MATOP_ASSEMBLY_BEGIN = 20
    MATOP_ASSEMBLY_END = 21
    MATOP_SET_OPTION = 22
    MATOP_ZERO_ENTRIES = 23
    MATOP_ZERO_ROWS = 24
    MATOP_LUFACTOR_SYMBOLIC = 25
    MATOP_LUFACTOR_NUMERIC = 26
    MATOP_CHOLESKY_FACTOR_SYMBOLIC = 27
    MATOP_CHOLESKY_FACTOR_NUMERIC = 28
    MATOP_SETUP_PREALLOCATION = 29
    MATOP_ILUFACTOR_SYMBOLIC = 30
    MATOP_ICCFACTOR_SYMBOLIC = 31
    MATOP_GET_DIAGONAL_BLOCK = 32
    MATOP_FREE_INTER_STRUCT = 33
    MATOP_DUPLICATE = 34
    MATOP_FORWARD_SOLVE = 35
    MATOP_BACKWARD_SOLVE = 36
    MATOP_ILUFACTOR = 37
    MATOP_ICCFACTOR = 38
    MATOP_AXPY = 39
    MATOP_CREATE_SUBMATRICES = 40
    MATOP_INCREASE_OVERLAP = 41
    MATOP_GET_VALUES = 42
    MATOP_COPY = 43
    MATOP_GET_ROW_MAX = 44
    MATOP_SCALE = 45
    MATOP_SHIFT = 46
    MATOP_DIAGONAL_SET = 47
    MATOP_ZERO_ROWS_COLUMNS = 48
    MATOP_SET_RANDOM = 49
    MATOP_GET_ROW_IJ = 50
    MATOP_RESTORE_ROW_IJ = 51
    MATOP_GET_COLUMN_IJ = 52
    MATOP_RESTORE_COLUMN_IJ = 53
    MATOP_FDCOLORING_CREATE = 54
    MATOP_COLORING_PATCH = 55
    MATOP_SET_UNFACTORED = 56
    MATOP_PERMUTE = 57
    MATOP_SET_VALUES_BLOCKED = 58
    MATOP_CREATE_SUBMATRIX = 59
    MATOP_DESTROY = 60
    MATOP_VIEW = 61
    MATOP_CONVERT_FROM = 62
    MATOP_MATMAT_MULT = 63
    MATOP_MATMAT_MULT_SYMBOLIC = 64
    MATOP_MATMAT_MULT_NUMERIC = 65
    MATOP_SET_LOCAL_TO_GLOBAL_MAP = 66
    MATOP_SET_VALUES_LOCAL = 67
    MATOP_ZERO_ROWS_LOCAL = 68
    MATOP_GET_ROW_MAX_ABS = 69
    MATOP_GET_ROW_MIN_ABS = 70
    MATOP_CONVERT = 71
    MATOP_SET_COLORING = 72
    MATOP_SET_VALUES_ADIFOR = 74
    MATOP_FD_COLORING_APPLY = 75
    MATOP_SET_FROM_OPTIONS = 76
    MATOP_MULT_CONSTRAINED = 77
    MATOP_MULT_TRANSPOSE_CONSTRAIN = 78
    MATOP_FIND_ZERO_DIAGONALS = 79
    MATOP_MULT_MULTIPLE = 80
    MATOP_SOLVE_MULTIPLE = 81
    MATOP_GET_INERTIA = 82
    MATOP_LOAD = 83
    MATOP_IS_SYMMETRIC = 84
    MATOP_IS_HERMITIAN = 85
    MATOP_IS_STRUCTURALLY_SYMMETRIC = 86
    MATOP_SET_VALUES_BLOCKEDLOCAL = 87
    MATOP_CREATE_VECS = 88
    MATOP_MAT_MULT = 89
    MATOP_MAT_MULT_SYMBOLIC = 90
    MATOP_MAT_MULT_NUMERIC = 91
    MATOP_PTAP = 92
    MATOP_PTAP_SYMBOLIC = 93
    MATOP_PTAP_NUMERIC = 94
    MATOP_MAT_TRANSPOSE_MULT = 95
    MATOP_MAT_TRANSPOSE_MULT_SYMBO = 96
    MATOP_MAT_TRANSPOSE_MULT_NUMER = 97
    MATOP_PRODUCTSETFROMOPTIONS = 99
    MATOP_PRODUCTSYMBOLIC = 100
    MATOP_PRODUCTNUMERIC = 101
    MATOP_CONJUGATE = 102
    MATOP_SET_VALUES_ROW = 104
    MATOP_REAL_PART = 105
    MATOP_IMAGINARY_PART = 106
    MATOP_GET_ROW_UPPER_TRIANGULAR = 107
    MATOP_RESTORE_ROW_UPPER_TRIANG = 108
    MATOP_MAT_SOLVE = 109
    MATOP_MAT_SOLVE_TRANSPOSE = 110
    MATOP_GET_ROW_MIN = 111
    MATOP_GET_COLUMN_VECTOR = 112
    MATOP_MISSING_DIAGONAL = 113
    MATOP_GET_SEQ_NONZERO_STRUCTUR = 114
    MATOP_CREATE = 115
    MATOP_GET_GHOSTS = 116
    MATOP_GET_LOCAL_SUB_MATRIX = 117
    MATOP_RESTORE_LOCALSUB_MATRIX = 118
    MATOP_MULT_DIAGONAL_BLOCK = 119
    MATOP_HERMITIAN_TRANSPOSE = 120
    MATOP_MULT_HERMITIAN_TRANSPOSE = 121
    MATOP_MULT_HERMITIAN_TRANS_ADD = 122
    MATOP_GET_MULTI_PROC_BLOCK = 123
    MATOP_FIND_NONZERO_ROWS = 124
    MATOP_GET_COLUMN_NORMS = 125
    MATOP_INVERT_BLOCK_DIAGONAL = 126
    MATOP_CREATE_SUB_MATRICES_MPI = 128
    MATOP_SET_VALUES_BATCH = 129
    MATOP_TRANSPOSE_MAT_MULT = 130
    MATOP_TRANSPOSE_MAT_MULT_SYMBO = 131
    MATOP_TRANSPOSE_MAT_MULT_NUMER = 132
    MATOP_TRANSPOSE_COLORING_CREAT = 133
    MATOP_TRANS_COLORING_APPLY_SPT = 134
    MATOP_TRANS_COLORING_APPLY_DEN = 135
    MATOP_RART = 136
    MATOP_RART_SYMBOLIC = 137
    MATOP_RART_NUMERIC = 138
    MATOP_SET_BLOCK_SIZES = 139
    MATOP_AYPX = 140
    MATOP_RESIDUAL = 141
    MATOP_FDCOLORING_SETUP = 142
    MATOP_MPICONCATENATESEQ = 144
    MATOP_DESTROYSUBMATRICES = 145
    MATOP_TRANSPOSE_SOLVE = 146
    MATOP_GET_VALUES_LOCAL = 147
end

function MatSetOperation(arg1, arg2, arg3)
    @ccall libpetsc.MatSetOperation(arg1::Mat, arg2::MatOperation, arg3::Ptr{Cvoid})::PetscErrorCode
end

function MatGetOperation(arg1, arg2, arg3)
    @ccall libpetsc.MatGetOperation(arg1::Mat, arg2::MatOperation, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function MatHasOperation(arg1, arg2, arg3)
    @ccall libpetsc.MatHasOperation(arg1::Mat, arg2::MatOperation, arg3::Ptr{PetscBool})::PetscErrorCode
end

function MatHasCongruentLayouts(arg1, arg2)
    @ccall libpetsc.MatHasCongruentLayouts(arg1::Mat, arg2::Ptr{PetscBool})::PetscErrorCode
end

function MatFreeIntermediateDataStructures(A)
    @ccall libpetsc.MatFreeIntermediateDataStructures(A::Mat)::PetscErrorCode
end

function MatShellSetOperation(arg1, arg2, arg3)
    @ccall libpetsc.MatShellSetOperation(arg1::Mat, arg2::MatOperation, arg3::Ptr{Cvoid})::PetscErrorCode
end

function MatShellGetOperation(arg1, arg2, arg3)
    @ccall libpetsc.MatShellGetOperation(arg1::Mat, arg2::MatOperation, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function MatShellSetContext(arg1, arg2)
    @ccall libpetsc.MatShellSetContext(arg1::Mat, arg2::Ptr{Cvoid})::PetscErrorCode
end

function MatShellSetVecType(arg1, arg2)
    @ccall libpetsc.MatShellSetVecType(arg1::Mat, arg2::VecType)::PetscErrorCode
end

function MatShellTestMult(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatShellTestMult(arg1::Mat, arg2::Ptr{Cvoid}, arg3::Vec, arg4::Ptr{Cvoid}, arg5::Ptr{PetscBool})::PetscErrorCode
end

function MatShellTestMultTranspose(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatShellTestMultTranspose(arg1::Mat, arg2::Ptr{Cvoid}, arg3::Vec, arg4::Ptr{Cvoid}, arg5::Ptr{PetscBool})::PetscErrorCode
end

function MatShellSetManageScalingShifts(arg1)
    @ccall libpetsc.MatShellSetManageScalingShifts(arg1::Mat)::PetscErrorCode
end

function MatShellSetMatProductOperation(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.MatShellSetMatProductOperation(arg1::Mat, arg2::MatProductType, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid}, arg6::MatType, arg7::MatType)::PetscErrorCode
end

function MatIsShell(arg1, arg2)
    @ccall libpetsc.MatIsShell(arg1::Mat, arg2::Ptr{PetscBool})::PetscErrorCode
end

function MatMPIBAIJSetHashTableFactor(arg1, arg2)
    @ccall libpetsc.MatMPIBAIJSetHashTableFactor(arg1::Mat, arg2::PetscReal)::PetscErrorCode
end

function MatISSetLocalMatType(arg1, arg2)
    @ccall libpetsc.MatISSetLocalMatType(arg1::Mat, arg2::MatType)::PetscErrorCode
end

function MatISSetPreallocation(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatISSetPreallocation(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::PetscInt, arg5::Ptr{PetscInt})::PetscErrorCode
end

function MatISStoreL2L(arg1, arg2)
    @ccall libpetsc.MatISStoreL2L(arg1::Mat, arg2::PetscBool)::PetscErrorCode
end

function MatISFixLocalEmpty(arg1, arg2)
    @ccall libpetsc.MatISFixLocalEmpty(arg1::Mat, arg2::PetscBool)::PetscErrorCode
end

function MatISGetLocalMat(arg1, arg2)
    @ccall libpetsc.MatISGetLocalMat(arg1::Mat, arg2::Ptr{Mat})::PetscErrorCode
end

function MatISRestoreLocalMat(arg1, arg2)
    @ccall libpetsc.MatISRestoreLocalMat(arg1::Mat, arg2::Ptr{Mat})::PetscErrorCode
end

function MatISSetLocalMat(arg1, arg2)
    @ccall libpetsc.MatISSetLocalMat(arg1::Mat, arg2::Mat)::PetscErrorCode
end

function MatISGetMPIXAIJ(arg1, arg2, arg3)
    @ccall libpetsc.MatISGetMPIXAIJ(arg1::Mat, arg2::MatReuse, arg3::Ptr{Mat})::PetscErrorCode
end

mutable struct _p_MatNullSpace end

const MatNullSpace = Ptr{_p_MatNullSpace}

function MatNullSpaceCreate(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatNullSpaceCreate(arg1::MPI_Comm, arg2::PetscBool, arg3::PetscInt, arg4::Ptr{Vec}, arg5::Ptr{MatNullSpace})::PetscErrorCode
end

function MatNullSpaceSetFunction(arg1, arg2, arg3)
    @ccall libpetsc.MatNullSpaceSetFunction(arg1::MatNullSpace, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function MatNullSpaceDestroy(arg1)
    @ccall libpetsc.MatNullSpaceDestroy(arg1::Ptr{MatNullSpace})::PetscErrorCode
end

function MatNullSpaceRemove(arg1, arg2)
    @ccall libpetsc.MatNullSpaceRemove(arg1::MatNullSpace, arg2::Vec)::PetscErrorCode
end

function MatGetNullSpace(arg1, arg2)
    @ccall libpetsc.MatGetNullSpace(arg1::Mat, arg2::Ptr{MatNullSpace})::PetscErrorCode
end

function MatGetTransposeNullSpace(arg1, arg2)
    @ccall libpetsc.MatGetTransposeNullSpace(arg1::Mat, arg2::Ptr{MatNullSpace})::PetscErrorCode
end

function MatSetTransposeNullSpace(arg1, arg2)
    @ccall libpetsc.MatSetTransposeNullSpace(arg1::Mat, arg2::MatNullSpace)::PetscErrorCode
end

function MatSetNullSpace(arg1, arg2)
    @ccall libpetsc.MatSetNullSpace(arg1::Mat, arg2::MatNullSpace)::PetscErrorCode
end

function MatSetNearNullSpace(arg1, arg2)
    @ccall libpetsc.MatSetNearNullSpace(arg1::Mat, arg2::MatNullSpace)::PetscErrorCode
end

function MatGetNearNullSpace(arg1, arg2)
    @ccall libpetsc.MatGetNearNullSpace(arg1::Mat, arg2::Ptr{MatNullSpace})::PetscErrorCode
end

function MatNullSpaceTest(arg1, arg2, arg3)
    @ccall libpetsc.MatNullSpaceTest(arg1::MatNullSpace, arg2::Mat, arg3::Ptr{PetscBool})::PetscErrorCode
end

function MatNullSpaceView(arg1, arg2)
    @ccall libpetsc.MatNullSpaceView(arg1::MatNullSpace, arg2::PetscViewer)::PetscErrorCode
end

function MatNullSpaceGetVecs(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatNullSpaceGetVecs(arg1::MatNullSpace, arg2::Ptr{PetscBool}, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{Vec}})::PetscErrorCode
end

function MatNullSpaceCreateRigidBody(arg1, arg2)
    @ccall libpetsc.MatNullSpaceCreateRigidBody(arg1::Vec, arg2::Ptr{MatNullSpace})::PetscErrorCode
end

function MatReorderingSeqSBAIJ(arg1, arg2)
    @ccall libpetsc.MatReorderingSeqSBAIJ(arg1::Mat, arg2::IS)::PetscErrorCode
end

function MatMPISBAIJSetHashTableFactor(arg1, arg2)
    @ccall libpetsc.MatMPISBAIJSetHashTableFactor(arg1::Mat, arg2::PetscReal)::PetscErrorCode
end

function MatSeqSBAIJSetColumnIndices(arg1, arg2)
    @ccall libpetsc.MatSeqSBAIJSetColumnIndices(arg1::Mat, arg2::Ptr{PetscInt})::PetscErrorCode
end

function MatCreateMAIJ(arg1, arg2, arg3)
    @ccall libpetsc.MatCreateMAIJ(arg1::Mat, arg2::PetscInt, arg3::Ptr{Mat})::PetscErrorCode
end

function MatMAIJRedimension(arg1, arg2, arg3)
    @ccall libpetsc.MatMAIJRedimension(arg1::Mat, arg2::PetscInt, arg3::Ptr{Mat})::PetscErrorCode
end

function MatMAIJGetAIJ(arg1, arg2)
    @ccall libpetsc.MatMAIJGetAIJ(arg1::Mat, arg2::Ptr{Mat})::PetscErrorCode
end

function MatComputeOperator(arg1, arg2, arg3)
    @ccall libpetsc.MatComputeOperator(arg1::Mat, arg2::MatType, arg3::Ptr{Mat})::PetscErrorCode
end

function MatComputeOperatorTranspose(arg1, arg2, arg3)
    @ccall libpetsc.MatComputeOperatorTranspose(arg1::Mat, arg2::MatType, arg3::Ptr{Mat})::PetscErrorCode
end

function MatComputeExplicitOperator(A, B)
    @ccall libpetsc.MatComputeExplicitOperator(A::Mat, B::Ptr{Mat})::PetscErrorCode
end

function MatComputeExplicitOperatorTranspose(A, B)
    @ccall libpetsc.MatComputeExplicitOperatorTranspose(A::Mat, B::Ptr{Mat})::PetscErrorCode
end

function MatCreateKAIJ(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatCreateKAIJ(arg1::Mat, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscScalar}, arg5::Ptr{PetscScalar}, arg6::Ptr{Mat})::PetscErrorCode
end

function MatKAIJGetAIJ(arg1, arg2)
    @ccall libpetsc.MatKAIJGetAIJ(arg1::Mat, arg2::Ptr{Mat})::PetscErrorCode
end

function MatKAIJGetS(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatKAIJGetS(arg1::Mat, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatKAIJGetSRead(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatKAIJGetSRead(arg1::Mat, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatKAIJRestoreS(arg1, arg2)
    @ccall libpetsc.MatKAIJRestoreS(arg1::Mat, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatKAIJRestoreSRead(arg1, arg2)
    @ccall libpetsc.MatKAIJRestoreSRead(arg1::Mat, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatKAIJGetT(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatKAIJGetT(arg1::Mat, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatKAIJGetTRead(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatKAIJGetTRead(arg1::Mat, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatKAIJRestoreT(arg1, arg2)
    @ccall libpetsc.MatKAIJRestoreT(arg1::Mat, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatKAIJRestoreTRead(arg1, arg2)
    @ccall libpetsc.MatKAIJRestoreTRead(arg1::Mat, arg2::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function MatKAIJSetAIJ(arg1, arg2)
    @ccall libpetsc.MatKAIJSetAIJ(arg1::Mat, arg2::Mat)::PetscErrorCode
end

function MatKAIJSetS(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatKAIJSetS(arg1::Mat, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscScalar})::PetscErrorCode
end

function MatKAIJSetT(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatKAIJSetT(arg1::Mat, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscScalar})::PetscErrorCode
end

function MatKAIJGetScaledIdentity(arg1, arg2)
    @ccall libpetsc.MatKAIJGetScaledIdentity(arg1::Mat, arg2::Ptr{PetscBool})::PetscErrorCode
end

function MatDiagonalScaleLocal(arg1, arg2)
    @ccall libpetsc.MatDiagonalScaleLocal(arg1::Mat, arg2::Vec)::PetscErrorCode
end

function MatCreateMFFD(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatCreateMFFD(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{Mat})::PetscErrorCode
end

function MatMFFDSetBase(arg1, arg2, arg3)
    @ccall libpetsc.MatMFFDSetBase(arg1::Mat, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function MatMFFDSetFunction(arg1, arg2, arg3)
    @ccall libpetsc.MatMFFDSetFunction(arg1::Mat, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function MatMFFDSetFunctioni(arg1, arg2)
    @ccall libpetsc.MatMFFDSetFunctioni(arg1::Mat, arg2::Ptr{Cvoid})::PetscErrorCode
end

function MatMFFDSetFunctioniBase(arg1, arg2)
    @ccall libpetsc.MatMFFDSetFunctioniBase(arg1::Mat, arg2::Ptr{Cvoid})::PetscErrorCode
end

function MatMFFDSetHHistory(arg1, arg2, arg3)
    @ccall libpetsc.MatMFFDSetHHistory(arg1::Mat, arg2::Ptr{PetscScalar}, arg3::PetscInt)::PetscErrorCode
end

function MatMFFDResetHHistory(arg1)
    @ccall libpetsc.MatMFFDResetHHistory(arg1::Mat)::PetscErrorCode
end

function MatMFFDSetFunctionError(arg1, arg2)
    @ccall libpetsc.MatMFFDSetFunctionError(arg1::Mat, arg2::PetscReal)::PetscErrorCode
end

function MatMFFDSetPeriod(arg1, arg2)
    @ccall libpetsc.MatMFFDSetPeriod(arg1::Mat, arg2::PetscInt)::PetscErrorCode
end

function MatMFFDGetH(arg1, arg2)
    @ccall libpetsc.MatMFFDGetH(arg1::Mat, arg2::Ptr{PetscScalar})::PetscErrorCode
end

function MatMFFDSetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.MatMFFDSetOptionsPrefix(arg1::Mat, arg2::Ptr{Cchar})::PetscErrorCode
end

function MatMFFDCheckPositivity(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatMFFDCheckPositivity(arg1::Ptr{Cvoid}, arg2::Vec, arg3::Vec, arg4::Ptr{PetscScalar})::PetscErrorCode
end

function MatMFFDSetCheckh(arg1, arg2, arg3)
    @ccall libpetsc.MatMFFDSetCheckh(arg1::Mat, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

mutable struct _p_MatMFFD end

const MatMFFD = Ptr{_p_MatMFFD}

const MatMFFDType = Ptr{Cchar}

function MatMFFDSetType(arg1, arg2)
    @ccall libpetsc.MatMFFDSetType(arg1::Mat, arg2::MatMFFDType)::PetscErrorCode
end

function MatMFFDRegister(arg1, arg2)
    @ccall libpetsc.MatMFFDRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function MatMFFDDSSetUmin(arg1, arg2)
    @ccall libpetsc.MatMFFDDSSetUmin(arg1::Mat, arg2::PetscReal)::PetscErrorCode
end

function MatMFFDWPSetComputeNormU(arg1, arg2)
    @ccall libpetsc.MatMFFDWPSetComputeNormU(arg1::Mat, arg2::PetscBool)::PetscErrorCode
end

function MatFDColoringSetType(arg1, arg2)
    @ccall libpetsc.MatFDColoringSetType(arg1::MatFDColoring, arg2::MatMFFDType)::PetscErrorCode
end

function PetscViewerMathematicaPutMatrix(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscViewerMathematicaPutMatrix(arg1::PetscViewer, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscReal})::PetscErrorCode
end

function PetscViewerMathematicaPutCSRMatrix(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscViewerMathematicaPutCSRMatrix(arg1::PetscViewer, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscReal})::PetscErrorCode
end

function MatBindToCPU(arg1, arg2)
    @ccall libpetsc.MatBindToCPU(arg1::Mat, arg2::PetscBool)::PetscErrorCode
end

function MatPinToCPU(A, flg)
    @ccall libpetsc.MatPinToCPU(A::Mat, flg::PetscBool)::PetscErrorCode
end

mutable struct _p_SplitCSRMat end

const PetscSplitCSRDataStructure = _p_SplitCSRMat

function MatCreateNest(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.MatCreateNest(arg1::MPI_Comm, arg2::PetscInt, arg3::Ptr{IS}, arg4::PetscInt, arg5::Ptr{IS}, arg6::Ptr{Mat}, arg7::Ptr{Mat})::PetscErrorCode
end

function MatNestGetSize(arg1, arg2, arg3)
    @ccall libpetsc.MatNestGetSize(arg1::Mat, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function MatNestGetISs(arg1, arg2, arg3)
    @ccall libpetsc.MatNestGetISs(arg1::Mat, arg2::Ptr{IS}, arg3::Ptr{IS})::PetscErrorCode
end

function MatNestGetLocalISs(arg1, arg2, arg3)
    @ccall libpetsc.MatNestGetLocalISs(arg1::Mat, arg2::Ptr{IS}, arg3::Ptr{IS})::PetscErrorCode
end

function MatNestGetSubMats(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatNestGetSubMats(arg1::Mat, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{Ptr{Mat}}})::PetscErrorCode
end

function MatNestGetSubMat(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatNestGetSubMat(arg1::Mat, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Mat})::PetscErrorCode
end

function MatNestSetVecType(arg1, arg2)
    @ccall libpetsc.MatNestSetVecType(arg1::Mat, arg2::VecType)::PetscErrorCode
end

function MatNestSetSubMats(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatNestSetSubMats(arg1::Mat, arg2::PetscInt, arg3::Ptr{IS}, arg4::PetscInt, arg5::Ptr{IS}, arg6::Ptr{Mat})::PetscErrorCode
end

function MatNestSetSubMat(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatNestSetSubMat(arg1::Mat, arg2::PetscInt, arg3::PetscInt, arg4::Mat)::PetscErrorCode
end

function MatChop(arg1, arg2)
    @ccall libpetsc.MatChop(arg1::Mat, arg2::PetscReal)::PetscErrorCode
end

function MatComputeBandwidth(arg1, arg2, arg3)
    @ccall libpetsc.MatComputeBandwidth(arg1::Mat, arg2::PetscReal, arg3::Ptr{PetscInt})::PetscErrorCode
end

function MatSubdomainsCreateCoalesce(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatSubdomainsCreateCoalesce(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{IS}})::PetscErrorCode
end

function MatPreallocatorPreallocate(arg1, arg2, arg3)
    @ccall libpetsc.MatPreallocatorPreallocate(arg1::Mat, arg2::PetscBool, arg3::Mat)::PetscErrorCode
end

function MatHeaderMerge(arg1, arg2)
    @ccall libpetsc.MatHeaderMerge(arg1::Mat, arg2::Ptr{Mat})::PetscErrorCode
end

function MatHeaderReplace(arg1, arg2)
    @ccall libpetsc.MatHeaderReplace(arg1::Mat, arg2::Ptr{Mat})::PetscErrorCode
end

mutable struct _p_DM end

const DM = Ptr{_p_DM}

@enum DMBoundaryType::UInt32 begin
    DM_BOUNDARY_NONE = 0
    DM_BOUNDARY_GHOSTED = 1
    DM_BOUNDARY_MIRROR = 2
    DM_BOUNDARY_PERIODIC = 3
    DM_BOUNDARY_TWIST = 4
end

@enum DMBoundaryConditionType::UInt32 begin
    DM_BC_ESSENTIAL = 1
    DM_BC_ESSENTIAL_FIELD = 5
    DM_BC_NATURAL = 2
    DM_BC_NATURAL_FIELD = 6
    DM_BC_ESSENTIAL_BD_FIELD = 9
    DM_BC_NATURAL_RIEMANN = 10
end

@enum DMPointLocationType::UInt32 begin
    DM_POINTLOCATION_NONE = 0
    DM_POINTLOCATION_NEAREST = 1
    DM_POINTLOCATION_REMOVE = 2
end

@enum DMAdaptationStrategy::UInt32 begin
    DM_ADAPTATION_INITIAL = 0
    DM_ADAPTATION_SEQUENTIAL = 1
    DM_ADAPTATION_MULTILEVEL = 2
end

@enum DMAdaptationCriterion::UInt32 begin
    DM_ADAPTATION_NONE = 0
    DM_ADAPTATION_REFINE = 1
    DM_ADAPTATION_LABEL = 2
    DM_ADAPTATION_METRIC = 3
end

@enum DMAdaptFlag::Int32 begin
    DM_ADAPT_DETERMINE = -1
    DM_ADAPT_KEEP = 0
    DM_ADAPT_REFINE = 1
    DM_ADAPT_COARSEN = 2
    DM_ADAPT_COARSEN_LAST = 3
    DM_ADAPT_RESERVED_COUNT = 4
end

@enum DMDirection::UInt32 begin
    DM_X = 0
    DM_Y = 1
    DM_Z = 2
end

@enum DMEnclosureType::UInt32 begin
    DM_ENC_EQUALITY = 0
    DM_ENC_SUPERMESH = 1
    DM_ENC_SUBMESH = 2
    DM_ENC_NONE = 3
    DM_ENC_UNKNOWN = 4
end

@enum DMPolytopeType::UInt32 begin
    DM_POLYTOPE_POINT = 0
    DM_POLYTOPE_SEGMENT = 1
    DM_POLYTOPE_POINT_PRISM_TENSOR = 2
    DM_POLYTOPE_TRIANGLE = 3
    DM_POLYTOPE_QUADRILATERAL = 4
    DM_POLYTOPE_SEG_PRISM_TENSOR = 5
    DM_POLYTOPE_TETRAHEDRON = 6
    DM_POLYTOPE_HEXAHEDRON = 7
    DM_POLYTOPE_TRI_PRISM = 8
    DM_POLYTOPE_TRI_PRISM_TENSOR = 9
    DM_POLYTOPE_QUAD_PRISM_TENSOR = 10
    DM_POLYTOPE_PYRAMID = 11
    DM_POLYTOPE_FV_GHOST = 12
    DM_POLYTOPE_INTERIOR_GHOST = 13
    DM_POLYTOPE_UNKNOWN = 14
    DM_NUM_POLYTOPES = 15
end

@enum PetscUnit::UInt32 begin
    PETSC_UNIT_LENGTH = 0
    PETSC_UNIT_MASS = 1
    PETSC_UNIT_TIME = 2
    PETSC_UNIT_CURRENT = 3
    PETSC_UNIT_TEMPERATURE = 4
    PETSC_UNIT_AMOUNT = 5
    PETSC_UNIT_LUMINOSITY = 6
    NUM_PETSC_UNITS = 7
end

mutable struct _p_DMField end

const DMField = Ptr{_p_DMField}

mutable struct _p_UniversalLabel end

const DMUniversalLabel = Ptr{_p_UniversalLabel}

mutable struct _p_PetscSpace end

const PetscSpace = Ptr{_p_PetscSpace}

@enum PetscSpacePolynomialType::UInt32 begin
    PETSCSPACE_POLYNOMIALTYPE_P = 0
    PETSCSPACE_POLYNOMIALTYPE_PMINUS_HDIV = 1
    PETSCSPACE_POLYNOMIALTYPE_PMINUS_HCURL = 2
end

mutable struct _p_PetscDualSpace end

const PetscDualSpace = Ptr{_p_PetscDualSpace}

@enum PetscDualSpaceReferenceCell::UInt32 begin
    PETSCDUALSPACE_REFCELL_SIMPLEX = 0
    PETSCDUALSPACE_REFCELL_TENSOR = 1
end

@enum PetscDualSpaceTransformType::UInt32 begin
    IDENTITY_TRANSFORM = 0
    COVARIANT_PIOLA_TRANSFORM = 1
    CONTRAVARIANT_PIOLA_TRANSFORM = 2
end

mutable struct _p_PetscFE end

const PetscFE = Ptr{_p_PetscFE}

@enum PetscFEJacobianType::UInt32 begin
    PETSCFE_JACOBIAN = 0
    PETSCFE_JACOBIAN_PRE = 1
    PETSCFE_JACOBIAN_DYN = 2
end

mutable struct _p_DMLabel end

const DMLabel = Ptr{_p_DMLabel}

function DMLabelCreate(arg1, arg2, arg3)
    @ccall libpetsc.DMLabelCreate(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{DMLabel})::PetscErrorCode
end

function DMLabelView(arg1, arg2)
    @ccall libpetsc.DMLabelView(arg1::DMLabel, arg2::PetscViewer)::PetscErrorCode
end

function DMLabelReset(arg1)
    @ccall libpetsc.DMLabelReset(arg1::DMLabel)::PetscErrorCode
end

function DMLabelDestroy(arg1)
    @ccall libpetsc.DMLabelDestroy(arg1::Ptr{DMLabel})::PetscErrorCode
end

function DMLabelGetDefaultValue(arg1, arg2)
    @ccall libpetsc.DMLabelGetDefaultValue(arg1::DMLabel, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMLabelSetDefaultValue(arg1, arg2)
    @ccall libpetsc.DMLabelSetDefaultValue(arg1::DMLabel, arg2::PetscInt)::PetscErrorCode
end

function DMLabelDuplicate(arg1, arg2)
    @ccall libpetsc.DMLabelDuplicate(arg1::DMLabel, arg2::Ptr{DMLabel})::PetscErrorCode
end

function DMLabelGetValue(arg1, arg2, arg3)
    @ccall libpetsc.DMLabelGetValue(arg1::DMLabel, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMLabelSetValue(arg1, arg2, arg3)
    @ccall libpetsc.DMLabelSetValue(arg1::DMLabel, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function DMLabelClearValue(arg1, arg2, arg3)
    @ccall libpetsc.DMLabelClearValue(arg1::DMLabel, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function DMLabelAddStratum(arg1, arg2)
    @ccall libpetsc.DMLabelAddStratum(arg1::DMLabel, arg2::PetscInt)::PetscErrorCode
end

function DMLabelAddStrata(arg1, arg2, arg3)
    @ccall libpetsc.DMLabelAddStrata(arg1::DMLabel, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMLabelAddStrataIS(arg1, arg2)
    @ccall libpetsc.DMLabelAddStrataIS(arg1::DMLabel, arg2::IS)::PetscErrorCode
end

function DMLabelInsertIS(arg1, arg2, arg3)
    @ccall libpetsc.DMLabelInsertIS(arg1::DMLabel, arg2::IS, arg3::PetscInt)::PetscErrorCode
end

function DMLabelGetNumValues(arg1, arg2)
    @ccall libpetsc.DMLabelGetNumValues(arg1::DMLabel, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMLabelGetStratumBounds(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMLabelGetStratumBounds(arg1::DMLabel, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMLabelGetValueIS(arg1, arg2)
    @ccall libpetsc.DMLabelGetValueIS(arg1::DMLabel, arg2::Ptr{IS})::PetscErrorCode
end

function DMLabelStratumHasPoint(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMLabelStratumHasPoint(arg1::DMLabel, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscBool})::PetscErrorCode
end

function DMLabelHasStratum(arg1, arg2, arg3)
    @ccall libpetsc.DMLabelHasStratum(arg1::DMLabel, arg2::PetscInt, arg3::Ptr{PetscBool})::PetscErrorCode
end

function DMLabelGetStratumSize(arg1, arg2, arg3)
    @ccall libpetsc.DMLabelGetStratumSize(arg1::DMLabel, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMLabelGetStratumIS(arg1, arg2, arg3)
    @ccall libpetsc.DMLabelGetStratumIS(arg1::DMLabel, arg2::PetscInt, arg3::Ptr{IS})::PetscErrorCode
end

function DMLabelSetStratumIS(arg1, arg2, arg3)
    @ccall libpetsc.DMLabelSetStratumIS(arg1::DMLabel, arg2::PetscInt, arg3::IS)::PetscErrorCode
end

function DMLabelSetStratumBounds(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMLabelSetStratumBounds(arg1::DMLabel, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt)::PetscErrorCode
end

function DMLabelClearStratum(arg1, arg2)
    @ccall libpetsc.DMLabelClearStratum(arg1::DMLabel, arg2::PetscInt)::PetscErrorCode
end

function DMLabelComputeIndex(arg1)
    @ccall libpetsc.DMLabelComputeIndex(arg1::DMLabel)::PetscErrorCode
end

function DMLabelCreateIndex(arg1, arg2, arg3)
    @ccall libpetsc.DMLabelCreateIndex(arg1::DMLabel, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function DMLabelDestroyIndex(arg1)
    @ccall libpetsc.DMLabelDestroyIndex(arg1::DMLabel)::PetscErrorCode
end

function DMLabelHasValue(arg1, arg2, arg3)
    @ccall libpetsc.DMLabelHasValue(arg1::DMLabel, arg2::PetscInt, arg3::Ptr{PetscBool})::PetscErrorCode
end

function DMLabelHasPoint(arg1, arg2, arg3)
    @ccall libpetsc.DMLabelHasPoint(arg1::DMLabel, arg2::PetscInt, arg3::Ptr{PetscBool})::PetscErrorCode
end

function DMLabelGetBounds(arg1, arg2, arg3)
    @ccall libpetsc.DMLabelGetBounds(arg1::DMLabel, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMLabelFilter(arg1, arg2, arg3)
    @ccall libpetsc.DMLabelFilter(arg1::DMLabel, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function DMLabelPermute(arg1, arg2, arg3)
    @ccall libpetsc.DMLabelPermute(arg1::DMLabel, arg2::IS, arg3::Ptr{DMLabel})::PetscErrorCode
end

function DMLabelDistribute(arg1, arg2, arg3)
    @ccall libpetsc.DMLabelDistribute(arg1::DMLabel, arg2::PetscSF, arg3::Ptr{DMLabel})::PetscErrorCode
end

function DMLabelGather(arg1, arg2, arg3)
    @ccall libpetsc.DMLabelGather(arg1::DMLabel, arg2::PetscSF, arg3::Ptr{DMLabel})::PetscErrorCode
end

function DMLabelConvertToSection(arg1, arg2, arg3)
    @ccall libpetsc.DMLabelConvertToSection(arg1::DMLabel, arg2::Ptr{PetscSection}, arg3::Ptr{IS})::PetscErrorCode
end

function PetscSectionCreateGlobalSectionLabel(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscSectionCreateGlobalSectionLabel(arg1::PetscSection, arg2::PetscSF, arg3::PetscBool, arg4::DMLabel, arg5::PetscInt, arg6::Ptr{PetscSection})::PetscErrorCode
end

function PetscSectionSymCreateLabel(arg1, arg2, arg3)
    @ccall libpetsc.PetscSectionSymCreateLabel(arg1::MPI_Comm, arg2::DMLabel, arg3::Ptr{PetscSectionSym})::PetscErrorCode
end

function PetscSectionSymLabelSetLabel(arg1, arg2)
    @ccall libpetsc.PetscSectionSymLabelSetLabel(arg1::PetscSectionSym, arg2::DMLabel)::PetscErrorCode
end

function PetscSectionSymLabelSetStratum(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscSectionSymLabelSetStratum(arg1::PetscSectionSym, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscCopyMode, arg7::Ptr{Ptr{PetscInt}}, arg8::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

mutable struct _p_PetscDS end

const PetscDS = Ptr{_p_PetscDS}

mutable struct _p_PetscWeakForm end

const PetscWeakForm = Ptr{_p_PetscWeakForm}

mutable struct _PetscHashFormKey
    label::DMLabel
    value::PetscInt
    field::PetscInt
    _PetscHashFormKey() = new()
end

const PetscHashFormKey = _PetscHashFormKey

function DMInitializePackage()
    @ccall libpetsc.DMInitializePackage()::PetscErrorCode
end

const DMType = Ptr{Cchar}

function DMCreate(arg1, arg2)
    @ccall libpetsc.DMCreate(arg1::MPI_Comm, arg2::Ptr{DM})::PetscErrorCode
end

function DMClone(arg1, arg2)
    @ccall libpetsc.DMClone(arg1::DM, arg2::Ptr{DM})::PetscErrorCode
end

function DMSetType(arg1, arg2)
    @ccall libpetsc.DMSetType(arg1::DM, arg2::DMType)::PetscErrorCode
end

function DMGetType(arg1, arg2)
    @ccall libpetsc.DMGetType(arg1::DM, arg2::Ptr{DMType})::PetscErrorCode
end

function DMRegister(arg1, arg2)
    @ccall libpetsc.DMRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMRegisterDestroy()
    @ccall libpetsc.DMRegisterDestroy()::PetscErrorCode
end

function DMView(arg1, arg2)
    @ccall libpetsc.DMView(arg1::DM, arg2::PetscViewer)::PetscErrorCode
end

function DMLoad(arg1, arg2)
    @ccall libpetsc.DMLoad(arg1::DM, arg2::PetscViewer)::PetscErrorCode
end

function DMDestroy(arg1)
    @ccall libpetsc.DMDestroy(arg1::Ptr{DM})::PetscErrorCode
end

function DMCreateGlobalVector(arg1, arg2)
    @ccall libpetsc.DMCreateGlobalVector(arg1::DM, arg2::Ptr{Vec})::PetscErrorCode
end

function DMCreateLocalVector(arg1, arg2)
    @ccall libpetsc.DMCreateLocalVector(arg1::DM, arg2::Ptr{Vec})::PetscErrorCode
end

function DMGetLocalVector(arg1, arg2)
    @ccall libpetsc.DMGetLocalVector(arg1::DM, arg2::Ptr{Vec})::PetscErrorCode
end

function DMRestoreLocalVector(arg1, arg2)
    @ccall libpetsc.DMRestoreLocalVector(arg1::DM, arg2::Ptr{Vec})::PetscErrorCode
end

function DMGetGlobalVector(arg1, arg2)
    @ccall libpetsc.DMGetGlobalVector(arg1::DM, arg2::Ptr{Vec})::PetscErrorCode
end

function DMRestoreGlobalVector(arg1, arg2)
    @ccall libpetsc.DMRestoreGlobalVector(arg1::DM, arg2::Ptr{Vec})::PetscErrorCode
end

function DMClearGlobalVectors(arg1)
    @ccall libpetsc.DMClearGlobalVectors(arg1::DM)::PetscErrorCode
end

function DMClearLocalVectors(arg1)
    @ccall libpetsc.DMClearLocalVectors(arg1::DM)::PetscErrorCode
end

function DMHasNamedGlobalVector(arg1, arg2, arg3)
    @ccall libpetsc.DMHasNamedGlobalVector(arg1::DM, arg2::Ptr{Cchar}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function DMGetNamedGlobalVector(arg1, arg2, arg3)
    @ccall libpetsc.DMGetNamedGlobalVector(arg1::DM, arg2::Ptr{Cchar}, arg3::Ptr{Vec})::PetscErrorCode
end

function DMRestoreNamedGlobalVector(arg1, arg2, arg3)
    @ccall libpetsc.DMRestoreNamedGlobalVector(arg1::DM, arg2::Ptr{Cchar}, arg3::Ptr{Vec})::PetscErrorCode
end

function DMHasNamedLocalVector(arg1, arg2, arg3)
    @ccall libpetsc.DMHasNamedLocalVector(arg1::DM, arg2::Ptr{Cchar}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function DMGetNamedLocalVector(arg1, arg2, arg3)
    @ccall libpetsc.DMGetNamedLocalVector(arg1::DM, arg2::Ptr{Cchar}, arg3::Ptr{Vec})::PetscErrorCode
end

function DMRestoreNamedLocalVector(arg1, arg2, arg3)
    @ccall libpetsc.DMRestoreNamedLocalVector(arg1::DM, arg2::Ptr{Cchar}, arg3::Ptr{Vec})::PetscErrorCode
end

function DMGetLocalToGlobalMapping(arg1, arg2)
    @ccall libpetsc.DMGetLocalToGlobalMapping(arg1::DM, arg2::Ptr{ISLocalToGlobalMapping})::PetscErrorCode
end

function DMCreateFieldIS(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMCreateFieldIS(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{Ptr{Cchar}}}, arg4::Ptr{Ptr{IS}})::PetscErrorCode
end

function DMGetBlockSize(arg1, arg2)
    @ccall libpetsc.DMGetBlockSize(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMCreateColoring(arg1, arg2, arg3)
    @ccall libpetsc.DMCreateColoring(arg1::DM, arg2::ISColoringType, arg3::Ptr{ISColoring})::PetscErrorCode
end

function DMCreateMatrix(arg1, arg2)
    @ccall libpetsc.DMCreateMatrix(arg1::DM, arg2::Ptr{Mat})::PetscErrorCode
end

function DMSetMatrixPreallocateOnly(arg1, arg2)
    @ccall libpetsc.DMSetMatrixPreallocateOnly(arg1::DM, arg2::PetscBool)::PetscErrorCode
end

function DMSetMatrixStructureOnly(arg1, arg2)
    @ccall libpetsc.DMSetMatrixStructureOnly(arg1::DM, arg2::PetscBool)::PetscErrorCode
end

function DMCreateInterpolation(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMCreateInterpolation(arg1::DM, arg2::DM, arg3::Ptr{Mat}, arg4::Ptr{Vec})::PetscErrorCode
end

function DMCreateRestriction(arg1, arg2, arg3)
    @ccall libpetsc.DMCreateRestriction(arg1::DM, arg2::DM, arg3::Ptr{Mat})::PetscErrorCode
end

function DMCreateInjection(arg1, arg2, arg3)
    @ccall libpetsc.DMCreateInjection(arg1::DM, arg2::DM, arg3::Ptr{Mat})::PetscErrorCode
end

function DMCreateMassMatrix(arg1, arg2, arg3)
    @ccall libpetsc.DMCreateMassMatrix(arg1::DM, arg2::DM, arg3::Ptr{Mat})::PetscErrorCode
end

function DMGetWorkArray(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMGetWorkArray(arg1::DM, arg2::PetscInt, arg3::MPI_Datatype, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMRestoreWorkArray(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMRestoreWorkArray(arg1::DM, arg2::PetscInt, arg3::MPI_Datatype, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMRefine(arg1, arg2, arg3)
    @ccall libpetsc.DMRefine(arg1::DM, arg2::MPI_Comm, arg3::Ptr{DM})::PetscErrorCode
end

function DMCoarsen(arg1, arg2, arg3)
    @ccall libpetsc.DMCoarsen(arg1::DM, arg2::MPI_Comm, arg3::Ptr{DM})::PetscErrorCode
end

function DMGetCoarseDM(arg1, arg2)
    @ccall libpetsc.DMGetCoarseDM(arg1::DM, arg2::Ptr{DM})::PetscErrorCode
end

function DMSetCoarseDM(arg1, arg2)
    @ccall libpetsc.DMSetCoarseDM(arg1::DM, arg2::DM)::PetscErrorCode
end

function DMGetFineDM(arg1, arg2)
    @ccall libpetsc.DMGetFineDM(arg1::DM, arg2::Ptr{DM})::PetscErrorCode
end

function DMSetFineDM(arg1, arg2)
    @ccall libpetsc.DMSetFineDM(arg1::DM, arg2::DM)::PetscErrorCode
end

function DMRefineHierarchy(arg1, arg2, arg3)
    @ccall libpetsc.DMRefineHierarchy(arg1::DM, arg2::PetscInt, arg3::Ptr{DM})::PetscErrorCode
end

function DMCoarsenHierarchy(arg1, arg2, arg3)
    @ccall libpetsc.DMCoarsenHierarchy(arg1::DM, arg2::PetscInt, arg3::Ptr{DM})::PetscErrorCode
end

function DMCoarsenHookAdd(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMCoarsenHookAdd(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMCoarsenHookRemove(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMCoarsenHookRemove(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMRefineHookAdd(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMRefineHookAdd(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMRefineHookRemove(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMRefineHookRemove(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMRestrict(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMRestrict(arg1::DM, arg2::Mat, arg3::Vec, arg4::Mat, arg5::DM)::PetscErrorCode
end

function DMInterpolate(arg1, arg2, arg3)
    @ccall libpetsc.DMInterpolate(arg1::DM, arg2::Mat, arg3::DM)::PetscErrorCode
end

function DMInterpolateSolution(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMInterpolateSolution(arg1::DM, arg2::DM, arg3::Mat, arg4::Vec, arg5::Vec)::PetscErrorCode
end

function DMSetFromOptions(arg1)
    @ccall libpetsc.DMSetFromOptions(arg1::DM)::PetscErrorCode
end

function DMViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.DMViewFromOptions(arg1::DM, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function DMAdaptLabel(arg1, arg2, arg3)
    @ccall libpetsc.DMAdaptLabel(arg1::DM, arg2::DMLabel, arg3::Ptr{DM})::PetscErrorCode
end

function DMAdaptMetric(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMAdaptMetric(arg1::DM, arg2::Vec, arg3::DMLabel, arg4::Ptr{DM})::PetscErrorCode
end

function DMSetUp(arg1)
    @ccall libpetsc.DMSetUp(arg1::DM)::PetscErrorCode
end

function DMCreateInterpolationScale(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMCreateInterpolationScale(arg1::DM, arg2::DM, arg3::Mat, arg4::Ptr{Vec})::PetscErrorCode
end

function DMCreateAggregates(arg1, arg2, arg3)
    @ccall libpetsc.DMCreateAggregates(arg1::DM, arg2::DM, arg3::Ptr{Mat})::PetscErrorCode
end

function DMGlobalToLocalHookAdd(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMGlobalToLocalHookAdd(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMLocalToGlobalHookAdd(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMLocalToGlobalHookAdd(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMGlobalToLocal(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMGlobalToLocal(arg1::DM, arg2::Vec, arg3::InsertMode, arg4::Vec)::PetscErrorCode
end

function DMGlobalToLocalBegin(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMGlobalToLocalBegin(arg1::DM, arg2::Vec, arg3::InsertMode, arg4::Vec)::PetscErrorCode
end

function DMGlobalToLocalEnd(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMGlobalToLocalEnd(arg1::DM, arg2::Vec, arg3::InsertMode, arg4::Vec)::PetscErrorCode
end

function DMLocalToGlobal(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMLocalToGlobal(arg1::DM, arg2::Vec, arg3::InsertMode, arg4::Vec)::PetscErrorCode
end

function DMLocalToGlobalBegin(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMLocalToGlobalBegin(arg1::DM, arg2::Vec, arg3::InsertMode, arg4::Vec)::PetscErrorCode
end

function DMLocalToGlobalEnd(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMLocalToGlobalEnd(arg1::DM, arg2::Vec, arg3::InsertMode, arg4::Vec)::PetscErrorCode
end

function DMLocalToLocalBegin(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMLocalToLocalBegin(arg1::DM, arg2::Vec, arg3::InsertMode, arg4::Vec)::PetscErrorCode
end

function DMLocalToLocalEnd(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMLocalToLocalEnd(arg1::DM, arg2::Vec, arg3::InsertMode, arg4::Vec)::PetscErrorCode
end

function DMConvert(arg1, arg2, arg3)
    @ccall libpetsc.DMConvert(arg1::DM, arg2::DMType, arg3::Ptr{DM})::PetscErrorCode
end

function DMGetDimension(arg1, arg2)
    @ccall libpetsc.DMGetDimension(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMSetDimension(arg1, arg2)
    @ccall libpetsc.DMSetDimension(arg1::DM, arg2::PetscInt)::PetscErrorCode
end

function DMGetDimPoints(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMGetDimPoints(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMGetUseNatural(arg1, arg2)
    @ccall libpetsc.DMGetUseNatural(arg1::DM, arg2::Ptr{PetscBool})::PetscErrorCode
end

function DMSetUseNatural(arg1, arg2)
    @ccall libpetsc.DMSetUseNatural(arg1::DM, arg2::PetscBool)::PetscErrorCode
end

function DMGetCoordinateDM(arg1, arg2)
    @ccall libpetsc.DMGetCoordinateDM(arg1::DM, arg2::Ptr{DM})::PetscErrorCode
end

function DMSetCoordinateDM(arg1, arg2)
    @ccall libpetsc.DMSetCoordinateDM(arg1::DM, arg2::DM)::PetscErrorCode
end

function DMGetCoordinateDim(arg1, arg2)
    @ccall libpetsc.DMGetCoordinateDim(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMSetCoordinateDim(arg1, arg2)
    @ccall libpetsc.DMSetCoordinateDim(arg1::DM, arg2::PetscInt)::PetscErrorCode
end

function DMGetCoordinateSection(arg1, arg2)
    @ccall libpetsc.DMGetCoordinateSection(arg1::DM, arg2::Ptr{PetscSection})::PetscErrorCode
end

function DMSetCoordinateSection(arg1, arg2, arg3)
    @ccall libpetsc.DMSetCoordinateSection(arg1::DM, arg2::PetscInt, arg3::PetscSection)::PetscErrorCode
end

function DMGetCoordinates(arg1, arg2)
    @ccall libpetsc.DMGetCoordinates(arg1::DM, arg2::Ptr{Vec})::PetscErrorCode
end

function DMSetCoordinates(arg1, arg2)
    @ccall libpetsc.DMSetCoordinates(arg1::DM, arg2::Vec)::PetscErrorCode
end

function DMGetCoordinatesLocal(arg1, arg2)
    @ccall libpetsc.DMGetCoordinatesLocal(arg1::DM, arg2::Ptr{Vec})::PetscErrorCode
end

function DMGetCoordinatesLocalSetUp(arg1)
    @ccall libpetsc.DMGetCoordinatesLocalSetUp(arg1::DM)::PetscErrorCode
end

function DMGetCoordinatesLocalNoncollective(arg1, arg2)
    @ccall libpetsc.DMGetCoordinatesLocalNoncollective(arg1::DM, arg2::Ptr{Vec})::PetscErrorCode
end

function DMGetCoordinatesLocalTuple(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMGetCoordinatesLocalTuple(arg1::DM, arg2::IS, arg3::Ptr{PetscSection}, arg4::Ptr{Vec})::PetscErrorCode
end

function DMSetCoordinatesLocal(arg1, arg2)
    @ccall libpetsc.DMSetCoordinatesLocal(arg1::DM, arg2::Vec)::PetscErrorCode
end

function DMLocatePoints(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMLocatePoints(arg1::DM, arg2::Vec, arg3::DMPointLocationType, arg4::Ptr{PetscSF})::PetscErrorCode
end

function DMGetPeriodicity(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMGetPeriodicity(arg1::DM, arg2::Ptr{PetscBool}, arg3::Ptr{Ptr{PetscReal}}, arg4::Ptr{Ptr{PetscReal}}, arg5::Ptr{Ptr{DMBoundaryType}})::PetscErrorCode
end

function DMSetPeriodicity(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMSetPeriodicity(arg1::DM, arg2::PetscBool, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{DMBoundaryType})::PetscErrorCode
end

function DMLocalizeCoordinate(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMLocalizeCoordinate(arg1::DM, arg2::Ptr{PetscScalar}, arg3::PetscBool, arg4::Ptr{PetscScalar})::PetscErrorCode
end

function DMLocalizeCoordinates(arg1)
    @ccall libpetsc.DMLocalizeCoordinates(arg1::DM)::PetscErrorCode
end

function DMGetCoordinatesLocalized(arg1, arg2)
    @ccall libpetsc.DMGetCoordinatesLocalized(arg1::DM, arg2::Ptr{PetscBool})::PetscErrorCode
end

function DMGetCoordinatesLocalizedLocal(arg1, arg2)
    @ccall libpetsc.DMGetCoordinatesLocalizedLocal(arg1::DM, arg2::Ptr{PetscBool})::PetscErrorCode
end

function DMGetNeighbors(arg1, arg2, arg3)
    @ccall libpetsc.DMGetNeighbors(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{PetscMPIInt}})::PetscErrorCode
end

function DMGetCoordinateField(arg1, arg2)
    @ccall libpetsc.DMGetCoordinateField(arg1::DM, arg2::Ptr{DMField})::PetscErrorCode
end

function DMSetCoordinateField(arg1, arg2)
    @ccall libpetsc.DMSetCoordinateField(arg1::DM, arg2::DMField)::PetscErrorCode
end

function DMGetBoundingBox(arg1, arg2, arg3)
    @ccall libpetsc.DMGetBoundingBox(arg1::DM, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal})::PetscErrorCode
end

function DMGetLocalBoundingBox(arg1, arg2, arg3)
    @ccall libpetsc.DMGetLocalBoundingBox(arg1::DM, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal})::PetscErrorCode
end

function DMProjectCoordinates(arg1, arg2)
    @ccall libpetsc.DMProjectCoordinates(arg1::DM, arg2::PetscFE)::PetscErrorCode
end

function DMSubDomainHookAdd(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMSubDomainHookAdd(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMSubDomainHookRemove(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMSubDomainHookRemove(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMSubDomainRestrict(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMSubDomainRestrict(arg1::DM, arg2::VecScatter, arg3::VecScatter, arg4::DM)::PetscErrorCode
end

function DMSetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.DMSetOptionsPrefix(arg1::DM, arg2::Ptr{Cchar})::PetscErrorCode
end

function DMAppendOptionsPrefix(arg1, arg2)
    @ccall libpetsc.DMAppendOptionsPrefix(arg1::DM, arg2::Ptr{Cchar})::PetscErrorCode
end

function DMGetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.DMGetOptionsPrefix(arg1::DM, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function DMSetVecType(arg1, arg2)
    @ccall libpetsc.DMSetVecType(arg1::DM, arg2::VecType)::PetscErrorCode
end

function DMGetVecType(arg1, arg2)
    @ccall libpetsc.DMGetVecType(arg1::DM, arg2::Ptr{VecType})::PetscErrorCode
end

function DMSetMatType(arg1, arg2)
    @ccall libpetsc.DMSetMatType(arg1::DM, arg2::MatType)::PetscErrorCode
end

function DMGetMatType(arg1, arg2)
    @ccall libpetsc.DMGetMatType(arg1::DM, arg2::Ptr{MatType})::PetscErrorCode
end

function DMSetISColoringType(arg1, arg2)
    @ccall libpetsc.DMSetISColoringType(arg1::DM, arg2::ISColoringType)::PetscErrorCode
end

function DMGetISColoringType(arg1, arg2)
    @ccall libpetsc.DMGetISColoringType(arg1::DM, arg2::Ptr{ISColoringType})::PetscErrorCode
end

function DMSetApplicationContext(arg1, arg2)
    @ccall libpetsc.DMSetApplicationContext(arg1::DM, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMSetApplicationContextDestroy(arg1, arg2)
    @ccall libpetsc.DMSetApplicationContextDestroy(arg1::DM, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMGetApplicationContext(arg1, arg2)
    @ccall libpetsc.DMGetApplicationContext(arg1::DM, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMSetVariableBounds(arg1, arg2)
    @ccall libpetsc.DMSetVariableBounds(arg1::DM, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMHasVariableBounds(arg1, arg2)
    @ccall libpetsc.DMHasVariableBounds(arg1::DM, arg2::Ptr{PetscBool})::PetscErrorCode
end

function DMHasColoring(arg1, arg2)
    @ccall libpetsc.DMHasColoring(arg1::DM, arg2::Ptr{PetscBool})::PetscErrorCode
end

function DMHasCreateRestriction(arg1, arg2)
    @ccall libpetsc.DMHasCreateRestriction(arg1::DM, arg2::Ptr{PetscBool})::PetscErrorCode
end

function DMHasCreateInjection(arg1, arg2)
    @ccall libpetsc.DMHasCreateInjection(arg1::DM, arg2::Ptr{PetscBool})::PetscErrorCode
end

function DMComputeVariableBounds(arg1, arg2, arg3)
    @ccall libpetsc.DMComputeVariableBounds(arg1::DM, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function DMCreateSubDM(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMCreateSubDM(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{IS}, arg5::Ptr{DM})::PetscErrorCode
end

function DMCreateSuperDM(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMCreateSuperDM(arg1::Ptr{DM}, arg2::PetscInt, arg3::Ptr{Ptr{IS}}, arg4::Ptr{DM})::PetscErrorCode
end

function DMCreateSectionSubDM(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMCreateSectionSubDM(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{IS}, arg5::Ptr{DM})::PetscErrorCode
end

function DMCreateSectionSuperDM(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMCreateSectionSuperDM(arg1::Ptr{DM}, arg2::PetscInt, arg3::Ptr{Ptr{IS}}, arg4::Ptr{DM})::PetscErrorCode
end

function DMCreateFieldDecomposition(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMCreateFieldDecomposition(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{Ptr{Cchar}}}, arg4::Ptr{Ptr{IS}}, arg5::Ptr{Ptr{DM}})::PetscErrorCode
end

function DMCreateDomainDecomposition(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMCreateDomainDecomposition(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{Ptr{Cchar}}}, arg4::Ptr{Ptr{IS}}, arg5::Ptr{Ptr{IS}}, arg6::Ptr{Ptr{DM}})::PetscErrorCode
end

function DMCreateDomainDecompositionScatters(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMCreateDomainDecompositionScatters(arg1::DM, arg2::PetscInt, arg3::Ptr{DM}, arg4::Ptr{Ptr{VecScatter}}, arg5::Ptr{Ptr{VecScatter}}, arg6::Ptr{Ptr{VecScatter}})::PetscErrorCode
end

function DMGetRefineLevel(arg1, arg2)
    @ccall libpetsc.DMGetRefineLevel(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMSetRefineLevel(arg1, arg2)
    @ccall libpetsc.DMSetRefineLevel(arg1::DM, arg2::PetscInt)::PetscErrorCode
end

function DMGetCoarsenLevel(arg1, arg2)
    @ccall libpetsc.DMGetCoarsenLevel(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMSetCoarsenLevel(arg1, arg2)
    @ccall libpetsc.DMSetCoarsenLevel(arg1::DM, arg2::PetscInt)::PetscErrorCode
end

function DMFinalizePackage()
    @ccall libpetsc.DMFinalizePackage()::PetscErrorCode
end

function VecGetDM(arg1, arg2)
    @ccall libpetsc.VecGetDM(arg1::Vec, arg2::Ptr{DM})::PetscErrorCode
end

function VecSetDM(arg1, arg2)
    @ccall libpetsc.VecSetDM(arg1::Vec, arg2::DM)::PetscErrorCode
end

function MatGetDM(arg1, arg2)
    @ccall libpetsc.MatGetDM(arg1::Mat, arg2::Ptr{DM})::PetscErrorCode
end

function MatSetDM(arg1, arg2)
    @ccall libpetsc.MatSetDM(arg1::Mat, arg2::DM)::PetscErrorCode
end

function MatFDColoringUseDM(arg1, arg2)
    @ccall libpetsc.MatFDColoringUseDM(arg1::Mat, arg2::MatFDColoring)::PetscErrorCode
end

mutable struct NLF_DAAD end

const NLF = Ptr{NLF_DAAD}

function DMPrintCellVector(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPrintCellVector(arg1::PetscInt, arg2::Ptr{Cchar}, arg3::PetscInt, arg4::Ptr{PetscScalar})::PetscErrorCode
end

function DMPrintCellMatrix(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPrintCellMatrix(arg1::PetscInt, arg2::Ptr{Cchar}, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscScalar})::PetscErrorCode
end

function DMPrintLocalVec(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPrintLocalVec(arg1::DM, arg2::Ptr{Cchar}, arg3::PetscReal, arg4::Vec)::PetscErrorCode
end

function DMSetNullSpaceConstructor(arg1, arg2, arg3)
    @ccall libpetsc.DMSetNullSpaceConstructor(arg1::DM, arg2::PetscInt, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMGetNullSpaceConstructor(arg1, arg2, arg3)
    @ccall libpetsc.DMGetNullSpaceConstructor(arg1::DM, arg2::PetscInt, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMSetNearNullSpaceConstructor(arg1, arg2, arg3)
    @ccall libpetsc.DMSetNearNullSpaceConstructor(arg1::DM, arg2::PetscInt, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMGetNearNullSpaceConstructor(arg1, arg2, arg3)
    @ccall libpetsc.DMGetNearNullSpaceConstructor(arg1::DM, arg2::PetscInt, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMGetSection(arg1, arg2)
    @ccall libpetsc.DMGetSection(arg1::DM, arg2::Ptr{PetscSection})::PetscErrorCode
end

function DMSetSection(arg1, arg2)
    @ccall libpetsc.DMSetSection(arg1::DM, arg2::PetscSection)::PetscErrorCode
end

function DMGetLocalSection(arg1, arg2)
    @ccall libpetsc.DMGetLocalSection(arg1::DM, arg2::Ptr{PetscSection})::PetscErrorCode
end

function DMSetLocalSection(arg1, arg2)
    @ccall libpetsc.DMSetLocalSection(arg1::DM, arg2::PetscSection)::PetscErrorCode
end

function DMGetGlobalSection(arg1, arg2)
    @ccall libpetsc.DMGetGlobalSection(arg1::DM, arg2::Ptr{PetscSection})::PetscErrorCode
end

function DMSetGlobalSection(arg1, arg2)
    @ccall libpetsc.DMSetGlobalSection(arg1::DM, arg2::PetscSection)::PetscErrorCode
end

function DMGetDefaultSection(dm, s)
    @ccall libpetsc.DMGetDefaultSection(dm::DM, s::Ptr{PetscSection})::PetscErrorCode
end

function DMSetDefaultSection(dm, s)
    @ccall libpetsc.DMSetDefaultSection(dm::DM, s::PetscSection)::PetscErrorCode
end

function DMGetDefaultGlobalSection(dm, s)
    @ccall libpetsc.DMGetDefaultGlobalSection(dm::DM, s::Ptr{PetscSection})::PetscErrorCode
end

function DMSetDefaultGlobalSection(dm, s)
    @ccall libpetsc.DMSetDefaultGlobalSection(dm::DM, s::PetscSection)::PetscErrorCode
end

function DMGetSectionSF(arg1, arg2)
    @ccall libpetsc.DMGetSectionSF(arg1::DM, arg2::Ptr{PetscSF})::PetscErrorCode
end

function DMSetSectionSF(arg1, arg2)
    @ccall libpetsc.DMSetSectionSF(arg1::DM, arg2::PetscSF)::PetscErrorCode
end

function DMCreateSectionSF(arg1, arg2, arg3)
    @ccall libpetsc.DMCreateSectionSF(arg1::DM, arg2::PetscSection, arg3::PetscSection)::PetscErrorCode
end

function DMGetDefaultSF(dm, s)
    @ccall libpetsc.DMGetDefaultSF(dm::DM, s::Ptr{PetscSF})::PetscErrorCode
end

function DMSetDefaultSF(dm, s)
    @ccall libpetsc.DMSetDefaultSF(dm::DM, s::PetscSF)::PetscErrorCode
end

function DMCreateDefaultSF(dm, l, g)
    @ccall libpetsc.DMCreateDefaultSF(dm::DM, l::PetscSection, g::PetscSection)::PetscErrorCode
end

function DMGetPointSF(arg1, arg2)
    @ccall libpetsc.DMGetPointSF(arg1::DM, arg2::Ptr{PetscSF})::PetscErrorCode
end

function DMSetPointSF(arg1, arg2)
    @ccall libpetsc.DMSetPointSF(arg1::DM, arg2::PetscSF)::PetscErrorCode
end

function DMGetDefaultConstraints(arg1, arg2, arg3)
    @ccall libpetsc.DMGetDefaultConstraints(arg1::DM, arg2::Ptr{PetscSection}, arg3::Ptr{Mat})::PetscErrorCode
end

function DMSetDefaultConstraints(arg1, arg2, arg3)
    @ccall libpetsc.DMSetDefaultConstraints(arg1::DM, arg2::PetscSection, arg3::Mat)::PetscErrorCode
end

function DMGetOutputDM(arg1, arg2)
    @ccall libpetsc.DMGetOutputDM(arg1::DM, arg2::Ptr{DM})::PetscErrorCode
end

function DMGetOutputSequenceNumber(arg1, arg2, arg3)
    @ccall libpetsc.DMGetOutputSequenceNumber(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscReal})::PetscErrorCode
end

function DMSetOutputSequenceNumber(arg1, arg2, arg3)
    @ccall libpetsc.DMSetOutputSequenceNumber(arg1::DM, arg2::PetscInt, arg3::PetscReal)::PetscErrorCode
end

function DMOutputSequenceLoad(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMOutputSequenceLoad(arg1::DM, arg2::PetscViewer, arg3::Ptr{Cchar}, arg4::PetscInt, arg5::Ptr{PetscReal})::PetscErrorCode
end

function DMGetNumFields(arg1, arg2)
    @ccall libpetsc.DMGetNumFields(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMSetNumFields(arg1, arg2)
    @ccall libpetsc.DMSetNumFields(arg1::DM, arg2::PetscInt)::PetscErrorCode
end

function DMGetField(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMGetField(arg1::DM, arg2::PetscInt, arg3::Ptr{DMLabel}, arg4::Ptr{PetscObject})::PetscErrorCode
end

function DMSetField(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMSetField(arg1::DM, arg2::PetscInt, arg3::DMLabel, arg4::PetscObject)::PetscErrorCode
end

function DMAddField(arg1, arg2, arg3)
    @ccall libpetsc.DMAddField(arg1::DM, arg2::DMLabel, arg3::PetscObject)::PetscErrorCode
end

function DMSetFieldAvoidTensor(arg1, arg2, arg3)
    @ccall libpetsc.DMSetFieldAvoidTensor(arg1::DM, arg2::PetscInt, arg3::PetscBool)::PetscErrorCode
end

function DMGetFieldAvoidTensor(arg1, arg2, arg3)
    @ccall libpetsc.DMGetFieldAvoidTensor(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscBool})::PetscErrorCode
end

function DMClearFields(arg1)
    @ccall libpetsc.DMClearFields(arg1::DM)::PetscErrorCode
end

function DMCopyFields(arg1, arg2)
    @ccall libpetsc.DMCopyFields(arg1::DM, arg2::DM)::PetscErrorCode
end

function DMGetAdjacency(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMGetAdjacency(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscBool}, arg4::Ptr{PetscBool})::PetscErrorCode
end

function DMSetAdjacency(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMSetAdjacency(arg1::DM, arg2::PetscInt, arg3::PetscBool, arg4::PetscBool)::PetscErrorCode
end

function DMGetBasicAdjacency(arg1, arg2, arg3)
    @ccall libpetsc.DMGetBasicAdjacency(arg1::DM, arg2::Ptr{PetscBool}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function DMSetBasicAdjacency(arg1, arg2, arg3)
    @ccall libpetsc.DMSetBasicAdjacency(arg1::DM, arg2::PetscBool, arg3::PetscBool)::PetscErrorCode
end

function DMGetNumDS(arg1, arg2)
    @ccall libpetsc.DMGetNumDS(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMGetDS(arg1, arg2)
    @ccall libpetsc.DMGetDS(arg1::DM, arg2::Ptr{PetscDS})::PetscErrorCode
end

function DMGetCellDS(arg1, arg2, arg3)
    @ccall libpetsc.DMGetCellDS(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscDS})::PetscErrorCode
end

function DMGetRegionDS(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMGetRegionDS(arg1::DM, arg2::DMLabel, arg3::Ptr{IS}, arg4::Ptr{PetscDS})::PetscErrorCode
end

function DMSetRegionDS(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMSetRegionDS(arg1::DM, arg2::DMLabel, arg3::IS, arg4::PetscDS)::PetscErrorCode
end

function DMGetRegionNumDS(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMGetRegionNumDS(arg1::DM, arg2::PetscInt, arg3::Ptr{DMLabel}, arg4::Ptr{IS}, arg5::Ptr{PetscDS})::PetscErrorCode
end

function DMSetRegionNumDS(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMSetRegionNumDS(arg1::DM, arg2::PetscInt, arg3::DMLabel, arg4::IS, arg5::PetscDS)::PetscErrorCode
end

function DMFindRegionNum(arg1, arg2, arg3)
    @ccall libpetsc.DMFindRegionNum(arg1::DM, arg2::PetscDS, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMCreateDS(arg1)
    @ccall libpetsc.DMCreateDS(arg1::DM)::PetscErrorCode
end

function DMClearDS(arg1)
    @ccall libpetsc.DMClearDS(arg1::DM)::PetscErrorCode
end

function DMCopyDS(arg1, arg2)
    @ccall libpetsc.DMCopyDS(arg1::DM, arg2::DM)::PetscErrorCode
end

function DMCopyDisc(arg1, arg2)
    @ccall libpetsc.DMCopyDisc(arg1::DM, arg2::DM)::PetscErrorCode
end

function DMComputeExactSolution(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMComputeExactSolution(arg1::DM, arg2::PetscReal, arg3::Vec, arg4::Vec)::PetscErrorCode
end

mutable struct _DMInterpolationInfo
    comm::MPI_Comm
    dim::PetscInt
    nInput::PetscInt
    points::Ptr{PetscReal}
    cells::Ptr{PetscInt}
    n::PetscInt
    coords::Vec
    dof::PetscInt
    _DMInterpolationInfo() = new()
end

const DMInterpolationInfo = Ptr{_DMInterpolationInfo}

function DMInterpolationCreate(arg1, arg2)
    @ccall libpetsc.DMInterpolationCreate(arg1::MPI_Comm, arg2::Ptr{DMInterpolationInfo})::PetscErrorCode
end

function DMInterpolationSetDim(arg1, arg2)
    @ccall libpetsc.DMInterpolationSetDim(arg1::DMInterpolationInfo, arg2::PetscInt)::PetscErrorCode
end

function DMInterpolationGetDim(arg1, arg2)
    @ccall libpetsc.DMInterpolationGetDim(arg1::DMInterpolationInfo, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMInterpolationSetDof(arg1, arg2)
    @ccall libpetsc.DMInterpolationSetDof(arg1::DMInterpolationInfo, arg2::PetscInt)::PetscErrorCode
end

function DMInterpolationGetDof(arg1, arg2)
    @ccall libpetsc.DMInterpolationGetDof(arg1::DMInterpolationInfo, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMInterpolationAddPoints(arg1, arg2, arg3)
    @ccall libpetsc.DMInterpolationAddPoints(arg1::DMInterpolationInfo, arg2::PetscInt, arg3::Ptr{PetscReal})::PetscErrorCode
end

function DMInterpolationSetUp(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMInterpolationSetUp(arg1::DMInterpolationInfo, arg2::DM, arg3::PetscBool, arg4::PetscBool)::PetscErrorCode
end

function DMInterpolationGetCoordinates(arg1, arg2)
    @ccall libpetsc.DMInterpolationGetCoordinates(arg1::DMInterpolationInfo, arg2::Ptr{Vec})::PetscErrorCode
end

function DMInterpolationGetVector(arg1, arg2)
    @ccall libpetsc.DMInterpolationGetVector(arg1::DMInterpolationInfo, arg2::Ptr{Vec})::PetscErrorCode
end

function DMInterpolationRestoreVector(arg1, arg2)
    @ccall libpetsc.DMInterpolationRestoreVector(arg1::DMInterpolationInfo, arg2::Ptr{Vec})::PetscErrorCode
end

function DMInterpolationEvaluate(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMInterpolationEvaluate(arg1::DMInterpolationInfo, arg2::DM, arg3::Vec, arg4::Vec)::PetscErrorCode
end

function DMInterpolationDestroy(arg1)
    @ccall libpetsc.DMInterpolationDestroy(arg1::Ptr{DMInterpolationInfo})::PetscErrorCode
end

function DMCreateLabel(arg1, arg2)
    @ccall libpetsc.DMCreateLabel(arg1::DM, arg2::Ptr{Cchar})::PetscErrorCode
end

function DMGetLabelValue(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMGetLabelValue(arg1::DM, arg2::Ptr{Cchar}, arg3::PetscInt, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMSetLabelValue(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMSetLabelValue(arg1::DM, arg2::Ptr{Cchar}, arg3::PetscInt, arg4::PetscInt)::PetscErrorCode
end

function DMClearLabelValue(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMClearLabelValue(arg1::DM, arg2::Ptr{Cchar}, arg3::PetscInt, arg4::PetscInt)::PetscErrorCode
end

function DMGetLabelSize(arg1, arg2, arg3)
    @ccall libpetsc.DMGetLabelSize(arg1::DM, arg2::Ptr{Cchar}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMGetLabelIdIS(arg1, arg2, arg3)
    @ccall libpetsc.DMGetLabelIdIS(arg1::DM, arg2::Ptr{Cchar}, arg3::Ptr{IS})::PetscErrorCode
end

function DMGetStratumSize(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMGetStratumSize(arg1::DM, arg2::Ptr{Cchar}, arg3::PetscInt, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMGetStratumIS(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMGetStratumIS(arg1::DM, arg2::Ptr{Cchar}, arg3::PetscInt, arg4::Ptr{IS})::PetscErrorCode
end

function DMSetStratumIS(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMSetStratumIS(arg1::DM, arg2::Ptr{Cchar}, arg3::PetscInt, arg4::IS)::PetscErrorCode
end

function DMClearLabelStratum(arg1, arg2, arg3)
    @ccall libpetsc.DMClearLabelStratum(arg1::DM, arg2::Ptr{Cchar}, arg3::PetscInt)::PetscErrorCode
end

function DMGetLabelOutput(arg1, arg2, arg3)
    @ccall libpetsc.DMGetLabelOutput(arg1::DM, arg2::Ptr{Cchar}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function DMSetLabelOutput(arg1, arg2, arg3)
    @ccall libpetsc.DMSetLabelOutput(arg1::DM, arg2::Ptr{Cchar}, arg3::PetscBool)::PetscErrorCode
end

function DMGetNumLabels(arg1, arg2)
    @ccall libpetsc.DMGetNumLabels(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMGetLabelName(arg1, arg2, arg3)
    @ccall libpetsc.DMGetLabelName(arg1::DM, arg2::PetscInt, arg3::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function DMHasLabel(arg1, arg2, arg3)
    @ccall libpetsc.DMHasLabel(arg1::DM, arg2::Ptr{Cchar}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function DMGetLabel(arg1, arg2, arg3)
    @ccall libpetsc.DMGetLabel(arg1::DM, arg2::Ptr{Cchar}, arg3::Ptr{DMLabel})::PetscErrorCode
end

function DMGetLabelByNum(arg1, arg2, arg3)
    @ccall libpetsc.DMGetLabelByNum(arg1::DM, arg2::PetscInt, arg3::Ptr{DMLabel})::PetscErrorCode
end

function DMAddLabel(arg1, arg2)
    @ccall libpetsc.DMAddLabel(arg1::DM, arg2::DMLabel)::PetscErrorCode
end

function DMRemoveLabel(arg1, arg2, arg3)
    @ccall libpetsc.DMRemoveLabel(arg1::DM, arg2::Ptr{Cchar}, arg3::Ptr{DMLabel})::PetscErrorCode
end

function DMRemoveLabelBySelf(arg1, arg2, arg3)
    @ccall libpetsc.DMRemoveLabelBySelf(arg1::DM, arg2::Ptr{DMLabel}, arg3::PetscBool)::PetscErrorCode
end

function DMCopyLabels(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMCopyLabels(arg1::DM, arg2::DM, arg3::PetscCopyMode, arg4::PetscBool)::PetscErrorCode
end

function DMAddBoundary(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12)
    @ccall libpetsc.DMAddBoundary(arg1::DM, arg2::DMBoundaryConditionType, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{PetscInt}, arg8::Ptr{Cvoid}, arg9::Ptr{Cvoid}, arg10::PetscInt, arg11::Ptr{PetscInt}, arg12::Ptr{Cvoid})::PetscErrorCode
end

function DMGetNumBoundary(arg1, arg2)
    @ccall libpetsc.DMGetNumBoundary(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMGetBoundary(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13)
    @ccall libpetsc.DMGetBoundary(arg1::DM, arg2::PetscInt, arg3::Ptr{DMBoundaryConditionType}, arg4::Ptr{Ptr{Cchar}}, arg5::Ptr{Ptr{Cchar}}, arg6::Ptr{PetscInt}, arg7::Ptr{PetscInt}, arg8::Ptr{Ptr{PetscInt}}, arg9::Ptr{Ptr{Cvoid}}, arg10::Ptr{Ptr{Cvoid}}, arg11::Ptr{PetscInt}, arg12::Ptr{Ptr{PetscInt}}, arg13::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMIsBoundaryPoint(arg1, arg2, arg3)
    @ccall libpetsc.DMIsBoundaryPoint(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscBool})::PetscErrorCode
end

function DMCopyBoundary(arg1, arg2)
    @ccall libpetsc.DMCopyBoundary(arg1::DM, arg2::DM)::PetscErrorCode
end

function DMProjectFunction(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMProjectFunction(arg1::DM, arg2::PetscReal, arg3::Ptr{Ptr{Cvoid}}, arg4::Ptr{Ptr{Cvoid}}, arg5::InsertMode, arg6::Vec)::PetscErrorCode
end

function DMProjectFunctionLocal(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMProjectFunctionLocal(arg1::DM, arg2::PetscReal, arg3::Ptr{Ptr{Cvoid}}, arg4::Ptr{Ptr{Cvoid}}, arg5::InsertMode, arg6::Vec)::PetscErrorCode
end

function DMProjectFunctionLabel(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11)
    @ccall libpetsc.DMProjectFunctionLabel(arg1::DM, arg2::PetscReal, arg3::DMLabel, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::PetscInt, arg7::Ptr{PetscInt}, arg8::Ptr{Ptr{Cvoid}}, arg9::Ptr{Ptr{Cvoid}}, arg10::InsertMode, arg11::Vec)::PetscErrorCode
end

function DMProjectFunctionLabelLocal(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11)
    @ccall libpetsc.DMProjectFunctionLabelLocal(arg1::DM, arg2::PetscReal, arg3::DMLabel, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::PetscInt, arg7::Ptr{PetscInt}, arg8::Ptr{Ptr{Cvoid}}, arg9::Ptr{Ptr{Cvoid}}, arg10::InsertMode, arg11::Vec)::PetscErrorCode
end

function DMProjectFieldLocal(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMProjectFieldLocal(arg1::DM, arg2::PetscReal, arg3::Vec, arg4::Ptr{Ptr{Cvoid}}, arg5::InsertMode, arg6::Vec)::PetscErrorCode
end

function DMProjectFieldLabelLocal(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11)
    @ccall libpetsc.DMProjectFieldLabelLocal(arg1::DM, arg2::PetscReal, arg3::DMLabel, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::PetscInt, arg7::Ptr{PetscInt}, arg8::Vec, arg9::Ptr{Ptr{Cvoid}}, arg10::InsertMode, arg11::Vec)::PetscErrorCode
end

function DMProjectBdFieldLabelLocal(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11)
    @ccall libpetsc.DMProjectBdFieldLabelLocal(arg1::DM, arg2::PetscReal, arg3::DMLabel, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::PetscInt, arg7::Ptr{PetscInt}, arg8::Vec, arg9::Ptr{Ptr{Cvoid}}, arg10::InsertMode, arg11::Vec)::PetscErrorCode
end

function DMComputeL2Diff(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMComputeL2Diff(arg1::DM, arg2::PetscReal, arg3::Ptr{Ptr{Cvoid}}, arg4::Ptr{Ptr{Cvoid}}, arg5::Vec, arg6::Ptr{PetscReal})::PetscErrorCode
end

function DMComputeL2GradientDiff(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMComputeL2GradientDiff(arg1::DM, arg2::PetscReal, arg3::Ptr{Ptr{Cvoid}}, arg4::Ptr{Ptr{Cvoid}}, arg5::Vec, arg6::Ptr{PetscReal}, arg7::Ptr{PetscReal})::PetscErrorCode
end

function DMComputeL2FieldDiff(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMComputeL2FieldDiff(arg1::DM, arg2::PetscReal, arg3::Ptr{Ptr{Cvoid}}, arg4::Ptr{Ptr{Cvoid}}, arg5::Vec, arg6::Ptr{PetscReal})::PetscErrorCode
end

function DMComputeError(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMComputeError(arg1::DM, arg2::Vec, arg3::Ptr{PetscReal}, arg4::Ptr{Vec})::PetscErrorCode
end

function DMHasBasisTransform(arg1, arg2)
    @ccall libpetsc.DMHasBasisTransform(arg1::DM, arg2::Ptr{PetscBool})::PetscErrorCode
end

function DMCopyTransform(arg1, arg2)
    @ccall libpetsc.DMCopyTransform(arg1::DM, arg2::DM)::PetscErrorCode
end

function DMGetCompatibility(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMGetCompatibility(arg1::DM, arg2::DM, arg3::Ptr{PetscBool}, arg4::Ptr{PetscBool})::PetscErrorCode
end

function DMMonitorSet(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMMonitorSet(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMMonitorCancel(arg1)
    @ccall libpetsc.DMMonitorCancel(arg1::DM)::PetscErrorCode
end

function DMMonitorSetFromOptions(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMMonitorSetFromOptions(arg1::DM, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{Cvoid}, arg6::Ptr{Cvoid}, arg7::Ptr{PetscBool})::PetscErrorCode
end

function DMMonitor(arg1)
    @ccall libpetsc.DMMonitor(arg1::DM)::PetscErrorCode
end

function DMPolytopeTypeGetDim(ct)
    @ccall libpetsc.DMPolytopeTypeGetDim(ct::DMPolytopeType)::PetscInt
end

function DMPolytopeTypeGetConeSize(ct)
    @ccall libpetsc.DMPolytopeTypeGetConeSize(ct::DMPolytopeType)::PetscInt
end

function DMPolytopeTypeGetNumVertices(ct)
    @ccall libpetsc.DMPolytopeTypeGetNumVertices(ct::DMPolytopeType)::PetscInt
end

@enum DMDAStencilType::UInt32 begin
    DMDA_STENCIL_STAR = 0
    DMDA_STENCIL_BOX = 1
end

@enum DMDAInterpolationType::UInt32 begin
    DMDA_Q0 = 0
    DMDA_Q1 = 1
end

@enum DMDAElementType::UInt32 begin
    DMDA_ELEMENT_P1 = 0
    DMDA_ELEMENT_Q1 = 1
end

mutable struct DMDALocalInfo
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
    DMDALocalInfo() = new()
end

const PFType = Ptr{Cchar}

mutable struct _p_PF end

const PF = Ptr{_p_PF}

function PFCreate(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PFCreate(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PF})::PetscErrorCode
end

function PFSetType(arg1, arg2, arg3)
    @ccall libpetsc.PFSetType(arg1::PF, arg2::PFType, arg3::Ptr{Cvoid})::PetscErrorCode
end

function PFSet(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PFSet(arg1::PF, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid}, arg6::Ptr{Cvoid})::PetscErrorCode
end

function PFApply(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PFApply(arg1::PF, arg2::PetscInt, arg3::Ptr{PetscScalar}, arg4::Ptr{PetscScalar})::PetscErrorCode
end

function PFApplyVec(arg1, arg2, arg3)
    @ccall libpetsc.PFApplyVec(arg1::PF, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function PFInitializePackage()
    @ccall libpetsc.PFInitializePackage()::PetscErrorCode
end

function PFRegister(arg1, arg2)
    @ccall libpetsc.PFRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PFDestroy(arg1)
    @ccall libpetsc.PFDestroy(arg1::Ptr{PF})::PetscErrorCode
end

function PFSetFromOptions(arg1)
    @ccall libpetsc.PFSetFromOptions(arg1::PF)::PetscErrorCode
end

function PFGetType(arg1, arg2)
    @ccall libpetsc.PFGetType(arg1::PF, arg2::Ptr{PFType})::PetscErrorCode
end

function PFView(arg1, arg2)
    @ccall libpetsc.PFView(arg1::PF, arg2::PetscViewer)::PetscErrorCode
end

function PFViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.PFViewFromOptions(arg1::PF, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

mutable struct _p_AO end

const AO = Ptr{_p_AO}

const AOType = Ptr{Cchar}

function AOInitializePackage()
    @ccall libpetsc.AOInitializePackage()::PetscErrorCode
end

function AOCreate(arg1, arg2)
    @ccall libpetsc.AOCreate(arg1::MPI_Comm, arg2::Ptr{AO})::PetscErrorCode
end

function AOSetIS(arg1, arg2, arg3)
    @ccall libpetsc.AOSetIS(arg1::AO, arg2::IS, arg3::IS)::PetscErrorCode
end

function AOSetFromOptions(arg1)
    @ccall libpetsc.AOSetFromOptions(arg1::AO)::PetscErrorCode
end

function AOCreateBasic(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.AOCreateBasic(arg1::MPI_Comm, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{AO})::PetscErrorCode
end

function AOCreateBasicIS(arg1, arg2, arg3)
    @ccall libpetsc.AOCreateBasicIS(arg1::IS, arg2::IS, arg3::Ptr{AO})::PetscErrorCode
end

function AOCreateMemoryScalable(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.AOCreateMemoryScalable(arg1::MPI_Comm, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{AO})::PetscErrorCode
end

function AOCreateMemoryScalableIS(arg1, arg2, arg3)
    @ccall libpetsc.AOCreateMemoryScalableIS(arg1::IS, arg2::IS, arg3::Ptr{AO})::PetscErrorCode
end

function AOCreateMapping(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.AOCreateMapping(arg1::MPI_Comm, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{AO})::PetscErrorCode
end

function AOCreateMappingIS(arg1, arg2, arg3)
    @ccall libpetsc.AOCreateMappingIS(arg1::IS, arg2::IS, arg3::Ptr{AO})::PetscErrorCode
end

function AOView(arg1, arg2)
    @ccall libpetsc.AOView(arg1::AO, arg2::PetscViewer)::PetscErrorCode
end

function AOViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.AOViewFromOptions(arg1::AO, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function AODestroy(arg1)
    @ccall libpetsc.AODestroy(arg1::Ptr{AO})::PetscErrorCode
end

function AOSetType(arg1, arg2)
    @ccall libpetsc.AOSetType(arg1::AO, arg2::AOType)::PetscErrorCode
end

function AOGetType(arg1, arg2)
    @ccall libpetsc.AOGetType(arg1::AO, arg2::Ptr{AOType})::PetscErrorCode
end

function AORegister(arg1, arg2)
    @ccall libpetsc.AORegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function AOPetscToApplication(arg1, arg2, arg3)
    @ccall libpetsc.AOPetscToApplication(arg1::AO, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function AOApplicationToPetsc(arg1, arg2, arg3)
    @ccall libpetsc.AOApplicationToPetsc(arg1::AO, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function AOPetscToApplicationIS(arg1, arg2)
    @ccall libpetsc.AOPetscToApplicationIS(arg1::AO, arg2::IS)::PetscErrorCode
end

function AOApplicationToPetscIS(arg1, arg2)
    @ccall libpetsc.AOApplicationToPetscIS(arg1::AO, arg2::IS)::PetscErrorCode
end

function AOPetscToApplicationPermuteInt(arg1, arg2, arg3)
    @ccall libpetsc.AOPetscToApplicationPermuteInt(arg1::AO, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function AOApplicationToPetscPermuteInt(arg1, arg2, arg3)
    @ccall libpetsc.AOApplicationToPetscPermuteInt(arg1::AO, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function AOPetscToApplicationPermuteReal(arg1, arg2, arg3)
    @ccall libpetsc.AOPetscToApplicationPermuteReal(arg1::AO, arg2::PetscInt, arg3::Ptr{PetscReal})::PetscErrorCode
end

function AOApplicationToPetscPermuteReal(arg1, arg2, arg3)
    @ccall libpetsc.AOApplicationToPetscPermuteReal(arg1::AO, arg2::PetscInt, arg3::Ptr{PetscReal})::PetscErrorCode
end

function AOMappingHasApplicationIndex(arg1, arg2, arg3)
    @ccall libpetsc.AOMappingHasApplicationIndex(arg1::AO, arg2::PetscInt, arg3::Ptr{PetscBool})::PetscErrorCode
end

function AOMappingHasPetscIndex(arg1, arg2, arg3)
    @ccall libpetsc.AOMappingHasPetscIndex(arg1::AO, arg2::PetscInt, arg3::Ptr{PetscBool})::PetscErrorCode
end

mutable struct _p_PetscQuadrature end

const PetscQuadrature = Ptr{_p_PetscQuadrature}

@enum PetscGaussLobattoLegendreCreateType::UInt32 begin
    PETSCGAUSSLOBATTOLEGENDRE_VIA_LINEAR_ALGEBRA = 0
    PETSCGAUSSLOBATTOLEGENDRE_VIA_NEWTON = 1
end

@enum PetscDTNodeType::Int32 begin
    PETSCDTNODES_DEFAULT = -1
    PETSCDTNODES_GAUSSJACOBI = 0
    PETSCDTNODES_EQUISPACED = 1
    PETSCDTNODES_TANHSINH = 2
end

function PetscQuadratureCreate(arg1, arg2)
    @ccall libpetsc.PetscQuadratureCreate(arg1::MPI_Comm, arg2::Ptr{PetscQuadrature})::PetscErrorCode
end

function PetscQuadratureDuplicate(arg1, arg2)
    @ccall libpetsc.PetscQuadratureDuplicate(arg1::PetscQuadrature, arg2::Ptr{PetscQuadrature})::PetscErrorCode
end

function PetscQuadratureGetOrder(arg1, arg2)
    @ccall libpetsc.PetscQuadratureGetOrder(arg1::PetscQuadrature, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscQuadratureSetOrder(arg1, arg2)
    @ccall libpetsc.PetscQuadratureSetOrder(arg1::PetscQuadrature, arg2::PetscInt)::PetscErrorCode
end

function PetscQuadratureGetNumComponents(arg1, arg2)
    @ccall libpetsc.PetscQuadratureGetNumComponents(arg1::PetscQuadrature, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscQuadratureSetNumComponents(arg1, arg2)
    @ccall libpetsc.PetscQuadratureSetNumComponents(arg1::PetscQuadrature, arg2::PetscInt)::PetscErrorCode
end

function PetscQuadratureGetData(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscQuadratureGetData(arg1::PetscQuadrature, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{Ptr{PetscReal}}, arg6::Ptr{Ptr{PetscReal}})::PetscErrorCode
end

function PetscQuadratureSetData(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscQuadratureSetData(arg1::PetscQuadrature, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal})::PetscErrorCode
end

function PetscQuadratureView(arg1, arg2)
    @ccall libpetsc.PetscQuadratureView(arg1::PetscQuadrature, arg2::PetscViewer)::PetscErrorCode
end

function PetscQuadratureDestroy(arg1)
    @ccall libpetsc.PetscQuadratureDestroy(arg1::Ptr{PetscQuadrature})::PetscErrorCode
end

function PetscQuadratureExpandComposite(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscQuadratureExpandComposite(arg1::PetscQuadrature, arg2::PetscInt, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscQuadrature})::PetscErrorCode
end

function PetscQuadraturePushForward(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscQuadraturePushForward(arg1::PetscQuadrature, arg2::PetscInt, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal}, arg6::PetscInt, arg7::Ptr{PetscQuadrature})::PetscErrorCode
end

function PetscDTLegendreEval(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscDTLegendreEval(arg1::PetscInt, arg2::Ptr{PetscReal}, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal}, arg7::Ptr{PetscReal})::PetscErrorCode
end

function PetscDTJacobiNorm(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDTJacobiNorm(arg1::PetscReal, arg2::PetscReal, arg3::PetscInt, arg4::Ptr{PetscReal})::PetscErrorCode
end

function PetscDTJacobiEval(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.PetscDTJacobiEval(arg1::PetscInt, arg2::PetscReal, arg3::PetscReal, arg4::Ptr{PetscReal}, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::Ptr{PetscReal}, arg8::Ptr{PetscReal}, arg9::Ptr{PetscReal})::PetscErrorCode
end

function PetscDTJacobiEvalJet(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscDTJacobiEvalJet(arg1::PetscReal, arg2::PetscReal, arg3::PetscInt, arg4::Ptr{PetscReal}, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{PetscReal})::PetscErrorCode
end

function PetscDTPKDEvalJet(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscDTPKDEvalJet(arg1::PetscInt, arg2::PetscInt, arg3::Ptr{PetscReal}, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{PetscReal})::PetscErrorCode
end

function PetscDTGaussQuadrature(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDTGaussQuadrature(arg1::PetscInt, arg2::PetscReal, arg3::PetscReal, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal})::PetscErrorCode
end

function PetscDTGaussJacobiQuadrature(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscDTGaussJacobiQuadrature(arg1::PetscInt, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal, arg6::Ptr{PetscReal}, arg7::Ptr{PetscReal})::PetscErrorCode
end

function PetscDTGaussLobattoJacobiQuadrature(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscDTGaussLobattoJacobiQuadrature(arg1::PetscInt, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal, arg6::Ptr{PetscReal}, arg7::Ptr{PetscReal})::PetscErrorCode
end

function PetscDTGaussLobattoLegendreQuadrature(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDTGaussLobattoLegendreQuadrature(arg1::PetscInt, arg2::PetscGaussLobattoLegendreCreateType, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal})::PetscErrorCode
end

function PetscDTReconstructPoly(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscDTReconstructPoly(arg1::PetscInt, arg2::PetscInt, arg3::Ptr{PetscReal}, arg4::PetscInt, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal})::PetscErrorCode
end

function PetscDTGaussTensorQuadrature(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscDTGaussTensorQuadrature(arg1::PetscInt, arg2::PetscInt, arg3::PetscInt, arg4::PetscReal, arg5::PetscReal, arg6::Ptr{PetscQuadrature})::PetscErrorCode
end

function PetscDTStroudConicalQuadrature(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscDTStroudConicalQuadrature(arg1::PetscInt, arg2::PetscInt, arg3::PetscInt, arg4::PetscReal, arg5::PetscReal, arg6::Ptr{PetscQuadrature})::PetscErrorCode
end

function PetscDTTanhSinhTensorQuadrature(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDTTanhSinhTensorQuadrature(arg1::PetscInt, arg2::PetscInt, arg3::PetscReal, arg4::PetscReal, arg5::Ptr{PetscQuadrature})::PetscErrorCode
end

function PetscDTTanhSinhIntegrate(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDTTanhSinhIntegrate(arg1::Ptr{Cvoid}, arg2::PetscReal, arg3::PetscReal, arg4::PetscInt, arg5::Ptr{PetscReal})::PetscErrorCode
end

function PetscDTTanhSinhIntegrateMPFR(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDTTanhSinhIntegrateMPFR(arg1::Ptr{Cvoid}, arg2::PetscReal, arg3::PetscReal, arg4::PetscInt, arg5::Ptr{PetscReal})::PetscErrorCode
end

function PetscGaussLobattoLegendreIntegrate(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscGaussLobattoLegendreIntegrate(arg1::PetscInt, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal})::PetscErrorCode
end

function PetscGaussLobattoLegendreElementLaplacianCreate(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscGaussLobattoLegendreElementLaplacianCreate(arg1::PetscInt, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{Ptr{Ptr{PetscReal}}})::PetscErrorCode
end

function PetscGaussLobattoLegendreElementLaplacianDestroy(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscGaussLobattoLegendreElementLaplacianDestroy(arg1::PetscInt, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{Ptr{Ptr{PetscReal}}})::PetscErrorCode
end

function PetscGaussLobattoLegendreElementGradientCreate(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscGaussLobattoLegendreElementGradientCreate(arg1::PetscInt, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{Ptr{Ptr{PetscReal}}}, arg5::Ptr{Ptr{Ptr{PetscReal}}})::PetscErrorCode
end

function PetscGaussLobattoLegendreElementGradientDestroy(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscGaussLobattoLegendreElementGradientDestroy(arg1::PetscInt, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{Ptr{Ptr{PetscReal}}}, arg5::Ptr{Ptr{Ptr{PetscReal}}})::PetscErrorCode
end

function PetscGaussLobattoLegendreElementAdvectionCreate(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscGaussLobattoLegendreElementAdvectionCreate(arg1::PetscInt, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{Ptr{Ptr{PetscReal}}})::PetscErrorCode
end

function PetscGaussLobattoLegendreElementAdvectionDestroy(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscGaussLobattoLegendreElementAdvectionDestroy(arg1::PetscInt, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{Ptr{Ptr{PetscReal}}})::PetscErrorCode
end

function PetscGaussLobattoLegendreElementMassCreate(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscGaussLobattoLegendreElementMassCreate(arg1::PetscInt, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{Ptr{Ptr{PetscReal}}})::PetscErrorCode
end

function PetscGaussLobattoLegendreElementMassDestroy(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscGaussLobattoLegendreElementMassDestroy(arg1::PetscInt, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{Ptr{Ptr{PetscReal}}})::PetscErrorCode
end

function PetscDTAltVApply(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDTAltVApply(arg1::PetscInt, arg2::PetscInt, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal})::PetscErrorCode
end

function PetscDTAltVWedge(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscDTAltVWedge(arg1::PetscInt, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal})::PetscErrorCode
end

function PetscDTAltVWedgeMatrix(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDTAltVWedgeMatrix(arg1::PetscInt, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal})::PetscErrorCode
end

function PetscDTAltVPullback(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscDTAltVPullback(arg1::PetscInt, arg2::PetscInt, arg3::Ptr{PetscReal}, arg4::PetscInt, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal})::PetscErrorCode
end

function PetscDTAltVPullbackMatrix(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDTAltVPullbackMatrix(arg1::PetscInt, arg2::PetscInt, arg3::Ptr{PetscReal}, arg4::PetscInt, arg5::Ptr{PetscReal})::PetscErrorCode
end

function PetscDTAltVInterior(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDTAltVInterior(arg1::PetscInt, arg2::PetscInt, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal})::PetscErrorCode
end

function PetscDTAltVInteriorMatrix(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDTAltVInteriorMatrix(arg1::PetscInt, arg2::PetscInt, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal})::PetscErrorCode
end

function PetscDTAltVInteriorPattern(arg1, arg2, arg3)
    @ccall libpetsc.PetscDTAltVInteriorPattern(arg1::PetscInt, arg2::PetscInt, arg3::Ptr{NTuple{3, PetscInt}})::PetscErrorCode
end

function PetscDTAltVStar(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDTAltVStar(arg1::PetscInt, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal})::PetscErrorCode
end

function PetscDTBaryToIndex(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDTBaryToIndex(arg1::PetscInt, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function PetscDTIndexToBary(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDTIndexToBary(arg1::PetscInt, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt})::PetscErrorCode
end

function PetscDTGradedOrderToIndex(arg1, arg2, arg3)
    @ccall libpetsc.PetscDTGradedOrderToIndex(arg1::PetscInt, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscDTIndexToGradedOrder(arg1, arg2, arg3)
    @ccall libpetsc.PetscDTIndexToGradedOrder(arg1::PetscInt, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscDTFactorial(n, factorial)
    @ccall libpetsc.PetscDTFactorial(n::PetscInt, factorial::Ptr{PetscReal})::PetscErrorCode
end

function PetscDTFactorialInt(n, factorial)
    @ccall libpetsc.PetscDTFactorialInt(n::PetscInt, factorial::Ptr{PetscInt})::PetscErrorCode
end

function PetscDTBinomial(n, k, binomial)
    @ccall libpetsc.PetscDTBinomial(n::PetscInt, k::PetscInt, binomial::Ptr{PetscReal})::PetscErrorCode
end

function PetscDTBinomialInt(n, k, binomial)
    @ccall libpetsc.PetscDTBinomialInt(n::PetscInt, k::PetscInt, binomial::Ptr{PetscInt})::PetscErrorCode
end

function PetscDTEnumPerm(n, k, perm, isOdd)
    @ccall libpetsc.PetscDTEnumPerm(n::PetscInt, k::PetscInt, perm::Ptr{PetscInt}, isOdd::Ptr{PetscBool})::PetscErrorCode
end

function PetscDTPermIndex(n, perm, k, isOdd)
    @ccall libpetsc.PetscDTPermIndex(n::PetscInt, perm::Ptr{PetscInt}, k::Ptr{PetscInt}, isOdd::Ptr{PetscBool})::PetscErrorCode
end

function PetscDTEnumSubset(n, k, j, subset)
    @ccall libpetsc.PetscDTEnumSubset(n::PetscInt, k::PetscInt, j::PetscInt, subset::Ptr{PetscInt})::PetscErrorCode
end

function PetscDTSubsetIndex(n, k, subset, index)
    @ccall libpetsc.PetscDTSubsetIndex(n::PetscInt, k::PetscInt, subset::Ptr{PetscInt}, index::Ptr{PetscInt})::PetscErrorCode
end

function PetscDTEnumSplit(n, k, j, perm, isOdd)
    @ccall libpetsc.PetscDTEnumSplit(n::PetscInt, k::PetscInt, j::PetscInt, perm::Ptr{PetscInt}, isOdd::Ptr{PetscBool})::PetscErrorCode
end

mutable struct _p_PetscTabulation
    K::PetscInt
    Nr::PetscInt
    Np::PetscInt
    Nb::PetscInt
    Nc::PetscInt
    cdim::PetscInt
    T::Ptr{Ptr{PetscReal}}
    _p_PetscTabulation() = new()
end

const PetscTabulation = Ptr{_p_PetscTabulation}

mutable struct _n_PetscFEGeom
    xi::Ptr{PetscReal}
    v::Ptr{PetscReal}
    J::Ptr{PetscReal}
    invJ::Ptr{PetscReal}
    detJ::Ptr{PetscReal}
    n::Ptr{PetscReal}
    face::Ptr{NTuple{2, PetscInt}}
    suppJ::NTuple{2, Ptr{PetscReal}}
    suppInvJ::NTuple{2, Ptr{PetscReal}}
    suppDetJ::NTuple{2, Ptr{PetscReal}}
    dim::PetscInt
    dimEmbed::PetscInt
    numCells::PetscInt
    numPoints::PetscInt
    isAffine::PetscBool
    _n_PetscFEGeom() = new()
end

const PetscFEGeom = _n_PetscFEGeom

function PetscFEInitializePackage()
    @ccall libpetsc.PetscFEInitializePackage()::PetscErrorCode
end

const PetscSpaceType = Ptr{Cchar}

function PetscSpaceCreate(arg1, arg2)
    @ccall libpetsc.PetscSpaceCreate(arg1::MPI_Comm, arg2::Ptr{PetscSpace})::PetscErrorCode
end

function PetscSpaceDestroy(arg1)
    @ccall libpetsc.PetscSpaceDestroy(arg1::Ptr{PetscSpace})::PetscErrorCode
end

function PetscSpaceSetType(arg1, arg2)
    @ccall libpetsc.PetscSpaceSetType(arg1::PetscSpace, arg2::PetscSpaceType)::PetscErrorCode
end

function PetscSpaceGetType(arg1, arg2)
    @ccall libpetsc.PetscSpaceGetType(arg1::PetscSpace, arg2::Ptr{PetscSpaceType})::PetscErrorCode
end

function PetscSpaceSetUp(arg1)
    @ccall libpetsc.PetscSpaceSetUp(arg1::PetscSpace)::PetscErrorCode
end

function PetscSpaceSetFromOptions(arg1)
    @ccall libpetsc.PetscSpaceSetFromOptions(arg1::PetscSpace)::PetscErrorCode
end

function PetscSpaceViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.PetscSpaceViewFromOptions(arg1::PetscSpace, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function PetscSpaceView(arg1, arg2)
    @ccall libpetsc.PetscSpaceView(arg1::PetscSpace, arg2::PetscViewer)::PetscErrorCode
end

function PetscSpaceRegister(arg1, arg2)
    @ccall libpetsc.PetscSpaceRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscSpaceRegisterDestroy()
    @ccall libpetsc.PetscSpaceRegisterDestroy()::PetscErrorCode
end

function PetscSpaceGetDimension(arg1, arg2)
    @ccall libpetsc.PetscSpaceGetDimension(arg1::PetscSpace, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscSpaceSetNumComponents(arg1, arg2)
    @ccall libpetsc.PetscSpaceSetNumComponents(arg1::PetscSpace, arg2::PetscInt)::PetscErrorCode
end

function PetscSpaceGetNumComponents(arg1, arg2)
    @ccall libpetsc.PetscSpaceGetNumComponents(arg1::PetscSpace, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscSpaceSetNumVariables(arg1, arg2)
    @ccall libpetsc.PetscSpaceSetNumVariables(arg1::PetscSpace, arg2::PetscInt)::PetscErrorCode
end

function PetscSpaceGetNumVariables(arg1, arg2)
    @ccall libpetsc.PetscSpaceGetNumVariables(arg1::PetscSpace, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscSpaceSetDegree(arg1, arg2, arg3)
    @ccall libpetsc.PetscSpaceSetDegree(arg1::PetscSpace, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function PetscSpaceGetDegree(arg1, arg2, arg3)
    @ccall libpetsc.PetscSpaceGetDegree(arg1::PetscSpace, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscSpaceEvaluate(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscSpaceEvaluate(arg1::PetscSpace, arg2::PetscInt, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal})::PetscErrorCode
end

function PetscSpaceGetHeightSubspace(arg1, arg2, arg3)
    @ccall libpetsc.PetscSpaceGetHeightSubspace(arg1::PetscSpace, arg2::PetscInt, arg3::Ptr{PetscSpace})::PetscErrorCode
end

function PetscSpacePolynomialSetSymmetric(arg1, arg2)
    @ccall libpetsc.PetscSpacePolynomialSetSymmetric(arg1::PetscSpace, arg2::PetscBool)::PetscErrorCode
end

function PetscSpacePolynomialGetSymmetric(arg1, arg2)
    @ccall libpetsc.PetscSpacePolynomialGetSymmetric(arg1::PetscSpace, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscSpacePolynomialSetTensor(arg1, arg2)
    @ccall libpetsc.PetscSpacePolynomialSetTensor(arg1::PetscSpace, arg2::PetscBool)::PetscErrorCode
end

function PetscSpacePolynomialGetTensor(arg1, arg2)
    @ccall libpetsc.PetscSpacePolynomialGetTensor(arg1::PetscSpace, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscSpaceTensorSetNumSubspaces(arg1, arg2)
    @ccall libpetsc.PetscSpaceTensorSetNumSubspaces(arg1::PetscSpace, arg2::PetscInt)::PetscErrorCode
end

function PetscSpaceTensorGetNumSubspaces(arg1, arg2)
    @ccall libpetsc.PetscSpaceTensorGetNumSubspaces(arg1::PetscSpace, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscSpaceTensorSetSubspace(arg1, arg2, arg3)
    @ccall libpetsc.PetscSpaceTensorSetSubspace(arg1::PetscSpace, arg2::PetscInt, arg3::PetscSpace)::PetscErrorCode
end

function PetscSpaceTensorGetSubspace(arg1, arg2, arg3)
    @ccall libpetsc.PetscSpaceTensorGetSubspace(arg1::PetscSpace, arg2::PetscInt, arg3::Ptr{PetscSpace})::PetscErrorCode
end

function PetscSpaceSumSetNumSubspaces(arg1, arg2)
    @ccall libpetsc.PetscSpaceSumSetNumSubspaces(arg1::PetscSpace, arg2::PetscInt)::PetscErrorCode
end

function PetscSpaceSumGetNumSubspaces(arg1, arg2)
    @ccall libpetsc.PetscSpaceSumGetNumSubspaces(arg1::PetscSpace, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscSpaceSumSetSubspace(arg1, arg2, arg3)
    @ccall libpetsc.PetscSpaceSumSetSubspace(arg1::PetscSpace, arg2::PetscInt, arg3::PetscSpace)::PetscErrorCode
end

function PetscSpaceSumGetSubspace(arg1, arg2, arg3)
    @ccall libpetsc.PetscSpaceSumGetSubspace(arg1::PetscSpace, arg2::PetscInt, arg3::Ptr{PetscSpace})::PetscErrorCode
end

function PetscSpaceSumSetConcatenate(arg1, arg2)
    @ccall libpetsc.PetscSpaceSumSetConcatenate(arg1::PetscSpace, arg2::PetscBool)::PetscErrorCode
end

function PetscSpaceSumGetConcatenate(arg1, arg2)
    @ccall libpetsc.PetscSpaceSumGetConcatenate(arg1::PetscSpace, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscSpaceCreateSum(numSubspaces, subspaces, concatenate, sumSpace)
    @ccall libpetsc.PetscSpaceCreateSum(numSubspaces::PetscInt, subspaces::Ptr{PetscSpace}, concatenate::PetscBool, sumSpace::Ptr{PetscSpace})::PetscErrorCode
end

function PetscSpacePointGetPoints(arg1, arg2)
    @ccall libpetsc.PetscSpacePointGetPoints(arg1::PetscSpace, arg2::Ptr{PetscQuadrature})::PetscErrorCode
end

function PetscSpacePointSetPoints(arg1, arg2)
    @ccall libpetsc.PetscSpacePointSetPoints(arg1::PetscSpace, arg2::PetscQuadrature)::PetscErrorCode
end

function PetscSpaceCreateSubspace(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscSpaceCreateSubspace(arg1::PetscSpace, arg2::PetscDualSpace, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal}, arg7::PetscCopyMode, arg8::Ptr{PetscSpace})::PetscErrorCode
end

const PetscDualSpaceType = Ptr{Cchar}

function PetscDualSpaceCreate(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceCreate(arg1::MPI_Comm, arg2::Ptr{PetscDualSpace})::PetscErrorCode
end

function PetscDualSpaceDestroy(arg1)
    @ccall libpetsc.PetscDualSpaceDestroy(arg1::Ptr{PetscDualSpace})::PetscErrorCode
end

function PetscDualSpaceDuplicate(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceDuplicate(arg1::PetscDualSpace, arg2::Ptr{PetscDualSpace})::PetscErrorCode
end

function PetscDualSpaceSetType(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceSetType(arg1::PetscDualSpace, arg2::PetscDualSpaceType)::PetscErrorCode
end

function PetscDualSpaceGetType(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceGetType(arg1::PetscDualSpace, arg2::Ptr{PetscDualSpaceType})::PetscErrorCode
end

function PetscDualSpaceGetUniform(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceGetUniform(arg1::PetscDualSpace, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscDualSpaceGetNumDof(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceGetNumDof(arg1::PetscDualSpace, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function PetscDualSpaceGetSection(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceGetSection(arg1::PetscDualSpace, arg2::Ptr{PetscSection})::PetscErrorCode
end

function PetscDualSpaceSetUp(arg1)
    @ccall libpetsc.PetscDualSpaceSetUp(arg1::PetscDualSpace)::PetscErrorCode
end

function PetscDualSpaceSetFromOptions(arg1)
    @ccall libpetsc.PetscDualSpaceSetFromOptions(arg1::PetscDualSpace)::PetscErrorCode
end

function PetscDualSpaceViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.PetscDualSpaceViewFromOptions(arg1::PetscDualSpace, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function PetscDualSpaceView(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceView(arg1::PetscDualSpace, arg2::PetscViewer)::PetscErrorCode
end

function PetscDualSpaceRegister(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscDualSpaceRegisterDestroy()
    @ccall libpetsc.PetscDualSpaceRegisterDestroy()::PetscErrorCode
end

function PetscDualSpaceGetDimension(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceGetDimension(arg1::PetscDualSpace, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscDualSpaceGetInteriorDimension(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceGetInteriorDimension(arg1::PetscDualSpace, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscDualSpaceSetNumComponents(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceSetNumComponents(arg1::PetscDualSpace, arg2::PetscInt)::PetscErrorCode
end

function PetscDualSpaceGetNumComponents(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceGetNumComponents(arg1::PetscDualSpace, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscDualSpaceSetOrder(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceSetOrder(arg1::PetscDualSpace, arg2::PetscInt)::PetscErrorCode
end

function PetscDualSpaceGetOrder(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceGetOrder(arg1::PetscDualSpace, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscDualSpaceSetDM(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceSetDM(arg1::PetscDualSpace, arg2::DM)::PetscErrorCode
end

function PetscDualSpaceGetDM(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceGetDM(arg1::PetscDualSpace, arg2::Ptr{DM})::PetscErrorCode
end

function PetscDualSpaceGetFunctional(arg1, arg2, arg3)
    @ccall libpetsc.PetscDualSpaceGetFunctional(arg1::PetscDualSpace, arg2::PetscInt, arg3::Ptr{PetscQuadrature})::PetscErrorCode
end

function PetscDualSpaceCreateReferenceCell(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDualSpaceCreateReferenceCell(arg1::PetscDualSpace, arg2::PetscInt, arg3::PetscBool, arg4::Ptr{DM})::PetscErrorCode
end

function PetscDualSpaceGetSymmetries(arg1, arg2, arg3)
    @ccall libpetsc.PetscDualSpaceGetSymmetries(arg1::PetscDualSpace, arg2::Ptr{Ptr{Ptr{Ptr{PetscInt}}}}, arg3::Ptr{Ptr{Ptr{Ptr{PetscScalar}}}})::PetscErrorCode
end

function PetscFEGeomCreate(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscFEGeomCreate(arg1::PetscQuadrature, arg2::PetscInt, arg3::PetscInt, arg4::PetscBool, arg5::Ptr{Ptr{PetscFEGeom}})::PetscErrorCode
end

function PetscFEGeomGetChunk(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscFEGeomGetChunk(arg1::Ptr{PetscFEGeom}, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Ptr{PetscFEGeom}})::PetscErrorCode
end

function PetscFEGeomRestoreChunk(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscFEGeomRestoreChunk(arg1::Ptr{PetscFEGeom}, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Ptr{PetscFEGeom}})::PetscErrorCode
end

function PetscFEGeomComplete(arg1)
    @ccall libpetsc.PetscFEGeomComplete(arg1::Ptr{PetscFEGeom})::PetscErrorCode
end

function PetscFEGeomDestroy(arg1)
    @ccall libpetsc.PetscFEGeomDestroy(arg1::Ptr{Ptr{PetscFEGeom}})::PetscErrorCode
end

function PetscDualSpaceApply(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscDualSpaceApply(arg1::PetscDualSpace, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscFEGeom}, arg5::PetscInt, arg6::Ptr{Cvoid}, arg7::Ptr{Cvoid}, arg8::Ptr{PetscScalar})::PetscErrorCode
end

function PetscDualSpaceApplyDefault(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscDualSpaceApplyDefault(arg1::PetscDualSpace, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscFEGeom}, arg5::PetscInt, arg6::Ptr{Cvoid}, arg7::Ptr{Cvoid}, arg8::Ptr{PetscScalar})::PetscErrorCode
end

function PetscDualSpaceGetAllData(arg1, arg2, arg3)
    @ccall libpetsc.PetscDualSpaceGetAllData(arg1::PetscDualSpace, arg2::Ptr{PetscQuadrature}, arg3::Ptr{Mat})::PetscErrorCode
end

function PetscDualSpaceCreateAllDataDefault(arg1, arg2, arg3)
    @ccall libpetsc.PetscDualSpaceCreateAllDataDefault(arg1::PetscDualSpace, arg2::Ptr{PetscQuadrature}, arg3::Ptr{Mat})::PetscErrorCode
end

function PetscDualSpaceGetInteriorData(arg1, arg2, arg3)
    @ccall libpetsc.PetscDualSpaceGetInteriorData(arg1::PetscDualSpace, arg2::Ptr{PetscQuadrature}, arg3::Ptr{Mat})::PetscErrorCode
end

function PetscDualSpaceCreateInteriorDataDefault(arg1, arg2, arg3)
    @ccall libpetsc.PetscDualSpaceCreateInteriorDataDefault(arg1::PetscDualSpace, arg2::Ptr{PetscQuadrature}, arg3::Ptr{Mat})::PetscErrorCode
end

function PetscDualSpaceApplyAll(arg1, arg2, arg3)
    @ccall libpetsc.PetscDualSpaceApplyAll(arg1::PetscDualSpace, arg2::Ptr{PetscScalar}, arg3::Ptr{PetscScalar})::PetscErrorCode
end

function PetscDualSpaceApplyAllDefault(arg1, arg2, arg3)
    @ccall libpetsc.PetscDualSpaceApplyAllDefault(arg1::PetscDualSpace, arg2::Ptr{PetscScalar}, arg3::Ptr{PetscScalar})::PetscErrorCode
end

function PetscDualSpaceApplyInterior(arg1, arg2, arg3)
    @ccall libpetsc.PetscDualSpaceApplyInterior(arg1::PetscDualSpace, arg2::Ptr{PetscScalar}, arg3::Ptr{PetscScalar})::PetscErrorCode
end

function PetscDualSpaceApplyInteriorDefault(arg1, arg2, arg3)
    @ccall libpetsc.PetscDualSpaceApplyInteriorDefault(arg1::PetscDualSpace, arg2::Ptr{PetscScalar}, arg3::Ptr{PetscScalar})::PetscErrorCode
end

function PetscDualSpaceGetFormDegree(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceGetFormDegree(arg1::PetscDualSpace, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscDualSpaceSetFormDegree(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceSetFormDegree(arg1::PetscDualSpace, arg2::PetscInt)::PetscErrorCode
end

function PetscDualSpaceGetDeRahm(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceGetDeRahm(arg1::PetscDualSpace, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscDualSpaceTransform(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscDualSpaceTransform(arg1::PetscDualSpace, arg2::PetscDualSpaceTransformType, arg3::PetscBool, arg4::Ptr{PetscFEGeom}, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{PetscScalar})::PetscErrorCode
end

function PetscDualSpaceTransformGradient(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscDualSpaceTransformGradient(arg1::PetscDualSpace, arg2::PetscDualSpaceTransformType, arg3::PetscBool, arg4::Ptr{PetscFEGeom}, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{PetscScalar})::PetscErrorCode
end

function PetscDualSpaceTransformHessian(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscDualSpaceTransformHessian(arg1::PetscDualSpace, arg2::PetscDualSpaceTransformType, arg3::PetscBool, arg4::Ptr{PetscFEGeom}, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{PetscScalar})::PetscErrorCode
end

function PetscDualSpacePullback(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDualSpacePullback(arg1::PetscDualSpace, arg2::Ptr{PetscFEGeom}, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscScalar})::PetscErrorCode
end

function PetscDualSpacePushforward(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDualSpacePushforward(arg1::PetscDualSpace, arg2::Ptr{PetscFEGeom}, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscScalar})::PetscErrorCode
end

function PetscDualSpacePushforwardGradient(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDualSpacePushforwardGradient(arg1::PetscDualSpace, arg2::Ptr{PetscFEGeom}, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscScalar})::PetscErrorCode
end

function PetscDualSpacePushforwardHessian(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscDualSpacePushforwardHessian(arg1::PetscDualSpace, arg2::Ptr{PetscFEGeom}, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscScalar})::PetscErrorCode
end

function PetscDualSpaceLagrangeGetContinuity(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceLagrangeGetContinuity(arg1::PetscDualSpace, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscDualSpaceLagrangeSetContinuity(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceLagrangeSetContinuity(arg1::PetscDualSpace, arg2::PetscBool)::PetscErrorCode
end

function PetscDualSpaceLagrangeGetTensor(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceLagrangeGetTensor(arg1::PetscDualSpace, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscDualSpaceLagrangeSetTensor(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceLagrangeSetTensor(arg1::PetscDualSpace, arg2::PetscBool)::PetscErrorCode
end

function PetscDualSpaceLagrangeGetTrimmed(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceLagrangeGetTrimmed(arg1::PetscDualSpace, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscDualSpaceLagrangeSetTrimmed(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceLagrangeSetTrimmed(arg1::PetscDualSpace, arg2::PetscBool)::PetscErrorCode
end

function PetscDualSpaceLagrangeGetNodeType(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDualSpaceLagrangeGetNodeType(arg1::PetscDualSpace, arg2::Ptr{PetscDTNodeType}, arg3::Ptr{PetscBool}, arg4::Ptr{PetscReal})::PetscErrorCode
end

function PetscDualSpaceLagrangeSetNodeType(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDualSpaceLagrangeSetNodeType(arg1::PetscDualSpace, arg2::PetscDTNodeType, arg3::PetscBool, arg4::PetscReal)::PetscErrorCode
end

function PetscDualSpaceLagrangeGetUseMoments(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceLagrangeGetUseMoments(arg1::PetscDualSpace, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscDualSpaceLagrangeSetUseMoments(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceLagrangeSetUseMoments(arg1::PetscDualSpace, arg2::PetscBool)::PetscErrorCode
end

function PetscDualSpaceLagrangeGetMomentOrder(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceLagrangeGetMomentOrder(arg1::PetscDualSpace, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscDualSpaceLagrangeSetMomentOrder(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceLagrangeSetMomentOrder(arg1::PetscDualSpace, arg2::PetscInt)::PetscErrorCode
end

function PetscDualSpaceGetHeightSubspace(arg1, arg2, arg3)
    @ccall libpetsc.PetscDualSpaceGetHeightSubspace(arg1::PetscDualSpace, arg2::PetscInt, arg3::Ptr{PetscDualSpace})::PetscErrorCode
end

function PetscDualSpaceGetPointSubspace(arg1, arg2, arg3)
    @ccall libpetsc.PetscDualSpaceGetPointSubspace(arg1::PetscDualSpace, arg2::PetscInt, arg3::Ptr{PetscDualSpace})::PetscErrorCode
end

function PetscDualSpaceSimpleSetDimension(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceSimpleSetDimension(arg1::PetscDualSpace, arg2::PetscInt)::PetscErrorCode
end

function PetscDualSpaceSimpleSetFunctional(arg1, arg2, arg3)
    @ccall libpetsc.PetscDualSpaceSimpleSetFunctional(arg1::PetscDualSpace, arg2::PetscInt, arg3::PetscQuadrature)::PetscErrorCode
end

function PetscDualSpaceRefinedSetCellSpaces(arg1, arg2)
    @ccall libpetsc.PetscDualSpaceRefinedSetCellSpaces(arg1::PetscDualSpace, arg2::Ptr{PetscDualSpace})::PetscErrorCode
end

const PetscFEType = Ptr{Cchar}

function PetscFECreate(arg1, arg2)
    @ccall libpetsc.PetscFECreate(arg1::MPI_Comm, arg2::Ptr{PetscFE})::PetscErrorCode
end

function PetscFEDestroy(arg1)
    @ccall libpetsc.PetscFEDestroy(arg1::Ptr{PetscFE})::PetscErrorCode
end

function PetscFESetType(arg1, arg2)
    @ccall libpetsc.PetscFESetType(arg1::PetscFE, arg2::PetscFEType)::PetscErrorCode
end

function PetscFEGetType(arg1, arg2)
    @ccall libpetsc.PetscFEGetType(arg1::PetscFE, arg2::Ptr{PetscFEType})::PetscErrorCode
end

function PetscFESetUp(arg1)
    @ccall libpetsc.PetscFESetUp(arg1::PetscFE)::PetscErrorCode
end

function PetscFESetFromOptions(arg1)
    @ccall libpetsc.PetscFESetFromOptions(arg1::PetscFE)::PetscErrorCode
end

function PetscFEViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.PetscFEViewFromOptions(arg1::PetscFE, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function PetscFESetName(arg1, arg2)
    @ccall libpetsc.PetscFESetName(arg1::PetscFE, arg2::Ptr{Cchar})::PetscErrorCode
end

function PetscFEView(arg1, arg2)
    @ccall libpetsc.PetscFEView(arg1::PetscFE, arg2::PetscViewer)::PetscErrorCode
end

function PetscFERegister(arg1, arg2)
    @ccall libpetsc.PetscFERegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscFERegisterDestroy()
    @ccall libpetsc.PetscFERegisterDestroy()::PetscErrorCode
end

function PetscFECreateDefault(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscFECreateDefault(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscBool, arg5::Ptr{Cchar}, arg6::PetscInt, arg7::Ptr{PetscFE})::PetscErrorCode
end

function PetscFECreateLagrange(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscFECreateLagrange(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscBool, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{PetscFE})::PetscErrorCode
end

function PetscFEGetDimension(arg1, arg2)
    @ccall libpetsc.PetscFEGetDimension(arg1::PetscFE, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscFEGetSpatialDimension(arg1, arg2)
    @ccall libpetsc.PetscFEGetSpatialDimension(arg1::PetscFE, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscFESetNumComponents(arg1, arg2)
    @ccall libpetsc.PetscFESetNumComponents(arg1::PetscFE, arg2::PetscInt)::PetscErrorCode
end

function PetscFEGetNumComponents(arg1, arg2)
    @ccall libpetsc.PetscFEGetNumComponents(arg1::PetscFE, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscFEGetTileSizes(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscFEGetTileSizes(arg1::PetscFE, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt})::PetscErrorCode
end

function PetscFESetTileSizes(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscFESetTileSizes(arg1::PetscFE, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt)::PetscErrorCode
end

function PetscFESetBasisSpace(arg1, arg2)
    @ccall libpetsc.PetscFESetBasisSpace(arg1::PetscFE, arg2::PetscSpace)::PetscErrorCode
end

function PetscFEGetBasisSpace(arg1, arg2)
    @ccall libpetsc.PetscFEGetBasisSpace(arg1::PetscFE, arg2::Ptr{PetscSpace})::PetscErrorCode
end

function PetscFESetDualSpace(arg1, arg2)
    @ccall libpetsc.PetscFESetDualSpace(arg1::PetscFE, arg2::PetscDualSpace)::PetscErrorCode
end

function PetscFEGetDualSpace(arg1, arg2)
    @ccall libpetsc.PetscFEGetDualSpace(arg1::PetscFE, arg2::Ptr{PetscDualSpace})::PetscErrorCode
end

function PetscFESetQuadrature(arg1, arg2)
    @ccall libpetsc.PetscFESetQuadrature(arg1::PetscFE, arg2::PetscQuadrature)::PetscErrorCode
end

function PetscFEGetQuadrature(arg1, arg2)
    @ccall libpetsc.PetscFEGetQuadrature(arg1::PetscFE, arg2::Ptr{PetscQuadrature})::PetscErrorCode
end

function PetscFESetFaceQuadrature(arg1, arg2)
    @ccall libpetsc.PetscFESetFaceQuadrature(arg1::PetscFE, arg2::PetscQuadrature)::PetscErrorCode
end

function PetscFEGetFaceQuadrature(arg1, arg2)
    @ccall libpetsc.PetscFEGetFaceQuadrature(arg1::PetscFE, arg2::Ptr{PetscQuadrature})::PetscErrorCode
end

function PetscFECopyQuadrature(arg1, arg2)
    @ccall libpetsc.PetscFECopyQuadrature(arg1::PetscFE, arg2::PetscFE)::PetscErrorCode
end

function PetscFEGetNumDof(arg1, arg2)
    @ccall libpetsc.PetscFEGetNumDof(arg1::PetscFE, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function PetscFEGetCellTabulation(arg1, arg2, arg3)
    @ccall libpetsc.PetscFEGetCellTabulation(arg1::PetscFE, arg2::PetscInt, arg3::Ptr{PetscTabulation})::PetscErrorCode
end

function PetscFEGetFaceTabulation(arg1, arg2, arg3)
    @ccall libpetsc.PetscFEGetFaceTabulation(arg1::PetscFE, arg2::PetscInt, arg3::Ptr{PetscTabulation})::PetscErrorCode
end

function PetscFEGetFaceCentroidTabulation(arg1, arg2)
    @ccall libpetsc.PetscFEGetFaceCentroidTabulation(arg1::PetscFE, arg2::Ptr{PetscTabulation})::PetscErrorCode
end

function PetscFECreateTabulation(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscFECreateTabulation(arg1::PetscFE, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscReal}, arg5::PetscInt, arg6::Ptr{PetscTabulation})::PetscErrorCode
end

function PetscFEComputeTabulation(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscFEComputeTabulation(arg1::PetscFE, arg2::PetscInt, arg3::Ptr{PetscReal}, arg4::PetscInt, arg5::PetscTabulation)::PetscErrorCode
end

function PetscTabulationDestroy(arg1)
    @ccall libpetsc.PetscTabulationDestroy(arg1::Ptr{PetscTabulation})::PetscErrorCode
end

function PetscFERefine(arg1, arg2)
    @ccall libpetsc.PetscFERefine(arg1::PetscFE, arg2::Ptr{PetscFE})::PetscErrorCode
end

function PetscFEGetHeightSubspace(arg1, arg2, arg3)
    @ccall libpetsc.PetscFEGetHeightSubspace(arg1::PetscFE, arg2::PetscInt, arg3::Ptr{PetscFE})::PetscErrorCode
end

function PetscFECreateCellGeometry(arg1, arg2, arg3)
    @ccall libpetsc.PetscFECreateCellGeometry(arg1::PetscFE, arg2::PetscQuadrature, arg3::Ptr{PetscFEGeom})::PetscErrorCode
end

function PetscFEDestroyCellGeometry(arg1, arg2)
    @ccall libpetsc.PetscFEDestroyCellGeometry(arg1::PetscFE, arg2::Ptr{PetscFEGeom})::PetscErrorCode
end

function PetscFEPushforward(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscFEPushforward(arg1::PetscFE, arg2::Ptr{PetscFEGeom}, arg3::PetscInt, arg4::Ptr{PetscScalar})::PetscErrorCode
end

function PetscFEPushforwardGradient(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscFEPushforwardGradient(arg1::PetscFE, arg2::Ptr{PetscFEGeom}, arg3::PetscInt, arg4::Ptr{PetscScalar})::PetscErrorCode
end

function PetscFEPushforwardHessian(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscFEPushforwardHessian(arg1::PetscFE, arg2::Ptr{PetscFEGeom}, arg3::PetscInt, arg4::Ptr{PetscScalar})::PetscErrorCode
end

function PetscFEIntegrate(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscFEIntegrate(arg1::PetscDS, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscFEGeom}, arg5::Ptr{PetscScalar}, arg6::PetscDS, arg7::Ptr{PetscScalar}, arg8::Ptr{PetscScalar})::PetscErrorCode
end

function PetscFEIntegrateBd(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.PetscFEIntegrateBd(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{Cvoid}, arg4::PetscInt, arg5::Ptr{PetscFEGeom}, arg6::Ptr{PetscScalar}, arg7::PetscDS, arg8::Ptr{PetscScalar}, arg9::Ptr{PetscScalar})::PetscErrorCode
end

function PetscFEIntegrateResidual(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.PetscFEIntegrateResidual(arg1::PetscDS, arg2::PetscHashFormKey, arg3::PetscInt, arg4::Ptr{PetscFEGeom}, arg5::Ptr{PetscScalar}, arg6::Ptr{PetscScalar}, arg7::PetscDS, arg8::Ptr{PetscScalar}, arg9::PetscReal, arg10::Ptr{PetscScalar})::PetscErrorCode
end

function PetscFEIntegrateBdResidual(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.PetscFEIntegrateBdResidual(arg1::PetscDS, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscFEGeom}, arg5::Ptr{PetscScalar}, arg6::Ptr{PetscScalar}, arg7::PetscDS, arg8::Ptr{PetscScalar}, arg9::PetscReal, arg10::Ptr{PetscScalar})::PetscErrorCode
end

function PetscFEIntegrateHybridResidual(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.PetscFEIntegrateHybridResidual(arg1::PetscDS, arg2::PetscHashFormKey, arg3::PetscInt, arg4::Ptr{PetscFEGeom}, arg5::Ptr{PetscScalar}, arg6::Ptr{PetscScalar}, arg7::PetscDS, arg8::Ptr{PetscScalar}, arg9::PetscReal, arg10::Ptr{PetscScalar})::PetscErrorCode
end

function PetscFEIntegrateJacobian(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12)
    @ccall libpetsc.PetscFEIntegrateJacobian(arg1::PetscDS, arg2::PetscFEJacobianType, arg3::PetscHashFormKey, arg4::PetscInt, arg5::Ptr{PetscFEGeom}, arg6::Ptr{PetscScalar}, arg7::Ptr{PetscScalar}, arg8::PetscDS, arg9::Ptr{PetscScalar}, arg10::PetscReal, arg11::PetscReal, arg12::Ptr{PetscScalar})::PetscErrorCode
end

function PetscFEIntegrateBdJacobian(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12)
    @ccall libpetsc.PetscFEIntegrateBdJacobian(arg1::PetscDS, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscFEGeom}, arg6::Ptr{PetscScalar}, arg7::Ptr{PetscScalar}, arg8::PetscDS, arg9::Ptr{PetscScalar}, arg10::PetscReal, arg11::PetscReal, arg12::Ptr{PetscScalar})::PetscErrorCode
end

function PetscFEIntegrateHybridJacobian(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13)
    @ccall libpetsc.PetscFEIntegrateHybridJacobian(arg1::PetscDS, arg2::PetscFEJacobianType, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{PetscFEGeom}, arg7::Ptr{PetscScalar}, arg8::Ptr{PetscScalar}, arg9::PetscDS, arg10::Ptr{PetscScalar}, arg11::PetscReal, arg12::PetscReal, arg13::Ptr{PetscScalar})::PetscErrorCode
end

function PetscFECompositeGetMapping(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscFECompositeGetMapping(arg1::PetscFE, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{PetscReal}}, arg4::Ptr{Ptr{PetscReal}}, arg5::Ptr{Ptr{PetscReal}})::PetscErrorCode
end

function PetscFECreateHeightTrace(arg1, arg2, arg3)
    @ccall libpetsc.PetscFECreateHeightTrace(arg1::PetscFE, arg2::PetscInt, arg3::Ptr{PetscFE})::PetscErrorCode
end

function PetscFECreatePointTrace(arg1, arg2, arg3)
    @ccall libpetsc.PetscFECreatePointTrace(arg1::PetscFE, arg2::PetscInt, arg3::Ptr{PetscFE})::PetscErrorCode
end

function PetscFEOpenCLSetRealType(arg1, arg2)
    @ccall libpetsc.PetscFEOpenCLSetRealType(arg1::PetscFE, arg2::PetscDataType)::PetscErrorCode
end

function PetscFEOpenCLGetRealType(arg1, arg2)
    @ccall libpetsc.PetscFEOpenCLGetRealType(arg1::PetscFE, arg2::Ptr{PetscDataType})::PetscErrorCode
end

function DMDASetInterpolationType(arg1, arg2)
    @ccall libpetsc.DMDASetInterpolationType(arg1::DM, arg2::DMDAInterpolationType)::PetscErrorCode
end

function DMDAGetInterpolationType(arg1, arg2)
    @ccall libpetsc.DMDAGetInterpolationType(arg1::DM, arg2::Ptr{DMDAInterpolationType})::PetscErrorCode
end

function DMDACreateAggregates(arg1, arg2, arg3)
    @ccall libpetsc.DMDACreateAggregates(arg1::DM, arg2::DM, arg3::Ptr{Mat})::PetscErrorCode
end

function DMDASetElementType(arg1, arg2)
    @ccall libpetsc.DMDASetElementType(arg1::DM, arg2::DMDAElementType)::PetscErrorCode
end

function DMDAGetElementType(arg1, arg2)
    @ccall libpetsc.DMDAGetElementType(arg1::DM, arg2::Ptr{DMDAElementType})::PetscErrorCode
end

function DMDAGetElements(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMDAGetElements(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMDARestoreElements(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMDARestoreElements(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMDAGetElementsSizes(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMDAGetElementsSizes(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMDAGetElementsCorners(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMDAGetElementsCorners(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMDAGetSubdomainCornersIS(arg1, arg2)
    @ccall libpetsc.DMDAGetSubdomainCornersIS(arg1::DM, arg2::Ptr{IS})::PetscErrorCode
end

function DMDARestoreSubdomainCornersIS(arg1, arg2)
    @ccall libpetsc.DMDARestoreSubdomainCornersIS(arg1::DM, arg2::Ptr{IS})::PetscErrorCode
end

function DMDACreate(arg1, arg2)
    @ccall libpetsc.DMDACreate(arg1::MPI_Comm, arg2::Ptr{DM})::PetscErrorCode
end

function DMDASetSizes(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMDASetSizes(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt)::PetscErrorCode
end

function DMDACreate1d(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMDACreate1d(arg1::MPI_Comm, arg2::DMBoundaryType, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::Ptr{DM})::PetscErrorCode
end

function DMDACreate2d(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13)
    @ccall libpetsc.DMDACreate2d(arg1::MPI_Comm, arg2::DMBoundaryType, arg3::DMBoundaryType, arg4::DMDAStencilType, arg5::PetscInt, arg6::PetscInt, arg7::PetscInt, arg8::PetscInt, arg9::PetscInt, arg10::PetscInt, arg11::Ptr{PetscInt}, arg12::Ptr{PetscInt}, arg13::Ptr{DM})::PetscErrorCode
end

function DMDACreate3d(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13, arg14, arg15, arg16, arg17)
    @ccall libpetsc.DMDACreate3d(arg1::MPI_Comm, arg2::DMBoundaryType, arg3::DMBoundaryType, arg4::DMBoundaryType, arg5::DMDAStencilType, arg6::PetscInt, arg7::PetscInt, arg8::PetscInt, arg9::PetscInt, arg10::PetscInt, arg11::PetscInt, arg12::PetscInt, arg13::PetscInt, arg14::Ptr{PetscInt}, arg15::Ptr{PetscInt}, arg16::Ptr{PetscInt}, arg17::Ptr{DM})::PetscErrorCode
end

function DMDAGlobalToNaturalBegin(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMDAGlobalToNaturalBegin(arg1::DM, arg2::Vec, arg3::InsertMode, arg4::Vec)::PetscErrorCode
end

function DMDAGlobalToNaturalEnd(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMDAGlobalToNaturalEnd(arg1::DM, arg2::Vec, arg3::InsertMode, arg4::Vec)::PetscErrorCode
end

function DMDANaturalToGlobalBegin(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMDANaturalToGlobalBegin(arg1::DM, arg2::Vec, arg3::InsertMode, arg4::Vec)::PetscErrorCode
end

function DMDANaturalToGlobalEnd(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMDANaturalToGlobalEnd(arg1::DM, arg2::Vec, arg3::InsertMode, arg4::Vec)::PetscErrorCode
end

function DMDALocalToLocalBegin(dm, g, mode, l)
    @ccall libpetsc.DMDALocalToLocalBegin(dm::DM, g::Vec, mode::InsertMode, l::Vec)::PetscErrorCode
end

function DMDALocalToLocalEnd(dm, g, mode, l)
    @ccall libpetsc.DMDALocalToLocalEnd(dm::DM, g::Vec, mode::InsertMode, l::Vec)::PetscErrorCode
end

function DMDACreateNaturalVector(arg1, arg2)
    @ccall libpetsc.DMDACreateNaturalVector(arg1::DM, arg2::Ptr{Vec})::PetscErrorCode
end

function DMDAGetCorners(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMDAGetCorners(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscInt}, arg7::Ptr{PetscInt})::PetscErrorCode
end

function DMDAGetGhostCorners(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMDAGetGhostCorners(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscInt}, arg7::Ptr{PetscInt})::PetscErrorCode
end

function DMDAGetInfo(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13, arg14)
    @ccall libpetsc.DMDAGetInfo(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscInt}, arg7::Ptr{PetscInt}, arg8::Ptr{PetscInt}, arg9::Ptr{PetscInt}, arg10::Ptr{PetscInt}, arg11::Ptr{DMBoundaryType}, arg12::Ptr{DMBoundaryType}, arg13::Ptr{DMBoundaryType}, arg14::Ptr{DMDAStencilType})::PetscErrorCode
end

function DMDAGetProcessorSubset(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMDAGetProcessorSubset(arg1::DM, arg2::DMDirection, arg3::PetscInt, arg4::Ptr{MPI_Comm})::PetscErrorCode
end

function DMDAGetProcessorSubsets(arg1, arg2, arg3)
    @ccall libpetsc.DMDAGetProcessorSubsets(arg1::DM, arg2::DMDirection, arg3::Ptr{MPI_Comm})::PetscErrorCode
end

function DMDAGetRay(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMDAGetRay(arg1::DM, arg2::DMDirection, arg3::PetscInt, arg4::Ptr{Vec}, arg5::Ptr{VecScatter})::PetscErrorCode
end

function DMDAGlobalToNaturalAllCreate(arg1, arg2)
    @ccall libpetsc.DMDAGlobalToNaturalAllCreate(arg1::DM, arg2::Ptr{VecScatter})::PetscErrorCode
end

function DMDANaturalAllToGlobalCreate(arg1, arg2)
    @ccall libpetsc.DMDANaturalAllToGlobalCreate(arg1::DM, arg2::Ptr{VecScatter})::PetscErrorCode
end

function DMDAGetScatter(arg1, arg2, arg3)
    @ccall libpetsc.DMDAGetScatter(arg1::DM, arg2::Ptr{VecScatter}, arg3::Ptr{VecScatter})::PetscErrorCode
end

function DMDAGetNeighbors(arg1, arg2)
    @ccall libpetsc.DMDAGetNeighbors(arg1::DM, arg2::Ptr{Ptr{PetscMPIInt}})::PetscErrorCode
end

function DMDASetAOType(arg1, arg2)
    @ccall libpetsc.DMDASetAOType(arg1::DM, arg2::AOType)::PetscErrorCode
end

function DMDAGetAO(arg1, arg2)
    @ccall libpetsc.DMDAGetAO(arg1::DM, arg2::Ptr{AO})::PetscErrorCode
end

function DMDASetUniformCoordinates(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMDASetUniformCoordinates(arg1::DM, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal, arg6::PetscReal, arg7::PetscReal)::PetscErrorCode
end

function DMDASetGLLCoordinates(arg1, arg2, arg3)
    @ccall libpetsc.DMDASetGLLCoordinates(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscReal})::PetscErrorCode
end

function DMDAGetCoordinateArray(arg1, arg2)
    @ccall libpetsc.DMDAGetCoordinateArray(arg1::DM, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMDARestoreCoordinateArray(arg1, arg2)
    @ccall libpetsc.DMDARestoreCoordinateArray(arg1::DM, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMDAGetLogicalCoordinate(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.DMDAGetLogicalCoordinate(arg1::DM, arg2::PetscScalar, arg3::PetscScalar, arg4::PetscScalar, arg5::Ptr{PetscInt}, arg6::Ptr{PetscInt}, arg7::Ptr{PetscInt}, arg8::Ptr{PetscScalar}, arg9::Ptr{PetscScalar}, arg10::Ptr{PetscScalar})::PetscErrorCode
end

function DMDAMapCoordsToPeriodicDomain(arg1, arg2, arg3)
    @ccall libpetsc.DMDAMapCoordsToPeriodicDomain(arg1::DM, arg2::Ptr{PetscScalar}, arg3::Ptr{PetscScalar})::PetscErrorCode
end

function DMDACreateCompatibleDMDA(arg1, arg2, arg3)
    @ccall libpetsc.DMDACreateCompatibleDMDA(arg1::DM, arg2::PetscInt, arg3::Ptr{DM})::PetscErrorCode
end

function DMDAGetReducedDMDA(arg1, arg2, arg3)
    @ccall libpetsc.DMDAGetReducedDMDA(arg1::DM, arg2::PetscInt, arg3::Ptr{DM})::PetscErrorCode
end

function DMDASetFieldName(arg1, arg2, arg3)
    @ccall libpetsc.DMDASetFieldName(arg1::DM, arg2::PetscInt, arg3::Ptr{Cchar})::PetscErrorCode
end

function DMDAGetFieldName(arg1, arg2, arg3)
    @ccall libpetsc.DMDAGetFieldName(arg1::DM, arg2::PetscInt, arg3::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function DMDASetFieldNames(arg1, arg2)
    @ccall libpetsc.DMDASetFieldNames(arg1::DM, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function DMDAGetFieldNames(arg1, arg2)
    @ccall libpetsc.DMDAGetFieldNames(arg1::DM, arg2::Ptr{Ptr{Ptr{Cchar}}})::PetscErrorCode
end

function DMDASetCoordinateName(arg1, arg2, arg3)
    @ccall libpetsc.DMDASetCoordinateName(arg1::DM, arg2::PetscInt, arg3::Ptr{Cchar})::PetscErrorCode
end

function DMDAGetCoordinateName(arg1, arg2, arg3)
    @ccall libpetsc.DMDAGetCoordinateName(arg1::DM, arg2::PetscInt, arg3::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function DMDASetBoundaryType(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMDASetBoundaryType(arg1::DM, arg2::DMBoundaryType, arg3::DMBoundaryType, arg4::DMBoundaryType)::PetscErrorCode
end

function DMDASetDof(arg1, arg2)
    @ccall libpetsc.DMDASetDof(arg1::DM, arg2::PetscInt)::PetscErrorCode
end

function DMDAGetDof(arg1, arg2)
    @ccall libpetsc.DMDAGetDof(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMDASetOverlap(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMDASetOverlap(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt)::PetscErrorCode
end

function DMDAGetOverlap(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMDAGetOverlap(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMDASetNumLocalSubDomains(arg1, arg2)
    @ccall libpetsc.DMDASetNumLocalSubDomains(arg1::DM, arg2::PetscInt)::PetscErrorCode
end

function DMDAGetNumLocalSubDomains(arg1, arg2)
    @ccall libpetsc.DMDAGetNumLocalSubDomains(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMDAGetOffset(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMDAGetOffset(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscInt}, arg7::Ptr{PetscInt})::PetscErrorCode
end

function DMDASetOffset(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMDASetOffset(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::PetscInt)::PetscErrorCode
end

function DMDAGetNonOverlappingRegion(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMDAGetNonOverlappingRegion(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscInt}, arg7::Ptr{PetscInt})::PetscErrorCode
end

function DMDASetNonOverlappingRegion(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMDASetNonOverlappingRegion(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::PetscInt)::PetscErrorCode
end

function DMDASetStencilWidth(arg1, arg2)
    @ccall libpetsc.DMDASetStencilWidth(arg1::DM, arg2::PetscInt)::PetscErrorCode
end

function DMDAGetStencilWidth(arg1, arg2)
    @ccall libpetsc.DMDAGetStencilWidth(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMDASetOwnershipRanges(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMDASetOwnershipRanges(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMDAGetOwnershipRanges(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMDAGetOwnershipRanges(arg1::DM, arg2::Ptr{Ptr{PetscInt}}, arg3::Ptr{Ptr{PetscInt}}, arg4::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMDASetNumProcs(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMDASetNumProcs(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt)::PetscErrorCode
end

function DMDASetStencilType(arg1, arg2)
    @ccall libpetsc.DMDASetStencilType(arg1::DM, arg2::DMDAStencilType)::PetscErrorCode
end

function DMDAGetStencilType(arg1, arg2)
    @ccall libpetsc.DMDAGetStencilType(arg1::DM, arg2::Ptr{DMDAStencilType})::PetscErrorCode
end

function DMDAVecGetArray(arg1, arg2, arg3)
    @ccall libpetsc.DMDAVecGetArray(arg1::DM, arg2::Vec, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMDAVecRestoreArray(arg1, arg2, arg3)
    @ccall libpetsc.DMDAVecRestoreArray(arg1::DM, arg2::Vec, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMDAVecGetArrayWrite(arg1, arg2, arg3)
    @ccall libpetsc.DMDAVecGetArrayWrite(arg1::DM, arg2::Vec, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMDAVecRestoreArrayWrite(arg1, arg2, arg3)
    @ccall libpetsc.DMDAVecRestoreArrayWrite(arg1::DM, arg2::Vec, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMDAVecGetArrayDOF(arg1, arg2, arg3)
    @ccall libpetsc.DMDAVecGetArrayDOF(arg1::DM, arg2::Vec, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMDAVecRestoreArrayDOF(arg1, arg2, arg3)
    @ccall libpetsc.DMDAVecRestoreArrayDOF(arg1::DM, arg2::Vec, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMDAVecGetArrayRead(arg1, arg2, arg3)
    @ccall libpetsc.DMDAVecGetArrayRead(arg1::DM, arg2::Vec, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMDAVecRestoreArrayRead(arg1, arg2, arg3)
    @ccall libpetsc.DMDAVecRestoreArrayRead(arg1::DM, arg2::Vec, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMDAVecGetArrayDOFRead(arg1, arg2, arg3)
    @ccall libpetsc.DMDAVecGetArrayDOFRead(arg1::DM, arg2::Vec, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMDAVecRestoreArrayDOFRead(arg1, arg2, arg3)
    @ccall libpetsc.DMDAVecRestoreArrayDOFRead(arg1::DM, arg2::Vec, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMDACreatePatchIS(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMDACreatePatchIS(arg1::DM, arg2::Ptr{MatStencil}, arg3::Ptr{MatStencil}, arg4::Ptr{IS}, arg5::PetscBool)::PetscErrorCode
end

mutable struct DMDACoor2d
    x::PetscScalar
    y::PetscScalar
    DMDACoor2d() = new()
end

mutable struct DMDACoor3d
    x::PetscScalar
    y::PetscScalar
    z::PetscScalar
    DMDACoor3d() = new()
end

function DMDAGetLocalInfo(arg1, arg2)
    @ccall libpetsc.DMDAGetLocalInfo(arg1::DM, arg2::Ptr{DMDALocalInfo})::PetscErrorCode
end

function MatRegisterDAAD()
    @ccall libpetsc.MatRegisterDAAD()::PetscErrorCode
end

function MatCreateDAAD(arg1, arg2)
    @ccall libpetsc.MatCreateDAAD(arg1::DM, arg2::Ptr{Mat})::PetscErrorCode
end

function MatCreateSeqUSFFT(arg1, arg2, arg3)
    @ccall libpetsc.MatCreateSeqUSFFT(arg1::Vec, arg2::DM, arg3::Ptr{Mat})::PetscErrorCode
end

function DMDASetGetMatrix(arg1, arg2)
    @ccall libpetsc.DMDASetGetMatrix(arg1::DM, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMDASetBlockFills(arg1, arg2, arg3)
    @ccall libpetsc.DMDASetBlockFills(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMDASetBlockFillsSparse(arg1, arg2, arg3)
    @ccall libpetsc.DMDASetBlockFillsSparse(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMDASetRefinementFactor(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMDASetRefinementFactor(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt)::PetscErrorCode
end

function DMDAGetRefinementFactor(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMDAGetRefinementFactor(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMDAGetArray(arg1, arg2, arg3)
    @ccall libpetsc.DMDAGetArray(arg1::DM, arg2::PetscBool, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMDARestoreArray(arg1, arg2, arg3)
    @ccall libpetsc.DMDARestoreArray(arg1::DM, arg2::PetscBool, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMDACreatePF(arg1, arg2)
    @ccall libpetsc.DMDACreatePF(arg1::DM, arg2::Ptr{PF})::PetscErrorCode
end

function DMDAGetNumCells(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMDAGetNumCells(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt})::PetscErrorCode
end

function DMDAGetCellPoint(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMDAGetCellPoint(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscInt})::PetscErrorCode
end

function DMDAGetNumVertices(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMDAGetNumVertices(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt})::PetscErrorCode
end

function DMDAGetNumFaces(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMDAGetNumFaces(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscInt}, arg7::Ptr{PetscInt})::PetscErrorCode
end

function DMDAGetHeightStratum(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMDAGetHeightStratum(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMDAGetDepthStratum(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMDAGetDepthStratum(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMDAComputeCellGeometryFEM(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMDAComputeCellGeometryFEM(arg1::DM, arg2::PetscInt, arg3::PetscQuadrature, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal}, arg7::Ptr{PetscReal})::PetscErrorCode
end

function DMDAGetTransitiveClosure(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMDAGetTransitiveClosure(arg1::DM, arg2::PetscInt, arg3::PetscBool, arg4::Ptr{PetscInt}, arg5::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMDARestoreTransitiveClosure(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMDARestoreTransitiveClosure(arg1::DM, arg2::PetscInt, arg3::PetscBool, arg4::Ptr{PetscInt}, arg5::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMDAConvertToCell(arg1, arg2, arg3)
    @ccall libpetsc.DMDAConvertToCell(arg1::DM, arg2::MatStencil, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMDASetVertexCoordinates(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMDASetVertexCoordinates(arg1::DM, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal, arg6::PetscReal, arg7::PetscReal)::PetscErrorCode
end

function DMDASetPreallocationCenterDimension(arg1, arg2)
    @ccall libpetsc.DMDASetPreallocationCenterDimension(arg1::DM, arg2::PetscInt)::PetscErrorCode
end

function DMDAGetPreallocationCenterDimension(arg1, arg2)
    @ccall libpetsc.DMDAGetPreallocationCenterDimension(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMCompositeCreate(arg1, arg2)
    @ccall libpetsc.DMCompositeCreate(arg1::MPI_Comm, arg2::Ptr{DM})::PetscErrorCode
end

function DMCompositeAddDM(arg1, arg2)
    @ccall libpetsc.DMCompositeAddDM(arg1::DM, arg2::DM)::PetscErrorCode
end

function DMCompositeSetCoupling(arg1, arg2)
    @ccall libpetsc.DMCompositeSetCoupling(arg1::DM, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMCompositeAddVecScatter(arg1, arg2)
    @ccall libpetsc.DMCompositeAddVecScatter(arg1::DM, arg2::VecScatter)::PetscErrorCode
end

function DMCompositeScatterArray(arg1, arg2, arg3)
    @ccall libpetsc.DMCompositeScatterArray(arg1::DM, arg2::Vec, arg3::Ptr{Vec})::PetscErrorCode
end

function DMCompositeGatherArray(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMCompositeGatherArray(arg1::DM, arg2::InsertMode, arg3::Vec, arg4::Ptr{Vec})::PetscErrorCode
end

function DMCompositeGetNumberDM(arg1, arg2)
    @ccall libpetsc.DMCompositeGetNumberDM(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMCompositeGetAccessArray(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMCompositeGetAccessArray(arg1::DM, arg2::Vec, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::Ptr{Vec})::PetscErrorCode
end

function DMCompositeRestoreAccessArray(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMCompositeRestoreAccessArray(arg1::DM, arg2::Vec, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::Ptr{Vec})::PetscErrorCode
end

function DMCompositeGetLocalAccessArray(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMCompositeGetLocalAccessArray(arg1::DM, arg2::Vec, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::Ptr{Vec})::PetscErrorCode
end

function DMCompositeRestoreLocalAccessArray(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMCompositeRestoreLocalAccessArray(arg1::DM, arg2::Vec, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::Ptr{Vec})::PetscErrorCode
end

function DMCompositeGetEntriesArray(arg1, arg2)
    @ccall libpetsc.DMCompositeGetEntriesArray(arg1::DM, arg2::Ptr{DM})::PetscErrorCode
end

function DMCompositeGetGlobalISs(arg1, arg2)
    @ccall libpetsc.DMCompositeGetGlobalISs(arg1::DM, arg2::Ptr{Ptr{IS}})::PetscErrorCode
end

function DMCompositeGetLocalISs(arg1, arg2)
    @ccall libpetsc.DMCompositeGetLocalISs(arg1::DM, arg2::Ptr{Ptr{IS}})::PetscErrorCode
end

function DMCompositeGetISLocalToGlobalMappings(arg1, arg2)
    @ccall libpetsc.DMCompositeGetISLocalToGlobalMappings(arg1::DM, arg2::Ptr{Ptr{ISLocalToGlobalMapping}})::PetscErrorCode
end

function DMPatchCreate(arg1, arg2)
    @ccall libpetsc.DMPatchCreate(arg1::MPI_Comm, arg2::Ptr{DM})::PetscErrorCode
end

function DMPatchZoom(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMPatchZoom(arg1::DM, arg2::MatStencil, arg3::MatStencil, arg4::MPI_Comm, arg5::Ptr{DM}, arg6::Ptr{PetscSF}, arg7::Ptr{PetscSF})::PetscErrorCode
end

function DMPatchSolve(arg1)
    @ccall libpetsc.DMPatchSolve(arg1::DM)::PetscErrorCode
end

function DMPatchGetPatchSize(arg1, arg2)
    @ccall libpetsc.DMPatchGetPatchSize(arg1::DM, arg2::Ptr{MatStencil})::PetscErrorCode
end

function DMPatchSetPatchSize(arg1, arg2)
    @ccall libpetsc.DMPatchSetPatchSize(arg1::DM, arg2::MatStencil)::PetscErrorCode
end

function DMPatchGetCommSize(arg1, arg2)
    @ccall libpetsc.DMPatchGetCommSize(arg1::DM, arg2::Ptr{MatStencil})::PetscErrorCode
end

function DMPatchSetCommSize(arg1, arg2)
    @ccall libpetsc.DMPatchSetCommSize(arg1::DM, arg2::MatStencil)::PetscErrorCode
end

function DMPatchGetCoarse(arg1, arg2)
    @ccall libpetsc.DMPatchGetCoarse(arg1::DM, arg2::Ptr{DM})::PetscErrorCode
end

function DMPatchCreateGrid(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMPatchCreateGrid(arg1::MPI_Comm, arg2::PetscInt, arg3::MatStencil, arg4::MatStencil, arg5::MatStencil, arg6::Ptr{DM})::PetscErrorCode
end

mutable struct _p_PetscPartitioner end

const PetscPartitioner = Ptr{_p_PetscPartitioner}

function PetscPartitionerInitializePackage()
    @ccall libpetsc.PetscPartitionerInitializePackage()::PetscErrorCode
end

function PetscPartitionerFinalizePackage()
    @ccall libpetsc.PetscPartitionerFinalizePackage()::PetscErrorCode
end

const PetscPartitionerType = Ptr{Cchar}

function PetscPartitionerRegister(arg1, arg2)
    @ccall libpetsc.PetscPartitionerRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscPartitionerCreate(arg1, arg2)
    @ccall libpetsc.PetscPartitionerCreate(arg1::MPI_Comm, arg2::Ptr{PetscPartitioner})::PetscErrorCode
end

function PetscPartitionerDestroy(arg1)
    @ccall libpetsc.PetscPartitionerDestroy(arg1::Ptr{PetscPartitioner})::PetscErrorCode
end

function PetscPartitionerSetType(arg1, arg2)
    @ccall libpetsc.PetscPartitionerSetType(arg1::PetscPartitioner, arg2::PetscPartitionerType)::PetscErrorCode
end

function PetscPartitionerGetType(arg1, arg2)
    @ccall libpetsc.PetscPartitionerGetType(arg1::PetscPartitioner, arg2::Ptr{PetscPartitionerType})::PetscErrorCode
end

function PetscPartitionerSetUp(arg1)
    @ccall libpetsc.PetscPartitionerSetUp(arg1::PetscPartitioner)::PetscErrorCode
end

function PetscPartitionerReset(arg1)
    @ccall libpetsc.PetscPartitionerReset(arg1::PetscPartitioner)::PetscErrorCode
end

function PetscPartitionerSetFromOptions(arg1)
    @ccall libpetsc.PetscPartitionerSetFromOptions(arg1::PetscPartitioner)::PetscErrorCode
end

function PetscPartitionerViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.PetscPartitionerViewFromOptions(arg1::PetscPartitioner, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function PetscPartitionerView(arg1, arg2)
    @ccall libpetsc.PetscPartitionerView(arg1::PetscPartitioner, arg2::PetscViewer)::PetscErrorCode
end

function PetscPartitionerPartition(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.PetscPartitionerPartition(arg1::PetscPartitioner, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::PetscSection, arg7::PetscSection, arg8::PetscSection, arg9::Ptr{IS})::PetscErrorCode
end

function PetscPartitionerShellSetPartition(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscPartitionerShellSetPartition(arg1::PetscPartitioner, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function PetscPartitionerShellSetRandom(arg1, arg2)
    @ccall libpetsc.PetscPartitionerShellSetRandom(arg1::PetscPartitioner, arg2::PetscBool)::PetscErrorCode
end

function PetscPartitionerShellGetRandom(arg1, arg2)
    @ccall libpetsc.PetscPartitionerShellGetRandom(arg1::PetscPartitioner, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscPartitionerMatPartitioningGetMatPartitioning(arg1, arg2)
    @ccall libpetsc.PetscPartitionerMatPartitioningGetMatPartitioning(arg1::PetscPartitioner, arg2::Ptr{MatPartitioning})::PetscErrorCode
end

mutable struct _p_DMPlexCellRefiner end

const DMPlexCellRefiner = Ptr{_p_DMPlexCellRefiner}

@enum DMPlexCellRefinerType::UInt32 begin
    DM_REFINER_REGULAR = 0
    DM_REFINER_TO_BOX = 1
    DM_REFINER_TO_SIMPLEX = 2
    DM_REFINER_ALFELD2D = 3
    DM_REFINER_ALFELD3D = 4
    DM_REFINER_POWELL_SABIN = 5
    DM_REFINER_BOUNDARYLAYER = 6
    DM_REFINER_SBR = 7
end

mutable struct _p_PetscLimiter end

const PetscLimiter = Ptr{_p_PetscLimiter}

mutable struct _p_PetscFV end

const PetscFV = Ptr{_p_PetscFV}

mutable struct PetscFVFaceGeom
    normal::NTuple{3, PetscReal}
    centroid::NTuple{3, PetscReal}
    grad::NTuple{2, NTuple{3, PetscScalar}}
    PetscFVFaceGeom() = new()
end

mutable struct PetscFVCellGeom
    centroid::NTuple{3, PetscReal}
    volume::PetscReal
    PetscFVCellGeom() = new()
end

const PetscLimiterType = Ptr{Cchar}

function PetscLimiterCreate(arg1, arg2)
    @ccall libpetsc.PetscLimiterCreate(arg1::MPI_Comm, arg2::Ptr{PetscLimiter})::PetscErrorCode
end

function PetscLimiterDestroy(arg1)
    @ccall libpetsc.PetscLimiterDestroy(arg1::Ptr{PetscLimiter})::PetscErrorCode
end

function PetscLimiterSetType(arg1, arg2)
    @ccall libpetsc.PetscLimiterSetType(arg1::PetscLimiter, arg2::PetscLimiterType)::PetscErrorCode
end

function PetscLimiterGetType(arg1, arg2)
    @ccall libpetsc.PetscLimiterGetType(arg1::PetscLimiter, arg2::Ptr{PetscLimiterType})::PetscErrorCode
end

function PetscLimiterSetUp(arg1)
    @ccall libpetsc.PetscLimiterSetUp(arg1::PetscLimiter)::PetscErrorCode
end

function PetscLimiterSetFromOptions(arg1)
    @ccall libpetsc.PetscLimiterSetFromOptions(arg1::PetscLimiter)::PetscErrorCode
end

function PetscLimiterViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.PetscLimiterViewFromOptions(arg1::PetscLimiter, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function PetscLimiterView(arg1, arg2)
    @ccall libpetsc.PetscLimiterView(arg1::PetscLimiter, arg2::PetscViewer)::PetscErrorCode
end

function PetscLimiterRegister(arg1, arg2)
    @ccall libpetsc.PetscLimiterRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscLimiterRegisterDestroy()
    @ccall libpetsc.PetscLimiterRegisterDestroy()::PetscErrorCode
end

function PetscLimiterLimit(arg1, arg2, arg3)
    @ccall libpetsc.PetscLimiterLimit(arg1::PetscLimiter, arg2::PetscReal, arg3::Ptr{PetscReal})::PetscErrorCode
end

function PetscFVInitializePackage()
    @ccall libpetsc.PetscFVInitializePackage()::PetscErrorCode
end

const PetscFVType = Ptr{Cchar}

function PetscFVCreate(arg1, arg2)
    @ccall libpetsc.PetscFVCreate(arg1::MPI_Comm, arg2::Ptr{PetscFV})::PetscErrorCode
end

function PetscFVDestroy(arg1)
    @ccall libpetsc.PetscFVDestroy(arg1::Ptr{PetscFV})::PetscErrorCode
end

function PetscFVSetType(arg1, arg2)
    @ccall libpetsc.PetscFVSetType(arg1::PetscFV, arg2::PetscFVType)::PetscErrorCode
end

function PetscFVGetType(arg1, arg2)
    @ccall libpetsc.PetscFVGetType(arg1::PetscFV, arg2::Ptr{PetscFVType})::PetscErrorCode
end

function PetscFVSetUp(arg1)
    @ccall libpetsc.PetscFVSetUp(arg1::PetscFV)::PetscErrorCode
end

function PetscFVSetFromOptions(arg1)
    @ccall libpetsc.PetscFVSetFromOptions(arg1::PetscFV)::PetscErrorCode
end

function PetscFVViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.PetscFVViewFromOptions(arg1::PetscFV, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function PetscFVView(arg1, arg2)
    @ccall libpetsc.PetscFVView(arg1::PetscFV, arg2::PetscViewer)::PetscErrorCode
end

function PetscFVRegister(arg1, arg2)
    @ccall libpetsc.PetscFVRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscFVRegisterDestroy()
    @ccall libpetsc.PetscFVRegisterDestroy()::PetscErrorCode
end

function PetscFVSetComponentName(arg1, arg2, arg3)
    @ccall libpetsc.PetscFVSetComponentName(arg1::PetscFV, arg2::PetscInt, arg3::Ptr{Cchar})::PetscErrorCode
end

function PetscFVGetComponentName(arg1, arg2, arg3)
    @ccall libpetsc.PetscFVGetComponentName(arg1::PetscFV, arg2::PetscInt, arg3::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PetscFVSetLimiter(arg1, arg2)
    @ccall libpetsc.PetscFVSetLimiter(arg1::PetscFV, arg2::PetscLimiter)::PetscErrorCode
end

function PetscFVGetLimiter(arg1, arg2)
    @ccall libpetsc.PetscFVGetLimiter(arg1::PetscFV, arg2::Ptr{PetscLimiter})::PetscErrorCode
end

function PetscFVSetNumComponents(arg1, arg2)
    @ccall libpetsc.PetscFVSetNumComponents(arg1::PetscFV, arg2::PetscInt)::PetscErrorCode
end

function PetscFVGetNumComponents(arg1, arg2)
    @ccall libpetsc.PetscFVGetNumComponents(arg1::PetscFV, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscFVSetSpatialDimension(arg1, arg2)
    @ccall libpetsc.PetscFVSetSpatialDimension(arg1::PetscFV, arg2::PetscInt)::PetscErrorCode
end

function PetscFVGetSpatialDimension(arg1, arg2)
    @ccall libpetsc.PetscFVGetSpatialDimension(arg1::PetscFV, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscFVSetComputeGradients(arg1, arg2)
    @ccall libpetsc.PetscFVSetComputeGradients(arg1::PetscFV, arg2::PetscBool)::PetscErrorCode
end

function PetscFVGetComputeGradients(arg1, arg2)
    @ccall libpetsc.PetscFVGetComputeGradients(arg1::PetscFV, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscFVSetQuadrature(arg1, arg2)
    @ccall libpetsc.PetscFVSetQuadrature(arg1::PetscFV, arg2::PetscQuadrature)::PetscErrorCode
end

function PetscFVGetQuadrature(arg1, arg2)
    @ccall libpetsc.PetscFVGetQuadrature(arg1::PetscFV, arg2::Ptr{PetscQuadrature})::PetscErrorCode
end

function PetscFVSetDualSpace(arg1, arg2)
    @ccall libpetsc.PetscFVSetDualSpace(arg1::PetscFV, arg2::PetscDualSpace)::PetscErrorCode
end

function PetscFVGetDualSpace(arg1, arg2)
    @ccall libpetsc.PetscFVGetDualSpace(arg1::PetscFV, arg2::Ptr{PetscDualSpace})::PetscErrorCode
end

function PetscFVRefine(arg1, arg2)
    @ccall libpetsc.PetscFVRefine(arg1::PetscFV, arg2::Ptr{PetscFV})::PetscErrorCode
end

function PetscFVGetCellTabulation(arg1, arg2)
    @ccall libpetsc.PetscFVGetCellTabulation(arg1::PetscFV, arg2::Ptr{PetscTabulation})::PetscErrorCode
end

function PetscFVCreateTabulation(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscFVCreateTabulation(arg1::PetscFV, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscReal}, arg5::PetscInt, arg6::Ptr{PetscTabulation})::PetscErrorCode
end

function PetscFVComputeGradient(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscFVComputeGradient(arg1::PetscFV, arg2::PetscInt, arg3::Ptr{PetscScalar}, arg4::Ptr{PetscScalar})::PetscErrorCode
end

function PetscFVIntegrateRHSFunction(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.PetscFVIntegrateRHSFunction(arg1::PetscFV, arg2::PetscDS, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscFVFaceGeom}, arg6::Ptr{PetscReal}, arg7::Ptr{PetscScalar}, arg8::Ptr{PetscScalar}, arg9::Ptr{PetscScalar}, arg10::Ptr{PetscScalar})::PetscErrorCode
end

function PetscFVLeastSquaresSetMaxFaces(arg1, arg2)
    @ccall libpetsc.PetscFVLeastSquaresSetMaxFaces(arg1::PetscFV, arg2::PetscInt)::PetscErrorCode
end

function PetscDualSpaceApplyFVM(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscDualSpaceApplyFVM(arg1::PetscDualSpace, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscFVCellGeom}, arg5::PetscInt, arg6::Ptr{Cvoid}, arg7::Ptr{Cvoid}, arg8::Ptr{PetscScalar})::PetscErrorCode
end

function DMFieldInitializePackage()
    @ccall libpetsc.DMFieldInitializePackage()::PetscErrorCode
end

function DMFieldFinalizePackage()
    @ccall libpetsc.DMFieldFinalizePackage()::PetscErrorCode
end

const DMFieldType = Ptr{Cchar}

function DMFieldSetType(arg1, arg2)
    @ccall libpetsc.DMFieldSetType(arg1::DMField, arg2::DMFieldType)::PetscErrorCode
end

function DMFieldGetType(arg1, arg2)
    @ccall libpetsc.DMFieldGetType(arg1::DMField, arg2::Ptr{DMFieldType})::PetscErrorCode
end

function DMFieldRegister(arg1, arg2)
    @ccall libpetsc.DMFieldRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

@enum DMFieldContinuity::UInt32 begin
    DMFIELD_VERTEX = 0
    DMFIELD_EDGE = 1
    DMFIELD_FACET = 2
    DMFIELD_CELL = 3
end

function DMFieldDestroy(arg1)
    @ccall libpetsc.DMFieldDestroy(arg1::Ptr{DMField})::PetscErrorCode
end

function DMFieldView(arg1, arg2)
    @ccall libpetsc.DMFieldView(arg1::DMField, arg2::PetscViewer)::PetscErrorCode
end

function DMFieldGetDM(arg1, arg2)
    @ccall libpetsc.DMFieldGetDM(arg1::DMField, arg2::Ptr{DM})::PetscErrorCode
end

function DMFieldGetNumComponents(arg1, arg2)
    @ccall libpetsc.DMFieldGetNumComponents(arg1::DMField, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMFieldGetContinuity(arg1, arg2)
    @ccall libpetsc.DMFieldGetContinuity(arg1::DMField, arg2::Ptr{DMFieldContinuity})::PetscErrorCode
end

function DMFieldEvaluate(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMFieldEvaluate(arg1::DMField, arg2::Vec, arg3::PetscDataType, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid}, arg6::Ptr{Cvoid})::PetscErrorCode
end

function DMFieldEvaluateFE(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMFieldEvaluateFE(arg1::DMField, arg2::IS, arg3::PetscQuadrature, arg4::PetscDataType, arg5::Ptr{Cvoid}, arg6::Ptr{Cvoid}, arg7::Ptr{Cvoid})::PetscErrorCode
end

function DMFieldEvaluateFV(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMFieldEvaluateFV(arg1::DMField, arg2::IS, arg3::PetscDataType, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid}, arg6::Ptr{Cvoid})::PetscErrorCode
end

function DMFieldCreateFEGeom(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMFieldCreateFEGeom(arg1::DMField, arg2::IS, arg3::PetscQuadrature, arg4::PetscBool, arg5::Ptr{Ptr{PetscFEGeom}})::PetscErrorCode
end

function DMFieldCreateDefaultQuadrature(arg1, arg2, arg3)
    @ccall libpetsc.DMFieldCreateDefaultQuadrature(arg1::DMField, arg2::IS, arg3::Ptr{PetscQuadrature})::PetscErrorCode
end

function DMFieldGetDegree(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMFieldGetDegree(arg1::DMField, arg2::IS, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMFieldCreateDA(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMFieldCreateDA(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscScalar}, arg4::Ptr{DMField})::PetscErrorCode
end

function DMFieldCreateDS(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMFieldCreateDS(arg1::DM, arg2::PetscInt, arg3::Vec, arg4::Ptr{DMField})::PetscErrorCode
end

function DMFieldCreateShell(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMFieldCreateShell(arg1::DM, arg2::PetscInt, arg3::DMFieldContinuity, arg4::Ptr{Cvoid}, arg5::Ptr{DMField})::PetscErrorCode
end

function DMFieldShellSetDestroy(arg1, arg2)
    @ccall libpetsc.DMFieldShellSetDestroy(arg1::DMField, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMFieldShellGetContext(arg1, arg2)
    @ccall libpetsc.DMFieldShellGetContext(arg1::DMField, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMFieldShellSetEvaluate(arg1, arg2)
    @ccall libpetsc.DMFieldShellSetEvaluate(arg1::DMField, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMFieldShellSetEvaluateFE(arg1, arg2)
    @ccall libpetsc.DMFieldShellSetEvaluateFE(arg1::DMField, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMFieldShellEvaluateFEDefault(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMFieldShellEvaluateFEDefault(arg1::DMField, arg2::IS, arg3::PetscQuadrature, arg4::PetscDataType, arg5::Ptr{Cvoid}, arg6::Ptr{Cvoid}, arg7::Ptr{Cvoid})::PetscErrorCode
end

function DMFieldShellSetEvaluateFV(arg1, arg2)
    @ccall libpetsc.DMFieldShellSetEvaluateFV(arg1::DMField, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMFieldShellEvaluateFVDefault(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMFieldShellEvaluateFVDefault(arg1::DMField, arg2::IS, arg3::PetscDataType, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid}, arg6::Ptr{Cvoid})::PetscErrorCode
end

function DMFieldShellSetGetDegree(arg1, arg2)
    @ccall libpetsc.DMFieldShellSetGetDegree(arg1::DMField, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMFieldShellSetCreateDefaultQuadrature(arg1, arg2)
    @ccall libpetsc.DMFieldShellSetCreateDefaultQuadrature(arg1::DMField, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscPartitionerDMPlexPartition(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscPartitionerDMPlexPartition(arg1::PetscPartitioner, arg2::DM, arg3::PetscSection, arg4::PetscSection, arg5::Ptr{IS})::PetscErrorCode
end

function DMPlexBuildFromCellList(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexBuildFromCellList(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexBuildFromCellListParallel(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMPlexBuildFromCellListParallel(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::Ptr{PetscSF})::PetscErrorCode
end

function DMPlexBuildCoordinatesFromCellList(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexBuildCoordinatesFromCellList(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscReal})::PetscErrorCode
end

function DMPlexBuildCoordinatesFromCellListParallel(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexBuildCoordinatesFromCellListParallel(arg1::DM, arg2::PetscInt, arg3::PetscSF, arg4::Ptr{PetscReal})::PetscErrorCode
end

function DMPlexCreate(arg1, arg2)
    @ccall libpetsc.DMPlexCreate(arg1::MPI_Comm, arg2::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateCohesiveSubmesh(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexCreateCohesiveSubmesh(arg1::DM, arg2::PetscBool, arg3::Ptr{Cchar}, arg4::PetscInt, arg5::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateFromCellListPetsc(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.DMPlexCreateFromCellListPetsc(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscBool, arg7::Ptr{PetscInt}, arg8::PetscInt, arg9::Ptr{PetscReal}, arg10::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateFromCellList(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.DMPlexCreateFromCellList(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscBool, arg7::Ptr{Cint}, arg8::PetscInt, arg9::Ptr{Cdouble}, arg10::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateFromCellListParallelPetsc(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12)
    @ccall libpetsc.DMPlexCreateFromCellListParallelPetsc(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::PetscBool, arg8::Ptr{PetscInt}, arg9::PetscInt, arg10::Ptr{PetscReal}, arg11::Ptr{PetscSF}, arg12::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateFromCellListParallel(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11)
    @ccall libpetsc.DMPlexCreateFromCellListParallel(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscBool, arg7::Ptr{Cint}, arg8::PetscInt, arg9::Ptr{PetscReal}, arg10::Ptr{PetscSF}, arg11::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateFromDAG(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMPlexCreateFromDAG(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscInt}, arg7::Ptr{PetscScalar})::PetscErrorCode
end

function DMPlexCreateReferenceCell(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexCreateReferenceCell(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscBool, arg4::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateReferenceCellByType(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexCreateReferenceCellByType(arg1::MPI_Comm, arg2::DMPolytopeType, arg3::Ptr{DM})::PetscErrorCode
end

function DMPlexSetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.DMPlexSetOptionsPrefix(arg1::DM, arg2::Ptr{Cchar})::PetscErrorCode
end

function DMPlexGetChart(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexGetChart(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexSetChart(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexSetChart(arg1::DM, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function DMPlexGetConeSize(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexGetConeSize(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexSetConeSize(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexSetConeSize(arg1::DM, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function DMPlexAddConeSize(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexAddConeSize(arg1::DM, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function DMPlexGetCone(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexGetCone(arg1::DM, arg2::PetscInt, arg3::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMPlexGetConeTuple(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexGetConeTuple(arg1::DM, arg2::IS, arg3::Ptr{PetscSection}, arg4::Ptr{IS})::PetscErrorCode
end

function DMPlexGetConeRecursive(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexGetConeRecursive(arg1::DM, arg2::IS, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{IS}}, arg5::Ptr{Ptr{PetscSection}})::PetscErrorCode
end

function DMPlexRestoreConeRecursive(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexRestoreConeRecursive(arg1::DM, arg2::IS, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{IS}}, arg5::Ptr{Ptr{PetscSection}})::PetscErrorCode
end

function DMPlexGetConeRecursiveVertices(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexGetConeRecursiveVertices(arg1::DM, arg2::IS, arg3::Ptr{IS})::PetscErrorCode
end

function DMPlexSetCone(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexSetCone(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexInsertCone(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexInsertCone(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt)::PetscErrorCode
end

function DMPlexInsertConeOrientation(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexInsertConeOrientation(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt)::PetscErrorCode
end

function DMPlexGetConeOrientation(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexGetConeOrientation(arg1::DM, arg2::PetscInt, arg3::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMPlexSetConeOrientation(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexSetConeOrientation(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexGetSupportSize(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexGetSupportSize(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexSetSupportSize(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexSetSupportSize(arg1::DM, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function DMPlexGetSupport(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexGetSupport(arg1::DM, arg2::PetscInt, arg3::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMPlexSetSupport(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexSetSupport(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexInsertSupport(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexInsertSupport(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt)::PetscErrorCode
end

function DMPlexGetConeSection(arg1, arg2)
    @ccall libpetsc.DMPlexGetConeSection(arg1::DM, arg2::Ptr{PetscSection})::PetscErrorCode
end

function DMPlexGetSupportSection(arg1, arg2)
    @ccall libpetsc.DMPlexGetSupportSection(arg1::DM, arg2::Ptr{PetscSection})::PetscErrorCode
end

function DMPlexGetCones(arg1, arg2)
    @ccall libpetsc.DMPlexGetCones(arg1::DM, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMPlexGetConeOrientations(arg1, arg2)
    @ccall libpetsc.DMPlexGetConeOrientations(arg1::DM, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMPlexGetMaxSizes(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexGetMaxSizes(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexSymmetrize(arg1)
    @ccall libpetsc.DMPlexSymmetrize(arg1::DM)::PetscErrorCode
end

function DMPlexStratify(arg1)
    @ccall libpetsc.DMPlexStratify(arg1::DM)::PetscErrorCode
end

function DMPlexEqual(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexEqual(arg1::DM, arg2::DM, arg3::Ptr{PetscBool})::PetscErrorCode
end

function DMPlexReverseCell(arg1, arg2)
    @ccall libpetsc.DMPlexReverseCell(arg1::DM, arg2::PetscInt)::PetscErrorCode
end

function DMPlexOrientCell(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexOrientCell(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexCompareOrientations(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMPlexCompareOrientations(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscBool})::PetscErrorCode
end

function DMPlexOrient(arg1)
    @ccall libpetsc.DMPlexOrient(arg1::DM)::PetscErrorCode
end

function DMPlexPreallocateOperator(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.DMPlexPreallocateOperator(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscInt}, arg7::Mat, arg8::PetscBool)::PetscErrorCode
end

function DMPlexGetPointLocal(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexGetPointLocal(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexPointLocalRead(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexPointLocalRead(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscScalar}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexPointLocalRef(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexPointLocalRef(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscScalar}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexGetPointLocalField(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexGetPointLocalField(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexPointLocalFieldRef(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexPointLocalFieldRef(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscScalar}, arg5::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexPointLocalFieldRead(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexPointLocalFieldRead(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscScalar}, arg5::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexGetPointGlobal(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexGetPointGlobal(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexPointGlobalRead(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexPointGlobalRead(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscScalar}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexPointGlobalRef(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexPointGlobalRef(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscScalar}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexGetPointGlobalField(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexGetPointGlobalField(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexPointGlobalFieldRef(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexPointGlobalFieldRef(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscScalar}, arg5::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexPointGlobalFieldRead(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexPointGlobalFieldRead(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscScalar}, arg5::Ptr{Cvoid})::PetscErrorCode
end

@enum DMPlexInterpolatedFlag::Int32 begin
    DMPLEX_INTERPOLATED_INVALID = -1
    DMPLEX_INTERPOLATED_NONE = 0
    DMPLEX_INTERPOLATED_PARTIAL = 1
    DMPLEX_INTERPOLATED_MIXED = 2
    DMPLEX_INTERPOLATED_FULL = 3
end

function DMPlexInterpolate(arg1, arg2)
    @ccall libpetsc.DMPlexInterpolate(arg1::DM, arg2::Ptr{DM})::PetscErrorCode
end

function DMPlexUninterpolate(arg1, arg2)
    @ccall libpetsc.DMPlexUninterpolate(arg1::DM, arg2::Ptr{DM})::PetscErrorCode
end

function DMPlexInterpolatePointSF(arg1, arg2)
    @ccall libpetsc.DMPlexInterpolatePointSF(arg1::DM, arg2::PetscSF)::PetscErrorCode
end

function DMPlexIsInterpolated(arg1, arg2)
    @ccall libpetsc.DMPlexIsInterpolated(arg1::DM, arg2::Ptr{DMPlexInterpolatedFlag})::PetscErrorCode
end

function DMPlexIsInterpolatedCollective(arg1, arg2)
    @ccall libpetsc.DMPlexIsInterpolatedCollective(arg1::DM, arg2::Ptr{DMPlexInterpolatedFlag})::PetscErrorCode
end

function DMPlexFilter(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexFilter(arg1::DM, arg2::DMLabel, arg3::PetscInt, arg4::Ptr{DM})::PetscErrorCode
end

function DMPlexGetCellNumbering(arg1, arg2)
    @ccall libpetsc.DMPlexGetCellNumbering(arg1::DM, arg2::Ptr{IS})::PetscErrorCode
end

function DMPlexGetVertexNumbering(arg1, arg2)
    @ccall libpetsc.DMPlexGetVertexNumbering(arg1::DM, arg2::Ptr{IS})::PetscErrorCode
end

function DMPlexCreatePointNumbering(arg1, arg2)
    @ccall libpetsc.DMPlexCreatePointNumbering(arg1::DM, arg2::Ptr{IS})::PetscErrorCode
end

function DMPlexCreateRankField(arg1, arg2)
    @ccall libpetsc.DMPlexCreateRankField(arg1::DM, arg2::Ptr{Vec})::PetscErrorCode
end

function DMPlexCreateLabelField(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexCreateLabelField(arg1::DM, arg2::DMLabel, arg3::Ptr{Vec})::PetscErrorCode
end

function DMPlexGetDepth(arg1, arg2)
    @ccall libpetsc.DMPlexGetDepth(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexGetDepthLabel(arg1, arg2)
    @ccall libpetsc.DMPlexGetDepthLabel(arg1::DM, arg2::Ptr{DMLabel})::PetscErrorCode
end

function DMPlexGetDepthStratum(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexGetDepthStratum(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexGetHeightStratum(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexGetHeightStratum(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexGetPointDepth(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexGetPointDepth(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexGetPointHeight(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexGetPointHeight(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexGetCellTypeLabel(arg1, arg2)
    @ccall libpetsc.DMPlexGetCellTypeLabel(arg1::DM, arg2::Ptr{DMLabel})::PetscErrorCode
end

function DMPlexGetCellType(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexGetCellType(arg1::DM, arg2::PetscInt, arg3::Ptr{DMPolytopeType})::PetscErrorCode
end

function DMPlexSetCellType(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexSetCellType(arg1::DM, arg2::PetscInt, arg3::DMPolytopeType)::PetscErrorCode
end

function DMPlexComputeCellTypes(arg1)
    @ccall libpetsc.DMPlexComputeCellTypes(arg1::DM)::PetscErrorCode
end

function DMPlexInvertCell(arg1, arg2)
    @ccall libpetsc.DMPlexInvertCell(arg1::DMPolytopeType, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexReorderCell(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexReorderCell(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexGetMeet(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexGetMeet(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMPlexGetFullMeet(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexGetFullMeet(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMPlexRestoreMeet(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexRestoreMeet(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMPlexGetJoin(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexGetJoin(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMPlexGetFullJoin(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexGetFullJoin(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMPlexRestoreJoin(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexRestoreJoin(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMPlexGetTransitiveClosure(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexGetTransitiveClosure(arg1::DM, arg2::PetscInt, arg3::PetscBool, arg4::Ptr{PetscInt}, arg5::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMPlexRestoreTransitiveClosure(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexRestoreTransitiveClosure(arg1::DM, arg2::PetscInt, arg3::PetscBool, arg4::Ptr{PetscInt}, arg5::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMPlexGenerate(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexGenerate(arg1::DM, arg2::Ptr{Cchar}, arg3::PetscBool, arg4::Ptr{DM})::PetscErrorCode
end

function DMPlexGenerateRegister(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexGenerateRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid}, arg5::PetscInt)::PetscErrorCode
end

function DMPlexGenerateRegisterAll()
    @ccall libpetsc.DMPlexGenerateRegisterAll()::PetscErrorCode
end

function DMPlexGenerateRegisterDestroy()
    @ccall libpetsc.DMPlexGenerateRegisterDestroy()::PetscErrorCode
end

function DMPlexCopyCoordinates(arg1, arg2)
    @ccall libpetsc.DMPlexCopyCoordinates(arg1::DM, arg2::DM)::PetscErrorCode
end

function DMPlexCreateDoublet(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMPlexCreateDoublet(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscBool, arg4::PetscBool, arg5::PetscReal, arg6::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateSquareBoundary(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexCreateSquareBoundary(arg1::DM, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexCreateCubeBoundary(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexCreateCubeBoundary(arg1::DM, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexCreateBoxMesh(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.DMPlexCreateBoxMesh(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscBool, arg4::Ptr{PetscInt}, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal}, arg7::Ptr{DMBoundaryType}, arg8::PetscBool, arg9::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateSphereMesh(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexCreateSphereMesh(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscBool, arg4::PetscReal, arg5::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateBallMesh(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexCreateBallMesh(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateHexCylinderMesh(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexCreateHexCylinderMesh(arg1::MPI_Comm, arg2::PetscInt, arg3::DMBoundaryType, arg4::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateWedgeCylinderMesh(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexCreateWedgeCylinderMesh(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscBool, arg4::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateWedgeBoxMesh(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.DMPlexCreateWedgeBoxMesh(arg1::MPI_Comm, arg2::Ptr{PetscInt}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{DMBoundaryType}, arg6::PetscBool, arg7::PetscBool, arg8::Ptr{DM})::PetscErrorCode
end

function DMPlexExtrude(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMPlexExtrude(arg1::DM, arg2::PetscInt, arg3::PetscReal, arg4::PetscBool, arg5::Ptr{PetscReal}, arg6::PetscBool, arg7::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateConeSection(arg1, arg2)
    @ccall libpetsc.DMPlexCreateConeSection(arg1::DM, arg2::Ptr{PetscSection})::PetscErrorCode
end

function DMPlexCheckSymmetry(arg1)
    @ccall libpetsc.DMPlexCheckSymmetry(arg1::DM)::PetscErrorCode
end

function DMPlexCheckSkeleton(arg1, arg2)
    @ccall libpetsc.DMPlexCheckSkeleton(arg1::DM, arg2::PetscInt)::PetscErrorCode
end

function DMPlexCheckFaces(arg1, arg2)
    @ccall libpetsc.DMPlexCheckFaces(arg1::DM, arg2::PetscInt)::PetscErrorCode
end

function DMPlexCheckGeometry(arg1)
    @ccall libpetsc.DMPlexCheckGeometry(arg1::DM)::PetscErrorCode
end

function DMPlexCheckPointSF(arg1)
    @ccall libpetsc.DMPlexCheckPointSF(arg1::DM)::PetscErrorCode
end

function DMPlexCheckInterfaceCones(arg1)
    @ccall libpetsc.DMPlexCheckInterfaceCones(arg1::DM)::PetscErrorCode
end

function DMPlexCheckCellShape(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexCheckCellShape(arg1::DM, arg2::PetscBool, arg3::PetscReal)::PetscErrorCode
end

function DMPlexComputeOrthogonalQuality(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexComputeOrthogonalQuality(arg1::DM, arg2::PetscFV, arg3::PetscReal, arg4::Ptr{Vec}, arg5::Ptr{DMLabel})::PetscErrorCode
end

function DMPlexTriangleSetOptions(arg1, arg2)
    @ccall libpetsc.DMPlexTriangleSetOptions(arg1::DM, arg2::Ptr{Cchar})::PetscErrorCode
end

function DMPlexTetgenSetOptions(arg1, arg2)
    @ccall libpetsc.DMPlexTetgenSetOptions(arg1::DM, arg2::Ptr{Cchar})::PetscErrorCode
end

function DMPlexCreateFromFile(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexCreateFromFile(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::PetscBool, arg4::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateExodus(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexCreateExodus(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscBool, arg4::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateExodusFromFile(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexCreateExodusFromFile(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::PetscBool, arg4::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateCGNS(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexCreateCGNS(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscBool, arg4::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateCGNSFromFile(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexCreateCGNSFromFile(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::PetscBool, arg4::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateGmsh(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexCreateGmsh(arg1::MPI_Comm, arg2::PetscViewer, arg3::PetscBool, arg4::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateGmshFromFile(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexCreateGmshFromFile(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::PetscBool, arg4::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateFluent(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexCreateFluent(arg1::MPI_Comm, arg2::PetscViewer, arg3::PetscBool, arg4::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateFluentFromFile(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexCreateFluentFromFile(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::PetscBool, arg4::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateMedFromFile(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexCreateMedFromFile(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::PetscBool, arg4::Ptr{DM})::PetscErrorCode
end

function DMPlexCreatePLYFromFile(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexCreatePLYFromFile(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::PetscBool, arg4::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateEGADSFromFile(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexCreateEGADSFromFile(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{DM})::PetscErrorCode
end

function PetscViewerExodusIIGetId(arg1, arg2)
    @ccall libpetsc.PetscViewerExodusIIGetId(arg1::PetscViewer, arg2::Ptr{Cint})::PetscErrorCode
end

function DMPlexCreateNeighborCSR(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexCreateNeighborCSR(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{PetscInt}}, arg5::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMPlexGetPartitioner(arg1, arg2)
    @ccall libpetsc.DMPlexGetPartitioner(arg1::DM, arg2::Ptr{PetscPartitioner})::PetscErrorCode
end

function DMPlexSetPartitioner(arg1, arg2)
    @ccall libpetsc.DMPlexSetPartitioner(arg1::DM, arg2::PetscPartitioner)::PetscErrorCode
end

function DMPlexCreatePartitionerGraph(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMPlexCreatePartitionerGraph(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{PetscInt}}, arg5::Ptr{Ptr{PetscInt}}, arg6::Ptr{IS})::PetscErrorCode
end

function DMPlexPartitionLabelInvert(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexPartitionLabelInvert(arg1::DM, arg2::DMLabel, arg3::PetscSF, arg4::DMLabel)::PetscErrorCode
end

function DMPlexPartitionLabelClosure(arg1, arg2)
    @ccall libpetsc.DMPlexPartitionLabelClosure(arg1::DM, arg2::DMLabel)::PetscErrorCode
end

function DMPlexPartitionLabelAdjacency(arg1, arg2)
    @ccall libpetsc.DMPlexPartitionLabelAdjacency(arg1::DM, arg2::DMLabel)::PetscErrorCode
end

function DMPlexPartitionLabelPropagate(arg1, arg2)
    @ccall libpetsc.DMPlexPartitionLabelPropagate(arg1::DM, arg2::DMLabel)::PetscErrorCode
end

function DMPlexPartitionLabelCreateSF(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexPartitionLabelCreateSF(arg1::DM, arg2::DMLabel, arg3::Ptr{PetscSF})::PetscErrorCode
end

function DMPlexSetPartitionBalance(arg1, arg2)
    @ccall libpetsc.DMPlexSetPartitionBalance(arg1::DM, arg2::PetscBool)::PetscErrorCode
end

function DMPlexGetPartitionBalance(arg1, arg2)
    @ccall libpetsc.DMPlexGetPartitionBalance(arg1::DM, arg2::Ptr{PetscBool})::PetscErrorCode
end

function DMPlexIsDistributed(arg1, arg2)
    @ccall libpetsc.DMPlexIsDistributed(arg1::DM, arg2::Ptr{PetscBool})::PetscErrorCode
end

function DMPlexDistribute(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexDistribute(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscSF}, arg4::Ptr{DM})::PetscErrorCode
end

function DMPlexDistributeOverlap(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexDistributeOverlap(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscSF}, arg4::Ptr{DM})::PetscErrorCode
end

function DMPlexGetOverlap(arg1, arg2)
    @ccall libpetsc.DMPlexGetOverlap(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexDistributeField(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMPlexDistributeField(arg1::DM, arg2::PetscSF, arg3::PetscSection, arg4::Vec, arg5::PetscSection, arg6::Vec)::PetscErrorCode
end

function DMPlexDistributeFieldIS(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMPlexDistributeFieldIS(arg1::DM, arg2::PetscSF, arg3::PetscSection, arg4::IS, arg5::PetscSection, arg6::Ptr{IS})::PetscErrorCode
end

function DMPlexDistributeData(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMPlexDistributeData(arg1::DM, arg2::PetscSF, arg3::PetscSection, arg4::MPI_Datatype, arg5::Ptr{Cvoid}, arg6::PetscSection, arg7::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMPlexRebalanceSharedPoints(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexRebalanceSharedPoints(arg1::DM, arg2::PetscInt, arg3::PetscBool, arg4::PetscBool, arg5::Ptr{PetscBool})::PetscErrorCode
end

function DMPlexMigrate(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexMigrate(arg1::DM, arg2::PetscSF, arg3::DM)::PetscErrorCode
end

function DMPlexGetGatherDM(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexGetGatherDM(arg1::DM, arg2::Ptr{PetscSF}, arg3::Ptr{DM})::PetscErrorCode
end

function DMPlexGetRedundantDM(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexGetRedundantDM(arg1::DM, arg2::Ptr{PetscSF}, arg3::Ptr{DM})::PetscErrorCode
end

function DMPlexSetAdjacencyUser(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexSetAdjacencyUser(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexGetAdjacencyUser(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexGetAdjacencyUser(arg1::DM, arg2::Ptr{Ptr{Cvoid}}, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMPlexSetAdjacencyUseAnchors(arg1, arg2)
    @ccall libpetsc.DMPlexSetAdjacencyUseAnchors(arg1::DM, arg2::PetscBool)::PetscErrorCode
end

function DMPlexGetAdjacencyUseAnchors(arg1, arg2)
    @ccall libpetsc.DMPlexGetAdjacencyUseAnchors(arg1::DM, arg2::Ptr{PetscBool})::PetscErrorCode
end

function DMPlexGetAdjacency(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexGetAdjacency(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMPlexSetMigrationSF(arg1, arg2)
    @ccall libpetsc.DMPlexSetMigrationSF(arg1::DM, arg2::PetscSF)::PetscErrorCode
end

function DMPlexGetMigrationSF(arg1, arg2)
    @ccall libpetsc.DMPlexGetMigrationSF(arg1::DM, arg2::Ptr{PetscSF})::PetscErrorCode
end

function DMPlexGetOrdering(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexGetOrdering(arg1::DM, arg2::MatOrderingType, arg3::DMLabel, arg4::Ptr{IS})::PetscErrorCode
end

function DMPlexPermute(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexPermute(arg1::DM, arg2::IS, arg3::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateProcessSF(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexCreateProcessSF(arg1::DM, arg2::PetscSF, arg3::Ptr{IS}, arg4::Ptr{PetscSF})::PetscErrorCode
end

function DMPlexCreateTwoSidedProcessSF(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.DMPlexCreateTwoSidedProcessSF(arg1::DM, arg2::PetscSF, arg3::PetscSection, arg4::IS, arg5::PetscSection, arg6::IS, arg7::Ptr{IS}, arg8::Ptr{PetscSF})::PetscErrorCode
end

function DMPlexDistributeOwnership(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexDistributeOwnership(arg1::DM, arg2::PetscSection, arg3::Ptr{IS}, arg4::PetscSection, arg5::Ptr{IS})::PetscErrorCode
end

function DMPlexCreateOverlapLabel(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMPlexCreateOverlapLabel(arg1::DM, arg2::PetscInt, arg3::PetscSection, arg4::IS, arg5::PetscSection, arg6::IS, arg7::Ptr{DMLabel})::PetscErrorCode
end

function DMPlexCreateOverlapMigrationSF(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexCreateOverlapMigrationSF(arg1::DM, arg2::PetscSF, arg3::Ptr{PetscSF})::PetscErrorCode
end

function DMPlexStratifyMigrationSF(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexStratifyMigrationSF(arg1::DM, arg2::PetscSF, arg3::Ptr{PetscSF})::PetscErrorCode
end

function DMPlexCreateSubmesh(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexCreateSubmesh(arg1::DM, arg2::DMLabel, arg3::PetscInt, arg4::PetscBool, arg5::Ptr{DM})::PetscErrorCode
end

function DMPlexCreateHybridMesh(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMPlexCreateHybridMesh(arg1::DM, arg2::DMLabel, arg3::DMLabel, arg4::Ptr{DMLabel}, arg5::Ptr{DMLabel}, arg6::Ptr{DM}, arg7::Ptr{DM})::PetscErrorCode
end

function DMPlexGetSubpointMap(arg1, arg2)
    @ccall libpetsc.DMPlexGetSubpointMap(arg1::DM, arg2::Ptr{DMLabel})::PetscErrorCode
end

function DMPlexSetSubpointMap(arg1, arg2)
    @ccall libpetsc.DMPlexSetSubpointMap(arg1::DM, arg2::DMLabel)::PetscErrorCode
end

function DMPlexGetSubpointIS(arg1, arg2)
    @ccall libpetsc.DMPlexGetSubpointIS(arg1::DM, arg2::Ptr{IS})::PetscErrorCode
end

function DMGetEnclosureRelation(arg1, arg2, arg3)
    @ccall libpetsc.DMGetEnclosureRelation(arg1::DM, arg2::DM, arg3::Ptr{DMEnclosureType})::PetscErrorCode
end

function DMGetEnclosurePoint(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMGetEnclosurePoint(arg1::DM, arg2::DM, arg3::DMEnclosureType, arg4::PetscInt, arg5::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexLabelComplete(arg1, arg2)
    @ccall libpetsc.DMPlexLabelComplete(arg1::DM, arg2::DMLabel)::PetscErrorCode
end

function DMPlexLabelCohesiveComplete(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexLabelCohesiveComplete(arg1::DM, arg2::DMLabel, arg3::DMLabel, arg4::PetscBool, arg5::DM)::PetscErrorCode
end

function DMPlexLabelAddCells(arg1, arg2)
    @ccall libpetsc.DMPlexLabelAddCells(arg1::DM, arg2::DMLabel)::PetscErrorCode
end

function DMPlexLabelAddFaceCells(arg1, arg2)
    @ccall libpetsc.DMPlexLabelAddFaceCells(arg1::DM, arg2::DMLabel)::PetscErrorCode
end

function DMPlexLabelClearCells(arg1, arg2)
    @ccall libpetsc.DMPlexLabelClearCells(arg1::DM, arg2::DMLabel)::PetscErrorCode
end

function DMPlexGetRefinementLimit(arg1, arg2)
    @ccall libpetsc.DMPlexGetRefinementLimit(arg1::DM, arg2::Ptr{PetscReal})::PetscErrorCode
end

function DMPlexSetRefinementLimit(arg1, arg2)
    @ccall libpetsc.DMPlexSetRefinementLimit(arg1::DM, arg2::PetscReal)::PetscErrorCode
end

function DMPlexGetRefinementUniform(arg1, arg2)
    @ccall libpetsc.DMPlexGetRefinementUniform(arg1::DM, arg2::Ptr{PetscBool})::PetscErrorCode
end

function DMPlexSetRefinementUniform(arg1, arg2)
    @ccall libpetsc.DMPlexSetRefinementUniform(arg1::DM, arg2::PetscBool)::PetscErrorCode
end

function DMPlexGetRefinementFunction(arg1, arg2)
    @ccall libpetsc.DMPlexGetRefinementFunction(arg1::DM, arg2::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMPlexSetRefinementFunction(arg1, arg2)
    @ccall libpetsc.DMPlexSetRefinementFunction(arg1::DM, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexCreateCoarsePointIS(arg1, arg2)
    @ccall libpetsc.DMPlexCreateCoarsePointIS(arg1::DM, arg2::Ptr{IS})::PetscErrorCode
end

function DMPlexGetRegularRefinement(arg1, arg2)
    @ccall libpetsc.DMPlexGetRegularRefinement(arg1::DM, arg2::Ptr{PetscBool})::PetscErrorCode
end

function DMPlexSetRegularRefinement(arg1, arg2)
    @ccall libpetsc.DMPlexSetRegularRefinement(arg1::DM, arg2::PetscBool)::PetscErrorCode
end

function DMPlexGetCellRefinerType(arg1, arg2)
    @ccall libpetsc.DMPlexGetCellRefinerType(arg1::DM, arg2::Ptr{DMPlexCellRefinerType})::PetscErrorCode
end

function DMPlexSetCellRefinerType(arg1, arg2)
    @ccall libpetsc.DMPlexSetCellRefinerType(arg1::DM, arg2::DMPlexCellRefinerType)::PetscErrorCode
end

function DMPlexGetNumFaceVertices(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexGetNumFaceVertices(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexGetOrientedFace(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.DMPlexGetOrientedFace(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::Ptr{PetscInt}, arg8::Ptr{PetscInt}, arg9::Ptr{PetscBool})::PetscErrorCode
end

function DMPlexGetMinRadius(arg1, arg2)
    @ccall libpetsc.DMPlexGetMinRadius(arg1::DM, arg2::Ptr{PetscReal})::PetscErrorCode
end

function DMPlexSetMinRadius(arg1, arg2)
    @ccall libpetsc.DMPlexSetMinRadius(arg1::DM, arg2::PetscReal)::PetscErrorCode
end

function DMPlexComputeProjection2Dto1D(arg1, arg2)
    @ccall libpetsc.DMPlexComputeProjection2Dto1D(arg1::Ptr{PetscScalar}, arg2::Ptr{PetscReal})::PetscErrorCode
end

function DMPlexComputeProjection3Dto1D(arg1, arg2)
    @ccall libpetsc.DMPlexComputeProjection3Dto1D(arg1::Ptr{PetscScalar}, arg2::Ptr{PetscReal})::PetscErrorCode
end

function DMPlexComputeProjection3Dto2D(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexComputeProjection3Dto2D(arg1::PetscInt, arg2::Ptr{PetscScalar}, arg3::Ptr{PetscReal})::PetscErrorCode
end

mutable struct _PetscGridHash end

const PetscGridHash = Ptr{_PetscGridHash}

function PetscGridHashCreate(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscGridHashCreate(arg1::MPI_Comm, arg2::PetscInt, arg3::Ptr{PetscScalar}, arg4::Ptr{PetscGridHash})::PetscErrorCode
end

function PetscGridHashEnlarge(arg1, arg2)
    @ccall libpetsc.PetscGridHashEnlarge(arg1::PetscGridHash, arg2::Ptr{PetscScalar})::PetscErrorCode
end

function PetscGridHashSetGrid(arg1, arg2, arg3)
    @ccall libpetsc.PetscGridHashSetGrid(arg1::PetscGridHash, arg2::Ptr{PetscInt}, arg3::Ptr{PetscReal})::PetscErrorCode
end

function PetscGridHashGetEnclosingBox(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscGridHashGetEnclosingBox(arg1::PetscGridHash, arg2::PetscInt, arg3::Ptr{PetscScalar}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt})::PetscErrorCode
end

function PetscGridHashDestroy(arg1)
    @ccall libpetsc.PetscGridHashDestroy(arg1::Ptr{PetscGridHash})::PetscErrorCode
end

function DMPlexFindVertices(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexFindVertices(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscReal}, arg4::PetscReal, arg5::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexComputeCellGeometryFVM(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexComputeCellGeometryFVM(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal})::PetscErrorCode
end

function DMPlexComputeGeometryFVM(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexComputeGeometryFVM(arg1::DM, arg2::Ptr{Vec}, arg3::Ptr{Vec})::PetscErrorCode
end

function DMPlexComputeGradientFVM(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexComputeGradientFVM(arg1::DM, arg2::PetscFV, arg3::Vec, arg4::Vec, arg5::Ptr{DM})::PetscErrorCode
end

function DMPlexGetDataFVM(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexGetDataFVM(arg1::DM, arg2::PetscFV, arg3::Ptr{Vec}, arg4::Ptr{Vec}, arg5::Ptr{DM})::PetscErrorCode
end

function DMPlexComputeGeometryFEM(arg1, arg2)
    @ccall libpetsc.DMPlexComputeGeometryFEM(arg1::DM, arg2::Ptr{Vec})::PetscErrorCode
end

function DMPlexGetGeometryFVM(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexGetGeometryFVM(arg1::DM, arg2::Ptr{Vec}, arg3::Ptr{Vec}, arg4::Ptr{PetscReal})::PetscErrorCode
end

function DMPlexGetGradientDM(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexGetGradientDM(arg1::DM, arg2::PetscFV, arg3::Ptr{DM})::PetscErrorCode
end

function DMPlexInsertBoundaryValues(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMPlexInsertBoundaryValues(arg1::DM, arg2::PetscBool, arg3::Vec, arg4::PetscReal, arg5::Vec, arg6::Vec, arg7::Vec)::PetscErrorCode
end

function DMPlexInsertTimeDerivativeBoundaryValues(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMPlexInsertTimeDerivativeBoundaryValues(arg1::DM, arg2::PetscBool, arg3::Vec, arg4::PetscReal, arg5::Vec, arg6::Vec, arg7::Vec)::PetscErrorCode
end

function DMPlexInsertBoundaryValuesEssential(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11)
    @ccall libpetsc.DMPlexInsertBoundaryValuesEssential(arg1::DM, arg2::PetscReal, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::DMLabel, arg7::PetscInt, arg8::Ptr{PetscInt}, arg9::Ptr{Cvoid}, arg10::Ptr{Cvoid}, arg11::Vec)::PetscErrorCode
end

function DMPlexInsertBoundaryValuesEssentialField(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12)
    @ccall libpetsc.DMPlexInsertBoundaryValuesEssentialField(arg1::DM, arg2::PetscReal, arg3::Vec, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::DMLabel, arg8::PetscInt, arg9::Ptr{PetscInt}, arg10::Ptr{Cvoid}, arg11::Ptr{Cvoid}, arg12::Vec)::PetscErrorCode
end

function DMPlexInsertBoundaryValuesEssentialBdField(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12)
    @ccall libpetsc.DMPlexInsertBoundaryValuesEssentialBdField(arg1::DM, arg2::PetscReal, arg3::Vec, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::DMLabel, arg8::PetscInt, arg9::Ptr{PetscInt}, arg10::Ptr{Cvoid}, arg11::Ptr{Cvoid}, arg12::Vec)::PetscErrorCode
end

function DMPlexInsertBoundaryValuesRiemann(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13, arg14)
    @ccall libpetsc.DMPlexInsertBoundaryValuesRiemann(arg1::DM, arg2::PetscReal, arg3::Vec, arg4::Vec, arg5::Vec, arg6::PetscInt, arg7::PetscInt, arg8::Ptr{PetscInt}, arg9::DMLabel, arg10::PetscInt, arg11::Ptr{PetscInt}, arg12::Ptr{Cvoid}, arg13::Ptr{Cvoid}, arg14::Vec)::PetscErrorCode
end

function DMPlexMarkBoundaryFaces(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexMarkBoundaryFaces(arg1::DM, arg2::PetscInt, arg3::DMLabel)::PetscErrorCode
end

function DMPlexCreateSection(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.DMPlexCreateSection(arg1::DM, arg2::Ptr{DMLabel}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::Ptr{IS}, arg8::Ptr{IS}, arg9::IS, arg10::Ptr{PetscSection})::PetscErrorCode
end

function DMPlexGetSubdomainSection(arg1, arg2)
    @ccall libpetsc.DMPlexGetSubdomainSection(arg1::DM, arg2::Ptr{PetscSection})::PetscErrorCode
end

function DMPlexComputeCellGeometryAffineFEM(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMPlexComputeCellGeometryAffineFEM(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal})::PetscErrorCode
end

function DMPlexComputeCellGeometryFEM(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMPlexComputeCellGeometryFEM(arg1::DM, arg2::PetscInt, arg3::PetscQuadrature, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal}, arg7::Ptr{PetscReal})::PetscErrorCode
end

function DMPlexCoordinatesToReference(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexCoordinatesToReference(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal})::PetscErrorCode
end

function DMPlexReferenceToCoordinates(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexReferenceToCoordinates(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal})::PetscErrorCode
end

function DMPlexShearGeometry(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexShearGeometry(arg1::DM, arg2::DMDirection, arg3::Ptr{PetscReal})::PetscErrorCode
end

function DMPlexRemapGeometry(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexRemapGeometry(arg1::DM, arg2::PetscReal, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexVecGetClosure(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMPlexVecGetClosure(arg1::DM, arg2::PetscSection, arg3::Vec, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function DMPlexVecRestoreClosure(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMPlexVecRestoreClosure(arg1::DM, arg2::PetscSection, arg3::Vec, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function DMPlexVecSetClosure(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMPlexVecSetClosure(arg1::DM, arg2::PetscSection, arg3::Vec, arg4::PetscInt, arg5::Ptr{PetscScalar}, arg6::InsertMode)::PetscErrorCode
end

function DMPlexMatSetClosure(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMPlexMatSetClosure(arg1::DM, arg2::PetscSection, arg3::PetscSection, arg4::Mat, arg5::PetscInt, arg6::Ptr{PetscScalar}, arg7::InsertMode)::PetscErrorCode
end

function DMPlexMatSetClosureGeneral(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.DMPlexMatSetClosureGeneral(arg1::DM, arg2::PetscSection, arg3::PetscSection, arg4::DM, arg5::PetscSection, arg6::PetscSection, arg7::Mat, arg8::PetscInt, arg9::Ptr{PetscScalar}, arg10::InsertMode)::PetscErrorCode
end

function DMPlexGetClosureIndices(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.DMPlexGetClosureIndices(arg1::DM, arg2::PetscSection, arg3::PetscSection, arg4::PetscInt, arg5::PetscBool, arg6::Ptr{PetscInt}, arg7::Ptr{Ptr{PetscInt}}, arg8::Ptr{PetscInt}, arg9::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function DMPlexRestoreClosureIndices(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.DMPlexRestoreClosureIndices(arg1::DM, arg2::PetscSection, arg3::PetscSection, arg4::PetscInt, arg5::PetscBool, arg6::Ptr{PetscInt}, arg7::Ptr{Ptr{PetscInt}}, arg8::Ptr{PetscInt}, arg9::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function DMPlexMatSetClosureRefined(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.DMPlexMatSetClosureRefined(arg1::DM, arg2::PetscSection, arg3::PetscSection, arg4::DM, arg5::PetscSection, arg6::PetscSection, arg7::Mat, arg8::PetscInt, arg9::Ptr{PetscScalar}, arg10::InsertMode)::PetscErrorCode
end

function DMPlexMatGetClosureIndicesRefined(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.DMPlexMatGetClosureIndicesRefined(arg1::DM, arg2::PetscSection, arg3::PetscSection, arg4::DM, arg5::PetscSection, arg6::PetscSection, arg7::PetscInt, arg8::Ptr{PetscInt}, arg9::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexCreateClosureIndex(arg1, arg2)
    @ccall libpetsc.DMPlexCreateClosureIndex(arg1::DM, arg2::PetscSection)::PetscErrorCode
end

function DMPlexSetClosurePermutationTensor(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexSetClosurePermutationTensor(arg1::DM, arg2::PetscInt, arg3::PetscSection)::PetscErrorCode
end

function DMPlexConstructGhostCells(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexConstructGhostCells(arg1::DM, arg2::Ptr{Cchar}, arg3::Ptr{PetscInt}, arg4::Ptr{DM})::PetscErrorCode
end

function DMPlexConstructCohesiveCells(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexConstructCohesiveCells(arg1::DM, arg2::DMLabel, arg3::DMLabel, arg4::Ptr{DM})::PetscErrorCode
end

function DMPlexGetVTKCellHeight(arg1, arg2)
    @ccall libpetsc.DMPlexGetVTKCellHeight(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexSetVTKCellHeight(arg1, arg2)
    @ccall libpetsc.DMPlexSetVTKCellHeight(arg1::DM, arg2::PetscInt)::PetscErrorCode
end

function DMPlexVTKWriteAll(arg1, arg2)
    @ccall libpetsc.DMPlexVTKWriteAll(arg1::PetscObject, arg2::PetscViewer)::PetscErrorCode
end

function DMPlexGetGhostCellStratum(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexGetGhostCellStratum(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexGetSimplexOrBoxCells(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexGetSimplexOrBoxCells(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexGetCellFields(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.DMPlexGetCellFields(arg1::DM, arg2::IS, arg3::Vec, arg4::Vec, arg5::Vec, arg6::Ptr{Ptr{PetscScalar}}, arg7::Ptr{Ptr{PetscScalar}}, arg8::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function DMPlexRestoreCellFields(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.DMPlexRestoreCellFields(arg1::DM, arg2::IS, arg3::Vec, arg4::Vec, arg5::Vec, arg6::Ptr{Ptr{PetscScalar}}, arg7::Ptr{Ptr{PetscScalar}}, arg8::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function DMPlexGetFaceFields(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11)
    @ccall libpetsc.DMPlexGetFaceFields(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::Vec, arg5::Vec, arg6::Vec, arg7::Vec, arg8::Vec, arg9::Ptr{PetscInt}, arg10::Ptr{Ptr{PetscScalar}}, arg11::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function DMPlexRestoreFaceFields(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11)
    @ccall libpetsc.DMPlexRestoreFaceFields(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::Vec, arg5::Vec, arg6::Vec, arg7::Vec, arg8::Vec, arg9::Ptr{PetscInt}, arg10::Ptr{Ptr{PetscScalar}}, arg11::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function DMPlexGetFaceGeometry(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.DMPlexGetFaceGeometry(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::Vec, arg5::Vec, arg6::Ptr{PetscInt}, arg7::Ptr{Ptr{PetscFVFaceGeom}}, arg8::Ptr{Ptr{PetscReal}})::PetscErrorCode
end

function DMPlexRestoreFaceGeometry(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.DMPlexRestoreFaceGeometry(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::Vec, arg5::Vec, arg6::Ptr{PetscInt}, arg7::Ptr{Ptr{PetscFVFaceGeom}}, arg8::Ptr{Ptr{PetscReal}})::PetscErrorCode
end

function DMPlexGetScale(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexGetScale(arg1::DM, arg2::PetscUnit, arg3::Ptr{PetscReal})::PetscErrorCode
end

function DMPlexSetScale(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexSetScale(arg1::DM, arg2::PetscUnit, arg3::PetscReal)::PetscErrorCode
end

mutable struct JacActionCtx
    dm::DM
    u::Vec
    J::Mat
    user::Ptr{Cvoid}
    JacActionCtx() = new()
end

function DMPlexSetMaxProjectionHeight(arg1, arg2)
    @ccall libpetsc.DMPlexSetMaxProjectionHeight(arg1::DM, arg2::PetscInt)::PetscErrorCode
end

function DMPlexGetMaxProjectionHeight(arg1, arg2)
    @ccall libpetsc.DMPlexGetMaxProjectionHeight(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexGetActivePoint(arg1, arg2)
    @ccall libpetsc.DMPlexGetActivePoint(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexSetActivePoint(arg1, arg2)
    @ccall libpetsc.DMPlexSetActivePoint(arg1::DM, arg2::PetscInt)::PetscErrorCode
end

function DMPlexProjectFieldLocal(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexProjectFieldLocal(arg1::DM, arg2::Vec, arg3::Ptr{Ptr{Cvoid}}, arg4::InsertMode, arg5::Vec)::PetscErrorCode
end

function DMPlexComputeL2DiffLocal(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMPlexComputeL2DiffLocal(arg1::DM, arg2::PetscReal, arg3::Ptr{Ptr{Cvoid}}, arg4::Ptr{Ptr{Cvoid}}, arg5::Vec, arg6::Ptr{PetscReal})::PetscErrorCode
end

function DMPlexComputeL2FieldDiff(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMPlexComputeL2FieldDiff(arg1::DM, arg2::PetscReal, arg3::Ptr{Ptr{Cvoid}}, arg4::Ptr{Ptr{Cvoid}}, arg5::Vec, arg6::Ptr{PetscReal})::PetscErrorCode
end

function DMPlexComputeL2DiffVec(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMPlexComputeL2DiffVec(arg1::DM, arg2::PetscReal, arg3::Ptr{Ptr{Cvoid}}, arg4::Ptr{Ptr{Cvoid}}, arg5::Vec, arg6::Vec)::PetscErrorCode
end

function DMPlexComputeCellwiseIntegralFEM(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexComputeCellwiseIntegralFEM(arg1::DM, arg2::Vec, arg3::Vec, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexComputeIntegralFEM(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexComputeIntegralFEM(arg1::DM, arg2::Vec, arg3::Ptr{PetscScalar}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexComputeBdIntegral(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.DMPlexComputeBdIntegral(arg1::DM, arg2::Vec, arg3::DMLabel, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{Cvoid}, arg7::Ptr{PetscScalar}, arg8::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexComputeInterpolatorNested(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexComputeInterpolatorNested(arg1::DM, arg2::DM, arg3::PetscBool, arg4::Mat, arg5::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexComputeInterpolatorGeneral(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexComputeInterpolatorGeneral(arg1::DM, arg2::DM, arg3::Mat, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexComputeGradientClementInterpolant(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexComputeGradientClementInterpolant(arg1::DM, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function DMPlexComputeInjectorFEM(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexComputeInjectorFEM(arg1::DM, arg2::DM, arg3::Ptr{VecScatter}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexComputeMassMatrixNested(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexComputeMassMatrixNested(arg1::DM, arg2::DM, arg3::Mat, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexComputeMassMatrixGeneral(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexComputeMassMatrixGeneral(arg1::DM, arg2::DM, arg3::Mat, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexCreateRigidBody(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexCreateRigidBody(arg1::DM, arg2::PetscInt, arg3::Ptr{MatNullSpace})::PetscErrorCode
end

function DMPlexCreateRigidBodies(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMPlexCreateRigidBodies(arg1::DM, arg2::PetscInt, arg3::DMLabel, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{MatNullSpace})::PetscErrorCode
end

function DMPlexSetSNESLocalFEM(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexSetSNESLocalFEM(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexSNESComputeBoundaryFEM(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexSNESComputeBoundaryFEM(arg1::DM, arg2::Vec, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexSNESComputeResidualFEM(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexSNESComputeResidualFEM(arg1::DM, arg2::Vec, arg3::Vec, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexSNESComputeJacobianFEM(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexSNESComputeJacobianFEM(arg1::DM, arg2::Vec, arg3::Mat, arg4::Mat, arg5::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexComputeJacobianAction(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.DMPlexComputeJacobianAction(arg1::DM, arg2::IS, arg3::PetscReal, arg4::PetscReal, arg5::Vec, arg6::Vec, arg7::Vec, arg8::Vec, arg9::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexComputeBdResidualSingle(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.DMPlexComputeBdResidualSingle(arg1::DM, arg2::PetscReal, arg3::DMLabel, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::PetscInt, arg7::Vec, arg8::Vec, arg9::Vec)::PetscErrorCode
end

function DMPlexComputeBdJacobianSingle(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11)
    @ccall libpetsc.DMPlexComputeBdJacobianSingle(arg1::DM, arg2::PetscReal, arg3::DMLabel, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::PetscInt, arg7::Vec, arg8::Vec, arg9::PetscReal, arg10::Mat, arg11::Mat)::PetscErrorCode
end

function DMPlexTSComputeBoundary(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexTSComputeBoundary(arg1::DM, arg2::PetscReal, arg3::Vec, arg4::Vec, arg5::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexTSComputeRHSFunctionFVM(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexTSComputeRHSFunctionFVM(arg1::DM, arg2::PetscReal, arg3::Vec, arg4::Vec, arg5::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexTSComputeIFunctionFEM(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMPlexTSComputeIFunctionFEM(arg1::DM, arg2::PetscReal, arg3::Vec, arg4::Vec, arg5::Vec, arg6::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexTSComputeIJacobianFEM(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.DMPlexTSComputeIJacobianFEM(arg1::DM, arg2::PetscReal, arg3::Vec, arg4::Vec, arg5::PetscReal, arg6::Mat, arg7::Mat, arg8::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexComputeRHSFunctionFVM(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMPlexComputeRHSFunctionFVM(arg1::DM, arg2::PetscReal, arg3::Vec, arg4::Vec, arg5::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexReconstructGradientsFVM(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexReconstructGradientsFVM(arg1::DM, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function DMPlexGetAnchors(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexGetAnchors(arg1::DM, arg2::Ptr{PetscSection}, arg3::Ptr{IS})::PetscErrorCode
end

function DMPlexSetAnchors(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexSetAnchors(arg1::DM, arg2::PetscSection, arg3::IS)::PetscErrorCode
end

function DMPlexSetReferenceTree(arg1, arg2)
    @ccall libpetsc.DMPlexSetReferenceTree(arg1::DM, arg2::DM)::PetscErrorCode
end

function DMPlexGetReferenceTree(arg1, arg2)
    @ccall libpetsc.DMPlexGetReferenceTree(arg1::DM, arg2::Ptr{DM})::PetscErrorCode
end

function DMPlexReferenceTreeGetChildSymmetry(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.DMPlexReferenceTreeGetChildSymmetry(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{PetscInt}, arg8::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexCreateDefaultReferenceTree(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexCreateDefaultReferenceTree(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscBool, arg4::Ptr{DM})::PetscErrorCode
end

function DMPlexSetTree(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexSetTree(arg1::DM, arg2::PetscSection, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexGetTree(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMPlexGetTree(arg1::DM, arg2::Ptr{PetscSection}, arg3::Ptr{Ptr{PetscInt}}, arg4::Ptr{Ptr{PetscInt}}, arg5::Ptr{PetscSection}, arg6::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMPlexGetTreeParent(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexGetTreeParent(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMPlexGetTreeChildren(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexGetTreeChildren(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMPlexTreeRefineCell(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexTreeRefineCell(arg1::DM, arg2::PetscInt, arg3::Ptr{DM})::PetscErrorCode
end

function DMPlexComputeInjectorReferenceTree(arg1, arg2)
    @ccall libpetsc.DMPlexComputeInjectorReferenceTree(arg1::DM, arg2::Ptr{Mat})::PetscErrorCode
end

function DMPlexTransferVecTree(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.DMPlexTransferVecTree(arg1::DM, arg2::Vec, arg3::DM, arg4::Vec, arg5::PetscSF, arg6::PetscSF, arg7::Ptr{PetscInt}, arg8::Ptr{PetscInt}, arg9::PetscBool, arg10::PetscReal)::PetscErrorCode
end

function DMPlexMonitorThroughput(arg1, arg2)
    @ccall libpetsc.DMPlexMonitorThroughput(arg1::DM, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexCreateGlobalToNaturalSF(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexCreateGlobalToNaturalSF(arg1::DM, arg2::PetscSection, arg3::PetscSF, arg4::Ptr{PetscSF})::PetscErrorCode
end

function DMPlexSetGlobalToNaturalSF(arg1, arg2)
    @ccall libpetsc.DMPlexSetGlobalToNaturalSF(arg1::DM, arg2::PetscSF)::PetscErrorCode
end

function DMPlexGetGlobalToNaturalSF(arg1, arg2)
    @ccall libpetsc.DMPlexGetGlobalToNaturalSF(arg1::DM, arg2::Ptr{PetscSF})::PetscErrorCode
end

function DMPlexGlobalToNaturalBegin(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexGlobalToNaturalBegin(arg1::DM, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function DMPlexGlobalToNaturalEnd(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexGlobalToNaturalEnd(arg1::DM, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function DMPlexNaturalToGlobalBegin(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexNaturalToGlobalBegin(arg1::DM, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function DMPlexNaturalToGlobalEnd(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexNaturalToGlobalEnd(arg1::DM, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function DMPlexAdapt(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexAdapt(arg1::DM, arg2::Vec, arg3::Ptr{Cchar}, arg4::Ptr{DM})::PetscErrorCode
end

function DMPlexSnapToGeomModel(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexSnapToGeomModel(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscScalar}, arg4::Ptr{PetscScalar})::PetscErrorCode
end

function DMPlexGlobalToLocalBasis(arg1, arg2)
    @ccall libpetsc.DMPlexGlobalToLocalBasis(arg1::DM, arg2::Vec)::PetscErrorCode
end

function DMPlexLocalToGlobalBasis(arg1, arg2)
    @ccall libpetsc.DMPlexLocalToGlobalBasis(arg1::DM, arg2::Vec)::PetscErrorCode
end

function DMPlexCreateBasisRotation(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexCreateBasisRotation(arg1::DM, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal)::PetscErrorCode
end

function DMPlexCellRefinerCreate(arg1, arg2)
    @ccall libpetsc.DMPlexCellRefinerCreate(arg1::DM, arg2::Ptr{DMPlexCellRefiner})::PetscErrorCode
end

function DMPlexCellRefinerSetUp(arg1)
    @ccall libpetsc.DMPlexCellRefinerSetUp(arg1::DMPlexCellRefiner)::PetscErrorCode
end

function DMPlexCellRefinerDestroy(arg1)
    @ccall libpetsc.DMPlexCellRefinerDestroy(arg1::Ptr{DMPlexCellRefiner})::PetscErrorCode
end

function DMPlexCellRefinerRefine(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.DMPlexCellRefinerRefine(arg1::DMPlexCellRefiner, arg2::DMPolytopeType, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{Ptr{DMPolytopeType}}, arg7::Ptr{Ptr{PetscInt}}, arg8::Ptr{Ptr{PetscInt}}, arg9::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMPlexCellRefinerGetAffineTransforms(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMPlexCellRefinerGetAffineTransforms(arg1::DMPlexCellRefiner, arg2::DMPolytopeType, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{PetscReal}}, arg5::Ptr{Ptr{PetscReal}}, arg6::Ptr{Ptr{PetscReal}})::PetscErrorCode
end

function DMPlexCellRefinerGetAffineFaceTransforms(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMPlexCellRefinerGetAffineFaceTransforms(arg1::DMPlexCellRefiner, arg2::DMPolytopeType, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{PetscReal}}, arg5::Ptr{Ptr{PetscReal}}, arg6::Ptr{Ptr{PetscReal}}, arg7::Ptr{Ptr{PetscReal}})::PetscErrorCode
end

function DMPlexRefineUniform(arg1, arg2, arg3)
    @ccall libpetsc.DMPlexRefineUniform(arg1::DM, arg2::DMPlexCellRefiner, arg3::Ptr{DM})::PetscErrorCode
end

function DMRedundantCreate(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMRedundantCreate(arg1::MPI_Comm, arg2::PetscMPIInt, arg3::PetscInt, arg4::Ptr{DM})::PetscErrorCode
end

function DMRedundantSetSize(arg1, arg2, arg3)
    @ccall libpetsc.DMRedundantSetSize(arg1::DM, arg2::PetscMPIInt, arg3::PetscInt)::PetscErrorCode
end

function DMRedundantGetSize(arg1, arg2, arg3)
    @ccall libpetsc.DMRedundantGetSize(arg1::DM, arg2::Ptr{PetscMPIInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMShellCreate(arg1, arg2)
    @ccall libpetsc.DMShellCreate(arg1::MPI_Comm, arg2::Ptr{DM})::PetscErrorCode
end

function DMShellSetContext(arg1, arg2)
    @ccall libpetsc.DMShellSetContext(arg1::DM, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMShellGetContext(arg1, arg2)
    @ccall libpetsc.DMShellGetContext(arg1::DM, arg2::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMShellSetMatrix(arg1, arg2)
    @ccall libpetsc.DMShellSetMatrix(arg1::DM, arg2::Mat)::PetscErrorCode
end

function DMShellSetGlobalVector(arg1, arg2)
    @ccall libpetsc.DMShellSetGlobalVector(arg1::DM, arg2::Vec)::PetscErrorCode
end

function DMShellGetGlobalVector(arg1, arg2)
    @ccall libpetsc.DMShellGetGlobalVector(arg1::DM, arg2::Ptr{Vec})::PetscErrorCode
end

function DMShellSetLocalVector(arg1, arg2)
    @ccall libpetsc.DMShellSetLocalVector(arg1::DM, arg2::Vec)::PetscErrorCode
end

function DMShellSetCreateGlobalVector(arg1, arg2)
    @ccall libpetsc.DMShellSetCreateGlobalVector(arg1::DM, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMShellSetCreateLocalVector(arg1, arg2)
    @ccall libpetsc.DMShellSetCreateLocalVector(arg1::DM, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMShellSetGlobalToLocal(arg1, arg2, arg3)
    @ccall libpetsc.DMShellSetGlobalToLocal(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMShellSetGlobalToLocalVecScatter(arg1, arg2)
    @ccall libpetsc.DMShellSetGlobalToLocalVecScatter(arg1::DM, arg2::VecScatter)::PetscErrorCode
end

function DMShellSetLocalToGlobal(arg1, arg2, arg3)
    @ccall libpetsc.DMShellSetLocalToGlobal(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMShellSetLocalToGlobalVecScatter(arg1, arg2)
    @ccall libpetsc.DMShellSetLocalToGlobalVecScatter(arg1::DM, arg2::VecScatter)::PetscErrorCode
end

function DMShellSetLocalToLocal(arg1, arg2, arg3)
    @ccall libpetsc.DMShellSetLocalToLocal(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMShellSetLocalToLocalVecScatter(arg1, arg2)
    @ccall libpetsc.DMShellSetLocalToLocalVecScatter(arg1::DM, arg2::VecScatter)::PetscErrorCode
end

function DMShellSetCreateMatrix(arg1, arg2)
    @ccall libpetsc.DMShellSetCreateMatrix(arg1::DM, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMShellSetCoarsen(arg1, arg2)
    @ccall libpetsc.DMShellSetCoarsen(arg1::DM, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMShellGetCoarsen(arg1, arg2)
    @ccall libpetsc.DMShellGetCoarsen(arg1::DM, arg2::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMShellSetRefine(arg1, arg2)
    @ccall libpetsc.DMShellSetRefine(arg1::DM, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMShellGetRefine(arg1, arg2)
    @ccall libpetsc.DMShellGetRefine(arg1::DM, arg2::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMShellSetCreateInterpolation(arg1, arg2)
    @ccall libpetsc.DMShellSetCreateInterpolation(arg1::DM, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMShellGetCreateInterpolation(arg1, arg2)
    @ccall libpetsc.DMShellGetCreateInterpolation(arg1::DM, arg2::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMShellSetCreateRestriction(arg1, arg2)
    @ccall libpetsc.DMShellSetCreateRestriction(arg1::DM, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMShellGetCreateRestriction(arg1, arg2)
    @ccall libpetsc.DMShellGetCreateRestriction(arg1::DM, arg2::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMShellSetCreateInjection(arg1, arg2)
    @ccall libpetsc.DMShellSetCreateInjection(arg1::DM, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMShellGetCreateInjection(arg1, arg2)
    @ccall libpetsc.DMShellGetCreateInjection(arg1::DM, arg2::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMShellSetCreateFieldDecomposition(arg1, arg2)
    @ccall libpetsc.DMShellSetCreateFieldDecomposition(arg1::DM, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMShellSetCreateDomainDecomposition(arg1, arg2)
    @ccall libpetsc.DMShellSetCreateDomainDecomposition(arg1::DM, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMShellSetCreateDomainDecompositionScatters(arg1, arg2)
    @ccall libpetsc.DMShellSetCreateDomainDecompositionScatters(arg1::DM, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMShellSetCreateSubDM(arg1, arg2)
    @ccall libpetsc.DMShellSetCreateSubDM(arg1::DM, arg2::Ptr{Cvoid})::PetscErrorCode
end

function DMShellGetCreateSubDM(arg1, arg2)
    @ccall libpetsc.DMShellGetCreateSubDM(arg1::DM, arg2::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMGlobalToLocalBeginDefaultShell(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMGlobalToLocalBeginDefaultShell(arg1::DM, arg2::Vec, arg3::InsertMode, arg4::Vec)::PetscErrorCode
end

function DMGlobalToLocalEndDefaultShell(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMGlobalToLocalEndDefaultShell(arg1::DM, arg2::Vec, arg3::InsertMode, arg4::Vec)::PetscErrorCode
end

function DMLocalToGlobalBeginDefaultShell(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMLocalToGlobalBeginDefaultShell(arg1::DM, arg2::Vec, arg3::InsertMode, arg4::Vec)::PetscErrorCode
end

function DMLocalToGlobalEndDefaultShell(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMLocalToGlobalEndDefaultShell(arg1::DM, arg2::Vec, arg3::InsertMode, arg4::Vec)::PetscErrorCode
end

function DMLocalToLocalBeginDefaultShell(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMLocalToLocalBeginDefaultShell(arg1::DM, arg2::Vec, arg3::InsertMode, arg4::Vec)::PetscErrorCode
end

function DMLocalToLocalEndDefaultShell(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMLocalToLocalEndDefaultShell(arg1::DM, arg2::Vec, arg3::InsertMode, arg4::Vec)::PetscErrorCode
end

function DMSlicedCreate(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.DMSlicedCreate(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{PetscInt}, arg7::Ptr{PetscInt}, arg8::Ptr{DM})::PetscErrorCode
end

function DMSlicedSetPreallocation(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMSlicedSetPreallocation(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::PetscInt, arg5::Ptr{PetscInt})::PetscErrorCode
end

function DMSlicedSetBlockFills(arg1, arg2, arg3)
    @ccall libpetsc.DMSlicedSetBlockFills(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMSlicedSetGhosts(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMSlicedSetGhosts(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscInt})::PetscErrorCode
end

@enum DMSwarmType::UInt32 begin
    DMSWARM_BASIC = 0
    DMSWARM_PIC = 1
end

@enum DMSwarmMigrateType::UInt32 begin
    DMSWARM_MIGRATE_BASIC = 0
    DMSWARM_MIGRATE_DMCELLNSCATTER = 1
    DMSWARM_MIGRATE_DMCELLEXACT = 2
    DMSWARM_MIGRATE_USER = 3
end

@enum DMSwarmCollectType::UInt32 begin
    DMSWARM_COLLECT_BASIC = 0
    DMSWARM_COLLECT_DMDABOUNDINGBOX = 1
    DMSWARM_COLLECT_GENERAL = 2
    DMSWARM_COLLECT_USER = 3
end

@enum DMSwarmPICLayoutType::UInt32 begin
    DMSWARMPIC_LAYOUT_REGULAR = 0
    DMSWARMPIC_LAYOUT_GAUSS = 1
    DMSWARMPIC_LAYOUT_SUBDIVISION = 2
end

function DMSwarmCreateGlobalVectorFromField(arg1, arg2, arg3)
    @ccall libpetsc.DMSwarmCreateGlobalVectorFromField(arg1::DM, arg2::Ptr{Cchar}, arg3::Ptr{Vec})::PetscErrorCode
end

function DMSwarmDestroyGlobalVectorFromField(arg1, arg2, arg3)
    @ccall libpetsc.DMSwarmDestroyGlobalVectorFromField(arg1::DM, arg2::Ptr{Cchar}, arg3::Ptr{Vec})::PetscErrorCode
end

function DMSwarmCreateLocalVectorFromField(arg1, arg2, arg3)
    @ccall libpetsc.DMSwarmCreateLocalVectorFromField(arg1::DM, arg2::Ptr{Cchar}, arg3::Ptr{Vec})::PetscErrorCode
end

function DMSwarmDestroyLocalVectorFromField(arg1, arg2, arg3)
    @ccall libpetsc.DMSwarmDestroyLocalVectorFromField(arg1::DM, arg2::Ptr{Cchar}, arg3::Ptr{Vec})::PetscErrorCode
end

function DMSwarmInitializeFieldRegister(arg1)
    @ccall libpetsc.DMSwarmInitializeFieldRegister(arg1::DM)::PetscErrorCode
end

function DMSwarmFinalizeFieldRegister(arg1)
    @ccall libpetsc.DMSwarmFinalizeFieldRegister(arg1::DM)::PetscErrorCode
end

function DMSwarmSetLocalSizes(arg1, arg2, arg3)
    @ccall libpetsc.DMSwarmSetLocalSizes(arg1::DM, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function DMSwarmRegisterPetscDatatypeField(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMSwarmRegisterPetscDatatypeField(arg1::DM, arg2::Ptr{Cchar}, arg3::PetscInt, arg4::PetscDataType)::PetscErrorCode
end

function DMSwarmRegisterUserStructField(arg1, arg2, arg3)
    @ccall libpetsc.DMSwarmRegisterUserStructField(arg1::DM, arg2::Ptr{Cchar}, arg3::Csize_t)::PetscErrorCode
end

function DMSwarmRegisterUserDatatypeField(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMSwarmRegisterUserDatatypeField(arg1::DM, arg2::Ptr{Cchar}, arg3::Csize_t, arg4::PetscInt)::PetscErrorCode
end

function DMSwarmGetField(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMSwarmGetField(arg1::DM, arg2::Ptr{Cchar}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscDataType}, arg5::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMSwarmRestoreField(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMSwarmRestoreField(arg1::DM, arg2::Ptr{Cchar}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscDataType}, arg5::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMSwarmVectorDefineField(arg1, arg2)
    @ccall libpetsc.DMSwarmVectorDefineField(arg1::DM, arg2::Ptr{Cchar})::PetscErrorCode
end

function DMSwarmAddPoint(arg1)
    @ccall libpetsc.DMSwarmAddPoint(arg1::DM)::PetscErrorCode
end

function DMSwarmAddNPoints(arg1, arg2)
    @ccall libpetsc.DMSwarmAddNPoints(arg1::DM, arg2::PetscInt)::PetscErrorCode
end

function DMSwarmRemovePoint(arg1)
    @ccall libpetsc.DMSwarmRemovePoint(arg1::DM)::PetscErrorCode
end

function DMSwarmRemovePointAtIndex(arg1, arg2)
    @ccall libpetsc.DMSwarmRemovePointAtIndex(arg1::DM, arg2::PetscInt)::PetscErrorCode
end

function DMSwarmCopyPoint(dm, arg2, arg3)
    @ccall libpetsc.DMSwarmCopyPoint(dm::DM, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function DMSwarmGetLocalSize(arg1, arg2)
    @ccall libpetsc.DMSwarmGetLocalSize(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMSwarmGetSize(arg1, arg2)
    @ccall libpetsc.DMSwarmGetSize(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMSwarmMigrate(arg1, arg2)
    @ccall libpetsc.DMSwarmMigrate(arg1::DM, arg2::PetscBool)::PetscErrorCode
end

function DMSwarmCollectViewCreate(arg1)
    @ccall libpetsc.DMSwarmCollectViewCreate(arg1::DM)::PetscErrorCode
end

function DMSwarmCollectViewDestroy(arg1)
    @ccall libpetsc.DMSwarmCollectViewDestroy(arg1::DM)::PetscErrorCode
end

function DMSwarmSetCellDM(arg1, arg2)
    @ccall libpetsc.DMSwarmSetCellDM(arg1::DM, arg2::DM)::PetscErrorCode
end

function DMSwarmGetCellDM(arg1, arg2)
    @ccall libpetsc.DMSwarmGetCellDM(arg1::DM, arg2::Ptr{DM})::PetscErrorCode
end

function DMSwarmSetType(arg1, arg2)
    @ccall libpetsc.DMSwarmSetType(arg1::DM, arg2::DMSwarmType)::PetscErrorCode
end

function DMSwarmSetPointsUniformCoordinates(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMSwarmSetPointsUniformCoordinates(arg1::DM, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscInt}, arg5::InsertMode)::PetscErrorCode
end

function DMSwarmSetPointCoordinates(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMSwarmSetPointCoordinates(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscReal}, arg4::PetscBool, arg5::InsertMode)::PetscErrorCode
end

function DMSwarmInsertPointsUsingCellDM(arg1, arg2, arg3)
    @ccall libpetsc.DMSwarmInsertPointsUsingCellDM(arg1::DM, arg2::DMSwarmPICLayoutType, arg3::PetscInt)::PetscErrorCode
end

function DMSwarmSetPointCoordinatesCellwise(arg1, arg2, arg3)
    @ccall libpetsc.DMSwarmSetPointCoordinatesCellwise(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscReal})::PetscErrorCode
end

function DMSwarmSetPointCoordinatesRandom(arg1, arg2)
    @ccall libpetsc.DMSwarmSetPointCoordinatesRandom(arg1::DM, arg2::PetscInt)::PetscErrorCode
end

function DMSwarmViewFieldsXDMF(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMSwarmViewFieldsXDMF(arg1::DM, arg2::Ptr{Cchar}, arg3::PetscInt, arg4::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function DMSwarmViewXDMF(arg1, arg2)
    @ccall libpetsc.DMSwarmViewXDMF(arg1::DM, arg2::Ptr{Cchar})::PetscErrorCode
end

function DMSwarmSortGetAccess(arg1)
    @ccall libpetsc.DMSwarmSortGetAccess(arg1::DM)::PetscErrorCode
end

function DMSwarmSortRestoreAccess(arg1)
    @ccall libpetsc.DMSwarmSortRestoreAccess(arg1::DM)::PetscErrorCode
end

function DMSwarmSortGetPointsPerCell(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMSwarmSortGetPointsPerCell(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMSwarmSortGetNumberOfPointsPerCell(arg1, arg2, arg3)
    @ccall libpetsc.DMSwarmSortGetNumberOfPointsPerCell(arg1::DM, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMSwarmSortGetIsValid(arg1, arg2)
    @ccall libpetsc.DMSwarmSortGetIsValid(arg1::DM, arg2::Ptr{PetscBool})::PetscErrorCode
end

function DMSwarmSortGetSizes(arg1, arg2, arg3)
    @ccall libpetsc.DMSwarmSortGetSizes(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMSwarmProjectFields(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMSwarmProjectFields(arg1::DM, arg2::PetscInt, arg3::Ptr{Ptr{Cchar}}, arg4::Ptr{Ptr{Vec}}, arg5::PetscBool)::PetscErrorCode
end

function DMSwarmCreateMassMatrixSquare(arg1, arg2, arg3)
    @ccall libpetsc.DMSwarmCreateMassMatrixSquare(arg1::DM, arg2::DM, arg3::Ptr{Mat})::PetscErrorCode
end

function DMCreate_Product(arg1)
    @ccall libpetsc.DMCreate_Product(arg1::DM)::PetscErrorCode
end

function DMProductGetDM(arg1, arg2, arg3)
    @ccall libpetsc.DMProductGetDM(arg1::DM, arg2::PetscInt, arg3::Ptr{DM})::PetscErrorCode
end

function DMProductSetDimensionIndex(arg1, arg2, arg3)
    @ccall libpetsc.DMProductSetDimensionIndex(arg1::DM, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function DMProductSetDM(arg1, arg2, arg3)
    @ccall libpetsc.DMProductSetDM(arg1::DM, arg2::PetscInt, arg3::DM)::PetscErrorCode
end

@enum DMStagStencilLocation::UInt32 begin
    DMSTAG_NULL_LOCATION = 0
    DMSTAG_BACK_DOWN_LEFT = 1
    DMSTAG_BACK_DOWN = 2
    DMSTAG_BACK_DOWN_RIGHT = 3
    DMSTAG_BACK_LEFT = 4
    DMSTAG_BACK = 5
    DMSTAG_BACK_RIGHT = 6
    DMSTAG_BACK_UP_LEFT = 7
    DMSTAG_BACK_UP = 8
    DMSTAG_BACK_UP_RIGHT = 9
    DMSTAG_DOWN_LEFT = 10
    DMSTAG_DOWN = 11
    DMSTAG_DOWN_RIGHT = 12
    DMSTAG_LEFT = 13
    DMSTAG_ELEMENT = 14
    DMSTAG_RIGHT = 15
    DMSTAG_UP_LEFT = 16
    DMSTAG_UP = 17
    DMSTAG_UP_RIGHT = 18
    DMSTAG_FRONT_DOWN_LEFT = 19
    DMSTAG_FRONT_DOWN = 20
    DMSTAG_FRONT_DOWN_RIGHT = 21
    DMSTAG_FRONT_LEFT = 22
    DMSTAG_FRONT = 23
    DMSTAG_FRONT_RIGHT = 24
    DMSTAG_FRONT_UP_LEFT = 25
    DMSTAG_FRONT_UP = 26
    DMSTAG_FRONT_UP_RIGHT = 27
end

mutable struct DMStagStencil
    loc::DMStagStencilLocation
    i::PetscInt
    j::PetscInt
    k::PetscInt
    c::PetscInt
    DMStagStencil() = new()
end

@enum DMStagStencilType::UInt32 begin
    DMSTAG_STENCIL_NONE = 0
    DMSTAG_STENCIL_STAR = 1
    DMSTAG_STENCIL_BOX = 2
end

function DMCreate_Stag(arg1)
    @ccall libpetsc.DMCreate_Stag(arg1::DM)::PetscErrorCode
end

function DMStagCreate1d(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.DMStagCreate1d(arg1::MPI_Comm, arg2::DMBoundaryType, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::DMStagStencilType, arg7::PetscInt, arg8::Ptr{PetscInt}, arg9::Ptr{DM})::PetscErrorCode
end

function DMStagCreate2d(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13, arg14, arg15)
    @ccall libpetsc.DMStagCreate2d(arg1::MPI_Comm, arg2::DMBoundaryType, arg3::DMBoundaryType, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::PetscInt, arg8::PetscInt, arg9::PetscInt, arg10::PetscInt, arg11::DMStagStencilType, arg12::PetscInt, arg13::Ptr{PetscInt}, arg14::Ptr{PetscInt}, arg15::Ptr{DM})::PetscErrorCode
end

function DMStagCreate3d(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13, arg14, arg15, arg16, arg17, arg18, arg19, arg20)
    @ccall libpetsc.DMStagCreate3d(arg1::MPI_Comm, arg2::DMBoundaryType, arg3::DMBoundaryType, arg4::DMBoundaryType, arg5::PetscInt, arg6::PetscInt, arg7::PetscInt, arg8::PetscInt, arg9::PetscInt, arg10::PetscInt, arg11::PetscInt, arg12::PetscInt, arg13::PetscInt, arg14::PetscInt, arg15::DMStagStencilType, arg16::PetscInt, arg17::Ptr{PetscInt}, arg18::Ptr{PetscInt}, arg19::Ptr{PetscInt}, arg20::Ptr{DM})::PetscErrorCode
end

function DMStagCreateCompatibleDMStag(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMStagCreateCompatibleDMStag(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{DM})::PetscErrorCode
end

function DMStagGetBoundaryTypes(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMStagGetBoundaryTypes(arg1::DM, arg2::Ptr{DMBoundaryType}, arg3::Ptr{DMBoundaryType}, arg4::Ptr{DMBoundaryType})::PetscErrorCode
end

function DMStagGetCorners(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.DMStagGetCorners(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscInt}, arg7::Ptr{PetscInt}, arg8::Ptr{PetscInt}, arg9::Ptr{PetscInt}, arg10::Ptr{PetscInt})::PetscErrorCode
end

function DMStagGetDOF(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMStagGetDOF(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt})::PetscErrorCode
end

function DMStagGetEntries(arg1, arg2)
    @ccall libpetsc.DMStagGetEntries(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMStagGetEntriesPerElement(arg1, arg2)
    @ccall libpetsc.DMStagGetEntriesPerElement(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMStagGetGhostCorners(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMStagGetGhostCorners(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscInt}, arg7::Ptr{PetscInt})::PetscErrorCode
end

function DMStagGetGlobalSizes(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMStagGetGlobalSizes(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMStagGetIsFirstRank(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMStagGetIsFirstRank(arg1::DM, arg2::Ptr{PetscBool}, arg3::Ptr{PetscBool}, arg4::Ptr{PetscBool})::PetscErrorCode
end

function DMStagGetIsLastRank(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMStagGetIsLastRank(arg1::DM, arg2::Ptr{PetscBool}, arg3::Ptr{PetscBool}, arg4::Ptr{PetscBool})::PetscErrorCode
end

function DMStagGetLocalSizes(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMStagGetLocalSizes(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMStagGetLocationDOF(arg1, arg2, arg3)
    @ccall libpetsc.DMStagGetLocationDOF(arg1::DM, arg2::DMStagStencilLocation, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMStagGetLocationSlot(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMStagGetLocationSlot(arg1::DM, arg2::DMStagStencilLocation, arg3::PetscInt, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMStagGetNumRanks(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMStagGetNumRanks(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMStagGetOwnershipRanges(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMStagGetOwnershipRanges(arg1::DM, arg2::Ptr{Ptr{PetscInt}}, arg3::Ptr{Ptr{PetscInt}}, arg4::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function DMStagGetProductCoordinateArrays(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMStagGetProductCoordinateArrays(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMStagGetProductCoordinateArraysRead(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMStagGetProductCoordinateArraysRead(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMStagGetProductCoordinateLocationSlot(arg1, arg2, arg3)
    @ccall libpetsc.DMStagGetProductCoordinateLocationSlot(arg1::DM, arg2::DMStagStencilLocation, arg3::Ptr{PetscInt})::PetscErrorCode
end

function DMStagGetStencilType(arg1, arg2)
    @ccall libpetsc.DMStagGetStencilType(arg1::DM, arg2::Ptr{DMStagStencilType})::PetscErrorCode
end

function DMStagGetStencilWidth(arg1, arg2)
    @ccall libpetsc.DMStagGetStencilWidth(arg1::DM, arg2::Ptr{PetscInt})::PetscErrorCode
end

function DMStagMatGetValuesStencil(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMStagMatGetValuesStencil(arg1::DM, arg2::Mat, arg3::PetscInt, arg4::Ptr{DMStagStencil}, arg5::PetscInt, arg6::Ptr{DMStagStencil}, arg7::Ptr{PetscScalar})::PetscErrorCode
end

function DMStagMatSetValuesStencil(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.DMStagMatSetValuesStencil(arg1::DM, arg2::Mat, arg3::PetscInt, arg4::Ptr{DMStagStencil}, arg5::PetscInt, arg6::Ptr{DMStagStencil}, arg7::Ptr{PetscScalar}, arg8::InsertMode)::PetscErrorCode
end

function DMStagMigrateVec(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMStagMigrateVec(arg1::DM, arg2::Vec, arg3::DM, arg4::Vec)::PetscErrorCode
end

function DMStagPopulateLocalToGlobalInjective(arg1)
    @ccall libpetsc.DMStagPopulateLocalToGlobalInjective(arg1::DM)::PetscErrorCode
end

function DMStagRestoreProductCoordinateArrays(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMStagRestoreProductCoordinateArrays(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMStagRestoreProductCoordinateArraysRead(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMStagRestoreProductCoordinateArraysRead(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMStagSetBoundaryTypes(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMStagSetBoundaryTypes(arg1::DM, arg2::DMBoundaryType, arg3::DMBoundaryType, arg4::DMBoundaryType)::PetscErrorCode
end

function DMStagSetCoordinateDMType(arg1, arg2)
    @ccall libpetsc.DMStagSetCoordinateDMType(arg1::DM, arg2::DMType)::PetscErrorCode
end

function DMStagSetDOF(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMStagSetDOF(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt)::PetscErrorCode
end

function DMStagSetGlobalSizes(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMStagSetGlobalSizes(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt)::PetscErrorCode
end

function DMStagSetNumRanks(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMStagSetNumRanks(arg1::DM, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt)::PetscErrorCode
end

function DMStagSetOwnershipRanges(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMStagSetOwnershipRanges(arg1::DM, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function DMStagSetStencilType(arg1, arg2)
    @ccall libpetsc.DMStagSetStencilType(arg1::DM, arg2::DMStagStencilType)::PetscErrorCode
end

function DMStagSetStencilWidth(arg1, arg2)
    @ccall libpetsc.DMStagSetStencilWidth(arg1::DM, arg2::PetscInt)::PetscErrorCode
end

function DMStagSetUniformCoordinates(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMStagSetUniformCoordinates(arg1::DM, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal, arg6::PetscReal, arg7::PetscReal)::PetscErrorCode
end

function DMStagSetUniformCoordinatesExplicit(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMStagSetUniformCoordinatesExplicit(arg1::DM, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal, arg6::PetscReal, arg7::PetscReal)::PetscErrorCode
end

function DMStagSetUniformCoordinatesProduct(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.DMStagSetUniformCoordinatesProduct(arg1::DM, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal, arg6::PetscReal, arg7::PetscReal)::PetscErrorCode
end

function DMStagVecGetArray(arg1, arg2, arg3)
    @ccall libpetsc.DMStagVecGetArray(arg1::DM, arg2::Vec, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMStagVecGetArrayRead(arg1, arg2, arg3)
    @ccall libpetsc.DMStagVecGetArrayRead(arg1::DM, arg2::Vec, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMStagVecGetValuesStencil(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMStagVecGetValuesStencil(arg1::DM, arg2::Vec, arg3::PetscInt, arg4::Ptr{DMStagStencil}, arg5::Ptr{PetscScalar})::PetscErrorCode
end

function DMStagVecRestoreArray(arg1, arg2, arg3)
    @ccall libpetsc.DMStagVecRestoreArray(arg1::DM, arg2::Vec, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMStagVecRestoreArrayRead(arg1, arg2, arg3)
    @ccall libpetsc.DMStagVecRestoreArrayRead(arg1::DM, arg2::Vec, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMStagVecSetValuesStencil(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMStagVecSetValuesStencil(arg1::DM, arg2::Vec, arg3::PetscInt, arg4::Ptr{DMStagStencil}, arg5::Ptr{PetscScalar}, arg6::InsertMode)::PetscErrorCode
end

function DMStagVecSplitToDMDA(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMStagVecSplitToDMDA(arg1::DM, arg2::Vec, arg3::DMStagStencilLocation, arg4::PetscInt, arg5::Ptr{DM}, arg6::Ptr{Vec})::PetscErrorCode
end

function DMStagGet1dCoordinateArraysDOFRead(dm, ax, ay, az)
    @ccall libpetsc.DMStagGet1dCoordinateArraysDOFRead(dm::DM, ax::Ptr{Cvoid}, ay::Ptr{Cvoid}, az::Ptr{Cvoid})::PetscErrorCode
end

function DMStagGet1dCoordinateLocationSlot(dm, loc, s)
    @ccall libpetsc.DMStagGet1dCoordinateLocationSlot(dm::DM, loc::DMStagStencilLocation, s::Ptr{PetscInt})::PetscErrorCode
end

function DMStagGetGhostType(dm, s)
    @ccall libpetsc.DMStagGetGhostType(dm::DM, s::Ptr{DMStagStencilType})::PetscErrorCode
end

function DMStagRestore1dCoordinateArraysDOFRead(dm, ax, ay, az)
    @ccall libpetsc.DMStagRestore1dCoordinateArraysDOFRead(dm::DM, ax::Ptr{Cvoid}, ay::Ptr{Cvoid}, az::Ptr{Cvoid})::PetscErrorCode
end

function DMStagSetGhostType(dm, s)
    @ccall libpetsc.DMStagSetGhostType(dm::DM, s::Ptr{DMStagStencilType})::PetscErrorCode
end

function DMStagVecGetArrayDOF(dm, v, a)
    @ccall libpetsc.DMStagVecGetArrayDOF(dm::DM, v::Vec, a::Ptr{Cvoid})::PetscErrorCode
end

function DMStagVecGetArrayDOFRead(dm, v, a)
    @ccall libpetsc.DMStagVecGetArrayDOFRead(dm::DM, v::Vec, a::Ptr{Cvoid})::PetscErrorCode
end

function DMStagVecRestoreArrayDOF(dm, v, a)
    @ccall libpetsc.DMStagVecRestoreArrayDOF(dm::DM, v::Vec, a::Ptr{Cvoid})::PetscErrorCode
end

function DMStagVecRestoreArrayDOFRead(dm, v, a)
    @ccall libpetsc.DMStagVecRestoreArrayDOFRead(dm::DM, v::Vec, a::Ptr{Cvoid})::PetscErrorCode
end

function PetscWeakFormCreate(arg1, arg2)
    @ccall libpetsc.PetscWeakFormCreate(arg1::MPI_Comm, arg2::Ptr{PetscWeakForm})::PetscErrorCode
end

function PetscWeakFormDestroy(arg1)
    @ccall libpetsc.PetscWeakFormDestroy(arg1::Ptr{PetscWeakForm})::PetscErrorCode
end

function PetscWeakFormView(arg1, arg2)
    @ccall libpetsc.PetscWeakFormView(arg1::PetscWeakForm, arg2::PetscViewer)::PetscErrorCode
end

function PetscWeakFormGetNumFields(arg1, arg2)
    @ccall libpetsc.PetscWeakFormGetNumFields(arg1::PetscWeakForm, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscWeakFormSetNumFields(arg1, arg2)
    @ccall libpetsc.PetscWeakFormSetNumFields(arg1::PetscWeakForm, arg2::PetscInt)::PetscErrorCode
end

function PetscHashFormKeySort(arg1, arg2)
    @ccall libpetsc.PetscHashFormKeySort(arg1::PetscInt, arg2::Ptr{PetscHashFormKey})::PetscErrorCode
end

function PetscWeakFormGetObjective(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscWeakFormGetObjective(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{Ptr{Ptr{Cvoid}}})::PetscErrorCode
end

function PetscWeakFormAddObjective(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscWeakFormAddObjective(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{Cvoid})::PetscErrorCode
end

function PetscWeakFormSetObjective(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscWeakFormSetObjective(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscWeakFormGetIndexObjective(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscWeakFormGetIndexObjective(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscWeakFormSetIndexObjective(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscWeakFormSetIndexObjective(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{Cvoid})::PetscErrorCode
end

function PetscWeakFormGetResidual(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscWeakFormGetResidual(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{Ptr{Ptr{Cvoid}}}, arg7::Ptr{PetscInt}, arg8::Ptr{Ptr{Ptr{Cvoid}}})::PetscErrorCode
end

function PetscWeakFormAddResidual(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscWeakFormAddResidual(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{Cvoid}, arg6::Ptr{Cvoid})::PetscErrorCode
end

function PetscWeakFormSetResidual(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscWeakFormSetResidual(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{Ptr{Cvoid}}, arg7::PetscInt, arg8::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscWeakFormSetIndexResidual(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscWeakFormSetIndexResidual(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{Cvoid}, arg7::PetscInt, arg8::Ptr{Cvoid})::PetscErrorCode
end

function PetscWeakFormHasJacobian(arg1, arg2)
    @ccall libpetsc.PetscWeakFormHasJacobian(arg1::PetscWeakForm, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscWeakFormGetJacobian(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13)
    @ccall libpetsc.PetscWeakFormGetJacobian(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::Ptr{Ptr{Ptr{Cvoid}}}, arg8::Ptr{PetscInt}, arg9::Ptr{Ptr{Ptr{Cvoid}}}, arg10::Ptr{PetscInt}, arg11::Ptr{Ptr{Ptr{Cvoid}}}, arg12::Ptr{PetscInt}, arg13::Ptr{Ptr{Ptr{Cvoid}}})::PetscErrorCode
end

function PetscWeakFormAddJacobian(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.PetscWeakFormAddJacobian(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{Cvoid}, arg7::Ptr{Cvoid}, arg8::Ptr{Cvoid}, arg9::Ptr{Cvoid})::PetscErrorCode
end

function PetscWeakFormSetJacobian(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13)
    @ccall libpetsc.PetscWeakFormSetJacobian(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{Ptr{Cvoid}}, arg8::PetscInt, arg9::Ptr{Ptr{Cvoid}}, arg10::PetscInt, arg11::Ptr{Ptr{Cvoid}}, arg12::PetscInt, arg13::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscWeakFormSetIndexJacobian(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13)
    @ccall libpetsc.PetscWeakFormSetIndexJacobian(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{Cvoid}, arg8::PetscInt, arg9::Ptr{Cvoid}, arg10::PetscInt, arg11::Ptr{Cvoid}, arg12::PetscInt, arg13::Ptr{Cvoid})::PetscErrorCode
end

function PetscWeakFormHasJacobianPreconditioner(arg1, arg2)
    @ccall libpetsc.PetscWeakFormHasJacobianPreconditioner(arg1::PetscWeakForm, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscWeakFormGetJacobianPreconditioner(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13)
    @ccall libpetsc.PetscWeakFormGetJacobianPreconditioner(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::Ptr{Ptr{Ptr{Cvoid}}}, arg8::Ptr{PetscInt}, arg9::Ptr{Ptr{Ptr{Cvoid}}}, arg10::Ptr{PetscInt}, arg11::Ptr{Ptr{Ptr{Cvoid}}}, arg12::Ptr{PetscInt}, arg13::Ptr{Ptr{Ptr{Cvoid}}})::PetscErrorCode
end

function PetscWeakFormAddJacobianPreconditioner(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.PetscWeakFormAddJacobianPreconditioner(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{Cvoid}, arg7::Ptr{Cvoid}, arg8::Ptr{Cvoid}, arg9::Ptr{Cvoid})::PetscErrorCode
end

function PetscWeakFormSetJacobianPreconditioner(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13)
    @ccall libpetsc.PetscWeakFormSetJacobianPreconditioner(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{Ptr{Cvoid}}, arg8::PetscInt, arg9::Ptr{Ptr{Cvoid}}, arg10::PetscInt, arg11::Ptr{Ptr{Cvoid}}, arg12::PetscInt, arg13::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscWeakFormSetIndexJacobianPreconditioner(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13)
    @ccall libpetsc.PetscWeakFormSetIndexJacobianPreconditioner(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{Cvoid}, arg8::PetscInt, arg9::Ptr{Cvoid}, arg10::PetscInt, arg11::Ptr{Cvoid}, arg12::PetscInt, arg13::Ptr{Cvoid})::PetscErrorCode
end

function PetscWeakFormHasDynamicJacobian(arg1, arg2)
    @ccall libpetsc.PetscWeakFormHasDynamicJacobian(arg1::PetscWeakForm, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscWeakFormGetDynamicJacobian(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13)
    @ccall libpetsc.PetscWeakFormGetDynamicJacobian(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::Ptr{Ptr{Ptr{Cvoid}}}, arg8::Ptr{PetscInt}, arg9::Ptr{Ptr{Ptr{Cvoid}}}, arg10::Ptr{PetscInt}, arg11::Ptr{Ptr{Ptr{Cvoid}}}, arg12::Ptr{PetscInt}, arg13::Ptr{Ptr{Ptr{Cvoid}}})::PetscErrorCode
end

function PetscWeakFormAddDynamicJacobian(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.PetscWeakFormAddDynamicJacobian(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{Cvoid}, arg7::Ptr{Cvoid}, arg8::Ptr{Cvoid}, arg9::Ptr{Cvoid})::PetscErrorCode
end

function PetscWeakFormSetDynamicJacobian(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13)
    @ccall libpetsc.PetscWeakFormSetDynamicJacobian(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{Ptr{Cvoid}}, arg8::PetscInt, arg9::Ptr{Ptr{Cvoid}}, arg10::PetscInt, arg11::Ptr{Ptr{Cvoid}}, arg12::PetscInt, arg13::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscWeakFormSetIndexDynamicJacobian(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13)
    @ccall libpetsc.PetscWeakFormSetIndexDynamicJacobian(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{Cvoid}, arg8::PetscInt, arg9::Ptr{Cvoid}, arg10::PetscInt, arg11::Ptr{Cvoid}, arg12::PetscInt, arg13::Ptr{Cvoid})::PetscErrorCode
end

function PetscWeakFormGetBdResidual(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscWeakFormGetBdResidual(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{Ptr{Ptr{Cvoid}}}, arg7::Ptr{PetscInt}, arg8::Ptr{Ptr{Ptr{Cvoid}}})::PetscErrorCode
end

function PetscWeakFormAddBdResidual(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscWeakFormAddBdResidual(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{Cvoid}, arg6::Ptr{Cvoid})::PetscErrorCode
end

function PetscWeakFormSetBdResidual(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscWeakFormSetBdResidual(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{Ptr{Cvoid}}, arg7::PetscInt, arg8::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscWeakFormSetIndexBdResidual(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.PetscWeakFormSetIndexBdResidual(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{Cvoid}, arg7::PetscInt, arg8::Ptr{Cvoid})::PetscErrorCode
end

function PetscWeakFormHasBdJacobian(arg1, arg2)
    @ccall libpetsc.PetscWeakFormHasBdJacobian(arg1::PetscWeakForm, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscWeakFormGetBdJacobian(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13)
    @ccall libpetsc.PetscWeakFormGetBdJacobian(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::Ptr{Ptr{Ptr{Cvoid}}}, arg8::Ptr{PetscInt}, arg9::Ptr{Ptr{Ptr{Cvoid}}}, arg10::Ptr{PetscInt}, arg11::Ptr{Ptr{Ptr{Cvoid}}}, arg12::Ptr{PetscInt}, arg13::Ptr{Ptr{Ptr{Cvoid}}})::PetscErrorCode
end

function PetscWeakFormAddBdJacobian(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.PetscWeakFormAddBdJacobian(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{Cvoid}, arg7::Ptr{Cvoid}, arg8::Ptr{Cvoid}, arg9::Ptr{Cvoid})::PetscErrorCode
end

function PetscWeakFormSetBdJacobian(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13)
    @ccall libpetsc.PetscWeakFormSetBdJacobian(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{Ptr{Cvoid}}, arg8::PetscInt, arg9::Ptr{Ptr{Cvoid}}, arg10::PetscInt, arg11::Ptr{Ptr{Cvoid}}, arg12::PetscInt, arg13::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscWeakFormSetIndexBdJacobian(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13)
    @ccall libpetsc.PetscWeakFormSetIndexBdJacobian(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{Cvoid}, arg8::PetscInt, arg9::Ptr{Cvoid}, arg10::PetscInt, arg11::Ptr{Cvoid}, arg12::PetscInt, arg13::Ptr{Cvoid})::PetscErrorCode
end

function PetscWeakFormHasBdJacobianPreconditioner(arg1, arg2)
    @ccall libpetsc.PetscWeakFormHasBdJacobianPreconditioner(arg1::PetscWeakForm, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscWeakFormGetBdJacobianPreconditioner(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13)
    @ccall libpetsc.PetscWeakFormGetBdJacobianPreconditioner(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::Ptr{Ptr{Ptr{Cvoid}}}, arg8::Ptr{PetscInt}, arg9::Ptr{Ptr{Ptr{Cvoid}}}, arg10::Ptr{PetscInt}, arg11::Ptr{Ptr{Ptr{Cvoid}}}, arg12::Ptr{PetscInt}, arg13::Ptr{Ptr{Ptr{Cvoid}}})::PetscErrorCode
end

function PetscWeakFormAddBdJacobianPreconditioner(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.PetscWeakFormAddBdJacobianPreconditioner(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{Cvoid}, arg7::Ptr{Cvoid}, arg8::Ptr{Cvoid}, arg9::Ptr{Cvoid})::PetscErrorCode
end

function PetscWeakFormSetBdJacobianPreconditioner(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13)
    @ccall libpetsc.PetscWeakFormSetBdJacobianPreconditioner(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{Ptr{Cvoid}}, arg8::PetscInt, arg9::Ptr{Ptr{Cvoid}}, arg10::PetscInt, arg11::Ptr{Ptr{Cvoid}}, arg12::PetscInt, arg13::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscWeakFormSetIndexBdJacobianPreconditioner(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13)
    @ccall libpetsc.PetscWeakFormSetIndexBdJacobianPreconditioner(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{Cvoid}, arg8::PetscInt, arg9::Ptr{Cvoid}, arg10::PetscInt, arg11::Ptr{Cvoid}, arg12::PetscInt, arg13::Ptr{Cvoid})::PetscErrorCode
end

function PetscWeakFormGetRiemannSolver(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscWeakFormGetRiemannSolver(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{Ptr{Ptr{Cvoid}}})::PetscErrorCode
end

function PetscWeakFormSetRiemannSolver(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscWeakFormSetRiemannSolver(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscWeakFormSetIndexRiemannSolver(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscWeakFormSetIndexRiemannSolver(arg1::PetscWeakForm, arg2::DMLabel, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{Cvoid})::PetscErrorCode
end

function PetscDSInitializePackage()
    @ccall libpetsc.PetscDSInitializePackage()::PetscErrorCode
end

const PetscDSType = Ptr{Cchar}

@enum PetscDiscType::UInt32 begin
    PETSC_DISC_NONE = 0
    PETSC_DISC_FE = 1
    PETSC_DISC_FV = 2
end

# typedef void ( * PetscPointFunc ) ( PetscInt , PetscInt , PetscInt , const PetscInt [ ] , const PetscInt [ ] , const PetscScalar [ ] , const PetscScalar [ ] , const PetscScalar [ ] , const PetscInt [ ] , const PetscInt [ ] , const PetscScalar [ ] , const PetscScalar [ ] , const PetscScalar [ ] , PetscReal , const PetscReal [ ] , PetscInt , const PetscScalar [ ] , PetscScalar [ ] )
const PetscPointFunc = Ptr{Cvoid}

# typedef void ( * PetscPointJac ) ( PetscInt , PetscInt , PetscInt , const PetscInt [ ] , const PetscInt [ ] , const PetscScalar [ ] , const PetscScalar [ ] , const PetscScalar [ ] , const PetscInt [ ] , const PetscInt [ ] , const PetscScalar [ ] , const PetscScalar [ ] , const PetscScalar [ ] , PetscReal , PetscReal , const PetscReal [ ] , PetscInt , const PetscScalar [ ] , PetscScalar [ ] )
const PetscPointJac = Ptr{Cvoid}

# typedef void ( * PetscBdPointFunc ) ( PetscInt , PetscInt , PetscInt , const PetscInt [ ] , const PetscInt [ ] , const PetscScalar [ ] , const PetscScalar [ ] , const PetscScalar [ ] , const PetscInt [ ] , const PetscInt [ ] , const PetscScalar [ ] , const PetscScalar [ ] , const PetscScalar [ ] , PetscReal , const PetscReal [ ] , const PetscReal [ ] , PetscInt , const PetscScalar [ ] , PetscScalar [ ] )
const PetscBdPointFunc = Ptr{Cvoid}

# typedef void ( * PetscBdPointJac ) ( PetscInt , PetscInt , PetscInt , const PetscInt [ ] , const PetscInt [ ] , const PetscScalar [ ] , const PetscScalar [ ] , const PetscScalar [ ] , const PetscInt [ ] , const PetscInt [ ] , const PetscScalar [ ] , const PetscScalar [ ] , const PetscScalar [ ] , PetscReal , PetscReal , const PetscReal [ ] , const PetscReal [ ] , PetscInt , const PetscScalar [ ] , PetscScalar [ ] )
const PetscBdPointJac = Ptr{Cvoid}

# typedef void ( * PetscRiemannFunc ) ( PetscInt , PetscInt , const PetscReal [ ] , const PetscReal [ ] , const PetscScalar [ ] , const PetscScalar [ ] , PetscInt , const PetscScalar [ ] , PetscScalar [ ] , void * )
const PetscRiemannFunc = Ptr{Cvoid}

# typedef PetscErrorCode ( * PetscSimplePointFunc ) ( PetscInt , PetscReal , const PetscReal [ ] , PetscInt , PetscScalar [ ] , void * )
const PetscSimplePointFunc = Ptr{Cvoid}

function PetscDSCreate(arg1, arg2)
    @ccall libpetsc.PetscDSCreate(arg1::MPI_Comm, arg2::Ptr{PetscDS})::PetscErrorCode
end

function PetscDSDestroy(arg1)
    @ccall libpetsc.PetscDSDestroy(arg1::Ptr{PetscDS})::PetscErrorCode
end

function PetscDSSetType(arg1, arg2)
    @ccall libpetsc.PetscDSSetType(arg1::PetscDS, arg2::PetscDSType)::PetscErrorCode
end

function PetscDSGetType(arg1, arg2)
    @ccall libpetsc.PetscDSGetType(arg1::PetscDS, arg2::Ptr{PetscDSType})::PetscErrorCode
end

function PetscDSSetUp(arg1)
    @ccall libpetsc.PetscDSSetUp(arg1::PetscDS)::PetscErrorCode
end

function PetscDSSetFromOptions(arg1)
    @ccall libpetsc.PetscDSSetFromOptions(arg1::PetscDS)::PetscErrorCode
end

function PetscDSViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.PetscDSViewFromOptions(arg1::PetscDS, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function PetscDSView(arg1, arg2)
    @ccall libpetsc.PetscDSView(arg1::PetscDS, arg2::PetscViewer)::PetscErrorCode
end

function PetscDSRegister(arg1, arg2)
    @ccall libpetsc.PetscDSRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PetscDSRegisterDestroy()
    @ccall libpetsc.PetscDSRegisterDestroy()::PetscErrorCode
end

function PetscDSGetHeightSubspace(arg1, arg2, arg3)
    @ccall libpetsc.PetscDSGetHeightSubspace(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{PetscDS})::PetscErrorCode
end

function PetscDSGetSpatialDimension(arg1, arg2)
    @ccall libpetsc.PetscDSGetSpatialDimension(arg1::PetscDS, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscDSGetCoordinateDimension(arg1, arg2)
    @ccall libpetsc.PetscDSGetCoordinateDimension(arg1::PetscDS, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscDSSetCoordinateDimension(arg1, arg2)
    @ccall libpetsc.PetscDSSetCoordinateDimension(arg1::PetscDS, arg2::PetscInt)::PetscErrorCode
end

function PetscDSGetHybrid(arg1, arg2)
    @ccall libpetsc.PetscDSGetHybrid(arg1::PetscDS, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscDSSetHybrid(arg1, arg2)
    @ccall libpetsc.PetscDSSetHybrid(arg1::PetscDS, arg2::PetscBool)::PetscErrorCode
end

function PetscDSGetNumFields(arg1, arg2)
    @ccall libpetsc.PetscDSGetNumFields(arg1::PetscDS, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscDSGetTotalDimension(arg1, arg2)
    @ccall libpetsc.PetscDSGetTotalDimension(arg1::PetscDS, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscDSGetTotalComponents(arg1, arg2)
    @ccall libpetsc.PetscDSGetTotalComponents(arg1::PetscDS, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscDSGetFieldIndex(arg1, arg2, arg3)
    @ccall libpetsc.PetscDSGetFieldIndex(arg1::PetscDS, arg2::PetscObject, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscDSGetFieldSize(arg1, arg2, arg3)
    @ccall libpetsc.PetscDSGetFieldSize(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscDSGetFieldOffset(arg1, arg2, arg3)
    @ccall libpetsc.PetscDSGetFieldOffset(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscDSGetDimensions(arg1, arg2)
    @ccall libpetsc.PetscDSGetDimensions(arg1::PetscDS, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function PetscDSGetComponents(arg1, arg2)
    @ccall libpetsc.PetscDSGetComponents(arg1::PetscDS, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function PetscDSGetComponentOffset(arg1, arg2, arg3)
    @ccall libpetsc.PetscDSGetComponentOffset(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscDSGetComponentOffsets(arg1, arg2)
    @ccall libpetsc.PetscDSGetComponentOffsets(arg1::PetscDS, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function PetscDSGetComponentDerivativeOffsets(arg1, arg2)
    @ccall libpetsc.PetscDSGetComponentDerivativeOffsets(arg1::PetscDS, arg2::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function PetscDSGetWeakForm(arg1, arg2)
    @ccall libpetsc.PetscDSGetWeakForm(arg1::PetscDS, arg2::Ptr{PetscWeakForm})::PetscErrorCode
end

function PetscDSSetWeakForm(arg1, arg2)
    @ccall libpetsc.PetscDSSetWeakForm(arg1::PetscDS, arg2::PetscWeakForm)::PetscErrorCode
end

function PetscDSGetDiscretization(arg1, arg2, arg3)
    @ccall libpetsc.PetscDSGetDiscretization(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{PetscObject})::PetscErrorCode
end

function PetscDSSetDiscretization(arg1, arg2, arg3)
    @ccall libpetsc.PetscDSSetDiscretization(arg1::PetscDS, arg2::PetscInt, arg3::PetscObject)::PetscErrorCode
end

function PetscDSAddDiscretization(arg1, arg2)
    @ccall libpetsc.PetscDSAddDiscretization(arg1::PetscDS, arg2::PetscObject)::PetscErrorCode
end

function PetscDSGetQuadrature(arg1, arg2)
    @ccall libpetsc.PetscDSGetQuadrature(arg1::PetscDS, arg2::Ptr{PetscQuadrature})::PetscErrorCode
end

function PetscDSGetImplicit(arg1, arg2, arg3)
    @ccall libpetsc.PetscDSGetImplicit(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{PetscBool})::PetscErrorCode
end

function PetscDSSetImplicit(arg1, arg2, arg3)
    @ccall libpetsc.PetscDSSetImplicit(arg1::PetscDS, arg2::PetscInt, arg3::PetscBool)::PetscErrorCode
end

function PetscDSGetJetDegree(arg1, arg2, arg3)
    @ccall libpetsc.PetscDSGetJetDegree(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PetscDSSetJetDegree(arg1, arg2, arg3)
    @ccall libpetsc.PetscDSSetJetDegree(arg1::PetscDS, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function PetscDSGetConstants(arg1, arg2, arg3)
    @ccall libpetsc.PetscDSGetConstants(arg1::PetscDS, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function PetscDSSetConstants(arg1, arg2, arg3)
    @ccall libpetsc.PetscDSSetConstants(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{PetscScalar})::PetscErrorCode
end

function PetscDSGetObjective(arg1, arg2, arg3)
    @ccall libpetsc.PetscDSGetObjective(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscDSSetObjective(arg1, arg2, arg3)
    @ccall libpetsc.PetscDSSetObjective(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{Cvoid})::PetscErrorCode
end

function PetscDSGetResidual(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDSGetResidual(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{Ptr{Cvoid}}, arg4::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscDSSetResidual(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDSSetResidual(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function PetscDSHasJacobian(arg1, arg2)
    @ccall libpetsc.PetscDSHasJacobian(arg1::PetscDS, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscDSGetJacobian(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscDSGetJacobian(arg1::PetscDS, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Ptr{Cvoid}}, arg5::Ptr{Ptr{Cvoid}}, arg6::Ptr{Ptr{Cvoid}}, arg7::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscDSSetJacobian(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscDSSetJacobian(arg1::PetscDS, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid}, arg6::Ptr{Cvoid}, arg7::Ptr{Cvoid})::PetscErrorCode
end

function PetscDSUseJacobianPreconditioner(arg1, arg2)
    @ccall libpetsc.PetscDSUseJacobianPreconditioner(arg1::PetscDS, arg2::PetscBool)::PetscErrorCode
end

function PetscDSHasJacobianPreconditioner(arg1, arg2)
    @ccall libpetsc.PetscDSHasJacobianPreconditioner(arg1::PetscDS, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscDSGetJacobianPreconditioner(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscDSGetJacobianPreconditioner(arg1::PetscDS, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Ptr{Cvoid}}, arg5::Ptr{Ptr{Cvoid}}, arg6::Ptr{Ptr{Cvoid}}, arg7::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscDSSetJacobianPreconditioner(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscDSSetJacobianPreconditioner(arg1::PetscDS, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid}, arg6::Ptr{Cvoid}, arg7::Ptr{Cvoid})::PetscErrorCode
end

function PetscDSHasDynamicJacobian(arg1, arg2)
    @ccall libpetsc.PetscDSHasDynamicJacobian(arg1::PetscDS, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscDSGetDynamicJacobian(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscDSGetDynamicJacobian(arg1::PetscDS, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Ptr{Cvoid}}, arg5::Ptr{Ptr{Cvoid}}, arg6::Ptr{Ptr{Cvoid}}, arg7::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscDSSetDynamicJacobian(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscDSSetDynamicJacobian(arg1::PetscDS, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid}, arg6::Ptr{Cvoid}, arg7::Ptr{Cvoid})::PetscErrorCode
end

function PetscDSGetRiemannSolver(arg1, arg2, arg3)
    @ccall libpetsc.PetscDSGetRiemannSolver(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscDSSetRiemannSolver(arg1, arg2, arg3)
    @ccall libpetsc.PetscDSSetRiemannSolver(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{Cvoid})::PetscErrorCode
end

function PetscDSGetUpdate(arg1, arg2, arg3)
    @ccall libpetsc.PetscDSGetUpdate(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscDSSetUpdate(arg1, arg2, arg3)
    @ccall libpetsc.PetscDSSetUpdate(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{Cvoid})::PetscErrorCode
end

function PetscDSGetContext(arg1, arg2, arg3)
    @ccall libpetsc.PetscDSGetContext(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscDSSetContext(arg1, arg2, arg3)
    @ccall libpetsc.PetscDSSetContext(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{Cvoid})::PetscErrorCode
end

function PetscDSGetBdResidual(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDSGetBdResidual(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{Ptr{Cvoid}}, arg4::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscDSSetBdResidual(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDSSetBdResidual(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function PetscDSHasBdJacobian(arg1, arg2)
    @ccall libpetsc.PetscDSHasBdJacobian(arg1::PetscDS, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscDSGetBdJacobian(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscDSGetBdJacobian(arg1::PetscDS, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Ptr{Cvoid}}, arg5::Ptr{Ptr{Cvoid}}, arg6::Ptr{Ptr{Cvoid}}, arg7::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscDSSetBdJacobian(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscDSSetBdJacobian(arg1::PetscDS, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid}, arg6::Ptr{Cvoid}, arg7::Ptr{Cvoid})::PetscErrorCode
end

function PetscDSHasBdJacobianPreconditioner(arg1, arg2)
    @ccall libpetsc.PetscDSHasBdJacobianPreconditioner(arg1::PetscDS, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PetscDSGetBdJacobianPreconditioner(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscDSGetBdJacobianPreconditioner(arg1::PetscDS, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Ptr{Cvoid}}, arg5::Ptr{Ptr{Cvoid}}, arg6::Ptr{Ptr{Cvoid}}, arg7::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscDSSetBdJacobianPreconditioner(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscDSSetBdJacobianPreconditioner(arg1::PetscDS, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid}, arg6::Ptr{Cvoid}, arg7::Ptr{Cvoid})::PetscErrorCode
end

function PetscDSGetExactSolution(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDSGetExactSolution(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{Ptr{Cvoid}}, arg4::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscDSSetExactSolution(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDSSetExactSolution(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function PetscDSGetExactSolutionTimeDerivative(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDSGetExactSolutionTimeDerivative(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{Ptr{Cvoid}}, arg4::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscDSSetExactSolutionTimeDerivative(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDSSetExactSolutionTimeDerivative(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function PetscDSGetTabulation(arg1, arg2)
    @ccall libpetsc.PetscDSGetTabulation(arg1::PetscDS, arg2::Ptr{Ptr{PetscTabulation}})::PetscErrorCode
end

function PetscDSGetFaceTabulation(arg1, arg2)
    @ccall libpetsc.PetscDSGetFaceTabulation(arg1::PetscDS, arg2::Ptr{Ptr{PetscTabulation}})::PetscErrorCode
end

function PetscDSGetEvaluationArrays(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDSGetEvaluationArrays(arg1::PetscDS, arg2::Ptr{Ptr{PetscScalar}}, arg3::Ptr{Ptr{PetscScalar}}, arg4::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function PetscDSGetWeakFormArrays(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.PetscDSGetWeakFormArrays(arg1::PetscDS, arg2::Ptr{Ptr{PetscScalar}}, arg3::Ptr{Ptr{PetscScalar}}, arg4::Ptr{Ptr{PetscScalar}}, arg5::Ptr{Ptr{PetscScalar}}, arg6::Ptr{Ptr{PetscScalar}}, arg7::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function PetscDSGetWorkspace(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PetscDSGetWorkspace(arg1::PetscDS, arg2::Ptr{Ptr{PetscReal}}, arg3::Ptr{Ptr{PetscScalar}}, arg4::Ptr{Ptr{PetscScalar}}, arg5::Ptr{Ptr{PetscScalar}}, arg6::Ptr{Ptr{PetscScalar}})::PetscErrorCode
end

function PetscDSCopyConstants(arg1, arg2)
    @ccall libpetsc.PetscDSCopyConstants(arg1::PetscDS, arg2::PetscDS)::PetscErrorCode
end

function PetscDSCopyEquations(arg1, arg2)
    @ccall libpetsc.PetscDSCopyEquations(arg1::PetscDS, arg2::PetscDS)::PetscErrorCode
end

function PetscDSSelectDiscretizations(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDSSelectDiscretizations(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::PetscDS)::PetscErrorCode
end

function PetscDSSelectEquations(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDSSelectEquations(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::PetscDS)::PetscErrorCode
end

function PetscDSAddBoundary(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12)
    @ccall libpetsc.PetscDSAddBoundary(arg1::PetscDS, arg2::DMBoundaryConditionType, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{PetscInt}, arg8::Ptr{Cvoid}, arg9::Ptr{Cvoid}, arg10::PetscInt, arg11::Ptr{PetscInt}, arg12::Ptr{Cvoid})::PetscErrorCode
end

function PetscDSUpdateBoundary(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13)
    @ccall libpetsc.PetscDSUpdateBoundary(arg1::PetscDS, arg2::PetscInt, arg3::DMBoundaryConditionType, arg4::Ptr{Cchar}, arg5::Ptr{Cchar}, arg6::PetscInt, arg7::PetscInt, arg8::Ptr{PetscInt}, arg9::Ptr{Cvoid}, arg10::Ptr{Cvoid}, arg11::PetscInt, arg12::Ptr{PetscInt}, arg13::Ptr{Cvoid})::PetscErrorCode
end

function PetscDSGetNumBoundary(arg1, arg2)
    @ccall libpetsc.PetscDSGetNumBoundary(arg1::PetscDS, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PetscDSGetBoundary(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13)
    @ccall libpetsc.PetscDSGetBoundary(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{DMBoundaryConditionType}, arg4::Ptr{Ptr{Cchar}}, arg5::Ptr{Ptr{Cchar}}, arg6::Ptr{PetscInt}, arg7::Ptr{PetscInt}, arg8::Ptr{Ptr{PetscInt}}, arg9::Ptr{Ptr{Cvoid}}, arg10::Ptr{Ptr{Cvoid}}, arg11::Ptr{PetscInt}, arg12::Ptr{Ptr{PetscInt}}, arg13::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PetscDSCopyBoundary(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscDSCopyBoundary(arg1::PetscDS, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::PetscDS)::PetscErrorCode
end

function CharacteristicInitializePackage()
    @ccall libpetsc.CharacteristicInitializePackage()::PetscErrorCode
end

mutable struct _p_Characteristic end

const Characteristic = Ptr{_p_Characteristic}

const CharacteristicType = Ptr{Cchar}

function CharacteristicCreate(arg1, arg2)
    @ccall libpetsc.CharacteristicCreate(arg1::MPI_Comm, arg2::Ptr{Characteristic})::PetscErrorCode
end

function CharacteristicSetType(arg1, arg2)
    @ccall libpetsc.CharacteristicSetType(arg1::Characteristic, arg2::CharacteristicType)::PetscErrorCode
end

function CharacteristicSetUp(arg1)
    @ccall libpetsc.CharacteristicSetUp(arg1::Characteristic)::PetscErrorCode
end

function CharacteristicSetVelocityInterpolation(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.CharacteristicSetVelocityInterpolation(arg1::Characteristic, arg2::DM, arg3::Vec, arg4::Vec, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::Ptr{Cvoid}, arg8::Ptr{Cvoid})::PetscErrorCode
end

function CharacteristicSetVelocityInterpolationLocal(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.CharacteristicSetVelocityInterpolationLocal(arg1::Characteristic, arg2::DM, arg3::Vec, arg4::Vec, arg5::PetscInt, arg6::Ptr{PetscInt}, arg7::Ptr{Cvoid}, arg8::Ptr{Cvoid})::PetscErrorCode
end

function CharacteristicSetFieldInterpolation(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.CharacteristicSetFieldInterpolation(arg1::Characteristic, arg2::DM, arg3::Vec, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{Cvoid}, arg7::Ptr{Cvoid})::PetscErrorCode
end

function CharacteristicSetFieldInterpolationLocal(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.CharacteristicSetFieldInterpolationLocal(arg1::Characteristic, arg2::DM, arg3::Vec, arg4::PetscInt, arg5::Ptr{PetscInt}, arg6::Ptr{Cvoid}, arg7::Ptr{Cvoid})::PetscErrorCode
end

function CharacteristicSolve(arg1, arg2, arg3)
    @ccall libpetsc.CharacteristicSolve(arg1::Characteristic, arg2::PetscReal, arg3::Vec)::PetscErrorCode
end

function CharacteristicDestroy(arg1)
    @ccall libpetsc.CharacteristicDestroy(arg1::Ptr{Characteristic})::PetscErrorCode
end

function CharacteristicRegister(arg1, arg2)
    @ccall libpetsc.CharacteristicRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

mutable struct _p_PC end

const PC = Ptr{_p_PC}

const PCType = Ptr{Cchar}

@enum PCSide::Int32 begin
    PC_SIDE_DEFAULT = -1
    PC_LEFT = 0
    PC_RIGHT = 1
    PC_SYMMETRIC = 2
end

@enum PCRichardsonConvergedReason::Int32 begin
    PCRICHARDSON_CONVERGED_RTOL = 2
    PCRICHARDSON_CONVERGED_ATOL = 3
    PCRICHARDSON_CONVERGED_ITS = 4
    PCRICHARDSON_DIVERGED_DTOL = -4
end

@enum PCJacobiType::UInt32 begin
    PC_JACOBI_DIAGONAL = 0
    PC_JACOBI_ROWMAX = 1
    PC_JACOBI_ROWSUM = 2
end

@enum PCASMType::UInt32 begin
    PC_ASM_BASIC = 3
    PC_ASM_RESTRICT = 1
    PC_ASM_INTERPOLATE = 2
    PC_ASM_NONE = 0
end

@enum PCGASMType::UInt32 begin
    PC_GASM_BASIC = 3
    PC_GASM_RESTRICT = 1
    PC_GASM_INTERPOLATE = 2
    PC_GASM_NONE = 0
end

@enum PCCompositeType::UInt32 begin
    PC_COMPOSITE_ADDITIVE = 0
    PC_COMPOSITE_MULTIPLICATIVE = 1
    PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE = 2
    PC_COMPOSITE_SPECIAL = 3
    PC_COMPOSITE_SCHUR = 4
    PC_COMPOSITE_GKB = 5
end

@enum PCFieldSplitSchurPreType::UInt32 begin
    PC_FIELDSPLIT_SCHUR_PRE_SELF = 0
    PC_FIELDSPLIT_SCHUR_PRE_SELFP = 1
    PC_FIELDSPLIT_SCHUR_PRE_A11 = 2
    PC_FIELDSPLIT_SCHUR_PRE_USER = 3
    PC_FIELDSPLIT_SCHUR_PRE_FULL = 4
end

@enum PCFieldSplitSchurFactType::UInt32 begin
    PC_FIELDSPLIT_SCHUR_FACT_DIAG = 0
    PC_FIELDSPLIT_SCHUR_FACT_LOWER = 1
    PC_FIELDSPLIT_SCHUR_FACT_UPPER = 2
    PC_FIELDSPLIT_SCHUR_FACT_FULL = 3
end

@enum PCPARMSGlobalType::UInt32 begin
    PC_PARMS_GLOBAL_RAS = 0
    PC_PARMS_GLOBAL_SCHUR = 1
    PC_PARMS_GLOBAL_BJ = 2
end

@enum PCPARMSLocalType::UInt32 begin
    PC_PARMS_LOCAL_ILU0 = 0
    PC_PARMS_LOCAL_ILUK = 1
    PC_PARMS_LOCAL_ILUT = 2
    PC_PARMS_LOCAL_ARMS = 3
end

const PCGAMGType = Ptr{Cchar}

const PCGAMGClassicalType = Ptr{Cchar}

@enum PCMGType::UInt32 begin
    PC_MG_MULTIPLICATIVE = 0
    PC_MG_ADDITIVE = 1
    PC_MG_FULL = 2
    PC_MG_KASKADE = 3
end

@enum PCMGCycleType::UInt32 begin
    PC_MG_CYCLE_V = 1
    PC_MG_CYCLE_W = 2
end

@enum PCMGGalerkinType::UInt32 begin
    PC_MG_GALERKIN_BOTH = 0
    PC_MG_GALERKIN_PMAT = 1
    PC_MG_GALERKIN_MAT = 2
    PC_MG_GALERKIN_NONE = 3
    PC_MG_GALERKIN_EXTERNAL = 4
end

@enum PCExoticType::UInt32 begin
    PC_EXOTIC_FACE = 0
    PC_EXOTIC_WIREBASKET = 1
end

@enum PCBDDCInterfaceExtType::UInt32 begin
    PC_BDDC_INTERFACE_EXT_DIRICHLET = 0
    PC_BDDC_INTERFACE_EXT_LUMP = 1
end

@enum PCMGCoarseSpaceType::UInt32 begin
    PCMG_POLYNOMIAL = 0
    PCMG_HARMONIC = 1
    PCMG_EIGENVECTOR = 2
    PCMG_GENERALIZED_EIGENVECTOR = 3
end

@enum PCPatchConstructType::UInt32 begin
    PC_PATCH_STAR = 0
    PC_PATCH_VANKA = 1
    PC_PATCH_PARDECOMP = 2
    PC_PATCH_USER = 3
    PC_PATCH_PYTHON = 4
end

@enum PCDeflationSpaceType::UInt32 begin
    PC_DEFLATION_SPACE_HAAR = 0
    PC_DEFLATION_SPACE_DB2 = 1
    PC_DEFLATION_SPACE_DB4 = 2
    PC_DEFLATION_SPACE_DB8 = 3
    PC_DEFLATION_SPACE_DB16 = 4
    PC_DEFLATION_SPACE_BIORTH22 = 5
    PC_DEFLATION_SPACE_MEYER = 6
    PC_DEFLATION_SPACE_AGGREGATION = 7
    PC_DEFLATION_SPACE_USER = 8
end

@enum PCHPDDMCoarseCorrectionType::UInt32 begin
    PC_HPDDM_COARSE_CORRECTION_DEFLATED = 0
    PC_HPDDM_COARSE_CORRECTION_ADDITIVE = 1
    PC_HPDDM_COARSE_CORRECTION_BALANCED = 2
end

@enum PCFailedReason::Int32 begin
    PC_SETUP_ERROR = -1
    PC_NOERROR = 0
    PC_FACTOR_STRUCT_ZEROPIVOT = 1
    PC_FACTOR_NUMERIC_ZEROPIVOT = 2
    PC_FACTOR_OUTMEMORY = 3
    PC_FACTOR_OTHER = 4
    PC_SUBPC_ERROR = 5
end

@enum PCGAMGLayoutType::UInt32 begin
    PCGAMG_LAYOUT_COMPACT = 0
    PCGAMG_LAYOUT_SPREAD = 1
end

function PCInitializePackage()
    @ccall libpetsc.PCInitializePackage()::PetscErrorCode
end

function PCCreate(arg1, arg2)
    @ccall libpetsc.PCCreate(arg1::MPI_Comm, arg2::Ptr{PC})::PetscErrorCode
end

function PCSetType(arg1, arg2)
    @ccall libpetsc.PCSetType(arg1::PC, arg2::PCType)::PetscErrorCode
end

function PCGetType(arg1, arg2)
    @ccall libpetsc.PCGetType(arg1::PC, arg2::Ptr{PCType})::PetscErrorCode
end

function PCSetUp(arg1)
    @ccall libpetsc.PCSetUp(arg1::PC)::PetscErrorCode
end

function PCSetFailedReason(arg1, arg2)
    @ccall libpetsc.PCSetFailedReason(arg1::PC, arg2::PCFailedReason)::PetscErrorCode
end

function PCGetFailedReason(arg1, arg2)
    @ccall libpetsc.PCGetFailedReason(arg1::PC, arg2::Ptr{PCFailedReason})::PetscErrorCode
end

function PCGetSetUpFailedReason(pc, reason)
    @ccall libpetsc.PCGetSetUpFailedReason(pc::PC, reason::Ptr{PCFailedReason})::PetscErrorCode
end

function PCGetFailedReasonRank(arg1, arg2)
    @ccall libpetsc.PCGetFailedReasonRank(arg1::PC, arg2::Ptr{PCFailedReason})::PetscErrorCode
end

function PCSetUpOnBlocks(arg1)
    @ccall libpetsc.PCSetUpOnBlocks(arg1::PC)::PetscErrorCode
end

function PCApply(arg1, arg2, arg3)
    @ccall libpetsc.PCApply(arg1::PC, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function PCMatApply(arg1, arg2, arg3)
    @ccall libpetsc.PCMatApply(arg1::PC, arg2::Mat, arg3::Mat)::PetscErrorCode
end

function PCApplySymmetricLeft(arg1, arg2, arg3)
    @ccall libpetsc.PCApplySymmetricLeft(arg1::PC, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function PCApplySymmetricRight(arg1, arg2, arg3)
    @ccall libpetsc.PCApplySymmetricRight(arg1::PC, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function PCApplyBAorAB(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PCApplyBAorAB(arg1::PC, arg2::PCSide, arg3::Vec, arg4::Vec, arg5::Vec)::PetscErrorCode
end

function PCApplyTranspose(arg1, arg2, arg3)
    @ccall libpetsc.PCApplyTranspose(arg1::PC, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function PCApplyTransposeExists(arg1, arg2)
    @ccall libpetsc.PCApplyTransposeExists(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCApplyBAorABTranspose(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PCApplyBAorABTranspose(arg1::PC, arg2::PCSide, arg3::Vec, arg4::Vec, arg5::Vec)::PetscErrorCode
end

function PCSetReusePreconditioner(arg1, arg2)
    @ccall libpetsc.PCSetReusePreconditioner(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCGetReusePreconditioner(arg1, arg2)
    @ccall libpetsc.PCGetReusePreconditioner(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCSetErrorIfFailure(arg1, arg2)
    @ccall libpetsc.PCSetErrorIfFailure(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCApplyRichardson(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11)
    @ccall libpetsc.PCApplyRichardson(arg1::PC, arg2::Vec, arg3::Vec, arg4::Vec, arg5::PetscReal, arg6::PetscReal, arg7::PetscReal, arg8::PetscInt, arg9::PetscBool, arg10::Ptr{PetscInt}, arg11::Ptr{PCRichardsonConvergedReason})::PetscErrorCode
end

function PCApplyRichardsonExists(arg1, arg2)
    @ccall libpetsc.PCApplyRichardsonExists(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCSetUseAmat(arg1, arg2)
    @ccall libpetsc.PCSetUseAmat(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCGetUseAmat(arg1, arg2)
    @ccall libpetsc.PCGetUseAmat(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCRegister(arg1, arg2)
    @ccall libpetsc.PCRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PCReset(arg1)
    @ccall libpetsc.PCReset(arg1::PC)::PetscErrorCode
end

function PCDestroy(arg1)
    @ccall libpetsc.PCDestroy(arg1::Ptr{PC})::PetscErrorCode
end

function PCSetFromOptions(arg1)
    @ccall libpetsc.PCSetFromOptions(arg1::PC)::PetscErrorCode
end

function PCFactorGetMatrix(arg1, arg2)
    @ccall libpetsc.PCFactorGetMatrix(arg1::PC, arg2::Ptr{Mat})::PetscErrorCode
end

function PCSetModifySubMatrices(arg1, arg2, arg3)
    @ccall libpetsc.PCSetModifySubMatrices(arg1::PC, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function PCModifySubMatrices(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PCModifySubMatrices(arg1::PC, arg2::PetscInt, arg3::Ptr{IS}, arg4::Ptr{IS}, arg5::Ptr{Mat}, arg6::Ptr{Cvoid})::PetscErrorCode
end

function PCSetOperators(arg1, arg2, arg3)
    @ccall libpetsc.PCSetOperators(arg1::PC, arg2::Mat, arg3::Mat)::PetscErrorCode
end

function PCGetOperators(arg1, arg2, arg3)
    @ccall libpetsc.PCGetOperators(arg1::PC, arg2::Ptr{Mat}, arg3::Ptr{Mat})::PetscErrorCode
end

function PCGetOperatorsSet(arg1, arg2, arg3)
    @ccall libpetsc.PCGetOperatorsSet(arg1::PC, arg2::Ptr{PetscBool}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function PCView(arg1, arg2)
    @ccall libpetsc.PCView(arg1::PC, arg2::PetscViewer)::PetscErrorCode
end

function PCLoad(arg1, arg2)
    @ccall libpetsc.PCLoad(arg1::PC, arg2::PetscViewer)::PetscErrorCode
end

function PCViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.PCViewFromOptions(arg1::PC, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function PCSetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.PCSetOptionsPrefix(arg1::PC, arg2::Ptr{Cchar})::PetscErrorCode
end

function PCAppendOptionsPrefix(arg1, arg2)
    @ccall libpetsc.PCAppendOptionsPrefix(arg1::PC, arg2::Ptr{Cchar})::PetscErrorCode
end

function PCGetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.PCGetOptionsPrefix(arg1::PC, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PCComputeOperator(arg1, arg2, arg3)
    @ccall libpetsc.PCComputeOperator(arg1::PC, arg2::MatType, arg3::Ptr{Mat})::PetscErrorCode
end

function PCComputeExplicitOperator(A, B)
    @ccall libpetsc.PCComputeExplicitOperator(A::PC, B::Ptr{Mat})::PetscErrorCode
end

function PCGetDiagonalScale(arg1, arg2)
    @ccall libpetsc.PCGetDiagonalScale(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCDiagonalScaleLeft(arg1, arg2, arg3)
    @ccall libpetsc.PCDiagonalScaleLeft(arg1::PC, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function PCDiagonalScaleRight(arg1, arg2, arg3)
    @ccall libpetsc.PCDiagonalScaleRight(arg1::PC, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function PCSetDiagonalScale(arg1, arg2)
    @ccall libpetsc.PCSetDiagonalScale(arg1::PC, arg2::Vec)::PetscErrorCode
end

function PCSetDM(arg1, arg2)
    @ccall libpetsc.PCSetDM(arg1::PC, arg2::DM)::PetscErrorCode
end

function PCGetDM(arg1, arg2)
    @ccall libpetsc.PCGetDM(arg1::PC, arg2::Ptr{DM})::PetscErrorCode
end

function PCGetInterpolations(arg1, arg2, arg3)
    @ccall libpetsc.PCGetInterpolations(arg1::PC, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{Mat}})::PetscErrorCode
end

function PCGetCoarseOperators(pc, arg2, arg3)
    @ccall libpetsc.PCGetCoarseOperators(pc::PC, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{Mat}})::PetscErrorCode
end

function PCSetCoordinates(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PCSetCoordinates(arg1::PC, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscReal})::PetscErrorCode
end

function PCSetApplicationContext(arg1, arg2)
    @ccall libpetsc.PCSetApplicationContext(arg1::PC, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PCGetApplicationContext(arg1, arg2)
    @ccall libpetsc.PCGetApplicationContext(arg1::PC, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PCJacobiSetType(arg1, arg2)
    @ccall libpetsc.PCJacobiSetType(arg1::PC, arg2::PCJacobiType)::PetscErrorCode
end

function PCJacobiGetType(arg1, arg2)
    @ccall libpetsc.PCJacobiGetType(arg1::PC, arg2::Ptr{PCJacobiType})::PetscErrorCode
end

function PCJacobiSetUseAbs(arg1, arg2)
    @ccall libpetsc.PCJacobiSetUseAbs(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCJacobiGetUseAbs(arg1, arg2)
    @ccall libpetsc.PCJacobiGetUseAbs(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCSORSetSymmetric(arg1, arg2)
    @ccall libpetsc.PCSORSetSymmetric(arg1::PC, arg2::MatSORType)::PetscErrorCode
end

function PCSORGetSymmetric(arg1, arg2)
    @ccall libpetsc.PCSORGetSymmetric(arg1::PC, arg2::Ptr{MatSORType})::PetscErrorCode
end

function PCSORSetOmega(arg1, arg2)
    @ccall libpetsc.PCSORSetOmega(arg1::PC, arg2::PetscReal)::PetscErrorCode
end

function PCSORGetOmega(arg1, arg2)
    @ccall libpetsc.PCSORGetOmega(arg1::PC, arg2::Ptr{PetscReal})::PetscErrorCode
end

function PCSORSetIterations(arg1, arg2, arg3)
    @ccall libpetsc.PCSORSetIterations(arg1::PC, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function PCSORGetIterations(arg1, arg2, arg3)
    @ccall libpetsc.PCSORGetIterations(arg1::PC, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PCEisenstatSetOmega(arg1, arg2)
    @ccall libpetsc.PCEisenstatSetOmega(arg1::PC, arg2::PetscReal)::PetscErrorCode
end

function PCEisenstatGetOmega(arg1, arg2)
    @ccall libpetsc.PCEisenstatGetOmega(arg1::PC, arg2::Ptr{PetscReal})::PetscErrorCode
end

function PCEisenstatSetNoDiagonalScaling(arg1, arg2)
    @ccall libpetsc.PCEisenstatSetNoDiagonalScaling(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCEisenstatGetNoDiagonalScaling(arg1, arg2)
    @ccall libpetsc.PCEisenstatGetNoDiagonalScaling(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCBJacobiSetTotalBlocks(arg1, arg2, arg3)
    @ccall libpetsc.PCBJacobiSetTotalBlocks(arg1::PC, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PCBJacobiGetTotalBlocks(arg1, arg2, arg3)
    @ccall libpetsc.PCBJacobiGetTotalBlocks(arg1::PC, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function PCBJacobiSetLocalBlocks(arg1, arg2, arg3)
    @ccall libpetsc.PCBJacobiSetLocalBlocks(arg1::PC, arg2::PetscInt, arg3::Ptr{PetscInt})::PetscErrorCode
end

function PCBJacobiGetLocalBlocks(arg1, arg2, arg3)
    @ccall libpetsc.PCBJacobiGetLocalBlocks(arg1::PC, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{PetscInt}})::PetscErrorCode
end

function PCShellSetApply(arg1, arg2)
    @ccall libpetsc.PCShellSetApply(arg1::PC, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PCShellSetMatApply(arg1, arg2)
    @ccall libpetsc.PCShellSetMatApply(arg1::PC, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PCShellSetApplySymmetricLeft(arg1, arg2)
    @ccall libpetsc.PCShellSetApplySymmetricLeft(arg1::PC, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PCShellSetApplySymmetricRight(arg1, arg2)
    @ccall libpetsc.PCShellSetApplySymmetricRight(arg1::PC, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PCShellSetApplyBA(arg1, arg2)
    @ccall libpetsc.PCShellSetApplyBA(arg1::PC, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PCShellSetApplyTranspose(arg1, arg2)
    @ccall libpetsc.PCShellSetApplyTranspose(arg1::PC, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PCShellSetSetUp(arg1, arg2)
    @ccall libpetsc.PCShellSetSetUp(arg1::PC, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PCShellSetApplyRichardson(arg1, arg2)
    @ccall libpetsc.PCShellSetApplyRichardson(arg1::PC, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PCShellSetView(arg1, arg2)
    @ccall libpetsc.PCShellSetView(arg1::PC, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PCShellSetDestroy(arg1, arg2)
    @ccall libpetsc.PCShellSetDestroy(arg1::PC, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PCShellSetContext(arg1, arg2)
    @ccall libpetsc.PCShellSetContext(arg1::PC, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PCShellGetContext(arg1, arg2)
    @ccall libpetsc.PCShellGetContext(arg1::PC, arg2::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PCShellSetName(arg1, arg2)
    @ccall libpetsc.PCShellSetName(arg1::PC, arg2::Ptr{Cchar})::PetscErrorCode
end

function PCShellGetName(arg1, arg2)
    @ccall libpetsc.PCShellGetName(arg1::PC, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PCFactorSetZeroPivot(arg1, arg2)
    @ccall libpetsc.PCFactorSetZeroPivot(arg1::PC, arg2::PetscReal)::PetscErrorCode
end

function PCFactorSetShiftType(arg1, arg2)
    @ccall libpetsc.PCFactorSetShiftType(arg1::PC, arg2::MatFactorShiftType)::PetscErrorCode
end

function PCFactorSetShiftAmount(arg1, arg2)
    @ccall libpetsc.PCFactorSetShiftAmount(arg1::PC, arg2::PetscReal)::PetscErrorCode
end

function PCFactorSetMatSolverType(arg1, arg2)
    @ccall libpetsc.PCFactorSetMatSolverType(arg1::PC, arg2::MatSolverType)::PetscErrorCode
end

function PCFactorGetMatSolverType(arg1, arg2)
    @ccall libpetsc.PCFactorGetMatSolverType(arg1::PC, arg2::Ptr{MatSolverType})::PetscErrorCode
end

function PCFactorSetUpMatSolverType(arg1)
    @ccall libpetsc.PCFactorSetUpMatSolverType(arg1::PC)::PetscErrorCode
end

function PCFactorSetMatSolverPackage(pc, stype)
    @ccall libpetsc.PCFactorSetMatSolverPackage(pc::PC, stype::MatSolverType)::PetscErrorCode
end

function PCFactorGetMatSolverPackage(pc, stype)
    @ccall libpetsc.PCFactorGetMatSolverPackage(pc::PC, stype::Ptr{MatSolverType})::PetscErrorCode
end

function PCFactorSetUpMatSolverPackage(pc)
    @ccall libpetsc.PCFactorSetUpMatSolverPackage(pc::PC)::PetscErrorCode
end

function PCFactorSetFill(arg1, arg2)
    @ccall libpetsc.PCFactorSetFill(arg1::PC, arg2::PetscReal)::PetscErrorCode
end

function PCFactorSetColumnPivot(arg1, arg2)
    @ccall libpetsc.PCFactorSetColumnPivot(arg1::PC, arg2::PetscReal)::PetscErrorCode
end

function PCFactorReorderForNonzeroDiagonal(arg1, arg2)
    @ccall libpetsc.PCFactorReorderForNonzeroDiagonal(arg1::PC, arg2::PetscReal)::PetscErrorCode
end

function PCFactorSetMatOrderingType(arg1, arg2)
    @ccall libpetsc.PCFactorSetMatOrderingType(arg1::PC, arg2::MatOrderingType)::PetscErrorCode
end

function PCFactorSetReuseOrdering(arg1, arg2)
    @ccall libpetsc.PCFactorSetReuseOrdering(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCFactorSetReuseFill(arg1, arg2)
    @ccall libpetsc.PCFactorSetReuseFill(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCFactorSetUseInPlace(arg1, arg2)
    @ccall libpetsc.PCFactorSetUseInPlace(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCFactorGetUseInPlace(arg1, arg2)
    @ccall libpetsc.PCFactorGetUseInPlace(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCFactorSetAllowDiagonalFill(arg1, arg2)
    @ccall libpetsc.PCFactorSetAllowDiagonalFill(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCFactorGetAllowDiagonalFill(arg1, arg2)
    @ccall libpetsc.PCFactorGetAllowDiagonalFill(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCFactorSetPivotInBlocks(arg1, arg2)
    @ccall libpetsc.PCFactorSetPivotInBlocks(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCFactorSetLevels(arg1, arg2)
    @ccall libpetsc.PCFactorSetLevels(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCFactorGetLevels(arg1, arg2)
    @ccall libpetsc.PCFactorGetLevels(arg1::PC, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PCFactorSetDropTolerance(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PCFactorSetDropTolerance(arg1::PC, arg2::PetscReal, arg3::PetscReal, arg4::PetscInt)::PetscErrorCode
end

function PCFactorGetZeroPivot(arg1, arg2)
    @ccall libpetsc.PCFactorGetZeroPivot(arg1::PC, arg2::Ptr{PetscReal})::PetscErrorCode
end

function PCFactorGetShiftAmount(arg1, arg2)
    @ccall libpetsc.PCFactorGetShiftAmount(arg1::PC, arg2::Ptr{PetscReal})::PetscErrorCode
end

function PCFactorGetShiftType(arg1, arg2)
    @ccall libpetsc.PCFactorGetShiftType(arg1::PC, arg2::Ptr{MatFactorShiftType})::PetscErrorCode
end

function PCASMSetLocalSubdomains(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PCASMSetLocalSubdomains(arg1::PC, arg2::PetscInt, arg3::Ptr{IS}, arg4::Ptr{IS})::PetscErrorCode
end

function PCASMSetTotalSubdomains(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PCASMSetTotalSubdomains(arg1::PC, arg2::PetscInt, arg3::Ptr{IS}, arg4::Ptr{IS})::PetscErrorCode
end

function PCASMSetOverlap(arg1, arg2)
    @ccall libpetsc.PCASMSetOverlap(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCASMSetDMSubdomains(arg1, arg2)
    @ccall libpetsc.PCASMSetDMSubdomains(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCASMGetDMSubdomains(arg1, arg2)
    @ccall libpetsc.PCASMGetDMSubdomains(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCASMSetSortIndices(arg1, arg2)
    @ccall libpetsc.PCASMSetSortIndices(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCASMSetType(arg1, arg2)
    @ccall libpetsc.PCASMSetType(arg1::PC, arg2::PCASMType)::PetscErrorCode
end

function PCASMGetType(arg1, arg2)
    @ccall libpetsc.PCASMGetType(arg1::PC, arg2::Ptr{PCASMType})::PetscErrorCode
end

function PCASMSetLocalType(arg1, arg2)
    @ccall libpetsc.PCASMSetLocalType(arg1::PC, arg2::PCCompositeType)::PetscErrorCode
end

function PCASMGetLocalType(arg1, arg2)
    @ccall libpetsc.PCASMGetLocalType(arg1::PC, arg2::Ptr{PCCompositeType})::PetscErrorCode
end

function PCASMCreateSubdomains(arg1, arg2, arg3)
    @ccall libpetsc.PCASMCreateSubdomains(arg1::Mat, arg2::PetscInt, arg3::Ptr{Ptr{IS}})::PetscErrorCode
end

function PCASMDestroySubdomains(arg1, arg2, arg3)
    @ccall libpetsc.PCASMDestroySubdomains(arg1::PetscInt, arg2::Ptr{IS}, arg3::Ptr{IS})::PetscErrorCode
end

function PCASMCreateSubdomains2D(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.PCASMCreateSubdomains2D(arg1::PetscInt, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::Ptr{PetscInt}, arg8::Ptr{Ptr{IS}}, arg9::Ptr{Ptr{IS}})::PetscErrorCode
end

function PCASMGetLocalSubdomains(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PCASMGetLocalSubdomains(arg1::PC, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{IS}}, arg4::Ptr{Ptr{IS}})::PetscErrorCode
end

function PCASMGetLocalSubmatrices(arg1, arg2, arg3)
    @ccall libpetsc.PCASMGetLocalSubmatrices(arg1::PC, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{Mat}})::PetscErrorCode
end

function PCASMGetSubMatType(arg1, arg2)
    @ccall libpetsc.PCASMGetSubMatType(arg1::PC, arg2::Ptr{MatType})::PetscErrorCode
end

function PCASMSetSubMatType(arg1, arg2)
    @ccall libpetsc.PCASMSetSubMatType(arg1::PC, arg2::MatType)::PetscErrorCode
end

function PCGASMSetTotalSubdomains(arg1, arg2)
    @ccall libpetsc.PCGASMSetTotalSubdomains(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCGASMSetSubdomains(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PCGASMSetSubdomains(arg1::PC, arg2::PetscInt, arg3::Ptr{IS}, arg4::Ptr{IS})::PetscErrorCode
end

function PCGASMSetOverlap(arg1, arg2)
    @ccall libpetsc.PCGASMSetOverlap(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCGASMSetUseDMSubdomains(arg1, arg2)
    @ccall libpetsc.PCGASMSetUseDMSubdomains(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCGASMGetUseDMSubdomains(arg1, arg2)
    @ccall libpetsc.PCGASMGetUseDMSubdomains(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCGASMSetSortIndices(arg1, arg2)
    @ccall libpetsc.PCGASMSetSortIndices(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCGASMSetType(arg1, arg2)
    @ccall libpetsc.PCGASMSetType(arg1::PC, arg2::PCGASMType)::PetscErrorCode
end

function PCGASMCreateSubdomains(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PCGASMCreateSubdomains(arg1::Mat, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{IS}})::PetscErrorCode
end

function PCGASMDestroySubdomains(arg1, arg2, arg3)
    @ccall libpetsc.PCGASMDestroySubdomains(arg1::PetscInt, arg2::Ptr{Ptr{IS}}, arg3::Ptr{Ptr{IS}})::PetscErrorCode
end

function PCGASMCreateSubdomains2D(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.PCGASMCreateSubdomains2D(arg1::PC, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::PetscInt, arg7::PetscInt, arg8::Ptr{PetscInt}, arg9::Ptr{Ptr{IS}}, arg10::Ptr{Ptr{IS}})::PetscErrorCode
end

function PCGASMGetSubdomains(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PCGASMGetSubdomains(arg1::PC, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{IS}}, arg4::Ptr{Ptr{IS}})::PetscErrorCode
end

function PCGASMGetSubmatrices(arg1, arg2, arg3)
    @ccall libpetsc.PCGASMGetSubmatrices(arg1::PC, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{Mat}})::PetscErrorCode
end

function PCCompositeSetType(arg1, arg2)
    @ccall libpetsc.PCCompositeSetType(arg1::PC, arg2::PCCompositeType)::PetscErrorCode
end

function PCCompositeGetType(arg1, arg2)
    @ccall libpetsc.PCCompositeGetType(arg1::PC, arg2::Ptr{PCCompositeType})::PetscErrorCode
end

function PCCompositeAddPCType(arg1, arg2)
    @ccall libpetsc.PCCompositeAddPCType(arg1::PC, arg2::PCType)::PetscErrorCode
end

function PCCompositeAddPC(arg1, arg2)
    @ccall libpetsc.PCCompositeAddPC(arg1::PC, arg2::PC)::PetscErrorCode
end

function PCCompositeGetNumberPC(arg1, arg2)
    @ccall libpetsc.PCCompositeGetNumberPC(arg1::PC, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PCCompositeGetPC(arg1, arg2, arg3)
    @ccall libpetsc.PCCompositeGetPC(arg1::PC, arg2::PetscInt, arg3::Ptr{PC})::PetscErrorCode
end

function PCCompositeSpecialSetAlpha(arg1, arg2)
    @ccall libpetsc.PCCompositeSpecialSetAlpha(arg1::PC, arg2::PetscScalar)::PetscErrorCode
end

function PCRedundantSetNumber(arg1, arg2)
    @ccall libpetsc.PCRedundantSetNumber(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCRedundantSetScatter(arg1, arg2, arg3)
    @ccall libpetsc.PCRedundantSetScatter(arg1::PC, arg2::VecScatter, arg3::VecScatter)::PetscErrorCode
end

function PCRedundantGetOperators(arg1, arg2, arg3)
    @ccall libpetsc.PCRedundantGetOperators(arg1::PC, arg2::Ptr{Mat}, arg3::Ptr{Mat})::PetscErrorCode
end

function PCSPAISetEpsilon(arg1, arg2)
    @ccall libpetsc.PCSPAISetEpsilon(arg1::PC, arg2::Cdouble)::PetscErrorCode
end

function PCSPAISetNBSteps(arg1, arg2)
    @ccall libpetsc.PCSPAISetNBSteps(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCSPAISetMax(arg1, arg2)
    @ccall libpetsc.PCSPAISetMax(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCSPAISetMaxNew(arg1, arg2)
    @ccall libpetsc.PCSPAISetMaxNew(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCSPAISetBlockSize(arg1, arg2)
    @ccall libpetsc.PCSPAISetBlockSize(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCSPAISetCacheSize(arg1, arg2)
    @ccall libpetsc.PCSPAISetCacheSize(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCSPAISetVerbose(arg1, arg2)
    @ccall libpetsc.PCSPAISetVerbose(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCSPAISetSp(arg1, arg2)
    @ccall libpetsc.PCSPAISetSp(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCHYPRESetType(arg1, arg2)
    @ccall libpetsc.PCHYPRESetType(arg1::PC, arg2::Ptr{Cchar})::PetscErrorCode
end

function PCHYPREGetType(arg1, arg2)
    @ccall libpetsc.PCHYPREGetType(arg1::PC, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function PCHYPRESetDiscreteGradient(arg1, arg2)
    @ccall libpetsc.PCHYPRESetDiscreteGradient(arg1::PC, arg2::Mat)::PetscErrorCode
end

function PCHYPRESetDiscreteCurl(arg1, arg2)
    @ccall libpetsc.PCHYPRESetDiscreteCurl(arg1::PC, arg2::Mat)::PetscErrorCode
end

function PCHYPRESetInterpolations(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PCHYPRESetInterpolations(arg1::PC, arg2::PetscInt, arg3::Mat, arg4::Ptr{Mat}, arg5::Mat, arg6::Ptr{Mat})::PetscErrorCode
end

function PCHYPRESetEdgeConstantVectors(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PCHYPRESetEdgeConstantVectors(arg1::PC, arg2::Vec, arg3::Vec, arg4::Vec)::PetscErrorCode
end

function PCHYPRESetAlphaPoissonMatrix(arg1, arg2)
    @ccall libpetsc.PCHYPRESetAlphaPoissonMatrix(arg1::PC, arg2::Mat)::PetscErrorCode
end

function PCHYPRESetBetaPoissonMatrix(arg1, arg2)
    @ccall libpetsc.PCHYPRESetBetaPoissonMatrix(arg1::PC, arg2::Mat)::PetscErrorCode
end

function PCFieldSplitSetFields(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PCFieldSplitSetFields(arg1::PC, arg2::Ptr{Cchar}, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt})::PetscErrorCode
end

function PCFieldSplitSetType(arg1, arg2)
    @ccall libpetsc.PCFieldSplitSetType(arg1::PC, arg2::PCCompositeType)::PetscErrorCode
end

function PCFieldSplitGetType(arg1, arg2)
    @ccall libpetsc.PCFieldSplitGetType(arg1::PC, arg2::Ptr{PCCompositeType})::PetscErrorCode
end

function PCFieldSplitSetBlockSize(arg1, arg2)
    @ccall libpetsc.PCFieldSplitSetBlockSize(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCFieldSplitSetIS(arg1, arg2, arg3)
    @ccall libpetsc.PCFieldSplitSetIS(arg1::PC, arg2::Ptr{Cchar}, arg3::IS)::PetscErrorCode
end

function PCFieldSplitGetIS(arg1, arg2, arg3)
    @ccall libpetsc.PCFieldSplitGetIS(arg1::PC, arg2::Ptr{Cchar}, arg3::Ptr{IS})::PetscErrorCode
end

function PCFieldSplitGetISByIndex(arg1, arg2, arg3)
    @ccall libpetsc.PCFieldSplitGetISByIndex(arg1::PC, arg2::PetscInt, arg3::Ptr{IS})::PetscErrorCode
end

function PCFieldSplitRestrictIS(arg1, arg2)
    @ccall libpetsc.PCFieldSplitRestrictIS(arg1::PC, arg2::IS)::PetscErrorCode
end

function PCFieldSplitSetDMSplits(arg1, arg2)
    @ccall libpetsc.PCFieldSplitSetDMSplits(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCFieldSplitGetDMSplits(arg1, arg2)
    @ccall libpetsc.PCFieldSplitGetDMSplits(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCFieldSplitSetDiagUseAmat(arg1, arg2)
    @ccall libpetsc.PCFieldSplitSetDiagUseAmat(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCFieldSplitGetDiagUseAmat(arg1, arg2)
    @ccall libpetsc.PCFieldSplitGetDiagUseAmat(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCFieldSplitSetOffDiagUseAmat(arg1, arg2)
    @ccall libpetsc.PCFieldSplitSetOffDiagUseAmat(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCFieldSplitGetOffDiagUseAmat(arg1, arg2)
    @ccall libpetsc.PCFieldSplitGetOffDiagUseAmat(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCFieldSplitSchurPrecondition(arg1, arg2, arg3)
    @ccall libpetsc.PCFieldSplitSchurPrecondition(arg1::PC, arg2::PCFieldSplitSchurPreType, arg3::Mat)::PetscErrorCode
end

function PCFieldSplitSetSchurPre(arg1, arg2, arg3)
    @ccall libpetsc.PCFieldSplitSetSchurPre(arg1::PC, arg2::PCFieldSplitSchurPreType, arg3::Mat)::PetscErrorCode
end

function PCFieldSplitGetSchurPre(arg1, arg2, arg3)
    @ccall libpetsc.PCFieldSplitGetSchurPre(arg1::PC, arg2::Ptr{PCFieldSplitSchurPreType}, arg3::Ptr{Mat})::PetscErrorCode
end

function PCFieldSplitSetSchurFactType(arg1, arg2)
    @ccall libpetsc.PCFieldSplitSetSchurFactType(arg1::PC, arg2::PCFieldSplitSchurFactType)::PetscErrorCode
end

function PCFieldSplitSetSchurScale(arg1, arg2)
    @ccall libpetsc.PCFieldSplitSetSchurScale(arg1::PC, arg2::PetscScalar)::PetscErrorCode
end

function PCFieldSplitGetSchurBlocks(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PCFieldSplitGetSchurBlocks(arg1::PC, arg2::Ptr{Mat}, arg3::Ptr{Mat}, arg4::Ptr{Mat}, arg5::Ptr{Mat})::PetscErrorCode
end

function PCFieldSplitSchurGetS(arg1, S)
    @ccall libpetsc.PCFieldSplitSchurGetS(arg1::PC, S::Ptr{Mat})::PetscErrorCode
end

function PCFieldSplitSchurRestoreS(arg1, S)
    @ccall libpetsc.PCFieldSplitSchurRestoreS(arg1::PC, S::Ptr{Mat})::PetscErrorCode
end

function PCFieldSplitGetDetectSaddlePoint(arg1, arg2)
    @ccall libpetsc.PCFieldSplitGetDetectSaddlePoint(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCFieldSplitSetDetectSaddlePoint(arg1, arg2)
    @ccall libpetsc.PCFieldSplitSetDetectSaddlePoint(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCFieldSplitSetGKBTol(arg1, arg2)
    @ccall libpetsc.PCFieldSplitSetGKBTol(arg1::PC, arg2::PetscReal)::PetscErrorCode
end

function PCFieldSplitSetGKBNu(arg1, arg2)
    @ccall libpetsc.PCFieldSplitSetGKBNu(arg1::PC, arg2::PetscReal)::PetscErrorCode
end

function PCFieldSplitSetGKBMaxit(arg1, arg2)
    @ccall libpetsc.PCFieldSplitSetGKBMaxit(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCFieldSplitSetGKBDelay(arg1, arg2)
    @ccall libpetsc.PCFieldSplitSetGKBDelay(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCGalerkinSetRestriction(arg1, arg2)
    @ccall libpetsc.PCGalerkinSetRestriction(arg1::PC, arg2::Mat)::PetscErrorCode
end

function PCGalerkinSetInterpolation(arg1, arg2)
    @ccall libpetsc.PCGalerkinSetInterpolation(arg1::PC, arg2::Mat)::PetscErrorCode
end

function PCGalerkinSetComputeSubmatrix(arg1, arg2, arg3)
    @ccall libpetsc.PCGalerkinSetComputeSubmatrix(arg1::PC, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function PCPythonSetType(arg1, arg2)
    @ccall libpetsc.PCPythonSetType(arg1::PC, arg2::Ptr{Cchar})::PetscErrorCode
end

function PCPARMSSetGlobal(arg1, arg2)
    @ccall libpetsc.PCPARMSSetGlobal(arg1::PC, arg2::PCPARMSGlobalType)::PetscErrorCode
end

function PCPARMSSetLocal(arg1, arg2)
    @ccall libpetsc.PCPARMSSetLocal(arg1::PC, arg2::PCPARMSLocalType)::PetscErrorCode
end

function PCPARMSSetSolveTolerances(arg1, arg2, arg3)
    @ccall libpetsc.PCPARMSSetSolveTolerances(arg1::PC, arg2::PetscReal, arg3::PetscInt)::PetscErrorCode
end

function PCPARMSSetSolveRestart(arg1, arg2)
    @ccall libpetsc.PCPARMSSetSolveRestart(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCPARMSSetNonsymPerm(arg1, arg2)
    @ccall libpetsc.PCPARMSSetNonsymPerm(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCPARMSSetFill(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PCPARMSSetFill(arg1::PC, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt)::PetscErrorCode
end

function PCGAMGSetType(arg1, arg2)
    @ccall libpetsc.PCGAMGSetType(arg1::PC, arg2::PCGAMGType)::PetscErrorCode
end

function PCGAMGGetType(arg1, arg2)
    @ccall libpetsc.PCGAMGGetType(arg1::PC, arg2::Ptr{PCGAMGType})::PetscErrorCode
end

function PCGAMGSetProcEqLim(arg1, arg2)
    @ccall libpetsc.PCGAMGSetProcEqLim(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCGAMGSetRepartition(arg1, arg2)
    @ccall libpetsc.PCGAMGSetRepartition(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCGAMGSetUseSAEstEig(arg1, arg2)
    @ccall libpetsc.PCGAMGSetUseSAEstEig(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCGAMGSetEstEigKSPMaxIt(arg1, arg2)
    @ccall libpetsc.PCGAMGSetEstEigKSPMaxIt(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCGAMGSetEstEigKSPType(arg1, arg2)
    @ccall libpetsc.PCGAMGSetEstEigKSPType(arg1::PC, arg2::Ptr{Cchar})::PetscErrorCode
end

function PCGAMGSetEigenvalues(arg1, arg2, arg3)
    @ccall libpetsc.PCGAMGSetEigenvalues(arg1::PC, arg2::PetscReal, arg3::PetscReal)::PetscErrorCode
end

function PCGAMGASMSetUseAggs(arg1, arg2)
    @ccall libpetsc.PCGAMGASMSetUseAggs(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCGAMGSetUseParallelCoarseGridSolve(arg1, arg2)
    @ccall libpetsc.PCGAMGSetUseParallelCoarseGridSolve(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCGAMGSetCpuPinCoarseGrids(arg1, arg2)
    @ccall libpetsc.PCGAMGSetCpuPinCoarseGrids(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCGAMGSetCoarseGridLayoutType(arg1, arg2)
    @ccall libpetsc.PCGAMGSetCoarseGridLayoutType(arg1::PC, arg2::PCGAMGLayoutType)::PetscErrorCode
end

function PCGAMGSetThreshold(arg1, arg2, arg3)
    @ccall libpetsc.PCGAMGSetThreshold(arg1::PC, arg2::Ptr{PetscReal}, arg3::PetscInt)::PetscErrorCode
end

function PCGAMGSetRankReductionFactors(arg1, arg2, arg3)
    @ccall libpetsc.PCGAMGSetRankReductionFactors(arg1::PC, arg2::Ptr{PetscInt}, arg3::PetscInt)::PetscErrorCode
end

function PCGAMGSetThresholdScale(arg1, arg2)
    @ccall libpetsc.PCGAMGSetThresholdScale(arg1::PC, arg2::PetscReal)::PetscErrorCode
end

function PCGAMGSetCoarseEqLim(arg1, arg2)
    @ccall libpetsc.PCGAMGSetCoarseEqLim(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCGAMGSetNlevels(arg1, arg2)
    @ccall libpetsc.PCGAMGSetNlevels(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCGAMGSetNSmooths(arg1, arg2)
    @ccall libpetsc.PCGAMGSetNSmooths(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCGAMGSetSymGraph(arg1, arg2)
    @ccall libpetsc.PCGAMGSetSymGraph(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCGAMGSetSquareGraph(arg1, arg2)
    @ccall libpetsc.PCGAMGSetSquareGraph(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCGAMGSetReuseInterpolation(arg1, arg2)
    @ccall libpetsc.PCGAMGSetReuseInterpolation(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCGAMGFinalizePackage()
    @ccall libpetsc.PCGAMGFinalizePackage()::PetscErrorCode
end

function PCGAMGInitializePackage()
    @ccall libpetsc.PCGAMGInitializePackage()::PetscErrorCode
end

function PCGAMGRegister(arg1, arg2)
    @ccall libpetsc.PCGAMGRegister(arg1::PCGAMGType, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PCGAMGClassicalSetType(arg1, arg2)
    @ccall libpetsc.PCGAMGClassicalSetType(arg1::PC, arg2::PCGAMGClassicalType)::PetscErrorCode
end

function PCGAMGClassicalGetType(arg1, arg2)
    @ccall libpetsc.PCGAMGClassicalGetType(arg1::PC, arg2::Ptr{PCGAMGClassicalType})::PetscErrorCode
end

function PCBDDCSetDiscreteGradient(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.PCBDDCSetDiscreteGradient(arg1::PC, arg2::Mat, arg3::PetscInt, arg4::PetscInt, arg5::PetscBool, arg6::PetscBool)::PetscErrorCode
end

function PCBDDCSetDivergenceMat(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PCBDDCSetDivergenceMat(arg1::PC, arg2::Mat, arg3::PetscBool, arg4::IS)::PetscErrorCode
end

function PCBDDCSetChangeOfBasisMat(arg1, arg2, arg3)
    @ccall libpetsc.PCBDDCSetChangeOfBasisMat(arg1::PC, arg2::Mat, arg3::PetscBool)::PetscErrorCode
end

function PCBDDCSetPrimalVerticesIS(arg1, arg2)
    @ccall libpetsc.PCBDDCSetPrimalVerticesIS(arg1::PC, arg2::IS)::PetscErrorCode
end

function PCBDDCSetPrimalVerticesLocalIS(arg1, arg2)
    @ccall libpetsc.PCBDDCSetPrimalVerticesLocalIS(arg1::PC, arg2::IS)::PetscErrorCode
end

function PCBDDCGetPrimalVerticesIS(arg1, arg2)
    @ccall libpetsc.PCBDDCGetPrimalVerticesIS(arg1::PC, arg2::Ptr{IS})::PetscErrorCode
end

function PCBDDCGetPrimalVerticesLocalIS(arg1, arg2)
    @ccall libpetsc.PCBDDCGetPrimalVerticesLocalIS(arg1::PC, arg2::Ptr{IS})::PetscErrorCode
end

function PCBDDCSetCoarseningRatio(arg1, arg2)
    @ccall libpetsc.PCBDDCSetCoarseningRatio(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCBDDCSetLevels(arg1, arg2)
    @ccall libpetsc.PCBDDCSetLevels(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCBDDCSetDirichletBoundaries(arg1, arg2)
    @ccall libpetsc.PCBDDCSetDirichletBoundaries(arg1::PC, arg2::IS)::PetscErrorCode
end

function PCBDDCSetDirichletBoundariesLocal(arg1, arg2)
    @ccall libpetsc.PCBDDCSetDirichletBoundariesLocal(arg1::PC, arg2::IS)::PetscErrorCode
end

function PCBDDCGetDirichletBoundaries(arg1, arg2)
    @ccall libpetsc.PCBDDCGetDirichletBoundaries(arg1::PC, arg2::Ptr{IS})::PetscErrorCode
end

function PCBDDCGetDirichletBoundariesLocal(arg1, arg2)
    @ccall libpetsc.PCBDDCGetDirichletBoundariesLocal(arg1::PC, arg2::Ptr{IS})::PetscErrorCode
end

function PCBDDCSetInterfaceExtType(arg1, arg2)
    @ccall libpetsc.PCBDDCSetInterfaceExtType(arg1::PC, arg2::PCBDDCInterfaceExtType)::PetscErrorCode
end

function PCBDDCSetNeumannBoundaries(arg1, arg2)
    @ccall libpetsc.PCBDDCSetNeumannBoundaries(arg1::PC, arg2::IS)::PetscErrorCode
end

function PCBDDCSetNeumannBoundariesLocal(arg1, arg2)
    @ccall libpetsc.PCBDDCSetNeumannBoundariesLocal(arg1::PC, arg2::IS)::PetscErrorCode
end

function PCBDDCGetNeumannBoundaries(arg1, arg2)
    @ccall libpetsc.PCBDDCGetNeumannBoundaries(arg1::PC, arg2::Ptr{IS})::PetscErrorCode
end

function PCBDDCGetNeumannBoundariesLocal(arg1, arg2)
    @ccall libpetsc.PCBDDCGetNeumannBoundariesLocal(arg1::PC, arg2::Ptr{IS})::PetscErrorCode
end

function PCBDDCSetDofsSplitting(arg1, arg2, arg3)
    @ccall libpetsc.PCBDDCSetDofsSplitting(arg1::PC, arg2::PetscInt, arg3::Ptr{IS})::PetscErrorCode
end

function PCBDDCSetDofsSplittingLocal(arg1, arg2, arg3)
    @ccall libpetsc.PCBDDCSetDofsSplittingLocal(arg1::PC, arg2::PetscInt, arg3::Ptr{IS})::PetscErrorCode
end

function PCBDDCSetLocalAdjacencyGraph(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PCBDDCSetLocalAdjacencyGraph(arg1::PC, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt}, arg5::PetscCopyMode)::PetscErrorCode
end

function PCBDDCCreateFETIDPOperators(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PCBDDCCreateFETIDPOperators(arg1::PC, arg2::PetscBool, arg3::Ptr{Cchar}, arg4::Ptr{Mat}, arg5::Ptr{PC})::PetscErrorCode
end

function PCBDDCMatFETIDPGetRHS(arg1, arg2, arg3)
    @ccall libpetsc.PCBDDCMatFETIDPGetRHS(arg1::Mat, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function PCBDDCMatFETIDPGetSolution(arg1, arg2, arg3)
    @ccall libpetsc.PCBDDCMatFETIDPGetSolution(arg1::Mat, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function PCBDDCFinalizePackage()
    @ccall libpetsc.PCBDDCFinalizePackage()::PetscErrorCode
end

function PCBDDCInitializePackage()
    @ccall libpetsc.PCBDDCInitializePackage()::PetscErrorCode
end

function PCISSetUseStiffnessScaling(arg1, arg2)
    @ccall libpetsc.PCISSetUseStiffnessScaling(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCISSetSubdomainScalingFactor(arg1, arg2)
    @ccall libpetsc.PCISSetSubdomainScalingFactor(arg1::PC, arg2::PetscScalar)::PetscErrorCode
end

function PCISSetSubdomainDiagonalScaling(arg1, arg2)
    @ccall libpetsc.PCISSetSubdomainDiagonalScaling(arg1::PC, arg2::Vec)::PetscErrorCode
end

function PCMGSetType(arg1, arg2)
    @ccall libpetsc.PCMGSetType(arg1::PC, arg2::PCMGType)::PetscErrorCode
end

function PCMGGetType(arg1, arg2)
    @ccall libpetsc.PCMGGetType(arg1::PC, arg2::Ptr{PCMGType})::PetscErrorCode
end

function PCMGSetLevels(arg1, arg2, arg3)
    @ccall libpetsc.PCMGSetLevels(arg1::PC, arg2::PetscInt, arg3::Ptr{MPI_Comm})::PetscErrorCode
end

function PCMGGetLevels(arg1, arg2)
    @ccall libpetsc.PCMGGetLevels(arg1::PC, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PCMGSetDistinctSmoothUp(arg1)
    @ccall libpetsc.PCMGSetDistinctSmoothUp(arg1::PC)::PetscErrorCode
end

function PCMGSetNumberSmooth(arg1, arg2)
    @ccall libpetsc.PCMGSetNumberSmooth(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCMGSetCycleType(arg1, arg2)
    @ccall libpetsc.PCMGSetCycleType(arg1::PC, arg2::PCMGCycleType)::PetscErrorCode
end

function PCMGSetCycleTypeOnLevel(arg1, arg2, arg3)
    @ccall libpetsc.PCMGSetCycleTypeOnLevel(arg1::PC, arg2::PetscInt, arg3::PCMGCycleType)::PetscErrorCode
end

function PCMGSetCyclesOnLevel(pc, l, t)
    @ccall libpetsc.PCMGSetCyclesOnLevel(pc::PC, l::PetscInt, t::PetscInt)::PetscErrorCode
end

function PCMGMultiplicativeSetCycles(arg1, arg2)
    @ccall libpetsc.PCMGMultiplicativeSetCycles(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCMGSetGalerkin(arg1, arg2)
    @ccall libpetsc.PCMGSetGalerkin(arg1::PC, arg2::PCMGGalerkinType)::PetscErrorCode
end

function PCMGGetGalerkin(arg1, arg2)
    @ccall libpetsc.PCMGGetGalerkin(arg1::PC, arg2::Ptr{PCMGGalerkinType})::PetscErrorCode
end

function PCMGSetAdaptInterpolation(arg1, arg2)
    @ccall libpetsc.PCMGSetAdaptInterpolation(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCMGGetAdaptInterpolation(arg1, arg2)
    @ccall libpetsc.PCMGGetAdaptInterpolation(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCMGSetAdaptCR(arg1, arg2)
    @ccall libpetsc.PCMGSetAdaptCR(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCMGGetAdaptCR(arg1, arg2)
    @ccall libpetsc.PCMGGetAdaptCR(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCMGSetRhs(arg1, arg2, arg3)
    @ccall libpetsc.PCMGSetRhs(arg1::PC, arg2::PetscInt, arg3::Vec)::PetscErrorCode
end

function PCMGSetX(arg1, arg2, arg3)
    @ccall libpetsc.PCMGSetX(arg1::PC, arg2::PetscInt, arg3::Vec)::PetscErrorCode
end

function PCMGSetR(arg1, arg2, arg3)
    @ccall libpetsc.PCMGSetR(arg1::PC, arg2::PetscInt, arg3::Vec)::PetscErrorCode
end

function PCMGSetRestriction(arg1, arg2, arg3)
    @ccall libpetsc.PCMGSetRestriction(arg1::PC, arg2::PetscInt, arg3::Mat)::PetscErrorCode
end

function PCMGGetRestriction(arg1, arg2, arg3)
    @ccall libpetsc.PCMGGetRestriction(arg1::PC, arg2::PetscInt, arg3::Ptr{Mat})::PetscErrorCode
end

function PCMGSetInjection(arg1, arg2, arg3)
    @ccall libpetsc.PCMGSetInjection(arg1::PC, arg2::PetscInt, arg3::Mat)::PetscErrorCode
end

function PCMGGetInjection(arg1, arg2, arg3)
    @ccall libpetsc.PCMGGetInjection(arg1::PC, arg2::PetscInt, arg3::Ptr{Mat})::PetscErrorCode
end

function PCMGSetInterpolation(arg1, arg2, arg3)
    @ccall libpetsc.PCMGSetInterpolation(arg1::PC, arg2::PetscInt, arg3::Mat)::PetscErrorCode
end

function PCMGSetOperators(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PCMGSetOperators(arg1::PC, arg2::PetscInt, arg3::Mat, arg4::Mat)::PetscErrorCode
end

function PCMGGetInterpolation(arg1, arg2, arg3)
    @ccall libpetsc.PCMGGetInterpolation(arg1::PC, arg2::PetscInt, arg3::Ptr{Mat})::PetscErrorCode
end

function PCMGSetRScale(arg1, arg2, arg3)
    @ccall libpetsc.PCMGSetRScale(arg1::PC, arg2::PetscInt, arg3::Vec)::PetscErrorCode
end

function PCMGGetRScale(arg1, arg2, arg3)
    @ccall libpetsc.PCMGGetRScale(arg1::PC, arg2::PetscInt, arg3::Ptr{Vec})::PetscErrorCode
end

function PCMGSetResidual(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PCMGSetResidual(arg1::PC, arg2::PetscInt, arg3::Ptr{Cvoid}, arg4::Mat)::PetscErrorCode
end

function PCMGSetResidualTranspose(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PCMGSetResidualTranspose(arg1::PC, arg2::PetscInt, arg3::Ptr{Cvoid}, arg4::Mat)::PetscErrorCode
end

function PCMGResidualDefault(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PCMGResidualDefault(arg1::Mat, arg2::Vec, arg3::Vec, arg4::Vec)::PetscErrorCode
end

function PCMGResidualTransposeDefault(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PCMGResidualTransposeDefault(arg1::Mat, arg2::Vec, arg3::Vec, arg4::Vec)::PetscErrorCode
end

function PCMGMatResidualDefault(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PCMGMatResidualDefault(arg1::Mat, arg2::Mat, arg3::Mat, arg4::Mat)::PetscErrorCode
end

function PCMGMatResidualTransposeDefault(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PCMGMatResidualTransposeDefault(arg1::Mat, arg2::Mat, arg3::Mat, arg4::Mat)::PetscErrorCode
end

function PCHMGSetReuseInterpolation(arg1, arg2)
    @ccall libpetsc.PCHMGSetReuseInterpolation(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCHMGSetUseSubspaceCoarsening(arg1, arg2)
    @ccall libpetsc.PCHMGSetUseSubspaceCoarsening(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCHMGSetInnerPCType(arg1, arg2)
    @ccall libpetsc.PCHMGSetInnerPCType(arg1::PC, arg2::PCType)::PetscErrorCode
end

function PCHMGSetCoarseningComponent(arg1, arg2)
    @ccall libpetsc.PCHMGSetCoarseningComponent(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCHMGUseMatMAIJ(arg1, arg2)
    @ccall libpetsc.PCHMGUseMatMAIJ(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCTelescopeGetSubcommType(arg1, arg2)
    @ccall libpetsc.PCTelescopeGetSubcommType(arg1::PC, arg2::Ptr{PetscSubcommType})::PetscErrorCode
end

function PCTelescopeSetSubcommType(arg1, arg2)
    @ccall libpetsc.PCTelescopeSetSubcommType(arg1::PC, arg2::PetscSubcommType)::PetscErrorCode
end

function PCTelescopeGetReductionFactor(arg1, arg2)
    @ccall libpetsc.PCTelescopeGetReductionFactor(arg1::PC, arg2::Ptr{PetscInt})::PetscErrorCode
end

function PCTelescopeSetReductionFactor(arg1, arg2)
    @ccall libpetsc.PCTelescopeSetReductionFactor(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCTelescopeGetIgnoreDM(arg1, arg2)
    @ccall libpetsc.PCTelescopeGetIgnoreDM(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCTelescopeSetIgnoreDM(arg1, arg2)
    @ccall libpetsc.PCTelescopeSetIgnoreDM(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCTelescopeGetUseCoarseDM(arg1, arg2)
    @ccall libpetsc.PCTelescopeGetUseCoarseDM(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCTelescopeSetUseCoarseDM(arg1, arg2)
    @ccall libpetsc.PCTelescopeSetUseCoarseDM(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCTelescopeGetIgnoreKSPComputeOperators(arg1, arg2)
    @ccall libpetsc.PCTelescopeGetIgnoreKSPComputeOperators(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCTelescopeSetIgnoreKSPComputeOperators(arg1, arg2)
    @ccall libpetsc.PCTelescopeSetIgnoreKSPComputeOperators(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCTelescopeGetDM(arg1, arg2)
    @ccall libpetsc.PCTelescopeGetDM(arg1::PC, arg2::Ptr{DM})::PetscErrorCode
end

function PCPatchSetSaveOperators(arg1, arg2)
    @ccall libpetsc.PCPatchSetSaveOperators(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCPatchGetSaveOperators(arg1, arg2)
    @ccall libpetsc.PCPatchGetSaveOperators(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCPatchSetPrecomputeElementTensors(arg1, arg2)
    @ccall libpetsc.PCPatchSetPrecomputeElementTensors(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCPatchGetPrecomputeElementTensors(arg1, arg2)
    @ccall libpetsc.PCPatchGetPrecomputeElementTensors(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCPatchSetPartitionOfUnity(arg1, arg2)
    @ccall libpetsc.PCPatchSetPartitionOfUnity(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCPatchGetPartitionOfUnity(arg1, arg2)
    @ccall libpetsc.PCPatchGetPartitionOfUnity(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCPatchSetSubMatType(arg1, arg2)
    @ccall libpetsc.PCPatchSetSubMatType(arg1::PC, arg2::MatType)::PetscErrorCode
end

function PCPatchGetSubMatType(arg1, arg2)
    @ccall libpetsc.PCPatchGetSubMatType(arg1::PC, arg2::Ptr{MatType})::PetscErrorCode
end

function PCPatchSetCellNumbering(arg1, arg2)
    @ccall libpetsc.PCPatchSetCellNumbering(arg1::PC, arg2::PetscSection)::PetscErrorCode
end

function PCPatchGetCellNumbering(arg1, arg2)
    @ccall libpetsc.PCPatchGetCellNumbering(arg1::PC, arg2::Ptr{PetscSection})::PetscErrorCode
end

function PCPatchSetConstructType(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PCPatchSetConstructType(arg1::PC, arg2::PCPatchConstructType, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function PCPatchGetConstructType(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PCPatchGetConstructType(arg1::PC, arg2::Ptr{PCPatchConstructType}, arg3::Ptr{Ptr{Cvoid}}, arg4::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function PCPatchSetDiscretisationInfo(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11)
    @ccall libpetsc.PCPatchSetDiscretisationInfo(arg1::PC, arg2::PetscInt, arg3::Ptr{DM}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{Ptr{PetscInt}}, arg7::Ptr{PetscInt}, arg8::PetscInt, arg9::Ptr{PetscInt}, arg10::PetscInt, arg11::Ptr{PetscInt})::PetscErrorCode
end

function PCPatchSetComputeOperator(arg1, arg2, arg3)
    @ccall libpetsc.PCPatchSetComputeOperator(arg1::PC, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function PCPatchSetComputeFunction(pc, func, ctx)
    @ccall libpetsc.PCPatchSetComputeFunction(pc::PC, func::Ptr{Cvoid}, ctx::Ptr{Cvoid})::PetscErrorCode
end

function PCPatchSetComputeOperatorInteriorFacets(arg1, arg2, arg3)
    @ccall libpetsc.PCPatchSetComputeOperatorInteriorFacets(arg1::PC, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function PCPatchSetComputeFunctionInteriorFacets(pc, func, ctx)
    @ccall libpetsc.PCPatchSetComputeFunctionInteriorFacets(pc::PC, func::Ptr{Cvoid}, ctx::Ptr{Cvoid})::PetscErrorCode
end

function PCLMVMSetMatLMVM(arg1, arg2)
    @ccall libpetsc.PCLMVMSetMatLMVM(arg1::PC, arg2::Mat)::PetscErrorCode
end

function PCLMVMGetMatLMVM(arg1, arg2)
    @ccall libpetsc.PCLMVMGetMatLMVM(arg1::PC, arg2::Ptr{Mat})::PetscErrorCode
end

function PCLMVMSetIS(arg1, arg2)
    @ccall libpetsc.PCLMVMSetIS(arg1::PC, arg2::IS)::PetscErrorCode
end

function PCLMVMClearIS(arg1)
    @ccall libpetsc.PCLMVMClearIS(arg1::PC)::PetscErrorCode
end

function PCExoticSetType(arg1, arg2)
    @ccall libpetsc.PCExoticSetType(arg1::PC, arg2::PCExoticType)::PetscErrorCode
end

function PCDeflationSetInitOnly(arg1, arg2)
    @ccall libpetsc.PCDeflationSetInitOnly(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCDeflationSetLevels(arg1, arg2)
    @ccall libpetsc.PCDeflationSetLevels(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCDeflationSetReductionFactor(arg1, arg2)
    @ccall libpetsc.PCDeflationSetReductionFactor(arg1::PC, arg2::PetscInt)::PetscErrorCode
end

function PCDeflationSetCorrectionFactor(arg1, arg2)
    @ccall libpetsc.PCDeflationSetCorrectionFactor(arg1::PC, arg2::PetscScalar)::PetscErrorCode
end

function PCDeflationSetSpaceToCompute(arg1, arg2, arg3)
    @ccall libpetsc.PCDeflationSetSpaceToCompute(arg1::PC, arg2::PCDeflationSpaceType, arg3::PetscInt)::PetscErrorCode
end

function PCDeflationSetSpace(arg1, arg2, arg3)
    @ccall libpetsc.PCDeflationSetSpace(arg1::PC, arg2::Mat, arg3::PetscBool)::PetscErrorCode
end

function PCDeflationSetProjectionNullSpaceMat(arg1, arg2)
    @ccall libpetsc.PCDeflationSetProjectionNullSpaceMat(arg1::PC, arg2::Mat)::PetscErrorCode
end

function PCDeflationSetCoarseMat(arg1, arg2)
    @ccall libpetsc.PCDeflationSetCoarseMat(arg1::PC, arg2::Mat)::PetscErrorCode
end

function PCDeflationGetPC(arg1, arg2)
    @ccall libpetsc.PCDeflationGetPC(arg1::PC, arg2::Ptr{PC})::PetscErrorCode
end

function PCHPDDMSetAuxiliaryMat(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PCHPDDMSetAuxiliaryMat(arg1::PC, arg2::IS, arg3::Mat, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid})::PetscErrorCode
end

function PCHPDDMSetRHSMat(arg1, arg2)
    @ccall libpetsc.PCHPDDMSetRHSMat(arg1::PC, arg2::Mat)::PetscErrorCode
end

function PCHPDDMHasNeumannMat(arg1, arg2)
    @ccall libpetsc.PCHPDDMHasNeumannMat(arg1::PC, arg2::PetscBool)::PetscErrorCode
end

function PCHPDDMSetCoarseCorrectionType(arg1, arg2)
    @ccall libpetsc.PCHPDDMSetCoarseCorrectionType(arg1::PC, arg2::PCHPDDMCoarseCorrectionType)::PetscErrorCode
end

function PCHPDDMGetCoarseCorrectionType(arg1, arg2)
    @ccall libpetsc.PCHPDDMGetCoarseCorrectionType(arg1::PC, arg2::Ptr{PCHPDDMCoarseCorrectionType})::PetscErrorCode
end

function PCHPDDMGetSTShareSubKSP(arg1, arg2)
    @ccall libpetsc.PCHPDDMGetSTShareSubKSP(arg1::PC, arg2::Ptr{PetscBool})::PetscErrorCode
end

function PCHPDDMFinalizePackage()
    @ccall libpetsc.PCHPDDMFinalizePackage()::PetscErrorCode
end

function PCHPDDMInitializePackage()
    @ccall libpetsc.PCHPDDMInitializePackage()::PetscErrorCode
end

function KSPInitializePackage()
    @ccall libpetsc.KSPInitializePackage()::PetscErrorCode
end

mutable struct _p_KSP end

const KSP = Ptr{_p_KSP}

const KSPType = Ptr{Cchar}

function KSPCreate(arg1, arg2)
    @ccall libpetsc.KSPCreate(arg1::MPI_Comm, arg2::Ptr{KSP})::PetscErrorCode
end

function KSPSetType(arg1, arg2)
    @ccall libpetsc.KSPSetType(arg1::KSP, arg2::KSPType)::PetscErrorCode
end

function KSPGetType(arg1, arg2)
    @ccall libpetsc.KSPGetType(arg1::KSP, arg2::Ptr{KSPType})::PetscErrorCode
end

function KSPSetUp(arg1)
    @ccall libpetsc.KSPSetUp(arg1::KSP)::PetscErrorCode
end

function KSPSetUpOnBlocks(arg1)
    @ccall libpetsc.KSPSetUpOnBlocks(arg1::KSP)::PetscErrorCode
end

function KSPSolve(arg1, arg2, arg3)
    @ccall libpetsc.KSPSolve(arg1::KSP, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function KSPSolveTranspose(arg1, arg2, arg3)
    @ccall libpetsc.KSPSolveTranspose(arg1::KSP, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function KSPSetUseExplicitTranspose(arg1, arg2)
    @ccall libpetsc.KSPSetUseExplicitTranspose(arg1::KSP, arg2::PetscBool)::PetscErrorCode
end

function KSPMatSolve(arg1, arg2, arg3)
    @ccall libpetsc.KSPMatSolve(arg1::KSP, arg2::Mat, arg3::Mat)::PetscErrorCode
end

function KSPSetMatSolveBatchSize(arg1, arg2)
    @ccall libpetsc.KSPSetMatSolveBatchSize(arg1::KSP, arg2::PetscInt)::PetscErrorCode
end

function KSPSetMatSolveBlockSize(ksp, n)
    @ccall libpetsc.KSPSetMatSolveBlockSize(ksp::KSP, n::PetscInt)::PetscErrorCode
end

function KSPGetMatSolveBatchSize(arg1, arg2)
    @ccall libpetsc.KSPGetMatSolveBatchSize(arg1::KSP, arg2::Ptr{PetscInt})::PetscErrorCode
end

function KSPGetMatSolveBlockSize(ksp, n)
    @ccall libpetsc.KSPGetMatSolveBlockSize(ksp::KSP, n::Ptr{PetscInt})::PetscErrorCode
end

function KSPReset(arg1)
    @ccall libpetsc.KSPReset(arg1::KSP)::PetscErrorCode
end

function KSPResetViewers(arg1)
    @ccall libpetsc.KSPResetViewers(arg1::KSP)::PetscErrorCode
end

function KSPDestroy(arg1)
    @ccall libpetsc.KSPDestroy(arg1::Ptr{KSP})::PetscErrorCode
end

function KSPSetReusePreconditioner(arg1, arg2)
    @ccall libpetsc.KSPSetReusePreconditioner(arg1::KSP, arg2::PetscBool)::PetscErrorCode
end

function KSPGetReusePreconditioner(arg1, arg2)
    @ccall libpetsc.KSPGetReusePreconditioner(arg1::KSP, arg2::Ptr{PetscBool})::PetscErrorCode
end

function KSPSetSkipPCSetFromOptions(arg1, arg2)
    @ccall libpetsc.KSPSetSkipPCSetFromOptions(arg1::KSP, arg2::PetscBool)::PetscErrorCode
end

function KSPCheckSolve(arg1, arg2, arg3)
    @ccall libpetsc.KSPCheckSolve(arg1::KSP, arg2::PC, arg3::Vec)::PetscErrorCode
end

function KSPRegister(arg1, arg2)
    @ccall libpetsc.KSPRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function KSPMonitorRegister(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.KSPMonitorRegister(arg1::Ptr{Cchar}, arg2::PetscViewerType, arg3::PetscViewerFormat, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid}, arg6::Ptr{Cvoid})::PetscErrorCode
end

function KSPSetPCSide(arg1, arg2)
    @ccall libpetsc.KSPSetPCSide(arg1::KSP, arg2::PCSide)::PetscErrorCode
end

function KSPGetPCSide(arg1, arg2)
    @ccall libpetsc.KSPGetPCSide(arg1::KSP, arg2::Ptr{PCSide})::PetscErrorCode
end

function KSPSetTolerances(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.KSPSetTolerances(arg1::KSP, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscInt)::PetscErrorCode
end

function KSPGetTolerances(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.KSPGetTolerances(arg1::KSP, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscInt})::PetscErrorCode
end

function KSPSetInitialGuessNonzero(arg1, arg2)
    @ccall libpetsc.KSPSetInitialGuessNonzero(arg1::KSP, arg2::PetscBool)::PetscErrorCode
end

function KSPGetInitialGuessNonzero(arg1, arg2)
    @ccall libpetsc.KSPGetInitialGuessNonzero(arg1::KSP, arg2::Ptr{PetscBool})::PetscErrorCode
end

function KSPSetErrorIfNotConverged(arg1, arg2)
    @ccall libpetsc.KSPSetErrorIfNotConverged(arg1::KSP, arg2::PetscBool)::PetscErrorCode
end

function KSPGetErrorIfNotConverged(arg1, arg2)
    @ccall libpetsc.KSPGetErrorIfNotConverged(arg1::KSP, arg2::Ptr{PetscBool})::PetscErrorCode
end

function KSPSetComputeEigenvalues(arg1, arg2)
    @ccall libpetsc.KSPSetComputeEigenvalues(arg1::KSP, arg2::PetscBool)::PetscErrorCode
end

function KSPSetComputeRitz(arg1, arg2)
    @ccall libpetsc.KSPSetComputeRitz(arg1::KSP, arg2::PetscBool)::PetscErrorCode
end

function KSPGetComputeEigenvalues(arg1, arg2)
    @ccall libpetsc.KSPGetComputeEigenvalues(arg1::KSP, arg2::Ptr{PetscBool})::PetscErrorCode
end

function KSPSetComputeSingularValues(arg1, arg2)
    @ccall libpetsc.KSPSetComputeSingularValues(arg1::KSP, arg2::PetscBool)::PetscErrorCode
end

function KSPGetComputeSingularValues(arg1, arg2)
    @ccall libpetsc.KSPGetComputeSingularValues(arg1::KSP, arg2::Ptr{PetscBool})::PetscErrorCode
end

function KSPGetRhs(arg1, arg2)
    @ccall libpetsc.KSPGetRhs(arg1::KSP, arg2::Ptr{Vec})::PetscErrorCode
end

function KSPGetSolution(arg1, arg2)
    @ccall libpetsc.KSPGetSolution(arg1::KSP, arg2::Ptr{Vec})::PetscErrorCode
end

function KSPGetResidualNorm(arg1, arg2)
    @ccall libpetsc.KSPGetResidualNorm(arg1::KSP, arg2::Ptr{PetscReal})::PetscErrorCode
end

function KSPGetIterationNumber(arg1, arg2)
    @ccall libpetsc.KSPGetIterationNumber(arg1::KSP, arg2::Ptr{PetscInt})::PetscErrorCode
end

function KSPGetTotalIterations(arg1, arg2)
    @ccall libpetsc.KSPGetTotalIterations(arg1::KSP, arg2::Ptr{PetscInt})::PetscErrorCode
end

function KSPCreateVecs(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.KSPCreateVecs(arg1::KSP, arg2::PetscInt, arg3::Ptr{Ptr{Vec}}, arg4::PetscInt, arg5::Ptr{Ptr{Vec}})::PetscErrorCode
end

function KSPGetVecs(ksp, n, x, m, y)
    @ccall libpetsc.KSPGetVecs(ksp::KSP, n::PetscInt, x::Ptr{Ptr{Vec}}, m::PetscInt, y::Ptr{Ptr{Vec}})::PetscErrorCode
end

function KSPSetPreSolve(arg1, arg2, arg3)
    @ccall libpetsc.KSPSetPreSolve(arg1::KSP, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function KSPSetPostSolve(arg1, arg2, arg3)
    @ccall libpetsc.KSPSetPostSolve(arg1::KSP, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function KSPSetPC(arg1, arg2)
    @ccall libpetsc.KSPSetPC(arg1::KSP, arg2::PC)::PetscErrorCode
end

function KSPGetPC(arg1, arg2)
    @ccall libpetsc.KSPGetPC(arg1::KSP, arg2::Ptr{PC})::PetscErrorCode
end

function KSPMonitor(arg1, arg2, arg3)
    @ccall libpetsc.KSPMonitor(arg1::KSP, arg2::PetscInt, arg3::PetscReal)::PetscErrorCode
end

function KSPMonitorSet(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorSet(arg1::KSP, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function KSPMonitorCancel(arg1)
    @ccall libpetsc.KSPMonitorCancel(arg1::KSP)::PetscErrorCode
end

function KSPGetMonitorContext(arg1, arg2)
    @ccall libpetsc.KSPGetMonitorContext(arg1::KSP, arg2::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function KSPGetResidualHistory(arg1, arg2, arg3)
    @ccall libpetsc.KSPGetResidualHistory(arg1::KSP, arg2::Ptr{Ptr{PetscReal}}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function KSPSetResidualHistory(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPSetResidualHistory(arg1::KSP, arg2::Ptr{PetscReal}, arg3::PetscInt, arg4::PetscBool)::PetscErrorCode
end

function KSPGetErrorHistory(arg1, arg2, arg3)
    @ccall libpetsc.KSPGetErrorHistory(arg1::KSP, arg2::Ptr{Ptr{PetscReal}}, arg3::Ptr{PetscInt})::PetscErrorCode
end

function KSPSetErrorHistory(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPSetErrorHistory(arg1::KSP, arg2::Ptr{PetscReal}, arg3::PetscInt, arg4::PetscBool)::PetscErrorCode
end

function KSPBuildSolutionDefault(arg1, arg2, arg3)
    @ccall libpetsc.KSPBuildSolutionDefault(arg1::KSP, arg2::Vec, arg3::Ptr{Vec})::PetscErrorCode
end

function KSPBuildResidualDefault(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPBuildResidualDefault(arg1::KSP, arg2::Vec, arg3::Vec, arg4::Ptr{Vec})::PetscErrorCode
end

function KSPDestroyDefault(arg1)
    @ccall libpetsc.KSPDestroyDefault(arg1::KSP)::PetscErrorCode
end

function KSPSetWorkVecs(arg1, arg2)
    @ccall libpetsc.KSPSetWorkVecs(arg1::KSP, arg2::PetscInt)::PetscErrorCode
end

function PCKSPGetKSP(arg1, arg2)
    @ccall libpetsc.PCKSPGetKSP(arg1::PC, arg2::Ptr{KSP})::PetscErrorCode
end

function PCKSPSetKSP(arg1, arg2)
    @ccall libpetsc.PCKSPSetKSP(arg1::PC, arg2::KSP)::PetscErrorCode
end

function PCBJacobiGetSubKSP(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PCBJacobiGetSubKSP(arg1::PC, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{KSP}})::PetscErrorCode
end

function PCASMGetSubKSP(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PCASMGetSubKSP(arg1::PC, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{KSP}})::PetscErrorCode
end

function PCGASMGetSubKSP(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PCGASMGetSubKSP(arg1::PC, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{Ptr{KSP}})::PetscErrorCode
end

function PCFieldSplitGetSubKSP(arg1, arg2, arg3)
    @ccall libpetsc.PCFieldSplitGetSubKSP(arg1::PC, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{KSP}})::PetscErrorCode
end

function PCFieldSplitSchurGetSubKSP(arg1, arg2, arg3)
    @ccall libpetsc.PCFieldSplitSchurGetSubKSP(arg1::PC, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{KSP}})::PetscErrorCode
end

function PCMGGetSmoother(arg1, arg2, arg3)
    @ccall libpetsc.PCMGGetSmoother(arg1::PC, arg2::PetscInt, arg3::Ptr{KSP})::PetscErrorCode
end

function PCMGGetSmootherDown(arg1, arg2, arg3)
    @ccall libpetsc.PCMGGetSmootherDown(arg1::PC, arg2::PetscInt, arg3::Ptr{KSP})::PetscErrorCode
end

function PCMGGetSmootherUp(arg1, arg2, arg3)
    @ccall libpetsc.PCMGGetSmootherUp(arg1::PC, arg2::PetscInt, arg3::Ptr{KSP})::PetscErrorCode
end

function PCMGGetCoarseSolve(arg1, arg2)
    @ccall libpetsc.PCMGGetCoarseSolve(arg1::PC, arg2::Ptr{KSP})::PetscErrorCode
end

function PCGalerkinGetKSP(arg1, arg2)
    @ccall libpetsc.PCGalerkinGetKSP(arg1::PC, arg2::Ptr{KSP})::PetscErrorCode
end

function PCDeflationGetCoarseKSP(arg1, arg2)
    @ccall libpetsc.PCDeflationGetCoarseKSP(arg1::PC, arg2::Ptr{KSP})::PetscErrorCode
end

function PCMGRegisterCoarseSpaceConstructor(arg1, arg2)
    @ccall libpetsc.PCMGRegisterCoarseSpaceConstructor(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PCMGGetCoarseSpaceConstructor(arg1, arg2)
    @ccall libpetsc.PCMGGetCoarseSpaceConstructor(arg1::Ptr{Cchar}, arg2::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function KSPBuildSolution(arg1, arg2, arg3)
    @ccall libpetsc.KSPBuildSolution(arg1::KSP, arg2::Vec, arg3::Ptr{Vec})::PetscErrorCode
end

function KSPBuildResidual(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPBuildResidual(arg1::KSP, arg2::Vec, arg3::Vec, arg4::Ptr{Vec})::PetscErrorCode
end

function KSPRichardsonSetScale(arg1, arg2)
    @ccall libpetsc.KSPRichardsonSetScale(arg1::KSP, arg2::PetscReal)::PetscErrorCode
end

function KSPRichardsonSetSelfScale(arg1, arg2)
    @ccall libpetsc.KSPRichardsonSetSelfScale(arg1::KSP, arg2::PetscBool)::PetscErrorCode
end

function KSPChebyshevSetEigenvalues(arg1, arg2, arg3)
    @ccall libpetsc.KSPChebyshevSetEigenvalues(arg1::KSP, arg2::PetscReal, arg3::PetscReal)::PetscErrorCode
end

function KSPChebyshevEstEigSet(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.KSPChebyshevEstEigSet(arg1::KSP, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal)::PetscErrorCode
end

function KSPChebyshevEstEigSetUseNoisy(arg1, arg2)
    @ccall libpetsc.KSPChebyshevEstEigSetUseNoisy(arg1::KSP, arg2::PetscBool)::PetscErrorCode
end

function KSPChebyshevEstEigGetKSP(arg1, arg2)
    @ccall libpetsc.KSPChebyshevEstEigGetKSP(arg1::KSP, arg2::Ptr{KSP})::PetscErrorCode
end

function KSPComputeExtremeSingularValues(arg1, arg2, arg3)
    @ccall libpetsc.KSPComputeExtremeSingularValues(arg1::KSP, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal})::PetscErrorCode
end

function KSPComputeEigenvalues(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.KSPComputeEigenvalues(arg1::KSP, arg2::PetscInt, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscInt})::PetscErrorCode
end

function KSPComputeEigenvaluesExplicitly(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPComputeEigenvaluesExplicitly(arg1::KSP, arg2::PetscInt, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal})::PetscErrorCode
end

function KSPComputeRitz(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.KSPComputeRitz(arg1::KSP, arg2::PetscBool, arg3::PetscBool, arg4::Ptr{PetscInt}, arg5::Ptr{Vec}, arg6::Ptr{PetscReal}, arg7::Ptr{PetscReal})::PetscErrorCode
end

@enum KSPFCDTruncationType::UInt32 begin
    KSP_FCD_TRUNC_TYPE_STANDARD = 0
    KSP_FCD_TRUNC_TYPE_NOTAY = 1
end

function KSPFCGSetMmax(arg1, arg2)
    @ccall libpetsc.KSPFCGSetMmax(arg1::KSP, arg2::PetscInt)::PetscErrorCode
end

function KSPFCGGetMmax(arg1, arg2)
    @ccall libpetsc.KSPFCGGetMmax(arg1::KSP, arg2::Ptr{PetscInt})::PetscErrorCode
end

function KSPFCGSetNprealloc(arg1, arg2)
    @ccall libpetsc.KSPFCGSetNprealloc(arg1::KSP, arg2::PetscInt)::PetscErrorCode
end

function KSPFCGGetNprealloc(arg1, arg2)
    @ccall libpetsc.KSPFCGGetNprealloc(arg1::KSP, arg2::Ptr{PetscInt})::PetscErrorCode
end

function KSPFCGSetTruncationType(arg1, arg2)
    @ccall libpetsc.KSPFCGSetTruncationType(arg1::KSP, arg2::KSPFCDTruncationType)::PetscErrorCode
end

function KSPFCGGetTruncationType(arg1, arg2)
    @ccall libpetsc.KSPFCGGetTruncationType(arg1::KSP, arg2::Ptr{KSPFCDTruncationType})::PetscErrorCode
end

function KSPPIPEFCGSetMmax(arg1, arg2)
    @ccall libpetsc.KSPPIPEFCGSetMmax(arg1::KSP, arg2::PetscInt)::PetscErrorCode
end

function KSPPIPEFCGGetMmax(arg1, arg2)
    @ccall libpetsc.KSPPIPEFCGGetMmax(arg1::KSP, arg2::Ptr{PetscInt})::PetscErrorCode
end

function KSPPIPEFCGSetNprealloc(arg1, arg2)
    @ccall libpetsc.KSPPIPEFCGSetNprealloc(arg1::KSP, arg2::PetscInt)::PetscErrorCode
end

function KSPPIPEFCGGetNprealloc(arg1, arg2)
    @ccall libpetsc.KSPPIPEFCGGetNprealloc(arg1::KSP, arg2::Ptr{PetscInt})::PetscErrorCode
end

function KSPPIPEFCGSetTruncationType(arg1, arg2)
    @ccall libpetsc.KSPPIPEFCGSetTruncationType(arg1::KSP, arg2::KSPFCDTruncationType)::PetscErrorCode
end

function KSPPIPEFCGGetTruncationType(arg1, arg2)
    @ccall libpetsc.KSPPIPEFCGGetTruncationType(arg1::KSP, arg2::Ptr{KSPFCDTruncationType})::PetscErrorCode
end

function KSPPIPEGCRSetMmax(arg1, arg2)
    @ccall libpetsc.KSPPIPEGCRSetMmax(arg1::KSP, arg2::PetscInt)::PetscErrorCode
end

function KSPPIPEGCRGetMmax(arg1, arg2)
    @ccall libpetsc.KSPPIPEGCRGetMmax(arg1::KSP, arg2::Ptr{PetscInt})::PetscErrorCode
end

function KSPPIPEGCRSetNprealloc(arg1, arg2)
    @ccall libpetsc.KSPPIPEGCRSetNprealloc(arg1::KSP, arg2::PetscInt)::PetscErrorCode
end

function KSPPIPEGCRGetNprealloc(arg1, arg2)
    @ccall libpetsc.KSPPIPEGCRGetNprealloc(arg1::KSP, arg2::Ptr{PetscInt})::PetscErrorCode
end

function KSPPIPEGCRSetTruncationType(arg1, arg2)
    @ccall libpetsc.KSPPIPEGCRSetTruncationType(arg1::KSP, arg2::KSPFCDTruncationType)::PetscErrorCode
end

function KSPPIPEGCRGetTruncationType(arg1, arg2)
    @ccall libpetsc.KSPPIPEGCRGetTruncationType(arg1::KSP, arg2::Ptr{KSPFCDTruncationType})::PetscErrorCode
end

function KSPPIPEGCRSetUnrollW(arg1, arg2)
    @ccall libpetsc.KSPPIPEGCRSetUnrollW(arg1::KSP, arg2::PetscBool)::PetscErrorCode
end

function KSPPIPEGCRGetUnrollW(arg1, arg2)
    @ccall libpetsc.KSPPIPEGCRGetUnrollW(arg1::KSP, arg2::Ptr{PetscBool})::PetscErrorCode
end

function KSPGMRESSetRestart(arg1, arg2)
    @ccall libpetsc.KSPGMRESSetRestart(arg1::KSP, arg2::PetscInt)::PetscErrorCode
end

function KSPGMRESGetRestart(arg1, arg2)
    @ccall libpetsc.KSPGMRESGetRestart(arg1::KSP, arg2::Ptr{PetscInt})::PetscErrorCode
end

function KSPGMRESSetHapTol(arg1, arg2)
    @ccall libpetsc.KSPGMRESSetHapTol(arg1::KSP, arg2::PetscReal)::PetscErrorCode
end

function KSPGMRESSetBreakdownTolerance(arg1, arg2)
    @ccall libpetsc.KSPGMRESSetBreakdownTolerance(arg1::KSP, arg2::PetscReal)::PetscErrorCode
end

function KSPGMRESSetPreAllocateVectors(arg1)
    @ccall libpetsc.KSPGMRESSetPreAllocateVectors(arg1::KSP)::PetscErrorCode
end

function KSPGMRESSetOrthogonalization(arg1, arg2)
    @ccall libpetsc.KSPGMRESSetOrthogonalization(arg1::KSP, arg2::Ptr{Cvoid})::PetscErrorCode
end

function KSPGMRESGetOrthogonalization(arg1, arg2)
    @ccall libpetsc.KSPGMRESGetOrthogonalization(arg1::KSP, arg2::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function KSPGMRESModifiedGramSchmidtOrthogonalization(arg1, arg2)
    @ccall libpetsc.KSPGMRESModifiedGramSchmidtOrthogonalization(arg1::KSP, arg2::PetscInt)::PetscErrorCode
end

function KSPGMRESClassicalGramSchmidtOrthogonalization(arg1, arg2)
    @ccall libpetsc.KSPGMRESClassicalGramSchmidtOrthogonalization(arg1::KSP, arg2::PetscInt)::PetscErrorCode
end

function KSPLGMRESSetAugDim(arg1, arg2)
    @ccall libpetsc.KSPLGMRESSetAugDim(arg1::KSP, arg2::PetscInt)::PetscErrorCode
end

function KSPLGMRESSetConstant(arg1)
    @ccall libpetsc.KSPLGMRESSetConstant(arg1::KSP)::PetscErrorCode
end

function KSPPIPEFGMRESSetShift(arg1, arg2)
    @ccall libpetsc.KSPPIPEFGMRESSetShift(arg1::KSP, arg2::PetscScalar)::PetscErrorCode
end

function KSPGCRSetRestart(arg1, arg2)
    @ccall libpetsc.KSPGCRSetRestart(arg1::KSP, arg2::PetscInt)::PetscErrorCode
end

function KSPGCRGetRestart(arg1, arg2)
    @ccall libpetsc.KSPGCRGetRestart(arg1::KSP, arg2::Ptr{PetscInt})::PetscErrorCode
end

function KSPGCRSetModifyPC(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPGCRSetModifyPC(arg1::KSP, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function KSPFETIDPGetInnerBDDC(arg1, arg2)
    @ccall libpetsc.KSPFETIDPGetInnerBDDC(arg1::KSP, arg2::Ptr{PC})::PetscErrorCode
end

function KSPFETIDPSetInnerBDDC(arg1, arg2)
    @ccall libpetsc.KSPFETIDPSetInnerBDDC(arg1::KSP, arg2::PC)::PetscErrorCode
end

function KSPFETIDPGetInnerKSP(arg1, arg2)
    @ccall libpetsc.KSPFETIDPGetInnerKSP(arg1::KSP, arg2::Ptr{KSP})::PetscErrorCode
end

function KSPFETIDPSetPressureOperator(arg1, arg2)
    @ccall libpetsc.KSPFETIDPSetPressureOperator(arg1::KSP, arg2::Mat)::PetscErrorCode
end

function KSPHPDDMSetDeflationSpace(arg1, arg2)
    @ccall libpetsc.KSPHPDDMSetDeflationSpace(arg1::KSP, arg2::Mat)::PetscErrorCode
end

function KSPHPDDMGetDeflationSpace(arg1, arg2)
    @ccall libpetsc.KSPHPDDMGetDeflationSpace(arg1::KSP, arg2::Ptr{Mat})::PetscErrorCode
end

function KSPHPDDMMatSolve(ksp, B, X)
    @ccall libpetsc.KSPHPDDMMatSolve(ksp::KSP, B::Mat, X::Mat)::PetscErrorCode
end

@enum KSPHPDDMType::UInt32 begin
    KSP_HPDDM_TYPE_GMRES = 0
    KSP_HPDDM_TYPE_BGMRES = 1
    KSP_HPDDM_TYPE_CG = 2
    KSP_HPDDM_TYPE_BCG = 3
    KSP_HPDDM_TYPE_GCRODR = 4
    KSP_HPDDM_TYPE_BGCRODR = 5
    KSP_HPDDM_TYPE_BFBCG = 6
    KSP_HPDDM_TYPE_PREONLY = 7
end

function KSPHPDDMSetType(arg1, arg2)
    @ccall libpetsc.KSPHPDDMSetType(arg1::KSP, arg2::KSPHPDDMType)::PetscErrorCode
end

function KSPHPDDMGetType(arg1, arg2)
    @ccall libpetsc.KSPHPDDMGetType(arg1::KSP, arg2::Ptr{KSPHPDDMType})::PetscErrorCode
end

@enum KSPGMRESCGSRefinementType::UInt32 begin
    KSP_GMRES_CGS_REFINE_NEVER = 0
    KSP_GMRES_CGS_REFINE_IFNEEDED = 1
    KSP_GMRES_CGS_REFINE_ALWAYS = 2
end

function KSPGMRESSetCGSRefinementType(arg1, arg2)
    @ccall libpetsc.KSPGMRESSetCGSRefinementType(arg1::KSP, arg2::KSPGMRESCGSRefinementType)::PetscErrorCode
end

function KSPGMRESGetCGSRefinementType(arg1, arg2)
    @ccall libpetsc.KSPGMRESGetCGSRefinementType(arg1::KSP, arg2::Ptr{KSPGMRESCGSRefinementType})::PetscErrorCode
end

function KSPFGMRESModifyPCNoChange(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.KSPFGMRESModifyPCNoChange(arg1::KSP, arg2::PetscInt, arg3::PetscInt, arg4::PetscReal, arg5::Ptr{Cvoid})::PetscErrorCode
end

function KSPFGMRESModifyPCKSP(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.KSPFGMRESModifyPCKSP(arg1::KSP, arg2::PetscInt, arg3::PetscInt, arg4::PetscReal, arg5::Ptr{Cvoid})::PetscErrorCode
end

function KSPFGMRESSetModifyPC(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPFGMRESSetModifyPC(arg1::KSP, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function KSPQCGSetTrustRegionRadius(arg1, arg2)
    @ccall libpetsc.KSPQCGSetTrustRegionRadius(arg1::KSP, arg2::PetscReal)::PetscErrorCode
end

function KSPQCGGetQuadratic(arg1, arg2)
    @ccall libpetsc.KSPQCGGetQuadratic(arg1::KSP, arg2::Ptr{PetscReal})::PetscErrorCode
end

function KSPQCGGetTrialStepNorm(arg1, arg2)
    @ccall libpetsc.KSPQCGGetTrialStepNorm(arg1::KSP, arg2::Ptr{PetscReal})::PetscErrorCode
end

function KSPBCGSLSetXRes(arg1, arg2)
    @ccall libpetsc.KSPBCGSLSetXRes(arg1::KSP, arg2::PetscReal)::PetscErrorCode
end

function KSPBCGSLSetPol(arg1, arg2)
    @ccall libpetsc.KSPBCGSLSetPol(arg1::KSP, arg2::PetscBool)::PetscErrorCode
end

function KSPBCGSLSetEll(arg1, arg2)
    @ccall libpetsc.KSPBCGSLSetEll(arg1::KSP, arg2::PetscInt)::PetscErrorCode
end

function KSPBCGSLSetUsePseudoinverse(arg1, arg2)
    @ccall libpetsc.KSPBCGSLSetUsePseudoinverse(arg1::KSP, arg2::PetscBool)::PetscErrorCode
end

function KSPSetFromOptions(arg1)
    @ccall libpetsc.KSPSetFromOptions(arg1::KSP)::PetscErrorCode
end

function KSPResetFromOptions(arg1)
    @ccall libpetsc.KSPResetFromOptions(arg1::KSP)::PetscErrorCode
end

function KSPAddOptionsChecker(arg1)
    @ccall libpetsc.KSPAddOptionsChecker(arg1::Ptr{Cvoid})::PetscErrorCode
end

function KSPMonitorSetFromOptions(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorSetFromOptions(arg1::KSP, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function KSPMonitorLGCreate(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11)
    @ccall libpetsc.KSPMonitorLGCreate(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::PetscInt, arg6::Ptr{Ptr{Cchar}}, arg7::Cint, arg8::Cint, arg9::Cint, arg10::Cint, arg11::Ptr{PetscDrawLG})::PetscErrorCode
end

function KSPMonitorResidual(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorResidual(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPMonitorResidualDraw(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorResidualDraw(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPMonitorResidualDrawLG(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorResidualDrawLG(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPMonitorResidualDrawLGCreate(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorResidualDrawLGCreate(arg1::PetscViewer, arg2::PetscViewerFormat, arg3::Ptr{Cvoid}, arg4::Ptr{Ptr{PetscViewerAndFormat}})::PetscErrorCode
end

function KSPMonitorResidualShort(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorResidualShort(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPMonitorResidualRange(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorResidualRange(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPMonitorTrueResidual(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorTrueResidual(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPMonitorTrueResidualDraw(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorTrueResidualDraw(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPMonitorTrueResidualDrawLG(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorTrueResidualDrawLG(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPMonitorTrueResidualDrawLGCreate(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorTrueResidualDrawLGCreate(arg1::PetscViewer, arg2::PetscViewerFormat, arg3::Ptr{Cvoid}, arg4::Ptr{Ptr{PetscViewerAndFormat}})::PetscErrorCode
end

function KSPMonitorTrueResidualMax(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorTrueResidualMax(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPMonitorError(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorError(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPMonitorErrorDraw(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorErrorDraw(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPMonitorErrorDrawLG(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorErrorDrawLG(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPMonitorErrorDrawLGCreate(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorErrorDrawLGCreate(arg1::PetscViewer, arg2::PetscViewerFormat, arg3::Ptr{Cvoid}, arg4::Ptr{Ptr{PetscViewerAndFormat}})::PetscErrorCode
end

function KSPMonitorSolution(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorSolution(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPMonitorSolutionDraw(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorSolutionDraw(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPMonitorSolutionDrawLG(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorSolutionDrawLG(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPMonitorSolutionDrawLGCreate(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorSolutionDrawLGCreate(arg1::PetscViewer, arg2::PetscViewerFormat, arg3::Ptr{Cvoid}, arg4::Ptr{Ptr{PetscViewerAndFormat}})::PetscErrorCode
end

function KSPMonitorSingularValue(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorSingularValue(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPMonitorSingularValueCreate(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorSingularValueCreate(arg1::PetscViewer, arg2::PetscViewerFormat, arg3::Ptr{Cvoid}, arg4::Ptr{Ptr{PetscViewerAndFormat}})::PetscErrorCode
end

function KSPMonitorDefault(ksp, n, rnorm, vf)
    @ccall libpetsc.KSPMonitorDefault(ksp::KSP, n::PetscInt, rnorm::PetscReal, vf::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPMonitorTrueResidualNorm(ksp, n, rnorm, vf)
    @ccall libpetsc.KSPMonitorTrueResidualNorm(ksp::KSP, n::PetscInt, rnorm::PetscReal, vf::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPMonitorTrueResidualMaxNorm(ksp, n, rnorm, vf)
    @ccall libpetsc.KSPMonitorTrueResidualMaxNorm(ksp::KSP, n::PetscInt, rnorm::PetscReal, vf::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPGMRESMonitorKrylov(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPGMRESMonitorKrylov(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{Cvoid})::PetscErrorCode
end

function KSPMonitorDynamicTolerance(ksp, its, fnorm, dummy)
    @ccall libpetsc.KSPMonitorDynamicTolerance(ksp::KSP, its::PetscInt, fnorm::PetscReal, dummy::Ptr{Cvoid})::PetscErrorCode
end

function KSPMonitorDynamicToleranceDestroy(dummy)
    @ccall libpetsc.KSPMonitorDynamicToleranceDestroy(dummy::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function KSPMonitorSAWs(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorSAWs(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{Cvoid})::PetscErrorCode
end

function KSPMonitorSAWsCreate(arg1, arg2)
    @ccall libpetsc.KSPMonitorSAWsCreate(arg1::KSP, arg2::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function KSPMonitorSAWsDestroy(arg1)
    @ccall libpetsc.KSPMonitorSAWsDestroy(arg1::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function KSPUnwindPreconditioner(arg1, arg2, arg3)
    @ccall libpetsc.KSPUnwindPreconditioner(arg1::KSP, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function KSPInitialResidual(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.KSPInitialResidual(arg1::KSP, arg2::Vec, arg3::Vec, arg4::Vec, arg5::Vec, arg6::Vec)::PetscErrorCode
end

function KSPSetOperators(arg1, arg2, arg3)
    @ccall libpetsc.KSPSetOperators(arg1::KSP, arg2::Mat, arg3::Mat)::PetscErrorCode
end

function KSPGetOperators(arg1, arg2, arg3)
    @ccall libpetsc.KSPGetOperators(arg1::KSP, arg2::Ptr{Mat}, arg3::Ptr{Mat})::PetscErrorCode
end

function KSPGetOperatorsSet(arg1, arg2, arg3)
    @ccall libpetsc.KSPGetOperatorsSet(arg1::KSP, arg2::Ptr{PetscBool}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function KSPSetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.KSPSetOptionsPrefix(arg1::KSP, arg2::Ptr{Cchar})::PetscErrorCode
end

function KSPAppendOptionsPrefix(arg1, arg2)
    @ccall libpetsc.KSPAppendOptionsPrefix(arg1::KSP, arg2::Ptr{Cchar})::PetscErrorCode
end

function KSPGetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.KSPGetOptionsPrefix(arg1::KSP, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function KSPSetDiagonalScale(arg1, arg2)
    @ccall libpetsc.KSPSetDiagonalScale(arg1::KSP, arg2::PetscBool)::PetscErrorCode
end

function KSPGetDiagonalScale(arg1, arg2)
    @ccall libpetsc.KSPGetDiagonalScale(arg1::KSP, arg2::Ptr{PetscBool})::PetscErrorCode
end

function KSPSetDiagonalScaleFix(arg1, arg2)
    @ccall libpetsc.KSPSetDiagonalScaleFix(arg1::KSP, arg2::PetscBool)::PetscErrorCode
end

function KSPGetDiagonalScaleFix(arg1, arg2)
    @ccall libpetsc.KSPGetDiagonalScaleFix(arg1::KSP, arg2::Ptr{PetscBool})::PetscErrorCode
end

function KSPView(arg1, arg2)
    @ccall libpetsc.KSPView(arg1::KSP, arg2::PetscViewer)::PetscErrorCode
end

function KSPLoad(arg1, arg2)
    @ccall libpetsc.KSPLoad(arg1::KSP, arg2::PetscViewer)::PetscErrorCode
end

function KSPViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.KSPViewFromOptions(arg1::KSP, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function KSPConvergedReasonView(arg1, arg2)
    @ccall libpetsc.KSPConvergedReasonView(arg1::KSP, arg2::PetscViewer)::PetscErrorCode
end

function KSPConvergedReasonViewSet(arg1, arg2, vctx, arg4)
    @ccall libpetsc.KSPConvergedReasonViewSet(arg1::KSP, arg2::Ptr{Cvoid}, vctx::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function KSPConvergedReasonViewFromOptions(arg1)
    @ccall libpetsc.KSPConvergedReasonViewFromOptions(arg1::KSP)::PetscErrorCode
end

function KSPConvergedReasonViewCancel(arg1)
    @ccall libpetsc.KSPConvergedReasonViewCancel(arg1::KSP)::PetscErrorCode
end

function KSPConvergedRateView(arg1, arg2)
    @ccall libpetsc.KSPConvergedRateView(arg1::KSP, arg2::PetscViewer)::PetscErrorCode
end

function KSPReasonView(ksp, v)
    @ccall libpetsc.KSPReasonView(ksp::KSP, v::PetscViewer)::PetscErrorCode
end

function KSPReasonViewFromOptions(ksp)
    @ccall libpetsc.KSPReasonViewFromOptions(ksp::KSP)::PetscErrorCode
end

function KSPLSQRSetExactMatNorm(arg1, arg2)
    @ccall libpetsc.KSPLSQRSetExactMatNorm(arg1::KSP, arg2::PetscBool)::PetscErrorCode
end

function KSPLSQRSetComputeStandardErrorVec(arg1, arg2)
    @ccall libpetsc.KSPLSQRSetComputeStandardErrorVec(arg1::KSP, arg2::PetscBool)::PetscErrorCode
end

function KSPLSQRGetStandardErrorVec(arg1, arg2)
    @ccall libpetsc.KSPLSQRGetStandardErrorVec(arg1::KSP, arg2::Ptr{Vec})::PetscErrorCode
end

function KSPLSQRGetNorms(arg1, arg2, arg3)
    @ccall libpetsc.KSPLSQRGetNorms(arg1::KSP, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal})::PetscErrorCode
end

function KSPLSQRMonitorResidual(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPLSQRMonitorResidual(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPLSQRMonitorResidualDrawLG(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPLSQRMonitorResidualDrawLG(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPLSQRMonitorResidualDrawLGCreate(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPLSQRMonitorResidualDrawLGCreate(arg1::PetscViewer, arg2::PetscViewerFormat, arg3::Ptr{Cvoid}, arg4::Ptr{Ptr{PetscViewerAndFormat}})::PetscErrorCode
end

function PCRedundantGetKSP(arg1, arg2)
    @ccall libpetsc.PCRedundantGetKSP(arg1::PC, arg2::Ptr{KSP})::PetscErrorCode
end

function PCRedistributeGetKSP(arg1, arg2)
    @ccall libpetsc.PCRedistributeGetKSP(arg1::PC, arg2::Ptr{KSP})::PetscErrorCode
end

function PCTelescopeGetKSP(arg1, arg2)
    @ccall libpetsc.PCTelescopeGetKSP(arg1::PC, arg2::Ptr{KSP})::PetscErrorCode
end

@enum KSPNormType::Int32 begin
    KSP_NORM_DEFAULT = -1
    KSP_NORM_NONE = 0
    KSP_NORM_PRECONDITIONED = 1
    KSP_NORM_UNPRECONDITIONED = 2
    KSP_NORM_NATURAL = 3
end

function KSPSetNormType(arg1, arg2)
    @ccall libpetsc.KSPSetNormType(arg1::KSP, arg2::KSPNormType)::PetscErrorCode
end

function KSPGetNormType(arg1, arg2)
    @ccall libpetsc.KSPGetNormType(arg1::KSP, arg2::Ptr{KSPNormType})::PetscErrorCode
end

function KSPSetSupportedNorm(ksp, arg2, arg3, arg4)
    @ccall libpetsc.KSPSetSupportedNorm(ksp::KSP, arg2::KSPNormType, arg3::PCSide, arg4::PetscInt)::PetscErrorCode
end

function KSPSetCheckNormIteration(arg1, arg2)
    @ccall libpetsc.KSPSetCheckNormIteration(arg1::KSP, arg2::PetscInt)::PetscErrorCode
end

function KSPSetLagNorm(arg1, arg2)
    @ccall libpetsc.KSPSetLagNorm(arg1::KSP, arg2::PetscBool)::PetscErrorCode
end

@enum KSPConvergedReason::Int32 begin
    KSP_CONVERGED_RTOL_NORMAL = 1
    KSP_CONVERGED_ATOL_NORMAL = 9
    KSP_CONVERGED_RTOL = 2
    KSP_CONVERGED_ATOL = 3
    KSP_CONVERGED_ITS = 4
    KSP_CONVERGED_CG_NEG_CURVE = 5
    KSP_CONVERGED_CG_CONSTRAINED = 6
    KSP_CONVERGED_STEP_LENGTH = 7
    KSP_CONVERGED_HAPPY_BREAKDOWN = 8
    KSP_DIVERGED_NULL = -2
    KSP_DIVERGED_ITS = -3
    KSP_DIVERGED_DTOL = -4
    KSP_DIVERGED_BREAKDOWN = -5
    KSP_DIVERGED_BREAKDOWN_BICG = -6
    KSP_DIVERGED_NONSYMMETRIC = -7
    KSP_DIVERGED_INDEFINITE_PC = -8
    KSP_DIVERGED_NANORINF = -9
    KSP_DIVERGED_INDEFINITE_MAT = -10
    KSP_DIVERGED_PC_FAILED = -11
    # KSP_DIVERGED_PCSETUP_FAILED = -11
    KSP_CONVERGED_ITERATING = 0
end

function KSPSetConvergenceTest(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPSetConvergenceTest(arg1::KSP, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function KSPGetConvergenceTest(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPGetConvergenceTest(arg1::KSP, arg2::Ptr{Ptr{Cvoid}}, arg3::Ptr{Ptr{Cvoid}}, arg4::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function KSPGetAndClearConvergenceTest(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPGetAndClearConvergenceTest(arg1::KSP, arg2::Ptr{Ptr{Cvoid}}, arg3::Ptr{Ptr{Cvoid}}, arg4::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function KSPGetConvergenceContext(arg1, arg2)
    @ccall libpetsc.KSPGetConvergenceContext(arg1::KSP, arg2::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function KSPConvergedDefault(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.KSPConvergedDefault(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{KSPConvergedReason}, arg5::Ptr{Cvoid})::PetscErrorCode
end

function KSPLSQRConvergedDefault(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.KSPLSQRConvergedDefault(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{KSPConvergedReason}, arg5::Ptr{Cvoid})::PetscErrorCode
end

function KSPConvergedDefaultDestroy(arg1)
    @ccall libpetsc.KSPConvergedDefaultDestroy(arg1::Ptr{Cvoid})::PetscErrorCode
end

function KSPConvergedDefaultCreate(arg1)
    @ccall libpetsc.KSPConvergedDefaultCreate(arg1::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function KSPConvergedDefaultSetUIRNorm(arg1)
    @ccall libpetsc.KSPConvergedDefaultSetUIRNorm(arg1::KSP)::PetscErrorCode
end

function KSPConvergedDefaultSetUMIRNorm(arg1)
    @ccall libpetsc.KSPConvergedDefaultSetUMIRNorm(arg1::KSP)::PetscErrorCode
end

function KSPConvergedDefaultSetConvergedMaxits(arg1, arg2)
    @ccall libpetsc.KSPConvergedDefaultSetConvergedMaxits(arg1::KSP, arg2::PetscBool)::PetscErrorCode
end

function KSPConvergedSkip(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.KSPConvergedSkip(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{KSPConvergedReason}, arg5::Ptr{Cvoid})::PetscErrorCode
end

function KSPGetConvergedReason(arg1, arg2)
    @ccall libpetsc.KSPGetConvergedReason(arg1::KSP, arg2::Ptr{KSPConvergedReason})::PetscErrorCode
end

function KSPGetConvergedReasonString(arg1, arg2)
    @ccall libpetsc.KSPGetConvergedReasonString(arg1::KSP, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function KSPComputeConvergenceRate(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.KSPComputeConvergenceRate(arg1::KSP, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal})::PetscErrorCode
end

function KSPComputeOperator(arg1, arg2, arg3)
    @ccall libpetsc.KSPComputeOperator(arg1::KSP, arg2::MatType, arg3::Ptr{Mat})::PetscErrorCode
end

function KSPComputeExplicitOperator(A, B)
    @ccall libpetsc.KSPComputeExplicitOperator(A::KSP, B::Ptr{Mat})::PetscErrorCode
end

@enum KSPCGType::UInt32 begin
    KSP_CG_SYMMETRIC = 0
    KSP_CG_HERMITIAN = 1
end

function KSPCGSetType(arg1, arg2)
    @ccall libpetsc.KSPCGSetType(arg1::KSP, arg2::KSPCGType)::PetscErrorCode
end

function KSPCGUseSingleReduction(arg1, arg2)
    @ccall libpetsc.KSPCGUseSingleReduction(arg1::KSP, arg2::PetscBool)::PetscErrorCode
end

function KSPCGSetRadius(arg1, arg2)
    @ccall libpetsc.KSPCGSetRadius(arg1::KSP, arg2::PetscReal)::PetscErrorCode
end

function KSPCGGetNormD(arg1, arg2)
    @ccall libpetsc.KSPCGGetNormD(arg1::KSP, arg2::Ptr{PetscReal})::PetscErrorCode
end

function KSPCGGetObjFcn(arg1, arg2)
    @ccall libpetsc.KSPCGGetObjFcn(arg1::KSP, arg2::Ptr{PetscReal})::PetscErrorCode
end

function KSPGLTRGetMinEig(arg1, arg2)
    @ccall libpetsc.KSPGLTRGetMinEig(arg1::KSP, arg2::Ptr{PetscReal})::PetscErrorCode
end

function KSPGLTRGetLambda(arg1, arg2)
    @ccall libpetsc.KSPGLTRGetLambda(arg1::KSP, arg2::Ptr{PetscReal})::PetscErrorCode
end

function KSPCGGLTRGetMinEig(ksp, x)
    @ccall libpetsc.KSPCGGLTRGetMinEig(ksp::KSP, x::Ptr{PetscReal})::PetscErrorCode
end

function KSPCGGLTRGetLambda(ksp, x)
    @ccall libpetsc.KSPCGGLTRGetLambda(ksp::KSP, x::Ptr{PetscReal})::PetscErrorCode
end

function KSPPythonSetType(arg1, arg2)
    @ccall libpetsc.KSPPythonSetType(arg1::KSP, arg2::Ptr{Cchar})::PetscErrorCode
end

function PCPreSolve(arg1, arg2)
    @ccall libpetsc.PCPreSolve(arg1::PC, arg2::KSP)::PetscErrorCode
end

function PCPostSolve(arg1, arg2)
    @ccall libpetsc.PCPostSolve(arg1::PC, arg2::KSP)::PetscErrorCode
end

function KSPMonitorLGRange(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorLGRange(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{Cvoid})::PetscErrorCode
end

function PCShellSetPreSolve(arg1, arg2)
    @ccall libpetsc.PCShellSetPreSolve(arg1::PC, arg2::Ptr{Cvoid})::PetscErrorCode
end

function PCShellSetPostSolve(arg1, arg2)
    @ccall libpetsc.PCShellSetPostSolve(arg1::PC, arg2::Ptr{Cvoid})::PetscErrorCode
end

mutable struct _p_KSPGuess end

const KSPGuess = Ptr{_p_KSPGuess}

const KSPGuessType = Ptr{Cchar}

function KSPGuessRegister(arg1, arg2)
    @ccall libpetsc.KSPGuessRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function KSPSetGuess(arg1, arg2)
    @ccall libpetsc.KSPSetGuess(arg1::KSP, arg2::KSPGuess)::PetscErrorCode
end

function KSPGetGuess(arg1, arg2)
    @ccall libpetsc.KSPGetGuess(arg1::KSP, arg2::Ptr{KSPGuess})::PetscErrorCode
end

function KSPGuessView(arg1, arg2)
    @ccall libpetsc.KSPGuessView(arg1::KSPGuess, arg2::PetscViewer)::PetscErrorCode
end

function KSPGuessDestroy(arg1)
    @ccall libpetsc.KSPGuessDestroy(arg1::Ptr{KSPGuess})::PetscErrorCode
end

function KSPGuessCreate(arg1, arg2)
    @ccall libpetsc.KSPGuessCreate(arg1::MPI_Comm, arg2::Ptr{KSPGuess})::PetscErrorCode
end

function KSPGuessSetType(arg1, arg2)
    @ccall libpetsc.KSPGuessSetType(arg1::KSPGuess, arg2::KSPGuessType)::PetscErrorCode
end

function KSPGuessGetType(arg1, arg2)
    @ccall libpetsc.KSPGuessGetType(arg1::KSPGuess, arg2::Ptr{KSPGuessType})::PetscErrorCode
end

function KSPGuessSetUp(arg1)
    @ccall libpetsc.KSPGuessSetUp(arg1::KSPGuess)::PetscErrorCode
end

function KSPGuessUpdate(arg1, arg2, arg3)
    @ccall libpetsc.KSPGuessUpdate(arg1::KSPGuess, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function KSPGuessFormGuess(arg1, arg2, arg3)
    @ccall libpetsc.KSPGuessFormGuess(arg1::KSPGuess, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function KSPGuessSetFromOptions(arg1)
    @ccall libpetsc.KSPGuessSetFromOptions(arg1::KSPGuess)::PetscErrorCode
end

function KSPGuessFischerSetModel(arg1, arg2, arg3)
    @ccall libpetsc.KSPGuessFischerSetModel(arg1::KSPGuess, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function KSPSetUseFischerGuess(arg1, arg2, arg3)
    @ccall libpetsc.KSPSetUseFischerGuess(arg1::KSP, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function KSPSetInitialGuessKnoll(arg1, arg2)
    @ccall libpetsc.KSPSetInitialGuessKnoll(arg1::KSP, arg2::PetscBool)::PetscErrorCode
end

function KSPGetInitialGuessKnoll(arg1, arg2)
    @ccall libpetsc.KSPGetInitialGuessKnoll(arg1::KSP, arg2::Ptr{PetscBool})::PetscErrorCode
end

@enum MatSchurComplementAinvType::UInt32 begin
    MAT_SCHUR_COMPLEMENT_AINV_DIAG = 0
    MAT_SCHUR_COMPLEMENT_AINV_LUMP = 1
    MAT_SCHUR_COMPLEMENT_AINV_BLOCK_DIAG = 2
end

@enum MatLMVMSymBroydenScaleType::UInt32 begin
    MAT_LMVM_SYMBROYDEN_SCALE_NONE = 0
    MAT_LMVM_SYMBROYDEN_SCALE_SCALAR = 1
    MAT_LMVM_SYMBROYDEN_SCALE_DIAGONAL = 2
    MAT_LMVM_SYMBROYDEN_SCALE_USER = 3
end

function MatCreateSchurComplement(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatCreateSchurComplement(arg1::Mat, arg2::Mat, arg3::Mat, arg4::Mat, arg5::Mat, arg6::Ptr{Mat})::PetscErrorCode
end

function MatSchurComplementGetKSP(arg1, arg2)
    @ccall libpetsc.MatSchurComplementGetKSP(arg1::Mat, arg2::Ptr{KSP})::PetscErrorCode
end

function MatSchurComplementSetKSP(arg1, arg2)
    @ccall libpetsc.MatSchurComplementSetKSP(arg1::Mat, arg2::KSP)::PetscErrorCode
end

function MatSchurComplementSetSubMatrices(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatSchurComplementSetSubMatrices(arg1::Mat, arg2::Mat, arg3::Mat, arg4::Mat, arg5::Mat, arg6::Mat)::PetscErrorCode
end

function MatSchurComplementUpdateSubMatrices(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatSchurComplementUpdateSubMatrices(arg1::Mat, arg2::Mat, arg3::Mat, arg4::Mat, arg5::Mat, arg6::Mat)::PetscErrorCode
end

function MatSchurComplementGetSubMatrices(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.MatSchurComplementGetSubMatrices(arg1::Mat, arg2::Ptr{Mat}, arg3::Ptr{Mat}, arg4::Ptr{Mat}, arg5::Ptr{Mat}, arg6::Ptr{Mat})::PetscErrorCode
end

function MatSchurComplementSetAinvType(arg1, arg2)
    @ccall libpetsc.MatSchurComplementSetAinvType(arg1::Mat, arg2::MatSchurComplementAinvType)::PetscErrorCode
end

function MatSchurComplementGetAinvType(arg1, arg2)
    @ccall libpetsc.MatSchurComplementGetAinvType(arg1::Mat, arg2::Ptr{MatSchurComplementAinvType})::PetscErrorCode
end

function MatSchurComplementGetPmat(arg1, arg2, arg3)
    @ccall libpetsc.MatSchurComplementGetPmat(arg1::Mat, arg2::MatReuse, arg3::Ptr{Mat})::PetscErrorCode
end

function MatSchurComplementComputeExplicitOperator(arg1, arg2)
    @ccall libpetsc.MatSchurComplementComputeExplicitOperator(arg1::Mat, arg2::Ptr{Mat})::PetscErrorCode
end

function MatGetSchurComplement(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.MatGetSchurComplement(arg1::Mat, arg2::IS, arg3::IS, arg4::IS, arg5::IS, arg6::MatReuse, arg7::Ptr{Mat}, arg8::MatSchurComplementAinvType, arg9::MatReuse, arg10::Ptr{Mat})::PetscErrorCode
end

function MatCreateSchurComplementPmat(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.MatCreateSchurComplementPmat(arg1::Mat, arg2::Mat, arg3::Mat, arg4::Mat, arg5::MatSchurComplementAinvType, arg6::MatReuse, arg7::Ptr{Mat})::PetscErrorCode
end

function MatCreateLMVMDFP(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatCreateLMVMDFP(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Mat})::PetscErrorCode
end

function MatCreateLMVMBFGS(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatCreateLMVMBFGS(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Mat})::PetscErrorCode
end

function MatCreateLMVMSR1(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatCreateLMVMSR1(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Mat})::PetscErrorCode
end

function MatCreateLMVMBroyden(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatCreateLMVMBroyden(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Mat})::PetscErrorCode
end

function MatCreateLMVMBadBroyden(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatCreateLMVMBadBroyden(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Mat})::PetscErrorCode
end

function MatCreateLMVMSymBroyden(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatCreateLMVMSymBroyden(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Mat})::PetscErrorCode
end

function MatCreateLMVMSymBadBroyden(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatCreateLMVMSymBadBroyden(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Mat})::PetscErrorCode
end

function MatCreateLMVMDiagBroyden(arg1, arg2, arg3, arg4)
    @ccall libpetsc.MatCreateLMVMDiagBroyden(arg1::MPI_Comm, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{Mat})::PetscErrorCode
end

function MatLMVMUpdate(arg1, arg2, arg3)
    @ccall libpetsc.MatLMVMUpdate(arg1::Mat, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function MatLMVMIsAllocated(arg1, arg2)
    @ccall libpetsc.MatLMVMIsAllocated(arg1::Mat, arg2::Ptr{PetscBool})::PetscErrorCode
end

function MatLMVMAllocate(arg1, arg2, arg3)
    @ccall libpetsc.MatLMVMAllocate(arg1::Mat, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function MatLMVMReset(arg1, arg2)
    @ccall libpetsc.MatLMVMReset(arg1::Mat, arg2::PetscBool)::PetscErrorCode
end

function MatLMVMResetShift(arg1)
    @ccall libpetsc.MatLMVMResetShift(arg1::Mat)::PetscErrorCode
end

function MatLMVMClearJ0(arg1)
    @ccall libpetsc.MatLMVMClearJ0(arg1::Mat)::PetscErrorCode
end

function MatLMVMSetJ0(arg1, arg2)
    @ccall libpetsc.MatLMVMSetJ0(arg1::Mat, arg2::Mat)::PetscErrorCode
end

function MatLMVMSetJ0Scale(arg1, arg2)
    @ccall libpetsc.MatLMVMSetJ0Scale(arg1::Mat, arg2::PetscReal)::PetscErrorCode
end

function MatLMVMSetJ0Diag(arg1, arg2)
    @ccall libpetsc.MatLMVMSetJ0Diag(arg1::Mat, arg2::Vec)::PetscErrorCode
end

function MatLMVMSetJ0PC(arg1, arg2)
    @ccall libpetsc.MatLMVMSetJ0PC(arg1::Mat, arg2::PC)::PetscErrorCode
end

function MatLMVMSetJ0KSP(arg1, arg2)
    @ccall libpetsc.MatLMVMSetJ0KSP(arg1::Mat, arg2::KSP)::PetscErrorCode
end

function MatLMVMApplyJ0Fwd(arg1, arg2, arg3)
    @ccall libpetsc.MatLMVMApplyJ0Fwd(arg1::Mat, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function MatLMVMApplyJ0Inv(arg1, arg2, arg3)
    @ccall libpetsc.MatLMVMApplyJ0Inv(arg1::Mat, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function MatLMVMGetJ0(arg1, arg2)
    @ccall libpetsc.MatLMVMGetJ0(arg1::Mat, arg2::Ptr{Mat})::PetscErrorCode
end

function MatLMVMGetJ0PC(arg1, arg2)
    @ccall libpetsc.MatLMVMGetJ0PC(arg1::Mat, arg2::Ptr{PC})::PetscErrorCode
end

function MatLMVMGetJ0KSP(arg1, arg2)
    @ccall libpetsc.MatLMVMGetJ0KSP(arg1::Mat, arg2::Ptr{KSP})::PetscErrorCode
end

function MatLMVMSetHistorySize(arg1, arg2)
    @ccall libpetsc.MatLMVMSetHistorySize(arg1::Mat, arg2::PetscInt)::PetscErrorCode
end

function MatLMVMGetUpdateCount(arg1, arg2)
    @ccall libpetsc.MatLMVMGetUpdateCount(arg1::Mat, arg2::Ptr{PetscInt})::PetscErrorCode
end

function MatLMVMGetRejectCount(arg1, arg2)
    @ccall libpetsc.MatLMVMGetRejectCount(arg1::Mat, arg2::Ptr{PetscInt})::PetscErrorCode
end

function MatLMVMSymBroydenSetDelta(arg1, arg2)
    @ccall libpetsc.MatLMVMSymBroydenSetDelta(arg1::Mat, arg2::PetscScalar)::PetscErrorCode
end

function MatLMVMSymBroydenSetScaleType(arg1, arg2)
    @ccall libpetsc.MatLMVMSymBroydenSetScaleType(arg1::Mat, arg2::MatLMVMSymBroydenScaleType)::PetscErrorCode
end

function KSPSetDM(arg1, arg2)
    @ccall libpetsc.KSPSetDM(arg1::KSP, arg2::DM)::PetscErrorCode
end

function KSPSetDMActive(arg1, arg2)
    @ccall libpetsc.KSPSetDMActive(arg1::KSP, arg2::PetscBool)::PetscErrorCode
end

function KSPGetDM(arg1, arg2)
    @ccall libpetsc.KSPGetDM(arg1::KSP, arg2::Ptr{DM})::PetscErrorCode
end

function KSPSetApplicationContext(arg1, arg2)
    @ccall libpetsc.KSPSetApplicationContext(arg1::KSP, arg2::Ptr{Cvoid})::PetscErrorCode
end

function KSPGetApplicationContext(arg1, arg2)
    @ccall libpetsc.KSPGetApplicationContext(arg1::KSP, arg2::Ptr{Cvoid})::PetscErrorCode
end

function KSPSetComputeRHS(arg1, func, arg3)
    @ccall libpetsc.KSPSetComputeRHS(arg1::KSP, func::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function KSPSetComputeOperators(arg1, arg2, arg3)
    @ccall libpetsc.KSPSetComputeOperators(arg1::KSP, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function KSPSetComputeInitialGuess(arg1, arg2, arg3)
    @ccall libpetsc.KSPSetComputeInitialGuess(arg1::KSP, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMKSPSetComputeOperators(arg1, arg2, arg3)
    @ccall libpetsc.DMKSPSetComputeOperators(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMKSPGetComputeOperators(arg1, arg2, arg3)
    @ccall libpetsc.DMKSPGetComputeOperators(arg1::DM, arg2::Ptr{Ptr{Cvoid}}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMKSPSetComputeRHS(arg1, arg2, arg3)
    @ccall libpetsc.DMKSPSetComputeRHS(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMKSPGetComputeRHS(arg1, arg2, arg3)
    @ccall libpetsc.DMKSPGetComputeRHS(arg1::DM, arg2::Ptr{Ptr{Cvoid}}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMKSPSetComputeInitialGuess(arg1, arg2, arg3)
    @ccall libpetsc.DMKSPSetComputeInitialGuess(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMKSPGetComputeInitialGuess(arg1, arg2, arg3)
    @ccall libpetsc.DMKSPGetComputeInitialGuess(arg1::DM, arg2::Ptr{Ptr{Cvoid}}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMGlobalToLocalSolve(arg1, arg2, arg3)
    @ccall libpetsc.DMGlobalToLocalSolve(arg1::DM, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function DMProjectField(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMProjectField(arg1::DM, arg2::PetscReal, arg3::Vec, arg4::Ptr{Ptr{Cvoid}}, arg5::InsertMode, arg6::Vec)::PetscErrorCode
end

function DMAdaptInterpolator(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.DMAdaptInterpolator(arg1::DM, arg2::DM, arg3::Mat, arg4::KSP, arg5::PetscInt, arg6::Ptr{Vec}, arg7::Ptr{Vec}, arg8::Ptr{Mat}, arg9::Ptr{Cvoid})::PetscErrorCode
end

function DMCheckInterpolator(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMCheckInterpolator(arg1::DM, arg2::Mat, arg3::PetscInt, arg4::Ptr{Vec}, arg5::Ptr{Vec}, arg6::PetscReal)::PetscErrorCode
end

mutable struct _p_SNES end

const SNES = Ptr{_p_SNES}

const SNESType = Ptr{Cchar}

function SNESInitializePackage()
    @ccall libpetsc.SNESInitializePackage()::PetscErrorCode
end

function SNESCreate(arg1, arg2)
    @ccall libpetsc.SNESCreate(arg1::MPI_Comm, arg2::Ptr{SNES})::PetscErrorCode
end

function SNESReset(arg1)
    @ccall libpetsc.SNESReset(arg1::SNES)::PetscErrorCode
end

function SNESDestroy(arg1)
    @ccall libpetsc.SNESDestroy(arg1::Ptr{SNES})::PetscErrorCode
end

function SNESSetType(arg1, arg2)
    @ccall libpetsc.SNESSetType(arg1::SNES, arg2::SNESType)::PetscErrorCode
end

function SNESMonitor(arg1, arg2, arg3)
    @ccall libpetsc.SNESMonitor(arg1::SNES, arg2::PetscInt, arg3::PetscReal)::PetscErrorCode
end

function SNESMonitorSet(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESMonitorSet(arg1::SNES, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function SNESMonitorSetFromOptions(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.SNESMonitorSetFromOptions(arg1::SNES, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{Cvoid}, arg6::Ptr{Cvoid})::PetscErrorCode
end

function SNESMonitorCancel(arg1)
    @ccall libpetsc.SNESMonitorCancel(arg1::SNES)::PetscErrorCode
end

function SNESMonitorSAWs(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESMonitorSAWs(arg1::SNES, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{Cvoid})::PetscErrorCode
end

function SNESMonitorSAWsCreate(arg1, arg2)
    @ccall libpetsc.SNESMonitorSAWsCreate(arg1::SNES, arg2::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function SNESMonitorSAWsDestroy(arg1)
    @ccall libpetsc.SNESMonitorSAWsDestroy(arg1::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function SNESSetConvergenceHistory(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.SNESSetConvergenceHistory(arg1::SNES, arg2::Ptr{PetscReal}, arg3::Ptr{PetscInt}, arg4::PetscInt, arg5::PetscBool)::PetscErrorCode
end

function SNESGetConvergenceHistory(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESGetConvergenceHistory(arg1::SNES, arg2::Ptr{Ptr{PetscReal}}, arg3::Ptr{Ptr{PetscInt}}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function SNESSetUp(arg1)
    @ccall libpetsc.SNESSetUp(arg1::SNES)::PetscErrorCode
end

function SNESSolve(arg1, arg2, arg3)
    @ccall libpetsc.SNESSolve(arg1::SNES, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function SNESSetErrorIfNotConverged(arg1, arg2)
    @ccall libpetsc.SNESSetErrorIfNotConverged(arg1::SNES, arg2::PetscBool)::PetscErrorCode
end

function SNESGetErrorIfNotConverged(arg1, arg2)
    @ccall libpetsc.SNESGetErrorIfNotConverged(arg1::SNES, arg2::Ptr{PetscBool})::PetscErrorCode
end

function SNESSetWorkVecs(arg1, arg2)
    @ccall libpetsc.SNESSetWorkVecs(arg1::SNES, arg2::PetscInt)::PetscErrorCode
end

function SNESAddOptionsChecker(arg1)
    @ccall libpetsc.SNESAddOptionsChecker(arg1::Ptr{Cvoid})::PetscErrorCode
end

function SNESSetUpdate(arg1, arg2)
    @ccall libpetsc.SNESSetUpdate(arg1::SNES, arg2::Ptr{Cvoid})::PetscErrorCode
end

function SNESRegister(arg1, arg2)
    @ccall libpetsc.SNESRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function SNESGetKSP(arg1, arg2)
    @ccall libpetsc.SNESGetKSP(arg1::SNES, arg2::Ptr{KSP})::PetscErrorCode
end

function SNESSetKSP(arg1, arg2)
    @ccall libpetsc.SNESSetKSP(arg1::SNES, arg2::KSP)::PetscErrorCode
end

function SNESSetSolution(arg1, arg2)
    @ccall libpetsc.SNESSetSolution(arg1::SNES, arg2::Vec)::PetscErrorCode
end

function SNESGetSolution(arg1, arg2)
    @ccall libpetsc.SNESGetSolution(arg1::SNES, arg2::Ptr{Vec})::PetscErrorCode
end

function SNESGetSolutionUpdate(arg1, arg2)
    @ccall libpetsc.SNESGetSolutionUpdate(arg1::SNES, arg2::Ptr{Vec})::PetscErrorCode
end

function SNESGetRhs(arg1, arg2)
    @ccall libpetsc.SNESGetRhs(arg1::SNES, arg2::Ptr{Vec})::PetscErrorCode
end

function SNESView(arg1, arg2)
    @ccall libpetsc.SNESView(arg1::SNES, arg2::PetscViewer)::PetscErrorCode
end

function SNESLoad(arg1, arg2)
    @ccall libpetsc.SNESLoad(arg1::SNES, arg2::PetscViewer)::PetscErrorCode
end

function SNESConvergedReasonViewSet(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESConvergedReasonViewSet(arg1::SNES, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function SNESViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.SNESViewFromOptions(arg1::SNES, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function SNESConvergedReasonView(arg1, arg2)
    @ccall libpetsc.SNESConvergedReasonView(arg1::SNES, arg2::PetscViewer)::PetscErrorCode
end

function SNESConvergedReasonViewFromOptions(arg1)
    @ccall libpetsc.SNESConvergedReasonViewFromOptions(arg1::SNES)::PetscErrorCode
end

function SNESConvergedReasonViewCancel(arg1)
    @ccall libpetsc.SNESConvergedReasonViewCancel(arg1::SNES)::PetscErrorCode
end

function SNESReasonView(snes, v)
    @ccall libpetsc.SNESReasonView(snes::SNES, v::PetscViewer)::PetscErrorCode
end

function SNESReasonViewFromOptions(snes)
    @ccall libpetsc.SNESReasonViewFromOptions(snes::SNES)::PetscErrorCode
end

function SNESSetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.SNESSetOptionsPrefix(arg1::SNES, arg2::Ptr{Cchar})::PetscErrorCode
end

function SNESAppendOptionsPrefix(arg1, arg2)
    @ccall libpetsc.SNESAppendOptionsPrefix(arg1::SNES, arg2::Ptr{Cchar})::PetscErrorCode
end

function SNESGetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.SNESGetOptionsPrefix(arg1::SNES, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function SNESSetFromOptions(arg1)
    @ccall libpetsc.SNESSetFromOptions(arg1::SNES)::PetscErrorCode
end

function SNESResetFromOptions(arg1)
    @ccall libpetsc.SNESResetFromOptions(arg1::SNES)::PetscErrorCode
end

function SNESSetUseMatrixFree(arg1, arg2, arg3)
    @ccall libpetsc.SNESSetUseMatrixFree(arg1::SNES, arg2::PetscBool, arg3::PetscBool)::PetscErrorCode
end

function SNESGetUseMatrixFree(arg1, arg2, arg3)
    @ccall libpetsc.SNESGetUseMatrixFree(arg1::SNES, arg2::Ptr{PetscBool}, arg3::Ptr{PetscBool})::PetscErrorCode
end

function MatCreateSNESMF(arg1, arg2)
    @ccall libpetsc.MatCreateSNESMF(arg1::SNES, arg2::Ptr{Mat})::PetscErrorCode
end

function MatSNESMFGetSNES(arg1, arg2)
    @ccall libpetsc.MatSNESMFGetSNES(arg1::Mat, arg2::Ptr{SNES})::PetscErrorCode
end

function MatSNESMFSetReuseBase(arg1, arg2)
    @ccall libpetsc.MatSNESMFSetReuseBase(arg1::Mat, arg2::PetscBool)::PetscErrorCode
end

function MatSNESMFGetReuseBase(arg1, arg2)
    @ccall libpetsc.MatSNESMFGetReuseBase(arg1::Mat, arg2::Ptr{PetscBool})::PetscErrorCode
end

function MatMFFDComputeJacobian(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.MatMFFDComputeJacobian(arg1::SNES, arg2::Vec, arg3::Mat, arg4::Mat, arg5::Ptr{Cvoid})::PetscErrorCode
end

function SNESGetType(arg1, arg2)
    @ccall libpetsc.SNESGetType(arg1::SNES, arg2::Ptr{SNESType})::PetscErrorCode
end

function SNESMonitorDefault(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESMonitorDefault(arg1::SNES, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function SNESMonitorScaling(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESMonitorScaling(arg1::SNES, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function SNESMonitorRange(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESMonitorRange(arg1::SNES, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function SNESMonitorRatio(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESMonitorRatio(arg1::SNES, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function SNESMonitorRatioSetUp(arg1, arg2)
    @ccall libpetsc.SNESMonitorRatioSetUp(arg1::SNES, arg2::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function SNESMonitorSolution(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESMonitorSolution(arg1::SNES, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function SNESMonitorResidual(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESMonitorResidual(arg1::SNES, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function SNESMonitorSolutionUpdate(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESMonitorSolutionUpdate(arg1::SNES, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function SNESMonitorDefaultShort(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESMonitorDefaultShort(arg1::SNES, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function SNESMonitorDefaultField(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESMonitorDefaultField(arg1::SNES, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function SNESMonitorJacUpdateSpectrum(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESMonitorJacUpdateSpectrum(arg1::SNES, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function SNESMonitorFields(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESMonitorFields(arg1::SNES, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPMonitorSNESResidual(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorSNESResidual(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPMonitorSNESResidualDrawLG(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorSNESResidualDrawLG(arg1::KSP, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function KSPMonitorSNESResidualDrawLGCreate(arg1, arg2, arg3, arg4)
    @ccall libpetsc.KSPMonitorSNESResidualDrawLGCreate(arg1::PetscViewer, arg2::PetscViewerFormat, arg3::Ptr{Cvoid}, arg4::Ptr{Ptr{PetscViewerAndFormat}})::PetscErrorCode
end

function SNESSetTolerances(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.SNESSetTolerances(arg1::SNES, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscInt, arg6::PetscInt)::PetscErrorCode
end

function SNESSetDivergenceTolerance(arg1, arg2)
    @ccall libpetsc.SNESSetDivergenceTolerance(arg1::SNES, arg2::PetscReal)::PetscErrorCode
end

function SNESGetTolerances(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.SNESGetTolerances(arg1::SNES, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscInt})::PetscErrorCode
end

function SNESGetDivergenceTolerance(arg1, arg2)
    @ccall libpetsc.SNESGetDivergenceTolerance(arg1::SNES, arg2::Ptr{PetscReal})::PetscErrorCode
end

function SNESSetTrustRegionTolerance(arg1, arg2)
    @ccall libpetsc.SNESSetTrustRegionTolerance(arg1::SNES, arg2::PetscReal)::PetscErrorCode
end

function SNESGetForceIteration(arg1, arg2)
    @ccall libpetsc.SNESGetForceIteration(arg1::SNES, arg2::Ptr{PetscBool})::PetscErrorCode
end

function SNESSetForceIteration(arg1, arg2)
    @ccall libpetsc.SNESSetForceIteration(arg1::SNES, arg2::PetscBool)::PetscErrorCode
end

function SNESGetIterationNumber(arg1, arg2)
    @ccall libpetsc.SNESGetIterationNumber(arg1::SNES, arg2::Ptr{PetscInt})::PetscErrorCode
end

function SNESSetIterationNumber(arg1, arg2)
    @ccall libpetsc.SNESSetIterationNumber(arg1::SNES, arg2::PetscInt)::PetscErrorCode
end

function SNESNewtonTRSetPreCheck(arg1, arg2, ctx)
    @ccall libpetsc.SNESNewtonTRSetPreCheck(arg1::SNES, arg2::Ptr{Cvoid}, ctx::Ptr{Cvoid})::PetscErrorCode
end

function SNESNewtonTRGetPreCheck(arg1, arg2, ctx)
    @ccall libpetsc.SNESNewtonTRGetPreCheck(arg1::SNES, arg2::Ptr{Ptr{Cvoid}}, ctx::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function SNESNewtonTRSetPostCheck(arg1, arg2, ctx)
    @ccall libpetsc.SNESNewtonTRSetPostCheck(arg1::SNES, arg2::Ptr{Cvoid}, ctx::Ptr{Cvoid})::PetscErrorCode
end

function SNESNewtonTRGetPostCheck(arg1, arg2, ctx)
    @ccall libpetsc.SNESNewtonTRGetPostCheck(arg1::SNES, arg2::Ptr{Ptr{Cvoid}}, ctx::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function SNESGetNonlinearStepFailures(arg1, arg2)
    @ccall libpetsc.SNESGetNonlinearStepFailures(arg1::SNES, arg2::Ptr{PetscInt})::PetscErrorCode
end

function SNESSetMaxNonlinearStepFailures(arg1, arg2)
    @ccall libpetsc.SNESSetMaxNonlinearStepFailures(arg1::SNES, arg2::PetscInt)::PetscErrorCode
end

function SNESGetMaxNonlinearStepFailures(arg1, arg2)
    @ccall libpetsc.SNESGetMaxNonlinearStepFailures(arg1::SNES, arg2::Ptr{PetscInt})::PetscErrorCode
end

function SNESGetNumberFunctionEvals(arg1, arg2)
    @ccall libpetsc.SNESGetNumberFunctionEvals(arg1::SNES, arg2::Ptr{PetscInt})::PetscErrorCode
end

function SNESSetLagPreconditioner(arg1, arg2)
    @ccall libpetsc.SNESSetLagPreconditioner(arg1::SNES, arg2::PetscInt)::PetscErrorCode
end

function SNESGetLagPreconditioner(arg1, arg2)
    @ccall libpetsc.SNESGetLagPreconditioner(arg1::SNES, arg2::Ptr{PetscInt})::PetscErrorCode
end

function SNESSetLagJacobian(arg1, arg2)
    @ccall libpetsc.SNESSetLagJacobian(arg1::SNES, arg2::PetscInt)::PetscErrorCode
end

function SNESGetLagJacobian(arg1, arg2)
    @ccall libpetsc.SNESGetLagJacobian(arg1::SNES, arg2::Ptr{PetscInt})::PetscErrorCode
end

function SNESSetLagPreconditionerPersists(arg1, arg2)
    @ccall libpetsc.SNESSetLagPreconditionerPersists(arg1::SNES, arg2::PetscBool)::PetscErrorCode
end

function SNESSetLagJacobianPersists(arg1, arg2)
    @ccall libpetsc.SNESSetLagJacobianPersists(arg1::SNES, arg2::PetscBool)::PetscErrorCode
end

function SNESSetGridSequence(arg1, arg2)
    @ccall libpetsc.SNESSetGridSequence(arg1::SNES, arg2::PetscInt)::PetscErrorCode
end

function SNESGetGridSequence(arg1, arg2)
    @ccall libpetsc.SNESGetGridSequence(arg1::SNES, arg2::Ptr{PetscInt})::PetscErrorCode
end

function SNESGetLinearSolveIterations(arg1, arg2)
    @ccall libpetsc.SNESGetLinearSolveIterations(arg1::SNES, arg2::Ptr{PetscInt})::PetscErrorCode
end

function SNESGetLinearSolveFailures(arg1, arg2)
    @ccall libpetsc.SNESGetLinearSolveFailures(arg1::SNES, arg2::Ptr{PetscInt})::PetscErrorCode
end

function SNESSetMaxLinearSolveFailures(arg1, arg2)
    @ccall libpetsc.SNESSetMaxLinearSolveFailures(arg1::SNES, arg2::PetscInt)::PetscErrorCode
end

function SNESGetMaxLinearSolveFailures(arg1, arg2)
    @ccall libpetsc.SNESGetMaxLinearSolveFailures(arg1::SNES, arg2::Ptr{PetscInt})::PetscErrorCode
end

function SNESSetCountersReset(arg1, arg2)
    @ccall libpetsc.SNESSetCountersReset(arg1::SNES, arg2::PetscBool)::PetscErrorCode
end

function SNESKSPSetUseEW(arg1, arg2)
    @ccall libpetsc.SNESKSPSetUseEW(arg1::SNES, arg2::PetscBool)::PetscErrorCode
end

function SNESKSPGetUseEW(arg1, arg2)
    @ccall libpetsc.SNESKSPGetUseEW(arg1::SNES, arg2::Ptr{PetscBool})::PetscErrorCode
end

function SNESKSPSetParametersEW(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.SNESKSPSetParametersEW(arg1::SNES, arg2::PetscInt, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal, arg6::PetscReal, arg7::PetscReal, arg8::PetscReal)::PetscErrorCode
end

function SNESKSPGetParametersEW(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.SNESKSPGetParametersEW(arg1::SNES, arg2::Ptr{PetscInt}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal}, arg7::Ptr{PetscReal}, arg8::Ptr{PetscReal})::PetscErrorCode
end

function SNESMonitorLGRange(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESMonitorLGRange(arg1::SNES, arg2::PetscInt, arg3::PetscReal, arg4::Ptr{Cvoid})::PetscErrorCode
end

function SNESSetApplicationContext(arg1, arg2)
    @ccall libpetsc.SNESSetApplicationContext(arg1::SNES, arg2::Ptr{Cvoid})::PetscErrorCode
end

function SNESGetApplicationContext(arg1, arg2)
    @ccall libpetsc.SNESGetApplicationContext(arg1::SNES, arg2::Ptr{Cvoid})::PetscErrorCode
end

function SNESSetComputeApplicationContext(arg1, arg2, arg3)
    @ccall libpetsc.SNESSetComputeApplicationContext(arg1::SNES, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function SNESPythonSetType(arg1, arg2)
    @ccall libpetsc.SNESPythonSetType(arg1::SNES, arg2::Ptr{Cchar})::PetscErrorCode
end

function SNESSetFunctionDomainError(arg1)
    @ccall libpetsc.SNESSetFunctionDomainError(arg1::SNES)::PetscErrorCode
end

function SNESGetFunctionDomainError(arg1, arg2)
    @ccall libpetsc.SNESGetFunctionDomainError(arg1::SNES, arg2::Ptr{PetscBool})::PetscErrorCode
end

function SNESGetJacobianDomainError(arg1, arg2)
    @ccall libpetsc.SNESGetJacobianDomainError(arg1::SNES, arg2::Ptr{PetscBool})::PetscErrorCode
end

function SNESSetJacobianDomainError(arg1)
    @ccall libpetsc.SNESSetJacobianDomainError(arg1::SNES)::PetscErrorCode
end

function SNESSetCheckJacobianDomainError(arg1, arg2)
    @ccall libpetsc.SNESSetCheckJacobianDomainError(arg1::SNES, arg2::PetscBool)::PetscErrorCode
end

function SNESGetCheckJacobianDomainError(arg1, arg2)
    @ccall libpetsc.SNESGetCheckJacobianDomainError(arg1::SNES, arg2::Ptr{PetscBool})::PetscErrorCode
end

@enum SNESConvergedReason::Int32 begin
    SNES_CONVERGED_FNORM_ABS = 2
    SNES_CONVERGED_FNORM_RELATIVE = 3
    SNES_CONVERGED_SNORM_RELATIVE = 4
    SNES_CONVERGED_ITS = 5
    SNES_DIVERGED_FUNCTION_DOMAIN = -1
    SNES_DIVERGED_FUNCTION_COUNT = -2
    SNES_DIVERGED_LINEAR_SOLVE = -3
    SNES_DIVERGED_FNORM_NAN = -4
    SNES_DIVERGED_MAX_IT = -5
    SNES_DIVERGED_LINE_SEARCH = -6
    SNES_DIVERGED_INNER = -7
    SNES_DIVERGED_LOCAL_MIN = -8
    SNES_DIVERGED_DTOL = -9
    SNES_DIVERGED_JACOBIAN_DOMAIN = -10
    SNES_DIVERGED_TR_DELTA = -11
    # SNES_CONVERGED_TR_DELTA = -11
    SNES_CONVERGED_ITERATING = 0
end

function SNESSetConvergenceTest(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESSetConvergenceTest(arg1::SNES, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function SNESConvergedDefault(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.SNESConvergedDefault(arg1::SNES, arg2::PetscInt, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal, arg6::Ptr{SNESConvergedReason}, arg7::Ptr{Cvoid})::PetscErrorCode
end

function SNESConvergedSkip(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.SNESConvergedSkip(arg1::SNES, arg2::PetscInt, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal, arg6::Ptr{SNESConvergedReason}, arg7::Ptr{Cvoid})::PetscErrorCode
end

function SNESConvergedCorrectPressure(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.SNESConvergedCorrectPressure(arg1::SNES, arg2::PetscInt, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal, arg6::Ptr{SNESConvergedReason}, arg7::Ptr{Cvoid})::PetscErrorCode
end

function SNESGetConvergedReason(arg1, arg2)
    @ccall libpetsc.SNESGetConvergedReason(arg1::SNES, arg2::Ptr{SNESConvergedReason})::PetscErrorCode
end

function SNESGetConvergedReasonString(arg1, arg2)
    @ccall libpetsc.SNESGetConvergedReasonString(arg1::SNES, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function SNESSetConvergedReason(arg1, arg2)
    @ccall libpetsc.SNESSetConvergedReason(arg1::SNES, arg2::SNESConvergedReason)::PetscErrorCode
end

function SNESSetFunction(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESSetFunction(arg1::SNES, arg2::Vec, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function SNESGetFunction(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESGetFunction(arg1::SNES, arg2::Ptr{Vec}, arg3::Ptr{Ptr{Cvoid}}, arg4::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function SNESComputeFunction(arg1, arg2, arg3)
    @ccall libpetsc.SNESComputeFunction(arg1::SNES, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function SNESSetInitialFunction(arg1, arg2)
    @ccall libpetsc.SNESSetInitialFunction(arg1::SNES, arg2::Vec)::PetscErrorCode
end

function SNESSetJacobian(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.SNESSetJacobian(arg1::SNES, arg2::Mat, arg3::Mat, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid})::PetscErrorCode
end

function SNESGetJacobian(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.SNESGetJacobian(arg1::SNES, arg2::Ptr{Mat}, arg3::Ptr{Mat}, arg4::Ptr{Ptr{Cvoid}}, arg5::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function SNESObjectiveComputeFunctionDefaultFD(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESObjectiveComputeFunctionDefaultFD(arg1::SNES, arg2::Vec, arg3::Vec, arg4::Ptr{Cvoid})::PetscErrorCode
end

function SNESComputeJacobianDefault(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.SNESComputeJacobianDefault(arg1::SNES, arg2::Vec, arg3::Mat, arg4::Mat, arg5::Ptr{Cvoid})::PetscErrorCode
end

function SNESComputeJacobianDefaultColor(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.SNESComputeJacobianDefaultColor(arg1::SNES, arg2::Vec, arg3::Mat, arg4::Mat, arg5::Ptr{Cvoid})::PetscErrorCode
end

function SNESSetComputeInitialGuess(arg1, arg2, arg3)
    @ccall libpetsc.SNESSetComputeInitialGuess(arg1::SNES, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function SNESSetPicard(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.SNESSetPicard(arg1::SNES, arg2::Vec, arg3::Ptr{Cvoid}, arg4::Mat, arg5::Mat, arg6::Ptr{Cvoid}, arg7::Ptr{Cvoid})::PetscErrorCode
end

function SNESGetPicard(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.SNESGetPicard(arg1::SNES, arg2::Ptr{Vec}, arg3::Ptr{Ptr{Cvoid}}, arg4::Ptr{Mat}, arg5::Ptr{Mat}, arg6::Ptr{Ptr{Cvoid}}, arg7::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function SNESPicardComputeFunction(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESPicardComputeFunction(arg1::SNES, arg2::Vec, arg3::Vec, arg4::Ptr{Cvoid})::PetscErrorCode
end

function SNESPicardComputeJacobian(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.SNESPicardComputeJacobian(arg1::SNES, arg2::Vec, arg3::Mat, arg4::Mat, arg5::Ptr{Cvoid})::PetscErrorCode
end

function SNESSetObjective(arg1, arg2, arg3)
    @ccall libpetsc.SNESSetObjective(arg1::SNES, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function SNESGetObjective(arg1, arg2, arg3)
    @ccall libpetsc.SNESGetObjective(arg1::SNES, arg2::Ptr{Ptr{Cvoid}}, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function SNESComputeObjective(arg1, arg2, arg3)
    @ccall libpetsc.SNESComputeObjective(arg1::SNES, arg2::Vec, arg3::Ptr{PetscReal})::PetscErrorCode
end

@enum SNESNormSchedule::Int32 begin
    SNES_NORM_DEFAULT = -1
    SNES_NORM_NONE = 0
    SNES_NORM_ALWAYS = 1
    SNES_NORM_INITIAL_ONLY = 2
    SNES_NORM_FINAL_ONLY = 3
    SNES_NORM_INITIAL_FINAL_ONLY = 4
end

function SNESSetNormSchedule(arg1, arg2)
    @ccall libpetsc.SNESSetNormSchedule(arg1::SNES, arg2::SNESNormSchedule)::PetscErrorCode
end

function SNESGetNormSchedule(arg1, arg2)
    @ccall libpetsc.SNESGetNormSchedule(arg1::SNES, arg2::Ptr{SNESNormSchedule})::PetscErrorCode
end

function SNESSetFunctionNorm(arg1, arg2)
    @ccall libpetsc.SNESSetFunctionNorm(arg1::SNES, arg2::PetscReal)::PetscErrorCode
end

function SNESGetFunctionNorm(arg1, arg2)
    @ccall libpetsc.SNESGetFunctionNorm(arg1::SNES, arg2::Ptr{PetscReal})::PetscErrorCode
end

function SNESGetUpdateNorm(arg1, arg2)
    @ccall libpetsc.SNESGetUpdateNorm(arg1::SNES, arg2::Ptr{PetscReal})::PetscErrorCode
end

function SNESGetSolutionNorm(arg1, arg2)
    @ccall libpetsc.SNESGetSolutionNorm(arg1::SNES, arg2::Ptr{PetscReal})::PetscErrorCode
end

@enum SNESFunctionType::Int32 begin
    SNES_FUNCTION_DEFAULT = -1
    SNES_FUNCTION_UNPRECONDITIONED = 0
    SNES_FUNCTION_PRECONDITIONED = 1
end

function SNESSetFunctionType(arg1, arg2)
    @ccall libpetsc.SNESSetFunctionType(arg1::SNES, arg2::SNESFunctionType)::PetscErrorCode
end

function SNESGetFunctionType(arg1, arg2)
    @ccall libpetsc.SNESGetFunctionType(arg1::SNES, arg2::Ptr{SNESFunctionType})::PetscErrorCode
end

function SNESSetNGS(arg1, arg2, arg3)
    @ccall libpetsc.SNESSetNGS(arg1::SNES, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function SNESGetNGS(arg1, arg2, arg3)
    @ccall libpetsc.SNESGetNGS(arg1::SNES, arg2::Ptr{Ptr{Cvoid}}, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function SNESComputeNGS(arg1, arg2, arg3)
    @ccall libpetsc.SNESComputeNGS(arg1::SNES, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function SNESNGSSetSweeps(arg1, arg2)
    @ccall libpetsc.SNESNGSSetSweeps(arg1::SNES, arg2::PetscInt)::PetscErrorCode
end

function SNESNGSGetSweeps(arg1, arg2)
    @ccall libpetsc.SNESNGSGetSweeps(arg1::SNES, arg2::Ptr{PetscInt})::PetscErrorCode
end

function SNESNGSSetTolerances(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.SNESNGSSetTolerances(arg1::SNES, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscInt)::PetscErrorCode
end

function SNESNGSGetTolerances(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.SNESNGSGetTolerances(arg1::SNES, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscInt})::PetscErrorCode
end

function SNESSetAlwaysComputesFinalResidual(arg1, arg2)
    @ccall libpetsc.SNESSetAlwaysComputesFinalResidual(arg1::SNES, arg2::PetscBool)::PetscErrorCode
end

function SNESGetAlwaysComputesFinalResidual(arg1, arg2)
    @ccall libpetsc.SNESGetAlwaysComputesFinalResidual(arg1::SNES, arg2::Ptr{PetscBool})::PetscErrorCode
end

function SNESShellGetContext(arg1, arg2)
    @ccall libpetsc.SNESShellGetContext(arg1::SNES, arg2::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function SNESShellSetContext(arg1, arg2)
    @ccall libpetsc.SNESShellSetContext(arg1::SNES, arg2::Ptr{Cvoid})::PetscErrorCode
end

function SNESShellSetSolve(arg1, arg2)
    @ccall libpetsc.SNESShellSetSolve(arg1::SNES, arg2::Ptr{Cvoid})::PetscErrorCode
end

mutable struct _p_LineSearch end

const SNESLineSearch = Ptr{_p_LineSearch}

const SNESLineSearchType = Ptr{Cchar}

# typedef PetscErrorCode ( * SNESLineSearchVIProjectFunc ) ( SNES , Vec )
const SNESLineSearchVIProjectFunc = Ptr{Cvoid}

# typedef PetscErrorCode ( * SNESLineSearchVINormFunc ) ( SNES , Vec , Vec , PetscReal * )
const SNESLineSearchVINormFunc = Ptr{Cvoid}

# typedef PetscErrorCode ( * SNESLineSearchApplyFunc ) ( SNESLineSearch )
const SNESLineSearchApplyFunc = Ptr{Cvoid}

# typedef PetscErrorCode ( * SNESLineSearchUserFunc ) ( SNESLineSearch , void * )
const SNESLineSearchUserFunc = Ptr{Cvoid}

function SNESLineSearchCreate(arg1, arg2)
    @ccall libpetsc.SNESLineSearchCreate(arg1::MPI_Comm, arg2::Ptr{SNESLineSearch})::PetscErrorCode
end

function SNESLineSearchReset(arg1)
    @ccall libpetsc.SNESLineSearchReset(arg1::SNESLineSearch)::PetscErrorCode
end

function SNESLineSearchView(arg1, arg2)
    @ccall libpetsc.SNESLineSearchView(arg1::SNESLineSearch, arg2::PetscViewer)::PetscErrorCode
end

function SNESLineSearchDestroy(arg1)
    @ccall libpetsc.SNESLineSearchDestroy(arg1::Ptr{SNESLineSearch})::PetscErrorCode
end

function SNESLineSearchGetType(arg1, arg2)
    @ccall libpetsc.SNESLineSearchGetType(arg1::SNESLineSearch, arg2::Ptr{SNESLineSearchType})::PetscErrorCode
end

function SNESLineSearchSetType(arg1, arg2)
    @ccall libpetsc.SNESLineSearchSetType(arg1::SNESLineSearch, arg2::SNESLineSearchType)::PetscErrorCode
end

function SNESLineSearchSetFromOptions(arg1)
    @ccall libpetsc.SNESLineSearchSetFromOptions(arg1::SNESLineSearch)::PetscErrorCode
end

function SNESLineSearchSetFunction(arg1, arg2)
    @ccall libpetsc.SNESLineSearchSetFunction(arg1::SNESLineSearch, arg2::Ptr{Cvoid})::PetscErrorCode
end

function SNESLineSearchSetUp(arg1)
    @ccall libpetsc.SNESLineSearchSetUp(arg1::SNESLineSearch)::PetscErrorCode
end

function SNESLineSearchApply(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.SNESLineSearchApply(arg1::SNESLineSearch, arg2::Vec, arg3::Vec, arg4::Ptr{PetscReal}, arg5::Vec)::PetscErrorCode
end

function SNESLineSearchPreCheck(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESLineSearchPreCheck(arg1::SNESLineSearch, arg2::Vec, arg3::Vec, arg4::Ptr{PetscBool})::PetscErrorCode
end

function SNESLineSearchPostCheck(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.SNESLineSearchPostCheck(arg1::SNESLineSearch, arg2::Vec, arg3::Vec, arg4::Vec, arg5::Ptr{PetscBool}, arg6::Ptr{PetscBool})::PetscErrorCode
end

function SNESLineSearchSetWorkVecs(arg1, arg2)
    @ccall libpetsc.SNESLineSearchSetWorkVecs(arg1::SNESLineSearch, arg2::PetscInt)::PetscErrorCode
end

function SNESLineSearchSetPreCheck(arg1, arg2, ctx)
    @ccall libpetsc.SNESLineSearchSetPreCheck(arg1::SNESLineSearch, arg2::Ptr{Cvoid}, ctx::Ptr{Cvoid})::PetscErrorCode
end

function SNESLineSearchSetPostCheck(arg1, arg2, ctx)
    @ccall libpetsc.SNESLineSearchSetPostCheck(arg1::SNESLineSearch, arg2::Ptr{Cvoid}, ctx::Ptr{Cvoid})::PetscErrorCode
end

function SNESLineSearchGetPreCheck(arg1, arg2, ctx)
    @ccall libpetsc.SNESLineSearchGetPreCheck(arg1::SNESLineSearch, arg2::Ptr{Ptr{Cvoid}}, ctx::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function SNESLineSearchGetPostCheck(arg1, arg2, ctx)
    @ccall libpetsc.SNESLineSearchGetPostCheck(arg1::SNESLineSearch, arg2::Ptr{Ptr{Cvoid}}, ctx::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function SNESLineSearchSetVIFunctions(arg1, arg2, arg3)
    @ccall libpetsc.SNESLineSearchSetVIFunctions(arg1::SNESLineSearch, arg2::SNESLineSearchVIProjectFunc, arg3::SNESLineSearchVINormFunc)::PetscErrorCode
end

function SNESLineSearchGetVIFunctions(arg1, arg2, arg3)
    @ccall libpetsc.SNESLineSearchGetVIFunctions(arg1::SNESLineSearch, arg2::Ptr{SNESLineSearchVIProjectFunc}, arg3::Ptr{SNESLineSearchVINormFunc})::PetscErrorCode
end

function SNESLineSearchSetSNES(arg1, arg2)
    @ccall libpetsc.SNESLineSearchSetSNES(arg1::SNESLineSearch, arg2::SNES)::PetscErrorCode
end

function SNESLineSearchGetSNES(arg1, arg2)
    @ccall libpetsc.SNESLineSearchGetSNES(arg1::SNESLineSearch, arg2::Ptr{SNES})::PetscErrorCode
end

function SNESLineSearchGetTolerances(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.SNESLineSearchGetTolerances(arg1::SNESLineSearch, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal}, arg7::Ptr{PetscInt})::PetscErrorCode
end

function SNESLineSearchSetTolerances(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.SNESLineSearchSetTolerances(arg1::SNESLineSearch, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal, arg6::PetscReal, arg7::PetscInt)::PetscErrorCode
end

function SNESLineSearchPreCheckPicard(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.SNESLineSearchPreCheckPicard(arg1::SNESLineSearch, arg2::Vec, arg3::Vec, arg4::Ptr{PetscBool}, arg5::Ptr{Cvoid})::PetscErrorCode
end

function SNESLineSearchGetLambda(arg1, arg2)
    @ccall libpetsc.SNESLineSearchGetLambda(arg1::SNESLineSearch, arg2::Ptr{PetscReal})::PetscErrorCode
end

function SNESLineSearchSetLambda(arg1, arg2)
    @ccall libpetsc.SNESLineSearchSetLambda(arg1::SNESLineSearch, arg2::PetscReal)::PetscErrorCode
end

function SNESLineSearchGetDamping(arg1, arg2)
    @ccall libpetsc.SNESLineSearchGetDamping(arg1::SNESLineSearch, arg2::Ptr{PetscReal})::PetscErrorCode
end

function SNESLineSearchSetDamping(arg1, arg2)
    @ccall libpetsc.SNESLineSearchSetDamping(arg1::SNESLineSearch, arg2::PetscReal)::PetscErrorCode
end

function SNESLineSearchGetOrder(arg1, order)
    @ccall libpetsc.SNESLineSearchGetOrder(arg1::SNESLineSearch, order::Ptr{PetscInt})::PetscErrorCode
end

function SNESLineSearchSetOrder(arg1, order)
    @ccall libpetsc.SNESLineSearchSetOrder(arg1::SNESLineSearch, order::PetscInt)::PetscErrorCode
end

@enum SNESLineSearchReason::UInt32 begin
    SNES_LINESEARCH_SUCCEEDED = 0
    SNES_LINESEARCH_FAILED_NANORINF = 1
    SNES_LINESEARCH_FAILED_DOMAIN = 2
    SNES_LINESEARCH_FAILED_REDUCT = 3
    SNES_LINESEARCH_FAILED_USER = 4
    SNES_LINESEARCH_FAILED_FUNCTION = 5
end

function SNESLineSearchGetReason(arg1, arg2)
    @ccall libpetsc.SNESLineSearchGetReason(arg1::SNESLineSearch, arg2::Ptr{SNESLineSearchReason})::PetscErrorCode
end

function SNESLineSearchSetReason(arg1, arg2)
    @ccall libpetsc.SNESLineSearchSetReason(arg1::SNESLineSearch, arg2::SNESLineSearchReason)::PetscErrorCode
end

function SNESLineSearchGetVecs(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.SNESLineSearchGetVecs(arg1::SNESLineSearch, arg2::Ptr{Vec}, arg3::Ptr{Vec}, arg4::Ptr{Vec}, arg5::Ptr{Vec}, arg6::Ptr{Vec})::PetscErrorCode
end

function SNESLineSearchSetVecs(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.SNESLineSearchSetVecs(arg1::SNESLineSearch, arg2::Vec, arg3::Vec, arg4::Vec, arg5::Vec, arg6::Vec)::PetscErrorCode
end

function SNESLineSearchGetNorms(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESLineSearchGetNorms(arg1::SNESLineSearch, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal})::PetscErrorCode
end

function SNESLineSearchSetNorms(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESLineSearchSetNorms(arg1::SNESLineSearch, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal)::PetscErrorCode
end

function SNESLineSearchComputeNorms(arg1)
    @ccall libpetsc.SNESLineSearchComputeNorms(arg1::SNESLineSearch)::PetscErrorCode
end

function SNESLineSearchSetComputeNorms(arg1, arg2)
    @ccall libpetsc.SNESLineSearchSetComputeNorms(arg1::SNESLineSearch, arg2::PetscBool)::PetscErrorCode
end

function SNESLineSearchMonitor(arg1)
    @ccall libpetsc.SNESLineSearchMonitor(arg1::SNESLineSearch)::PetscErrorCode
end

function SNESLineSearchMonitorSet(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESLineSearchMonitorSet(arg1::SNESLineSearch, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function SNESLineSearchMonitorSetFromOptions(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.SNESLineSearchMonitorSetFromOptions(arg1::SNESLineSearch, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{Cvoid}, arg6::Ptr{Cvoid})::PetscErrorCode
end

function SNESLineSearchMonitorCancel(arg1)
    @ccall libpetsc.SNESLineSearchMonitorCancel(arg1::SNESLineSearch)::PetscErrorCode
end

function SNESLineSearchMonitorUpdate(arg1, arg2)
    @ccall libpetsc.SNESLineSearchMonitorUpdate(arg1::SNESLineSearch, arg2::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function SNESLineSearchSetDefaultMonitor(arg1, arg2)
    @ccall libpetsc.SNESLineSearchSetDefaultMonitor(arg1::SNESLineSearch, arg2::PetscViewer)::PetscErrorCode
end

function SNESLineSearchGetDefaultMonitor(arg1, arg2)
    @ccall libpetsc.SNESLineSearchGetDefaultMonitor(arg1::SNESLineSearch, arg2::Ptr{PetscViewer})::PetscErrorCode
end

function SNESLineSearchMonitorSolutionUpdate(arg1, arg2)
    @ccall libpetsc.SNESLineSearchMonitorSolutionUpdate(arg1::SNESLineSearch, arg2::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function SNESLineSearchAppendOptionsPrefix(arg1, prefix)
    @ccall libpetsc.SNESLineSearchAppendOptionsPrefix(arg1::SNESLineSearch, prefix::Ptr{Cchar})::PetscErrorCode
end

function SNESLineSearchGetOptionsPrefix(arg1, prefix)
    @ccall libpetsc.SNESLineSearchGetOptionsPrefix(arg1::SNESLineSearch, prefix::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function SNESLineSearchShellSetUserFunc(arg1, arg2, arg3)
    @ccall libpetsc.SNESLineSearchShellSetUserFunc(arg1::SNESLineSearch, arg2::SNESLineSearchUserFunc, arg3::Ptr{Cvoid})::PetscErrorCode
end

function SNESLineSearchShellGetUserFunc(arg1, arg2, arg3)
    @ccall libpetsc.SNESLineSearchShellGetUserFunc(arg1::SNESLineSearch, arg2::Ptr{SNESLineSearchUserFunc}, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function SNESLineSearchBTSetAlpha(arg1, arg2)
    @ccall libpetsc.SNESLineSearchBTSetAlpha(arg1::SNESLineSearch, arg2::PetscReal)::PetscErrorCode
end

function SNESLineSearchBTGetAlpha(arg1, arg2)
    @ccall libpetsc.SNESLineSearchBTGetAlpha(arg1::SNESLineSearch, arg2::Ptr{PetscReal})::PetscErrorCode
end

function SNESLineSearchRegister(arg1, arg2)
    @ccall libpetsc.SNESLineSearchRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function SNESVISetVariableBounds(arg1, arg2, arg3)
    @ccall libpetsc.SNESVISetVariableBounds(arg1::SNES, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function SNESVISetComputeVariableBounds(arg1, arg2)
    @ccall libpetsc.SNESVISetComputeVariableBounds(arg1::SNES, arg2::Ptr{Cvoid})::PetscErrorCode
end

function SNESVIGetInactiveSet(arg1, arg2)
    @ccall libpetsc.SNESVIGetInactiveSet(arg1::SNES, arg2::Ptr{IS})::PetscErrorCode
end

function SNESVIGetActiveSetIS(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESVIGetActiveSetIS(arg1::SNES, arg2::Vec, arg3::Vec, arg4::Ptr{IS})::PetscErrorCode
end

function SNESVIComputeInactiveSetFnorm(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESVIComputeInactiveSetFnorm(arg1::SNES, arg2::Vec, arg3::Vec, arg4::Ptr{PetscReal})::PetscErrorCode
end

function SNESVISetRedundancyCheck(arg1, arg2, arg3)
    @ccall libpetsc.SNESVISetRedundancyCheck(arg1::SNES, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function SNESTestLocalMin(arg1)
    @ccall libpetsc.SNESTestLocalMin(arg1::SNES)::PetscErrorCode
end

function SNESComputeJacobian(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESComputeJacobian(arg1::SNES, arg2::Vec, arg3::Mat, arg4::Mat)::PetscErrorCode
end

function SNESTestJacobian(arg1)
    @ccall libpetsc.SNESTestJacobian(arg1::SNES)::PetscErrorCode
end

function SNESSetDM(arg1, arg2)
    @ccall libpetsc.SNESSetDM(arg1::SNES, arg2::DM)::PetscErrorCode
end

function SNESGetDM(arg1, arg2)
    @ccall libpetsc.SNESGetDM(arg1::SNES, arg2::Ptr{DM})::PetscErrorCode
end

function SNESSetNPC(arg1, arg2)
    @ccall libpetsc.SNESSetNPC(arg1::SNES, arg2::SNES)::PetscErrorCode
end

function SNESGetNPC(arg1, arg2)
    @ccall libpetsc.SNESGetNPC(arg1::SNES, arg2::Ptr{SNES})::PetscErrorCode
end

function SNESHasNPC(arg1, arg2)
    @ccall libpetsc.SNESHasNPC(arg1::SNES, arg2::Ptr{PetscBool})::PetscErrorCode
end

function SNESApplyNPC(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESApplyNPC(arg1::SNES, arg2::Vec, arg3::Vec, arg4::Vec)::PetscErrorCode
end

function SNESGetNPCFunction(arg1, arg2, arg3)
    @ccall libpetsc.SNESGetNPCFunction(arg1::SNES, arg2::Vec, arg3::Ptr{PetscReal})::PetscErrorCode
end

function SNESComputeFunctionDefaultNPC(arg1, arg2, arg3)
    @ccall libpetsc.SNESComputeFunctionDefaultNPC(arg1::SNES, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function SNESSetNPCSide(arg1, arg2)
    @ccall libpetsc.SNESSetNPCSide(arg1::SNES, arg2::PCSide)::PetscErrorCode
end

function SNESGetNPCSide(arg1, arg2)
    @ccall libpetsc.SNESGetNPCSide(arg1::SNES, arg2::Ptr{PCSide})::PetscErrorCode
end

function SNESSetLineSearch(arg1, arg2)
    @ccall libpetsc.SNESSetLineSearch(arg1::SNES, arg2::SNESLineSearch)::PetscErrorCode
end

function SNESGetLineSearch(arg1, arg2)
    @ccall libpetsc.SNESGetLineSearch(arg1::SNES, arg2::Ptr{SNESLineSearch})::PetscErrorCode
end

function SNESRestrictHookAdd(arg1, arg2, arg3)
    @ccall libpetsc.SNESRestrictHookAdd(arg1::SNES, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function SNESGetSNESLineSearch(snes, ls)
    @ccall libpetsc.SNESGetSNESLineSearch(snes::SNES, ls::Ptr{SNESLineSearch})::PetscErrorCode
end

function SNESSetSNESLineSearch(snes, ls)
    @ccall libpetsc.SNESSetSNESLineSearch(snes::SNES, ls::SNESLineSearch)::PetscErrorCode
end

function SNESSetUpMatrices(arg1)
    @ccall libpetsc.SNESSetUpMatrices(arg1::SNES)::PetscErrorCode
end

function DMSNESSetFunction(arg1, arg2, arg3)
    @ccall libpetsc.DMSNESSetFunction(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMSNESGetFunction(arg1, arg2, arg3)
    @ccall libpetsc.DMSNESGetFunction(arg1::DM, arg2::Ptr{Ptr{Cvoid}}, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMSNESSetNGS(arg1, arg2, arg3)
    @ccall libpetsc.DMSNESSetNGS(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMSNESGetNGS(arg1, arg2, arg3)
    @ccall libpetsc.DMSNESGetNGS(arg1::DM, arg2::Ptr{Ptr{Cvoid}}, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMSNESSetJacobian(arg1, arg2, arg3)
    @ccall libpetsc.DMSNESSetJacobian(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMSNESGetJacobian(arg1, arg2, arg3)
    @ccall libpetsc.DMSNESGetJacobian(arg1::DM, arg2::Ptr{Ptr{Cvoid}}, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMSNESSetPicard(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMSNESSetPicard(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMSNESGetPicard(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMSNESGetPicard(arg1::DM, arg2::Ptr{Ptr{Cvoid}}, arg3::Ptr{Ptr{Cvoid}}, arg4::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMSNESSetObjective(arg1, arg2, arg3)
    @ccall libpetsc.DMSNESSetObjective(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMSNESGetObjective(arg1, arg2, arg3)
    @ccall libpetsc.DMSNESGetObjective(arg1::DM, arg2::Ptr{Ptr{Cvoid}}, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMCopyDMSNES(arg1, arg2)
    @ccall libpetsc.DMCopyDMSNES(arg1::DM, arg2::DM)::PetscErrorCode
end

# typedef PetscErrorCode ( * DMDASNESFunction ) ( DMDALocalInfo * , void * , void * , void * )
const DMDASNESFunction = Ptr{Cvoid}

# typedef PetscErrorCode ( * DMDASNESJacobian ) ( DMDALocalInfo * , void * , Mat , Mat , void * )
const DMDASNESJacobian = Ptr{Cvoid}

# typedef PetscErrorCode ( * DMDASNESObjective ) ( DMDALocalInfo * , void * , PetscReal * , void * )
const DMDASNESObjective = Ptr{Cvoid}

function DMDASNESSetFunctionLocal(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMDASNESSetFunctionLocal(arg1::DM, arg2::InsertMode, arg3::DMDASNESFunction, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMDASNESSetJacobianLocal(arg1, arg2, arg3)
    @ccall libpetsc.DMDASNESSetJacobianLocal(arg1::DM, arg2::DMDASNESJacobian, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMDASNESSetObjectiveLocal(arg1, arg2, arg3)
    @ccall libpetsc.DMDASNESSetObjectiveLocal(arg1::DM, arg2::DMDASNESObjective, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMDASNESSetPicardLocal(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMDASNESSetPicardLocal(arg1::DM, arg2::InsertMode, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid})::PetscErrorCode
end

function DMSNESSetBoundaryLocal(arg1, arg2, arg3)
    @ccall libpetsc.DMSNESSetBoundaryLocal(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMSNESSetFunctionLocal(arg1, arg2, arg3)
    @ccall libpetsc.DMSNESSetFunctionLocal(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMSNESSetJacobianLocal(arg1, arg2, arg3)
    @ccall libpetsc.DMSNESSetJacobianLocal(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMSNESGetBoundaryLocal(arg1, arg2, arg3)
    @ccall libpetsc.DMSNESGetBoundaryLocal(arg1::DM, arg2::Ptr{Ptr{Cvoid}}, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMSNESGetFunctionLocal(arg1, arg2, arg3)
    @ccall libpetsc.DMSNESGetFunctionLocal(arg1::DM, arg2::Ptr{Ptr{Cvoid}}, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMSNESGetJacobianLocal(arg1, arg2, arg3)
    @ccall libpetsc.DMSNESGetJacobianLocal(arg1::DM, arg2::Ptr{Ptr{Cvoid}}, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function SNESMultiblockSetFields(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESMultiblockSetFields(arg1::SNES, arg2::Ptr{Cchar}, arg3::PetscInt, arg4::Ptr{PetscInt})::PetscErrorCode
end

function SNESMultiblockSetIS(arg1, arg2, arg3)
    @ccall libpetsc.SNESMultiblockSetIS(arg1::SNES, arg2::Ptr{Cchar}, arg3::IS)::PetscErrorCode
end

function SNESMultiblockSetBlockSize(arg1, arg2)
    @ccall libpetsc.SNESMultiblockSetBlockSize(arg1::SNES, arg2::PetscInt)::PetscErrorCode
end

function SNESMultiblockSetType(arg1, arg2)
    @ccall libpetsc.SNESMultiblockSetType(arg1::SNES, arg2::PCCompositeType)::PetscErrorCode
end

const SNESMSType = Ptr{Cchar}

function SNESMSRegister(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.SNESMSRegister(arg1::SNESMSType, arg2::PetscInt, arg3::PetscInt, arg4::PetscReal, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal}, arg7::Ptr{PetscReal})::PetscErrorCode
end

function SNESMSGetType(arg1, arg2)
    @ccall libpetsc.SNESMSGetType(arg1::SNES, arg2::Ptr{SNESMSType})::PetscErrorCode
end

function SNESMSSetType(arg1, arg2)
    @ccall libpetsc.SNESMSSetType(arg1::SNES, arg2::SNESMSType)::PetscErrorCode
end

function SNESMSGetDamping(arg1, arg2)
    @ccall libpetsc.SNESMSGetDamping(arg1::SNES, arg2::Ptr{PetscReal})::PetscErrorCode
end

function SNESMSSetDamping(arg1, arg2)
    @ccall libpetsc.SNESMSSetDamping(arg1::SNES, arg2::PetscReal)::PetscErrorCode
end

function SNESMSFinalizePackage()
    @ccall libpetsc.SNESMSFinalizePackage()::PetscErrorCode
end

function SNESMSInitializePackage()
    @ccall libpetsc.SNESMSInitializePackage()::PetscErrorCode
end

function SNESMSRegisterDestroy()
    @ccall libpetsc.SNESMSRegisterDestroy()::PetscErrorCode
end

@enum SNESNGMRESRestartType::UInt32 begin
    SNES_NGMRES_RESTART_NONE = 0
    SNES_NGMRES_RESTART_PERIODIC = 1
    SNES_NGMRES_RESTART_DIFFERENCE = 2
end

@enum SNESNGMRESSelectType::UInt32 begin
    SNES_NGMRES_SELECT_NONE = 0
    SNES_NGMRES_SELECT_DIFFERENCE = 1
    SNES_NGMRES_SELECT_LINESEARCH = 2
end

function SNESNGMRESSetRestartType(arg1, arg2)
    @ccall libpetsc.SNESNGMRESSetRestartType(arg1::SNES, arg2::SNESNGMRESRestartType)::PetscErrorCode
end

function SNESNGMRESSetSelectType(arg1, arg2)
    @ccall libpetsc.SNESNGMRESSetSelectType(arg1::SNES, arg2::SNESNGMRESSelectType)::PetscErrorCode
end

function SNESNGMRESSetRestartFmRise(arg1, arg2)
    @ccall libpetsc.SNESNGMRESSetRestartFmRise(arg1::SNES, arg2::PetscBool)::PetscErrorCode
end

function SNESNGMRESGetRestartFmRise(arg1, arg2)
    @ccall libpetsc.SNESNGMRESGetRestartFmRise(arg1::SNES, arg2::Ptr{PetscBool})::PetscErrorCode
end

@enum SNESNCGType::UInt32 begin
    SNES_NCG_FR = 0
    SNES_NCG_PRP = 1
    SNES_NCG_HS = 2
    SNES_NCG_DY = 3
    SNES_NCG_CD = 4
end

function SNESNCGSetType(arg1, arg2)
    @ccall libpetsc.SNESNCGSetType(arg1::SNES, arg2::SNESNCGType)::PetscErrorCode
end

@enum SNESQNScaleType::UInt32 begin
    SNES_QN_SCALE_DEFAULT = 0
    SNES_QN_SCALE_NONE = 1
    SNES_QN_SCALE_SCALAR = 2
    SNES_QN_SCALE_DIAGONAL = 3
    SNES_QN_SCALE_JACOBIAN = 4
end

@enum SNESQNRestartType::UInt32 begin
    SNES_QN_RESTART_DEFAULT = 0
    SNES_QN_RESTART_NONE = 1
    SNES_QN_RESTART_POWELL = 2
    SNES_QN_RESTART_PERIODIC = 3
end

@enum SNESQNType::UInt32 begin
    SNES_QN_LBFGS = 0
    SNES_QN_BROYDEN = 1
    SNES_QN_BADBROYDEN = 2
end

function SNESQNSetType(arg1, arg2)
    @ccall libpetsc.SNESQNSetType(arg1::SNES, arg2::SNESQNType)::PetscErrorCode
end

function SNESQNSetScaleType(arg1, arg2)
    @ccall libpetsc.SNESQNSetScaleType(arg1::SNES, arg2::SNESQNScaleType)::PetscErrorCode
end

function SNESQNSetRestartType(arg1, arg2)
    @ccall libpetsc.SNESQNSetRestartType(arg1::SNES, arg2::SNESQNRestartType)::PetscErrorCode
end

function SNESNASMGetType(arg1, arg2)
    @ccall libpetsc.SNESNASMGetType(arg1::SNES, arg2::Ptr{PCASMType})::PetscErrorCode
end

function SNESNASMSetType(arg1, arg2)
    @ccall libpetsc.SNESNASMSetType(arg1::SNES, arg2::PCASMType)::PetscErrorCode
end

function SNESNASMGetSubdomains(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.SNESNASMGetSubdomains(arg1::SNES, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{SNES}}, arg4::Ptr{Ptr{VecScatter}}, arg5::Ptr{Ptr{VecScatter}}, arg6::Ptr{Ptr{VecScatter}})::PetscErrorCode
end

function SNESNASMSetSubdomains(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.SNESNASMSetSubdomains(arg1::SNES, arg2::PetscInt, arg3::Ptr{SNES}, arg4::Ptr{VecScatter}, arg5::Ptr{VecScatter}, arg6::Ptr{VecScatter})::PetscErrorCode
end

function SNESNASMSetDamping(arg1, arg2)
    @ccall libpetsc.SNESNASMSetDamping(arg1::SNES, arg2::PetscReal)::PetscErrorCode
end

function SNESNASMGetDamping(arg1, arg2)
    @ccall libpetsc.SNESNASMGetDamping(arg1::SNES, arg2::Ptr{PetscReal})::PetscErrorCode
end

function SNESNASMGetSubdomainVecs(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.SNESNASMGetSubdomainVecs(arg1::SNES, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{Vec}}, arg4::Ptr{Ptr{Vec}}, arg5::Ptr{Ptr{Vec}}, arg6::Ptr{Ptr{Vec}})::PetscErrorCode
end

function SNESNASMSetComputeFinalJacobian(arg1, arg2)
    @ccall libpetsc.SNESNASMSetComputeFinalJacobian(arg1::SNES, arg2::PetscBool)::PetscErrorCode
end

function SNESNASMGetSNES(arg1, arg2, arg3)
    @ccall libpetsc.SNESNASMGetSNES(arg1::SNES, arg2::PetscInt, arg3::Ptr{SNES})::PetscErrorCode
end

function SNESNASMGetNumber(arg1, arg2)
    @ccall libpetsc.SNESNASMGetNumber(arg1::SNES, arg2::Ptr{PetscInt})::PetscErrorCode
end

function SNESNASMSetWeight(arg1, arg2)
    @ccall libpetsc.SNESNASMSetWeight(arg1::SNES, arg2::Vec)::PetscErrorCode
end

@enum SNESCompositeType::UInt32 begin
    SNES_COMPOSITE_ADDITIVE = 0
    SNES_COMPOSITE_MULTIPLICATIVE = 1
    SNES_COMPOSITE_ADDITIVEOPTIMAL = 2
end

function SNESCompositeSetType(arg1, arg2)
    @ccall libpetsc.SNESCompositeSetType(arg1::SNES, arg2::SNESCompositeType)::PetscErrorCode
end

function SNESCompositeAddSNES(arg1, arg2)
    @ccall libpetsc.SNESCompositeAddSNES(arg1::SNES, arg2::SNESType)::PetscErrorCode
end

function SNESCompositeGetSNES(arg1, arg2, arg3)
    @ccall libpetsc.SNESCompositeGetSNES(arg1::SNES, arg2::PetscInt, arg3::Ptr{SNES})::PetscErrorCode
end

function SNESCompositeGetNumber(arg1, arg2)
    @ccall libpetsc.SNESCompositeGetNumber(arg1::SNES, arg2::Ptr{PetscInt})::PetscErrorCode
end

function SNESCompositeSetDamping(arg1, arg2, arg3)
    @ccall libpetsc.SNESCompositeSetDamping(arg1::SNES, arg2::PetscInt, arg3::PetscReal)::PetscErrorCode
end

function SNESPatchSetDiscretisationInfo(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11)
    @ccall libpetsc.SNESPatchSetDiscretisationInfo(arg1::SNES, arg2::PetscInt, arg3::Ptr{DM}, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{Ptr{PetscInt}}, arg7::Ptr{PetscInt}, arg8::PetscInt, arg9::Ptr{PetscInt}, arg10::PetscInt, arg11::Ptr{PetscInt})::PetscErrorCode
end

function SNESPatchSetComputeOperator(arg1, func, arg3)
    @ccall libpetsc.SNESPatchSetComputeOperator(arg1::SNES, func::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function SNESPatchSetComputeFunction(arg1, func, arg3)
    @ccall libpetsc.SNESPatchSetComputeFunction(arg1::SNES, func::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function SNESPatchSetConstructType(arg1, arg2, func, arg4)
    @ccall libpetsc.SNESPatchSetConstructType(arg1::SNES, arg2::PCPatchConstructType, func::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function SNESPatchSetCellNumbering(arg1, arg2)
    @ccall libpetsc.SNESPatchSetCellNumbering(arg1::SNES, arg2::PetscSection)::PetscErrorCode
end

@enum SNESFASType::UInt32 begin
    SNES_FAS_MULTIPLICATIVE = 0
    SNES_FAS_ADDITIVE = 1
    SNES_FAS_FULL = 2
    SNES_FAS_KASKADE = 3
end

function SNESFASSetType(arg1, arg2)
    @ccall libpetsc.SNESFASSetType(arg1::SNES, arg2::SNESFASType)::PetscErrorCode
end

function SNESFASGetType(arg1, arg2)
    @ccall libpetsc.SNESFASGetType(arg1::SNES, arg2::Ptr{SNESFASType})::PetscErrorCode
end

function SNESFASSetLevels(arg1, arg2, arg3)
    @ccall libpetsc.SNESFASSetLevels(arg1::SNES, arg2::PetscInt, arg3::Ptr{MPI_Comm})::PetscErrorCode
end

function SNESFASGetLevels(arg1, arg2)
    @ccall libpetsc.SNESFASGetLevels(arg1::SNES, arg2::Ptr{PetscInt})::PetscErrorCode
end

function SNESFASGetCycleSNES(arg1, arg2, arg3)
    @ccall libpetsc.SNESFASGetCycleSNES(arg1::SNES, arg2::PetscInt, arg3::Ptr{SNES})::PetscErrorCode
end

function SNESFASSetNumberSmoothUp(arg1, arg2)
    @ccall libpetsc.SNESFASSetNumberSmoothUp(arg1::SNES, arg2::PetscInt)::PetscErrorCode
end

function SNESFASSetNumberSmoothDown(arg1, arg2)
    @ccall libpetsc.SNESFASSetNumberSmoothDown(arg1::SNES, arg2::PetscInt)::PetscErrorCode
end

function SNESFASSetCycles(arg1, arg2)
    @ccall libpetsc.SNESFASSetCycles(arg1::SNES, arg2::PetscInt)::PetscErrorCode
end

function SNESFASSetMonitor(arg1, arg2, arg3)
    @ccall libpetsc.SNESFASSetMonitor(arg1::SNES, arg2::Ptr{PetscViewerAndFormat}, arg3::PetscBool)::PetscErrorCode
end

function SNESFASSetLog(arg1, arg2)
    @ccall libpetsc.SNESFASSetLog(arg1::SNES, arg2::PetscBool)::PetscErrorCode
end

function SNESFASSetGalerkin(arg1, arg2)
    @ccall libpetsc.SNESFASSetGalerkin(arg1::SNES, arg2::PetscBool)::PetscErrorCode
end

function SNESFASGetGalerkin(arg1, arg2)
    @ccall libpetsc.SNESFASGetGalerkin(arg1::SNES, arg2::Ptr{PetscBool})::PetscErrorCode
end

function SNESFASGalerkinFunctionDefault(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESFASGalerkinFunctionDefault(arg1::SNES, arg2::Vec, arg3::Vec, arg4::Ptr{Cvoid})::PetscErrorCode
end

function SNESFASCycleGetSmoother(arg1, arg2)
    @ccall libpetsc.SNESFASCycleGetSmoother(arg1::SNES, arg2::Ptr{SNES})::PetscErrorCode
end

function SNESFASCycleGetSmootherUp(arg1, arg2)
    @ccall libpetsc.SNESFASCycleGetSmootherUp(arg1::SNES, arg2::Ptr{SNES})::PetscErrorCode
end

function SNESFASCycleGetSmootherDown(arg1, arg2)
    @ccall libpetsc.SNESFASCycleGetSmootherDown(arg1::SNES, arg2::Ptr{SNES})::PetscErrorCode
end

function SNESFASCycleGetCorrection(arg1, arg2)
    @ccall libpetsc.SNESFASCycleGetCorrection(arg1::SNES, arg2::Ptr{SNES})::PetscErrorCode
end

function SNESFASCycleGetInterpolation(arg1, arg2)
    @ccall libpetsc.SNESFASCycleGetInterpolation(arg1::SNES, arg2::Ptr{Mat})::PetscErrorCode
end

function SNESFASCycleGetRestriction(arg1, arg2)
    @ccall libpetsc.SNESFASCycleGetRestriction(arg1::SNES, arg2::Ptr{Mat})::PetscErrorCode
end

function SNESFASCycleGetInjection(arg1, arg2)
    @ccall libpetsc.SNESFASCycleGetInjection(arg1::SNES, arg2::Ptr{Mat})::PetscErrorCode
end

function SNESFASCycleGetRScale(arg1, arg2)
    @ccall libpetsc.SNESFASCycleGetRScale(arg1::SNES, arg2::Ptr{Vec})::PetscErrorCode
end

function SNESFASCycleSetCycles(arg1, arg2)
    @ccall libpetsc.SNESFASCycleSetCycles(arg1::SNES, arg2::PetscInt)::PetscErrorCode
end

function SNESFASCycleIsFine(arg1, arg2)
    @ccall libpetsc.SNESFASCycleIsFine(arg1::SNES, arg2::Ptr{PetscBool})::PetscErrorCode
end

function SNESFASSetInterpolation(arg1, arg2, arg3)
    @ccall libpetsc.SNESFASSetInterpolation(arg1::SNES, arg2::PetscInt, arg3::Mat)::PetscErrorCode
end

function SNESFASGetInterpolation(arg1, arg2, arg3)
    @ccall libpetsc.SNESFASGetInterpolation(arg1::SNES, arg2::PetscInt, arg3::Ptr{Mat})::PetscErrorCode
end

function SNESFASSetRestriction(arg1, arg2, arg3)
    @ccall libpetsc.SNESFASSetRestriction(arg1::SNES, arg2::PetscInt, arg3::Mat)::PetscErrorCode
end

function SNESFASGetRestriction(arg1, arg2, arg3)
    @ccall libpetsc.SNESFASGetRestriction(arg1::SNES, arg2::PetscInt, arg3::Ptr{Mat})::PetscErrorCode
end

function SNESFASSetInjection(arg1, arg2, arg3)
    @ccall libpetsc.SNESFASSetInjection(arg1::SNES, arg2::PetscInt, arg3::Mat)::PetscErrorCode
end

function SNESFASGetInjection(arg1, arg2, arg3)
    @ccall libpetsc.SNESFASGetInjection(arg1::SNES, arg2::PetscInt, arg3::Ptr{Mat})::PetscErrorCode
end

function SNESFASSetRScale(arg1, arg2, arg3)
    @ccall libpetsc.SNESFASSetRScale(arg1::SNES, arg2::PetscInt, arg3::Vec)::PetscErrorCode
end

function SNESFASGetRScale(arg1, arg2, arg3)
    @ccall libpetsc.SNESFASGetRScale(arg1::SNES, arg2::PetscInt, arg3::Ptr{Vec})::PetscErrorCode
end

function SNESFASSetContinuation(arg1, arg2)
    @ccall libpetsc.SNESFASSetContinuation(arg1::SNES, arg2::PetscBool)::PetscErrorCode
end

function SNESFASGetSmoother(arg1, arg2, arg3)
    @ccall libpetsc.SNESFASGetSmoother(arg1::SNES, arg2::PetscInt, arg3::Ptr{SNES})::PetscErrorCode
end

function SNESFASGetSmootherUp(arg1, arg2, arg3)
    @ccall libpetsc.SNESFASGetSmootherUp(arg1::SNES, arg2::PetscInt, arg3::Ptr{SNES})::PetscErrorCode
end

function SNESFASGetSmootherDown(arg1, arg2, arg3)
    @ccall libpetsc.SNESFASGetSmootherDown(arg1::SNES, arg2::PetscInt, arg3::Ptr{SNES})::PetscErrorCode
end

function SNESFASGetCoarseSolve(arg1, arg2)
    @ccall libpetsc.SNESFASGetCoarseSolve(arg1::SNES, arg2::Ptr{SNES})::PetscErrorCode
end

function SNESFASFullSetDownSweep(arg1, arg2)
    @ccall libpetsc.SNESFASFullSetDownSweep(arg1::SNES, arg2::PetscBool)::PetscErrorCode
end

function SNESFASCreateCoarseVec(arg1, arg2)
    @ccall libpetsc.SNESFASCreateCoarseVec(arg1::SNES, arg2::Ptr{Vec})::PetscErrorCode
end

function SNESFASRestrict(arg1, arg2, arg3)
    @ccall libpetsc.SNESFASRestrict(arg1::SNES, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function SNESFASFullSetTotal(arg1, arg2)
    @ccall libpetsc.SNESFASFullSetTotal(arg1::SNES, arg2::PetscBool)::PetscErrorCode
end

function SNESFASFullGetTotal(arg1, arg2)
    @ccall libpetsc.SNESFASFullGetTotal(arg1::SNES, arg2::Ptr{PetscBool})::PetscErrorCode
end

function DMSNESCheckDiscretization(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMSNESCheckDiscretization(arg1::SNES, arg2::DM, arg3::PetscReal, arg4::Vec, arg5::PetscReal, arg6::Ptr{PetscReal})::PetscErrorCode
end

function DMSNESCheckResidual(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.DMSNESCheckResidual(arg1::SNES, arg2::DM, arg3::Vec, arg4::PetscReal, arg5::Ptr{PetscReal})::PetscErrorCode
end

function DMSNESCheckJacobian(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.DMSNESCheckJacobian(arg1::SNES, arg2::DM, arg3::Vec, arg4::PetscReal, arg5::Ptr{PetscBool}, arg6::Ptr{PetscReal})::PetscErrorCode
end

function DMSNESCheckFromOptions(arg1, arg2)
    @ccall libpetsc.DMSNESCheckFromOptions(arg1::SNES, arg2::Vec)::PetscErrorCode
end

mutable struct _p_PetscConvEst end

const PetscConvEst = Ptr{_p_PetscConvEst}

function PetscConvEstCreate(arg1, arg2)
    @ccall libpetsc.PetscConvEstCreate(arg1::MPI_Comm, arg2::Ptr{PetscConvEst})::PetscErrorCode
end

function PetscConvEstDestroy(arg1)
    @ccall libpetsc.PetscConvEstDestroy(arg1::Ptr{PetscConvEst})::PetscErrorCode
end

function PetscConvEstView(arg1, arg2)
    @ccall libpetsc.PetscConvEstView(arg1::PetscConvEst, arg2::PetscViewer)::PetscErrorCode
end

function PetscConvEstSetFromOptions(arg1)
    @ccall libpetsc.PetscConvEstSetFromOptions(arg1::PetscConvEst)::PetscErrorCode
end

function PetscConvEstGetSolver(arg1, arg2)
    @ccall libpetsc.PetscConvEstGetSolver(arg1::PetscConvEst, arg2::Ptr{PetscObject})::PetscErrorCode
end

function PetscConvEstSetSolver(arg1, arg2)
    @ccall libpetsc.PetscConvEstSetSolver(arg1::PetscConvEst, arg2::PetscObject)::PetscErrorCode
end

function PetscConvEstSetUp(arg1)
    @ccall libpetsc.PetscConvEstSetUp(arg1::PetscConvEst)::PetscErrorCode
end

function PetscConvEstComputeInitialGuess(arg1, arg2, arg3, arg4)
    @ccall libpetsc.PetscConvEstComputeInitialGuess(arg1::PetscConvEst, arg2::PetscInt, arg3::DM, arg4::Vec)::PetscErrorCode
end

function PetscConvEstComputeError(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.PetscConvEstComputeError(arg1::PetscConvEst, arg2::PetscInt, arg3::DM, arg4::Vec, arg5::Ptr{PetscReal})::PetscErrorCode
end

function PetscConvEstGetConvRate(arg1, arg2)
    @ccall libpetsc.PetscConvEstGetConvRate(arg1::PetscConvEst, arg2::Ptr{PetscReal})::PetscErrorCode
end

function PetscConvEstMonitorDefault(arg1, arg2)
    @ccall libpetsc.PetscConvEstMonitorDefault(arg1::PetscConvEst, arg2::PetscInt)::PetscErrorCode
end

function PetscConvEstRateView(arg1, arg2, arg3)
    @ccall libpetsc.PetscConvEstRateView(arg1::PetscConvEst, arg2::Ptr{PetscReal}, arg3::PetscViewer)::PetscErrorCode
end

mutable struct _p_TS end

const TS = Ptr{_p_TS}

const TSType = Ptr{Cchar}

@enum TSProblemType::UInt32 begin
    TS_LINEAR = 0
    TS_NONLINEAR = 1
end

@enum TSEquationType::Int32 begin
    TS_EQ_UNSPECIFIED = -1
    TS_EQ_EXPLICIT = 0
    TS_EQ_ODE_EXPLICIT = 1
    TS_EQ_DAE_SEMI_EXPLICIT_INDEX1 = 100
    TS_EQ_DAE_SEMI_EXPLICIT_INDEX2 = 200
    TS_EQ_DAE_SEMI_EXPLICIT_INDEX3 = 300
    TS_EQ_DAE_SEMI_EXPLICIT_INDEXHI = 500
    TS_EQ_IMPLICIT = 1000
    TS_EQ_ODE_IMPLICIT = 1001
    TS_EQ_DAE_IMPLICIT_INDEX1 = 1100
    TS_EQ_DAE_IMPLICIT_INDEX2 = 1200
    TS_EQ_DAE_IMPLICIT_INDEX3 = 1300
    TS_EQ_DAE_IMPLICIT_INDEXHI = 1500
end

@enum TSConvergedReason::Int32 begin
    TS_CONVERGED_ITERATING = 0
    TS_CONVERGED_TIME = 1
    TS_CONVERGED_ITS = 2
    TS_CONVERGED_USER = 3
    TS_CONVERGED_EVENT = 4
    TS_CONVERGED_PSEUDO_FATOL = 5
    TS_CONVERGED_PSEUDO_FRTOL = 6
    TS_DIVERGED_NONLINEAR_SOLVE = -1
    TS_DIVERGED_STEP_REJECTED = -2
    TSFORWARD_DIVERGED_LINEAR_SOLVE = -3
    TSADJOINT_DIVERGED_LINEAR_SOLVE = -4
end

@enum TSExactFinalTimeOption::UInt32 begin
    TS_EXACTFINALTIME_UNSPECIFIED = 0
    TS_EXACTFINALTIME_STEPOVER = 1
    TS_EXACTFINALTIME_INTERPOLATE = 2
    TS_EXACTFINALTIME_MATCHSTEP = 3
end

function TSInitializePackage()
    @ccall libpetsc.TSInitializePackage()::PetscErrorCode
end

function TSFinalizePackage()
    @ccall libpetsc.TSFinalizePackage()::PetscErrorCode
end

function TSCreate(arg1, arg2)
    @ccall libpetsc.TSCreate(arg1::MPI_Comm, arg2::Ptr{TS})::PetscErrorCode
end

function TSClone(arg1, arg2)
    @ccall libpetsc.TSClone(arg1::TS, arg2::Ptr{TS})::PetscErrorCode
end

function TSDestroy(arg1)
    @ccall libpetsc.TSDestroy(arg1::Ptr{TS})::PetscErrorCode
end

function TSSetProblemType(arg1, arg2)
    @ccall libpetsc.TSSetProblemType(arg1::TS, arg2::TSProblemType)::PetscErrorCode
end

function TSGetProblemType(arg1, arg2)
    @ccall libpetsc.TSGetProblemType(arg1::TS, arg2::Ptr{TSProblemType})::PetscErrorCode
end

function TSMonitor(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSMonitor(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec)::PetscErrorCode
end

function TSMonitorSet(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSMonitorSet(arg1::TS, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function TSMonitorSetFromOptions(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TSMonitorSetFromOptions(arg1::TS, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{Cvoid}, arg6::Ptr{Cvoid})::PetscErrorCode
end

function TSMonitorCancel(arg1)
    @ccall libpetsc.TSMonitorCancel(arg1::TS)::PetscErrorCode
end

function TSSetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.TSSetOptionsPrefix(arg1::TS, arg2::Ptr{Cchar})::PetscErrorCode
end

function TSAppendOptionsPrefix(arg1, arg2)
    @ccall libpetsc.TSAppendOptionsPrefix(arg1::TS, arg2::Ptr{Cchar})::PetscErrorCode
end

function TSGetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.TSGetOptionsPrefix(arg1::TS, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function TSSetFromOptions(arg1)
    @ccall libpetsc.TSSetFromOptions(arg1::TS)::PetscErrorCode
end

function TSSetUp(arg1)
    @ccall libpetsc.TSSetUp(arg1::TS)::PetscErrorCode
end

function TSReset(arg1)
    @ccall libpetsc.TSReset(arg1::TS)::PetscErrorCode
end

function TSSetSolution(arg1, arg2)
    @ccall libpetsc.TSSetSolution(arg1::TS, arg2::Vec)::PetscErrorCode
end

function TSGetSolution(arg1, arg2)
    @ccall libpetsc.TSGetSolution(arg1::TS, arg2::Ptr{Vec})::PetscErrorCode
end

function TS2SetSolution(arg1, arg2, arg3)
    @ccall libpetsc.TS2SetSolution(arg1::TS, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function TS2GetSolution(arg1, arg2, arg3)
    @ccall libpetsc.TS2GetSolution(arg1::TS, arg2::Ptr{Vec}, arg3::Ptr{Vec})::PetscErrorCode
end

function TSGetSolutionComponents(arg1, arg2, arg3)
    @ccall libpetsc.TSGetSolutionComponents(arg1::TS, arg2::Ptr{PetscInt}, arg3::Ptr{Vec})::PetscErrorCode
end

function TSGetAuxSolution(arg1, arg2)
    @ccall libpetsc.TSGetAuxSolution(arg1::TS, arg2::Ptr{Vec})::PetscErrorCode
end

function TSGetTimeError(arg1, arg2, arg3)
    @ccall libpetsc.TSGetTimeError(arg1::TS, arg2::PetscInt, arg3::Ptr{Vec})::PetscErrorCode
end

function TSSetTimeError(arg1, arg2)
    @ccall libpetsc.TSSetTimeError(arg1::TS, arg2::Vec)::PetscErrorCode
end

function TSSetRHSJacobianP(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSSetRHSJacobianP(arg1::TS, arg2::Mat, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function TSGetRHSJacobianP(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSGetRHSJacobianP(arg1::TS, arg2::Ptr{Mat}, arg3::Ptr{Ptr{Cvoid}}, arg4::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function TSComputeRHSJacobianP(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSComputeRHSJacobianP(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Mat)::PetscErrorCode
end

function TSSetIJacobianP(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSSetIJacobianP(arg1::TS, arg2::Mat, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function TSComputeIJacobianP(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.TSComputeIJacobianP(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Vec, arg5::PetscReal, arg6::Mat, arg7::PetscBool)::PetscErrorCode
end

function TSComputeDRDPFunction(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSComputeDRDPFunction(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Ptr{Vec})::PetscErrorCode
end

function TSComputeDRDUFunction(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSComputeDRDUFunction(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Ptr{Vec})::PetscErrorCode
end

function TSSetIHessianProduct(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.TSSetIHessianProduct(arg1::TS, arg2::Ptr{Vec}, arg3::Ptr{Cvoid}, arg4::Ptr{Vec}, arg5::Ptr{Cvoid}, arg6::Ptr{Vec}, arg7::Ptr{Cvoid}, arg8::Ptr{Vec}, arg9::Ptr{Cvoid}, arg10::Ptr{Cvoid})::PetscErrorCode
end

function TSComputeIHessianProductFunctionUU(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TSComputeIHessianProductFunctionUU(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Ptr{Vec}, arg5::Vec, arg6::Ptr{Vec})::PetscErrorCode
end

function TSComputeIHessianProductFunctionUP(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TSComputeIHessianProductFunctionUP(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Ptr{Vec}, arg5::Vec, arg6::Ptr{Vec})::PetscErrorCode
end

function TSComputeIHessianProductFunctionPU(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TSComputeIHessianProductFunctionPU(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Ptr{Vec}, arg5::Vec, arg6::Ptr{Vec})::PetscErrorCode
end

function TSComputeIHessianProductFunctionPP(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TSComputeIHessianProductFunctionPP(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Ptr{Vec}, arg5::Vec, arg6::Ptr{Vec})::PetscErrorCode
end

function TSSetRHSHessianProduct(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    @ccall libpetsc.TSSetRHSHessianProduct(arg1::TS, arg2::Ptr{Vec}, arg3::Ptr{Cvoid}, arg4::Ptr{Vec}, arg5::Ptr{Cvoid}, arg6::Ptr{Vec}, arg7::Ptr{Cvoid}, arg8::Ptr{Vec}, arg9::Ptr{Cvoid}, arg10::Ptr{Cvoid})::PetscErrorCode
end

function TSComputeRHSHessianProductFunctionUU(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TSComputeRHSHessianProductFunctionUU(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Ptr{Vec}, arg5::Vec, arg6::Ptr{Vec})::PetscErrorCode
end

function TSComputeRHSHessianProductFunctionUP(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TSComputeRHSHessianProductFunctionUP(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Ptr{Vec}, arg5::Vec, arg6::Ptr{Vec})::PetscErrorCode
end

function TSComputeRHSHessianProductFunctionPU(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TSComputeRHSHessianProductFunctionPU(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Ptr{Vec}, arg5::Vec, arg6::Ptr{Vec})::PetscErrorCode
end

function TSComputeRHSHessianProductFunctionPP(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TSComputeRHSHessianProductFunctionPP(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Ptr{Vec}, arg5::Vec, arg6::Ptr{Vec})::PetscErrorCode
end

function TSSetCostHessianProducts(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSSetCostHessianProducts(arg1::TS, arg2::PetscInt, arg3::Ptr{Vec}, arg4::Ptr{Vec}, arg5::Vec)::PetscErrorCode
end

function TSGetCostHessianProducts(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSGetCostHessianProducts(arg1::TS, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{Vec}}, arg4::Ptr{Ptr{Vec}}, arg5::Ptr{Vec})::PetscErrorCode
end

function TSComputeSNESJacobian(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSComputeSNESJacobian(arg1::TS, arg2::Vec, arg3::Mat, arg4::Mat)::PetscErrorCode
end

mutable struct _p_TSTrajectory end

const TSTrajectory = Ptr{_p_TSTrajectory}

const TSTrajectoryType = Ptr{Cchar}

function TSSetSaveTrajectory(arg1)
    @ccall libpetsc.TSSetSaveTrajectory(arg1::TS)::PetscErrorCode
end

function TSResetTrajectory(arg1)
    @ccall libpetsc.TSResetTrajectory(arg1::TS)::PetscErrorCode
end

function TSTrajectoryCreate(arg1, arg2)
    @ccall libpetsc.TSTrajectoryCreate(arg1::MPI_Comm, arg2::Ptr{TSTrajectory})::PetscErrorCode
end

function TSTrajectoryReset(arg1)
    @ccall libpetsc.TSTrajectoryReset(arg1::TSTrajectory)::PetscErrorCode
end

function TSTrajectoryDestroy(arg1)
    @ccall libpetsc.TSTrajectoryDestroy(arg1::Ptr{TSTrajectory})::PetscErrorCode
end

function TSTrajectoryView(arg1, arg2)
    @ccall libpetsc.TSTrajectoryView(arg1::TSTrajectory, arg2::PetscViewer)::PetscErrorCode
end

function TSTrajectorySetType(arg1, arg2, arg3)
    @ccall libpetsc.TSTrajectorySetType(arg1::TSTrajectory, arg2::TS, arg3::TSTrajectoryType)::PetscErrorCode
end

function TSTrajectoryGetType(arg1, arg2, arg3)
    @ccall libpetsc.TSTrajectoryGetType(arg1::TSTrajectory, arg2::TS, arg3::Ptr{TSTrajectoryType})::PetscErrorCode
end

function TSTrajectorySet(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSTrajectorySet(arg1::TSTrajectory, arg2::TS, arg3::PetscInt, arg4::PetscReal, arg5::Vec)::PetscErrorCode
end

function TSTrajectoryGet(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSTrajectoryGet(arg1::TSTrajectory, arg2::TS, arg3::PetscInt, arg4::Ptr{PetscReal})::PetscErrorCode
end

function TSTrajectoryGetVecs(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TSTrajectoryGetVecs(arg1::TSTrajectory, arg2::TS, arg3::PetscInt, arg4::Ptr{PetscReal}, arg5::Vec, arg6::Vec)::PetscErrorCode
end

function TSTrajectoryGetUpdatedHistoryVecs(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSTrajectoryGetUpdatedHistoryVecs(arg1::TSTrajectory, arg2::TS, arg3::PetscReal, arg4::Ptr{Vec}, arg5::Ptr{Vec})::PetscErrorCode
end

function TSTrajectoryGetNumSteps(arg1, arg2)
    @ccall libpetsc.TSTrajectoryGetNumSteps(arg1::TSTrajectory, arg2::Ptr{PetscInt})::PetscErrorCode
end

function TSTrajectoryRestoreUpdatedHistoryVecs(arg1, arg2, arg3)
    @ccall libpetsc.TSTrajectoryRestoreUpdatedHistoryVecs(arg1::TSTrajectory, arg2::Ptr{Vec}, arg3::Ptr{Vec})::PetscErrorCode
end

function TSTrajectorySetFromOptions(arg1, arg2)
    @ccall libpetsc.TSTrajectorySetFromOptions(arg1::TSTrajectory, arg2::TS)::PetscErrorCode
end

function TSTrajectoryRegister(arg1, arg2)
    @ccall libpetsc.TSTrajectoryRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TSTrajectoryRegisterAll()
    @ccall libpetsc.TSTrajectoryRegisterAll()::PetscErrorCode
end

function TSTrajectorySetUp(arg1, arg2)
    @ccall libpetsc.TSTrajectorySetUp(arg1::TSTrajectory, arg2::TS)::PetscErrorCode
end

function TSTrajectorySetUseHistory(arg1, arg2)
    @ccall libpetsc.TSTrajectorySetUseHistory(arg1::TSTrajectory, arg2::PetscBool)::PetscErrorCode
end

function TSTrajectorySetMonitor(arg1, arg2)
    @ccall libpetsc.TSTrajectorySetMonitor(arg1::TSTrajectory, arg2::PetscBool)::PetscErrorCode
end

function TSTrajectorySetVariableNames(arg1, arg2)
    @ccall libpetsc.TSTrajectorySetVariableNames(arg1::TSTrajectory, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function TSTrajectorySetTransform(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSTrajectorySetTransform(arg1::TSTrajectory, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function TSTrajectorySetSolutionOnly(arg1, arg2)
    @ccall libpetsc.TSTrajectorySetSolutionOnly(arg1::TSTrajectory, arg2::PetscBool)::PetscErrorCode
end

function TSTrajectoryGetSolutionOnly(arg1, arg2)
    @ccall libpetsc.TSTrajectoryGetSolutionOnly(arg1::TSTrajectory, arg2::Ptr{PetscBool})::PetscErrorCode
end

function TSTrajectorySetKeepFiles(arg1, arg2)
    @ccall libpetsc.TSTrajectorySetKeepFiles(arg1::TSTrajectory, arg2::PetscBool)::PetscErrorCode
end

function TSTrajectorySetDirname(arg1, arg2)
    @ccall libpetsc.TSTrajectorySetDirname(arg1::TSTrajectory, arg2::Ptr{Cchar})::PetscErrorCode
end

function TSTrajectorySetFiletemplate(arg1, arg2)
    @ccall libpetsc.TSTrajectorySetFiletemplate(arg1::TSTrajectory, arg2::Ptr{Cchar})::PetscErrorCode
end

function TSGetTrajectory(arg1, arg2)
    @ccall libpetsc.TSGetTrajectory(arg1::TS, arg2::Ptr{TSTrajectory})::PetscErrorCode
end

function TSSetCostGradients(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSSetCostGradients(arg1::TS, arg2::PetscInt, arg3::Ptr{Vec}, arg4::Ptr{Vec})::PetscErrorCode
end

function TSGetCostGradients(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSGetCostGradients(arg1::TS, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{Vec}}, arg4::Ptr{Ptr{Vec}})::PetscErrorCode
end

function TSSetCostIntegrand(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.TSSetCostIntegrand(arg1::TS, arg2::PetscInt, arg3::Vec, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid}, arg6::Ptr{Cvoid}, arg7::PetscBool, arg8::Ptr{Cvoid})::PetscErrorCode
end

function TSGetCostIntegral(arg1, arg2)
    @ccall libpetsc.TSGetCostIntegral(arg1::TS, arg2::Ptr{Vec})::PetscErrorCode
end

function TSComputeCostIntegrand(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSComputeCostIntegrand(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Vec)::PetscErrorCode
end

function TSCreateQuadratureTS(arg1, arg2, arg3)
    @ccall libpetsc.TSCreateQuadratureTS(arg1::TS, arg2::PetscBool, arg3::Ptr{TS})::PetscErrorCode
end

function TSGetQuadratureTS(arg1, arg2, arg3)
    @ccall libpetsc.TSGetQuadratureTS(arg1::TS, arg2::Ptr{PetscBool}, arg3::Ptr{TS})::PetscErrorCode
end

function TSAdjointSetFromOptions(arg1, arg2)
    @ccall libpetsc.TSAdjointSetFromOptions(arg1::Ptr{PetscOptionItems}, arg2::TS)::PetscErrorCode
end

function TSAdjointMonitor(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.TSAdjointMonitor(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec, arg5::PetscInt, arg6::Ptr{Vec}, arg7::Ptr{Vec})::PetscErrorCode
end

function TSAdjointMonitorSet(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSAdjointMonitorSet(arg1::TS, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function TSAdjointMonitorCancel(arg1)
    @ccall libpetsc.TSAdjointMonitorCancel(arg1::TS)::PetscErrorCode
end

function TSAdjointMonitorSetFromOptions(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TSAdjointMonitorSetFromOptions(arg1::TS, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Ptr{Cchar}, arg5::Ptr{Cvoid}, arg6::Ptr{Cvoid})::PetscErrorCode
end

function TSAdjointSetRHSJacobian(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSAdjointSetRHSJacobian(arg1::TS, arg2::Mat, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function TSAdjointComputeRHSJacobian(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSAdjointComputeRHSJacobian(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Mat)::PetscErrorCode
end

function TSAdjointComputeDRDPFunction(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSAdjointComputeDRDPFunction(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Ptr{Vec})::PetscErrorCode
end

function TSAdjointComputeDRDYFunction(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSAdjointComputeDRDYFunction(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Ptr{Vec})::PetscErrorCode
end

function TSAdjointSolve(arg1)
    @ccall libpetsc.TSAdjointSolve(arg1::TS)::PetscErrorCode
end

function TSAdjointSetSteps(arg1, arg2)
    @ccall libpetsc.TSAdjointSetSteps(arg1::TS, arg2::PetscInt)::PetscErrorCode
end

function TSAdjointStep(arg1)
    @ccall libpetsc.TSAdjointStep(arg1::TS)::PetscErrorCode
end

function TSAdjointSetUp(arg1)
    @ccall libpetsc.TSAdjointSetUp(arg1::TS)::PetscErrorCode
end

function TSAdjointReset(arg1)
    @ccall libpetsc.TSAdjointReset(arg1::TS)::PetscErrorCode
end

function TSAdjointCostIntegral(arg1)
    @ccall libpetsc.TSAdjointCostIntegral(arg1::TS)::PetscErrorCode
end

function TSAdjointSetForward(arg1, arg2)
    @ccall libpetsc.TSAdjointSetForward(arg1::TS, arg2::Mat)::PetscErrorCode
end

function TSAdjointResetForward(arg1)
    @ccall libpetsc.TSAdjointResetForward(arg1::TS)::PetscErrorCode
end

function TSForwardSetSensitivities(arg1, arg2, arg3)
    @ccall libpetsc.TSForwardSetSensitivities(arg1::TS, arg2::PetscInt, arg3::Mat)::PetscErrorCode
end

function TSForwardGetSensitivities(arg1, arg2, arg3)
    @ccall libpetsc.TSForwardGetSensitivities(arg1::TS, arg2::Ptr{PetscInt}, arg3::Ptr{Mat})::PetscErrorCode
end

function TSForwardSetIntegralGradients(arg1, arg2, arg3)
    @ccall libpetsc.TSForwardSetIntegralGradients(arg1::TS, arg2::PetscInt, arg3::Ptr{Vec})::PetscErrorCode
end

function TSForwardGetIntegralGradients(arg1, arg2, arg3)
    @ccall libpetsc.TSForwardGetIntegralGradients(arg1::TS, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{Vec}})::PetscErrorCode
end

function TSForwardSetUp(arg1)
    @ccall libpetsc.TSForwardSetUp(arg1::TS)::PetscErrorCode
end

function TSForwardReset(arg1)
    @ccall libpetsc.TSForwardReset(arg1::TS)::PetscErrorCode
end

function TSForwardCostIntegral(arg1)
    @ccall libpetsc.TSForwardCostIntegral(arg1::TS)::PetscErrorCode
end

function TSForwardStep(arg1)
    @ccall libpetsc.TSForwardStep(arg1::TS)::PetscErrorCode
end

function TSForwardSetInitialSensitivities(arg1, arg2)
    @ccall libpetsc.TSForwardSetInitialSensitivities(arg1::TS, arg2::Mat)::PetscErrorCode
end

function TSForwardGetStages(arg1, arg2, arg3)
    @ccall libpetsc.TSForwardGetStages(arg1::TS, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{Mat}})::PetscErrorCode
end

function TSSetMaxSteps(arg1, arg2)
    @ccall libpetsc.TSSetMaxSteps(arg1::TS, arg2::PetscInt)::PetscErrorCode
end

function TSGetMaxSteps(arg1, arg2)
    @ccall libpetsc.TSGetMaxSteps(arg1::TS, arg2::Ptr{PetscInt})::PetscErrorCode
end

function TSSetMaxTime(arg1, arg2)
    @ccall libpetsc.TSSetMaxTime(arg1::TS, arg2::PetscReal)::PetscErrorCode
end

function TSGetMaxTime(arg1, arg2)
    @ccall libpetsc.TSGetMaxTime(arg1::TS, arg2::Ptr{PetscReal})::PetscErrorCode
end

function TSSetExactFinalTime(arg1, arg2)
    @ccall libpetsc.TSSetExactFinalTime(arg1::TS, arg2::TSExactFinalTimeOption)::PetscErrorCode
end

function TSGetExactFinalTime(arg1, arg2)
    @ccall libpetsc.TSGetExactFinalTime(arg1::TS, arg2::Ptr{TSExactFinalTimeOption})::PetscErrorCode
end

function TSSetInitialTimeStep(arg1, arg2, arg3)
    @ccall libpetsc.TSSetInitialTimeStep(arg1::TS, arg2::PetscReal, arg3::PetscReal)::PetscErrorCode
end

function TSSetDuration(arg1, arg2, arg3)
    @ccall libpetsc.TSSetDuration(arg1::TS, arg2::PetscInt, arg3::PetscReal)::PetscErrorCode
end

function TSGetDuration(arg1, arg2, arg3)
    @ccall libpetsc.TSGetDuration(arg1::TS, arg2::Ptr{PetscInt}, arg3::Ptr{PetscReal})::PetscErrorCode
end

function TSGetTimeStepNumber(arg1, arg2)
    @ccall libpetsc.TSGetTimeStepNumber(arg1::TS, arg2::Ptr{PetscInt})::PetscErrorCode
end

function TSGetTotalSteps(arg1, arg2)
    @ccall libpetsc.TSGetTotalSteps(arg1::TS, arg2::Ptr{PetscInt})::PetscErrorCode
end

function TSMonitorDefault(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSMonitorDefault(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec, arg5::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function TSMonitorExtreme(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSMonitorExtreme(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec, arg5::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

mutable struct _n_TSMonitorDrawCtx end

const TSMonitorDrawCtx = Ptr{_n_TSMonitorDrawCtx}

function TSMonitorDrawCtxCreate(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.TSMonitorDrawCtxCreate(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Cint, arg5::Cint, arg6::Cint, arg7::Cint, arg8::PetscInt, arg9::Ptr{TSMonitorDrawCtx})::PetscErrorCode
end

function TSMonitorDrawCtxDestroy(arg1)
    @ccall libpetsc.TSMonitorDrawCtxDestroy(arg1::Ptr{TSMonitorDrawCtx})::PetscErrorCode
end

function TSMonitorDrawSolution(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSMonitorDrawSolution(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TSMonitorDrawSolutionPhase(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSMonitorDrawSolutionPhase(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TSMonitorDrawError(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSMonitorDrawError(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TSMonitorDrawSolutionFunction(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSMonitorDrawSolutionFunction(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TSAdjointMonitorDefault(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.TSAdjointMonitorDefault(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec, arg5::PetscInt, arg6::Ptr{Vec}, arg7::Ptr{Vec}, arg8::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function TSAdjointMonitorDrawSensi(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.TSAdjointMonitorDrawSensi(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec, arg5::PetscInt, arg6::Ptr{Vec}, arg7::Ptr{Vec}, arg8::Ptr{Cvoid})::PetscErrorCode
end

function TSMonitorSolution(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSMonitorSolution(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec, arg5::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

function TSMonitorSolutionVTK(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSMonitorSolutionVTK(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TSMonitorSolutionVTKDestroy(arg1)
    @ccall libpetsc.TSMonitorSolutionVTKDestroy(arg1::Ptr{Cvoid})::PetscErrorCode
end

function TSStep(arg1)
    @ccall libpetsc.TSStep(arg1::TS)::PetscErrorCode
end

function TSEvaluateWLTE(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSEvaluateWLTE(arg1::TS, arg2::NormType, arg3::Ptr{PetscInt}, arg4::Ptr{PetscReal})::PetscErrorCode
end

function TSEvaluateStep(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSEvaluateStep(arg1::TS, arg2::PetscInt, arg3::Vec, arg4::Ptr{PetscBool})::PetscErrorCode
end

function TSSolve(arg1, arg2)
    @ccall libpetsc.TSSolve(arg1::TS, arg2::Vec)::PetscErrorCode
end

function TSGetEquationType(arg1, arg2)
    @ccall libpetsc.TSGetEquationType(arg1::TS, arg2::Ptr{TSEquationType})::PetscErrorCode
end

function TSSetEquationType(arg1, arg2)
    @ccall libpetsc.TSSetEquationType(arg1::TS, arg2::TSEquationType)::PetscErrorCode
end

function TSGetConvergedReason(arg1, arg2)
    @ccall libpetsc.TSGetConvergedReason(arg1::TS, arg2::Ptr{TSConvergedReason})::PetscErrorCode
end

function TSSetConvergedReason(arg1, arg2)
    @ccall libpetsc.TSSetConvergedReason(arg1::TS, arg2::TSConvergedReason)::PetscErrorCode
end

function TSGetSolveTime(arg1, arg2)
    @ccall libpetsc.TSGetSolveTime(arg1::TS, arg2::Ptr{PetscReal})::PetscErrorCode
end

function TSGetSNESIterations(arg1, arg2)
    @ccall libpetsc.TSGetSNESIterations(arg1::TS, arg2::Ptr{PetscInt})::PetscErrorCode
end

function TSGetKSPIterations(arg1, arg2)
    @ccall libpetsc.TSGetKSPIterations(arg1::TS, arg2::Ptr{PetscInt})::PetscErrorCode
end

function TSGetStepRejections(arg1, arg2)
    @ccall libpetsc.TSGetStepRejections(arg1::TS, arg2::Ptr{PetscInt})::PetscErrorCode
end

function TSSetMaxStepRejections(arg1, arg2)
    @ccall libpetsc.TSSetMaxStepRejections(arg1::TS, arg2::PetscInt)::PetscErrorCode
end

function TSGetSNESFailures(arg1, arg2)
    @ccall libpetsc.TSGetSNESFailures(arg1::TS, arg2::Ptr{PetscInt})::PetscErrorCode
end

function TSSetMaxSNESFailures(arg1, arg2)
    @ccall libpetsc.TSSetMaxSNESFailures(arg1::TS, arg2::PetscInt)::PetscErrorCode
end

function TSSetErrorIfStepFails(arg1, arg2)
    @ccall libpetsc.TSSetErrorIfStepFails(arg1::TS, arg2::PetscBool)::PetscErrorCode
end

function TSRestartStep(arg1)
    @ccall libpetsc.TSRestartStep(arg1::TS)::PetscErrorCode
end

function TSRollBack(arg1)
    @ccall libpetsc.TSRollBack(arg1::TS)::PetscErrorCode
end

function TSGetStages(arg1, arg2, arg3)
    @ccall libpetsc.TSGetStages(arg1::TS, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{Vec}})::PetscErrorCode
end

function TSGetTime(arg1, arg2)
    @ccall libpetsc.TSGetTime(arg1::TS, arg2::Ptr{PetscReal})::PetscErrorCode
end

function TSSetTime(arg1, arg2)
    @ccall libpetsc.TSSetTime(arg1::TS, arg2::PetscReal)::PetscErrorCode
end

function TSGetPrevTime(arg1, arg2)
    @ccall libpetsc.TSGetPrevTime(arg1::TS, arg2::Ptr{PetscReal})::PetscErrorCode
end

function TSGetTimeStep(arg1, arg2)
    @ccall libpetsc.TSGetTimeStep(arg1::TS, arg2::Ptr{PetscReal})::PetscErrorCode
end

function TSSetTimeStep(arg1, arg2)
    @ccall libpetsc.TSSetTimeStep(arg1::TS, arg2::PetscReal)::PetscErrorCode
end

function TSGetStepNumber(arg1, arg2)
    @ccall libpetsc.TSGetStepNumber(arg1::TS, arg2::Ptr{PetscInt})::PetscErrorCode
end

function TSSetStepNumber(arg1, arg2)
    @ccall libpetsc.TSSetStepNumber(arg1::TS, arg2::PetscInt)::PetscErrorCode
end

# typedef PetscErrorCode ( * TSRHSFunction ) ( TS , PetscReal , Vec , Vec , void * )
const TSRHSFunction = Ptr{Cvoid}

# typedef PetscErrorCode ( * TSRHSJacobian ) ( TS , PetscReal , Vec , Mat , Mat , void * )
const TSRHSJacobian = Ptr{Cvoid}

# typedef PetscErrorCode ( * TSRHSJacobianP ) ( TS , PetscReal , Vec , Mat , void * )
const TSRHSJacobianP = Ptr{Cvoid}

function TSSetRHSFunction(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSSetRHSFunction(arg1::TS, arg2::Vec, arg3::TSRHSFunction, arg4::Ptr{Cvoid})::PetscErrorCode
end

function TSGetRHSFunction(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSGetRHSFunction(arg1::TS, arg2::Ptr{Vec}, arg3::Ptr{TSRHSFunction}, arg4::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function TSSetRHSJacobian(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSSetRHSJacobian(arg1::TS, arg2::Mat, arg3::Mat, arg4::TSRHSJacobian, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TSGetRHSJacobian(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSGetRHSJacobian(arg1::TS, arg2::Ptr{Mat}, arg3::Ptr{Mat}, arg4::Ptr{TSRHSJacobian}, arg5::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function TSRHSJacobianSetReuse(arg1, arg2)
    @ccall libpetsc.TSRHSJacobianSetReuse(arg1::TS, arg2::PetscBool)::PetscErrorCode
end

# typedef PetscErrorCode ( * TSSolutionFunction ) ( TS , PetscReal , Vec , void * )
const TSSolutionFunction = Ptr{Cvoid}

function TSSetSolutionFunction(arg1, arg2, arg3)
    @ccall libpetsc.TSSetSolutionFunction(arg1::TS, arg2::TSSolutionFunction, arg3::Ptr{Cvoid})::PetscErrorCode
end

# typedef PetscErrorCode ( * TSForcingFunction ) ( TS , PetscReal , Vec , void * )
const TSForcingFunction = Ptr{Cvoid}

function TSSetForcingFunction(arg1, arg2, arg3)
    @ccall libpetsc.TSSetForcingFunction(arg1::TS, arg2::TSForcingFunction, arg3::Ptr{Cvoid})::PetscErrorCode
end

# typedef PetscErrorCode ( * TSIFunction ) ( TS , PetscReal , Vec , Vec , Vec , void * )
const TSIFunction = Ptr{Cvoid}

# typedef PetscErrorCode ( * TSIJacobian ) ( TS , PetscReal , Vec , Vec , PetscReal , Mat , Mat , void * )
const TSIJacobian = Ptr{Cvoid}

function TSSetIFunction(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSSetIFunction(arg1::TS, arg2::Vec, arg3::TSIFunction, arg4::Ptr{Cvoid})::PetscErrorCode
end

function TSGetIFunction(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSGetIFunction(arg1::TS, arg2::Ptr{Vec}, arg3::Ptr{TSIFunction}, arg4::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function TSSetIJacobian(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSSetIJacobian(arg1::TS, arg2::Mat, arg3::Mat, arg4::TSIJacobian, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TSGetIJacobian(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSGetIJacobian(arg1::TS, arg2::Ptr{Mat}, arg3::Ptr{Mat}, arg4::Ptr{TSIJacobian}, arg5::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

# typedef PetscErrorCode ( * TSI2Function ) ( TS , PetscReal , Vec , Vec , Vec , Vec , void * )
const TSI2Function = Ptr{Cvoid}

# typedef PetscErrorCode ( * TSI2Jacobian ) ( TS , PetscReal , Vec , Vec , Vec , PetscReal , PetscReal , Mat , Mat , void * )
const TSI2Jacobian = Ptr{Cvoid}

function TSSetI2Function(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSSetI2Function(arg1::TS, arg2::Vec, arg3::TSI2Function, arg4::Ptr{Cvoid})::PetscErrorCode
end

function TSGetI2Function(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSGetI2Function(arg1::TS, arg2::Ptr{Vec}, arg3::Ptr{TSI2Function}, arg4::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function TSSetI2Jacobian(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSSetI2Jacobian(arg1::TS, arg2::Mat, arg3::Mat, arg4::TSI2Jacobian, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TSGetI2Jacobian(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSGetI2Jacobian(arg1::TS, arg2::Ptr{Mat}, arg3::Ptr{Mat}, arg4::Ptr{TSI2Jacobian}, arg5::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function TSRHSSplitSetIS(arg1, arg2, arg3)
    @ccall libpetsc.TSRHSSplitSetIS(arg1::TS, arg2::Ptr{Cchar}, arg3::IS)::PetscErrorCode
end

function TSRHSSplitGetIS(arg1, arg2, arg3)
    @ccall libpetsc.TSRHSSplitGetIS(arg1::TS, arg2::Ptr{Cchar}, arg3::Ptr{IS})::PetscErrorCode
end

function TSRHSSplitSetRHSFunction(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSRHSSplitSetRHSFunction(arg1::TS, arg2::Ptr{Cchar}, arg3::Vec, arg4::TSRHSFunction, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TSRHSSplitGetSubTS(arg1, arg2, arg3)
    @ccall libpetsc.TSRHSSplitGetSubTS(arg1::TS, arg2::Ptr{Cchar}, arg3::Ptr{TS})::PetscErrorCode
end

function TSRHSSplitGetSubTSs(arg1, arg2, arg3)
    @ccall libpetsc.TSRHSSplitGetSubTSs(arg1::TS, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{TS}})::PetscErrorCode
end

function TSSetUseSplitRHSFunction(arg1, arg2)
    @ccall libpetsc.TSSetUseSplitRHSFunction(arg1::TS, arg2::PetscBool)::PetscErrorCode
end

function TSGetUseSplitRHSFunction(arg1, arg2)
    @ccall libpetsc.TSGetUseSplitRHSFunction(arg1::TS, arg2::Ptr{PetscBool})::PetscErrorCode
end

function TSComputeRHSFunctionLinear(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSComputeRHSFunctionLinear(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Vec, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TSComputeRHSJacobianConstant(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TSComputeRHSJacobianConstant(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Mat, arg5::Mat, arg6::Ptr{Cvoid})::PetscErrorCode
end

function TSComputeIFunctionLinear(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TSComputeIFunctionLinear(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Vec, arg5::Vec, arg6::Ptr{Cvoid})::PetscErrorCode
end

function TSComputeIJacobianConstant(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.TSComputeIJacobianConstant(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Vec, arg5::PetscReal, arg6::Mat, arg7::Mat, arg8::Ptr{Cvoid})::PetscErrorCode
end

function TSComputeSolutionFunction(arg1, arg2, arg3)
    @ccall libpetsc.TSComputeSolutionFunction(arg1::TS, arg2::PetscReal, arg3::Vec)::PetscErrorCode
end

function TSComputeForcingFunction(arg1, arg2, arg3)
    @ccall libpetsc.TSComputeForcingFunction(arg1::TS, arg2::PetscReal, arg3::Vec)::PetscErrorCode
end

function TSComputeIJacobianDefaultColor(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.TSComputeIJacobianDefaultColor(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Vec, arg5::PetscReal, arg6::Mat, arg7::Mat, arg8::Ptr{Cvoid})::PetscErrorCode
end

function TSSetPreStep(arg1, arg2)
    @ccall libpetsc.TSSetPreStep(arg1::TS, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TSSetPreStage(arg1, arg2)
    @ccall libpetsc.TSSetPreStage(arg1::TS, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TSSetPostStage(arg1, arg2)
    @ccall libpetsc.TSSetPostStage(arg1::TS, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TSSetPostEvaluate(arg1, arg2)
    @ccall libpetsc.TSSetPostEvaluate(arg1::TS, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TSSetPostStep(arg1, arg2)
    @ccall libpetsc.TSSetPostStep(arg1::TS, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TSPreStep(arg1)
    @ccall libpetsc.TSPreStep(arg1::TS)::PetscErrorCode
end

function TSPreStage(arg1, arg2)
    @ccall libpetsc.TSPreStage(arg1::TS, arg2::PetscReal)::PetscErrorCode
end

function TSPostStage(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSPostStage(arg1::TS, arg2::PetscReal, arg3::PetscInt, arg4::Ptr{Vec})::PetscErrorCode
end

function TSPostEvaluate(arg1)
    @ccall libpetsc.TSPostEvaluate(arg1::TS)::PetscErrorCode
end

function TSPostStep(arg1)
    @ccall libpetsc.TSPostStep(arg1::TS)::PetscErrorCode
end

function TSInterpolate(arg1, arg2, arg3)
    @ccall libpetsc.TSInterpolate(arg1::TS, arg2::PetscReal, arg3::Vec)::PetscErrorCode
end

function TSSetTolerances(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSSetTolerances(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::PetscReal, arg5::Vec)::PetscErrorCode
end

function TSGetTolerances(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSGetTolerances(arg1::TS, arg2::Ptr{PetscReal}, arg3::Ptr{Vec}, arg4::Ptr{PetscReal}, arg5::Ptr{Vec})::PetscErrorCode
end

function TSErrorWeightedNormInfinity(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TSErrorWeightedNormInfinity(arg1::TS, arg2::Vec, arg3::Vec, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal})::PetscErrorCode
end

function TSErrorWeightedNorm2(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TSErrorWeightedNorm2(arg1::TS, arg2::Vec, arg3::Vec, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal})::PetscErrorCode
end

function TSErrorWeightedNorm(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.TSErrorWeightedNorm(arg1::TS, arg2::Vec, arg3::Vec, arg4::NormType, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal}, arg7::Ptr{PetscReal})::PetscErrorCode
end

function TSErrorWeightedENormInfinity(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.TSErrorWeightedENormInfinity(arg1::TS, arg2::Vec, arg3::Vec, arg4::Vec, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal}, arg7::Ptr{PetscReal})::PetscErrorCode
end

function TSErrorWeightedENorm2(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.TSErrorWeightedENorm2(arg1::TS, arg2::Vec, arg3::Vec, arg4::Vec, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal}, arg7::Ptr{PetscReal})::PetscErrorCode
end

function TSErrorWeightedENorm(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.TSErrorWeightedENorm(arg1::TS, arg2::Vec, arg3::Vec, arg4::Vec, arg5::NormType, arg6::Ptr{PetscReal}, arg7::Ptr{PetscReal}, arg8::Ptr{PetscReal})::PetscErrorCode
end

function TSSetCFLTimeLocal(arg1, arg2)
    @ccall libpetsc.TSSetCFLTimeLocal(arg1::TS, arg2::PetscReal)::PetscErrorCode
end

function TSGetCFLTime(arg1, arg2)
    @ccall libpetsc.TSGetCFLTime(arg1::TS, arg2::Ptr{PetscReal})::PetscErrorCode
end

function TSSetFunctionDomainError(arg1, arg2)
    @ccall libpetsc.TSSetFunctionDomainError(arg1::TS, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TSFunctionDomainError(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSFunctionDomainError(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Ptr{PetscBool})::PetscErrorCode
end

function TSPseudoSetTimeStep(arg1, arg2, arg3)
    @ccall libpetsc.TSPseudoSetTimeStep(arg1::TS, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function TSPseudoTimeStepDefault(arg1, arg2, arg3)
    @ccall libpetsc.TSPseudoTimeStepDefault(arg1::TS, arg2::Ptr{PetscReal}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function TSPseudoComputeTimeStep(arg1, arg2)
    @ccall libpetsc.TSPseudoComputeTimeStep(arg1::TS, arg2::Ptr{PetscReal})::PetscErrorCode
end

function TSPseudoSetMaxTimeStep(arg1, arg2)
    @ccall libpetsc.TSPseudoSetMaxTimeStep(arg1::TS, arg2::PetscReal)::PetscErrorCode
end

function TSPseudoSetVerifyTimeStep(arg1, arg2, arg3)
    @ccall libpetsc.TSPseudoSetVerifyTimeStep(arg1::TS, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function TSPseudoVerifyTimeStepDefault(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSPseudoVerifyTimeStepDefault(arg1::TS, arg2::Vec, arg3::Ptr{Cvoid}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscBool})::PetscErrorCode
end

function TSPseudoVerifyTimeStep(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSPseudoVerifyTimeStep(arg1::TS, arg2::Vec, arg3::Ptr{PetscReal}, arg4::Ptr{PetscBool})::PetscErrorCode
end

function TSPseudoSetTimeStepIncrement(arg1, arg2)
    @ccall libpetsc.TSPseudoSetTimeStepIncrement(arg1::TS, arg2::PetscReal)::PetscErrorCode
end

function TSPseudoIncrementDtFromInitialDt(arg1)
    @ccall libpetsc.TSPseudoIncrementDtFromInitialDt(arg1::TS)::PetscErrorCode
end

function TSPythonSetType(arg1, arg2)
    @ccall libpetsc.TSPythonSetType(arg1::TS, arg2::Ptr{Cchar})::PetscErrorCode
end

function TSComputeRHSFunction(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSComputeRHSFunction(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Vec)::PetscErrorCode
end

function TSComputeRHSJacobian(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSComputeRHSJacobian(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Mat, arg5::Mat)::PetscErrorCode
end

function TSComputeIFunction(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TSComputeIFunction(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Vec, arg5::Vec, arg6::PetscBool)::PetscErrorCode
end

function TSComputeIJacobian(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libpetsc.TSComputeIJacobian(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Vec, arg5::PetscReal, arg6::Mat, arg7::Mat, arg8::PetscBool)::PetscErrorCode
end

function TSComputeI2Function(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TSComputeI2Function(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Vec, arg5::Vec, arg6::Vec)::PetscErrorCode
end

function TSComputeI2Jacobian(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.TSComputeI2Jacobian(arg1::TS, arg2::PetscReal, arg3::Vec, arg4::Vec, arg5::Vec, arg6::PetscReal, arg7::PetscReal, arg8::Mat, arg9::Mat)::PetscErrorCode
end

function TSComputeLinearStability(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSComputeLinearStability(arg1::TS, arg2::PetscReal, arg3::PetscReal, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal})::PetscErrorCode
end

function TSVISetVariableBounds(arg1, arg2, arg3)
    @ccall libpetsc.TSVISetVariableBounds(arg1::TS, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function DMTSSetBoundaryLocal(arg1, arg2, arg3)
    @ccall libpetsc.DMTSSetBoundaryLocal(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMTSSetRHSFunction(arg1, arg2, arg3)
    @ccall libpetsc.DMTSSetRHSFunction(arg1::DM, arg2::TSRHSFunction, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMTSGetRHSFunction(arg1, arg2, arg3)
    @ccall libpetsc.DMTSGetRHSFunction(arg1::DM, arg2::Ptr{TSRHSFunction}, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMTSSetRHSJacobian(arg1, arg2, arg3)
    @ccall libpetsc.DMTSSetRHSJacobian(arg1::DM, arg2::TSRHSJacobian, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMTSGetRHSJacobian(arg1, arg2, arg3)
    @ccall libpetsc.DMTSGetRHSJacobian(arg1::DM, arg2::Ptr{TSRHSJacobian}, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMTSSetIFunction(arg1, arg2, arg3)
    @ccall libpetsc.DMTSSetIFunction(arg1::DM, arg2::TSIFunction, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMTSGetIFunction(arg1, arg2, arg3)
    @ccall libpetsc.DMTSGetIFunction(arg1::DM, arg2::Ptr{TSIFunction}, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMTSSetIJacobian(arg1, arg2, arg3)
    @ccall libpetsc.DMTSSetIJacobian(arg1::DM, arg2::TSIJacobian, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMTSGetIJacobian(arg1, arg2, arg3)
    @ccall libpetsc.DMTSGetIJacobian(arg1::DM, arg2::Ptr{TSIJacobian}, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMTSSetI2Function(arg1, arg2, arg3)
    @ccall libpetsc.DMTSSetI2Function(arg1::DM, arg2::TSI2Function, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMTSGetI2Function(arg1, arg2, arg3)
    @ccall libpetsc.DMTSGetI2Function(arg1::DM, arg2::Ptr{TSI2Function}, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMTSSetI2Jacobian(arg1, arg2, arg3)
    @ccall libpetsc.DMTSSetI2Jacobian(arg1::DM, arg2::TSI2Jacobian, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMTSGetI2Jacobian(arg1, arg2, arg3)
    @ccall libpetsc.DMTSGetI2Jacobian(arg1::DM, arg2::Ptr{TSI2Jacobian}, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

# typedef PetscErrorCode ( * TSTransientVariable ) ( TS , Vec , Vec , void * )
const TSTransientVariable = Ptr{Cvoid}

function TSSetTransientVariable(arg1, arg2, arg3)
    @ccall libpetsc.TSSetTransientVariable(arg1::TS, arg2::TSTransientVariable, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMTSSetTransientVariable(arg1, arg2, arg3)
    @ccall libpetsc.DMTSSetTransientVariable(arg1::DM, arg2::TSTransientVariable, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMTSGetTransientVariable(arg1, arg2, arg3)
    @ccall libpetsc.DMTSGetTransientVariable(arg1::DM, arg2::Ptr{TSTransientVariable}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function TSComputeTransientVariable(arg1, arg2, arg3)
    @ccall libpetsc.TSComputeTransientVariable(arg1::TS, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function TSHasTransientVariable(arg1, arg2)
    @ccall libpetsc.TSHasTransientVariable(arg1::TS, arg2::Ptr{PetscBool})::PetscErrorCode
end

function DMTSSetSolutionFunction(arg1, arg2, arg3)
    @ccall libpetsc.DMTSSetSolutionFunction(arg1::DM, arg2::TSSolutionFunction, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMTSGetSolutionFunction(arg1, arg2, arg3)
    @ccall libpetsc.DMTSGetSolutionFunction(arg1::DM, arg2::Ptr{TSSolutionFunction}, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMTSSetForcingFunction(arg1, arg2, arg3)
    @ccall libpetsc.DMTSSetForcingFunction(arg1::DM, arg2::TSForcingFunction, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMTSGetForcingFunction(arg1, arg2, arg3)
    @ccall libpetsc.DMTSGetForcingFunction(arg1::DM, arg2::Ptr{TSForcingFunction}, arg3::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function DMTSCheckFromOptions(arg1, arg2)
    @ccall libpetsc.DMTSCheckFromOptions(arg1::TS, arg2::Vec)::PetscErrorCode
end

function DMTSSetIFunctionLocal(arg1, arg2, arg3)
    @ccall libpetsc.DMTSSetIFunctionLocal(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMTSSetIJacobianLocal(arg1, arg2, arg3)
    @ccall libpetsc.DMTSSetIJacobianLocal(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMTSSetRHSFunctionLocal(arg1, arg2, arg3)
    @ccall libpetsc.DMTSSetRHSFunctionLocal(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMTSSetIFunctionSerialize(arg1, arg2, arg3)
    @ccall libpetsc.DMTSSetIFunctionSerialize(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMTSSetIJacobianSerialize(arg1, arg2, arg3)
    @ccall libpetsc.DMTSSetIJacobianSerialize(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

# typedef PetscErrorCode ( * DMDATSRHSFunctionLocal ) ( DMDALocalInfo * , PetscReal , void * , void * , void * )
const DMDATSRHSFunctionLocal = Ptr{Cvoid}

# typedef PetscErrorCode ( * DMDATSRHSJacobianLocal ) ( DMDALocalInfo * , PetscReal , void * , Mat , Mat , void * )
const DMDATSRHSJacobianLocal = Ptr{Cvoid}

# typedef PetscErrorCode ( * DMDATSIFunctionLocal ) ( DMDALocalInfo * , PetscReal , void * , void * , void * , void * )
const DMDATSIFunctionLocal = Ptr{Cvoid}

# typedef PetscErrorCode ( * DMDATSIJacobianLocal ) ( DMDALocalInfo * , PetscReal , void * , void * , PetscReal , Mat , Mat , void * )
const DMDATSIJacobianLocal = Ptr{Cvoid}

function DMDATSSetRHSFunctionLocal(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMDATSSetRHSFunctionLocal(arg1::DM, arg2::InsertMode, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMDATSSetRHSJacobianLocal(arg1, arg2, arg3)
    @ccall libpetsc.DMDATSSetRHSJacobianLocal(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMDATSSetIFunctionLocal(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMDATSSetIFunctionLocal(arg1::DM, arg2::InsertMode, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function DMDATSSetIJacobianLocal(arg1, arg2, arg3)
    @ccall libpetsc.DMDATSSetIJacobianLocal(arg1::DM, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function DMPlexTSGetGeometryFVM(arg1, arg2, arg3, arg4)
    @ccall libpetsc.DMPlexTSGetGeometryFVM(arg1::DM, arg2::Ptr{Vec}, arg3::Ptr{Vec}, arg4::Ptr{PetscReal})::PetscErrorCode
end

mutable struct _n_TSMonitorLGCtx end

const TSMonitorLGCtx = Ptr{_n_TSMonitorLGCtx}

mutable struct TSMonitorDMDARayCtx
    ray::Vec
    scatter::VecScatter
    viewer::PetscViewer
    lgctx::TSMonitorLGCtx
    TSMonitorDMDARayCtx() = new()
end

function TSMonitorDMDARayDestroy(arg1)
    @ccall libpetsc.TSMonitorDMDARayDestroy(arg1::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function TSMonitorDMDARay(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSMonitorDMDARay(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TSMonitorLGDMDARay(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSMonitorLGDMDARay(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TSGetType(arg1, arg2)
    @ccall libpetsc.TSGetType(arg1::TS, arg2::Ptr{TSType})::PetscErrorCode
end

function TSSetType(arg1, arg2)
    @ccall libpetsc.TSSetType(arg1::TS, arg2::TSType)::PetscErrorCode
end

function TSRegister(arg1, arg2)
    @ccall libpetsc.TSRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TSGetSNES(arg1, arg2)
    @ccall libpetsc.TSGetSNES(arg1::TS, arg2::Ptr{SNES})::PetscErrorCode
end

function TSSetSNES(arg1, arg2)
    @ccall libpetsc.TSSetSNES(arg1::TS, arg2::SNES)::PetscErrorCode
end

function TSGetKSP(arg1, arg2)
    @ccall libpetsc.TSGetKSP(arg1::TS, arg2::Ptr{KSP})::PetscErrorCode
end

function TSView(arg1, arg2)
    @ccall libpetsc.TSView(arg1::TS, arg2::PetscViewer)::PetscErrorCode
end

function TSLoad(arg1, arg2)
    @ccall libpetsc.TSLoad(arg1::TS, arg2::PetscViewer)::PetscErrorCode
end

function TSViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.TSViewFromOptions(arg1::TS, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function TSTrajectoryViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.TSTrajectoryViewFromOptions(arg1::TSTrajectory, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function TSSetApplicationContext(arg1, arg2)
    @ccall libpetsc.TSSetApplicationContext(arg1::TS, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TSGetApplicationContext(arg1, arg2)
    @ccall libpetsc.TSGetApplicationContext(arg1::TS, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TSMonitorLGCtxCreate(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.TSMonitorLGCtxCreate(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Cint, arg5::Cint, arg6::Cint, arg7::Cint, arg8::PetscInt, arg9::Ptr{TSMonitorLGCtx})::PetscErrorCode
end

function TSMonitorLGCtxDestroy(arg1)
    @ccall libpetsc.TSMonitorLGCtxDestroy(arg1::Ptr{TSMonitorLGCtx})::PetscErrorCode
end

function TSMonitorLGTimeStep(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSMonitorLGTimeStep(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TSMonitorLGSolution(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSMonitorLGSolution(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TSMonitorLGSetVariableNames(arg1, arg2)
    @ccall libpetsc.TSMonitorLGSetVariableNames(arg1::TS, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function TSMonitorLGGetVariableNames(arg1, arg2)
    @ccall libpetsc.TSMonitorLGGetVariableNames(arg1::TS, arg2::Ptr{Ptr{Ptr{Cchar}}})::PetscErrorCode
end

function TSMonitorLGCtxSetVariableNames(arg1, arg2)
    @ccall libpetsc.TSMonitorLGCtxSetVariableNames(arg1::TSMonitorLGCtx, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function TSMonitorLGSetDisplayVariables(arg1, arg2)
    @ccall libpetsc.TSMonitorLGSetDisplayVariables(arg1::TS, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function TSMonitorLGCtxSetDisplayVariables(arg1, arg2)
    @ccall libpetsc.TSMonitorLGCtxSetDisplayVariables(arg1::TSMonitorLGCtx, arg2::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function TSMonitorLGSetTransform(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSMonitorLGSetTransform(arg1::TS, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function TSMonitorLGCtxSetTransform(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSMonitorLGCtxSetTransform(arg1::TSMonitorLGCtx, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function TSMonitorLGError(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSMonitorLGError(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TSMonitorLGSNESIterations(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSMonitorLGSNESIterations(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TSMonitorLGKSPIterations(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSMonitorLGKSPIterations(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TSMonitorError(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSMonitorError(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec, arg5::Ptr{PetscViewerAndFormat})::PetscErrorCode
end

mutable struct _n_TSMonitorLGCtxNetwork
    nlg::PetscInt
    lg::Ptr{PetscDrawLG}
    semilogy::PetscBool
    howoften::PetscInt
    _n_TSMonitorLGCtxNetwork() = new()
end

const TSMonitorLGCtxNetwork = Ptr{_n_TSMonitorLGCtxNetwork}

function TSMonitorLGCtxNetworkDestroy(arg1)
    @ccall libpetsc.TSMonitorLGCtxNetworkDestroy(arg1::Ptr{TSMonitorLGCtxNetwork})::PetscErrorCode
end

function TSMonitorLGCtxNetworkCreate(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.TSMonitorLGCtxNetworkCreate(arg1::TS, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Cint, arg5::Cint, arg6::Cint, arg7::Cint, arg8::PetscInt, arg9::Ptr{TSMonitorLGCtxNetwork})::PetscErrorCode
end

function TSMonitorLGCtxNetworkSolution(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSMonitorLGCtxNetworkSolution(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec, arg5::Ptr{Cvoid})::PetscErrorCode
end

mutable struct _n_TSMonitorEnvelopeCtx end

const TSMonitorEnvelopeCtx = Ptr{_n_TSMonitorEnvelopeCtx}

function TSMonitorEnvelopeCtxCreate(arg1, arg2)
    @ccall libpetsc.TSMonitorEnvelopeCtxCreate(arg1::TS, arg2::Ptr{TSMonitorEnvelopeCtx})::PetscErrorCode
end

function TSMonitorEnvelope(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSMonitorEnvelope(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TSMonitorEnvelopeGetBounds(arg1, arg2, arg3)
    @ccall libpetsc.TSMonitorEnvelopeGetBounds(arg1::TS, arg2::Ptr{Vec}, arg3::Ptr{Vec})::PetscErrorCode
end

function TSMonitorEnvelopeCtxDestroy(arg1)
    @ccall libpetsc.TSMonitorEnvelopeCtxDestroy(arg1::Ptr{TSMonitorEnvelopeCtx})::PetscErrorCode
end

mutable struct _n_TSMonitorSPEigCtx end

const TSMonitorSPEigCtx = Ptr{_n_TSMonitorSPEigCtx}

function TSMonitorSPEigCtxCreate(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.TSMonitorSPEigCtxCreate(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Cint, arg5::Cint, arg6::Cint, arg7::Cint, arg8::PetscInt, arg9::Ptr{TSMonitorSPEigCtx})::PetscErrorCode
end

function TSMonitorSPEigCtxDestroy(arg1)
    @ccall libpetsc.TSMonitorSPEigCtxDestroy(arg1::Ptr{TSMonitorSPEigCtx})::PetscErrorCode
end

function TSMonitorSPEig(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSMonitorSPEig(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec, arg5::Ptr{Cvoid})::PetscErrorCode
end

mutable struct _n_TSMonitorSPCtx end

const TSMonitorSPCtx = Ptr{_n_TSMonitorSPCtx}

function TSMonitorSPCtxCreate(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.TSMonitorSPCtxCreate(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Cint, arg5::Cint, arg6::Cint, arg7::Cint, arg8::PetscInt, arg9::Ptr{TSMonitorSPCtx})::PetscErrorCode
end

function TSMonitorSPCtxDestroy(arg1)
    @ccall libpetsc.TSMonitorSPCtxDestroy(arg1::Ptr{TSMonitorSPCtx})::PetscErrorCode
end

function TSMonitorSPSwarmSolution(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSMonitorSPSwarmSolution(arg1::TS, arg2::PetscInt, arg3::PetscReal, arg4::Vec, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TSSetEventHandler(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.TSSetEventHandler(arg1::TS, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscBool}, arg5::Ptr{Cvoid}, arg6::Ptr{Cvoid}, arg7::Ptr{Cvoid})::PetscErrorCode
end

function TSSetPostEventIntervalStep(arg1, arg2)
    @ccall libpetsc.TSSetPostEventIntervalStep(arg1::TS, arg2::PetscReal)::PetscErrorCode
end

function TSSetEventTolerances(arg1, arg2, arg3)
    @ccall libpetsc.TSSetEventTolerances(arg1::TS, arg2::PetscReal, arg3::Ptr{PetscReal})::PetscErrorCode
end

function TSGetNumEvents(arg1, arg2)
    @ccall libpetsc.TSGetNumEvents(arg1::TS, arg2::Ptr{PetscInt})::PetscErrorCode
end

const TSSSPType = Ptr{Cchar}

function TSSSPSetType(arg1, arg2)
    @ccall libpetsc.TSSSPSetType(arg1::TS, arg2::TSSSPType)::PetscErrorCode
end

function TSSSPGetType(arg1, arg2)
    @ccall libpetsc.TSSSPGetType(arg1::TS, arg2::Ptr{TSSSPType})::PetscErrorCode
end

function TSSSPSetNumStages(arg1, arg2)
    @ccall libpetsc.TSSSPSetNumStages(arg1::TS, arg2::PetscInt)::PetscErrorCode
end

function TSSSPGetNumStages(arg1, arg2)
    @ccall libpetsc.TSSSPGetNumStages(arg1::TS, arg2::Ptr{PetscInt})::PetscErrorCode
end

function TSSSPInitializePackage()
    @ccall libpetsc.TSSSPInitializePackage()::PetscErrorCode
end

function TSSSPFinalizePackage()
    @ccall libpetsc.TSSSPFinalizePackage()::PetscErrorCode
end

mutable struct _p_TSAdapt end

const TSAdapt = Ptr{_p_TSAdapt}

const TSAdaptType = Ptr{Cchar}

function TSGetAdapt(arg1, arg2)
    @ccall libpetsc.TSGetAdapt(arg1::TS, arg2::Ptr{TSAdapt})::PetscErrorCode
end

function TSAdaptRegister(arg1, arg2)
    @ccall libpetsc.TSAdaptRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TSAdaptInitializePackage()
    @ccall libpetsc.TSAdaptInitializePackage()::PetscErrorCode
end

function TSAdaptFinalizePackage()
    @ccall libpetsc.TSAdaptFinalizePackage()::PetscErrorCode
end

function TSAdaptCreate(arg1, arg2)
    @ccall libpetsc.TSAdaptCreate(arg1::MPI_Comm, arg2::Ptr{TSAdapt})::PetscErrorCode
end

function TSAdaptSetType(arg1, arg2)
    @ccall libpetsc.TSAdaptSetType(arg1::TSAdapt, arg2::TSAdaptType)::PetscErrorCode
end

function TSAdaptGetType(arg1, arg2)
    @ccall libpetsc.TSAdaptGetType(arg1::TSAdapt, arg2::Ptr{TSAdaptType})::PetscErrorCode
end

function TSAdaptSetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.TSAdaptSetOptionsPrefix(arg1::TSAdapt, arg2::Ptr{Cchar})::PetscErrorCode
end

function TSAdaptCandidatesClear(arg1)
    @ccall libpetsc.TSAdaptCandidatesClear(arg1::TSAdapt)::PetscErrorCode
end

function TSAdaptCandidateAdd(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.TSAdaptCandidateAdd(arg1::TSAdapt, arg2::Ptr{Cchar}, arg3::PetscInt, arg4::PetscInt, arg5::PetscReal, arg6::PetscReal, arg7::PetscBool)::PetscErrorCode
end

function TSAdaptCandidatesGet(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TSAdaptCandidatesGet(arg1::TSAdapt, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{PetscInt}}, arg4::Ptr{Ptr{PetscInt}}, arg5::Ptr{Ptr{PetscReal}}, arg6::Ptr{Ptr{PetscReal}})::PetscErrorCode
end

function TSAdaptChoose(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TSAdaptChoose(arg1::TSAdapt, arg2::TS, arg3::PetscReal, arg4::Ptr{PetscInt}, arg5::Ptr{PetscReal}, arg6::Ptr{PetscBool})::PetscErrorCode
end

function TSAdaptCheckStage(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSAdaptCheckStage(arg1::TSAdapt, arg2::TS, arg3::PetscReal, arg4::Vec, arg5::Ptr{PetscBool})::PetscErrorCode
end

function TSAdaptView(arg1, arg2)
    @ccall libpetsc.TSAdaptView(arg1::TSAdapt, arg2::PetscViewer)::PetscErrorCode
end

function TSAdaptLoad(arg1, arg2)
    @ccall libpetsc.TSAdaptLoad(arg1::TSAdapt, arg2::PetscViewer)::PetscErrorCode
end

function TSAdaptSetFromOptions(arg1, arg2)
    @ccall libpetsc.TSAdaptSetFromOptions(arg1::Ptr{PetscOptionItems}, arg2::TSAdapt)::PetscErrorCode
end

function TSAdaptReset(arg1)
    @ccall libpetsc.TSAdaptReset(arg1::TSAdapt)::PetscErrorCode
end

function TSAdaptDestroy(arg1)
    @ccall libpetsc.TSAdaptDestroy(arg1::Ptr{TSAdapt})::PetscErrorCode
end

function TSAdaptSetMonitor(arg1, arg2)
    @ccall libpetsc.TSAdaptSetMonitor(arg1::TSAdapt, arg2::PetscBool)::PetscErrorCode
end

function TSAdaptSetAlwaysAccept(arg1, arg2)
    @ccall libpetsc.TSAdaptSetAlwaysAccept(arg1::TSAdapt, arg2::PetscBool)::PetscErrorCode
end

function TSAdaptSetSafety(arg1, arg2, arg3)
    @ccall libpetsc.TSAdaptSetSafety(arg1::TSAdapt, arg2::PetscReal, arg3::PetscReal)::PetscErrorCode
end

function TSAdaptGetSafety(arg1, arg2, arg3)
    @ccall libpetsc.TSAdaptGetSafety(arg1::TSAdapt, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal})::PetscErrorCode
end

function TSAdaptSetMaxIgnore(arg1, arg2)
    @ccall libpetsc.TSAdaptSetMaxIgnore(arg1::TSAdapt, arg2::PetscReal)::PetscErrorCode
end

function TSAdaptGetMaxIgnore(arg1, arg2)
    @ccall libpetsc.TSAdaptGetMaxIgnore(arg1::TSAdapt, arg2::Ptr{PetscReal})::PetscErrorCode
end

function TSAdaptSetClip(arg1, arg2, arg3)
    @ccall libpetsc.TSAdaptSetClip(arg1::TSAdapt, arg2::PetscReal, arg3::PetscReal)::PetscErrorCode
end

function TSAdaptGetClip(arg1, arg2, arg3)
    @ccall libpetsc.TSAdaptGetClip(arg1::TSAdapt, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal})::PetscErrorCode
end

function TSAdaptSetScaleSolveFailed(arg1, arg2)
    @ccall libpetsc.TSAdaptSetScaleSolveFailed(arg1::TSAdapt, arg2::PetscReal)::PetscErrorCode
end

function TSAdaptGetScaleSolveFailed(arg1, arg2)
    @ccall libpetsc.TSAdaptGetScaleSolveFailed(arg1::TSAdapt, arg2::Ptr{PetscReal})::PetscErrorCode
end

function TSAdaptSetStepLimits(arg1, arg2, arg3)
    @ccall libpetsc.TSAdaptSetStepLimits(arg1::TSAdapt, arg2::PetscReal, arg3::PetscReal)::PetscErrorCode
end

function TSAdaptGetStepLimits(arg1, arg2, arg3)
    @ccall libpetsc.TSAdaptGetStepLimits(arg1::TSAdapt, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal})::PetscErrorCode
end

function TSAdaptSetCheckStage(arg1, arg2)
    @ccall libpetsc.TSAdaptSetCheckStage(arg1::TSAdapt, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TSAdaptHistorySetHistory(arg1, n, hist, arg4)
    @ccall libpetsc.TSAdaptHistorySetHistory(arg1::TSAdapt, n::PetscInt, hist::Ptr{PetscReal}, arg4::PetscBool)::PetscErrorCode
end

function TSAdaptHistorySetTrajectory(arg1, arg2, arg3)
    @ccall libpetsc.TSAdaptHistorySetTrajectory(arg1::TSAdapt, arg2::TSTrajectory, arg3::PetscBool)::PetscErrorCode
end

function TSAdaptHistoryGetStep(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSAdaptHistoryGetStep(arg1::TSAdapt, arg2::PetscInt, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal})::PetscErrorCode
end

function TSAdaptSetTimeStepIncreaseDelay(arg1, arg2)
    @ccall libpetsc.TSAdaptSetTimeStepIncreaseDelay(arg1::TSAdapt, arg2::PetscInt)::PetscErrorCode
end

function TSAdaptDSPSetFilter(arg1, arg2)
    @ccall libpetsc.TSAdaptDSPSetFilter(arg1::TSAdapt, arg2::Ptr{Cchar})::PetscErrorCode
end

function TSAdaptDSPSetPID(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSAdaptDSPSetPID(arg1::TSAdapt, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal)::PetscErrorCode
end

mutable struct _p_TSGLLEAdapt end

const TSGLLEAdapt = Ptr{_p_TSGLLEAdapt}

const TSGLLEAdaptType = Ptr{Cchar}

function TSGLLEAdaptRegister(arg1, arg2)
    @ccall libpetsc.TSGLLEAdaptRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TSGLLEAdaptInitializePackage()
    @ccall libpetsc.TSGLLEAdaptInitializePackage()::PetscErrorCode
end

function TSGLLEAdaptFinalizePackage()
    @ccall libpetsc.TSGLLEAdaptFinalizePackage()::PetscErrorCode
end

function TSGLLEAdaptCreate(arg1, arg2)
    @ccall libpetsc.TSGLLEAdaptCreate(arg1::MPI_Comm, arg2::Ptr{TSGLLEAdapt})::PetscErrorCode
end

function TSGLLEAdaptSetType(arg1, arg2)
    @ccall libpetsc.TSGLLEAdaptSetType(arg1::TSGLLEAdapt, arg2::TSGLLEAdaptType)::PetscErrorCode
end

function TSGLLEAdaptSetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.TSGLLEAdaptSetOptionsPrefix(arg1::TSGLLEAdapt, arg2::Ptr{Cchar})::PetscErrorCode
end

function TSGLLEAdaptChoose(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11)
    @ccall libpetsc.TSGLLEAdaptChoose(arg1::TSGLLEAdapt, arg2::PetscInt, arg3::Ptr{PetscInt}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal}, arg6::PetscInt, arg7::PetscReal, arg8::PetscReal, arg9::Ptr{PetscInt}, arg10::Ptr{PetscReal}, arg11::Ptr{PetscBool})::PetscErrorCode
end

function TSGLLEAdaptView(arg1, arg2)
    @ccall libpetsc.TSGLLEAdaptView(arg1::TSGLLEAdapt, arg2::PetscViewer)::PetscErrorCode
end

function TSGLLEAdaptSetFromOptions(arg1, arg2)
    @ccall libpetsc.TSGLLEAdaptSetFromOptions(arg1::Ptr{PetscOptionItems}, arg2::TSGLLEAdapt)::PetscErrorCode
end

function TSGLLEAdaptDestroy(arg1)
    @ccall libpetsc.TSGLLEAdaptDestroy(arg1::Ptr{TSGLLEAdapt})::PetscErrorCode
end

const TSGLLEAcceptType = Ptr{Cchar}

# typedef PetscErrorCode ( * TSGLLEAcceptFunction ) ( TS , PetscReal , PetscReal , const PetscReal [ ] , PetscBool * )
const TSGLLEAcceptFunction = Ptr{Cvoid}

function TSGLLEAcceptRegister(arg1, arg2)
    @ccall libpetsc.TSGLLEAcceptRegister(arg1::Ptr{Cchar}, arg2::TSGLLEAcceptFunction)::PetscErrorCode
end

const TSGLLEType = Ptr{Cchar}

function TSGLLERegister(arg1, arg2)
    @ccall libpetsc.TSGLLERegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TSGLLEInitializePackage()
    @ccall libpetsc.TSGLLEInitializePackage()::PetscErrorCode
end

function TSGLLEFinalizePackage()
    @ccall libpetsc.TSGLLEFinalizePackage()::PetscErrorCode
end

function TSGLLESetType(arg1, arg2)
    @ccall libpetsc.TSGLLESetType(arg1::TS, arg2::TSGLLEType)::PetscErrorCode
end

function TSGLLEGetAdapt(arg1, arg2)
    @ccall libpetsc.TSGLLEGetAdapt(arg1::TS, arg2::Ptr{TSGLLEAdapt})::PetscErrorCode
end

function TSGLLESetAcceptType(arg1, arg2)
    @ccall libpetsc.TSGLLESetAcceptType(arg1::TS, arg2::TSGLLEAcceptType)::PetscErrorCode
end

function TSEIMEXSetMaxRows(ts, arg2)
    @ccall libpetsc.TSEIMEXSetMaxRows(ts::TS, arg2::PetscInt)::PetscErrorCode
end

function TSEIMEXSetRowCol(ts, arg2, arg3)
    @ccall libpetsc.TSEIMEXSetRowCol(ts::TS, arg2::PetscInt, arg3::PetscInt)::PetscErrorCode
end

function TSEIMEXSetOrdAdapt(arg1, arg2)
    @ccall libpetsc.TSEIMEXSetOrdAdapt(arg1::TS, arg2::PetscBool)::PetscErrorCode
end

const TSRKType = Ptr{Cchar}

function TSRKGetOrder(arg1, arg2)
    @ccall libpetsc.TSRKGetOrder(arg1::TS, arg2::Ptr{PetscInt})::PetscErrorCode
end

function TSRKGetType(arg1, arg2)
    @ccall libpetsc.TSRKGetType(arg1::TS, arg2::Ptr{TSRKType})::PetscErrorCode
end

function TSRKSetType(arg1, arg2)
    @ccall libpetsc.TSRKSetType(arg1::TS, arg2::TSRKType)::PetscErrorCode
end

function TSRKGetTableau(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.TSRKGetTableau(arg1::TS, arg2::Ptr{PetscInt}, arg3::Ptr{Ptr{PetscReal}}, arg4::Ptr{Ptr{PetscReal}}, arg5::Ptr{Ptr{PetscReal}}, arg6::Ptr{Ptr{PetscReal}}, arg7::Ptr{PetscInt}, arg8::Ptr{Ptr{PetscReal}}, arg9::Ptr{PetscBool})::PetscErrorCode
end

function TSRKSetMultirate(arg1, arg2)
    @ccall libpetsc.TSRKSetMultirate(arg1::TS, arg2::PetscBool)::PetscErrorCode
end

function TSRKGetMultirate(arg1, arg2)
    @ccall libpetsc.TSRKGetMultirate(arg1::TS, arg2::Ptr{PetscBool})::PetscErrorCode
end

function TSRKRegister(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.TSRKRegister(arg1::TSRKType, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal}, arg7::Ptr{PetscReal}, arg8::PetscInt, arg9::Ptr{PetscReal})::PetscErrorCode
end

function TSRKInitializePackage()
    @ccall libpetsc.TSRKInitializePackage()::PetscErrorCode
end

function TSRKFinalizePackage()
    @ccall libpetsc.TSRKFinalizePackage()::PetscErrorCode
end

function TSRKRegisterDestroy()
    @ccall libpetsc.TSRKRegisterDestroy()::PetscErrorCode
end

const TSMPRKType = Ptr{Cchar}

function TSMPRKGetType(ts, arg2)
    @ccall libpetsc.TSMPRKGetType(ts::TS, arg2::Ptr{TSMPRKType})::PetscErrorCode
end

function TSMPRKSetType(ts, arg2)
    @ccall libpetsc.TSMPRKSetType(ts::TS, arg2::TSMPRKType)::PetscErrorCode
end

function TSMPRKRegister(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13, arg14, arg15, arg16)
    @ccall libpetsc.TSMPRKRegister(arg1::TSMPRKType, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscInt, arg6::Ptr{PetscReal}, arg7::Ptr{PetscReal}, arg8::Ptr{PetscReal}, arg9::Ptr{PetscInt}, arg10::Ptr{PetscReal}, arg11::Ptr{PetscReal}, arg12::Ptr{PetscReal}, arg13::Ptr{PetscInt}, arg14::Ptr{PetscReal}, arg15::Ptr{PetscReal}, arg16::Ptr{PetscReal})::PetscErrorCode
end

function TSMPRKInitializePackage()
    @ccall libpetsc.TSMPRKInitializePackage()::PetscErrorCode
end

function TSMPRKFinalizePackage()
    @ccall libpetsc.TSMPRKFinalizePackage()::PetscErrorCode
end

function TSMPRKRegisterDestroy()
    @ccall libpetsc.TSMPRKRegisterDestroy()::PetscErrorCode
end

const TSGLEEType = Ptr{Cchar}

function TSGLEEGetType(ts, arg2)
    @ccall libpetsc.TSGLEEGetType(ts::TS, arg2::Ptr{TSGLEEType})::PetscErrorCode
end

function TSGLEESetType(ts, arg2)
    @ccall libpetsc.TSGLEESetType(ts::TS, arg2::TSGLEEType)::PetscErrorCode
end

function TSGLEERegister(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13, arg14, arg15, arg16, arg17)
    @ccall libpetsc.TSGLEERegister(arg1::TSGLEEType, arg2::PetscInt, arg3::PetscInt, arg4::PetscInt, arg5::PetscReal, arg6::Ptr{PetscReal}, arg7::Ptr{PetscReal}, arg8::Ptr{PetscReal}, arg9::Ptr{PetscReal}, arg10::Ptr{PetscReal}, arg11::Ptr{PetscReal}, arg12::Ptr{PetscReal}, arg13::Ptr{PetscReal}, arg14::Ptr{PetscReal}, arg15::Ptr{PetscReal}, arg16::PetscInt, arg17::Ptr{PetscReal})::PetscErrorCode
end

function TSGLEEFinalizePackage()
    @ccall libpetsc.TSGLEEFinalizePackage()::PetscErrorCode
end

function TSGLEEInitializePackage()
    @ccall libpetsc.TSGLEEInitializePackage()::PetscErrorCode
end

function TSGLEERegisterDestroy()
    @ccall libpetsc.TSGLEERegisterDestroy()::PetscErrorCode
end

const TSARKIMEXType = Ptr{Cchar}

function TSARKIMEXGetType(ts, arg2)
    @ccall libpetsc.TSARKIMEXGetType(ts::TS, arg2::Ptr{TSARKIMEXType})::PetscErrorCode
end

function TSARKIMEXSetType(ts, arg2)
    @ccall libpetsc.TSARKIMEXSetType(ts::TS, arg2::TSARKIMEXType)::PetscErrorCode
end

function TSARKIMEXSetFullyImplicit(arg1, arg2)
    @ccall libpetsc.TSARKIMEXSetFullyImplicit(arg1::TS, arg2::PetscBool)::PetscErrorCode
end

function TSARKIMEXGetFullyImplicit(arg1, arg2)
    @ccall libpetsc.TSARKIMEXGetFullyImplicit(arg1::TS, arg2::Ptr{PetscBool})::PetscErrorCode
end

function TSARKIMEXRegister(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13, arg14)
    @ccall libpetsc.TSARKIMEXRegister(arg1::TSARKIMEXType, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal}, arg7::Ptr{PetscReal}, arg8::Ptr{PetscReal}, arg9::Ptr{PetscReal}, arg10::Ptr{PetscReal}, arg11::Ptr{PetscReal}, arg12::PetscInt, arg13::Ptr{PetscReal}, arg14::Ptr{PetscReal})::PetscErrorCode
end

function TSARKIMEXInitializePackage()
    @ccall libpetsc.TSARKIMEXInitializePackage()::PetscErrorCode
end

function TSARKIMEXFinalizePackage()
    @ccall libpetsc.TSARKIMEXFinalizePackage()::PetscErrorCode
end

function TSARKIMEXRegisterDestroy()
    @ccall libpetsc.TSARKIMEXRegisterDestroy()::PetscErrorCode
end

const TSRosWType = Ptr{Cchar}

function TSRosWGetType(arg1, arg2)
    @ccall libpetsc.TSRosWGetType(arg1::TS, arg2::Ptr{TSRosWType})::PetscErrorCode
end

function TSRosWSetType(arg1, arg2)
    @ccall libpetsc.TSRosWSetType(arg1::TS, arg2::TSRosWType)::PetscErrorCode
end

function TSRosWSetRecomputeJacobian(arg1, arg2)
    @ccall libpetsc.TSRosWSetRecomputeJacobian(arg1::TS, arg2::PetscBool)::PetscErrorCode
end

function TSRosWRegister(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.TSRosWRegister(arg1::TSRosWType, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal}, arg7::Ptr{PetscReal}, arg8::PetscInt, arg9::Ptr{PetscReal})::PetscErrorCode
end

function TSRosWRegisterRos4(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TSRosWRegisterRos4(arg1::TSRosWType, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal, arg6::PetscReal)::PetscErrorCode
end

function TSRosWInitializePackage()
    @ccall libpetsc.TSRosWInitializePackage()::PetscErrorCode
end

function TSRosWFinalizePackage()
    @ccall libpetsc.TSRosWFinalizePackage()::PetscErrorCode
end

function TSRosWRegisterDestroy()
    @ccall libpetsc.TSRosWRegisterDestroy()::PetscErrorCode
end

function TSBDFSetOrder(arg1, arg2)
    @ccall libpetsc.TSBDFSetOrder(arg1::TS, arg2::PetscInt)::PetscErrorCode
end

function TSBDFGetOrder(arg1, arg2)
    @ccall libpetsc.TSBDFGetOrder(arg1::TS, arg2::Ptr{PetscInt})::PetscErrorCode
end

const TSBasicSymplecticType = Ptr{Cchar}

function TSBasicSymplecticSetType(arg1, arg2)
    @ccall libpetsc.TSBasicSymplecticSetType(arg1::TS, arg2::TSBasicSymplecticType)::PetscErrorCode
end

function TSBasicSymplecticGetType(arg1, arg2)
    @ccall libpetsc.TSBasicSymplecticGetType(arg1::TS, arg2::Ptr{TSBasicSymplecticType})::PetscErrorCode
end

function TSBasicSymplecticRegister(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSBasicSymplecticRegister(arg1::TSBasicSymplecticType, arg2::PetscInt, arg3::PetscInt, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal})::PetscErrorCode
end

function TSBasicSymplecticInitializePackage()
    @ccall libpetsc.TSBasicSymplecticInitializePackage()::PetscErrorCode
end

function TSBasicSymplecticFinalizePackage()
    @ccall libpetsc.TSBasicSymplecticFinalizePackage()::PetscErrorCode
end

function TSBasicSymplecticRegisterDestroy()
    @ccall libpetsc.TSBasicSymplecticRegisterDestroy()::PetscErrorCode
end

function TSDiscGradSetFormulation(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSDiscGradSetFormulation(arg1::TS, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TSThetaSetTheta(arg1, arg2)
    @ccall libpetsc.TSThetaSetTheta(arg1::TS, arg2::PetscReal)::PetscErrorCode
end

function TSThetaGetTheta(arg1, arg2)
    @ccall libpetsc.TSThetaGetTheta(arg1::TS, arg2::Ptr{PetscReal})::PetscErrorCode
end

function TSThetaGetEndpoint(arg1, arg2)
    @ccall libpetsc.TSThetaGetEndpoint(arg1::TS, arg2::Ptr{PetscBool})::PetscErrorCode
end

function TSThetaSetEndpoint(arg1, arg2)
    @ccall libpetsc.TSThetaSetEndpoint(arg1::TS, arg2::PetscBool)::PetscErrorCode
end

function TSAlphaSetRadius(arg1, arg2)
    @ccall libpetsc.TSAlphaSetRadius(arg1::TS, arg2::PetscReal)::PetscErrorCode
end

function TSAlphaSetParams(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSAlphaSetParams(arg1::TS, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal)::PetscErrorCode
end

function TSAlphaGetParams(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TSAlphaGetParams(arg1::TS, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal})::PetscErrorCode
end

function TSAlpha2SetRadius(arg1, arg2)
    @ccall libpetsc.TSAlpha2SetRadius(arg1::TS, arg2::PetscReal)::PetscErrorCode
end

function TSAlpha2SetParams(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSAlpha2SetParams(arg1::TS, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal)::PetscErrorCode
end

function TSAlpha2GetParams(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TSAlpha2GetParams(arg1::TS, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal})::PetscErrorCode
end

function TSSetDM(arg1, arg2)
    @ccall libpetsc.TSSetDM(arg1::TS, arg2::DM)::PetscErrorCode
end

function TSGetDM(arg1, arg2)
    @ccall libpetsc.TSGetDM(arg1::TS, arg2::Ptr{DM})::PetscErrorCode
end

function SNESTSFormFunction(arg1, arg2, arg3, arg4)
    @ccall libpetsc.SNESTSFormFunction(arg1::SNES, arg2::Vec, arg3::Vec, arg4::Ptr{Cvoid})::PetscErrorCode
end

function SNESTSFormJacobian(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.SNESTSFormJacobian(arg1::SNES, arg2::Vec, arg3::Mat, arg4::Mat, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TSRHSJacobianTest(arg1, arg2)
    @ccall libpetsc.TSRHSJacobianTest(arg1::TS, arg2::Ptr{PetscBool})::PetscErrorCode
end

function TSRHSJacobianTestTranspose(arg1, arg2)
    @ccall libpetsc.TSRHSJacobianTestTranspose(arg1::TS, arg2::Ptr{PetscBool})::PetscErrorCode
end

function TSGetComputeInitialCondition(arg1, arg2)
    @ccall libpetsc.TSGetComputeInitialCondition(arg1::TS, arg2::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function TSSetComputeInitialCondition(arg1, arg2)
    @ccall libpetsc.TSSetComputeInitialCondition(arg1::TS, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TSComputeInitialCondition(arg1, arg2)
    @ccall libpetsc.TSComputeInitialCondition(arg1::TS, arg2::Vec)::PetscErrorCode
end

function TSGetComputeExactError(arg1, arg2)
    @ccall libpetsc.TSGetComputeExactError(arg1::TS, arg2::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function TSSetComputeExactError(arg1, arg2)
    @ccall libpetsc.TSSetComputeExactError(arg1::TS, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TSComputeExactError(arg1, arg2, arg3)
    @ccall libpetsc.TSComputeExactError(arg1::TS, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function PetscConvEstUseTS(arg1, arg2)
    @ccall libpetsc.PetscConvEstUseTS(arg1::PetscConvEst, arg2::PetscBool)::PetscErrorCode
end

function TSSetMatStructure(arg1, arg2)
    @ccall libpetsc.TSSetMatStructure(arg1::TS, arg2::MatStructure)::PetscErrorCode
end

function VecFischer(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.VecFischer(arg1::Vec, arg2::Vec, arg3::Vec, arg4::Vec, arg5::Vec)::PetscErrorCode
end

function VecSFischer(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.VecSFischer(arg1::Vec, arg2::Vec, arg3::Vec, arg4::Vec, arg5::PetscReal, arg6::Vec)::PetscErrorCode
end

function MatDFischer(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.MatDFischer(arg1::Mat, arg2::Vec, arg3::Vec, arg4::Vec, arg5::Vec, arg6::Vec, arg7::Vec, arg8::Vec, arg9::Vec)::PetscErrorCode
end

function MatDSFischer(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11)
    @ccall libpetsc.MatDSFischer(arg1::Mat, arg2::Vec, arg3::Vec, arg4::Vec, arg5::Vec, arg6::PetscReal, arg7::Vec, arg8::Vec, arg9::Vec, arg10::Vec, arg11::Vec)::PetscErrorCode
end

function TaoSoftThreshold(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TaoSoftThreshold(arg1::Vec, arg2::PetscReal, arg3::PetscReal, arg4::Vec)::PetscErrorCode
end

@enum TaoSubsetType::UInt32 begin
    TAO_SUBSET_SUBVEC = 0
    TAO_SUBSET_MASK = 1
    TAO_SUBSET_MATRIXFREE = 2
end

@enum TaoADMMUpdateType::UInt32 begin
    TAO_ADMM_UPDATE_BASIC = 0
    TAO_ADMM_UPDATE_ADAPTIVE = 1
    TAO_ADMM_UPDATE_ADAPTIVE_RELAXED = 2
end

@enum TaoADMMRegularizerType::UInt32 begin
    TAO_ADMM_REGULARIZER_USER = 0
    TAO_ADMM_REGULARIZER_SOFT_THRESH = 1
end

@enum TaoALMMType::UInt32 begin
    TAO_ALMM_CLASSIC = 0
    TAO_ALMM_PHR = 1
end

mutable struct _p_Tao end

const Tao = Ptr{_p_Tao}

const TaoType = Ptr{Cchar}

@enum TaoConvergedReason::Int32 begin
    TAO_CONVERGED_GATOL = 3
    TAO_CONVERGED_GRTOL = 4
    TAO_CONVERGED_GTTOL = 5
    TAO_CONVERGED_STEPTOL = 6
    TAO_CONVERGED_MINF = 7
    TAO_CONVERGED_USER = 8
    TAO_DIVERGED_MAXITS = -2
    TAO_DIVERGED_NAN = -4
    TAO_DIVERGED_MAXFCN = -5
    TAO_DIVERGED_LS_FAILURE = -6
    TAO_DIVERGED_TR_REDUCTION = -7
    TAO_DIVERGED_USER = -8
    TAO_CONTINUE_ITERATING = 0
end

function TaoInitializePackage()
    @ccall libpetsc.TaoInitializePackage()::PetscErrorCode
end

function TaoFinalizePackage()
    @ccall libpetsc.TaoFinalizePackage()::PetscErrorCode
end

function TaoCreate(arg1, arg2)
    @ccall libpetsc.TaoCreate(arg1::MPI_Comm, arg2::Ptr{Tao})::PetscErrorCode
end

function TaoSetFromOptions(arg1)
    @ccall libpetsc.TaoSetFromOptions(arg1::Tao)::PetscErrorCode
end

function TaoSetUp(arg1)
    @ccall libpetsc.TaoSetUp(arg1::Tao)::PetscErrorCode
end

function TaoSetType(arg1, arg2)
    @ccall libpetsc.TaoSetType(arg1::Tao, arg2::TaoType)::PetscErrorCode
end

function TaoGetType(arg1, arg2)
    @ccall libpetsc.TaoGetType(arg1::Tao, arg2::Ptr{TaoType})::PetscErrorCode
end

function TaoSetApplicationContext(arg1, arg2)
    @ccall libpetsc.TaoSetApplicationContext(arg1::Tao, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TaoGetApplicationContext(arg1, arg2)
    @ccall libpetsc.TaoGetApplicationContext(arg1::Tao, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TaoDestroy(arg1)
    @ccall libpetsc.TaoDestroy(arg1::Ptr{Tao})::PetscErrorCode
end

function TaoSetOptionsPrefix(arg1, arg2)
    @ccall libpetsc.TaoSetOptionsPrefix(arg1::Tao, arg2::Ptr{Cchar})::PetscErrorCode
end

function TaoView(arg1, arg2)
    @ccall libpetsc.TaoView(arg1::Tao, arg2::PetscViewer)::PetscErrorCode
end

function TaoViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.TaoViewFromOptions(arg1::Tao, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function TaoSolve(arg1)
    @ccall libpetsc.TaoSolve(arg1::Tao)::PetscErrorCode
end

function TaoRegister(arg1, arg2)
    @ccall libpetsc.TaoRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TaoRegisterDestroy()
    @ccall libpetsc.TaoRegisterDestroy()::PetscErrorCode
end

function TaoGetConvergedReason(arg1, arg2)
    @ccall libpetsc.TaoGetConvergedReason(arg1::Tao, arg2::Ptr{TaoConvergedReason})::PetscErrorCode
end

function TaoGetSolutionStatus(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.TaoGetSolutionStatus(arg1::Tao, arg2::Ptr{PetscInt}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscReal}, arg6::Ptr{PetscReal}, arg7::Ptr{TaoConvergedReason})::PetscErrorCode
end

function TaoSetConvergedReason(arg1, arg2)
    @ccall libpetsc.TaoSetConvergedReason(arg1::Tao, arg2::TaoConvergedReason)::PetscErrorCode
end

function TaoSetInitialVector(arg1, arg2)
    @ccall libpetsc.TaoSetInitialVector(arg1::Tao, arg2::Vec)::PetscErrorCode
end

function TaoGetSolutionVector(arg1, arg2)
    @ccall libpetsc.TaoGetSolutionVector(arg1::Tao, arg2::Ptr{Vec})::PetscErrorCode
end

function TaoGetGradientVector(arg1, arg2)
    @ccall libpetsc.TaoGetGradientVector(arg1::Tao, arg2::Ptr{Vec})::PetscErrorCode
end

function TaoSetGradientNorm(arg1, arg2)
    @ccall libpetsc.TaoSetGradientNorm(arg1::Tao, arg2::Mat)::PetscErrorCode
end

function TaoGetGradientNorm(arg1, arg2)
    @ccall libpetsc.TaoGetGradientNorm(arg1::Tao, arg2::Ptr{Mat})::PetscErrorCode
end

function TaoSetLMVMMatrix(arg1, arg2)
    @ccall libpetsc.TaoSetLMVMMatrix(arg1::Tao, arg2::Mat)::PetscErrorCode
end

function TaoGetLMVMMatrix(arg1, arg2)
    @ccall libpetsc.TaoGetLMVMMatrix(arg1::Tao, arg2::Ptr{Mat})::PetscErrorCode
end

function TaoSetRecycleHistory(arg1, arg2)
    @ccall libpetsc.TaoSetRecycleHistory(arg1::Tao, arg2::PetscBool)::PetscErrorCode
end

function TaoGetRecycleHistory(arg1, arg2)
    @ccall libpetsc.TaoGetRecycleHistory(arg1::Tao, arg2::Ptr{PetscBool})::PetscErrorCode
end

function TaoLMVMSetH0(arg1, arg2)
    @ccall libpetsc.TaoLMVMSetH0(arg1::Tao, arg2::Mat)::PetscErrorCode
end

function TaoLMVMGetH0(arg1, arg2)
    @ccall libpetsc.TaoLMVMGetH0(arg1::Tao, arg2::Ptr{Mat})::PetscErrorCode
end

function TaoLMVMGetH0KSP(arg1, arg2)
    @ccall libpetsc.TaoLMVMGetH0KSP(arg1::Tao, arg2::Ptr{KSP})::PetscErrorCode
end

function TaoLMVMRecycle(arg1, arg2)
    @ccall libpetsc.TaoLMVMRecycle(arg1::Tao, arg2::PetscBool)::PetscErrorCode
end

function TaoSetObjectiveRoutine(arg1, arg2, arg3)
    @ccall libpetsc.TaoSetObjectiveRoutine(arg1::Tao, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function TaoSetGradientRoutine(arg1, arg2, arg3)
    @ccall libpetsc.TaoSetGradientRoutine(arg1::Tao, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function TaoSetObjectiveAndGradientRoutine(arg1, arg2, arg3)
    @ccall libpetsc.TaoSetObjectiveAndGradientRoutine(arg1::Tao, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function TaoSetHessianRoutine(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TaoSetHessianRoutine(arg1::Tao, arg2::Mat, arg3::Mat, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TaoSetResidualRoutine(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TaoSetResidualRoutine(arg1::Tao, arg2::Vec, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function TaoSetResidualWeights(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TaoSetResidualWeights(arg1::Tao, arg2::Vec, arg3::PetscInt, arg4::Ptr{PetscInt}, arg5::Ptr{PetscInt}, arg6::Ptr{PetscReal})::PetscErrorCode
end

function TaoSetConstraintsRoutine(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TaoSetConstraintsRoutine(arg1::Tao, arg2::Vec, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function TaoSetInequalityConstraintsRoutine(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TaoSetInequalityConstraintsRoutine(arg1::Tao, arg2::Vec, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function TaoSetEqualityConstraintsRoutine(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TaoSetEqualityConstraintsRoutine(arg1::Tao, arg2::Vec, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function TaoSetJacobianResidualRoutine(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TaoSetJacobianResidualRoutine(arg1::Tao, arg2::Mat, arg3::Mat, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TaoSetJacobianRoutine(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TaoSetJacobianRoutine(arg1::Tao, arg2::Mat, arg3::Mat, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TaoSetJacobianStateRoutine(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TaoSetJacobianStateRoutine(arg1::Tao, arg2::Mat, arg3::Mat, arg4::Mat, arg5::Ptr{Cvoid}, arg6::Ptr{Cvoid})::PetscErrorCode
end

function TaoSetJacobianDesignRoutine(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TaoSetJacobianDesignRoutine(arg1::Tao, arg2::Mat, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function TaoSetJacobianInequalityRoutine(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TaoSetJacobianInequalityRoutine(arg1::Tao, arg2::Mat, arg3::Mat, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TaoSetJacobianEqualityRoutine(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TaoSetJacobianEqualityRoutine(arg1::Tao, arg2::Mat, arg3::Mat, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TaoShellSetSolve(arg1, arg2)
    @ccall libpetsc.TaoShellSetSolve(arg1::Tao, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TaoShellSetContext(arg1, arg2)
    @ccall libpetsc.TaoShellSetContext(arg1::Tao, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TaoShellGetContext(arg1, arg2)
    @ccall libpetsc.TaoShellGetContext(arg1::Tao, arg2::Ptr{Ptr{Cvoid}})::PetscErrorCode
end

function TaoSetSeparableObjectiveRoutine(tao, res, func, ctx)
    @ccall libpetsc.TaoSetSeparableObjectiveRoutine(tao::Tao, res::Vec, func::Ptr{Cvoid}, ctx::Ptr{Cvoid})::PetscErrorCode
end

function TaoSetSeparableObjectiveWeights(tao, sigma_v, n, rows, cols, vals)
    @ccall libpetsc.TaoSetSeparableObjectiveWeights(tao::Tao, sigma_v::Vec, n::PetscInt, rows::Ptr{PetscInt}, cols::Ptr{PetscInt}, vals::Ptr{PetscReal})::PetscErrorCode
end

function TaoSetStateDesignIS(arg1, arg2, arg3)
    @ccall libpetsc.TaoSetStateDesignIS(arg1::Tao, arg2::IS, arg3::IS)::PetscErrorCode
end

function TaoComputeObjective(arg1, arg2, arg3)
    @ccall libpetsc.TaoComputeObjective(arg1::Tao, arg2::Vec, arg3::Ptr{PetscReal})::PetscErrorCode
end

function TaoComputeResidual(arg1, arg2, arg3)
    @ccall libpetsc.TaoComputeResidual(arg1::Tao, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function TaoTestGradient(arg1, arg2, arg3)
    @ccall libpetsc.TaoTestGradient(arg1::Tao, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function TaoComputeGradient(arg1, arg2, arg3)
    @ccall libpetsc.TaoComputeGradient(arg1::Tao, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function TaoComputeObjectiveAndGradient(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TaoComputeObjectiveAndGradient(arg1::Tao, arg2::Vec, arg3::Ptr{PetscReal}, arg4::Vec)::PetscErrorCode
end

function TaoComputeConstraints(arg1, arg2, arg3)
    @ccall libpetsc.TaoComputeConstraints(arg1::Tao, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function TaoComputeInequalityConstraints(arg1, arg2, arg3)
    @ccall libpetsc.TaoComputeInequalityConstraints(arg1::Tao, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function TaoComputeEqualityConstraints(arg1, arg2, arg3)
    @ccall libpetsc.TaoComputeEqualityConstraints(arg1::Tao, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function TaoDefaultComputeGradient(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TaoDefaultComputeGradient(arg1::Tao, arg2::Vec, arg3::Vec, arg4::Ptr{Cvoid})::PetscErrorCode
end

function TaoIsObjectiveDefined(arg1, arg2)
    @ccall libpetsc.TaoIsObjectiveDefined(arg1::Tao, arg2::Ptr{PetscBool})::PetscErrorCode
end

function TaoIsGradientDefined(arg1, arg2)
    @ccall libpetsc.TaoIsGradientDefined(arg1::Tao, arg2::Ptr{PetscBool})::PetscErrorCode
end

function TaoIsObjectiveAndGradientDefined(arg1, arg2)
    @ccall libpetsc.TaoIsObjectiveAndGradientDefined(arg1::Tao, arg2::Ptr{PetscBool})::PetscErrorCode
end

function TaoComputeSeparableObjective(tao, X, F)
    @ccall libpetsc.TaoComputeSeparableObjective(tao::Tao, X::Vec, F::Vec)::PetscErrorCode
end

function TaoTestHessian(arg1)
    @ccall libpetsc.TaoTestHessian(arg1::Tao)::PetscErrorCode
end

function TaoComputeHessian(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TaoComputeHessian(arg1::Tao, arg2::Vec, arg3::Mat, arg4::Mat)::PetscErrorCode
end

function TaoComputeResidualJacobian(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TaoComputeResidualJacobian(arg1::Tao, arg2::Vec, arg3::Mat, arg4::Mat)::PetscErrorCode
end

function TaoComputeJacobian(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TaoComputeJacobian(arg1::Tao, arg2::Vec, arg3::Mat, arg4::Mat)::PetscErrorCode
end

function TaoComputeJacobianState(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TaoComputeJacobianState(arg1::Tao, arg2::Vec, arg3::Mat, arg4::Mat, arg5::Mat)::PetscErrorCode
end

function TaoComputeJacobianEquality(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TaoComputeJacobianEquality(arg1::Tao, arg2::Vec, arg3::Mat, arg4::Mat)::PetscErrorCode
end

function TaoComputeJacobianInequality(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TaoComputeJacobianInequality(arg1::Tao, arg2::Vec, arg3::Mat, arg4::Mat)::PetscErrorCode
end

function TaoComputeJacobianDesign(arg1, arg2, arg3)
    @ccall libpetsc.TaoComputeJacobianDesign(arg1::Tao, arg2::Vec, arg3::Mat)::PetscErrorCode
end

function TaoDefaultComputeHessian(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TaoDefaultComputeHessian(arg1::Tao, arg2::Vec, arg3::Mat, arg4::Mat, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TaoDefaultComputeHessianColor(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TaoDefaultComputeHessianColor(arg1::Tao, arg2::Vec, arg3::Mat, arg4::Mat, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TaoDefaultComputeHessianMFFD(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TaoDefaultComputeHessianMFFD(arg1::Tao, arg2::Vec, arg3::Mat, arg4::Mat, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TaoComputeDualVariables(arg1, arg2, arg3)
    @ccall libpetsc.TaoComputeDualVariables(arg1::Tao, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function TaoSetVariableBounds(arg1, arg2, arg3)
    @ccall libpetsc.TaoSetVariableBounds(arg1::Tao, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function TaoGetVariableBounds(arg1, arg2, arg3)
    @ccall libpetsc.TaoGetVariableBounds(arg1::Tao, arg2::Ptr{Vec}, arg3::Ptr{Vec})::PetscErrorCode
end

function TaoGetDualVariables(arg1, arg2, arg3)
    @ccall libpetsc.TaoGetDualVariables(arg1::Tao, arg2::Ptr{Vec}, arg3::Ptr{Vec})::PetscErrorCode
end

function TaoSetInequalityBounds(arg1, arg2, arg3)
    @ccall libpetsc.TaoSetInequalityBounds(arg1::Tao, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function TaoGetInequalityBounds(arg1, arg2, arg3)
    @ccall libpetsc.TaoGetInequalityBounds(arg1::Tao, arg2::Ptr{Vec}, arg3::Ptr{Vec})::PetscErrorCode
end

function TaoSetVariableBoundsRoutine(arg1, arg2, arg3)
    @ccall libpetsc.TaoSetVariableBoundsRoutine(arg1::Tao, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function TaoComputeVariableBounds(arg1)
    @ccall libpetsc.TaoComputeVariableBounds(arg1::Tao)::PetscErrorCode
end

function TaoGetTolerances(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TaoGetTolerances(arg1::Tao, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal})::PetscErrorCode
end

function TaoSetTolerances(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TaoSetTolerances(arg1::Tao, arg2::PetscReal, arg3::PetscReal, arg4::PetscReal)::PetscErrorCode
end

function TaoGetConstraintTolerances(arg1, arg2, arg3)
    @ccall libpetsc.TaoGetConstraintTolerances(arg1::Tao, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal})::PetscErrorCode
end

function TaoSetConstraintTolerances(arg1, arg2, arg3)
    @ccall libpetsc.TaoSetConstraintTolerances(arg1::Tao, arg2::PetscReal, arg3::PetscReal)::PetscErrorCode
end

function TaoSetFunctionLowerBound(arg1, arg2)
    @ccall libpetsc.TaoSetFunctionLowerBound(arg1::Tao, arg2::PetscReal)::PetscErrorCode
end

function TaoSetInitialTrustRegionRadius(arg1, arg2)
    @ccall libpetsc.TaoSetInitialTrustRegionRadius(arg1::Tao, arg2::PetscReal)::PetscErrorCode
end

function TaoSetMaximumIterations(arg1, arg2)
    @ccall libpetsc.TaoSetMaximumIterations(arg1::Tao, arg2::PetscInt)::PetscErrorCode
end

function TaoSetMaximumFunctionEvaluations(arg1, arg2)
    @ccall libpetsc.TaoSetMaximumFunctionEvaluations(arg1::Tao, arg2::PetscInt)::PetscErrorCode
end

function TaoGetFunctionLowerBound(arg1, arg2)
    @ccall libpetsc.TaoGetFunctionLowerBound(arg1::Tao, arg2::Ptr{PetscReal})::PetscErrorCode
end

function TaoGetInitialTrustRegionRadius(arg1, arg2)
    @ccall libpetsc.TaoGetInitialTrustRegionRadius(arg1::Tao, arg2::Ptr{PetscReal})::PetscErrorCode
end

function TaoGetCurrentTrustRegionRadius(arg1, arg2)
    @ccall libpetsc.TaoGetCurrentTrustRegionRadius(arg1::Tao, arg2::Ptr{PetscReal})::PetscErrorCode
end

function TaoGetMaximumIterations(arg1, arg2)
    @ccall libpetsc.TaoGetMaximumIterations(arg1::Tao, arg2::Ptr{PetscInt})::PetscErrorCode
end

function TaoGetCurrentFunctionEvaluations(arg1, arg2)
    @ccall libpetsc.TaoGetCurrentFunctionEvaluations(arg1::Tao, arg2::Ptr{PetscInt})::PetscErrorCode
end

function TaoGetMaximumFunctionEvaluations(arg1, arg2)
    @ccall libpetsc.TaoGetMaximumFunctionEvaluations(arg1::Tao, arg2::Ptr{PetscInt})::PetscErrorCode
end

function TaoGetIterationNumber(arg1, arg2)
    @ccall libpetsc.TaoGetIterationNumber(arg1::Tao, arg2::Ptr{PetscInt})::PetscErrorCode
end

function TaoSetIterationNumber(arg1, arg2)
    @ccall libpetsc.TaoSetIterationNumber(arg1::Tao, arg2::PetscInt)::PetscErrorCode
end

function TaoGetTotalIterationNumber(arg1, arg2)
    @ccall libpetsc.TaoGetTotalIterationNumber(arg1::Tao, arg2::Ptr{PetscInt})::PetscErrorCode
end

function TaoSetTotalIterationNumber(arg1, arg2)
    @ccall libpetsc.TaoSetTotalIterationNumber(arg1::Tao, arg2::PetscInt)::PetscErrorCode
end

function TaoGetResidualNorm(arg1, arg2)
    @ccall libpetsc.TaoGetResidualNorm(arg1::Tao, arg2::Ptr{PetscReal})::PetscErrorCode
end

function TaoGetObjective(arg1, arg2)
    @ccall libpetsc.TaoGetObjective(arg1::Tao, arg2::Ptr{PetscReal})::PetscErrorCode
end

function TaoAppendOptionsPrefix(arg1, p)
    @ccall libpetsc.TaoAppendOptionsPrefix(arg1::Tao, p::Ptr{Cchar})::PetscErrorCode
end

function TaoGetOptionsPrefix(arg1, p)
    @ccall libpetsc.TaoGetOptionsPrefix(arg1::Tao, p::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function TaoResetStatistics(arg1)
    @ccall libpetsc.TaoResetStatistics(arg1::Tao)::PetscErrorCode
end

function TaoSetUpdate(arg1, arg2, arg3)
    @ccall libpetsc.TaoSetUpdate(arg1::Tao, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function TaoGetKSP(arg1, arg2)
    @ccall libpetsc.TaoGetKSP(arg1::Tao, arg2::Ptr{KSP})::PetscErrorCode
end

function TaoGetLinearSolveIterations(arg1, arg2)
    @ccall libpetsc.TaoGetLinearSolveIterations(arg1::Tao, arg2::Ptr{PetscInt})::PetscErrorCode
end

mutable struct _p_TaoLineSearch end

const TaoLineSearch = Ptr{_p_TaoLineSearch}

@enum TaoLineSearchConvergedReason::Int32 begin
    TAOLINESEARCH_FAILED_INFORNAN = -1
    TAOLINESEARCH_FAILED_BADPARAMETER = -2
    TAOLINESEARCH_FAILED_ASCENT = -3
    TAOLINESEARCH_CONTINUE_ITERATING = 0
    TAOLINESEARCH_SUCCESS = 1
    TAOLINESEARCH_SUCCESS_USER = 2
    TAOLINESEARCH_HALTED_OTHER = 3
    TAOLINESEARCH_HALTED_MAXFCN = 4
    TAOLINESEARCH_HALTED_UPPERBOUND = 5
    TAOLINESEARCH_HALTED_LOWERBOUND = 6
    TAOLINESEARCH_HALTED_RTOL = 7
    TAOLINESEARCH_HALTED_USER = 8
end

const TaoLineSearchType = Ptr{Cchar}

function TaoLineSearchCreate(arg1, arg2)
    @ccall libpetsc.TaoLineSearchCreate(arg1::MPI_Comm, arg2::Ptr{TaoLineSearch})::PetscErrorCode
end

function TaoLineSearchSetFromOptions(arg1)
    @ccall libpetsc.TaoLineSearchSetFromOptions(arg1::TaoLineSearch)::PetscErrorCode
end

function TaoLineSearchSetUp(arg1)
    @ccall libpetsc.TaoLineSearchSetUp(arg1::TaoLineSearch)::PetscErrorCode
end

function TaoLineSearchDestroy(arg1)
    @ccall libpetsc.TaoLineSearchDestroy(arg1::Ptr{TaoLineSearch})::PetscErrorCode
end

function TaoLineSearchMonitor(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TaoLineSearchMonitor(arg1::TaoLineSearch, arg2::PetscInt, arg3::PetscReal, arg4::PetscReal)::PetscErrorCode
end

function TaoLineSearchView(arg1, arg2)
    @ccall libpetsc.TaoLineSearchView(arg1::TaoLineSearch, arg2::PetscViewer)::PetscErrorCode
end

function TaoLineSearchViewFromOptions(arg1, arg2, arg3)
    @ccall libpetsc.TaoLineSearchViewFromOptions(arg1::TaoLineSearch, arg2::PetscObject, arg3::Ptr{Cchar})::PetscErrorCode
end

function TaoLineSearchSetOptionsPrefix(arg1, prefix)
    @ccall libpetsc.TaoLineSearchSetOptionsPrefix(arg1::TaoLineSearch, prefix::Ptr{Cchar})::PetscErrorCode
end

function TaoLineSearchReset(arg1)
    @ccall libpetsc.TaoLineSearchReset(arg1::TaoLineSearch)::PetscErrorCode
end

function TaoLineSearchAppendOptionsPrefix(arg1, prefix)
    @ccall libpetsc.TaoLineSearchAppendOptionsPrefix(arg1::TaoLineSearch, prefix::Ptr{Cchar})::PetscErrorCode
end

function TaoLineSearchGetOptionsPrefix(arg1, prefix)
    @ccall libpetsc.TaoLineSearchGetOptionsPrefix(arg1::TaoLineSearch, prefix::Ptr{Ptr{Cchar}})::PetscErrorCode
end

function TaoLineSearchApply(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.TaoLineSearchApply(arg1::TaoLineSearch, arg2::Vec, arg3::Ptr{PetscReal}, arg4::Vec, arg5::Vec, arg6::Ptr{PetscReal}, arg7::Ptr{TaoLineSearchConvergedReason})::PetscErrorCode
end

function TaoLineSearchGetStepLength(arg1, arg2)
    @ccall libpetsc.TaoLineSearchGetStepLength(arg1::TaoLineSearch, arg2::Ptr{PetscReal})::PetscErrorCode
end

function TaoLineSearchGetStartingVector(arg1, arg2)
    @ccall libpetsc.TaoLineSearchGetStartingVector(arg1::TaoLineSearch, arg2::Ptr{Vec})::PetscErrorCode
end

function TaoLineSearchGetStepDirection(arg1, arg2)
    @ccall libpetsc.TaoLineSearchGetStepDirection(arg1::TaoLineSearch, arg2::Ptr{Vec})::PetscErrorCode
end

function TaoLineSearchSetInitialStepLength(arg1, arg2)
    @ccall libpetsc.TaoLineSearchSetInitialStepLength(arg1::TaoLineSearch, arg2::PetscReal)::PetscErrorCode
end

function TaoLineSearchGetSolution(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TaoLineSearchGetSolution(arg1::TaoLineSearch, arg2::Vec, arg3::Ptr{PetscReal}, arg4::Vec, arg5::Ptr{PetscReal}, arg6::Ptr{TaoLineSearchConvergedReason})::PetscErrorCode
end

function TaoLineSearchGetFullStepObjective(arg1, arg2)
    @ccall libpetsc.TaoLineSearchGetFullStepObjective(arg1::TaoLineSearch, arg2::Ptr{PetscReal})::PetscErrorCode
end

function TaoLineSearchGetNumberFunctionEvaluations(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TaoLineSearchGetNumberFunctionEvaluations(arg1::TaoLineSearch, arg2::Ptr{PetscInt}, arg3::Ptr{PetscInt}, arg4::Ptr{PetscInt})::PetscErrorCode
end

function TaoLineSearchGetType(arg1, arg2)
    @ccall libpetsc.TaoLineSearchGetType(arg1::TaoLineSearch, arg2::Ptr{TaoLineSearchType})::PetscErrorCode
end

function TaoLineSearchSetType(arg1, arg2)
    @ccall libpetsc.TaoLineSearchSetType(arg1::TaoLineSearch, arg2::TaoLineSearchType)::PetscErrorCode
end

function TaoLineSearchIsUsingTaoRoutines(arg1, arg2)
    @ccall libpetsc.TaoLineSearchIsUsingTaoRoutines(arg1::TaoLineSearch, arg2::Ptr{PetscBool})::PetscErrorCode
end

function TaoLineSearchSetObjectiveAndGTSRoutine(arg1, arg2, arg3)
    @ccall libpetsc.TaoLineSearchSetObjectiveAndGTSRoutine(arg1::TaoLineSearch, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function TaoLineSearchSetObjectiveRoutine(arg1, arg2, arg3)
    @ccall libpetsc.TaoLineSearchSetObjectiveRoutine(arg1::TaoLineSearch, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function TaoLineSearchSetGradientRoutine(arg1, arg2, arg3)
    @ccall libpetsc.TaoLineSearchSetGradientRoutine(arg1::TaoLineSearch, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function TaoLineSearchSetObjectiveAndGradientRoutine(arg1, arg2, arg3)
    @ccall libpetsc.TaoLineSearchSetObjectiveAndGradientRoutine(arg1::TaoLineSearch, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function TaoLineSearchComputeObjective(arg1, arg2, arg3)
    @ccall libpetsc.TaoLineSearchComputeObjective(arg1::TaoLineSearch, arg2::Vec, arg3::Ptr{PetscReal})::PetscErrorCode
end

function TaoLineSearchComputeGradient(arg1, arg2, arg3)
    @ccall libpetsc.TaoLineSearchComputeGradient(arg1::TaoLineSearch, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function TaoLineSearchComputeObjectiveAndGradient(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TaoLineSearchComputeObjectiveAndGradient(arg1::TaoLineSearch, arg2::Vec, arg3::Ptr{PetscReal}, arg4::Vec)::PetscErrorCode
end

function TaoLineSearchComputeObjectiveAndGTS(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TaoLineSearchComputeObjectiveAndGTS(arg1::TaoLineSearch, arg2::Vec, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal})::PetscErrorCode
end

function TaoLineSearchSetVariableBounds(arg1, arg2, arg3)
    @ccall libpetsc.TaoLineSearchSetVariableBounds(arg1::TaoLineSearch, arg2::Vec, arg3::Vec)::PetscErrorCode
end

function TaoLineSearchInitializePackage()
    @ccall libpetsc.TaoLineSearchInitializePackage()::PetscErrorCode
end

function TaoLineSearchFinalizePackage()
    @ccall libpetsc.TaoLineSearchFinalizePackage()::PetscErrorCode
end

function TaoLineSearchRegister(arg1, arg2)
    @ccall libpetsc.TaoLineSearchRegister(arg1::Ptr{Cchar}, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TaoLineSearchUseTaoRoutines(arg1, arg2)
    @ccall libpetsc.TaoLineSearchUseTaoRoutines(arg1::TaoLineSearch, arg2::Tao)::PetscErrorCode
end

function TaoGetLineSearch(arg1, arg2)
    @ccall libpetsc.TaoGetLineSearch(arg1::Tao, arg2::Ptr{TaoLineSearch})::PetscErrorCode
end

function TaoSetConvergenceHistory(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libpetsc.TaoSetConvergenceHistory(arg1::Tao, arg2::Ptr{PetscReal}, arg3::Ptr{PetscReal}, arg4::Ptr{PetscReal}, arg5::Ptr{PetscInt}, arg6::PetscInt, arg7::PetscBool)::PetscErrorCode
end

function TaoGetConvergenceHistory(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TaoGetConvergenceHistory(arg1::Tao, arg2::Ptr{Ptr{PetscReal}}, arg3::Ptr{Ptr{PetscReal}}, arg4::Ptr{Ptr{PetscReal}}, arg5::Ptr{Ptr{PetscInt}}, arg6::Ptr{PetscInt})::PetscErrorCode
end

function TaoSetMonitor(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TaoSetMonitor(arg1::Tao, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function TaoCancelMonitors(arg1)
    @ccall libpetsc.TaoCancelMonitors(arg1::Tao)::PetscErrorCode
end

function TaoMonitorDefault(arg1, arg2)
    @ccall libpetsc.TaoMonitorDefault(arg1::Tao, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TaoDefaultMonitor(tao, ctx)
    @ccall libpetsc.TaoDefaultMonitor(tao::Tao, ctx::Ptr{Cvoid})::PetscErrorCode
end

function TaoDefaultGMonitor(arg1, arg2)
    @ccall libpetsc.TaoDefaultGMonitor(arg1::Tao, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TaoDefaultSMonitor(arg1, arg2)
    @ccall libpetsc.TaoDefaultSMonitor(arg1::Tao, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TaoDefaultCMonitor(arg1, arg2)
    @ccall libpetsc.TaoDefaultCMonitor(arg1::Tao, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TaoSolutionMonitor(arg1, arg2)
    @ccall libpetsc.TaoSolutionMonitor(arg1::Tao, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TaoResidualMonitor(arg1, arg2)
    @ccall libpetsc.TaoResidualMonitor(arg1::Tao, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TaoGradientMonitor(arg1, arg2)
    @ccall libpetsc.TaoGradientMonitor(arg1::Tao, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TaoStepDirectionMonitor(arg1, arg2)
    @ccall libpetsc.TaoStepDirectionMonitor(arg1::Tao, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TaoDrawSolutionMonitor(arg1, arg2)
    @ccall libpetsc.TaoDrawSolutionMonitor(arg1::Tao, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TaoDrawStepMonitor(arg1, arg2)
    @ccall libpetsc.TaoDrawStepMonitor(arg1::Tao, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TaoDrawGradientMonitor(arg1, arg2)
    @ccall libpetsc.TaoDrawGradientMonitor(arg1::Tao, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TaoAddLineSearchCounts(arg1)
    @ccall libpetsc.TaoAddLineSearchCounts(arg1::Tao)::PetscErrorCode
end

function TaoDefaultConvergenceTest(arg1, arg2)
    @ccall libpetsc.TaoDefaultConvergenceTest(arg1::Tao, arg2::Ptr{Cvoid})::PetscErrorCode
end

function TaoSetConvergenceTest(arg1, arg2, arg3)
    @ccall libpetsc.TaoSetConvergenceTest(arg1::Tao, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function TaoLCLSetStateDesignIS(arg1, arg2, arg3)
    @ccall libpetsc.TaoLCLSetStateDesignIS(arg1::Tao, arg2::IS, arg3::IS)::PetscErrorCode
end

function TaoMonitor(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libpetsc.TaoMonitor(arg1::Tao, arg2::PetscInt, arg3::PetscReal, arg4::PetscReal, arg5::PetscReal, arg6::PetscReal)::PetscErrorCode
end

mutable struct _n_TaoMonitorDrawCtx end

const TaoMonitorDrawCtx = Ptr{_n_TaoMonitorDrawCtx}

function TaoMonitorDrawCtxCreate(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    @ccall libpetsc.TaoMonitorDrawCtxCreate(arg1::MPI_Comm, arg2::Ptr{Cchar}, arg3::Ptr{Cchar}, arg4::Cint, arg5::Cint, arg6::Cint, arg7::Cint, arg8::PetscInt, arg9::Ptr{TaoMonitorDrawCtx})::PetscErrorCode
end

function TaoMonitorDrawCtxDestroy(arg1)
    @ccall libpetsc.TaoMonitorDrawCtxDestroy(arg1::Ptr{TaoMonitorDrawCtx})::PetscErrorCode
end

function TaoBRGNGetSubsolver(arg1, arg2)
    @ccall libpetsc.TaoBRGNGetSubsolver(arg1::Tao, arg2::Ptr{Tao})::PetscErrorCode
end

function TaoBRGNSetRegularizerObjectiveAndGradientRoutine(arg1, arg2, arg3)
    @ccall libpetsc.TaoBRGNSetRegularizerObjectiveAndGradientRoutine(arg1::Tao, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function TaoBRGNSetRegularizerHessianRoutine(arg1, arg2, arg3, arg4)
    @ccall libpetsc.TaoBRGNSetRegularizerHessianRoutine(arg1::Tao, arg2::Mat, arg3::Ptr{Cvoid}, arg4::Ptr{Cvoid})::PetscErrorCode
end

function TaoBRGNSetRegularizerWeight(arg1, arg2)
    @ccall libpetsc.TaoBRGNSetRegularizerWeight(arg1::Tao, arg2::PetscReal)::PetscErrorCode
end

function TaoBRGNSetL1SmoothEpsilon(arg1, arg2)
    @ccall libpetsc.TaoBRGNSetL1SmoothEpsilon(arg1::Tao, arg2::PetscReal)::PetscErrorCode
end

function TaoBRGNSetDictionaryMatrix(arg1, arg2)
    @ccall libpetsc.TaoBRGNSetDictionaryMatrix(arg1::Tao, arg2::Mat)::PetscErrorCode
end

function TaoBRGNGetDampingVector(arg1, arg2)
    @ccall libpetsc.TaoBRGNGetDampingVector(arg1::Tao, arg2::Ptr{Vec})::PetscErrorCode
end

function TaoADMMGetMisfitSubsolver(arg1, arg2)
    @ccall libpetsc.TaoADMMGetMisfitSubsolver(arg1::Tao, arg2::Ptr{Tao})::PetscErrorCode
end

function TaoADMMGetRegularizationSubsolver(arg1, arg2)
    @ccall libpetsc.TaoADMMGetRegularizationSubsolver(arg1::Tao, arg2::Ptr{Tao})::PetscErrorCode
end

function TaoADMMGetDualVector(arg1, arg2)
    @ccall libpetsc.TaoADMMGetDualVector(arg1::Tao, arg2::Ptr{Vec})::PetscErrorCode
end

function TaoADMMGetSpectralPenalty(arg1, arg2)
    @ccall libpetsc.TaoADMMGetSpectralPenalty(arg1::Tao, arg2::Ptr{PetscReal})::PetscErrorCode
end

function TaoADMMSetSpectralPenalty(arg1, arg2)
    @ccall libpetsc.TaoADMMSetSpectralPenalty(arg1::Tao, arg2::PetscReal)::PetscErrorCode
end

function TaoGetADMMParentTao(arg1, arg2)
    @ccall libpetsc.TaoGetADMMParentTao(arg1::Tao, arg2::Ptr{Tao})::PetscErrorCode
end

function TaoADMMSetConstraintVectorRHS(arg1, arg2)
    @ccall libpetsc.TaoADMMSetConstraintVectorRHS(arg1::Tao, arg2::Vec)::PetscErrorCode
end

function TaoADMMSetRegularizerCoefficient(arg1, arg2)
    @ccall libpetsc.TaoADMMSetRegularizerCoefficient(arg1::Tao, arg2::PetscReal)::PetscErrorCode
end

function TaoADMMSetMisfitConstraintJacobian(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TaoADMMSetMisfitConstraintJacobian(arg1::Tao, arg2::Mat, arg3::Mat, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TaoADMMSetRegularizerConstraintJacobian(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TaoADMMSetRegularizerConstraintJacobian(arg1::Tao, arg2::Mat, arg3::Mat, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TaoADMMSetRegularizerHessianRoutine(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TaoADMMSetRegularizerHessianRoutine(arg1::Tao, arg2::Mat, arg3::Mat, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TaoADMMSetRegularizerObjectiveAndGradientRoutine(arg1, arg2, arg3)
    @ccall libpetsc.TaoADMMSetRegularizerObjectiveAndGradientRoutine(arg1::Tao, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function TaoADMMSetMisfitHessianRoutine(arg1, arg2, arg3, arg4, arg5)
    @ccall libpetsc.TaoADMMSetMisfitHessianRoutine(arg1::Tao, arg2::Mat, arg3::Mat, arg4::Ptr{Cvoid}, arg5::Ptr{Cvoid})::PetscErrorCode
end

function TaoADMMSetMisfitObjectiveAndGradientRoutine(arg1, arg2, arg3)
    @ccall libpetsc.TaoADMMSetMisfitObjectiveAndGradientRoutine(arg1::Tao, arg2::Ptr{Cvoid}, arg3::Ptr{Cvoid})::PetscErrorCode
end

function TaoADMMSetMisfitHessianChangeStatus(arg1, arg2)
    @ccall libpetsc.TaoADMMSetMisfitHessianChangeStatus(arg1::Tao, arg2::PetscBool)::PetscErrorCode
end

function TaoADMMSetRegHessianChangeStatus(arg1, arg2)
    @ccall libpetsc.TaoADMMSetRegHessianChangeStatus(arg1::Tao, arg2::PetscBool)::PetscErrorCode
end

function TaoADMMSetMinimumSpectralPenalty(arg1, arg2)
    @ccall libpetsc.TaoADMMSetMinimumSpectralPenalty(arg1::Tao, arg2::PetscReal)::PetscErrorCode
end

function TaoADMMSetRegularizerType(arg1, arg2)
    @ccall libpetsc.TaoADMMSetRegularizerType(arg1::Tao, arg2::TaoADMMRegularizerType)::PetscErrorCode
end

function TaoADMMGetRegularizerType(arg1, arg2)
    @ccall libpetsc.TaoADMMGetRegularizerType(arg1::Tao, arg2::Ptr{TaoADMMRegularizerType})::PetscErrorCode
end

function TaoADMMSetUpdateType(arg1, arg2)
    @ccall libpetsc.TaoADMMSetUpdateType(arg1::Tao, arg2::TaoADMMUpdateType)::PetscErrorCode
end

function TaoADMMGetUpdateType(arg1, arg2)
    @ccall libpetsc.TaoADMMGetUpdateType(arg1::Tao, arg2::Ptr{TaoADMMUpdateType})::PetscErrorCode
end

function TaoALMMGetType(arg1, arg2)
    @ccall libpetsc.TaoALMMGetType(arg1::Tao, arg2::Ptr{TaoALMMType})::PetscErrorCode
end

function TaoALMMSetType(arg1, arg2)
    @ccall libpetsc.TaoALMMSetType(arg1::Tao, arg2::TaoALMMType)::PetscErrorCode
end

function TaoALMMGetSubsolver(arg1, arg2)
    @ccall libpetsc.TaoALMMGetSubsolver(arg1::Tao, arg2::Ptr{Tao})::PetscErrorCode
end

function TaoALMMSetSubsolver(arg1, arg2)
    @ccall libpetsc.TaoALMMSetSubsolver(arg1::Tao, arg2::Tao)::PetscErrorCode
end

function TaoALMMGetMultipliers(arg1, arg2)
    @ccall libpetsc.TaoALMMGetMultipliers(arg1::Tao, arg2::Ptr{Vec})::PetscErrorCode
end

function TaoALMMSetMultipliers(arg1, arg2)
    @ccall libpetsc.TaoALMMSetMultipliers(arg1::Tao, arg2::Vec)::PetscErrorCode
end

function TaoALMMGetPrimalIS(arg1, arg2, arg3)
    @ccall libpetsc.TaoALMMGetPrimalIS(arg1::Tao, arg2::Ptr{IS}, arg3::Ptr{IS})::PetscErrorCode
end

function TaoALMMGetDualIS(arg1, arg2, arg3)
    @ccall libpetsc.TaoALMMGetDualIS(arg1::Tao, arg2::Ptr{IS}, arg3::Ptr{IS})::PetscErrorCode
end

const PETSC_ARCH = ""

const PETSC_BLASLAPACK_UNDERSCORE = 1

const PETSC_CLANGUAGE_C = 1

# Skipping MacroDefinition: PETSC_CXX_INLINE inline

const PETSC_CXX_RESTRICT = __restrict

# Skipping MacroDefinition: PETSC_C_INLINE inline

const PETSC_C_RESTRICT = __restrict

const PETSC_DIR = "/workspace/destdir"

const PETSC_DIR_SEPARATOR = Cchar('/')

const PETSC_DO_NOT_SWAP_CHILD_FOR_DEBUGGER = 1

const PETSC_FORTRAN_CHARLEN_T = size_t

# Skipping MacroDefinition: PETSC_FORTRAN_TYPE_INITIALIZE = - 2

# Skipping MacroDefinition: PETSC_FUNCTION_NAME_C __func__

# Skipping MacroDefinition: PETSC_FUNCTION_NAME_CXX __func__

const PETSC_HAVE_ACCESS = 1

const PETSC_HAVE_ATOLL = 1

const PETSC_HAVE_ATTRIBUTEALIGNED = 1

const PETSC_HAVE_BUILTIN_EXPECT = 1

const PETSC_HAVE_BZERO = 1

const PETSC_HAVE_C99_COMPLEX = 1

const PETSC_HAVE_CLOCK = 1

const PETSC_HAVE_CLOSURE = 1

const PETSC_HAVE_CXX = 1

const PETSC_HAVE_CXX_COMPLEX = 1

const PETSC_HAVE_CXX_DIALECT_CXX11 = 1

const PETSC_HAVE_CXX_DIALECT_CXX14 = 1

const PETSC_HAVE_DLADDR = 1

const PETSC_HAVE_DLCLOSE = 1

const PETSC_HAVE_DLERROR = 1

const PETSC_HAVE_DLOPEN = 1

const PETSC_HAVE_DLSYM = 1

const PETSC_HAVE_DOUBLE_ALIGN_MALLOC = 1

const PETSC_HAVE_DRAND48 = 1

const PETSC_HAVE_DYNAMIC_LIBRARIES = 1

const PETSC_HAVE_ERF = 1

const PETSC_HAVE_FORK = 1

const PETSC_HAVE_FORTRAN = 1

const PETSC_HAVE_FORTRAN_FLUSH = 1

const PETSC_HAVE_FORTRAN_GET_COMMAND_ARGUMENT = 1

const PETSC_HAVE_FORTRAN_TYPE_STAR = 1

const PETSC_HAVE_FORTRAN_UNDERSCORE = 1

const PETSC_HAVE_GETCWD = 1

const PETSC_HAVE_GETDOMAINNAME = 1

const PETSC_HAVE_GETHOSTBYNAME = 1

const PETSC_HAVE_GETHOSTNAME = 1

const PETSC_HAVE_GETPAGESIZE = 1

const PETSC_HAVE_GETRUSAGE = 1

const PETSC_HAVE_GETWD = 1

const PETSC_HAVE_ISINF = 1

const PETSC_HAVE_ISNAN = 1

const PETSC_HAVE_ISNORMAL = 1

const PETSC_HAVE_LGAMMA = 1

const PETSC_HAVE_LOG2 = 1

const PETSC_HAVE_LSEEK = 1

const PETSC_HAVE_MEMMOVE = 1

const PETSC_HAVE_MMAP = 1

const PETSC_HAVE_MPICH_NUMVERSION = 30402300

const PETSC_HAVE_MPIIO = 1

const PETSC_HAVE_MPI_COMBINER_CONTIGUOUS = 1

const PETSC_HAVE_MPI_COMBINER_DUP = 1

const PETSC_HAVE_MPI_COMBINER_NAMED = 1

const PETSC_HAVE_MPI_EXSCAN = 1

const PETSC_HAVE_MPI_F90MODULE = 1

const PETSC_HAVE_MPI_F90MODULE_VISIBILITY = 1

const PETSC_HAVE_MPI_FINALIZED = 1

const PETSC_HAVE_MPI_GET_ACCUMULATE = 1

const PETSC_HAVE_MPI_GET_LIBRARY_VERSION = 1

const PETSC_HAVE_MPI_IALLREDUCE = 1

const PETSC_HAVE_MPI_IBARRIER = 1

const PETSC_HAVE_MPI_INIT_THREAD = 1

const PETSC_HAVE_MPI_INT64_T = 1

const PETSC_HAVE_MPI_IN_PLACE = 1

const PETSC_HAVE_MPI_LONG_DOUBLE = 1

const PETSC_HAVE_MPI_NONBLOCKING_COLLECTIVES = 1

const PETSC_HAVE_MPI_ONE_SIDED = 1

const PETSC_HAVE_MPI_PROCESS_SHARED_MEMORY = 1

const PETSC_HAVE_MPI_REDUCE_LOCAL = 1

const PETSC_HAVE_MPI_REDUCE_SCATTER = 1

const PETSC_HAVE_MPI_REDUCE_SCATTER_BLOCK = 1

const PETSC_HAVE_MPI_RGET = 1

const PETSC_HAVE_MPI_TYPE_DUP = 1

const PETSC_HAVE_MPI_TYPE_GET_ENVELOPE = 1

const PETSC_HAVE_MPI_WIN_CREATE = 1

const PETSC_HAVE_NANOSLEEP = 1

const PETSC_HAVE_PACKAGES = ":blaslapack:mathlib:mpi:pthread:regex:"

const PETSC_HAVE_POPEN = 1

const PETSC_HAVE_PTHREAD = 1

const PETSC_HAVE_RAND = 1

const PETSC_HAVE_READLINK = 1

const PETSC_HAVE_REALPATH = 1

const PETSC_HAVE_REGEX = 1

const PETSC_HAVE_RTLD_GLOBAL = 1

const PETSC_HAVE_RTLD_LAZY = 1

const PETSC_HAVE_RTLD_LOCAL = 1

const PETSC_HAVE_RTLD_NOW = 1

const PETSC_HAVE_SLEEP = 1

const PETSC_HAVE_SNPRINTF = 1

const PETSC_HAVE_SOCKET = 1

const PETSC_HAVE_SO_REUSEADDR = 1

const PETSC_HAVE_STRCASECMP = 1

const PETSC_HAVE_STRUCT_SIGACTION = 1

const PETSC_HAVE_TGAMMA = 1

const PETSC_HAVE_TIME = 1

const PETSC_HAVE_UNAME = 1

const PETSC_HAVE_USLEEP = 1

const PETSC_HAVE_VA_COPY = 1

const PETSC_HAVE_VSNPRINTF = 1

const PETSC_IS_COLORING_MAX = USHRT_MAX

const PETSC_IS_COLORING_VALUE_TYPE = Cshort

const PETSC_IS_COLORING_VALUE_TYPE_F = integer2

const PETSC_LEVEL1_DCACHE_LINESIZE = 32

const PETSC_LIB_DIR = "/workspace/destdir/lib"

const PETSC_MAX_PATH_LEN = 1024

const PETSC_MEMALIGN = 16

const PETSC_MPICC_SHOW = "Unavailable"

const PETSC_MPIU_IS_COLORING_VALUE_TYPE = MPI_UNSIGNED_SHORT

const PETSC_PREFETCH_HINT_NTA = _MM_HINT_NTA

const PETSC_PREFETCH_HINT_T0 = _MM_HINT_T0

const PETSC_PREFETCH_HINT_T1 = _MM_HINT_T1

const PETSC_PREFETCH_HINT_T2 = _MM_HINT_T2

const PETSC_PYTHON_EXE = "/usr/bin/python3"

const PETSC_REPLACE_DIR_SEPARATOR = Cchar('\\')

const PETSC_RTLD_DEFAULT = 1

const PETSC_SIZEOF_ENUM = 4

const PETSC_SIZEOF_INT = 4

const PETSC_SIZEOF_LONG = 8

const PETSC_SIZEOF_LONG_LONG = 8

const PETSC_SIZEOF_SHORT = 2

const PETSC_SIZEOF_SIZE_T = 8

const PETSC_SIZEOF_VOID_P = 8

const PETSC_SLSUFFIX = "dylib"

const PETSC_UINTPTR_T = uintptr_t

const PETSC_UNUSED = __attribute(unused)

const PETSC_USE_64BIT_INDICES = 1

const PETSC_USE_AVX512_KERNELS = 1

const PETSC_USE_BACKWARD_LOOP = 1

const PETSC_USE_COMPLEX = 1

const PETSC_USE_CTABLE = 1

const PETSC_USE_DEBUGGER = "gdb"

const PETSC_USE_INFO = 1

const PETSC_USE_LOG = 1

const PETSC_USE_MALLOC_COALESCED = 1

const PETSC_USE_REAL_SINGLE = 1

const PETSC_USE_SHARED_LIBRARIES = 1

const PETSC_USE_SINGLE_LIBRARY = 1

const PETSC_USE_SOCKET_VIEWER = 1

const PETSC_USE_VISIBILITY_C = 1

const PETSC_USE_VISIBILITY_CXX = 1

const PETSC_USING_64BIT_PTR = 1

const PETSC_USING_F2003 = 1

const PETSC_USING_F90FREEFORM = 1

const PETSC__BSD_SOURCE = 1

const PETSC__DEFAULT_SOURCE = 1

const PETSC__GNU_SOURCE = 1

const PETSC_HAVE_COMPLEX = 1

const PETSC_REAL = PETSC_FLOAT

const PETSC_SCALAR = PETSC_COMPLEX

const PETSC_FORTRANADDR = PETSC_LONG

const PETSC_BINARY_INT_SIZE = 32 ÷ 8

const PETSC_BINARY_FLOAT_SIZE = 32 ÷ 8

const PETSC_BINARY_CHAR_SIZE = 8 ÷ 8

const PETSC_BINARY_SHORT_SIZE = 16 ÷ 8

const PETSC_BINARY_DOUBLE_SIZE = 64 ÷ 8

# Skipping MacroDefinition: PETSC_BINARY_SCALAR_SIZE sizeof ( PetscScalar )

const PETSC_FUNCTION_NAME = PETSC_FUNCTION_NAME_C

const PETSC_RESTRICT = PETSC_C_RESTRICT

const PETSC_INLINE = PETSC_C_INLINE

# Skipping MacroDefinition: PETSC_STATIC_INLINE static PETSC_INLINE

# Skipping MacroDefinition: PETSC_DLLEXPORT __attribute__ ( ( visibility ( "default" ) ) )

# Skipping MacroDefinition: PETSC_DLLIMPORT __attribute__ ( ( visibility ( "default" ) ) )

# Skipping MacroDefinition: PETSC_VISIBILITY_INTERNAL __attribute__ ( ( visibility ( "hidden" ) ) )

const PETSC_VISIBILITY_PUBLIC = PETSC_DLLIMPORT

# Skipping MacroDefinition: PETSC_EXTERN extern PETSC_VISIBILITY_PUBLIC

# Skipping MacroDefinition: PETSC_INTERN extern PETSC_VISIBILITY_INTERNAL

const PETSC_VERSION_RELEASE = 1

const PETSC_VERSION_MAJOR = 3

const PETSC_VERSION_MINOR = 15

const PETSC_VERSION_SUBMINOR = 2

const PETSC_VERSION_PATCH = 0

const PETSC_RELEASE_DATE = "Mar 30, 2021"

const PETSC_VERSION_DATE = "Jul 10, 2021"

const PETSC_VERSION_GIT = "v3.15.2"

const PETSC_VERSION_DATE_GIT = "2021-07-10 11:22:33 -0500"

const PETSC_VERSION_ = PETSC_VERSION_EQ

const PETSC_AUTHOR_INFO = "       The PETSc Team\n    petsc-maint@mcs.anl.gov\n https://www.mcs.anl.gov/petsc/\n"

const MPICH_SKIP_MPICXX = 1

const OMPI_SKIP_MPICXX = 1

const PetscDefined_arg_1 = $(Expr(:incomplete, "incomplete: premature end of input"))

const PetscDefined_arg_ = $(Expr(:incomplete, "incomplete: premature end of input"))

const MPIU_INT64 = MPI_INT64_T

const PetscInt64_FMT = PRId64

const MPIU_INT = MPIU_INT64

const PetscInt_FMT = PetscInt64_FMT

const MPIU_REAL = MPI_FLOAT

const MPIU_C_COMPLEX = (MPI_C_COMPLEX(PETSC_DEPRECATED_MACRO))("GCC warning \"MPIU_C_COMPLEX macro is deprecated use MPI_C_COMPLEX (since version 3.15)\"")

const MPIU_C_DOUBLE_COMPLEX = (MPI_C_DOUBLE_COMPLEX(PETSC_DEPRECATED_MACRO))("GCC warning \"MPIU_C_DOUBLE_COMPLEX macro is deprecated use MPI_C_DOUBLE_COMPLEX (since version 3.15)\"")

const MPIU_COMPLEX = MPI_C_COMPLEX

const MPIU_SCALAR = MPIU_COMPLEX

const PETSC_PI = PetscRealConstant(3.141592653589793)

const PETSC_PHI = PetscRealConstant(1.618033988749895)

const PETSC_SQRT2 = PetscRealConstant(1.4142135623730951)

const PETSC_MAX_INT = Clong(9223372036854775807)

const PETSC_MIN_INT = -PETSC_MAX_INT - 1

const PETSC_MAX_UINT16 = 65535

const PETSC_MAX_REAL = Float32(3.4028234663852886e38)

const PETSC_MIN_REAL = -PETSC_MAX_REAL

const PETSC_MACHINE_EPSILON = Float32(1.1920929e-7)

const PETSC_SQRT_MACHINE_EPSILON = Float32(0.000345266983)

const PETSC_SMALL = Float32(1.0e-5)

const PETSC_INFINITY = PETSC_MAX_REAL ÷ 4

const PETSC_NINFINITY = -PETSC_INFINITY

const MPIU_MATSCALAR = MPIU_SCALAR

const PETSC_IGNORE = NULL

const PETSC_NULL = NULL

const PETSC_DECIDE = -1

const PETSC_DETERMINE = PETSC_DECIDE

const PETSC_DEFAULT = -2

const PETSC_COMM_SELF = MPI_COMM_SELF

const MPIU_PETSCLOGDOUBLE = MPI_DOUBLE

const MPIU_2PETSCLOGDOUBLE = MPI_2DOUBLE_PRECISION

const MPIU_SUM = MPI_SUM

const MPIU_MAX = MPI_MAX

const MPIU_MIN = MPI_MIN

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

const PETSC_ERR_MPI_LIB_INCOMP = 87

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

const PETSC_ERR_OPT_OVERWRITE = 93

const PETSC_ERR_WRONG_MPI_SIZE = 94

const PETSC_ERR_USER_INPUT = 95

const PETSC_ERR_GPU_RESOURCE = 96

const PETSC_ERR_GPU = 97

const PETSC_ERR_MPI = 98

const PETSC_ERR_MAX_VALUE = 99

# Skipping MacroDefinition: CHKMEMQ do { PetscErrorCode _7_ierr = PetscMallocValidate ( __LINE__ , PETSC_FUNCTION_NAME , __FILE__ ) ; CHKERRQ ( _7_ierr ) ; } while ( 0 )

const CHKMEMA = PetscMallocValidate(__LINE__, PETSC_FUNCTION_NAME, __FILE__)

const PETSCSTACKSIZE = 64

# Skipping MacroDefinition: PetscStackPopNoCheck do { } while ( 0 )

const PetscStackPop = CHKMEMQ

const PETSC_SMALLEST_CLASSID = 1211211

const PETSC_MAX_OPTION_NAME = 512

const PETSC_EVENT = 1311311

const PETSC_FLOPS_PER_OP = 4.0

const PETSC_MPI_INT_MAX = 2147483647

const PETSC_MPI_INT_MIN = -2147483647

const PETSC_BLAS_INT_MAX = 2147483647

const PETSC_BLAS_INT_MIN = -2147483647

const PETSC_BITS_PER_BYTE = CHAR_BIT

const PETSCRAND = "rand"

const PETSCRAND48 = "rand48"

const PETSCSPRNG = "sprng"

const PETSCRANDER48 = "rander48"

const PETSCRANDOM123 = "random123"

const PETSCCURAND = "curand"

const PETSC_BAG_FILE_CLASSID = 1211219

const PETSC_DRAW_X = "x"

const PETSC_DRAW_NULL = "null"

const PETSC_DRAW_WIN32 = "win32"

const PETSC_DRAW_TIKZ = "tikz"

const PETSC_DRAW_IMAGE = "image"

const PETSCVIEWERSOCKET = "socket"

const PETSCVIEWERASCII = "ascii"

const PETSCVIEWERBINARY = "binary"

const PETSCVIEWERSTRING = "string"

const PETSCVIEWERDRAW = "draw"

const PETSCVIEWERVU = "vu"

const PETSCVIEWERMATHEMATICA = "mathematica"

const PETSCVIEWERHDF5 = "hdf5"

const PETSCVIEWERVTK = "vtk"

const PETSCVIEWERMATLAB = "matlab"

const PETSCVIEWERSAWS = "saws"

const PETSCVIEWERGLVIS = "glvis"

const PETSCVIEWERADIOS = "adios"

const PETSCVIEWERADIOS2 = "adios2"

const PETSCVIEWEREXODUSII = "exodusii"

const PETSC_VIEWER_ASCII_VTK_ATTR = (PETSC_VIEWER_ASCII_VTK(PETSC_DEPRECATED_ENUM))("Legacy VTK deprecated; use PetscViewerVTKOpen() with XML (.vtr .vts .vtu) format (since 3.14)")

const PETSC_VIEWER_ASCII_VTK_CELL_ATTR = (PETSC_VIEWER_ASCII_VTK_CELL(PETSC_DEPRECATED_ENUM))("Legacy VTK deprecated; use PetscViewerVTKOpen() with XML (.vtr .vts .vtu) format (since 3.14)")

const PETSC_VIEWER_ASCII_VTK_COORDS_ATTR = (PETSC_VIEWER_ASCII_VTK_COORDS(PETSC_DEPRECATED_ENUM))("Legacy VTK deprecated; use PetscViewerVTKOpen() with XML (.vtr .vts .vtu) format (since 3.14)")

const PETSC_VIEWER_STDERR_SELF = PETSC_VIEWER_STDERR_(PETSC_COMM_SELF)

const PETSC_VIEWER_STDERR_WORLD = PETSC_VIEWER_STDERR_(PETSC_COMM_WORLD)

const PETSC_VIEWER_STDOUT_WORLD = PETSC_VIEWER_STDOUT_(PETSC_COMM_WORLD)

const PETSC_VIEWER_STDOUT_SELF = PETSC_VIEWER_STDOUT_(PETSC_COMM_SELF)

const PETSC_VIEWER_DRAW_WORLD = PETSC_VIEWER_DRAW_(PETSC_COMM_WORLD)

const PETSC_VIEWER_DRAW_SELF = PETSC_VIEWER_DRAW_(PETSC_COMM_SELF)

const PETSC_VIEWER_SOCKET_WORLD = PETSC_VIEWER_SOCKET_(PETSC_COMM_WORLD)

const PETSC_VIEWER_SOCKET_SELF = PETSC_VIEWER_SOCKET_(PETSC_COMM_SELF)

const PETSC_VIEWER_BINARY_WORLD = PETSC_VIEWER_BINARY_(PETSC_COMM_WORLD)

const PETSC_VIEWER_BINARY_SELF = PETSC_VIEWER_BINARY_(PETSC_COMM_SELF)

const PETSC_VIEWER_MATLAB_WORLD = PETSC_VIEWER_MATLAB_(PETSC_COMM_WORLD)

const PETSC_VIEWER_MATLAB_SELF = PETSC_VIEWER_MATLAB_(PETSC_COMM_SELF)

const PETSC_VIEWER_MATHEMATICA_WORLD = (PetscViewerInitializeMathematicaWorld_Private(), PETSC_VIEWER_MATHEMATICA_WORLD_PRIVATE)

const PETSC_MATLAB_ENGINE_WORLD = PETSC_MATLAB_ENGINE_(PETSC_COMM_WORLD)

const PETSC_MATLAB_ENGINE_SELF = PETSC_MATLAB_ENGINE_(PETSC_COMM_SELF)

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

const PETSC_DRAW_MAXCOLOR = 256

const PETSC_DRAW_FULL_SIZE = -3

const PETSC_DRAW_HALF_SIZE = -4

const PETSC_DRAW_THIRD_SIZE = -5

const PETSC_DRAW_QUARTER_SIZE = -6

const IS_FILE_CLASSID = 1211218

const ISGENERAL = "general"

const ISSTRIDE = "stride"

const ISBLOCK = "block"

const ISLOCALTOGLOBALMAPPINGBASIC = "basic"

const ISLOCALTOGLOBALMAPPINGHASH = "hash"

const IS_COLORING_MAX = PETSC_IS_COLORING_MAX

const MPIU_COLORING_VALUE = PETSC_MPIU_IS_COLORING_VALUE_TYPE

const VECSEQ = "seq"

const VECMPI = "mpi"

const VECSTANDARD = "standard"

const VECSHARED = "shared"

const VECSEQVIENNACL = "seqviennacl"

const VECMPIVIENNACL = "mpiviennacl"

const VECVIENNACL = "viennacl"

const VECSEQCUDA = "seqcuda"

const VECMPICUDA = "mpicuda"

const VECCUDA = "cuda"

const VECSEQHIP = "seqhip"

const VECMPIHIP = "mpihip"

const VECHIP = "hip"

const VECNEST = "nest"

const VECSEQKOKKOS = "seqkokkos"

const VECMPIKOKKOS = "mpikokkos"

const VECKOKKOS = "kokkos"

const REAL_FILE_CLASSID = 1211213

const VEC_FILE_CLASSID = 1211214

const NORM_MAX = NORM_INFINITY

const VECTAGGERABSOLUTE = "absolute"

const VECTAGGERRELATIVE = "relative"

const VECTAGGERCDF = "cdf"

const VECTAGGEROR = "or"

const VECTAGGERAND = "and"

const PETSCSFBASIC = "basic"

const PETSCSFNEIGHBOR = "neighbor"

const PETSCSFALLGATHERV = "allgatherv"

const PETSCSFALLGATHER = "allgather"

const PETSCSFGATHERV = "gatherv"

const PETSCSFGATHER = "gather"

const PETSCSFALLTOALL = "alltoall"

const PETSCSFWINDOW = "window"

const MPIU_REPLACE = (MPI_REPLACE(PETSC_DEPRECATED_MACRO))("GCC warning \"MPIU_REPLACE macro is deprecated use MPI_REPLACE (since version 3.15)\"")

const MATSAME = "same"

const MATMAIJ = "maij"

const MATSEQMAIJ = "seqmaij"

const MATMPIMAIJ = "mpimaij"

const MATKAIJ = "kaij"

const MATSEQKAIJ = "seqkaij"

const MATMPIKAIJ = "mpikaij"

const MATIS = "is"

const MATAIJ = "aij"

const MATSEQAIJ = "seqaij"

const MATMPIAIJ = "mpiaij"

const MATAIJCRL = "aijcrl"

const MATSEQAIJCRL = "seqaijcrl"

const MATMPIAIJCRL = "mpiaijcrl"

const MATAIJCUSPARSE = "aijcusparse"

const MATSEQAIJCUSPARSE = "seqaijcusparse"

const MATMPIAIJCUSPARSE = "mpiaijcusparse"

const MATAIJKOKKOS = "aijkokkos"

const MATSEQAIJKOKKOS = "seqaijkokkos"

const MATMPIAIJKOKKOS = "mpiaijkokkos"

const MATAIJVIENNACL = "aijviennacl"

const MATSEQAIJVIENNACL = "seqaijviennacl"

const MATMPIAIJVIENNACL = "mpiaijviennacl"

const MATAIJPERM = "aijperm"

const MATSEQAIJPERM = "seqaijperm"

const MATMPIAIJPERM = "mpiaijperm"

const MATAIJSELL = "aijsell"

const MATSEQAIJSELL = "seqaijsell"

const MATMPIAIJSELL = "mpiaijsell"

const MATAIJMKL = "aijmkl"

const MATSEQAIJMKL = "seqaijmkl"

const MATMPIAIJMKL = "mpiaijmkl"

const MATBAIJMKL = "baijmkl"

const MATSEQBAIJMKL = "seqbaijmkl"

const MATMPIBAIJMKL = "mpibaijmkl"

const MATSHELL = "shell"

const MATDENSE = "dense"

const MATDENSECUDA = "densecuda"

const MATSEQDENSE = "seqdense"

const MATSEQDENSECUDA = "seqdensecuda"

const MATMPIDENSE = "mpidense"

const MATMPIDENSECUDA = "mpidensecuda"

const MATELEMENTAL = "elemental"

const MATSCALAPACK = "scalapack"

const MATBAIJ = "baij"

const MATSEQBAIJ = "seqbaij"

const MATMPIBAIJ = "mpibaij"

const MATMPIADJ = "mpiadj"

const MATSBAIJ = "sbaij"

const MATSEQSBAIJ = "seqsbaij"

const MATMPISBAIJ = "mpisbaij"

const MATMFFD = "mffd"

const MATNORMAL = "normal"

const MATNORMALHERMITIAN = "normalh"

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

const MATHYPRE = "hypre"

const MATHYPRESTRUCT = "hyprestruct"

const MATHYPRESSTRUCT = "hypresstruct"

const MATSUBMATRIX = "submatrix"

const MATLOCALREF = "localref"

const MATNEST = "nest"

const MATPREALLOCATOR = "preallocator"

const MATSELL = "sell"

const MATSEQSELL = "seqsell"

const MATMPISELL = "mpisell"

const MATDUMMY = "dummy"

const MATLMVM = "lmvm"

const MATLMVMDFP = "lmvmdfp"

const MATLMVMBFGS = "lmvmbfgs"

const MATLMVMSR1 = "lmvmsr1"

const MATLMVMBROYDEN = "lmvmbroyden"

const MATLMVMBADBROYDEN = "lmvmbadbroyden"

const MATLMVMSYMBROYDEN = "lmvmsymbroyden"

const MATLMVMSYMBADBROYDEN = "lmvmsymbadbroyden"

const MATLMVMDIAGBROYDEN = "lmvmdiagbroyden"

const MATCONSTANTDIAGONAL = "constantdiagonal"

const MATHARA = "hara"

const MATSOLVERSUPERLU = "superlu"

const MATSOLVERSUPERLU_DIST = "superlu_dist"

const MATSOLVERSTRUMPACK = "strumpack"

const MATSOLVERUMFPACK = "umfpack"

const MATSOLVERCHOLMOD = "cholmod"

const MATSOLVERKLU = "klu"

const MATSOLVERSPARSEELEMENTAL = "sparseelemental"

const MATSOLVERELEMENTAL = "elemental"

const MATSOLVERSCALAPACK = "scalapack"

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

const MATSOLVERCUSPARSEBAND = "cusparseband"

const MATSOLVERCUDA = "cuda"

const MATSOLVERKOKKOS = "kokkos"

const MATSOLVERKOKKOSDEVICE = "kokkosdevice"

const MATPRODUCTALGORITHM_DEFAULT = "default"

const MAT_FILE_CLASSID = 1211216

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

const MATORDERINGNATURAL_OR_ND = "natural_or_nd"

const MATORDERINGEXTERNAL = "external"

const MATCOLORINGJP = "jp"

const MATCOLORINGPOWER = "power"

const MATCOLORINGNATURAL = "natural"

const MATCOLORINGSL = "sl"

const MATCOLORINGLF = "lf"

const MATCOLORINGID = "id"

const MATCOLORINGGREEDY = "greedy"

const MATPARTITIONINGCURRENT = "current"

const MATPARTITIONINGAVERAGE = "average"

const MATPARTITIONINGSQUARE = "square"

const MATPARTITIONINGPARMETIS = "parmetis"

const MATPARTITIONINGCHACO = "chaco"

const MATPARTITIONINGPARTY = "party"

const MATPARTITIONINGPTSCOTCH = "ptscotch"

const MATPARTITIONINGHIERARCH = "hierarch"

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

const MATRIX_BINARY_FORMAT_DENSE = -1

const MATMFFD_DS = "ds"

const MATMFFD_WP = "wp"

const PETSCSECTIONSYMLABEL = "label"

const DMLOCATEPOINT_POINT_NOT_FOUND = -367

const DMDA = "da"

const DMCOMPOSITE = "composite"

const DMSLICED = "sliced"

const DMSHELL = "shell"

const DMPLEX = "plex"

const DMREDUNDANT = "redundant"

const DMPATCH = "patch"

const DMMOAB = "moab"

const DMNETWORK = "network"

const DMFOREST = "forest"

const DMP4EST = "p4est"

const DMP8EST = "p8est"

const DMSWARM = "swarm"

const DMPRODUCT = "product"

const DMSTAG = "stag"

const DM_FILE_CLASSID = 1211221

const PFCONSTANT = "constant"

const PFMAT = "mat"

const PFSTRING = "string"

const PFQUICK = "quick"

const PFIDENTITY = "identity"

const PFMATLAB = "matlab"

const AOBASIC = "basic"

const AOADVANCED = "advanced"

const AOMAPPING = "mapping"

const AOMEMORYSCALABLE = "memoryscalable"

const PETSC_FACTORIAL_MAX = 20

const PETSC_BINOMIAL_MAX = 61

const PETSCSPACEPOLYNOMIAL = "poly"

const PETSCSPACETENSOR = "tensor"

const PETSCSPACESUM = "sum"

const PETSCSPACEPOINT = "point"

const PETSCSPACESUBSPACE = "subspace"

const PETSCDUALSPACELAGRANGE = "lagrange"

const PETSCDUALSPACESIMPLE = "simple"

const PETSCDUALSPACEREFINED = "refined"

const PETSCDUALSPACEBDM = "bdm"

const PETSCFEBASIC = "basic"

const PETSCFEOPENCL = "opencl"

const PETSCFECOMPOSITE = "composite"

const MATSEQUSFFT = "sequsfft"

const PETSCPARTITIONERPARMETIS = "parmetis"

const PETSCPARTITIONERPTSCOTCH = "ptscotch"

const PETSCPARTITIONERCHACO = "chaco"

const PETSCPARTITIONERSIMPLE = "simple"

const PETSCPARTITIONERSHELL = "shell"

const PETSCPARTITIONERGATHER = "gather"

const PETSCPARTITIONERMATPARTITIONING = "matpartitioning"

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

const DMFIELDDA = "da"

const DMFIELDDS = "ds"

const DMFIELDSHELL = "shell"

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

const PCVPBJACOBI = "vpbjacobi"

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

const PCCHOWILUVIENNACL = "chowiluviennacl"

const PCROWSCALINGVIENNACL = "rowscalingviennacl"

const PCSAVIENNACL = "saviennacl"

const PCBDDC = "bddc"

const PCKACZMARZ = "kaczmarz"

const PCTELESCOPE = "telescope"

const PCPATCH = "patch"

const PCLMVM = "lmvm"

const PCHMG = "hmg"

const PCDEFLATION = "deflation"

const PCHPDDM = "hpddm"

const PCHARA = "hara"

const PC_SIDE_MAX = PC_SYMMETRIC + 1

const PCGAMGAGG = "agg"

const PCGAMGGEO = "geo"

const PCGAMGCLASSICAL = "classical"

const PCGAMGCLASSICALDIRECT = "direct"

const PCGAMGCLASSICALSTANDARD = "standard"

const PC_MG_CASCADE = $(Expr(:toplevel, :PC_MG_KASKADE))

const PC_FILE_CLASSID = 1211222

const KSPRICHARDSON = "richardson"

const KSPCHEBYSHEV = "chebyshev"

const KSPCG = "cg"

const KSPGROPPCG = "groppcg"

const KSPPIPECG = "pipecg"

const KSPPIPECGRR = "pipecgrr"

const KSPPIPELCG = "pipelcg"

const KSPPIPEPRCG = "pipeprcg"

const KSPPIPECG2 = "pipecg2"

const KSPCGNE = "cgne"

const KSPNASH = "nash"

const KSPSTCG = "stcg"

const KSPGLTR = "gltr"

const KSPCGNASH = (PETSC_DEPRECATED_MACRO("GCC warning \"KSPCGNASH macro is deprecated use KSPNASH (since version 3.11)\""))("nash")

const KSPCGSTCG = (PETSC_DEPRECATED_MACRO("GCC warning \"KSPCGSTCG macro is deprecated use KSPSTCG (since version 3.11)\""))("stcg")

const KSPCGGLTR = (PETSC_DEPRECATED_MACRO("GCC warning \"KSPCGGLTR macro is deprecated use KSPSGLTR (since version 3.11)\""))("gltr")

const KSPFCG = "fcg"

const KSPPIPEFCG = "pipefcg"

const KSPGMRES = "gmres"

const KSPPIPEFGMRES = "pipefgmres"

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

const KSPPIPEBCGS = "pipebcgs"

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

const KSPPIPEGCR = "pipegcr"

const KSPTSIRM = "tsirm"

const KSPCGLS = "cgls"

const KSPFETIDP = "fetidp"

const KSPHPDDM = "hpddm"

const KSP_FILE_CLASSID = 1211223

const KSP_NORM_MAX = KSP_NORM_NATURAL + 1

const KSP_DIVERGED_PCSETUP_FAILED_DEPRECATED = (KSP_DIVERGED_PCSETUP_FAILED(PETSC_DEPRECATED_ENUM))("Use KSP_DIVERGED_PC_FAILED (since version 3.11)")

const KSPDefaultConverged = (KSPDefaultConverged, KSPConvergedDefault)

const KSPDefaultConvergedDestroy = (KSPDefaultConvergedDestroy, KSPConvergedDefaultDestroy)

const KSPDefaultConvergedCreate = (KSPDefaultConvergedCreate, KSPConvergedDefaultCreate)

const KSPDefaultConvergedSetUIRNorm = (KSPDefaultConvergedSetUIRNorm, KSPConvergedDefaultSetUIRNorm)

const KSPDefaultConvergedSetUMIRNorm = (KSPDefaultConvergedSetUMIRNorm, KSPConvergedDefaultSetUMIRNorm)

const KSPSkipConverged = (KSPSkipConverged, KSPConvergedSkip)

const KSPGUESSFISCHER = "fischer"

const KSPGUESSPOD = "pod"

const SNESNEWTONLS = "newtonls"

const SNESNEWTONTR = "newtontr"

const SNESPYTHON = "python"

const SNESNRICHARDSON = "nrichardson"

const SNESKSPONLY = "ksponly"

const SNESKSPTRANSPOSEONLY = "ksptransposeonly"

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

const SNESPATCH = "patch"

const SNES_FILE_CLASSID = 1211224

const SNES_CONVERGED_TR_DELTA_DEPRECATED = (SNES_CONVERGED_TR_DELTA(PETSC_DEPRECATED_ENUM))("Use SNES_DIVERGED_TR_DELTA (since version 3.12)")

const SNESSkipConverged = (SNESSkipConverged, SNESConvergedSkip)

const SNESLINESEARCHBT = "bt"

const SNESLINESEARCHNLEQERR = "nleqerr"

const SNESLINESEARCHBASIC = "basic"

const SNESLINESEARCHL2 = "l2"

const SNESLINESEARCHCP = "cp"

const SNESLINESEARCHSHELL = "shell"

const SNESLINESEARCHNCGLINEAR = "ncglinear"

const SNES_LINESEARCH_ORDER_LINEAR = 1

const SNES_LINESEARCH_ORDER_QUADRATIC = 2

const SNES_LINESEARCH_ORDER_CUBIC = 3

const SNESMSM62 = "m62"

const SNESMSEULER = "euler"

const SNESMSJAMESON83 = "jameson83"

const SNESMSVLTP11 = "vltp11"

const SNESMSVLTP21 = "vltp21"

const SNESMSVLTP31 = "vltp31"

const SNESMSVLTP41 = "vltp41"

const SNESMSVLTP51 = "vltp51"

const SNESMSVLTP61 = "vltp61"

const TSEULER = "euler"

const TSBEULER = "beuler"

const TSBASICSYMPLECTIC = "basicsymplectic"

const TSPSEUDO = "pseudo"

const TSCN = "cn"

const TSSUNDIALS = "sundials"

const TSRK = "rk"

const TSPYTHON = "python"

const TSTHETA = "theta"

const TSALPHA = "alpha"

const TSALPHA2 = "alpha2"

const TSGLLE = "glle"

const TSGLEE = "glee"

const TSSSP = "ssp"

const TSARKIMEX = "arkimex"

const TSROSW = "rosw"

const TSEIMEX = "eimex"

const TSMIMEX = "mimex"

const TSBDF = "bdf"

const TSRADAU5 = "radau5"

const TSMPRK = "mprk"

const TSDISCGRAD = "discgrad"

const TSTRAJECTORYBASIC = "basic"

const TSTRAJECTORYSINGLEFILE = "singlefile"

const TSTRAJECTORYMEMORY = "memory"

const TSTRAJECTORYVISUALIZATION = "visualization"

const TS_FILE_CLASSID = 1211225

const TSSSPRKS2 = "rks2"

const TSSSPRKS3 = "rks3"

const TSSSPRK104 = "rk104"

const TSADAPTNONE = "none"

const TSADAPTBASIC = "basic"

const TSADAPTDSP = "dsp"

const TSADAPTCFL = "cfl"

const TSADAPTGLEE = "glee"

const TSADAPTHISTORY = "history"

const TSGLLEADAPT_NONE = "none"

const TSGLLEADAPT_SIZE = "size"

const TSGLLEADAPT_BOTH = "both"

const TSGLLEACCEPT_ALWAYS = "always"

const TSGLLE_IRKS = "irks"

const TSEIMEXType = $(Expr(:incomplete, "incomplete: premature end of input"))

const TSRK1FE = "1fe"

const TSRK2A = "2a"

const TSRK3 = "3"

const TSRK3BS = "3bs"

const TSRK4 = "4"

const TSRK5F = "5f"

const TSRK5DP = "5dp"

const TSRK5BS = "5bs"

const TSRK6VR = "6vr"

const TSRK7VR = "7vr"

const TSRK8VR = "8vr"

const TSMPRK2A22 = "2a22"

const TSMPRK2A23 = "2a23"

const TSMPRK2A32 = "2a32"

const TSMPRK2A33 = "2a33"

const TSMPRKP2 = "p2"

const TSMPRKP3 = "p3"

const TSGLEEi1 = "BE1"

const TSGLEE23 = "23"

const TSGLEE24 = "24"

const TSGLEE25I = "25i"

const TSGLEE35 = "35"

const TSGLEEEXRK2A = "exrk2a"

const TSGLEERK32G1 = "rk32g1"

const TSGLEERK285EX = "rk285ex"

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

const TSBASICSYMPLECTICSIEULER = "1"

const TSBASICSYMPLECTICVELVERLET = "2"

const TSBASICSYMPLECTIC3 = "3"

const TSBASICSYMPLECTIC4 = "4"

const TAOLMVM = "lmvm"

const TAONLS = "nls"

const TAONTR = "ntr"

const TAONTL = "ntl"

const TAOCG = "cg"

const TAOTRON = "tron"

const TAOOWLQN = "owlqn"

const TAOBMRM = "bmrm"

const TAOBLMVM = "blmvm"

const TAOBQNLS = "bqnls"

const TAOBNCG = "bncg"

const TAOBNLS = "bnls"

const TAOBNTR = "bntr"

const TAOBNTL = "bntl"

const TAOBQNKLS = "bqnkls"

const TAOBQNKTR = "bqnktr"

const TAOBQNKTL = "bqnktl"

const TAOBQPIP = "bqpip"

const TAOGPCG = "gpcg"

const TAONM = "nm"

const TAOPOUNDERS = "pounders"

const TAOBRGN = "brgn"

const TAOLCL = "lcl"

const TAOSSILS = "ssils"

const TAOSSFLS = "ssfls"

const TAOASILS = "asils"

const TAOASFLS = "asfls"

const TAOIPM = "ipm"

const TAOPDIPM = "pdipm"

const TAOSHELL = "shell"

const TAOADMM = "admm"

const TAOALMM = "almm"

const TAOLINESEARCHUNIT = "unit"

const TAOLINESEARCHMT = "more-thuente"

const TAOLINESEARCHGPCG = "gpcg"

const TAOLINESEARCHARMIJO = "armijo"

const TAOLINESEARCHOWARMIJO = "owarmijo"

const TAOLINESEARCHIPM = "ipm"

