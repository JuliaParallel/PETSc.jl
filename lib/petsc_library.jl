using MPI
const MPI_Comm = MPI.MPI_Comm
const MPI_Datatype = MPI.MPI_Datatype
const MPI_File = MPI.MPI_File
const MPI_Aint = MPI.MPI_Aint
const MPI_Info = MPI.MPI_Info
const MPI_Win = MPI.MPI_Win
const MPI_Offset = MPI.MPI_Offset
const MPI_Op = MPI.MPI_Op
const MPI_UNSIGNED_SHORT = MPI.MPI_UNSIGNED_SHORT
const MPI_INT64_T = MPI.MPI_INT64_T
const MPI_FLOAT = MPI.MPI_FLOAT
const MPI_COMM_SELF = MPI.MPI_COMM_SELF
const MPI_DOUBLE = MPI.MPI_DOUBLE
const MPI_SUM = MPI.MPI_SUM
const MPI_MAX = MPI.MPI_MAX
const MPI_MIN = MPI.MPI_MIN
const MPI_REPLACE = MPI.MPI_REPLACE

# We know these will be Cvoid, so just set them to be that
const PetscOptions = Ptr{Cvoid}
const PetscViewer = Ptr{Cvoid}
const PetscObject = Ptr{Cvoid}
const Vec = Ptr{Cvoid}
const Mat = Ptr{Cvoid}
const KSP = Ptr{Cvoid}
const VecType = Cstring
const KSPType = Cstring

const __darwin_off_t = Int64

mutable struct ADIOI_FileD end

const off_t = __darwin_off_t

const PetscClassId = Cint

const PetscMPIInt = Cint

@enum PetscEnum::UInt32 begin
    ENUM_DUMMY = 0
end

const PetscShort = Cshort

const PetscChar = Cchar

const PetscInt64 = Int64

const PetscBLASInt = Cint

@enum PetscCopyMode::UInt32 begin
    PETSC_COPY_VALUES = 0
    PETSC_OWN_POINTER = 1
    PETSC_USE_POINTER = 2
end

const PetscLogDouble = Cdouble

mutable struct _p_PetscToken end

const PetscToken = Ptr{_p_PetscToken}

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

@for_petsc function PetscSignReal(::$UnionPetscLib, a)
    ccall((:PetscSignReal, $petsc_library), $PetscReal, ($PetscReal,), a)
end

@for_petsc function PetscCMPLX(::$UnionPetscLib, x, y)
    ccall(
        (:PetscCMPLX, $petsc_library),
        $PetscComplex,
        ($PetscReal, $PetscReal),
        x,
        y,
    )
end

@enum PetscScalarPrecision::UInt32 begin
    PETSC_SCALAR_DOUBLE = 0
    PETSC_SCALAR_SINGLE = 1
    PETSC_SCALAR_LONG_DOUBLE = 2
    PETSC_SCALAR_HALF = 3
end

@for_petsc function PetscIsInfReal(::$UnionPetscLib, arg1)
    ccall((:PetscIsInfReal, $petsc_library), PetscBool, ($PetscReal,), arg1)
end

@for_petsc function PetscIsNanReal(::$UnionPetscLib, arg1)
    ccall((:PetscIsNanReal, $petsc_library), PetscBool, ($PetscReal,), arg1)
end

@for_petsc function PetscIsNormalReal(::$UnionPetscLib, arg1)
    ccall((:PetscIsNormalReal, $petsc_library), PetscBool, ($PetscReal,), arg1)
end

@for_petsc function PetscIsInfOrNanReal(::$UnionPetscLib, v)
    ccall((:PetscIsInfOrNanReal, $petsc_library), PetscBool, ($PetscReal,), v)
end

@for_petsc function PetscIsInfScalar(::$UnionPetscLib, v)
    ccall((:PetscIsInfScalar, $petsc_library), PetscBool, ($PetscScalar,), v)
end

@for_petsc function PetscIsNanScalar(::$UnionPetscLib, v)
    ccall((:PetscIsNanScalar, $petsc_library), PetscBool, ($PetscScalar,), v)
end

@for_petsc function PetscIsInfOrNanScalar(::$UnionPetscLib, v)
    ccall(
        (:PetscIsInfOrNanScalar, $petsc_library),
        PetscBool,
        ($PetscScalar,),
        v,
    )
end

@for_petsc function PetscIsNormalScalar(::$UnionPetscLib, v)
    ccall((:PetscIsNormalScalar, $petsc_library), PetscBool, ($PetscScalar,), v)
end

@for_petsc function PetscIsCloseAtTol(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    ccall(
        (:PetscIsCloseAtTol, $petsc_library),
        PetscBool,
        ($PetscReal, $PetscReal, $PetscReal, $PetscReal),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscEqualReal(::$UnionPetscLib, arg1, arg2)
    ccall(
        (:PetscEqualReal, $petsc_library),
        PetscBool,
        ($PetscReal, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function PetscEqualScalar(::$UnionPetscLib, arg1, arg2)
    ccall(
        (:PetscEqualScalar, $petsc_library),
        PetscBool,
        ($PetscScalar, $PetscScalar),
        arg1,
        arg2,
    )
end

@for_petsc function PetscPowInt(::$UnionPetscLib, base, power)
    ccall(
        (:PetscPowInt, $petsc_library),
        $PetscInt,
        ($PetscInt, $PetscInt),
        base,
        power,
    )
end

@for_petsc function PetscPowInt64(::$UnionPetscLib, base, power)
    ccall(
        (:PetscPowInt64, $petsc_library),
        PetscInt64,
        ($PetscInt, $PetscInt),
        base,
        power,
    )
end

@for_petsc function PetscPowRealInt(::$UnionPetscLib, base, power)
    ccall(
        (:PetscPowRealInt, $petsc_library),
        $PetscReal,
        ($PetscReal, $PetscInt),
        base,
        power,
    )
end

@for_petsc function PetscPowScalarInt(::$UnionPetscLib, base, power)
    ccall(
        (:PetscPowScalarInt, $petsc_library),
        $PetscScalar,
        ($PetscScalar, $PetscInt),
        base,
        power,
    )
end

@for_petsc function PetscPowScalarReal(::$UnionPetscLib, base, power)
    ccall(
        (:PetscPowScalarReal, $petsc_library),
        $PetscScalar,
        ($PetscScalar, $PetscReal),
        base,
        power,
    )
end

@for_petsc function PetscLinearRegression(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscLinearRegression, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscSetHelpVersionFunctions(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSetHelpVersionFunctions, $petsc_library),
        PetscErrorCode,
        (Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscCommDuplicate(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscCommDuplicate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{MPI_Comm}, Ptr{Cint}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscCommDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscCommDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{MPI_Comm},),
        arg1,
    )
end

@for_petsc function PetscMallocSetCoalesce(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscMallocSetCoalesce, $petsc_library),
        PetscErrorCode,
        (PetscBool,),
        arg1,
    )
end

@for_petsc function PetscMallocSet(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscMallocSet, $petsc_library),
        PetscErrorCode,
        (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscMallocClear(::$UnionPetscLib)
    @chk ccall((:PetscMallocClear, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscMallocSetDRAM(::$UnionPetscLib)
    @chk ccall((:PetscMallocSetDRAM, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscMallocResetDRAM(::$UnionPetscLib)
    @chk ccall((:PetscMallocResetDRAM, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscMallocDump(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscMallocDump, $petsc_library),
        PetscErrorCode,
        (Ptr{Libc.FILE},),
        arg1,
    )
end

@for_petsc function PetscMallocView(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscMallocView, $petsc_library),
        PetscErrorCode,
        (Ptr{Libc.FILE},),
        arg1,
    )
end

@for_petsc function PetscMallocGetCurrentUsage(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscMallocGetCurrentUsage, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscLogDouble},),
        arg1,
    )
end

@for_petsc function PetscMallocGetMaximumUsage(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscMallocGetMaximumUsage, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscLogDouble},),
        arg1,
    )
end

@for_petsc function PetscMallocPushMaximumUsage(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscMallocPushMaximumUsage, $petsc_library),
        PetscErrorCode,
        (Cint,),
        arg1,
    )
end

@for_petsc function PetscMallocPopMaximumUsage(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscMallocPopMaximumUsage, $petsc_library),
        PetscErrorCode,
        (Cint, Ptr{PetscLogDouble}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscMallocSetDebug(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscMallocSetDebug, $petsc_library),
        PetscErrorCode,
        (PetscBool, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscMallocGetDebug(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscMallocGetDebug, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscBool}, Ptr{PetscBool}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscMallocValidate(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscMallocValidate, $petsc_library),
        PetscErrorCode,
        (Cint, Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscMallocViewSet(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscMallocViewSet, $petsc_library),
        PetscErrorCode,
        (PetscLogDouble,),
        arg1,
    )
end

@for_petsc function PetscMallocViewGet(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscMallocViewGet, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscBool},),
        arg1,
    )
end

@for_petsc function PetscMallocLogRequestedSizeSet(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscMallocLogRequestedSizeSet, $petsc_library),
        PetscErrorCode,
        (PetscBool,),
        arg1,
    )
end

@for_petsc function PetscMallocLogRequestedSizeGet(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscMallocLogRequestedSizeGet, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscBool},),
        arg1,
    )
end

@for_petsc function PetscDataTypeToMPIDataType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDataTypeToMPIDataType, $petsc_library),
        PetscErrorCode,
        (PetscDataType, Ptr{MPI_Datatype}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscMPIDataTypeToPetscDataType(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscMPIDataTypeToPetscDataType, $petsc_library),
        PetscErrorCode,
        (MPI_Datatype, Ptr{PetscDataType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDataTypeGetSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDataTypeGetSize, $petsc_library),
        PetscErrorCode,
        (PetscDataType, Ptr{Csize_t}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDataTypeFromString(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDataTypeFromString, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{PetscDataType}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscMemcmp(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscMemcmp, $petsc_library),
        PetscErrorCode,
        (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscStrlen(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscStrlen, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Csize_t}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscStrToArray(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscStrToArray, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Cchar, Ptr{Cint}, Ptr{Ptr{Ptr{Cchar}}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscStrToArrayDestroy(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscStrToArrayDestroy, $petsc_library),
        PetscErrorCode,
        (Cint, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscStrcmp(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscStrcmp, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscStrgrt(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscStrgrt, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscStrcasecmp(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscStrcasecmp, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscStrncmp(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscStrncmp, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Csize_t, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscStrcpy(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscStrcpy, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscStrcat(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscStrcat, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscStrlcat(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscStrlcat, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscStrncpy(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscStrncpy, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscStrchr(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscStrchr, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Cchar, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscStrtolower(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscStrtolower, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar},),
        arg1,
    )
end

@for_petsc function PetscStrtoupper(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscStrtoupper, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar},),
        arg1,
    )
end

@for_petsc function PetscStrrchr(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscStrrchr, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Cchar, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscStrstr(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscStrstr, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscStrrstr(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscStrrstr, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscStrendswith(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscStrendswith, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscStrbeginswith(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscStrbeginswith, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscStrendswithwhich(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscStrendswithwhich, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Ptr{Cchar}}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscStrallocpy(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscStrallocpy, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscStrArrayallocpy(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscStrArrayallocpy, $petsc_library),
        PetscErrorCode,
        (Ptr{Ptr{Cchar}}, Ptr{Ptr{Ptr{Cchar}}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscStrArrayDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscStrArrayDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{Ptr{Ptr{Cchar}}},),
        arg1,
    )
end

@for_petsc function PetscStrNArrayallocpy(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscStrNArrayallocpy, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{Ptr{Cchar}}, Ptr{Ptr{Ptr{Cchar}}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscStrNArrayDestroy(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscStrNArrayDestroy, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{Ptr{Ptr{Cchar}}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscStrreplace(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscStrreplace, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscStrcmpNoError(::$UnionPetscLib, arg1, arg2, arg3)
    ccall(
        (:PetscStrcmpNoError, $petsc_library),
        Cvoid,
        (Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscTokenCreate(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscTokenCreate, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Cchar, Ptr{PetscToken}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscTokenFind(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscTokenFind, $petsc_library),
        PetscErrorCode,
        (PetscToken, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscTokenDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscTokenDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscToken},),
        arg1,
    )
end

@for_petsc function PetscStrInList(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscStrInList, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Cchar, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscEListFind(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscEListFind, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            Ptr{Ptr{Cchar}},
            Ptr{Cchar},
            Ptr{$PetscInt},
            Ptr{PetscBool},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscEnumFind(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscEnumFind, $petsc_library),
        PetscErrorCode,
        (Ptr{Ptr{Cchar}}, Ptr{Cchar}, Ptr{PetscEnum}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscMaxSum(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscMaxSum, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MPIULong_Send(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MPIULong_Send, $petsc_library),
        PetscErrorCode,
        (
            Ptr{Cvoid},
            $PetscInt,
            MPI_Datatype,
            PetscMPIInt,
            PetscMPIInt,
            MPI_Comm,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MPIULong_Recv(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MPIULong_Recv, $petsc_library),
        PetscErrorCode,
        (
            Ptr{Cvoid},
            $PetscInt,
            MPI_Datatype,
            PetscMPIInt,
            PetscMPIInt,
            MPI_Comm,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@enum PetscErrorType::UInt32 begin
    PETSC_ERROR_INITIAL = 0
    PETSC_ERROR_REPEAT = 1
    PETSC_ERROR_IN_CXX = 2
end

@for_petsc function PetscErrorPrintfInitialize(::$UnionPetscLib)
    @chk ccall(
        (:PetscErrorPrintfInitialize, $petsc_library),
        PetscErrorCode,
        (),
    )
end

@for_petsc function PetscErrorMessage(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscErrorMessage, $petsc_library),
        PetscErrorCode,
        (Cint, Ptr{Ptr{Cchar}}, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscTraceBackErrorHandler(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscTraceBackErrorHandler, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            Cint,
            Ptr{Cchar},
            Ptr{Cchar},
            PetscErrorCode,
            PetscErrorType,
            Ptr{Cchar},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function PetscIgnoreErrorHandler(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscIgnoreErrorHandler, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            Cint,
            Ptr{Cchar},
            Ptr{Cchar},
            PetscErrorCode,
            PetscErrorType,
            Ptr{Cchar},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function PetscEmacsClientErrorHandler(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscEmacsClientErrorHandler, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            Cint,
            Ptr{Cchar},
            Ptr{Cchar},
            PetscErrorCode,
            PetscErrorType,
            Ptr{Cchar},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function PetscMPIAbortErrorHandler(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscMPIAbortErrorHandler, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            Cint,
            Ptr{Cchar},
            Ptr{Cchar},
            PetscErrorCode,
            PetscErrorType,
            Ptr{Cchar},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function PetscAbortErrorHandler(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscAbortErrorHandler, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            Cint,
            Ptr{Cchar},
            Ptr{Cchar},
            PetscErrorCode,
            PetscErrorType,
            Ptr{Cchar},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function PetscAttachDebuggerErrorHandler(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscAttachDebuggerErrorHandler, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            Cint,
            Ptr{Cchar},
            Ptr{Cchar},
            PetscErrorCode,
            PetscErrorType,
            Ptr{Cchar},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function PetscReturnErrorHandler(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscReturnErrorHandler, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            Cint,
            Ptr{Cchar},
            Ptr{Cchar},
            PetscErrorCode,
            PetscErrorType,
            Ptr{Cchar},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function PetscPushErrorHandler(::$UnionPetscLib, handler, arg2)
    @chk ccall(
        (:PetscPushErrorHandler, $petsc_library),
        PetscErrorCode,
        (Ptr{Cvoid}, Ptr{Cvoid}),
        handler,
        arg2,
    )
end

@for_petsc function PetscPopErrorHandler(::$UnionPetscLib)
    @chk ccall((:PetscPopErrorHandler, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscSignalHandlerDefault(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSignalHandlerDefault, $petsc_library),
        PetscErrorCode,
        (Cint, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscPushSignalHandler(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscPushSignalHandler, $petsc_library),
        PetscErrorCode,
        (Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscPopSignalHandler(::$UnionPetscLib)
    @chk ccall((:PetscPopSignalHandler, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscCheckPointerSetIntensity(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscCheckPointerSetIntensity, $petsc_library),
        PetscErrorCode,
        ($PetscInt,),
        arg1,
    )
end

@for_petsc function PetscSignalSegvCheckPointerOrMpi(::$UnionPetscLib)
    ccall((:PetscSignalSegvCheckPointerOrMpi, $petsc_library), Cvoid, ())
end

@for_petsc function PetscSignalSegvCheckPointer(::$UnionPetscLib)
    ccall((:PetscSignalSegvCheckPointer, $petsc_library), Cvoid, ())
end

@enum PetscFPTrap::UInt32 begin
    PETSC_FP_TRAP_OFF = 0
    PETSC_FP_TRAP_ON = 1
end

@for_petsc function PetscSetFPTrap(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscSetFPTrap, $petsc_library),
        PetscErrorCode,
        (PetscFPTrap,),
        arg1,
    )
end

@for_petsc function PetscFPTrapPush(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscFPTrapPush, $petsc_library),
        PetscErrorCode,
        (PetscFPTrap,),
        arg1,
    )
end

@for_petsc function PetscFPTrapPop(::$UnionPetscLib)
    @chk ccall((:PetscFPTrapPop, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscDetermineInitialFPTrap(::$UnionPetscLib)
    @chk ccall(
        (:PetscDetermineInitialFPTrap, $petsc_library),
        PetscErrorCode,
        (),
    )
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

@for_petsc function PetscStackCopy(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscStackCopy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscStack}, Ptr{PetscStack}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscStackPrint(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscStackPrint, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscStack}, Ptr{Libc.FILE}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscStackActive(::$UnionPetscLib)
    ccall((:PetscStackActive, $petsc_library), PetscBool, ())
end

@for_petsc function PetscStackCreate(::$UnionPetscLib)
    @chk ccall((:PetscStackCreate, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscStackView(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscStackView, $petsc_library),
        PetscErrorCode,
        (Ptr{Libc.FILE},),
        arg1,
    )
end

@for_petsc function PetscStackDestroy(::$UnionPetscLib)
    @chk ccall((:PetscStackDestroy, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscClassIdRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscClassIdRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{PetscClassId}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectGetId(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectGetId, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{PetscObjectId}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectCompareId(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscObjectCompareId, $petsc_library),
        PetscErrorCode,
        (PetscObject, PetscObjectId, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscMemoryGetCurrentUsage(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscMemoryGetCurrentUsage, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscLogDouble},),
        arg1,
    )
end

@for_petsc function PetscMemoryGetMaximumUsage(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscMemoryGetMaximumUsage, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscLogDouble},),
        arg1,
    )
end

@for_petsc function PetscMemorySetGetMaximumUsage(::$UnionPetscLib)
    @chk ccall(
        (:PetscMemorySetGetMaximumUsage, $petsc_library),
        PetscErrorCode,
        (),
    )
end

@for_petsc function PetscMemoryTrace(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscMemoryTrace, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar},),
        arg1,
    )
end

@for_petsc function PetscSleep(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscSleep, $petsc_library),
        PetscErrorCode,
        ($PetscReal,),
        arg1,
    )
end

@for_petsc function PetscInitialize(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscInitialize, $petsc_library),
        PetscErrorCode,
        (Ptr{Cint}, Ptr{Ptr{Ptr{Cchar}}}, Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscInitializeNoPointers(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscInitializeNoPointers, $petsc_library),
        PetscErrorCode,
        (Cint, Ptr{Ptr{Cchar}}, Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscInitializeNoArguments(::$UnionPetscLib)
    @chk ccall(
        (:PetscInitializeNoArguments, $petsc_library),
        PetscErrorCode,
        (),
    )
end

@for_petsc function PetscInitialized(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscInitialized, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscBool},),
        arg1,
    )
end

@for_petsc function PetscFinalized(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscFinalized, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscBool},),
        arg1,
    )
end

@for_petsc function PetscFinalize(::$UnionPetscLib)
    @chk ccall((:PetscFinalize, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscInitializeFortran(::$UnionPetscLib)
    @chk ccall((:PetscInitializeFortran, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscGetArgs(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscGetArgs, $petsc_library),
        PetscErrorCode,
        (Ptr{Cint}, Ptr{Ptr{Ptr{Cchar}}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscGetArguments(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscGetArguments, $petsc_library),
        PetscErrorCode,
        (Ptr{Ptr{Ptr{Cchar}}},),
        arg1,
    )
end

@for_petsc function PetscFreeArguments(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscFreeArguments, $petsc_library),
        PetscErrorCode,
        (Ptr{Ptr{Cchar}},),
        arg1,
    )
end

@for_petsc function PetscEnd(::$UnionPetscLib)
    @chk ccall((:PetscEnd, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscSysInitializePackage(::$UnionPetscLib)
    @chk ccall((:PetscSysInitializePackage, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscPythonInitialize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscPythonInitialize, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscPythonFinalize(::$UnionPetscLib)
    @chk ccall((:PetscPythonFinalize, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscPythonPrintError(::$UnionPetscLib)
    @chk ccall((:PetscPythonPrintError, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscPythonMonitorSet(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscPythonMonitorSet, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscMonitorCompare(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscMonitorCompare, $petsc_library),
        PetscErrorCode,
        (
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{PetscBool},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

# typedef void ( * * PetscVoidStarFunction ) ( void )
const PetscVoidStarFunction = Ptr{Ptr{Cvoid}}

# typedef void ( * PetscVoidFunction ) ( void )
const PetscVoidFunction = Ptr{Cvoid}

# typedef PetscErrorCode ( * PetscErrorCodeFunction ) ( void )
const PetscErrorCodeFunction = Ptr{Cvoid}

@for_petsc function PetscObjectDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscObjectDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscObject},),
        arg1,
    )
end

@for_petsc function PetscObjectGetComm(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectGetComm, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{MPI_Comm}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectGetClassId(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectGetClassId, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{PetscClassId}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectGetClassName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectGetClassName, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectSetType, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectGetType, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectSetName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectSetName, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectGetName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectGetName, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectSetTabLevel(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectSetTabLevel, $petsc_library),
        PetscErrorCode,
        (PetscObject, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectGetTabLevel(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectGetTabLevel, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectIncrementTabLevel(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscObjectIncrementTabLevel, $petsc_library),
        PetscErrorCode,
        (PetscObject, PetscObject, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscObjectReference(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscObjectReference, $petsc_library),
        PetscErrorCode,
        (PetscObject,),
        arg1,
    )
end

@for_petsc function PetscObjectGetReference(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectGetReference, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectDereference(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscObjectDereference, $petsc_library),
        PetscErrorCode,
        (PetscObject,),
        arg1,
    )
end

@for_petsc function PetscObjectGetNewTag(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectGetNewTag, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{PetscMPIInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectCompose(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscObjectCompose, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{Cchar}, PetscObject),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscObjectRemoveReference(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectRemoveReference, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectQuery(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscObjectQuery, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{Cchar}, Ptr{PetscObject}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscObjectSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscObjectSetFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscObject,),
        arg1,
    )
end

@for_petsc function PetscObjectSetUp(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscObjectSetUp, $petsc_library),
        PetscErrorCode,
        (PetscObject,),
        arg1,
    )
end

@for_petsc function PetscObjectSetPrintedOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscObjectSetPrintedOptions, $petsc_library),
        PetscErrorCode,
        (PetscObject,),
        arg1,
    )
end

@for_petsc function PetscObjectInheritPrintedOptions(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscObjectInheritPrintedOptions, $petsc_library),
        PetscErrorCode,
        (PetscObject, PetscObject),
        arg1,
        arg2,
    )
end

@for_petsc function PetscCommGetNewTag(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscCommGetNewTag, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscMPIInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscOptionsCreate(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscOptionsCreate, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscOptions},),
        arg1,
    )
end

@for_petsc function PetscOptionsPush(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscOptionsPush, $petsc_library),
        PetscErrorCode,
        (PetscOptions,),
        arg1,
    )
end

@for_petsc function PetscOptionsPop(::$UnionPetscLib)
    @chk ccall((:PetscOptionsPop, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscOptionsDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscOptionsDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscOptions},),
        arg1,
    )
end

@for_petsc function PetscOptionsCreateDefault(::$UnionPetscLib)
    @chk ccall((:PetscOptionsCreateDefault, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscOptionsDestroyDefault(::$UnionPetscLib)
    @chk ccall(
        (:PetscOptionsDestroyDefault, $petsc_library),
        PetscErrorCode,
        (),
    )
end

@for_petsc function PetscOptionsHasHelp(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscOptionsHasHelp, $petsc_library),
        PetscErrorCode,
        (PetscOptions, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscOptionsHasName(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscOptionsHasName, $petsc_library),
        PetscErrorCode,
        (PetscOptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscOptionsGetBool(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscOptionsGetBool, $petsc_library),
        PetscErrorCode,
        (PetscOptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscOptionsGetInt(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscOptionsGetInt, $petsc_library),
        PetscErrorCode,
        (PetscOptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{$PetscInt}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscOptionsGetEnum(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscOptionsGetEnum, $petsc_library),
        PetscErrorCode,
        (
            PetscOptions,
            Ptr{Cchar},
            Ptr{Cchar},
            Ptr{Ptr{Cchar}},
            Ptr{PetscEnum},
            Ptr{PetscBool},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscOptionsGetEList(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscOptionsGetEList, $petsc_library),
        PetscErrorCode,
        (
            PetscOptions,
            Ptr{Cchar},
            Ptr{Cchar},
            Ptr{Ptr{Cchar}},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{PetscBool},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscOptionsGetReal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscOptionsGetReal, $petsc_library),
        PetscErrorCode,
        (PetscOptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{$PetscReal}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscOptionsGetScalar(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscOptionsGetScalar, $petsc_library),
        PetscErrorCode,
        (
            PetscOptions,
            Ptr{Cchar},
            Ptr{Cchar},
            Ptr{$PetscScalar},
            Ptr{PetscBool},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscOptionsGetString(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscOptionsGetString, $petsc_library),
        PetscErrorCode,
        (
            PetscOptions,
            Ptr{Cchar},
            Ptr{Cchar},
            Ptr{Cchar},
            Csize_t,
            Ptr{PetscBool},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscOptionsGetBoolArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscOptionsGetBoolArray, $petsc_library),
        PetscErrorCode,
        (
            PetscOptions,
            Ptr{Cchar},
            Ptr{Cchar},
            Ptr{PetscBool},
            Ptr{$PetscInt},
            Ptr{PetscBool},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscOptionsGetEnumArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscOptionsGetEnumArray, $petsc_library),
        PetscErrorCode,
        (
            PetscOptions,
            Ptr{Cchar},
            Ptr{Cchar},
            Ptr{Ptr{Cchar}},
            Ptr{PetscEnum},
            Ptr{$PetscInt},
            Ptr{PetscBool},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscOptionsGetIntArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscOptionsGetIntArray, $petsc_library),
        PetscErrorCode,
        (
            PetscOptions,
            Ptr{Cchar},
            Ptr{Cchar},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{PetscBool},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscOptionsGetRealArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscOptionsGetRealArray, $petsc_library),
        PetscErrorCode,
        (
            PetscOptions,
            Ptr{Cchar},
            Ptr{Cchar},
            Ptr{$PetscReal},
            Ptr{$PetscInt},
            Ptr{PetscBool},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscOptionsGetScalarArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscOptionsGetScalarArray, $petsc_library),
        PetscErrorCode,
        (
            PetscOptions,
            Ptr{Cchar},
            Ptr{Cchar},
            Ptr{$PetscScalar},
            Ptr{$PetscInt},
            Ptr{PetscBool},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscOptionsGetStringArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscOptionsGetStringArray, $petsc_library),
        PetscErrorCode,
        (
            PetscOptions,
            Ptr{Cchar},
            Ptr{Cchar},
            Ptr{Ptr{Cchar}},
            Ptr{$PetscInt},
            Ptr{PetscBool},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscOptionsValidKey(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscOptionsValidKey, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscOptionsSetAlias(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscOptionsSetAlias, $petsc_library),
        PetscErrorCode,
        (PetscOptions, Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscOptionsSetValue(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscOptionsSetValue, $petsc_library),
        PetscErrorCode,
        (PetscOptions, Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscOptionsClearValue(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscOptionsClearValue, $petsc_library),
        PetscErrorCode,
        (PetscOptions, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscOptionsFindPair(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscOptionsFindPair, $petsc_library),
        PetscErrorCode,
        (PetscOptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{Ptr{Cchar}}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscOptionsGetAll(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscOptionsGetAll, $petsc_library),
        PetscErrorCode,
        (PetscOptions, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscOptionsAllUsed(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscOptionsAllUsed, $petsc_library),
        PetscErrorCode,
        (PetscOptions, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscOptionsUsed(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscOptionsUsed, $petsc_library),
        PetscErrorCode,
        (PetscOptions, Ptr{Cchar}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscOptionsLeft(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscOptionsLeft, $petsc_library),
        PetscErrorCode,
        (PetscOptions,),
        arg1,
    )
end

@for_petsc function PetscOptionsLeftGet(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscOptionsLeftGet, $petsc_library),
        PetscErrorCode,
        (
            PetscOptions,
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cchar}}},
            Ptr{Ptr{Ptr{Cchar}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscOptionsLeftRestore(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscOptionsLeftRestore, $petsc_library),
        PetscErrorCode,
        (
            PetscOptions,
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cchar}}},
            Ptr{Ptr{Ptr{Cchar}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscOptionsView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscOptionsView, $petsc_library),
        PetscErrorCode,
        (PetscOptions, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscOptionsReject(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscOptionsReject, $petsc_library),
        PetscErrorCode,
        (PetscOptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscOptionsInsert(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscOptionsInsert, $petsc_library),
        PetscErrorCode,
        (PetscOptions, Ptr{Cint}, Ptr{Ptr{Ptr{Cchar}}}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscOptionsInsertFile(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscOptionsInsertFile, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, PetscOptions, Ptr{Cchar}, PetscBool),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscOptionsInsertFileYAML(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscOptionsInsertFileYAML, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, PetscOptions, Ptr{Cchar}, PetscBool),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscOptionsInsertString(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscOptionsInsertString, $petsc_library),
        PetscErrorCode,
        (PetscOptions, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscOptionsInsertStringYAML(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscOptionsInsertStringYAML, $petsc_library),
        PetscErrorCode,
        (PetscOptions, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscOptionsInsertArgs(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscOptionsInsertArgs, $petsc_library),
        PetscErrorCode,
        (PetscOptions, Cint, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscOptionsClear(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscOptionsClear, $petsc_library),
        PetscErrorCode,
        (PetscOptions,),
        arg1,
    )
end

@for_petsc function PetscOptionsPrefixPush(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscOptionsPrefixPush, $petsc_library),
        PetscErrorCode,
        (PetscOptions, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscOptionsPrefixPop(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscOptionsPrefixPop, $petsc_library),
        PetscErrorCode,
        (PetscOptions,),
        arg1,
    )
end

@for_petsc function PetscOptionsGetenv(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscOptionsGetenv, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Csize_t, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscOptionsStringToBool(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscOptionsStringToBool, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscOptionsStringToInt(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscOptionsStringToInt, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscOptionsStringToReal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscOptionsStringToReal, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscOptionsStringToScalar(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscOptionsStringToScalar, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{$PetscScalar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscOptionsMonitorSet(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscOptionsMonitorSet, $petsc_library),
        PetscErrorCode,
        (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscOptionsMonitorDefault(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscOptionsMonitorDefault, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscObjectSetOptions(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectSetOptions, $petsc_library),
        PetscErrorCode,
        (PetscObject, PetscOptions),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectGetOptions(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectGetOptions, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{PetscOptions}),
        arg1,
        arg2,
    )
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

mutable struct __JL__n_PetscOptionItem end

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

@for_petsc function PetscOptionsSAWsDestroy(::$UnionPetscLib)
    @chk ccall((:PetscOptionsSAWsDestroy, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscObjectAddOptionsHandler(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscObjectAddOptionsHandler, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscObjectDestroyOptionsHandlers(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscObjectDestroyOptionsHandlers, $petsc_library),
        PetscErrorCode,
        (PetscObject,),
        arg1,
    )
end

@for_petsc function PetscMallocTraceSet(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscMallocTraceSet, $petsc_library),
        PetscErrorCode,
        (PetscViewer, PetscBool, PetscLogDouble),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscMallocTraceGet(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscMallocTraceGet, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscBool},),
        arg1,
    )
end

@for_petsc function PetscObjectsListGetGlobalNumbering(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscObjectsListGetGlobalNumbering, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, Ptr{PetscObject}, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscMemoryView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscMemoryView, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectPrintClassNamePrefixType(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscObjectPrintClassNamePrefixType, $petsc_library),
        PetscErrorCode,
        (PetscObject, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectView, $petsc_library),
        PetscErrorCode,
        (PetscObject, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectSetOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectSetOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectAppendOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectAppendOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectPrependOptionsPrefix(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscObjectPrependOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectGetOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectGetOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectChangeTypeName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectChangeTypeName, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectRegisterDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscObjectRegisterDestroy, $petsc_library),
        PetscErrorCode,
        (PetscObject,),
        arg1,
    )
end

@for_petsc function PetscObjectRegisterDestroyAll(::$UnionPetscLib)
    @chk ccall(
        (:PetscObjectRegisterDestroyAll, $petsc_library),
        PetscErrorCode,
        (),
    )
end

@for_petsc function PetscObjectViewFromOptions(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscObjectViewFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscObject, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscObjectName(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscObjectName, $petsc_library),
        PetscErrorCode,
        (PetscObject,),
        arg1,
    )
end

@for_petsc function PetscObjectTypeCompare(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscObjectTypeCompare, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{Cchar}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscObjectBaseTypeCompare(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscObjectBaseTypeCompare, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{Cchar}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscRegisterFinalize(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscRegisterFinalize, $petsc_library),
        PetscErrorCode,
        (Ptr{Cvoid},),
        arg1,
    )
end

@for_petsc function PetscRegisterFinalizeAll(::$UnionPetscLib)
    @chk ccall((:PetscRegisterFinalizeAll, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscDLOpen(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDLOpen, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, PetscDLMode, Ptr{PetscDLHandle}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDLClose(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDLClose, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscDLHandle},),
        arg1,
    )
end

@for_petsc function PetscDLSym(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDLSym, $petsc_library),
        PetscErrorCode,
        (PetscDLHandle, Ptr{Cchar}, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDLAddr(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDLAddr, $petsc_library),
        PetscErrorCode,
        (Ptr{Cvoid}, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectsDump(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectsDump, $petsc_library),
        PetscErrorCode,
        (Ptr{Libc.FILE}, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectListDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscObjectListDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscObjectList},),
        arg1,
    )
end

@for_petsc function PetscObjectListFind(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscObjectListFind, $petsc_library),
        PetscErrorCode,
        (PetscObjectList, Ptr{Cchar}, Ptr{PetscObject}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscObjectListReverseFind(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscObjectListReverseFind, $petsc_library),
        PetscErrorCode,
        (PetscObjectList, PetscObject, Ptr{Ptr{Cchar}}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscObjectListAdd(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscObjectListAdd, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscObjectList}, Ptr{Cchar}, PetscObject),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscObjectListRemoveReference(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectListRemoveReference, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscObjectList}, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscObjectListDuplicate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectListDuplicate, $petsc_library),
        PetscErrorCode,
        (PetscObjectList, Ptr{PetscObjectList}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFunctionListDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscFunctionListDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscFunctionList},),
        arg1,
    )
end

@for_petsc function PetscFunctionListPrintTypes(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:PetscFunctionListPrintTypes, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            Ptr{Libc.FILE},
            Ptr{Cchar},
            Ptr{Cchar},
            Ptr{Cchar},
            Ptr{Cchar},
            PetscFunctionList,
            Ptr{Cchar},
            Ptr{Cchar},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function PetscFunctionListDuplicate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFunctionListDuplicate, $petsc_library),
        PetscErrorCode,
        (PetscFunctionList, Ptr{PetscFunctionList}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFunctionListView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFunctionListView, $petsc_library),
        PetscErrorCode,
        (PetscFunctionList, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFunctionListGet(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscFunctionListGet, $petsc_library),
        PetscErrorCode,
        (PetscFunctionList, Ptr{Ptr{Ptr{Cchar}}}, Ptr{Cint}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDLLibraryAppend(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDLLibraryAppend, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscDLLibrary}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDLLibraryPrepend(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDLLibraryPrepend, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscDLLibrary}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDLLibrarySym(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDLLibrarySym, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            Ptr{PetscDLLibrary},
            Ptr{Cchar},
            Ptr{Cchar},
            Ptr{Ptr{Cvoid}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDLLibraryPrintPath(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDLLibraryPrintPath, $petsc_library),
        PetscErrorCode,
        (PetscDLLibrary,),
        arg1,
    )
end

@for_petsc function PetscDLLibraryRetrieve(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDLLibraryRetrieve, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Csize_t, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDLLibraryOpen(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDLLibraryOpen, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{PetscDLLibrary}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDLLibraryClose(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDLLibraryClose, $petsc_library),
        PetscErrorCode,
        (PetscDLLibrary,),
        arg1,
    )
end

@for_petsc function PetscSplitOwnership(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSplitOwnership, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSplitOwnershipBlock(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSplitOwnershipBlock, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSplitOwnershipEqual(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSplitOwnershipEqual, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSequentialPhaseBegin(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSequentialPhaseBegin, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, PetscMPIInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSequentialPhaseEnd(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSequentialPhaseEnd, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, PetscMPIInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscBarrier(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscBarrier, $petsc_library),
        PetscErrorCode,
        (PetscObject,),
        arg1,
    )
end

@for_petsc function PetscMPIDump(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscMPIDump, $petsc_library),
        PetscErrorCode,
        (Ptr{Libc.FILE},),
        arg1,
    )
end

@for_petsc function PetscGlobalMinMaxInt(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscGlobalMinMaxInt, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscGlobalMinMaxReal(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscGlobalMinMaxReal, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscGetCPUTime(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscGetCPUTime, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscLogDouble},),
        arg1,
    )
end

@for_petsc function PetscTime(::$UnionPetscLib, v)
    @chk ccall(
        (:PetscTime, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscLogDouble},),
        v,
    )
end

@for_petsc function PetscTimeSubtract(::$UnionPetscLib, v)
    @chk ccall(
        (:PetscTimeSubtract, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscLogDouble},),
        v,
    )
end

@for_petsc function PetscTimeAdd(::$UnionPetscLib, v)
    @chk ccall(
        (:PetscTimeAdd, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscLogDouble},),
        v,
    )
end

@enum PetscInfoCommFlag::Int32 begin
    PETSC_INFO_COMM_ALL = -1
    PETSC_INFO_COMM_NO_SELF = 0
    PETSC_INFO_COMM_ONLY_SELF = 1
end

@for_petsc function PetscInfoDeactivateClass(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscInfoDeactivateClass, $petsc_library),
        PetscErrorCode,
        (PetscClassId,),
        arg1,
    )
end

@for_petsc function PetscInfoActivateClass(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscInfoActivateClass, $petsc_library),
        PetscErrorCode,
        (PetscClassId,),
        arg1,
    )
end

@for_petsc function PetscInfoEnabled(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscInfoEnabled, $petsc_library),
        PetscErrorCode,
        (PetscClassId, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscInfoAllow(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscInfoAllow, $petsc_library),
        PetscErrorCode,
        (PetscBool,),
        arg1,
    )
end

@for_petsc function PetscInfoSetFile(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscInfoSetFile, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscInfoGetFile(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscInfoGetFile, $petsc_library),
        PetscErrorCode,
        (Ptr{Ptr{Cchar}}, Ptr{Ptr{Libc.FILE}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscInfoSetClasses(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscInfoSetClasses, $petsc_library),
        PetscErrorCode,
        (PetscBool, $PetscInt, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscInfoGetClass(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscInfoGetClass, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscInfoGetInfo(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscInfoGetInfo, $petsc_library),
        PetscErrorCode,
        (
            Ptr{PetscBool},
            Ptr{PetscBool},
            Ptr{PetscBool},
            Ptr{PetscBool},
            Ptr{PetscInfoCommFlag},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscInfoProcessClass(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscInfoProcessClass, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, $PetscInt, Ptr{PetscClassId}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscInfoSetFilterCommSelf(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscInfoSetFilterCommSelf, $petsc_library),
        PetscErrorCode,
        (PetscInfoCommFlag,),
        arg1,
    )
end

@for_petsc function PetscInfoSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscInfoSetFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscOptions,),
        arg1,
    )
end

@for_petsc function PetscInfoDestroy(::$UnionPetscLib)
    @chk ccall((:PetscInfoDestroy, $petsc_library), PetscErrorCode, ())
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

@for_petsc function PetscLogObjectParent(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogObjectParent, $petsc_library),
        PetscErrorCode,
        (PetscObject, PetscObject),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogObjectMemory(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogObjectMemory, $petsc_library),
        PetscErrorCode,
        (PetscObject, PetscLogDouble),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogGetStageLog(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogGetStageLog, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscStageLog},),
        arg1,
    )
end

@for_petsc function PetscStageLogGetCurrent(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscStageLogGetCurrent, $petsc_library),
        PetscErrorCode,
        (PetscStageLog, Ptr{Cint}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscStageLogGetEventPerfLog(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscStageLogGetEventPerfLog, $petsc_library),
        PetscErrorCode,
        (PetscStageLog, Cint, Ptr{PetscEventPerfLog}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogFlops(::$UnionPetscLib, n)
    @chk ccall(
        (:PetscLogFlops, $petsc_library),
        PetscErrorCode,
        (PetscLogDouble,),
        n,
    )
end

@for_petsc function PetscGetFlops(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscGetFlops, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscLogDouble},),
        arg1,
    )
end

@for_petsc function PetscLogDefaultBegin(::$UnionPetscLib)
    @chk ccall((:PetscLogDefaultBegin, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscLogAllBegin(::$UnionPetscLib)
    @chk ccall((:PetscLogAllBegin, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscLogNestedBegin(::$UnionPetscLib)
    @chk ccall((:PetscLogNestedBegin, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscLogTraceBegin(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogTraceBegin, $petsc_library),
        PetscErrorCode,
        (Ptr{Libc.FILE},),
        arg1,
    )
end

@for_petsc function PetscLogActions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogActions, $petsc_library),
        PetscErrorCode,
        (PetscBool,),
        arg1,
    )
end

@for_petsc function PetscLogObjects(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogObjects, $petsc_library),
        PetscErrorCode,
        (PetscBool,),
        arg1,
    )
end

@for_petsc function PetscLogSetThreshold(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogSetThreshold, $petsc_library),
        PetscErrorCode,
        (PetscLogDouble, Ptr{PetscLogDouble}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogSet(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogSet, $petsc_library),
        PetscErrorCode,
        (Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogView(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogView, $petsc_library),
        PetscErrorCode,
        (PetscViewer,),
        arg1,
    )
end

@for_petsc function PetscLogViewFromOptions(::$UnionPetscLib)
    @chk ccall((:PetscLogViewFromOptions, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscLogDump(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogDump, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar},),
        arg1,
    )
end

@for_petsc function PetscLogStageRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogStageRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{PetscLogStage}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogStagePush(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogStagePush, $petsc_library),
        PetscErrorCode,
        (PetscLogStage,),
        arg1,
    )
end

@for_petsc function PetscLogStagePop(::$UnionPetscLib)
    @chk ccall((:PetscLogStagePop, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscLogStageSetActive(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogStageSetActive, $petsc_library),
        PetscErrorCode,
        (PetscLogStage, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogStageGetActive(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogStageGetActive, $petsc_library),
        PetscErrorCode,
        (PetscLogStage, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogStageSetVisible(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogStageSetVisible, $petsc_library),
        PetscErrorCode,
        (PetscLogStage, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogStageGetVisible(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogStageGetVisible, $petsc_library),
        PetscErrorCode,
        (PetscLogStage, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogStageGetId(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogStageGetId, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{PetscLogStage}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogEventRegister(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscLogEventRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, PetscClassId, Ptr{PetscLogEvent}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogEventSetCollective(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogEventSetCollective, $petsc_library),
        PetscErrorCode,
        (PetscLogEvent, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogEventIncludeClass(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogEventIncludeClass, $petsc_library),
        PetscErrorCode,
        (PetscClassId,),
        arg1,
    )
end

@for_petsc function PetscLogEventExcludeClass(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogEventExcludeClass, $petsc_library),
        PetscErrorCode,
        (PetscClassId,),
        arg1,
    )
end

@for_petsc function PetscLogEventActivate(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogEventActivate, $petsc_library),
        PetscErrorCode,
        (PetscLogEvent,),
        arg1,
    )
end

@for_petsc function PetscLogEventDeactivate(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogEventDeactivate, $petsc_library),
        PetscErrorCode,
        (PetscLogEvent,),
        arg1,
    )
end

@for_petsc function PetscLogEventDeactivatePush(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogEventDeactivatePush, $petsc_library),
        PetscErrorCode,
        (PetscLogEvent,),
        arg1,
    )
end

@for_petsc function PetscLogEventDeactivatePop(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogEventDeactivatePop, $petsc_library),
        PetscErrorCode,
        (PetscLogEvent,),
        arg1,
    )
end

@for_petsc function PetscLogEventSetActiveAll(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogEventSetActiveAll, $petsc_library),
        PetscErrorCode,
        (PetscLogEvent, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogEventActivateClass(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogEventActivateClass, $petsc_library),
        PetscErrorCode,
        (PetscClassId,),
        arg1,
    )
end

@for_petsc function PetscLogEventDeactivateClass(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogEventDeactivateClass, $petsc_library),
        PetscErrorCode,
        (PetscClassId,),
        arg1,
    )
end

@for_petsc function PetscLogEventGetId(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogEventGetId, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{PetscLogEvent}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogEventGetPerfInfo(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscLogEventGetPerfInfo, $petsc_library),
        PetscErrorCode,
        (Cint, PetscLogEvent, Ptr{PetscEventPerfInfo}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogEventSetDof(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscLogEventSetDof, $petsc_library),
        PetscErrorCode,
        (PetscLogEvent, $PetscInt, PetscLogDouble),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogEventSetError(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscLogEventSetError, $petsc_library),
        PetscErrorCode,
        (PetscLogEvent, $PetscInt, PetscLogDouble),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogEventSynchronize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogEventSynchronize, $petsc_library),
        PetscErrorCode,
        (PetscLogEvent, MPI_Comm),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogEventGetFlops(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogEventGetFlops, $petsc_library),
        PetscErrorCode,
        (PetscLogEvent, Ptr{PetscLogDouble}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogEventZeroFlops(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogEventZeroFlops, $petsc_library),
        PetscErrorCode,
        (PetscLogEvent,),
        arg1,
    )
end

@for_petsc function PetscMPITypeSize(::$UnionPetscLib, count, type, length)
    @chk ccall(
        (:PetscMPITypeSize, $petsc_library),
        PetscErrorCode,
        ($PetscInt, MPI_Datatype, Ptr{PetscLogDouble}),
        count,
        type,
        length,
    )
end

@for_petsc function PetscMPITypeSizeComm(
    ::$UnionPetscLib,
    comm,
    counts,
    type,
    length,
)
    @chk ccall(
        (:PetscMPITypeSizeComm, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscMPIInt}, MPI_Datatype, Ptr{PetscLogDouble}),
        comm,
        counts,
        type,
        length,
    )
end

@for_petsc function PetscMPITypeSizeCount(
    ::$UnionPetscLib,
    n,
    counts,
    type,
    length,
)
    @chk ccall(
        (:PetscMPITypeSizeCount, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{PetscMPIInt}, MPI_Datatype, Ptr{PetscLogDouble}),
        n,
        counts,
        type,
        length,
    )
end

@for_petsc function PetscMPIParallelComm(::$UnionPetscLib, comm)
    ccall((:PetscMPIParallelComm, $petsc_library), Cint, (MPI_Comm,), comm)
end

@for_petsc function PetscFixFilename(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFixFilename, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFOpen(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscFOpen, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Ptr{Ptr{Libc.FILE}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscFClose(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFClose, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Libc.FILE}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFormatRealArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscFormatRealArray, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Csize_t, Ptr{Cchar}, $PetscInt, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscFormatConvertGetSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFormatConvertGetSize, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Csize_t}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFormatConvert(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFormatConvert, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscPOpen(::$UnionPetscLib, arg1, arg2, arg3, arg4, arg5)
    @chk ccall(
        (:PetscPOpen, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Ptr{Libc.FILE}}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscPClose(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscPClose, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Libc.FILE}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscPOpenSetMachine(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscPOpenSetMachine, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar},),
        arg1,
    )
end

@for_petsc function PetscSynchronizedFlush(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSynchronizedFlush, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Libc.FILE}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSynchronizedFGets(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSynchronizedFGets, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Libc.FILE}, Csize_t, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscStartMatlab(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscStartMatlab, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Ptr{Ptr{Libc.FILE}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscStartJava(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscStartJava, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Ptr{Ptr{Libc.FILE}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscGetPetscDir(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscGetPetscDir, $petsc_library),
        PetscErrorCode,
        (Ptr{Ptr{Cchar}},),
        arg1,
    )
end

@for_petsc function PetscContainerGetPointer(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscContainerGetPointer, $petsc_library),
        PetscErrorCode,
        (PetscContainer, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscContainerSetPointer(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscContainerSetPointer, $petsc_library),
        PetscErrorCode,
        (PetscContainer, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscContainerDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscContainerDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscContainer},),
        arg1,
    )
end

@for_petsc function PetscContainerCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscContainerCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscContainer}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscContainerSetUserDestroy(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscContainerSetUserDestroy, $petsc_library),
        PetscErrorCode,
        (PetscContainer, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscContainerUserDestroyDefault(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscContainerUserDestroyDefault, $petsc_library),
        PetscErrorCode,
        (Ptr{Cvoid},),
        arg1,
    )
end

@for_petsc function PetscIntView(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscIntView, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscInt}, PetscViewer),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscRealView(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscRealView, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscReal}, PetscViewer),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscScalarView(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscScalarView, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscScalar}, PetscViewer),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscMemmove(::$UnionPetscLib, a, b, n)
    @chk ccall(
        (:PetscMemmove, $petsc_library),
        PetscErrorCode,
        (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t),
        a,
        b,
        n,
    )
end

@for_petsc function PetscMemcpy(::$UnionPetscLib, a, b, n)
    @chk ccall(
        (:PetscMemcpy, $petsc_library),
        PetscErrorCode,
        (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t),
        a,
        b,
        n,
    )
end

@for_petsc function PetscMemzero(::$UnionPetscLib, a, n)
    @chk ccall(
        (:PetscMemzero, $petsc_library),
        PetscErrorCode,
        (Ptr{Cvoid}, Csize_t),
        a,
        n,
    )
end

@for_petsc function PetscIntCast(::$UnionPetscLib, a, b)
    @chk ccall(
        (:PetscIntCast, $petsc_library),
        PetscErrorCode,
        (PetscInt64, Ptr{$PetscInt}),
        a,
        b,
    )
end

@for_petsc function PetscBLASIntCast(::$UnionPetscLib, a, b)
    @chk ccall(
        (:PetscBLASIntCast, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{PetscBLASInt}),
        a,
        b,
    )
end

@for_petsc function PetscMPIIntCast(::$UnionPetscLib, a, b)
    @chk ccall(
        (:PetscMPIIntCast, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{PetscMPIInt}),
        a,
        b,
    )
end

@for_petsc function PetscRealIntMultTruncate(::$UnionPetscLib, a, b)
    ccall(
        (:PetscRealIntMultTruncate, $petsc_library),
        $PetscInt,
        ($PetscReal, $PetscInt),
        a,
        b,
    )
end

@for_petsc function PetscIntMultTruncate(::$UnionPetscLib, a, b)
    ccall(
        (:PetscIntMultTruncate, $petsc_library),
        $PetscInt,
        ($PetscInt, $PetscInt),
        a,
        b,
    )
end

@for_petsc function PetscIntSumTruncate(::$UnionPetscLib, a, b)
    ccall(
        (:PetscIntSumTruncate, $petsc_library),
        $PetscInt,
        ($PetscInt, $PetscInt),
        a,
        b,
    )
end

@for_petsc function PetscIntMultError(::$UnionPetscLib, a, b, result)
    @chk ccall(
        (:PetscIntMultError, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, Ptr{$PetscInt}),
        a,
        b,
        result,
    )
end

@for_petsc function PetscIntSumError(::$UnionPetscLib, a, b, result)
    @chk ccall(
        (:PetscIntSumError, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, Ptr{$PetscInt}),
        a,
        b,
        result,
    )
end

@for_petsc function PetscGetArchType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscGetArchType, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
    )
end

@for_petsc function PetscGetHostName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscGetHostName, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
    )
end

@for_petsc function PetscGetUserName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscGetUserName, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
    )
end

@for_petsc function PetscGetProgramName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscGetProgramName, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSetProgramName(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscSetProgramName, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar},),
        arg1,
    )
end

@for_petsc function PetscGetDate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscGetDate, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
    )
end

@for_petsc function PetscGetVersion(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscGetVersion, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
    )
end

@for_petsc function PetscGetVersionNumber(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscGetVersionNumber, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSortedInt(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSortedInt, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscInt}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSortedMPIInt(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSortedMPIInt, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{PetscMPIInt}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSortedReal(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSortedReal, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscReal}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSortInt(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSortInt, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSortReverseInt(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSortReverseInt, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSortedRemoveDupsInt(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSortedRemoveDupsInt, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSortRemoveDupsInt(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSortRemoveDupsInt, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscCheckDupsInt(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscCheckDupsInt, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscInt}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscFindInt(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscFindInt, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscFindMPIInt(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscFindMPIInt, $petsc_library),
        PetscErrorCode,
        (PetscMPIInt, $PetscInt, Ptr{PetscMPIInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSortIntWithPermutation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSortIntWithPermutation, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSortStrWithPermutation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSortStrWithPermutation, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{Ptr{Cchar}}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSortIntWithArray(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSortIntWithArray, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSortIntWithArrayPair(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSortIntWithArrayPair, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSortMPIInt(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSortMPIInt, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{PetscMPIInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSortRemoveDupsMPIInt(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSortRemoveDupsMPIInt, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscInt}, Ptr{PetscMPIInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSortMPIIntWithArray(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSortMPIIntWithArray, $petsc_library),
        PetscErrorCode,
        (PetscMPIInt, Ptr{PetscMPIInt}, Ptr{PetscMPIInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSortMPIIntWithIntArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSortMPIIntWithIntArray, $petsc_library),
        PetscErrorCode,
        (PetscMPIInt, Ptr{PetscMPIInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSortIntWithScalarArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSortIntWithScalarArray, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSortIntWithDataArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscSortIntWithDataArray, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscInt}, Ptr{Cvoid}, Csize_t, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscSortReal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSortReal, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSortRealWithArrayInt(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSortRealWithArrayInt, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscReal}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSortRealWithPermutation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSortRealWithPermutation, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscReal}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSortRemoveDupsReal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSortRemoveDupsReal, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscInt}, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFindReal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscFindReal, $petsc_library),
        PetscErrorCode,
        ($PetscReal, $PetscInt, Ptr{$PetscReal}, $PetscReal, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscSortSplit(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscSortSplit, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSortSplitReal(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscSortSplitReal, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscProcessTree(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscProcessTree, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            Ptr{PetscBool},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscInt}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function PetscMergeIntArrayPair(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:PetscMergeIntArrayPair, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscInt}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function PetscMergeIntArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscMergeIntArray, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscMergeMPIIntArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscMergeMPIIntArray, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            Ptr{PetscMPIInt},
            $PetscInt,
            Ptr{PetscMPIInt},
            Ptr{$PetscInt},
            Ptr{Ptr{PetscMPIInt}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscParallelSortedInt(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscParallelSortedInt, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, Ptr{$PetscInt}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscTimSort(::$UnionPetscLib, arg1, arg2, arg3, arg4, arg5)
    @chk ccall(
        (:PetscTimSort, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{Cvoid}, Csize_t, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscIntSortSemiOrdered(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscIntSortSemiOrdered, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscMPIIntSortSemiOrdered(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscMPIIntSortSemiOrdered, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{PetscMPIInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscRealSortSemiOrdered(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscRealSortSemiOrdered, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscTimSortWithArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscTimSortWithArray, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            Ptr{Cvoid},
            Csize_t,
            Ptr{Cvoid},
            Csize_t,
            Ptr{Cvoid},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscIntSortSemiOrderedWithArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscIntSortSemiOrderedWithArray, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscMPIIntSortSemiOrderedWithArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscMPIIntSortSemiOrderedWithArray, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{PetscMPIInt}, Ptr{PetscMPIInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscRealSortSemiOrderedWithArrayInt(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscRealSortSemiOrderedWithArrayInt, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscReal}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSetDisplay(::$UnionPetscLib)
    @chk ccall((:PetscSetDisplay, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscGetDisplay(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscGetDisplay, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
    )
end

const PetscRandomType = Ptr{Cchar}

@for_petsc function PetscRandomInitializePackage(::$UnionPetscLib)
    @chk ccall(
        (:PetscRandomInitializePackage, $petsc_library),
        PetscErrorCode,
        (),
    )
end

@for_petsc function PetscRandomRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscRandomRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscRandomSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscRandomSetType, $petsc_library),
        PetscErrorCode,
        (PetscRandom, PetscRandomType),
        arg1,
        arg2,
    )
end

@for_petsc function PetscRandomSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscRandomSetFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscRandom,),
        arg1,
    )
end

@for_petsc function PetscRandomGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscRandomGetType, $petsc_library),
        PetscErrorCode,
        (PetscRandom, Ptr{PetscRandomType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscRandomViewFromOptions(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscRandomViewFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscRandom, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscRandomView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscRandomView, $petsc_library),
        PetscErrorCode,
        (PetscRandom, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscRandomCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscRandomCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscRandom}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscRandomGetValue(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscRandomGetValue, $petsc_library),
        PetscErrorCode,
        (PetscRandom, Ptr{$PetscScalar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscRandomGetValueReal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscRandomGetValueReal, $petsc_library),
        PetscErrorCode,
        (PetscRandom, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscRandomGetValues(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscRandomGetValues, $petsc_library),
        PetscErrorCode,
        (PetscRandom, $PetscInt, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscRandomGetValuesReal(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscRandomGetValuesReal, $petsc_library),
        PetscErrorCode,
        (PetscRandom, $PetscInt, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscRandomGetInterval(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscRandomGetInterval, $petsc_library),
        PetscErrorCode,
        (PetscRandom, Ptr{$PetscScalar}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscRandomSetInterval(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscRandomSetInterval, $petsc_library),
        PetscErrorCode,
        (PetscRandom, $PetscScalar, $PetscScalar),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscRandomSetSeed(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscRandomSetSeed, $petsc_library),
        PetscErrorCode,
        (PetscRandom, Culong),
        arg1,
        arg2,
    )
end

@for_petsc function PetscRandomGetSeed(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscRandomGetSeed, $petsc_library),
        PetscErrorCode,
        (PetscRandom, Ptr{Culong}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscRandomSeed(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscRandomSeed, $petsc_library),
        PetscErrorCode,
        (PetscRandom,),
        arg1,
    )
end

@for_petsc function PetscRandomDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscRandomDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscRandom},),
        arg1,
    )
end

@for_petsc function PetscGetFullPath(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscGetFullPath, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscGetRelativePath(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscGetRelativePath, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscGetWorkingDirectory(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscGetWorkingDirectory, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
    )
end

@for_petsc function PetscGetRealPath(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscGetRealPath, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscGetHomeDirectory(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscGetHomeDirectory, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
    )
end

@for_petsc function PetscTestFile(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscTestFile, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Cchar, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscTestDirectory(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscTestDirectory, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Cchar, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscMkdir(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscMkdir, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar},),
        arg1,
    )
end

@for_petsc function PetscMkdtemp(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscMkdtemp, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar},),
        arg1,
    )
end

@for_petsc function PetscRMTree(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscRMTree, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar},),
        arg1,
    )
end

@for_petsc function PetscBinaryBigEndian(::$UnionPetscLib)
    ccall((:PetscBinaryBigEndian, $petsc_library), PetscBool, ())
end

@for_petsc function PetscBinaryRead(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscBinaryRead, $petsc_library),
        PetscErrorCode,
        (Cint, Ptr{Cvoid}, $PetscInt, Ptr{$PetscInt}, PetscDataType),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscBinarySynchronizedRead(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscBinarySynchronizedRead, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Cint, Ptr{Cvoid}, $PetscInt, Ptr{$PetscInt}, PetscDataType),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscBinaryWrite(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscBinaryWrite, $petsc_library),
        PetscErrorCode,
        (Cint, Ptr{Cvoid}, $PetscInt, PetscDataType),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscBinarySynchronizedWrite(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscBinarySynchronizedWrite, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Cint, Ptr{Cvoid}, $PetscInt, PetscDataType),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscBinaryOpen(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscBinaryOpen, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, PetscFileMode, Ptr{Cint}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscBinaryClose(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscBinaryClose, $petsc_library),
        PetscErrorCode,
        (Cint,),
        arg1,
    )
end

@for_petsc function PetscSharedTmp(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSharedTmp, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSharedWorkingDirectory(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSharedWorkingDirectory, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscGetTmp(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscGetTmp, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscFileRetrieve(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscFileRetrieve, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Csize_t, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscLs(::$UnionPetscLib, arg1, arg2, arg3, arg4, arg5)
    @chk ccall(
        (:PetscLs, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Csize_t, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscOpenSocket(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscOpenSocket, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Cint, Ptr{Cint}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscBinarySeek(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscBinarySeek, $petsc_library),
        PetscErrorCode,
        (Cint, off_t, PetscBinarySeekType, Ptr{off_t}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscBinarySynchronizedSeek(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscBinarySynchronizedSeek, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Cint, off_t, PetscBinarySeekType, Ptr{off_t}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscByteSwap(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscByteSwap, $petsc_library),
        PetscErrorCode,
        (Ptr{Cvoid}, PetscDataType, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSetDebugTerminal(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscSetDebugTerminal, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar},),
        arg1,
    )
end

@for_petsc function PetscSetDebugger(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSetDebugger, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSetDefaultDebugger(::$UnionPetscLib)
    @chk ccall((:PetscSetDefaultDebugger, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscSetDebuggerFromString(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscSetDebuggerFromString, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar},),
        arg1,
    )
end

@for_petsc function PetscAttachDebugger(::$UnionPetscLib)
    @chk ccall((:PetscAttachDebugger, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscStopForDebugger(::$UnionPetscLib)
    @chk ccall((:PetscStopForDebugger, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscWaitOnError(::$UnionPetscLib)
    @chk ccall((:PetscWaitOnError, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscGatherNumberOfMessages(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscGatherNumberOfMessages, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscMPIInt}, Ptr{PetscMPIInt}, Ptr{PetscMPIInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscGatherMessageLengths(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscGatherMessageLengths, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            PetscMPIInt,
            PetscMPIInt,
            Ptr{PetscMPIInt},
            Ptr{Ptr{PetscMPIInt}},
            Ptr{Ptr{PetscMPIInt}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscGatherMessageLengths2(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscGatherMessageLengths2, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            PetscMPIInt,
            PetscMPIInt,
            Ptr{PetscMPIInt},
            Ptr{PetscMPIInt},
            Ptr{Ptr{PetscMPIInt}},
            Ptr{Ptr{PetscMPIInt}},
            Ptr{Ptr{PetscMPIInt}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function PetscCommBuildTwoSided(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:PetscCommBuildTwoSided, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            PetscMPIInt,
            MPI_Datatype,
            PetscMPIInt,
            Ptr{PetscMPIInt},
            Ptr{Cvoid},
            Ptr{PetscMPIInt},
            Ptr{Ptr{PetscMPIInt}},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function PetscCommBuildTwoSidedF(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    send,
    recv,
    ctx,
)
    @chk ccall(
        (:PetscCommBuildTwoSidedF, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            PetscMPIInt,
            MPI_Datatype,
            PetscMPIInt,
            Ptr{PetscMPIInt},
            Ptr{Cvoid},
            Ptr{PetscMPIInt},
            Ptr{Ptr{PetscMPIInt}},
            Ptr{Cvoid},
            PetscMPIInt,
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        send,
        recv,
        ctx,
    )
end

@for_petsc function PetscCommBuildTwoSidedSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscCommBuildTwoSidedSetType, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, PetscBuildTwoSidedType),
        arg1,
        arg2,
    )
end

@for_petsc function PetscCommBuildTwoSidedGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscCommBuildTwoSidedGetType, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscBuildTwoSidedType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSSEIsEnabled(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSSEIsEnabled, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscBool}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscObjectComm(::$UnionPetscLib, arg1)
    ccall((:PetscObjectComm, $petsc_library), MPI_Comm, (PetscObject,), arg1)
end

@for_petsc function PetscSubcommParent(::$UnionPetscLib, scomm)
    ccall(
        (:PetscSubcommParent, $petsc_library),
        MPI_Comm,
        (PetscSubcomm,),
        scomm,
    )
end

@for_petsc function PetscSubcommChild(::$UnionPetscLib, scomm)
    ccall(
        (:PetscSubcommChild, $petsc_library),
        MPI_Comm,
        (PetscSubcomm,),
        scomm,
    )
end

@for_petsc function PetscSubcommContiguousParent(::$UnionPetscLib, scomm)
    ccall(
        (:PetscSubcommContiguousParent, $petsc_library),
        MPI_Comm,
        (PetscSubcomm,),
        scomm,
    )
end

@for_petsc function PetscSubcommCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSubcommCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscSubcomm}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSubcommDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscSubcommDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscSubcomm},),
        arg1,
    )
end

@for_petsc function PetscSubcommSetNumber(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSubcommSetNumber, $petsc_library),
        PetscErrorCode,
        (PetscSubcomm, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSubcommSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSubcommSetType, $petsc_library),
        PetscErrorCode,
        (PetscSubcomm, PetscSubcommType),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSubcommSetTypeGeneral(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSubcommSetTypeGeneral, $petsc_library),
        PetscErrorCode,
        (PetscSubcomm, PetscMPIInt, PetscMPIInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSubcommView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSubcommView, $petsc_library),
        PetscErrorCode,
        (PetscSubcomm, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSubcommSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscSubcommSetFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscSubcomm,),
        arg1,
    )
end

@for_petsc function PetscSubcommSetOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSubcommSetOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (PetscSubcomm, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSubcommGetParent(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSubcommGetParent, $petsc_library),
        PetscErrorCode,
        (PetscSubcomm, Ptr{MPI_Comm}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSubcommGetContiguousParent(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscSubcommGetContiguousParent, $petsc_library),
        PetscErrorCode,
        (PetscSubcomm, Ptr{MPI_Comm}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSubcommGetChild(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSubcommGetChild, $petsc_library),
        PetscErrorCode,
        (PetscSubcomm, Ptr{MPI_Comm}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscHeapCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscHeapCreate, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{PetscHeap}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscHeapAdd(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscHeapAdd, $petsc_library),
        PetscErrorCode,
        (PetscHeap, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscHeapPop(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscHeapPop, $petsc_library),
        PetscErrorCode,
        (PetscHeap, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscHeapPeek(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscHeapPeek, $petsc_library),
        PetscErrorCode,
        (PetscHeap, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscHeapStash(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscHeapStash, $petsc_library),
        PetscErrorCode,
        (PetscHeap, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscHeapUnstash(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscHeapUnstash, $petsc_library),
        PetscErrorCode,
        (PetscHeap,),
        arg1,
    )
end

@for_petsc function PetscHeapDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscHeapDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscHeap},),
        arg1,
    )
end

@for_petsc function PetscHeapView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscHeapView, $petsc_library),
        PetscErrorCode,
        (PetscHeap, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscProcessPlacementView(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscProcessPlacementView, $petsc_library),
        PetscErrorCode,
        (PetscViewer,),
        arg1,
    )
end

@for_petsc function PetscShmCommGet(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscShmCommGet, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscShmComm}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscShmCommGlobalToLocal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscShmCommGlobalToLocal, $petsc_library),
        PetscErrorCode,
        (PetscShmComm, PetscMPIInt, Ptr{PetscMPIInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscShmCommLocalToGlobal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscShmCommLocalToGlobal, $petsc_library),
        PetscErrorCode,
        (PetscShmComm, PetscMPIInt, Ptr{PetscMPIInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscShmCommGetMpiShmComm(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscShmCommGetMpiShmComm, $petsc_library),
        PetscErrorCode,
        (PetscShmComm, Ptr{MPI_Comm}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscOmpCtrlCreate(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscOmpCtrlCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, Ptr{PetscOmpCtrl}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscOmpCtrlGetOmpComms(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscOmpCtrlGetOmpComms, $petsc_library),
        PetscErrorCode,
        (PetscOmpCtrl, Ptr{MPI_Comm}, Ptr{MPI_Comm}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscOmpCtrlDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscOmpCtrlDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscOmpCtrl},),
        arg1,
    )
end

@for_petsc function PetscOmpCtrlBarrier(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscOmpCtrlBarrier, $petsc_library),
        PetscErrorCode,
        (PetscOmpCtrl,),
        arg1,
    )
end

@for_petsc function PetscOmpCtrlOmpRegionOnMasterBegin(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscOmpCtrlOmpRegionOnMasterBegin, $petsc_library),
        PetscErrorCode,
        (PetscOmpCtrl,),
        arg1,
    )
end

@for_petsc function PetscOmpCtrlOmpRegionOnMasterEnd(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscOmpCtrlOmpRegionOnMasterEnd, $petsc_library),
        PetscErrorCode,
        (PetscOmpCtrl,),
        arg1,
    )
end

@for_petsc function PetscSegBufferCreate(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSegBufferCreate, $petsc_library),
        PetscErrorCode,
        (Csize_t, Csize_t, Ptr{PetscSegBuffer}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSegBufferDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscSegBufferDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscSegBuffer},),
        arg1,
    )
end

@for_petsc function PetscSegBufferGet(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSegBufferGet, $petsc_library),
        PetscErrorCode,
        (PetscSegBuffer, Csize_t, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSegBufferExtractAlloc(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSegBufferExtractAlloc, $petsc_library),
        PetscErrorCode,
        (PetscSegBuffer, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSegBufferExtractTo(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSegBufferExtractTo, $petsc_library),
        PetscErrorCode,
        (PetscSegBuffer, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSegBufferExtractInPlace(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSegBufferExtractInPlace, $petsc_library),
        PetscErrorCode,
        (PetscSegBuffer, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSegBufferGetSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSegBufferGetSize, $petsc_library),
        PetscErrorCode,
        (PetscSegBuffer, Ptr{Csize_t}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSegBufferUnuse(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSegBufferUnuse, $petsc_library),
        PetscErrorCode,
        (PetscSegBuffer, Csize_t),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSegBufferGetInts(::$UnionPetscLib, seg, count, slot)
    @chk ccall(
        (:PetscSegBufferGetInts, $petsc_library),
        PetscErrorCode,
        (PetscSegBuffer, Csize_t, Ptr{Ptr{$PetscInt}}),
        seg,
        count,
        slot,
    )
end

@for_petsc function PetscOptionsHelpPrintedDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscOptionsHelpPrintedDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscOptionsHelpPrinted},),
        arg1,
    )
end

@for_petsc function PetscOptionsHelpPrintedCreate(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscOptionsHelpPrintedCreate, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscOptionsHelpPrinted},),
        arg1,
    )
end

@for_petsc function PetscOptionsHelpPrintedCheck(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscOptionsHelpPrintedCheck, $petsc_library),
        PetscErrorCode,
        (PetscOptionsHelpPrinted, Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscCitationsRegister(::$UnionPetscLib, cit, set)
    @chk ccall(
        (:PetscCitationsRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{PetscBool}),
        cit,
        set,
    )
end

@for_petsc function PetscURLShorten(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscURLShorten, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscGoogleDriveAuthorize(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscGoogleDriveAuthorize, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscGoogleDriveRefresh(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscGoogleDriveRefresh, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscGoogleDriveUpload(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscGoogleDriveUpload, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscBoxAuthorize(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscBoxAuthorize, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscBoxRefresh(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscBoxRefresh, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscGlobusGetTransfers(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscGlobusGetTransfers, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscTextBelt(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscTextBelt, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscTellMyCell(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscTellMyCell, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscPullJSONValue(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscPullJSONValue, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Csize_t, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscPushJSONValue(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscPushJSONValue, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MPIU_Win_allocate_shared(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MPIU_Win_allocate_shared, $petsc_library),
        PetscErrorCode,
        (MPI_Aint, PetscMPIInt, MPI_Info, MPI_Comm, Ptr{Cvoid}, Ptr{MPI_Win}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MPIU_Win_shared_query(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MPIU_Win_shared_query, $petsc_library),
        PetscErrorCode,
        (MPI_Win, PetscMPIInt, Ptr{MPI_Aint}, Ptr{PetscMPIInt}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscHasExternalPackage(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscHasExternalPackage, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

mutable struct _n_PetscBag end

const PetscBag = Ptr{_n_PetscBag}

mutable struct _n_PetscBagItem end

const PetscBagItem = Ptr{_n_PetscBagItem}

@for_petsc function PetscBagCreate(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscBagCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Csize_t, Ptr{PetscBag}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscBagDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscBagDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscBag},),
        arg1,
    )
end

@for_petsc function PetscBagGetData(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscBagGetData, $petsc_library),
        PetscErrorCode,
        (PetscBag, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscBagRegisterReal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscBagRegisterReal, $petsc_library),
        PetscErrorCode,
        (PetscBag, Ptr{Cvoid}, $PetscReal, Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscBagRegisterRealArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscBagRegisterRealArray, $petsc_library),
        PetscErrorCode,
        (PetscBag, Ptr{Cvoid}, $PetscInt, Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscBagRegisterString(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscBagRegisterString, $petsc_library),
        PetscErrorCode,
        (PetscBag, Ptr{Cvoid}, $PetscInt, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscBagRegisterScalar(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscBagRegisterScalar, $petsc_library),
        PetscErrorCode,
        (PetscBag, Ptr{Cvoid}, $PetscScalar, Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscBagRegisterInt(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscBagRegisterInt, $petsc_library),
        PetscErrorCode,
        (PetscBag, Ptr{Cvoid}, $PetscInt, Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscBagRegisterInt64(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscBagRegisterInt64, $petsc_library),
        PetscErrorCode,
        (PetscBag, Ptr{Cvoid}, PetscInt64, Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscBagRegisterIntArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscBagRegisterIntArray, $petsc_library),
        PetscErrorCode,
        (PetscBag, Ptr{Cvoid}, $PetscInt, Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscBagRegisterEnum(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscBagRegisterEnum, $petsc_library),
        PetscErrorCode,
        (
            PetscBag,
            Ptr{Cvoid},
            Ptr{Ptr{Cchar}},
            PetscEnum,
            Ptr{Cchar},
            Ptr{Cchar},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscBagRegisterBool(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscBagRegisterBool, $petsc_library),
        PetscErrorCode,
        (PetscBag, Ptr{Cvoid}, PetscBool, Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscBagRegisterBoolArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscBagRegisterBoolArray, $petsc_library),
        PetscErrorCode,
        (PetscBag, Ptr{Cvoid}, $PetscInt, Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscBagGetNames(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscBagGetNames, $petsc_library),
        PetscErrorCode,
        (PetscBag, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscBagSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscBagSetFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscBag,),
        arg1,
    )
end

@for_petsc function PetscBagGetName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscBagGetName, $petsc_library),
        PetscErrorCode,
        (PetscBag, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscBagSetName(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscBagSetName, $petsc_library),
        PetscErrorCode,
        (PetscBag, Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscBagSetOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscBagSetOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (PetscBag, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscBagView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscBagView, $petsc_library),
        PetscErrorCode,
        (PetscBag, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscBagLoad(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscBagLoad, $petsc_library),
        PetscErrorCode,
        (PetscViewer, PetscBag),
        arg1,
        arg2,
    )
end

@for_petsc function PetscBagViewFromOptions(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscBagViewFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscBag, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscBagSetViewer(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscBagSetViewer, $petsc_library),
        PetscErrorCode,
        (PetscBag, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscBagSetLoader(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscBagSetLoader, $petsc_library),
        PetscErrorCode,
        (PetscBag, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscBagSetDestroy(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscBagSetDestroy, $petsc_library),
        PetscErrorCode,
        (PetscBag, Ptr{Cvoid}),
        arg1,
        arg2,
    )
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

@for_petsc function PetscViewerInitializePackage(::$UnionPetscLib)
    @chk ccall(
        (:PetscViewerInitializePackage, $petsc_library),
        PetscErrorCode,
        (),
    )
end

@for_petsc function PetscViewerRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscViewer}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscViewerSetFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscViewer,),
        arg1,
    )
end

@for_petsc function PetscViewerASCIIOpenWithFILE(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscViewerASCIIOpenWithFILE, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Libc.FILE}, Ptr{PetscViewer}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscViewerASCIIOpen(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscViewerASCIIOpen, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{PetscViewer}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscViewerASCIISetFILE(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerASCIISetFILE, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Libc.FILE}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerBinaryOpen(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscViewerBinaryOpen, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, PetscFileMode, Ptr{PetscViewer}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscViewerADIOSOpen(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscViewerADIOSOpen, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, PetscFileMode, Ptr{PetscViewer}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscViewerADIOS2Open(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscViewerADIOS2Open, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, PetscFileMode, Ptr{PetscViewer}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscViewerBinaryGetFlowControl(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscViewerBinaryGetFlowControl, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerBinarySetFlowControl(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscViewerBinarySetFlowControl, $petsc_library),
        PetscErrorCode,
        (PetscViewer, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerBinarySetUseMPIIO(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerBinarySetUseMPIIO, $petsc_library),
        PetscErrorCode,
        (PetscViewer, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerBinaryGetUseMPIIO(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerBinaryGetUseMPIIO, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerBinaryGetMPIIODescriptor(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscViewerBinaryGetMPIIODescriptor, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{MPI_File}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerBinaryGetMPIIOOffset(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscViewerBinaryGetMPIIOOffset, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{MPI_Offset}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerBinaryAddMPIIOOffset(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscViewerBinaryAddMPIIOOffset, $petsc_library),
        PetscErrorCode,
        (PetscViewer, MPI_Offset),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerSocketOpen(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscViewerSocketOpen, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Cint, Ptr{PetscViewer}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscViewerStringOpen(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscViewerStringOpen, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Csize_t, Ptr{PetscViewer}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscViewerDrawOpen(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscViewerDrawOpen, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            Ptr{Cchar},
            Ptr{Cchar},
            Cint,
            Cint,
            Cint,
            Cint,
            Ptr{PetscViewer},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function PetscViewerDrawSetDrawType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerDrawSetDrawType, $petsc_library),
        PetscErrorCode,
        (PetscViewer, PetscDrawType),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerDrawGetDrawType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerDrawGetDrawType, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{PetscDrawType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerDrawSetTitle(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerDrawSetTitle, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerDrawGetTitle(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerDrawGetTitle, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerDrawGetDraw(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscViewerDrawGetDraw, $petsc_library),
        PetscErrorCode,
        (PetscViewer, $PetscInt, Ptr{PetscDraw}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscViewerDrawBaseAdd(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerDrawBaseAdd, $petsc_library),
        PetscErrorCode,
        (PetscViewer, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerDrawBaseSet(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerDrawBaseSet, $petsc_library),
        PetscErrorCode,
        (PetscViewer, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerDrawGetDrawLG(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscViewerDrawGetDrawLG, $petsc_library),
        PetscErrorCode,
        (PetscViewer, $PetscInt, Ptr{PetscDrawLG}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscViewerDrawGetDrawAxis(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscViewerDrawGetDrawAxis, $petsc_library),
        PetscErrorCode,
        (PetscViewer, $PetscInt, Ptr{PetscDrawAxis}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscViewerMathematicaOpen(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscViewerMathematicaOpen, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Cint, Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscViewer}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscViewerSiloOpen(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscViewerSiloOpen, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{PetscViewer}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscViewerMatlabOpen(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscViewerMatlabOpen, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, PetscFileMode, Ptr{PetscViewer}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@enum PetscViewerGLVisType::UInt32 begin
    PETSC_VIEWER_GLVIS_DUMP = 0
    PETSC_VIEWER_GLVIS_SOCKET = 1
end

@for_petsc function PetscViewerGLVisOpen(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscViewerGLVisOpen, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            PetscViewerGLVisType,
            Ptr{Cchar},
            $PetscInt,
            Ptr{PetscViewer},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscViewerGLVisSetPrecision(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerGLVisSetPrecision, $petsc_library),
        PetscErrorCode,
        (PetscViewer, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerGLVisSetSnapId(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerGLVisSetSnapId, $petsc_library),
        PetscErrorCode,
        (PetscViewer, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerGLVisSetFields(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscViewerGLVisSetFields, $petsc_library),
        PetscErrorCode,
        (
            PetscViewer,
            $PetscInt,
            Ptr{Ptr{Cchar}},
            Ptr{$PetscInt},
            Ptr{Cvoid},
            Ptr{PetscObject},
            Ptr{Cvoid},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function PetscViewerGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerGetType, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{PetscViewerType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerSetType, $petsc_library),
        PetscErrorCode,
        (PetscViewer, PetscViewerType),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscViewerDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscViewer},),
        arg1,
    )
end

@for_petsc function PetscViewerGetSubViewer(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscViewerGetSubViewer, $petsc_library),
        PetscErrorCode,
        (PetscViewer, MPI_Comm, Ptr{PetscViewer}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscViewerRestoreSubViewer(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscViewerRestoreSubViewer, $petsc_library),
        PetscErrorCode,
        (PetscViewer, MPI_Comm, Ptr{PetscViewer}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscViewerSetUp(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscViewerSetUp, $petsc_library),
        PetscErrorCode,
        (PetscViewer,),
        arg1,
    )
end

@for_petsc function PetscViewerView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerView, $petsc_library),
        PetscErrorCode,
        (PetscViewer, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerViewFromOptions(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscViewerViewFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscViewer, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscViewerSetOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerSetOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerAppendOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerAppendOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerGetOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerGetOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerReadable(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerReadable, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerWritable(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerWritable, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerCheckReadable(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscViewerCheckReadable, $petsc_library),
        PetscErrorCode,
        (PetscViewer,),
        arg1,
    )
end

@for_petsc function PetscViewerCheckWritable(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscViewerCheckWritable, $petsc_library),
        PetscErrorCode,
        (PetscViewer,),
        arg1,
    )
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

@for_petsc function PetscViewerSetFormat(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerSetFormat, $petsc_library),
        PetscErrorCode,
        (PetscViewer, PetscViewerFormat),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerPushFormat(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerPushFormat, $petsc_library),
        PetscErrorCode,
        (PetscViewer, PetscViewerFormat),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerPopFormat(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscViewerPopFormat, $petsc_library),
        PetscErrorCode,
        (PetscViewer,),
        arg1,
    )
end

@for_petsc function PetscViewerGetFormat(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerGetFormat, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{PetscViewerFormat}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerFlush(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscViewerFlush, $petsc_library),
        PetscErrorCode,
        (PetscViewer,),
        arg1,
    )
end

@for_petsc function PetscOptionsPushGetViewerOff(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscOptionsPushGetViewerOff, $petsc_library),
        PetscErrorCode,
        (PetscBool,),
        arg1,
    )
end

@for_petsc function PetscOptionsPopGetViewerOff(::$UnionPetscLib)
    @chk ccall(
        (:PetscOptionsPopGetViewerOff, $petsc_library),
        PetscErrorCode,
        (),
    )
end

@for_petsc function PetscOptionsGetViewerOff(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscOptionsGetViewerOff, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscBool},),
        arg1,
    )
end

@for_petsc function PetscOptionsGetViewer(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscOptionsGetViewer, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            PetscOptions,
            Ptr{Cchar},
            Ptr{Cchar},
            Ptr{PetscViewer},
            Ptr{PetscViewerFormat},
            Ptr{PetscBool},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

mutable struct PetscViewerAndFormat
    viewer::PetscViewer
    format::PetscViewerFormat
    lg::PetscDrawLG
    data::Ptr{Cvoid}
    PetscViewerAndFormat() = new()
end

@for_petsc function PetscViewerAndFormatCreate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscViewerAndFormatCreate, $petsc_library),
        PetscErrorCode,
        (PetscViewer, PetscViewerFormat, Ptr{Ptr{PetscViewerAndFormat}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscViewerAndFormatDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscViewerAndFormatDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{Ptr{PetscViewerAndFormat}},),
        arg1,
    )
end

@for_petsc function PetscViewerASCIIGetPointer(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerASCIIGetPointer, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Ptr{Libc.FILE}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerFileGetMode(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerFileGetMode, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{PetscFileMode}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerFileSetMode(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerFileSetMode, $petsc_library),
        PetscErrorCode,
        (PetscViewer, PetscFileMode),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerRead(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscViewerRead, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Cvoid}, $PetscInt, Ptr{$PetscInt}, PetscDataType),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscViewerASCIIPushSynchronized(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscViewerASCIIPushSynchronized, $petsc_library),
        PetscErrorCode,
        (PetscViewer,),
        arg1,
    )
end

@for_petsc function PetscViewerASCIIPopSynchronized(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscViewerASCIIPopSynchronized, $petsc_library),
        PetscErrorCode,
        (PetscViewer,),
        arg1,
    )
end

@for_petsc function PetscViewerASCIIPushTab(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscViewerASCIIPushTab, $petsc_library),
        PetscErrorCode,
        (PetscViewer,),
        arg1,
    )
end

@for_petsc function PetscViewerASCIIPopTab(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscViewerASCIIPopTab, $petsc_library),
        PetscErrorCode,
        (PetscViewer,),
        arg1,
    )
end

@for_petsc function PetscViewerASCIIUseTabs(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerASCIIUseTabs, $petsc_library),
        PetscErrorCode,
        (PetscViewer, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerASCIISetTab(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerASCIISetTab, $petsc_library),
        PetscErrorCode,
        (PetscViewer, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerASCIIGetTab(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerASCIIGetTab, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerASCIIAddTab(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerASCIIAddTab, $petsc_library),
        PetscErrorCode,
        (PetscViewer, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerASCIISubtractTab(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerASCIISubtractTab, $petsc_library),
        PetscErrorCode,
        (PetscViewer, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerASCIIRead(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscViewerASCIIRead, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Cvoid}, $PetscInt, Ptr{$PetscInt}, PetscDataType),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscViewerBinaryGetDescriptor(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerBinaryGetDescriptor, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Cint}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerBinaryGetInfoPointer(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscViewerBinaryGetInfoPointer, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Ptr{Libc.FILE}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerBinaryRead(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscViewerBinaryRead, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Cvoid}, $PetscInt, Ptr{$PetscInt}, PetscDataType),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscViewerBinaryWrite(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscViewerBinaryWrite, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Cvoid}, $PetscInt, PetscDataType),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscViewerBinaryReadAll(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscViewerBinaryReadAll, $petsc_library),
        PetscErrorCode,
        (
            PetscViewer,
            Ptr{Cvoid},
            $PetscInt,
            $PetscInt,
            $PetscInt,
            PetscDataType,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscViewerBinaryWriteAll(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscViewerBinaryWriteAll, $petsc_library),
        PetscErrorCode,
        (
            PetscViewer,
            Ptr{Cvoid},
            $PetscInt,
            $PetscInt,
            $PetscInt,
            PetscDataType,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscViewerStringSetString(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscViewerStringSetString, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscViewerStringGetStringRead(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscViewerStringGetStringRead, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Ptr{Cchar}}, Ptr{Csize_t}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscViewerStringSetOwnString(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscViewerStringSetOwnString, $petsc_library),
        PetscErrorCode,
        (PetscViewer,),
        arg1,
    )
end

@for_petsc function PetscViewerDrawClear(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscViewerDrawClear, $petsc_library),
        PetscErrorCode,
        (PetscViewer,),
        arg1,
    )
end

@for_petsc function PetscViewerDrawSetHold(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerDrawSetHold, $petsc_library),
        PetscErrorCode,
        (PetscViewer, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerDrawGetHold(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerDrawGetHold, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerDrawSetPause(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerDrawSetPause, $petsc_library),
        PetscErrorCode,
        (PetscViewer, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerDrawGetPause(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerDrawGetPause, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerDrawSetInfo(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscViewerDrawSetInfo, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Cchar}, Ptr{Cchar}, Cint, Cint, Cint, Cint),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscViewerDrawResize(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscViewerDrawResize, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Cint, Cint),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscViewerDrawSetBounds(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscViewerDrawSetBounds, $petsc_library),
        PetscErrorCode,
        (PetscViewer, $PetscInt, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscViewerDrawGetBounds(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscViewerDrawGetBounds, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{$PetscInt}, Ptr{Ptr{$PetscReal}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscViewerSocketSetConnection(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscViewerSocketSetConnection, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Cchar}, Cint),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscViewerBinarySkipInfo(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscViewerBinarySkipInfo, $petsc_library),
        PetscErrorCode,
        (PetscViewer,),
        arg1,
    )
end

@for_petsc function PetscViewerBinarySetSkipInfo(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerBinarySetSkipInfo, $petsc_library),
        PetscErrorCode,
        (PetscViewer, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerBinaryGetSkipInfo(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerBinaryGetSkipInfo, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerBinarySetSkipOptions(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscViewerBinarySetSkipOptions, $petsc_library),
        PetscErrorCode,
        (PetscViewer, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerBinaryGetSkipOptions(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscViewerBinaryGetSkipOptions, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerBinarySetSkipHeader(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerBinarySetSkipHeader, $petsc_library),
        PetscErrorCode,
        (PetscViewer, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerBinaryGetSkipHeader(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerBinaryGetSkipHeader, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerBinaryReadStringArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscViewerBinaryReadStringArray, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Ptr{Ptr{Cchar}}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerBinaryWriteStringArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscViewerBinaryWriteStringArray, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerFileSetName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerFileSetName, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerFileGetName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerFileGetName, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerVUGetPointer(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerVUGetPointer, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Ptr{Libc.FILE}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerVUSetVecSeen(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerVUSetVecSeen, $petsc_library),
        PetscErrorCode,
        (PetscViewer, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerVUGetVecSeen(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerVUGetVecSeen, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerVUFlushDeferred(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscViewerVUFlushDeferred, $petsc_library),
        PetscErrorCode,
        (PetscViewer,),
        arg1,
    )
end

@for_petsc function PetscViewerVUSetMode(::$UnionPetscLib, viewer, mode)
    @chk ccall(
        (:PetscViewerVUSetMode, $petsc_library),
        PetscErrorCode,
        (PetscViewer, PetscFileMode),
        viewer,
        mode,
    )
end

@for_petsc function PetscViewerMathematicaInitializePackage(::$UnionPetscLib)
    @chk ccall(
        (:PetscViewerMathematicaInitializePackage, $petsc_library),
        PetscErrorCode,
        (),
    )
end

@for_petsc function PetscViewerMathematicaFinalizePackage(::$UnionPetscLib)
    @chk ccall(
        (:PetscViewerMathematicaFinalizePackage, $petsc_library),
        PetscErrorCode,
        (),
    )
end

@for_petsc function PetscViewerMathematicaGetName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerMathematicaGetName, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerMathematicaSetName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerMathematicaSetName, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerMathematicaClearName(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscViewerMathematicaClearName, $petsc_library),
        PetscErrorCode,
        (PetscViewer,),
        arg1,
    )
end

@for_petsc function PetscViewerMathematicaSkipPackets(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscViewerMathematicaSkipPackets, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Cint),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerSiloGetName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerSiloGetName, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerSiloSetName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerSiloSetName, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerSiloClearName(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscViewerSiloClearName, $petsc_library),
        PetscErrorCode,
        (PetscViewer,),
        arg1,
    )
end

@for_petsc function PetscViewerSiloGetMeshName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerSiloGetMeshName, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerSiloSetMeshName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerSiloSetMeshName, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerSiloClearMeshName(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscViewerSiloClearMeshName, $petsc_library),
        PetscErrorCode,
        (PetscViewer,),
        arg1,
    )
end

@enum PetscViewerVTKFieldType::UInt32 begin
    PETSC_VTK_INVALID = 0
    PETSC_VTK_POINT_FIELD = 1
    PETSC_VTK_POINT_VECTOR_FIELD = 2
    PETSC_VTK_CELL_FIELD = 3
    PETSC_VTK_CELL_VECTOR_FIELD = 4
end

@for_petsc function PetscViewerVTKAddField(
    ::$UnionPetscLib,
    arg1,
    arg2,
    PetscViewerVTKWriteFunction,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscViewerVTKAddField, $petsc_library),
        PetscErrorCode,
        (
            PetscViewer,
            PetscObject,
            Ptr{Cvoid},
            $PetscInt,
            PetscViewerVTKFieldType,
            PetscBool,
            PetscObject,
        ),
        arg1,
        arg2,
        PetscViewerVTKWriteFunction,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscViewerVTKGetDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerVTKGetDM, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{PetscObject}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerVTKOpen(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscViewerVTKOpen, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, PetscFileMode, Ptr{PetscViewer}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PETSC_VIEWER_STDOUT_(::$UnionPetscLib, arg1)
    ccall(
        (:PETSC_VIEWER_STDOUT_, $petsc_library),
        PetscViewer,
        (MPI_Comm,),
        arg1,
    )
end

@for_petsc function PetscViewerASCIIGetStdout(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerASCIIGetStdout, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscViewer}),
        arg1,
        arg2,
    )
end

@for_petsc function PETSC_VIEWER_STDERR_(::$UnionPetscLib, arg1)
    ccall(
        (:PETSC_VIEWER_STDERR_, $petsc_library),
        PetscViewer,
        (MPI_Comm,),
        arg1,
    )
end

@for_petsc function PetscViewerASCIIGetStderr(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerASCIIGetStderr, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscViewer}),
        arg1,
        arg2,
    )
end

@for_petsc function PETSC_VIEWER_DRAW_(::$UnionPetscLib, arg1)
    ccall((:PETSC_VIEWER_DRAW_, $petsc_library), PetscViewer, (MPI_Comm,), arg1)
end

@for_petsc function PETSC_VIEWER_SOCKET_(::$UnionPetscLib, arg1)
    ccall(
        (:PETSC_VIEWER_SOCKET_, $petsc_library),
        PetscViewer,
        (MPI_Comm,),
        arg1,
    )
end

@for_petsc function PETSC_VIEWER_BINARY_(::$UnionPetscLib, arg1)
    ccall(
        (:PETSC_VIEWER_BINARY_, $petsc_library),
        PetscViewer,
        (MPI_Comm,),
        arg1,
    )
end

@for_petsc function PETSC_VIEWER_MATLAB_(::$UnionPetscLib, arg1)
    ccall(
        (:PETSC_VIEWER_MATLAB_, $petsc_library),
        PetscViewer,
        (MPI_Comm,),
        arg1,
    )
end

@for_petsc function PETSC_VIEWER_HDF5_(::$UnionPetscLib, arg1)
    ccall((:PETSC_VIEWER_HDF5_, $petsc_library), PetscViewer, (MPI_Comm,), arg1)
end

@for_petsc function PETSC_VIEWER_GLVIS_(::$UnionPetscLib, arg1)
    ccall(
        (:PETSC_VIEWER_GLVIS_, $petsc_library),
        PetscViewer,
        (MPI_Comm,),
        arg1,
    )
end

@for_petsc function PETSC_VIEWER_EXODUSII_(::$UnionPetscLib, arg1)
    ccall(
        (:PETSC_VIEWER_EXODUSII_, $petsc_library),
        PetscViewer,
        (MPI_Comm,),
        arg1,
    )
end

@for_petsc function PetscViewerFlowControlStart(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscViewerFlowControlStart, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscViewerFlowControlStepMain(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscViewerFlowControlStepMain, $petsc_library),
        PetscErrorCode,
        (PetscViewer, $PetscInt, Ptr{$PetscInt}, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscViewerFlowControlEndMain(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerFlowControlEndMain, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerFlowControlStepWorker(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscViewerFlowControlStepWorker, $petsc_library),
        PetscErrorCode,
        (PetscViewer, PetscMPIInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscViewerFlowControlEndWorker(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscViewerFlowControlEndWorker, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerFlowControlStepMaster(
    ::$UnionPetscLib,
    viewer,
    i,
    mcnt,
    cnt,
)
    @chk ccall(
        (:PetscViewerFlowControlStepMaster, $petsc_library),
        PetscErrorCode,
        (PetscViewer, $PetscInt, Ptr{$PetscInt}, $PetscInt),
        viewer,
        i,
        mcnt,
        cnt,
    )
end

@for_petsc function PetscViewerFlowControlEndMaster(
    ::$UnionPetscLib,
    viewer,
    mcnt,
)
    @chk ccall(
        (:PetscViewerFlowControlEndMaster, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{$PetscInt}),
        viewer,
        mcnt,
    )
end

@for_petsc function PetscViewerMatlabPutArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscViewerMatlabPutArray, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Cint, Cint, Ptr{$PetscScalar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscViewerMatlabGetArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscViewerMatlabGetArray, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Cint, Cint, Ptr{$PetscScalar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscViewerMatlabPutVariable(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscViewerMatlabPutVariable, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

mutable struct _n_PetscViewers end

const PetscViewers = Ptr{_n_PetscViewers}

@for_petsc function PetscViewersCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewersCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscViewers}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewersDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscViewersDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscViewers},),
        arg1,
    )
end

@for_petsc function PetscViewersGetViewer(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscViewersGetViewer, $petsc_library),
        PetscErrorCode,
        (PetscViewers, $PetscInt, Ptr{PetscViewer}),
        arg1,
        arg2,
        arg3,
    )
end

const PetscBT = Ptr{Cchar}

@for_petsc function PetscBTLength(::$UnionPetscLib, m)
    ccall((:PetscBTLength, $petsc_library), $PetscInt, ($PetscInt,), m)
end

@for_petsc function PetscBTMemzero(::$UnionPetscLib, m, array)
    @chk ccall(
        (:PetscBTMemzero, $petsc_library),
        PetscErrorCode,
        ($PetscInt, PetscBT),
        m,
        array,
    )
end

@for_petsc function PetscBTDestroy(::$UnionPetscLib, array)
    @chk ccall(
        (:PetscBTDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscBT},),
        array,
    )
end

@for_petsc function PetscBTLookup(::$UnionPetscLib, array, index)
    ccall(
        (:PetscBTLookup, $petsc_library),
        Cchar,
        (PetscBT, $PetscInt),
        array,
        index,
    )
end

@for_petsc function PetscBTView(::$UnionPetscLib, m, bt, viewer)
    @chk ccall(
        (:PetscBTView, $petsc_library),
        PetscErrorCode,
        ($PetscInt, PetscBT, PetscViewer),
        m,
        bt,
        viewer,
    )
end

@for_petsc function PetscBTCreate(::$UnionPetscLib, m, array)
    @chk ccall(
        (:PetscBTCreate, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{PetscBT}),
        m,
        array,
    )
end

@for_petsc function PetscBTLookupSet(::$UnionPetscLib, array, index)
    ccall(
        (:PetscBTLookupSet, $petsc_library),
        Cchar,
        (PetscBT, $PetscInt),
        array,
        index,
    )
end

@for_petsc function PetscBTSet(::$UnionPetscLib, array, index)
    @chk ccall(
        (:PetscBTSet, $petsc_library),
        PetscErrorCode,
        (PetscBT, $PetscInt),
        array,
        index,
    )
end

@for_petsc function PetscBTNegate(::$UnionPetscLib, array, index)
    @chk ccall(
        (:PetscBTNegate, $petsc_library),
        PetscErrorCode,
        (PetscBT, $PetscInt),
        array,
        index,
    )
end

@for_petsc function PetscBTLookupClear(::$UnionPetscLib, array, index)
    ccall(
        (:PetscBTLookupClear, $petsc_library),
        Cchar,
        (PetscBT, $PetscInt),
        array,
        index,
    )
end

@for_petsc function PetscBTClear(::$UnionPetscLib, array, index)
    @chk ccall(
        (:PetscBTClear, $petsc_library),
        PetscErrorCode,
        (PetscBT, $PetscInt),
        array,
        index,
    )
end

mutable struct _p_PetscMatlabEngine end

const PetscMatlabEngine = Ptr{_p_PetscMatlabEngine}

@for_petsc function PetscMatlabEngineCreate(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscMatlabEngineCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{PetscMatlabEngine}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscMatlabEngineDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscMatlabEngineDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscMatlabEngine},),
        arg1,
    )
end

@for_petsc function PetscMatlabEngineGetOutput(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscMatlabEngineGetOutput, $petsc_library),
        PetscErrorCode,
        (PetscMatlabEngine, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscMatlabEnginePrintOutput(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscMatlabEnginePrintOutput, $petsc_library),
        PetscErrorCode,
        (PetscMatlabEngine, Ptr{Libc.FILE}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscMatlabEnginePut(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscMatlabEnginePut, $petsc_library),
        PetscErrorCode,
        (PetscMatlabEngine, PetscObject),
        arg1,
        arg2,
    )
end

@for_petsc function PetscMatlabEngineGet(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscMatlabEngineGet, $petsc_library),
        PetscErrorCode,
        (PetscMatlabEngine, PetscObject),
        arg1,
        arg2,
    )
end

@for_petsc function PetscMatlabEnginePutArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscMatlabEnginePutArray, $petsc_library),
        PetscErrorCode,
        (PetscMatlabEngine, Cint, Cint, Ptr{$PetscScalar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscMatlabEngineGetArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscMatlabEngineGetArray, $petsc_library),
        PetscErrorCode,
        (PetscMatlabEngine, Cint, Cint, Ptr{$PetscScalar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PETSC_MATLAB_ENGINE_(::$UnionPetscLib, arg1)
    ccall(
        (:PETSC_MATLAB_ENGINE_, $petsc_library),
        PetscMatlabEngine,
        (MPI_Comm,),
        arg1,
    )
end

@for_petsc function PetscDrawInitializePackage(::$UnionPetscLib)
    @chk ccall(
        (:PetscDrawInitializePackage, $petsc_library),
        PetscErrorCode,
        (),
    )
end

@for_petsc function PetscDrawFinalizePackage(::$UnionPetscLib)
    @chk ccall((:PetscDrawFinalizePackage, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscDrawRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawGetType, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Ptr{PetscDrawType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawSetType, $petsc_library),
        PetscErrorCode,
        (PetscDraw, PetscDrawType),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawCreate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscDrawCreate, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            Ptr{Cchar},
            Ptr{Cchar},
            Cint,
            Cint,
            Cint,
            Cint,
            Ptr{PetscDraw},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function PetscDrawSetOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawSetOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawSetFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscDraw,),
        arg1,
    )
end

@for_petsc function PetscDrawSetSave(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawSetSave, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawSetSaveMovie(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawSetSaveMovie, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawSetSaveFinalImage(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawSetSaveFinalImage, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawView, $petsc_library),
        PetscErrorCode,
        (PetscDraw, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawViewFromOptions(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDrawViewFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscDraw, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDrawRealToColor(::$UnionPetscLib, value, min, max)
    ccall(
        (:PetscDrawRealToColor, $petsc_library),
        Cint,
        ($PetscReal, $PetscReal, $PetscReal),
        value,
        min,
        max,
    )
end

@for_petsc function PetscDrawOpenX(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscDrawOpenX, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            Ptr{Cchar},
            Ptr{Cchar},
            Cint,
            Cint,
            Cint,
            Cint,
            Ptr{PetscDraw},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function PetscDrawOpenImage(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDrawOpenImage, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Cint, Cint, Ptr{PetscDraw}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDrawOpenNull(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawOpenNull, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscDraw}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscDraw},),
        arg1,
    )
end

@for_petsc function PetscDrawIsNull(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawIsNull, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawGetPopup(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawGetPopup, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Ptr{PetscDraw}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawScalePopup(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDrawScalePopup, $petsc_library),
        PetscErrorCode,
        (PetscDraw, $PetscReal, $PetscReal),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDrawCheckResizedWindow(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawCheckResizedWindow, $petsc_library),
        PetscErrorCode,
        (PetscDraw,),
        arg1,
    )
end

@for_petsc function PetscDrawResizeWindow(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDrawResizeWindow, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Cint, Cint),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDrawGetWindowSize(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDrawGetWindowSize, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Ptr{Cint}, Ptr{Cint}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDrawPixelToCoordinate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDrawPixelToCoordinate, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Cint, Cint, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDrawCoordinateToPixel(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDrawCoordinateToPixel, $petsc_library),
        PetscErrorCode,
        (PetscDraw, $PetscReal, $PetscReal, Ptr{Cint}, Ptr{Cint}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDrawIndicatorFunction(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscDrawIndicatorFunction, $petsc_library),
        PetscErrorCode,
        (
            PetscDraw,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            Cint,
            Ptr{Cvoid},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function PetscDrawLine(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscDrawLine, $petsc_library),
        PetscErrorCode,
        (PetscDraw, $PetscReal, $PetscReal, $PetscReal, $PetscReal, Cint),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscDrawArrow(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscDrawArrow, $petsc_library),
        PetscErrorCode,
        (PetscDraw, $PetscReal, $PetscReal, $PetscReal, $PetscReal, Cint),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscDrawLineSetWidth(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawLineSetWidth, $petsc_library),
        PetscErrorCode,
        (PetscDraw, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawLineGetWidth(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawLineGetWidth, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@enum PetscDrawMarkerType::UInt32 begin
    PETSC_DRAW_MARKER_CROSS = 0
    PETSC_DRAW_MARKER_POINT = 1
    PETSC_DRAW_MARKER_PLUS = 2
    PETSC_DRAW_MARKER_CIRCLE = 3
end

@for_petsc function PetscDrawMarker(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscDrawMarker, $petsc_library),
        PetscErrorCode,
        (PetscDraw, $PetscReal, $PetscReal, Cint),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDrawSetMarkerType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawSetMarkerType, $petsc_library),
        PetscErrorCode,
        (PetscDraw, PetscDrawMarkerType),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawGetMarkerType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawGetMarkerType, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Ptr{PetscDrawMarkerType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawPoint(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscDrawPoint, $petsc_library),
        PetscErrorCode,
        (PetscDraw, $PetscReal, $PetscReal, Cint),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDrawPointPixel(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDrawPointPixel, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Cint, Cint, Cint),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDrawPointSetSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawPointSetSize, $petsc_library),
        PetscErrorCode,
        (PetscDraw, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawRectangle(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:PetscDrawRectangle, $petsc_library),
        PetscErrorCode,
        (
            PetscDraw,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            Cint,
            Cint,
            Cint,
            Cint,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function PetscDrawTriangle(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
)
    @chk ccall(
        (:PetscDrawTriangle, $petsc_library),
        PetscErrorCode,
        (
            PetscDraw,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            Cint,
            Cint,
            Cint,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
    )
end

@for_petsc function PetscDrawEllipse(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscDrawEllipse, $petsc_library),
        PetscErrorCode,
        (PetscDraw, $PetscReal, $PetscReal, $PetscReal, $PetscReal, Cint),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscDrawTensorContourPatch(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscDrawTensorContourPatch, $petsc_library),
        PetscErrorCode,
        (
            PetscDraw,
            Cint,
            Cint,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            $PetscReal,
            $PetscReal,
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function PetscDrawTensorContour(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscDrawTensorContour, $petsc_library),
        PetscErrorCode,
        (
            PetscDraw,
            Cint,
            Cint,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscDrawString(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDrawString, $petsc_library),
        PetscErrorCode,
        (PetscDraw, $PetscReal, $PetscReal, Cint, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDrawStringCentered(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDrawStringCentered, $petsc_library),
        PetscErrorCode,
        (PetscDraw, $PetscReal, $PetscReal, Cint, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDrawStringBoxed(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscDrawStringBoxed, $petsc_library),
        PetscErrorCode,
        (
            PetscDraw,
            $PetscReal,
            $PetscReal,
            Cint,
            Cint,
            Ptr{Cchar},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function PetscDrawStringVertical(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDrawStringVertical, $petsc_library),
        PetscErrorCode,
        (PetscDraw, $PetscReal, $PetscReal, Cint, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDrawStringSetSize(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDrawStringSetSize, $petsc_library),
        PetscErrorCode,
        (PetscDraw, $PetscReal, $PetscReal),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDrawStringGetSize(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDrawStringGetSize, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDrawSetViewPort(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDrawSetViewPort, $petsc_library),
        PetscErrorCode,
        (PetscDraw, $PetscReal, $PetscReal, $PetscReal, $PetscReal),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDrawGetViewPort(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDrawGetViewPort, $petsc_library),
        PetscErrorCode,
        (
            PetscDraw,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDrawSplitViewPort(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawSplitViewPort, $petsc_library),
        PetscErrorCode,
        (PetscDraw,),
        arg1,
    )
end

@for_petsc function PetscDrawSetCoordinates(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDrawSetCoordinates, $petsc_library),
        PetscErrorCode,
        (PetscDraw, $PetscReal, $PetscReal, $PetscReal, $PetscReal),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDrawGetCoordinates(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDrawGetCoordinates, $petsc_library),
        PetscErrorCode,
        (
            PetscDraw,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDrawSetTitle(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawSetTitle, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawAppendTitle(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawAppendTitle, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawGetTitle(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawGetTitle, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawSetPause(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawSetPause, $petsc_library),
        PetscErrorCode,
        (PetscDraw, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawGetPause(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawGetPause, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawPause(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawPause, $petsc_library),
        PetscErrorCode,
        (PetscDraw,),
        arg1,
    )
end

@for_petsc function PetscDrawSetDoubleBuffer(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawSetDoubleBuffer, $petsc_library),
        PetscErrorCode,
        (PetscDraw,),
        arg1,
    )
end

@for_petsc function PetscDrawClear(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawClear, $petsc_library),
        PetscErrorCode,
        (PetscDraw,),
        arg1,
    )
end

@for_petsc function PetscDrawFlush(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawFlush, $petsc_library),
        PetscErrorCode,
        (PetscDraw,),
        arg1,
    )
end

@for_petsc function PetscDrawSave(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawSave, $petsc_library),
        PetscErrorCode,
        (PetscDraw,),
        arg1,
    )
end

@for_petsc function PetscDrawSaveMovie(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawSaveMovie, $petsc_library),
        PetscErrorCode,
        (PetscDraw,),
        arg1,
    )
end

@for_petsc function PetscDrawBOP(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawBOP, $petsc_library),
        PetscErrorCode,
        (PetscDraw,),
        arg1,
    )
end

@for_petsc function PetscDrawEOP(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawEOP, $petsc_library),
        PetscErrorCode,
        (PetscDraw,),
        arg1,
    )
end

@for_petsc function PetscDrawSetDisplay(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawSetDisplay, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawGetSingleton(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawGetSingleton, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Ptr{PetscDraw}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawRestoreSingleton(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawRestoreSingleton, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Ptr{PetscDraw}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawGetCurrentPoint(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDrawGetCurrentPoint, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDrawSetCurrentPoint(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDrawSetCurrentPoint, $petsc_library),
        PetscErrorCode,
        (PetscDraw, $PetscReal, $PetscReal),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDrawPushCurrentPoint(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscDrawPushCurrentPoint, $petsc_library),
        PetscErrorCode,
        (PetscDraw, $PetscReal, $PetscReal),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDrawPopCurrentPoint(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawPopCurrentPoint, $petsc_library),
        PetscErrorCode,
        (PetscDraw,),
        arg1,
    )
end

@for_petsc function PetscDrawGetBoundingBox(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDrawGetBoundingBox, $petsc_library),
        PetscErrorCode,
        (
            PetscDraw,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
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

@for_petsc function PetscDrawGetMouseButton(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscDrawGetMouseButton, $petsc_library),
        PetscErrorCode,
        (
            PetscDraw,
            Ptr{PetscDrawButton},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscDrawZoom(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDrawZoom, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDrawAxisCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawAxisCreate, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Ptr{PetscDrawAxis}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawAxisDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawAxisDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscDrawAxis},),
        arg1,
    )
end

@for_petsc function PetscDrawAxisDraw(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawAxisDraw, $petsc_library),
        PetscErrorCode,
        (PetscDrawAxis,),
        arg1,
    )
end

@for_petsc function PetscDrawAxisSetLimits(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDrawAxisSetLimits, $petsc_library),
        PetscErrorCode,
        (PetscDrawAxis, $PetscReal, $PetscReal, $PetscReal, $PetscReal),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDrawAxisGetLimits(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDrawAxisGetLimits, $petsc_library),
        PetscErrorCode,
        (
            PetscDrawAxis,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDrawAxisSetHoldLimits(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawAxisSetHoldLimits, $petsc_library),
        PetscErrorCode,
        (PetscDrawAxis, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawAxisSetColors(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDrawAxisSetColors, $petsc_library),
        PetscErrorCode,
        (PetscDrawAxis, Cint, Cint, Cint),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDrawAxisSetLabels(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDrawAxisSetLabels, $petsc_library),
        PetscErrorCode,
        (PetscDrawAxis, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDrawLGCreate(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDrawLGCreate, $petsc_library),
        PetscErrorCode,
        (PetscDraw, $PetscInt, Ptr{PetscDrawLG}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDrawLGDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawLGDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscDrawLG},),
        arg1,
    )
end

@for_petsc function PetscDrawLGAddPoint(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDrawLGAddPoint, $petsc_library),
        PetscErrorCode,
        (PetscDrawLG, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDrawLGAddCommonPoint(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscDrawLGAddCommonPoint, $petsc_library),
        PetscErrorCode,
        (PetscDrawLG, $PetscReal, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDrawLGAddPoints(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDrawLGAddPoints, $petsc_library),
        PetscErrorCode,
        (PetscDrawLG, $PetscInt, Ptr{Ptr{$PetscReal}}, Ptr{Ptr{$PetscReal}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDrawLGDraw(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawLGDraw, $petsc_library),
        PetscErrorCode,
        (PetscDrawLG,),
        arg1,
    )
end

@for_petsc function PetscDrawLGSave(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawLGSave, $petsc_library),
        PetscErrorCode,
        (PetscDrawLG,),
        arg1,
    )
end

@for_petsc function PetscDrawLGView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawLGView, $petsc_library),
        PetscErrorCode,
        (PetscDrawLG, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawLGReset(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawLGReset, $petsc_library),
        PetscErrorCode,
        (PetscDrawLG,),
        arg1,
    )
end

@for_petsc function PetscDrawLGSetDimension(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawLGSetDimension, $petsc_library),
        PetscErrorCode,
        (PetscDrawLG, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawLGGetDimension(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawLGGetDimension, $petsc_library),
        PetscErrorCode,
        (PetscDrawLG, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawLGSetLegend(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawLGSetLegend, $petsc_library),
        PetscErrorCode,
        (PetscDrawLG, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawLGGetAxis(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawLGGetAxis, $petsc_library),
        PetscErrorCode,
        (PetscDrawLG, Ptr{PetscDrawAxis}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawLGGetDraw(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawLGGetDraw, $petsc_library),
        PetscErrorCode,
        (PetscDrawLG, Ptr{PetscDraw}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawLGSetUseMarkers(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawLGSetUseMarkers, $petsc_library),
        PetscErrorCode,
        (PetscDrawLG, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawLGSetLimits(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDrawLGSetLimits, $petsc_library),
        PetscErrorCode,
        (PetscDrawLG, $PetscReal, $PetscReal, $PetscReal, $PetscReal),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDrawLGSetColors(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawLGSetColors, $petsc_library),
        PetscErrorCode,
        (PetscDrawLG, Ptr{Cint}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawLGSetOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawLGSetOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (PetscDrawLG, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawLGSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawLGSetFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscDrawLG,),
        arg1,
    )
end

@for_petsc function PetscDrawSPCreate(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDrawSPCreate, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Cint, Ptr{PetscDrawSP}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDrawSPDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawSPDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscDrawSP},),
        arg1,
    )
end

@for_petsc function PetscDrawSPAddPoint(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDrawSPAddPoint, $petsc_library),
        PetscErrorCode,
        (PetscDrawSP, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDrawSPAddPoints(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDrawSPAddPoints, $petsc_library),
        PetscErrorCode,
        (PetscDrawSP, Cint, Ptr{Ptr{$PetscReal}}, Ptr{Ptr{$PetscReal}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDrawSPDraw(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawSPDraw, $petsc_library),
        PetscErrorCode,
        (PetscDrawSP, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawSPSave(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawSPSave, $petsc_library),
        PetscErrorCode,
        (PetscDrawSP,),
        arg1,
    )
end

@for_petsc function PetscDrawSPReset(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawSPReset, $petsc_library),
        PetscErrorCode,
        (PetscDrawSP,),
        arg1,
    )
end

@for_petsc function PetscDrawSPSetDimension(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawSPSetDimension, $petsc_library),
        PetscErrorCode,
        (PetscDrawSP, Cint),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawSPGetAxis(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawSPGetAxis, $petsc_library),
        PetscErrorCode,
        (PetscDrawSP, Ptr{PetscDrawAxis}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawSPGetDraw(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawSPGetDraw, $petsc_library),
        PetscErrorCode,
        (PetscDrawSP, Ptr{PetscDraw}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawSPSetLimits(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDrawSPSetLimits, $petsc_library),
        PetscErrorCode,
        (PetscDrawSP, $PetscReal, $PetscReal, $PetscReal, $PetscReal),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDrawLGSPDraw(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawLGSPDraw, $petsc_library),
        PetscErrorCode,
        (PetscDrawLG, PetscDrawSP),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawHGCreate(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDrawHGCreate, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Cint, Ptr{PetscDrawHG}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDrawHGDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawHGDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscDrawHG},),
        arg1,
    )
end

@for_petsc function PetscDrawHGAddValue(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawHGAddValue, $petsc_library),
        PetscErrorCode,
        (PetscDrawHG, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawHGDraw(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawHGDraw, $petsc_library),
        PetscErrorCode,
        (PetscDrawHG,),
        arg1,
    )
end

@for_petsc function PetscDrawHGSave(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawHGSave, $petsc_library),
        PetscErrorCode,
        (PetscDrawHG,),
        arg1,
    )
end

@for_petsc function PetscDrawHGView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawHGView, $petsc_library),
        PetscErrorCode,
        (PetscDrawHG, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawHGReset(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawHGReset, $petsc_library),
        PetscErrorCode,
        (PetscDrawHG,),
        arg1,
    )
end

@for_petsc function PetscDrawHGGetAxis(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawHGGetAxis, $petsc_library),
        PetscErrorCode,
        (PetscDrawHG, Ptr{PetscDrawAxis}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawHGGetDraw(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawHGGetDraw, $petsc_library),
        PetscErrorCode,
        (PetscDrawHG, Ptr{PetscDraw}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawHGSetLimits(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDrawHGSetLimits, $petsc_library),
        PetscErrorCode,
        (PetscDrawHG, $PetscReal, $PetscReal, Cint, Cint),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDrawHGSetNumberBins(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawHGSetNumberBins, $petsc_library),
        PetscErrorCode,
        (PetscDrawHG, Cint),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawHGSetColor(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawHGSetColor, $petsc_library),
        PetscErrorCode,
        (PetscDrawHG, Cint),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawHGCalcStats(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawHGCalcStats, $petsc_library),
        PetscErrorCode,
        (PetscDrawHG, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawHGIntegerBins(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawHGIntegerBins, $petsc_library),
        PetscErrorCode,
        (PetscDrawHG, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawBarCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawBarCreate, $petsc_library),
        PetscErrorCode,
        (PetscDraw, Ptr{PetscDrawBar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawBarSetData(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDrawBarSetData, $petsc_library),
        PetscErrorCode,
        (PetscDrawBar, $PetscInt, Ptr{$PetscReal}, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDrawBarDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawBarDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscDrawBar},),
        arg1,
    )
end

@for_petsc function PetscDrawBarDraw(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawBarDraw, $petsc_library),
        PetscErrorCode,
        (PetscDrawBar,),
        arg1,
    )
end

@for_petsc function PetscDrawBarSave(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawBarSave, $petsc_library),
        PetscErrorCode,
        (PetscDrawBar,),
        arg1,
    )
end

@for_petsc function PetscDrawBarSetColor(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawBarSetColor, $petsc_library),
        PetscErrorCode,
        (PetscDrawBar, Cint),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawBarSetLimits(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDrawBarSetLimits, $petsc_library),
        PetscErrorCode,
        (PetscDrawBar, $PetscReal, $PetscReal),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDrawBarSort(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDrawBarSort, $petsc_library),
        PetscErrorCode,
        (PetscDrawBar, PetscBool, $PetscReal),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDrawBarSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawBarSetFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscDrawBar,),
        arg1,
    )
end

@for_petsc function PetscDrawBarGetAxis(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawBarGetAxis, $petsc_library),
        PetscErrorCode,
        (PetscDrawBar, Ptr{PetscDrawAxis}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawBarGetDraw(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawBarGetDraw, $petsc_library),
        PetscErrorCode,
        (PetscDrawBar, Ptr{PetscDraw}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDrawUtilitySetCmap(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDrawUtilitySetCmap, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Cint, Ptr{Cuchar}, Ptr{Cuchar}, Ptr{Cuchar}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDrawUtilitySetGamma(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDrawUtilitySetGamma, $petsc_library),
        PetscErrorCode,
        ($PetscReal,),
        arg1,
    )
end

mutable struct _p_PetscSF end

const PetscSF = Ptr{_p_PetscSF}

const PetscSFType = Ptr{Cchar}

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

@for_petsc function ISInitializePackage(::$UnionPetscLib)
    @chk ccall((:ISInitializePackage, $petsc_library), PetscErrorCode, ())
end

const ISType = Ptr{Cchar}

@for_petsc function ISSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISSetType, $petsc_library),
        PetscErrorCode,
        (IS, ISType),
        arg1,
        arg2,
    )
end

@for_petsc function ISGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISGetType, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{ISType}),
        arg1,
        arg2,
    )
end

@for_petsc function ISRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function ISCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{IS}),
        arg1,
        arg2,
    )
end

@for_petsc function ISCreateGeneral(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:ISCreateGeneral, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, Ptr{$PetscInt}, PetscCopyMode, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function ISGeneralSetIndices(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:ISGeneralSetIndices, $petsc_library),
        PetscErrorCode,
        (IS, $PetscInt, Ptr{$PetscInt}, PetscCopyMode),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function ISCreateBlock(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:ISCreateBlock, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            PetscCopyMode,
            Ptr{IS},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function ISBlockSetIndices(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:ISBlockSetIndices, $petsc_library),
        PetscErrorCode,
        (IS, $PetscInt, $PetscInt, Ptr{$PetscInt}, PetscCopyMode),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function ISCreateStride(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:ISCreateStride, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function ISStrideSetStride(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:ISStrideSetStride, $petsc_library),
        PetscErrorCode,
        (IS, $PetscInt, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function ISDestroy(::$UnionPetscLib, arg1)
    @chk ccall((:ISDestroy, $petsc_library), PetscErrorCode, (Ptr{IS},), arg1)
end

@for_petsc function ISSetPermutation(::$UnionPetscLib, arg1)
    @chk ccall((:ISSetPermutation, $petsc_library), PetscErrorCode, (IS,), arg1)
end

@for_petsc function ISPermutation(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISPermutation, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function ISSetIdentity(::$UnionPetscLib, arg1)
    @chk ccall((:ISSetIdentity, $petsc_library), PetscErrorCode, (IS,), arg1)
end

@for_petsc function ISIdentity(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISIdentity, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function ISContiguousLocal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:ISContiguousLocal, $petsc_library),
        PetscErrorCode,
        (IS, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
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

@for_petsc function ISSetInfo(::$UnionPetscLib, arg1, arg2, arg3, arg4, arg5)
    @chk ccall(
        (:ISSetInfo, $petsc_library),
        PetscErrorCode,
        (IS, ISInfo, ISInfoType, PetscBool, PetscBool),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function ISGetInfo(::$UnionPetscLib, arg1, arg2, arg3, arg4, arg5)
    @chk ccall(
        (:ISGetInfo, $petsc_library),
        PetscErrorCode,
        (IS, ISInfo, ISInfoType, PetscBool, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function ISClearInfoCache(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISClearInfoCache, $petsc_library),
        PetscErrorCode,
        (IS, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function ISGetIndices(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISGetIndices, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function ISRestoreIndices(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISRestoreIndices, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function ISGetTotalIndices(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISGetTotalIndices, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function ISRestoreTotalIndices(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISRestoreTotalIndices, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function ISGetNonlocalIndices(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISGetNonlocalIndices, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function ISRestoreNonlocalIndices(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISRestoreNonlocalIndices, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function ISGetNonlocalIS(::$UnionPetscLib, arg1, is)
    @chk ccall(
        (:ISGetNonlocalIS, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{IS}),
        arg1,
        is,
    )
end

@for_petsc function ISRestoreNonlocalIS(::$UnionPetscLib, arg1, is)
    @chk ccall(
        (:ISRestoreNonlocalIS, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{IS}),
        arg1,
        is,
    )
end

@for_petsc function ISGetSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISGetSize, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function ISGetLocalSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISGetLocalSize, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function ISInvertPermutation(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:ISInvertPermutation, $petsc_library),
        PetscErrorCode,
        (IS, $PetscInt, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function ISView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISView, $petsc_library),
        PetscErrorCode,
        (IS, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function ISViewFromOptions(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:ISViewFromOptions, $petsc_library),
        PetscErrorCode,
        (IS, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function ISLoad(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISLoad, $petsc_library),
        PetscErrorCode,
        (IS, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function ISEqual(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:ISEqual, $petsc_library),
        PetscErrorCode,
        (IS, IS, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function ISEqualUnsorted(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:ISEqualUnsorted, $petsc_library),
        PetscErrorCode,
        (IS, IS, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function ISSort(::$UnionPetscLib, arg1)
    @chk ccall((:ISSort, $petsc_library), PetscErrorCode, (IS,), arg1)
end

@for_petsc function ISSortRemoveDups(::$UnionPetscLib, arg1)
    @chk ccall((:ISSortRemoveDups, $petsc_library), PetscErrorCode, (IS,), arg1)
end

@for_petsc function ISSorted(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISSorted, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function ISDifference(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:ISDifference, $petsc_library),
        PetscErrorCode,
        (IS, IS, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function ISSum(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:ISSum, $petsc_library),
        PetscErrorCode,
        (IS, IS, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function ISExpand(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:ISExpand, $petsc_library),
        PetscErrorCode,
        (IS, IS, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function ISIntersect(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:ISIntersect, $petsc_library),
        PetscErrorCode,
        (IS, IS, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function ISGetMinMax(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:ISGetMinMax, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function ISLocate(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:ISLocate, $petsc_library),
        PetscErrorCode,
        (IS, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function ISGetPointRange(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:ISGetPointRange, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function ISRestorePointRange(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:ISRestorePointRange, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function ISGetPointSubrange(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:ISGetPointSubrange, $petsc_library),
        PetscErrorCode,
        (IS, $PetscInt, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function ISBlockGetIndices(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISBlockGetIndices, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function ISBlockRestoreIndices(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISBlockRestoreIndices, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function ISBlockGetLocalSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISBlockGetLocalSize, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function ISBlockGetSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISBlockGetSize, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function ISGetBlockSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISGetBlockSize, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function ISSetBlockSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISSetBlockSize, $petsc_library),
        PetscErrorCode,
        (IS, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function ISStrideGetInfo(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:ISStrideGetInfo, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function ISToGeneral(::$UnionPetscLib, arg1)
    @chk ccall((:ISToGeneral, $petsc_library), PetscErrorCode, (IS,), arg1)
end

@for_petsc function ISDuplicate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISDuplicate, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{IS}),
        arg1,
        arg2,
    )
end

@for_petsc function ISCopy(::$UnionPetscLib, arg1, arg2)
    @chk ccall((:ISCopy, $petsc_library), PetscErrorCode, (IS, IS), arg1, arg2)
end

@for_petsc function ISAllGather(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISAllGather, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{IS}),
        arg1,
        arg2,
    )
end

@for_petsc function ISComplement(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:ISComplement, $petsc_library),
        PetscErrorCode,
        (IS, $PetscInt, $PetscInt, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function ISConcatenate(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:ISConcatenate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, Ptr{IS}, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function ISListToPair(::$UnionPetscLib, arg1, arg2, arg3, arg4, arg5)
    @chk ccall(
        (:ISListToPair, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, Ptr{IS}, Ptr{IS}, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function ISPairToList(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:ISPairToList, $petsc_library),
        PetscErrorCode,
        (IS, IS, Ptr{$PetscInt}, Ptr{Ptr{IS}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function ISEmbed(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:ISEmbed, $petsc_library),
        PetscErrorCode,
        (IS, IS, PetscBool, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function ISSortPermutation(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:ISSortPermutation, $petsc_library),
        PetscErrorCode,
        (IS, PetscBool, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function ISOnComm(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:ISOnComm, $petsc_library),
        PetscErrorCode,
        (IS, MPI_Comm, PetscCopyMode, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function ISRenumber(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:ISRenumber, $petsc_library),
        PetscErrorCode,
        (IS, IS, Ptr{$PetscInt}, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function ISCreateSubIS(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:ISCreateSubIS, $petsc_library),
        PetscErrorCode,
        (IS, IS, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function ISGeneralFilter(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:ISGeneralFilter, $petsc_library),
        PetscErrorCode,
        (IS, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@enum ISGlobalToLocalMappingMode::UInt32 begin
    IS_GTOLM_MASK = 0
    IS_GTOLM_DROP = 1
end

const ISLocalToGlobalMappingType = Ptr{Cchar}

@for_petsc function ISLocalToGlobalMappingSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISLocalToGlobalMappingSetType, $petsc_library),
        PetscErrorCode,
        (ISLocalToGlobalMapping, ISLocalToGlobalMappingType),
        arg1,
        arg2,
    )
end

@for_petsc function ISLocalToGlobalMappingRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISLocalToGlobalMappingRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function ISLocalToGlobalMappingRegisterAll(::$UnionPetscLib)
    @chk ccall(
        (:ISLocalToGlobalMappingRegisterAll, $petsc_library),
        PetscErrorCode,
        (),
    )
end

@for_petsc function ISLocalToGlobalMappingCreate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:ISLocalToGlobalMappingCreate, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            PetscCopyMode,
            Ptr{ISLocalToGlobalMapping},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function ISLocalToGlobalMappingCreateIS(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISLocalToGlobalMappingCreateIS, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{ISLocalToGlobalMapping}),
        arg1,
        arg2,
    )
end

@for_petsc function ISLocalToGlobalMappingCreateSF(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:ISLocalToGlobalMappingCreateSF, $petsc_library),
        PetscErrorCode,
        (PetscSF, $PetscInt, Ptr{ISLocalToGlobalMapping}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function ISLocalToGlobalMappingSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:ISLocalToGlobalMappingSetFromOptions, $petsc_library),
        PetscErrorCode,
        (ISLocalToGlobalMapping,),
        arg1,
    )
end

@for_petsc function ISLocalToGlobalMappingSetUp(::$UnionPetscLib, arg1)
    @chk ccall(
        (:ISLocalToGlobalMappingSetUp, $petsc_library),
        PetscErrorCode,
        (ISLocalToGlobalMapping,),
        arg1,
    )
end

@for_petsc function ISLocalToGlobalMappingView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISLocalToGlobalMappingView, $petsc_library),
        PetscErrorCode,
        (ISLocalToGlobalMapping, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function ISLocalToGlobalMappingViewFromOptions(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:ISLocalToGlobalMappingViewFromOptions, $petsc_library),
        PetscErrorCode,
        (ISLocalToGlobalMapping, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function ISLocalToGlobalMappingDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:ISLocalToGlobalMappingDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{ISLocalToGlobalMapping},),
        arg1,
    )
end

@for_petsc function ISLocalToGlobalMappingApply(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:ISLocalToGlobalMappingApply, $petsc_library),
        PetscErrorCode,
        (ISLocalToGlobalMapping, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function ISLocalToGlobalMappingApplyBlock(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:ISLocalToGlobalMappingApplyBlock, $petsc_library),
        PetscErrorCode,
        (ISLocalToGlobalMapping, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function ISLocalToGlobalMappingApplyIS(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:ISLocalToGlobalMappingApplyIS, $petsc_library),
        PetscErrorCode,
        (ISLocalToGlobalMapping, IS, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function ISGlobalToLocalMappingApply(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:ISGlobalToLocalMappingApply, $petsc_library),
        PetscErrorCode,
        (
            ISLocalToGlobalMapping,
            ISGlobalToLocalMappingMode,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function ISGlobalToLocalMappingApplyBlock(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:ISGlobalToLocalMappingApplyBlock, $petsc_library),
        PetscErrorCode,
        (
            ISLocalToGlobalMapping,
            ISGlobalToLocalMappingMode,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function ISGlobalToLocalMappingApplyIS(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:ISGlobalToLocalMappingApplyIS, $petsc_library),
        PetscErrorCode,
        (ISLocalToGlobalMapping, ISGlobalToLocalMappingMode, IS, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function ISLocalToGlobalMappingGetSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISLocalToGlobalMappingGetSize, $petsc_library),
        PetscErrorCode,
        (ISLocalToGlobalMapping, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function ISLocalToGlobalMappingGetNodeInfo(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:ISLocalToGlobalMappingGetNodeInfo, $petsc_library),
        PetscErrorCode,
        (
            ISLocalToGlobalMapping,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{Ptr{$PetscInt}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function ISLocalToGlobalMappingRestoreNodeInfo(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:ISLocalToGlobalMappingRestoreNodeInfo, $petsc_library),
        PetscErrorCode,
        (
            ISLocalToGlobalMapping,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{Ptr{$PetscInt}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function ISLocalToGlobalMappingGetInfo(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:ISLocalToGlobalMappingGetInfo, $petsc_library),
        PetscErrorCode,
        (
            ISLocalToGlobalMapping,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{Ptr{$PetscInt}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function ISLocalToGlobalMappingRestoreInfo(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:ISLocalToGlobalMappingRestoreInfo, $petsc_library),
        PetscErrorCode,
        (
            ISLocalToGlobalMapping,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{Ptr{$PetscInt}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function ISLocalToGlobalMappingGetBlockInfo(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:ISLocalToGlobalMappingGetBlockInfo, $petsc_library),
        PetscErrorCode,
        (
            ISLocalToGlobalMapping,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{Ptr{$PetscInt}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function ISLocalToGlobalMappingRestoreBlockInfo(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:ISLocalToGlobalMappingRestoreBlockInfo, $petsc_library),
        PetscErrorCode,
        (
            ISLocalToGlobalMapping,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{Ptr{$PetscInt}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function ISLocalToGlobalMappingGetIndices(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:ISLocalToGlobalMappingGetIndices, $petsc_library),
        PetscErrorCode,
        (ISLocalToGlobalMapping, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function ISLocalToGlobalMappingRestoreIndices(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:ISLocalToGlobalMappingRestoreIndices, $petsc_library),
        PetscErrorCode,
        (ISLocalToGlobalMapping, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function ISLocalToGlobalMappingGetBlockIndices(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:ISLocalToGlobalMappingGetBlockIndices, $petsc_library),
        PetscErrorCode,
        (ISLocalToGlobalMapping, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function ISLocalToGlobalMappingRestoreBlockIndices(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:ISLocalToGlobalMappingRestoreBlockIndices, $petsc_library),
        PetscErrorCode,
        (ISLocalToGlobalMapping, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function ISLocalToGlobalMappingConcatenate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:ISLocalToGlobalMappingConcatenate, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            Ptr{ISLocalToGlobalMapping},
            Ptr{ISLocalToGlobalMapping},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function ISLocalToGlobalMappingGetBlockSize(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:ISLocalToGlobalMappingGetBlockSize, $petsc_library),
        PetscErrorCode,
        (ISLocalToGlobalMapping, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function ISLocalToGlobalMappingSetBlockSize(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:ISLocalToGlobalMappingSetBlockSize, $petsc_library),
        PetscErrorCode,
        (ISLocalToGlobalMapping, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function ISLocalToGlobalMappingDuplicate(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:ISLocalToGlobalMappingDuplicate, $petsc_library),
        PetscErrorCode,
        (ISLocalToGlobalMapping, Ptr{ISLocalToGlobalMapping}),
        arg1,
        arg2,
    )
end

@enum ISColoringType::UInt32 begin
    IS_COLORING_GLOBAL = 0
    IS_COLORING_LOCAL = 1
end

const ISColoringValue = Cushort

@for_petsc function ISAllGatherColors(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:ISAllGatherColors, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            Ptr{ISColoringValue},
            Ptr{$PetscInt},
            Ptr{Ptr{ISColoringValue}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function ISColoringCreate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:ISColoringCreate, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            Ptr{ISColoringValue},
            PetscCopyMode,
            Ptr{ISColoring},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function ISColoringDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:ISColoringDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{ISColoring},),
        arg1,
    )
end

@for_petsc function ISColoringView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISColoringView, $petsc_library),
        PetscErrorCode,
        (ISColoring, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function ISColoringViewFromOptions(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:ISColoringViewFromOptions, $petsc_library),
        PetscErrorCode,
        (ISColoring, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function ISColoringGetIS(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:ISColoringGetIS, $petsc_library),
        PetscErrorCode,
        (ISColoring, PetscCopyMode, Ptr{$PetscInt}, Ptr{Ptr{IS}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function ISColoringRestoreIS(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:ISColoringRestoreIS, $petsc_library),
        PetscErrorCode,
        (ISColoring, PetscCopyMode, Ptr{Ptr{IS}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function ISColoringReference(::$UnionPetscLib, arg1)
    @chk ccall(
        (:ISColoringReference, $petsc_library),
        PetscErrorCode,
        (ISColoring,),
        arg1,
    )
end

@for_petsc function ISColoringSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISColoringSetType, $petsc_library),
        PetscErrorCode,
        (ISColoring, ISColoringType),
        arg1,
        arg2,
    )
end

@for_petsc function ISColoringGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISColoringGetType, $petsc_library),
        PetscErrorCode,
        (ISColoring, Ptr{ISColoringType}),
        arg1,
        arg2,
    )
end

@for_petsc function ISColoringGetColors(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:ISColoringGetColors, $petsc_library),
        PetscErrorCode,
        (ISColoring, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{ISColoringValue}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function ISBuildTwoSided(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:ISBuildTwoSided, $petsc_library),
        PetscErrorCode,
        (IS, IS, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function ISPartitioningToNumbering(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISPartitioningToNumbering, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{IS}),
        arg1,
        arg2,
    )
end

@for_petsc function ISPartitioningCount(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:ISPartitioningCount, $petsc_library),
        PetscErrorCode,
        (IS, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function ISCompressIndicesGeneral(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:ISCompressIndicesGeneral, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{IS}, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function ISCompressIndicesSorted(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:ISCompressIndicesSorted, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, $PetscInt, Ptr{IS}, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function ISExpandIndicesGeneral(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:ISExpandIndicesGeneral, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{IS}, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@enum ScatterMode::UInt32 begin
    SCATTER_FORWARD = 0
    SCATTER_REVERSE = 1
    SCATTER_FORWARD_LOCAL = 2
    SCATTER_REVERSE_LOCAL = 3
    # SCATTER_LOCAL = 2
end

@for_petsc function VecScatterSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecScatterSetType, $petsc_library),
        PetscErrorCode,
        (VecScatter, VecScatterType),
        arg1,
        arg2,
    )
end

@for_petsc function VecScatterGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecScatterGetType, $petsc_library),
        PetscErrorCode,
        (VecScatter, Ptr{VecScatterType}),
        arg1,
        arg2,
    )
end

@for_petsc function VecScatterSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:VecScatterSetFromOptions, $petsc_library),
        PetscErrorCode,
        (VecScatter,),
        arg1,
    )
end

@for_petsc function VecScatterRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecScatterRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function VecScatterCreate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:VecScatterCreate, $petsc_library),
        PetscErrorCode,
        (Vec, IS, Vec, IS, Ptr{VecScatter}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function VecInitializePackage(::$UnionPetscLib)
    @chk ccall((:VecInitializePackage, $petsc_library), PetscErrorCode, ())
end

@for_petsc function VecFinalizePackage(::$UnionPetscLib)
    @chk ccall((:VecFinalizePackage, $petsc_library), PetscErrorCode, ())
end

@for_petsc function VecCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function VecCreateSeq(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecCreateSeq, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecCreateMPI(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecCreateMPI, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, $PetscInt, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecCreateSeqWithArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:VecCreateSeqWithArray, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function VecCreateMPIWithArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:VecCreateMPIWithArray, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscScalar},
            Ptr{Vec},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function VecCreateShared(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecCreateShared, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, $PetscInt, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecCreateNode(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecCreateNode, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, $PetscInt, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:VecSetFromOptions, $petsc_library),
        PetscErrorCode,
        (Vec,),
        arg1,
    )
end

@for_petsc function VecViewFromOptions(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecViewFromOptions, $petsc_library),
        PetscErrorCode,
        (Vec, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecSetUp(::$UnionPetscLib, arg1)
    @chk ccall((:VecSetUp, $petsc_library), PetscErrorCode, (Vec,), arg1)
end

@for_petsc function VecDestroy(::$UnionPetscLib, arg1)
    @chk ccall((:VecDestroy, $petsc_library), PetscErrorCode, (Ptr{Vec},), arg1)
end

@for_petsc function VecZeroEntries(::$UnionPetscLib, arg1)
    @chk ccall((:VecZeroEntries, $petsc_library), PetscErrorCode, (Vec,), arg1)
end

@for_petsc function VecSetOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecSetOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function VecAppendOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecAppendOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function VecGetOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecGetOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function VecSetSizes(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecSetSizes, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecDotNorm2(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecDotNorm2, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Ptr{$PetscScalar}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecDot(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecDot, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecDotRealPart(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecDotRealPart, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecTDot(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecTDot, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecMDot(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecMDot, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Ptr{Vec}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecMTDot(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecMTDot, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Ptr{Vec}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecGetSubVector(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecGetSubVector, $petsc_library),
        PetscErrorCode,
        (Vec, IS, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecRestoreSubVector(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecRestoreSubVector, $petsc_library),
        PetscErrorCode,
        (Vec, IS, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecConcatenate(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecConcatenate, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{Vec}, Ptr{Vec}, Ptr{Ptr{IS}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@enum NormType::UInt32 begin
    NORM_1 = 0
    NORM_2 = 1
    NORM_FROBENIUS = 2
    NORM_INFINITY = 3
    NORM_1_AND_2 = 4
end

@for_petsc function VecNorm(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecNorm, $petsc_library),
        PetscErrorCode,
        (Vec, NormType, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecNormAvailable(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecNormAvailable, $petsc_library),
        PetscErrorCode,
        (Vec, NormType, Ptr{PetscBool}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecNormalize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecNormalize, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function VecSum(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecSum, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscScalar}),
        arg1,
        arg2,
    )
end

@for_petsc function VecMax(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecMax, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscInt}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecMin(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecMin, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscInt}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecScale(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecScale, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscScalar),
        arg1,
        arg2,
    )
end

@for_petsc function VecCopy(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecCopy, $petsc_library),
        PetscErrorCode,
        (Vec, Vec),
        arg1,
        arg2,
    )
end

@for_petsc function VecSetRandom(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecSetRandom, $petsc_library),
        PetscErrorCode,
        (Vec, PetscRandom),
        arg1,
        arg2,
    )
end

@for_petsc function VecSet(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecSet, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscScalar),
        arg1,
        arg2,
    )
end

@for_petsc function VecSetInf(::$UnionPetscLib, arg1)
    @chk ccall((:VecSetInf, $petsc_library), PetscErrorCode, (Vec,), arg1)
end

@for_petsc function VecSwap(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecSwap, $petsc_library),
        PetscErrorCode,
        (Vec, Vec),
        arg1,
        arg2,
    )
end

@for_petsc function VecAXPY(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecAXPY, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscScalar, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecAXPBY(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecAXPBY, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscScalar, $PetscScalar, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecMAXPY(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecMAXPY, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Ptr{$PetscScalar}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecAYPX(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecAYPX, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscScalar, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecWAXPY(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecWAXPY, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscScalar, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecAXPBYPCZ(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:VecAXPBYPCZ, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscScalar, $PetscScalar, $PetscScalar, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function VecPointwiseMax(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecPointwiseMax, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecPointwiseMaxAbs(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecPointwiseMaxAbs, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecPointwiseMin(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecPointwiseMin, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecPointwiseMult(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecPointwiseMult, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecPointwiseDivide(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecPointwiseDivide, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecMaxPointwiseDivide(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecMaxPointwiseDivide, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecShift(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecShift, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscScalar),
        arg1,
        arg2,
    )
end

@for_petsc function VecReciprocal(::$UnionPetscLib, arg1)
    @chk ccall((:VecReciprocal, $petsc_library), PetscErrorCode, (Vec,), arg1)
end

@for_petsc function VecPermute(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecPermute, $petsc_library),
        PetscErrorCode,
        (Vec, IS, PetscBool),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecSqrtAbs(::$UnionPetscLib, arg1)
    @chk ccall((:VecSqrtAbs, $petsc_library), PetscErrorCode, (Vec,), arg1)
end

@for_petsc function VecLog(::$UnionPetscLib, arg1)
    @chk ccall((:VecLog, $petsc_library), PetscErrorCode, (Vec,), arg1)
end

@for_petsc function VecExp(::$UnionPetscLib, arg1)
    @chk ccall((:VecExp, $petsc_library), PetscErrorCode, (Vec,), arg1)
end

@for_petsc function VecAbs(::$UnionPetscLib, arg1)
    @chk ccall((:VecAbs, $petsc_library), PetscErrorCode, (Vec,), arg1)
end

@for_petsc function VecDuplicate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecDuplicate, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function VecDuplicateVecs(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecDuplicateVecs, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Ptr{Ptr{Vec}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecDestroyVecs(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecDestroyVecs, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{Ptr{Vec}}),
        arg1,
        arg2,
    )
end

@for_petsc function VecStrideNormAll(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecStrideNormAll, $petsc_library),
        PetscErrorCode,
        (Vec, NormType, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecStrideMaxAll(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecStrideMaxAll, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscInt}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecStrideMinAll(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecStrideMinAll, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscInt}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecStrideScaleAll(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecStrideScaleAll, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscScalar}),
        arg1,
        arg2,
    )
end

@for_petsc function VecUniqueEntries(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecUniqueEntries, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscInt}, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecStrideNorm(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecStrideNorm, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, NormType, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecStrideMax(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecStrideMax, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecStrideMin(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecStrideMin, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecStrideScale(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecStrideScale, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, $PetscScalar),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecStrideSet(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecStrideSet, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, $PetscScalar),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecStrideGather(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecStrideGather, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Vec, InsertMode),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecStrideScatter(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecStrideScatter, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Vec, InsertMode),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecStrideGatherAll(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecStrideGatherAll, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Vec}, InsertMode),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecStrideScatterAll(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecStrideScatterAll, $petsc_library),
        PetscErrorCode,
        (Ptr{Vec}, Vec, InsertMode),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecStrideSubSetScatter(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:VecStrideSubSetScatter, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Vec, InsertMode),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function VecStrideSubSetGather(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:VecStrideSubSetGather, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Vec, InsertMode),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function VecSetValues(::$UnionPetscLib, arg1, arg2, arg3, arg4, arg5)
    @chk ccall(
        (:VecSetValues, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, InsertMode),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function VecGetValues(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecGetValues, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecAssemblyBegin(::$UnionPetscLib, arg1)
    @chk ccall(
        (:VecAssemblyBegin, $petsc_library),
        PetscErrorCode,
        (Vec,),
        arg1,
    )
end

@for_petsc function VecAssemblyEnd(::$UnionPetscLib, arg1)
    @chk ccall((:VecAssemblyEnd, $petsc_library), PetscErrorCode, (Vec,), arg1)
end

@for_petsc function VecStashSetInitialSize(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecStashSetInitialSize, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecStashView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecStashView, $petsc_library),
        PetscErrorCode,
        (Vec, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function VecStashViewFromOptions(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecStashViewFromOptions, $petsc_library),
        PetscErrorCode,
        (Vec, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecStashGetInfo(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:VecStashGetInfo, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function VecSetValue(::$UnionPetscLib, v, i, va, mode)
    @chk ccall(
        (:VecSetValue, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, $PetscScalar, InsertMode),
        v,
        i,
        va,
        mode,
    )
end

@for_petsc function VecSetBlockSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecSetBlockSize, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function VecGetBlockSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecGetBlockSize, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function VecSetValuesBlocked(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:VecSetValuesBlocked, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, InsertMode),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function VecSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecSetType, $petsc_library),
        PetscErrorCode,
        (Vec, VecType),
        arg1,
        arg2,
    )
end

@for_petsc function VecGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecGetType, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{VecType}),
        arg1,
        arg2,
    )
end

@for_petsc function VecRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function VecScatterBegin(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:VecScatterBegin, $petsc_library),
        PetscErrorCode,
        (VecScatter, Vec, Vec, InsertMode, ScatterMode),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function VecScatterEnd(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:VecScatterEnd, $petsc_library),
        PetscErrorCode,
        (VecScatter, Vec, Vec, InsertMode, ScatterMode),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function VecScatterDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:VecScatterDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{VecScatter},),
        arg1,
    )
end

@for_petsc function VecScatterSetUp(::$UnionPetscLib, arg1)
    @chk ccall(
        (:VecScatterSetUp, $petsc_library),
        PetscErrorCode,
        (VecScatter,),
        arg1,
    )
end

@for_petsc function VecScatterCopy(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecScatterCopy, $petsc_library),
        PetscErrorCode,
        (VecScatter, Ptr{VecScatter}),
        arg1,
        arg2,
    )
end

@for_petsc function VecScatterView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecScatterView, $petsc_library),
        PetscErrorCode,
        (VecScatter, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function VecScatterViewFromOptions(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:VecScatterViewFromOptions, $petsc_library),
        PetscErrorCode,
        (VecScatter, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecScatterRemap(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecScatterRemap, $petsc_library),
        PetscErrorCode,
        (VecScatter, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecScatterGetMerged(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecScatterGetMerged, $petsc_library),
        PetscErrorCode,
        (VecScatter, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function VecGetArray4d(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
)
    @chk ccall(
        (:VecGetArray4d, $petsc_library),
        PetscErrorCode,
        (
            Vec,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
    )
end

@for_petsc function VecRestoreArray4d(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
)
    @chk ccall(
        (:VecRestoreArray4d, $petsc_library),
        PetscErrorCode,
        (
            Vec,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
    )
end

@for_petsc function VecGetArray3d(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:VecGetArray3d, $petsc_library),
        PetscErrorCode,
        (
            Vec,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function VecRestoreArray3d(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:VecRestoreArray3d, $petsc_library),
        PetscErrorCode,
        (
            Vec,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function VecGetArray2d(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:VecGetArray2d, $petsc_library),
        PetscErrorCode,
        (
            Vec,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Ptr{$PetscScalar}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function VecRestoreArray2d(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:VecRestoreArray2d, $petsc_library),
        PetscErrorCode,
        (
            Vec,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Ptr{$PetscScalar}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function VecGetArray1d(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecGetArray1d, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, $PetscInt, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecRestoreArray1d(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecRestoreArray1d, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, $PetscInt, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecGetArray4dWrite(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
)
    @chk ccall(
        (:VecGetArray4dWrite, $petsc_library),
        PetscErrorCode,
        (
            Vec,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
    )
end

@for_petsc function VecRestoreArray4dWrite(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
)
    @chk ccall(
        (:VecRestoreArray4dWrite, $petsc_library),
        PetscErrorCode,
        (
            Vec,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
    )
end

@for_petsc function VecGetArray3dWrite(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:VecGetArray3dWrite, $petsc_library),
        PetscErrorCode,
        (
            Vec,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function VecRestoreArray3dWrite(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:VecRestoreArray3dWrite, $petsc_library),
        PetscErrorCode,
        (
            Vec,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function VecGetArray2dWrite(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:VecGetArray2dWrite, $petsc_library),
        PetscErrorCode,
        (
            Vec,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Ptr{$PetscScalar}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function VecRestoreArray2dWrite(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:VecRestoreArray2dWrite, $petsc_library),
        PetscErrorCode,
        (
            Vec,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Ptr{$PetscScalar}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function VecGetArray1dWrite(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecGetArray1dWrite, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, $PetscInt, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecRestoreArray1dWrite(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:VecRestoreArray1dWrite, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, $PetscInt, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecGetArray4dRead(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
)
    @chk ccall(
        (:VecGetArray4dRead, $petsc_library),
        PetscErrorCode,
        (
            Vec,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
    )
end

@for_petsc function VecRestoreArray4dRead(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
)
    @chk ccall(
        (:VecRestoreArray4dRead, $petsc_library),
        PetscErrorCode,
        (
            Vec,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
    )
end

@for_petsc function VecGetArray3dRead(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:VecGetArray3dRead, $petsc_library),
        PetscErrorCode,
        (
            Vec,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function VecRestoreArray3dRead(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:VecRestoreArray3dRead, $petsc_library),
        PetscErrorCode,
        (
            Vec,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function VecGetArray2dRead(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:VecGetArray2dRead, $petsc_library),
        PetscErrorCode,
        (
            Vec,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Ptr{$PetscScalar}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function VecRestoreArray2dRead(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:VecRestoreArray2dRead, $petsc_library),
        PetscErrorCode,
        (
            Vec,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Ptr{$PetscScalar}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function VecGetArray1dRead(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecGetArray1dRead, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, $PetscInt, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecRestoreArray1dRead(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:VecRestoreArray1dRead, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, $PetscInt, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecPlaceArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecPlaceArray, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscScalar}),
        arg1,
        arg2,
    )
end

@for_petsc function VecResetArray(::$UnionPetscLib, arg1)
    @chk ccall((:VecResetArray, $petsc_library), PetscErrorCode, (Vec,), arg1)
end

@for_petsc function VecReplaceArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecReplaceArray, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscScalar}),
        arg1,
        arg2,
    )
end

@for_petsc function VecGetArrays(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecGetArrays, $petsc_library),
        PetscErrorCode,
        (Ptr{Vec}, $PetscInt, Ptr{Ptr{Ptr{$PetscScalar}}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecRestoreArrays(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecRestoreArrays, $petsc_library),
        PetscErrorCode,
        (Ptr{Vec}, $PetscInt, Ptr{Ptr{Ptr{$PetscScalar}}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecView, $petsc_library),
        PetscErrorCode,
        (Vec, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function VecEqual(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecEqual, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecLoad(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecLoad, $petsc_library),
        PetscErrorCode,
        (Vec, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function VecGetSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecGetSize, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function VecGetLocalSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecGetLocalSize, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function VecGetOwnershipRange(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecGetOwnershipRange, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecGetOwnershipRanges(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecGetOwnershipRanges, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function VecSetLocalToGlobalMapping(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecSetLocalToGlobalMapping, $petsc_library),
        PetscErrorCode,
        (Vec, ISLocalToGlobalMapping),
        arg1,
        arg2,
    )
end

@for_petsc function VecSetValuesLocal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:VecSetValuesLocal, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, InsertMode),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function VecCUDAGetArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecCUDAGetArray, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function VecCUDARestoreArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecCUDARestoreArray, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function VecCUDAGetArrayRead(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecCUDAGetArrayRead, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function VecCUDARestoreArrayRead(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecCUDARestoreArrayRead, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function VecCUDAGetArrayWrite(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecCUDAGetArrayWrite, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function VecCUDARestoreArrayWrite(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecCUDARestoreArrayWrite, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function VecCUDAPlaceArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecCUDAPlaceArray, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscScalar}),
        arg1,
        arg2,
    )
end

@for_petsc function VecCUDAReplaceArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecCUDAReplaceArray, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscScalar}),
        arg1,
        arg2,
    )
end

@for_petsc function VecCUDAResetArray(::$UnionPetscLib, arg1)
    @chk ccall(
        (:VecCUDAResetArray, $petsc_library),
        PetscErrorCode,
        (Vec,),
        arg1,
    )
end

@for_petsc function VecHIPGetArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecHIPGetArray, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function VecHIPRestoreArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecHIPRestoreArray, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function VecHIPGetArrayRead(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecHIPGetArrayRead, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function VecHIPRestoreArrayRead(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecHIPRestoreArrayRead, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function VecHIPGetArrayWrite(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecHIPGetArrayWrite, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function VecHIPRestoreArrayWrite(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecHIPRestoreArrayWrite, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function VecHIPPlaceArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecHIPPlaceArray, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscScalar}),
        arg1,
        arg2,
    )
end

@for_petsc function VecHIPReplaceArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecHIPReplaceArray, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscScalar}),
        arg1,
        arg2,
    )
end

@for_petsc function VecHIPResetArray(::$UnionPetscLib, arg1)
    @chk ccall(
        (:VecHIPResetArray, $petsc_library),
        PetscErrorCode,
        (Vec,),
        arg1,
    )
end

@for_petsc function VecViennaCLGetCLContext(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecViennaCLGetCLContext, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Csize_t}),
        arg1,
        arg2,
    )
end

@for_petsc function VecViennaCLGetCLQueue(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecViennaCLGetCLQueue, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Csize_t}),
        arg1,
        arg2,
    )
end

@for_petsc function VecViennaCLGetCLMemRead(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecViennaCLGetCLMemRead, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Csize_t}),
        arg1,
        arg2,
    )
end

@for_petsc function VecViennaCLGetCLMemWrite(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecViennaCLGetCLMemWrite, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Csize_t}),
        arg1,
        arg2,
    )
end

@for_petsc function VecViennaCLRestoreCLMemWrite(::$UnionPetscLib, arg1)
    @chk ccall(
        (:VecViennaCLRestoreCLMemWrite, $petsc_library),
        PetscErrorCode,
        (Vec,),
        arg1,
    )
end

@for_petsc function VecViennaCLGetCLMem(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecViennaCLGetCLMem, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Csize_t}),
        arg1,
        arg2,
    )
end

@for_petsc function VecViennaCLRestoreCLMem(::$UnionPetscLib, arg1)
    @chk ccall(
        (:VecViennaCLRestoreCLMem, $petsc_library),
        PetscErrorCode,
        (Vec,),
        arg1,
    )
end

@for_petsc function VecSetValueLocal(::$UnionPetscLib, v, i, va, mode)
    @chk ccall(
        (:VecSetValueLocal, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, $PetscScalar, InsertMode),
        v,
        i,
        va,
        mode,
    )
end

@for_petsc function VecSetValuesBlockedLocal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:VecSetValuesBlockedLocal, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, InsertMode),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function VecGetLocalToGlobalMapping(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecGetLocalToGlobalMapping, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{ISLocalToGlobalMapping}),
        arg1,
        arg2,
    )
end

@for_petsc function VecDotBegin(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecDotBegin, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecDotEnd(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecDotEnd, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecTDotBegin(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecTDotBegin, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecTDotEnd(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecTDotEnd, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecNormBegin(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecNormBegin, $petsc_library),
        PetscErrorCode,
        (Vec, NormType, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecNormEnd(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecNormEnd, $petsc_library),
        PetscErrorCode,
        (Vec, NormType, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecMDotBegin(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecMDotBegin, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Ptr{Vec}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecMDotEnd(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecMDotEnd, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Ptr{Vec}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecMTDotBegin(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecMTDotBegin, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Ptr{Vec}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecMTDotEnd(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecMTDotEnd, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Ptr{Vec}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscCommSplitReductionBegin(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscCommSplitReductionBegin, $petsc_library),
        PetscErrorCode,
        (MPI_Comm,),
        arg1,
    )
end

@for_petsc function VecBindToCPU(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecBindToCPU, $petsc_library),
        PetscErrorCode,
        (Vec, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function VecPinToCPU(::$UnionPetscLib, v, flg)
    @chk ccall(
        (:VecPinToCPU, $petsc_library),
        PetscErrorCode,
        (Vec, PetscBool),
        v,
        flg,
    )
end

@for_petsc function VecSetPinnedMemoryMin(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecSetPinnedMemoryMin, $petsc_library),
        PetscErrorCode,
        (Vec, Csize_t),
        arg1,
        arg2,
    )
end

@for_petsc function VecGetPinnedMemoryMin(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecGetPinnedMemoryMin, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Csize_t}),
        arg1,
        arg2,
    )
end

@enum PetscOffloadMask::UInt32 begin
    PETSC_OFFLOAD_UNALLOCATED = 0
    PETSC_OFFLOAD_CPU = 1
    PETSC_OFFLOAD_GPU = 2
    PETSC_OFFLOAD_BOTH = 3
    PETSC_OFFLOAD_VECKOKKOS = 256
end

@for_petsc function VecGetOffloadMask(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecGetOffloadMask, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{PetscOffloadMask}),
        arg1,
        arg2,
    )
end

@enum VecOption::UInt32 begin
    VEC_IGNORE_OFF_PROC_ENTRIES = 0
    VEC_IGNORE_NEGATIVE_INDICES = 1
    VEC_SUBSET_OFF_PROC_ENTRIES = 2
end

@for_petsc function VecSetOption(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecSetOption, $petsc_library),
        PetscErrorCode,
        (Vec, VecOption, PetscBool),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecGetArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecGetArray, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function VecGetArrayWrite(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecGetArrayWrite, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function VecGetArrayRead(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecGetArrayRead, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function VecRestoreArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecRestoreArray, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function VecRestoreArrayWrite(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecRestoreArrayWrite, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function VecRestoreArrayRead(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecRestoreArrayRead, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function VecGetLocalVector(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecGetLocalVector, $petsc_library),
        PetscErrorCode,
        (Vec, Vec),
        arg1,
        arg2,
    )
end

@for_petsc function VecRestoreLocalVector(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecRestoreLocalVector, $petsc_library),
        PetscErrorCode,
        (Vec, Vec),
        arg1,
        arg2,
    )
end

@for_petsc function VecGetLocalVectorRead(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecGetLocalVectorRead, $petsc_library),
        PetscErrorCode,
        (Vec, Vec),
        arg1,
        arg2,
    )
end

@for_petsc function VecRestoreLocalVectorRead(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecRestoreLocalVectorRead, $petsc_library),
        PetscErrorCode,
        (Vec, Vec),
        arg1,
        arg2,
    )
end

@for_petsc function VecGetArrayAndMemType(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecGetArrayAndMemType, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}, Ptr{PetscMemType}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecRestoreArrayAndMemType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecRestoreArrayAndMemType, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function VecGetArrayReadAndMemType(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:VecGetArrayReadAndMemType, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}, Ptr{PetscMemType}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecRestoreArrayReadAndMemType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecRestoreArrayReadAndMemType, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function VecGetArrayPair(::$UnionPetscLib, x, y, xv, yv)
    @chk ccall(
        (:VecGetArrayPair, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Ptr{Ptr{$PetscScalar}}, Ptr{Ptr{$PetscScalar}}),
        x,
        y,
        xv,
        yv,
    )
end

@for_petsc function VecRestoreArrayPair(::$UnionPetscLib, x, y, xv, yv)
    @chk ccall(
        (:VecRestoreArrayPair, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Ptr{Ptr{$PetscScalar}}, Ptr{Ptr{$PetscScalar}}),
        x,
        y,
        xv,
        yv,
    )
end

@for_petsc function VecValidValues(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecValidValues, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, PetscBool),
        arg1,
        arg2,
        arg3,
    )
end

@enum VecOperation::UInt32 begin
    VECOP_DUPLICATE = 0
    VECOP_VIEW = 33
    VECOP_LOAD = 41
    VECOP_VIEWNATIVE = 68
    VECOP_LOADNATIVE = 69
end

@for_petsc function VecSetOperation(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecSetOperation, $petsc_library),
        PetscErrorCode,
        (Vec, VecOperation, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecMPISetGhost(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecMPISetGhost, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecCreateGhost(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:VecCreateGhost, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function VecCreateGhostWithArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:VecCreateGhostWithArray, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscScalar},
            Ptr{Vec},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function VecCreateGhostBlock(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:VecCreateGhostBlock, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Vec},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function VecCreateGhostBlockWithArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:VecCreateGhostBlockWithArray, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscScalar},
            Ptr{Vec},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function VecGhostGetLocalForm(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecGhostGetLocalForm, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function VecGhostRestoreLocalForm(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecGhostRestoreLocalForm, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function VecGhostIsLocalForm(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecGhostIsLocalForm, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecGhostUpdateBegin(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecGhostUpdateBegin, $petsc_library),
        PetscErrorCode,
        (Vec, InsertMode, ScatterMode),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecGhostUpdateEnd(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecGhostUpdateEnd, $petsc_library),
        PetscErrorCode,
        (Vec, InsertMode, ScatterMode),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecConjugate(::$UnionPetscLib, arg1)
    @chk ccall((:VecConjugate, $petsc_library), PetscErrorCode, (Vec,), arg1)
end

@for_petsc function VecImaginaryPart(::$UnionPetscLib, arg1)
    @chk ccall(
        (:VecImaginaryPart, $petsc_library),
        PetscErrorCode,
        (Vec,),
        arg1,
    )
end

@for_petsc function VecRealPart(::$UnionPetscLib, arg1)
    @chk ccall((:VecRealPart, $petsc_library), PetscErrorCode, (Vec,), arg1)
end

@for_petsc function VecScatterCreateToAll(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecScatterCreateToAll, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{VecScatter}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecScatterCreateToZero(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecScatterCreateToZero, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{VecScatter}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function ISComplementVec(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:ISComplementVec, $petsc_library),
        PetscErrorCode,
        (IS, Vec, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecPow(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecPow, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscScalar),
        arg1,
        arg2,
    )
end

@for_petsc function VecMedian(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecMedian, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecWhichInactive(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:VecWhichInactive, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Vec, Vec, PetscBool, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function VecWhichBetween(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecWhichBetween, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Vec, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecWhichBetweenOrEqual(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:VecWhichBetweenOrEqual, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Vec, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecWhichGreaterThan(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecWhichGreaterThan, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecWhichLessThan(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecWhichLessThan, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecWhichEqual(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecWhichEqual, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecISAXPY(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecISAXPY, $petsc_library),
        PetscErrorCode,
        (Vec, IS, $PetscScalar, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecISCopy(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecISCopy, $petsc_library),
        PetscErrorCode,
        (Vec, IS, ScatterMode, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecISSet(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecISSet, $petsc_library),
        PetscErrorCode,
        (Vec, IS, $PetscScalar),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecBoundGradientProjection(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:VecBoundGradientProjection, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Vec, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function VecStepBoundInfo(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:VecStepBoundInfo, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Vec, Vec, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function VecStepMax(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecStepMax, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecStepMaxBounded(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:VecStepMaxBounded, $petsc_library),
        PetscErrorCode,
        (Vec, Vec, Vec, Vec, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscViewerMathematicaGetVector(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscViewerMathematicaGetVector, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Vec),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerMathematicaPutVector(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscViewerMathematicaPutVector, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Vec),
        arg1,
        arg2,
    )
end

@for_petsc function VecNestGetSubVecs(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecNestGetSubVecs, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscInt}, Ptr{Ptr{Vec}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecNestGetSubVec(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecNestGetSubVec, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecNestSetSubVecs(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecNestSetSubVecs, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Ptr{$PetscInt}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecNestSetSubVec(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecNestSetSubVec, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecCreateNest(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:VecCreateNest, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, Ptr{IS}, Ptr{Vec}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function VecNestGetSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecNestGetSize, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscOptionsGetVec(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscOptionsGetVec, $petsc_library),
        PetscErrorCode,
        (PetscOptions, Ptr{Cchar}, Ptr{Cchar}, Vec, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function VecChop(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecChop, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionVecView(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSectionVecView, $petsc_library),
        PetscErrorCode,
        (PetscSection, Vec, PetscViewer),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecGetValuesSection(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:VecGetValuesSection, $petsc_library),
        PetscErrorCode,
        (Vec, PetscSection, $PetscInt, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecSetValuesSection(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:VecSetValuesSection, $petsc_library),
        PetscErrorCode,
        (Vec, PetscSection, $PetscInt, Ptr{$PetscScalar}, InsertMode),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscSectionVecNorm(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscSectionVecNorm, $petsc_library),
        PetscErrorCode,
        (PetscSection, PetscSection, Vec, NormType, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

mutable struct _p_VecTagger end

const VecTagger = Ptr{_p_VecTagger}

const VecTaggerType = Ptr{Cchar}

@for_petsc function VecTaggerRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecTaggerRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function VecTaggerCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecTaggerCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{VecTagger}),
        arg1,
        arg2,
    )
end

@for_petsc function VecTaggerSetBlockSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecTaggerSetBlockSize, $petsc_library),
        PetscErrorCode,
        (VecTagger, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function VecTaggerGetBlockSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecTaggerGetBlockSize, $petsc_library),
        PetscErrorCode,
        (VecTagger, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function VecTaggerSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecTaggerSetType, $petsc_library),
        PetscErrorCode,
        (VecTagger, VecTaggerType),
        arg1,
        arg2,
    )
end

@for_petsc function VecTaggerGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecTaggerGetType, $petsc_library),
        PetscErrorCode,
        (VecTagger, Ptr{VecTaggerType}),
        arg1,
        arg2,
    )
end

@for_petsc function VecTaggerSetInvert(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecTaggerSetInvert, $petsc_library),
        PetscErrorCode,
        (VecTagger, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function VecTaggerGetInvert(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecTaggerGetInvert, $petsc_library),
        PetscErrorCode,
        (VecTagger, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function VecTaggerSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:VecTaggerSetFromOptions, $petsc_library),
        PetscErrorCode,
        (VecTagger,),
        arg1,
    )
end

@for_petsc function VecTaggerSetUp(::$UnionPetscLib, arg1)
    @chk ccall(
        (:VecTaggerSetUp, $petsc_library),
        PetscErrorCode,
        (VecTagger,),
        arg1,
    )
end

@for_petsc function VecTaggerView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecTaggerView, $petsc_library),
        PetscErrorCode,
        (VecTagger, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function VecTaggerComputeIS(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecTaggerComputeIS, $petsc_library),
        PetscErrorCode,
        (VecTagger, Vec, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecTaggerDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:VecTaggerDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{VecTagger},),
        arg1,
    )
end

@enum VecTaggerCDFMethod::UInt32 begin
    VECTAGGER_CDF_GATHER = 0
    VECTAGGER_CDF_ITERATIVE = 1
    VECTAGGER_CDF_NUM_METHODS = 2
end

@for_petsc function VecTaggerCDFSetMethod(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecTaggerCDFSetMethod, $petsc_library),
        PetscErrorCode,
        (VecTagger, VecTaggerCDFMethod),
        arg1,
        arg2,
    )
end

@for_petsc function VecTaggerCDFGetMethod(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecTaggerCDFGetMethod, $petsc_library),
        PetscErrorCode,
        (VecTagger, Ptr{VecTaggerCDFMethod}),
        arg1,
        arg2,
    )
end

@for_petsc function VecTaggerCDFIterativeSetTolerances(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:VecTaggerCDFIterativeSetTolerances, $petsc_library),
        PetscErrorCode,
        (VecTagger, $PetscInt, $PetscReal, $PetscReal),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecTaggerCDFIterativeGetTolerances(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:VecTaggerCDFIterativeGetTolerances, $petsc_library),
        PetscErrorCode,
        (VecTagger, Ptr{$PetscInt}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecTaggerOrSetSubs(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecTaggerOrSetSubs, $petsc_library),
        PetscErrorCode,
        (VecTagger, $PetscInt, Ptr{VecTagger}, PetscCopyMode),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecTaggerOrGetSubs(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecTaggerOrGetSubs, $petsc_library),
        PetscErrorCode,
        (VecTagger, Ptr{$PetscInt}, Ptr{Ptr{VecTagger}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecTaggerAndSetSubs(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:VecTaggerAndSetSubs, $petsc_library),
        PetscErrorCode,
        (VecTagger, $PetscInt, Ptr{VecTagger}, PetscCopyMode),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function VecTaggerAndGetSubs(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecTaggerAndGetSubs, $petsc_library),
        PetscErrorCode,
        (VecTagger, Ptr{$PetscInt}, Ptr{Ptr{VecTagger}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecTaggerInitializePackage(::$UnionPetscLib)
    @chk ccall(
        (:VecTaggerInitializePackage, $petsc_library),
        PetscErrorCode,
        (),
    )
end

@for_petsc function VecTaggerFinalizePackage(::$UnionPetscLib)
    @chk ccall((:VecTaggerFinalizePackage, $petsc_library), PetscErrorCode, ())
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

@for_petsc function PetscSFRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSFRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSFInitializePackage(::$UnionPetscLib)
    @chk ccall((:PetscSFInitializePackage, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscSFFinalizePackage(::$UnionPetscLib)
    @chk ccall((:PetscSFFinalizePackage, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscSFCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSFCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscSF}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSFDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscSFDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscSF},),
        arg1,
    )
end

@for_petsc function PetscSFSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSFSetType, $petsc_library),
        PetscErrorCode,
        (PetscSF, PetscSFType),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSFGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSFGetType, $petsc_library),
        PetscErrorCode,
        (PetscSF, Ptr{PetscSFType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSFView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSFView, $petsc_library),
        PetscErrorCode,
        (PetscSF, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSFViewFromOptions(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSFViewFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscSF, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSFSetUp(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscSFSetUp, $petsc_library),
        PetscErrorCode,
        (PetscSF,),
        arg1,
    )
end

@for_petsc function PetscSFSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscSFSetFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscSF,),
        arg1,
    )
end

@for_petsc function PetscSFDuplicate(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSFDuplicate, $petsc_library),
        PetscErrorCode,
        (PetscSF, PetscSFDuplicateOption, Ptr{PetscSF}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSFWindowSetSyncType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSFWindowSetSyncType, $petsc_library),
        PetscErrorCode,
        (PetscSF, PetscSFWindowSyncType),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSFWindowGetSyncType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSFWindowGetSyncType, $petsc_library),
        PetscErrorCode,
        (PetscSF, Ptr{PetscSFWindowSyncType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSFWindowSetFlavorType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSFWindowSetFlavorType, $petsc_library),
        PetscErrorCode,
        (PetscSF, PetscSFWindowFlavorType),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSFWindowGetFlavorType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSFWindowGetFlavorType, $petsc_library),
        PetscErrorCode,
        (PetscSF, Ptr{PetscSFWindowFlavorType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSFWindowSetInfo(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSFWindowSetInfo, $petsc_library),
        PetscErrorCode,
        (PetscSF, MPI_Info),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSFWindowGetInfo(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSFWindowGetInfo, $petsc_library),
        PetscErrorCode,
        (PetscSF, Ptr{MPI_Info}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSFSetRankOrder(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSFSetRankOrder, $petsc_library),
        PetscErrorCode,
        (PetscSF, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSFGetLeafRange(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSFGetLeafRange, $petsc_library),
        PetscErrorCode,
        (PetscSF, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSFCreateEmbeddedRootSF(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSFCreateEmbeddedRootSF, $petsc_library),
        PetscErrorCode,
        (PetscSF, $PetscInt, Ptr{$PetscInt}, Ptr{PetscSF}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSFCreateEmbeddedLeafSF(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSFCreateEmbeddedLeafSF, $petsc_library),
        PetscErrorCode,
        (PetscSF, $PetscInt, Ptr{$PetscInt}, Ptr{PetscSF}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSFReset(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscSFReset, $petsc_library),
        PetscErrorCode,
        (PetscSF,),
        arg1,
    )
end

@for_petsc function PetscSFGetRootRanks(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscSFGetRootRanks, $petsc_library),
        PetscErrorCode,
        (
            PetscSF,
            Ptr{$PetscInt},
            Ptr{Ptr{PetscMPIInt}},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscInt}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscSFGetLeafRanks(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscSFGetLeafRanks, $petsc_library),
        PetscErrorCode,
        (
            PetscSF,
            Ptr{$PetscInt},
            Ptr{Ptr{PetscMPIInt}},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscInt}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscSFGetMultiSF(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSFGetMultiSF, $petsc_library),
        PetscErrorCode,
        (PetscSF, Ptr{PetscSF}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSFCreateInverseSF(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSFCreateInverseSF, $petsc_library),
        PetscErrorCode,
        (PetscSF, Ptr{PetscSF}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSFSetGraphSection(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSFSetGraphSection, $petsc_library),
        PetscErrorCode,
        (PetscSF, PetscSection, PetscSection),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSFCreateRemoteOffsets(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSFCreateRemoteOffsets, $petsc_library),
        PetscErrorCode,
        (PetscSF, PetscSection, PetscSection, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSFDistributeSection(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSFDistributeSection, $petsc_library),
        PetscErrorCode,
        (PetscSF, PetscSection, Ptr{Ptr{$PetscInt}}, PetscSection),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSFCreateSectionSF(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscSFCreateSectionSF, $petsc_library),
        PetscErrorCode,
        (PetscSF, PetscSection, Ptr{$PetscInt}, PetscSection, Ptr{PetscSF}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscSFBcastBegin(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscSFBcastBegin, $petsc_library),
        PetscErrorCode,
        (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}, MPI_Op),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscSFBcastEnd(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscSFBcastEnd, $petsc_library),
        PetscErrorCode,
        (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}, MPI_Op),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscSFBcastWithMemTypeBegin(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscSFBcastWithMemTypeBegin, $petsc_library),
        PetscErrorCode,
        (
            PetscSF,
            MPI_Datatype,
            PetscMemType,
            Ptr{Cvoid},
            PetscMemType,
            Ptr{Cvoid},
            MPI_Op,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscSFReduceBegin(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscSFReduceBegin, $petsc_library),
        PetscErrorCode,
        (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}, MPI_Op),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscSFReduceEnd(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscSFReduceEnd, $petsc_library),
        PetscErrorCode,
        (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}, MPI_Op),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscSFReduceWithMemTypeBegin(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscSFReduceWithMemTypeBegin, $petsc_library),
        PetscErrorCode,
        (
            PetscSF,
            MPI_Datatype,
            PetscMemType,
            Ptr{Cvoid},
            PetscMemType,
            Ptr{Cvoid},
            MPI_Op,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscSFFetchAndOpBegin(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscSFFetchAndOpBegin, $petsc_library),
        PetscErrorCode,
        (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, MPI_Op),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscSFFetchAndOpEnd(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscSFFetchAndOpEnd, $petsc_library),
        PetscErrorCode,
        (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, MPI_Op),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscSFComputeDegreeBegin(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSFComputeDegreeBegin, $petsc_library),
        PetscErrorCode,
        (PetscSF, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSFComputeDegreeEnd(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSFComputeDegreeEnd, $petsc_library),
        PetscErrorCode,
        (PetscSF, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSFComputeMultiRootOriginalNumbering(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSFComputeMultiRootOriginalNumbering, $petsc_library),
        PetscErrorCode,
        (PetscSF, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSFGatherBegin(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscSFGatherBegin, $petsc_library),
        PetscErrorCode,
        (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSFGatherEnd(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscSFGatherEnd, $petsc_library),
        PetscErrorCode,
        (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSFScatterBegin(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSFScatterBegin, $petsc_library),
        PetscErrorCode,
        (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSFScatterEnd(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscSFScatterEnd, $petsc_library),
        PetscErrorCode,
        (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSFCompose(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSFCompose, $petsc_library),
        PetscErrorCode,
        (PetscSF, PetscSF, Ptr{PetscSF}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSFComposeInverse(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSFComposeInverse, $petsc_library),
        PetscErrorCode,
        (PetscSF, PetscSF, Ptr{PetscSF}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSFGetRanks(
    ::$UnionPetscLib,
    sf,
    nranks,
    ranks,
    roffset,
    rmine,
    rremote,
)
    @chk ccall(
        (:PetscSFGetRanks, $petsc_library),
        PetscErrorCode,
        (
            PetscSF,
            Ptr{$PetscInt},
            Ptr{Ptr{PetscMPIInt}},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscInt}},
        ),
        sf,
        nranks,
        ranks,
        roffset,
        rmine,
        rremote,
    )
end

@for_petsc function PetscSFCreateEmbeddedSF(
    ::$UnionPetscLib,
    sf,
    nselected,
    selected,
    esf,
)
    @chk ccall(
        (:PetscSFCreateEmbeddedSF, $petsc_library),
        PetscErrorCode,
        (PetscSF, $PetscInt, Ptr{$PetscInt}, Ptr{PetscSF}),
        sf,
        nselected,
        selected,
        esf,
    )
end

@for_petsc function PetscSFBcastAndOpBegin(
    ::$UnionPetscLib,
    sf,
    unit,
    rootdata,
    leafdata,
    op,
)
    @chk ccall(
        (:PetscSFBcastAndOpBegin, $petsc_library),
        PetscErrorCode,
        (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}, MPI_Op),
        sf,
        unit,
        rootdata,
        leafdata,
        op,
    )
end

@for_petsc function PetscSFBcastAndOpEnd(
    ::$UnionPetscLib,
    sf,
    unit,
    rootdata,
    leafdata,
    op,
)
    @chk ccall(
        (:PetscSFBcastAndOpEnd, $petsc_library),
        PetscErrorCode,
        (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}, MPI_Op),
        sf,
        unit,
        rootdata,
        leafdata,
        op,
    )
end

@for_petsc function PetscSFBcastAndOpWithMemtypeBegin(
    ::$UnionPetscLib,
    sf,
    unit,
    rootmtype,
    rootdata,
    leafmtype,
    leafdata,
    op,
)
    @chk ccall(
        (:PetscSFBcastAndOpWithMemtypeBegin, $petsc_library),
        PetscErrorCode,
        (
            PetscSF,
            MPI_Datatype,
            PetscMemType,
            Ptr{Cvoid},
            PetscMemType,
            Ptr{Cvoid},
            MPI_Op,
        ),
        sf,
        unit,
        rootmtype,
        rootdata,
        leafmtype,
        leafdata,
        op,
    )
end

@for_petsc function PetscSectionCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscSection}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionClone(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionClone, $petsc_library),
        PetscErrorCode,
        (PetscSection, Ptr{PetscSection}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscSectionSetFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscSection,),
        arg1,
    )
end

@for_petsc function PetscSectionCopy(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionCopy, $petsc_library),
        PetscErrorCode,
        (PetscSection, PetscSection),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionCompare(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSectionCompare, $petsc_library),
        PetscErrorCode,
        (PetscSection, PetscSection, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionGetNumFields(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionGetNumFields, $petsc_library),
        PetscErrorCode,
        (PetscSection, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionSetNumFields(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionSetNumFields, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionGetFieldName(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSectionGetFieldName, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionSetFieldName(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSectionSetFieldName, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionGetComponentName(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSectionGetComponentName, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, $PetscInt, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSectionSetComponentName(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSectionSetComponentName, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, $PetscInt, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSectionGetFieldComponents(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSectionGetFieldComponents, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionSetFieldComponents(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSectionSetFieldComponents, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionGetChart(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSectionGetChart, $petsc_library),
        PetscErrorCode,
        (PetscSection, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionSetChart(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSectionSetChart, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionGetPermutation(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionGetPermutation, $petsc_library),
        PetscErrorCode,
        (PetscSection, Ptr{IS}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionSetPermutation(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionSetPermutation, $petsc_library),
        PetscErrorCode,
        (PetscSection, IS),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionGetPointMajor(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionGetPointMajor, $petsc_library),
        PetscErrorCode,
        (PetscSection, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionSetPointMajor(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionSetPointMajor, $petsc_library),
        PetscErrorCode,
        (PetscSection, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionGetDof(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSectionGetDof, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionSetDof(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSectionSetDof, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionAddDof(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSectionAddDof, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionGetFieldDof(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSectionGetFieldDof, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSectionSetFieldDof(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSectionSetFieldDof, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSectionAddFieldDof(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSectionAddFieldDof, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSectionHasConstraints(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionHasConstraints, $petsc_library),
        PetscErrorCode,
        (PetscSection, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionGetConstraintDof(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSectionGetConstraintDof, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionSetConstraintDof(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSectionSetConstraintDof, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionAddConstraintDof(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSectionAddConstraintDof, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionGetFieldConstraintDof(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSectionGetFieldConstraintDof, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSectionSetFieldConstraintDof(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSectionSetFieldConstraintDof, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSectionAddFieldConstraintDof(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSectionAddFieldConstraintDof, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSectionGetConstraintIndices(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSectionGetConstraintIndices, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionSetConstraintIndices(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSectionSetConstraintIndices, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionGetFieldConstraintIndices(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSectionGetFieldConstraintIndices, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, $PetscInt, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSectionSetFieldConstraintIndices(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSectionSetFieldConstraintIndices, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSectionSetUpBC(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscSectionSetUpBC, $petsc_library),
        PetscErrorCode,
        (PetscSection,),
        arg1,
    )
end

@for_petsc function PetscSectionSetUp(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscSectionSetUp, $petsc_library),
        PetscErrorCode,
        (PetscSection,),
        arg1,
    )
end

@for_petsc function PetscSectionGetMaxDof(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionGetMaxDof, $petsc_library),
        PetscErrorCode,
        (PetscSection, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionGetStorageSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionGetStorageSize, $petsc_library),
        PetscErrorCode,
        (PetscSection, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionGetConstrainedStorageSize(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscSectionGetConstrainedStorageSize, $petsc_library),
        PetscErrorCode,
        (PetscSection, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionGetOffset(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSectionGetOffset, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionSetOffset(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSectionSetOffset, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionGetFieldOffset(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSectionGetFieldOffset, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSectionSetFieldOffset(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSectionSetFieldOffset, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSectionGetFieldPointOffset(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSectionGetFieldPointOffset, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSectionGetOffsetRange(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSectionGetOffsetRange, $petsc_library),
        PetscErrorCode,
        (PetscSection, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionView, $petsc_library),
        PetscErrorCode,
        (PetscSection, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionViewFromOptions(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSectionViewFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscSection, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionReset(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscSectionReset, $petsc_library),
        PetscErrorCode,
        (PetscSection,),
        arg1,
    )
end

@for_petsc function PetscSectionDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscSectionDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscSection},),
        arg1,
    )
end

@for_petsc function PetscSectionCreateGlobalSection(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscSectionCreateGlobalSection, $petsc_library),
        PetscErrorCode,
        (PetscSection, PetscSF, PetscBool, PetscBool, Ptr{PetscSection}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscSectionCreateGlobalSectionCensored(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscSectionCreateGlobalSectionCensored, $petsc_library),
        PetscErrorCode,
        (
            PetscSection,
            PetscSF,
            PetscBool,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{PetscSection},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscSectionCreateSubsection(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSectionCreateSubsection, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, Ptr{$PetscInt}, Ptr{PetscSection}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSectionCreateSupersection(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSectionCreateSupersection, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscSection}, $PetscInt, Ptr{PetscSection}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionCreateSubmeshSection(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSectionCreateSubmeshSection, $petsc_library),
        PetscErrorCode,
        (PetscSection, IS, Ptr{PetscSection}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionPermute(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSectionPermute, $petsc_library),
        PetscErrorCode,
        (PetscSection, IS, Ptr{PetscSection}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionGetField(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSectionGetField, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, Ptr{PetscSection}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionSetUseFieldOffsets(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionSetUseFieldOffsets, $petsc_library),
        PetscErrorCode,
        (PetscSection, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionGetUseFieldOffsets(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionGetUseFieldOffsets, $petsc_library),
        PetscErrorCode,
        (PetscSection, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionExtractDofsFromArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscSectionExtractDofsFromArray, $petsc_library),
        PetscErrorCode,
        (
            PetscSection,
            MPI_Datatype,
            Ptr{Cvoid},
            IS,
            Ptr{PetscSection},
            Ptr{Ptr{Cvoid}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscSectionSetClosureIndex(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSectionSetClosureIndex, $petsc_library),
        PetscErrorCode,
        (PetscSection, PetscObject, PetscSection, IS),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSectionGetClosureIndex(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSectionGetClosureIndex, $petsc_library),
        PetscErrorCode,
        (PetscSection, PetscObject, Ptr{PetscSection}, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSectionSetClosurePermutation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSectionSetClosurePermutation, $petsc_library),
        PetscErrorCode,
        (PetscSection, PetscObject, $PetscInt, IS),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscSectionGetClosurePermutation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscSectionGetClosurePermutation, $petsc_library),
        PetscErrorCode,
        (PetscSection, PetscObject, $PetscInt, $PetscInt, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscSectionGetClosureInversePermutation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscSectionGetClosureInversePermutation, $petsc_library),
        PetscErrorCode,
        (PetscSection, PetscObject, $PetscInt, $PetscInt, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscSectionSymSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionSymSetType, $petsc_library),
        PetscErrorCode,
        (PetscSectionSym, PetscSectionSymType),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionSymGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionSymGetType, $petsc_library),
        PetscErrorCode,
        (PetscSectionSym, Ptr{PetscSectionSymType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionSymRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionSymRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionSymCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionSymCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscSectionSym}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionSymDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscSectionSymDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscSectionSym},),
        arg1,
    )
end

@for_petsc function PetscSectionSymView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionSymView, $petsc_library),
        PetscErrorCode,
        (PetscSectionSym, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionSetSym(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionSetSym, $petsc_library),
        PetscErrorCode,
        (PetscSection, PetscSectionSym),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionGetSym(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionGetSym, $petsc_library),
        PetscErrorCode,
        (PetscSection, Ptr{PetscSectionSym}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionSetFieldSym(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSectionSetFieldSym, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, PetscSectionSym),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionGetFieldSym(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSectionGetFieldSym, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, Ptr{PetscSectionSym}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionGetPointSyms(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscSectionGetPointSyms, $petsc_library),
        PetscErrorCode,
        (
            PetscSection,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{$PetscInt}}},
            Ptr{Ptr{Ptr{$PetscScalar}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscSectionRestorePointSyms(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscSectionRestorePointSyms, $petsc_library),
        PetscErrorCode,
        (
            PetscSection,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{$PetscInt}}},
            Ptr{Ptr{Ptr{$PetscScalar}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscSectionGetFieldPointSyms(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscSectionGetFieldPointSyms, $petsc_library),
        PetscErrorCode,
        (
            PetscSection,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{$PetscInt}}},
            Ptr{Ptr{Ptr{$PetscScalar}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscSectionRestoreFieldPointSyms(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscSectionRestoreFieldPointSyms, $petsc_library),
        PetscErrorCode,
        (
            PetscSection,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{$PetscInt}}},
            Ptr{Ptr{Ptr{$PetscScalar}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

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

@for_petsc function MatGetFactor(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatGetFactor, $petsc_library),
        PetscErrorCode,
        (Mat, MatSolverType, MatFactorType, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatGetFactorAvailable(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatGetFactorAvailable, $petsc_library),
        PetscErrorCode,
        (Mat, MatSolverType, MatFactorType, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatFactorGetUseOrdering(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatFactorGetUseOrdering, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function MatFactorGetSolverType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatFactorGetSolverType, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{MatSolverType}),
        arg1,
        arg2,
    )
end

@for_petsc function MatGetFactorType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetFactorType, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{MatFactorType}),
        arg1,
        arg2,
    )
end

@for_petsc function MatSetFactorType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSetFactorType, $petsc_library),
        PetscErrorCode,
        (Mat, MatFactorType),
        arg1,
        arg2,
    )
end

# typedef PetscErrorCode ( * MatSolverFunction ) ( Mat , MatFactorType , Mat * )
const MatSolverFunction = Ptr{Cvoid}

@for_petsc function MatSolverTypeRegister(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatSolverTypeRegister, $petsc_library),
        PetscErrorCode,
        (MatSolverType, MatType, MatFactorType, MatSolverFunction),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatSolverTypeGet(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatSolverTypeGet, $petsc_library),
        PetscErrorCode,
        (
            MatSolverType,
            MatType,
            MatFactorType,
            Ptr{PetscBool},
            Ptr{PetscBool},
            Ptr{MatSolverFunction},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

const MatSolverPackage = MatSolverType

@for_petsc function MatSolverPackageRegister(
    ::$UnionPetscLib,
    stype,
    mtype,
    ftype,
    f,
)
    @chk ccall(
        (:MatSolverPackageRegister, $petsc_library),
        PetscErrorCode,
        (MatSolverType, MatType, MatFactorType, MatSolverFunction),
        stype,
        mtype,
        ftype,
        f,
    )
end

@for_petsc function MatSolverPackageGet(
    ::$UnionPetscLib,
    stype,
    mtype,
    ftype,
    foundmtype,
    foundstype,
    f,
)
    @chk ccall(
        (:MatSolverPackageGet, $petsc_library),
        PetscErrorCode,
        (
            MatSolverType,
            MatType,
            MatFactorType,
            Ptr{PetscBool},
            Ptr{PetscBool},
            Ptr{MatSolverFunction},
        ),
        stype,
        mtype,
        ftype,
        foundmtype,
        foundstype,
        f,
    )
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

@for_petsc function MatProductCreate(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatProductCreate, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, Mat, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatProductCreateWithMat(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatProductCreateWithMat, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, Mat, Mat),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatProductSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatProductSetType, $petsc_library),
        PetscErrorCode,
        (Mat, MatProductType),
        arg1,
        arg2,
    )
end

@for_petsc function MatProductSetAlgorithm(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatProductSetAlgorithm, $petsc_library),
        PetscErrorCode,
        (Mat, MatProductAlgorithm),
        arg1,
        arg2,
    )
end

@for_petsc function MatProductSetFill(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatProductSetFill, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function MatProductSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatProductSetFromOptions, $petsc_library),
        PetscErrorCode,
        (Mat,),
        arg1,
    )
end

@for_petsc function MatProductSymbolic(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatProductSymbolic, $petsc_library),
        PetscErrorCode,
        (Mat,),
        arg1,
    )
end

@for_petsc function MatProductNumeric(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatProductNumeric, $petsc_library),
        PetscErrorCode,
        (Mat,),
        arg1,
    )
end

@for_petsc function MatProductReplaceMats(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatProductReplaceMats, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, Mat, Mat),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatProductClear(::$UnionPetscLib, arg1)
    @chk ccall((:MatProductClear, $petsc_library), PetscErrorCode, (Mat,), arg1)
end

@for_petsc function MatProductView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatProductView, $petsc_library),
        PetscErrorCode,
        (Mat, PetscViewer),
        arg1,
        arg2,
    )
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

@for_petsc function MatInitializePackage(::$UnionPetscLib)
    @chk ccall((:MatInitializePackage, $petsc_library), PetscErrorCode, ())
end

@for_petsc function MatCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function MatSetSizes(::$UnionPetscLib, arg1, arg2, arg3, arg4, arg5)
    @chk ccall(
        (:MatSetSizes, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSetType, $petsc_library),
        PetscErrorCode,
        (Mat, MatType),
        arg1,
        arg2,
    )
end

@for_petsc function MatGetVecType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetVecType, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{VecType}),
        arg1,
        arg2,
    )
end

@for_petsc function MatSetVecType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSetVecType, $petsc_library),
        PetscErrorCode,
        (Mat, VecType),
        arg1,
        arg2,
    )
end

@for_petsc function MatSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatSetFromOptions, $petsc_library),
        PetscErrorCode,
        (Mat,),
        arg1,
    )
end

@for_petsc function MatViewFromOptions(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatViewFromOptions, $petsc_library),
        PetscErrorCode,
        (Mat, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function MatRegisterRootName(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatRegisterRootName, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatSetOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSetOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function MatAppendOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatAppendOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function MatGetOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatSetErrorIfFailure(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSetErrorIfFailure, $petsc_library),
        PetscErrorCode,
        (Mat, PetscBool),
        arg1,
        arg2,
    )
end

@enum MatStructure::UInt32 begin
    DIFFERENT_NONZERO_PATTERN = 0
    SUBSET_NONZERO_PATTERN = 1
    SAME_NONZERO_PATTERN = 2
    UNKNOWN_NONZERO_PATTERN = 3
end

@for_petsc function MatCreateSeqSELL(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatCreateSeqSELL, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MatCreateSELL(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
)
    @chk ccall(
        (:MatCreateSELL, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Mat},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
    )
end

@for_petsc function MatSeqSELLSetPreallocation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatSeqSELLSetPreallocation, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMPISELLSetPreallocation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatMPISELLSetPreallocation, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatCreateSeqDense(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatCreateSeqDense, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatCreateDense(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:MatCreateDense, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscScalar},
            Ptr{Mat},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function MatCreateSeqAIJ(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatCreateSeqAIJ, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MatCreateAIJ(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
)
    @chk ccall(
        (:MatCreateAIJ, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Mat},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
    )
end

@for_petsc function MatCreateMPIAIJWithArrays(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:MatCreateMPIAIJWithArrays, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscScalar},
            Ptr{Mat},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function MatUpdateMPIAIJWithArrays(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:MatUpdateMPIAIJWithArrays, $petsc_library),
        PetscErrorCode,
        (
            Mat,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscScalar},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function MatCreateMPIAIJWithSplitArrays(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
)
    @chk ccall(
        (:MatCreateMPIAIJWithSplitArrays, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscScalar},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscScalar},
            Ptr{Mat},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
    )
end

@for_petsc function MatCreateMPIAIJWithSeqAIJ(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatCreateMPIAIJWithSeqAIJ, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Mat, Mat, Ptr{$PetscInt}, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatCreateSeqBAIJ(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:MatCreateSeqBAIJ, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Mat},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function MatCreateBAIJ(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
)
    @chk ccall(
        (:MatCreateBAIJ, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Mat},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
    )
end

@for_petsc function MatCreateMPIBAIJWithArrays(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
)
    @chk ccall(
        (:MatCreateMPIBAIJWithArrays, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscScalar},
            Ptr{Mat},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
    )
end

@for_petsc function MatSetPreallocationCOO(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatSetPreallocationCOO, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatSetValuesCOO(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatSetValuesCOO, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscScalar}, InsertMode),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatCreateMPIAdj(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:MatCreateMPIAdj, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{Mat},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function MatCreateSeqSBAIJ(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:MatCreateSeqSBAIJ, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Mat},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function MatCreateSBAIJ(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
)
    @chk ccall(
        (:MatCreateSBAIJ, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Mat},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
    )
end

@for_petsc function MatCreateMPISBAIJWithArrays(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
)
    @chk ccall(
        (:MatCreateMPISBAIJWithArrays, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscScalar},
            Ptr{Mat},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
    )
end

@for_petsc function MatSeqSBAIJSetPreallocationCSR(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatSeqSBAIJSetPreallocationCSR, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatMPISBAIJSetPreallocationCSR(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatMPISBAIJSetPreallocationCSR, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatXAIJSetPreallocation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatXAIJSetPreallocation, $petsc_library),
        PetscErrorCode,
        (
            Mat,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MatCreateShell(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:MatCreateShell, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Cvoid},
            Ptr{Mat},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function MatCreateNormal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCreateNormal, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function MatCreateNormalHermitian(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCreateNormalHermitian, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function MatCreateLRC(::$UnionPetscLib, arg1, arg2, arg3, arg4, arg5)
    @chk ccall(
        (:MatCreateLRC, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, Vec, Mat, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatLRCGetMats(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatLRCGetMats, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}, Ptr{Mat}, Ptr{Vec}, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatCreateIS(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:MatCreateIS, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            ISLocalToGlobalMapping,
            ISLocalToGlobalMapping,
            Ptr{Mat},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function MatCreateSeqAIJCRL(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatCreateSeqAIJCRL, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MatCreateMPIAIJCRL(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:MatCreateMPIAIJCRL, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Mat},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function MatCreateScatter(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatCreateScatter, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, VecScatter, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatScatterSetVecScatter(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatScatterSetVecScatter, $petsc_library),
        PetscErrorCode,
        (Mat, VecScatter),
        arg1,
        arg2,
    )
end

@for_petsc function MatScatterGetVecScatter(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatScatterGetVecScatter, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{VecScatter}),
        arg1,
        arg2,
    )
end

@for_petsc function MatCreateBlockMat(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:MatCreateBlockMat, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Mat},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function MatCompositeAddMat(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCompositeAddMat, $petsc_library),
        PetscErrorCode,
        (Mat, Mat),
        arg1,
        arg2,
    )
end

@for_petsc function MatCompositeMerge(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatCompositeMerge, $petsc_library),
        PetscErrorCode,
        (Mat,),
        arg1,
    )
end

@enum MatCompositeMergeType::UInt32 begin
    MAT_COMPOSITE_MERGE_RIGHT = 0
    MAT_COMPOSITE_MERGE_LEFT = 1
end

@for_petsc function MatCompositeSetMergeType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCompositeSetMergeType, $petsc_library),
        PetscErrorCode,
        (Mat, MatCompositeMergeType),
        arg1,
        arg2,
    )
end

@for_petsc function MatCreateComposite(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatCreateComposite, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, Ptr{Mat}, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@enum MatCompositeType::UInt32 begin
    MAT_COMPOSITE_ADDITIVE = 0
    MAT_COMPOSITE_MULTIPLICATIVE = 1
end

@for_petsc function MatCompositeSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCompositeSetType, $petsc_library),
        PetscErrorCode,
        (Mat, MatCompositeType),
        arg1,
        arg2,
    )
end

@for_petsc function MatCompositeGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCompositeGetType, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{MatCompositeType}),
        arg1,
        arg2,
    )
end

@for_petsc function MatCompositeSetMatStructure(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCompositeSetMatStructure, $petsc_library),
        PetscErrorCode,
        (Mat, MatStructure),
        arg1,
        arg2,
    )
end

@for_petsc function MatCompositeGetMatStructure(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCompositeGetMatStructure, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{MatStructure}),
        arg1,
        arg2,
    )
end

@for_petsc function MatCompositeGetNumberMat(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCompositeGetNumberMat, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function MatCompositeGetMat(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatCompositeGetMat, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatCompositeSetScalings(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCompositeSetScalings, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscScalar}),
        arg1,
        arg2,
    )
end

@for_petsc function MatCreateFFT(::$UnionPetscLib, arg1, arg2, arg3, arg4, arg5)
    @chk ccall(
        (:MatCreateFFT, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, Ptr{$PetscInt}, MatType, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatCreateSeqCUFFT(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatCreateSeqCUFFT, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, Ptr{$PetscInt}, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatCreateTranspose(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCreateTranspose, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function MatTransposeGetMat(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatTransposeGetMat, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function MatCreateHermitianTranspose(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCreateHermitianTranspose, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function MatHermitianTransposeGetMat(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatHermitianTransposeGetMat, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function MatCreateSubMatrixVirtual(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatCreateSubMatrixVirtual, $petsc_library),
        PetscErrorCode,
        (Mat, IS, IS, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatSubMatrixVirtualUpdate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatSubMatrixVirtualUpdate, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, IS, IS),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatCreateLocalRef(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatCreateLocalRef, $petsc_library),
        PetscErrorCode,
        (Mat, IS, IS, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatCreateConstantDiagonal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:MatCreateConstantDiagonal, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscScalar,
            Ptr{Mat},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function MatPythonSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatPythonSetType, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function MatResetPreallocation(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatResetPreallocation, $petsc_library),
        PetscErrorCode,
        (Mat,),
        arg1,
    )
end

@for_petsc function MatSetUp(::$UnionPetscLib, arg1)
    @chk ccall((:MatSetUp, $petsc_library), PetscErrorCode, (Mat,), arg1)
end

@for_petsc function MatDestroy(::$UnionPetscLib, arg1)
    @chk ccall((:MatDestroy, $petsc_library), PetscErrorCode, (Ptr{Mat},), arg1)
end

@for_petsc function MatGetNonzeroState(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetNonzeroState, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{PetscObjectState}),
        arg1,
        arg2,
    )
end

@for_petsc function MatConjugate(::$UnionPetscLib, arg1)
    @chk ccall((:MatConjugate, $petsc_library), PetscErrorCode, (Mat,), arg1)
end

@for_petsc function MatRealPart(::$UnionPetscLib, arg1)
    @chk ccall((:MatRealPart, $petsc_library), PetscErrorCode, (Mat,), arg1)
end

@for_petsc function MatImaginaryPart(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatImaginaryPart, $petsc_library),
        PetscErrorCode,
        (Mat,),
        arg1,
    )
end

@for_petsc function MatGetDiagonalBlock(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetDiagonalBlock, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function MatGetTrace(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetTrace, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscScalar}),
        arg1,
        arg2,
    )
end

@for_petsc function MatInvertBlockDiagonal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatInvertBlockDiagonal, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatInvertVariableBlockDiagonal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatInvertVariableBlockDiagonal, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatInvertBlockDiagonalMat(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatInvertBlockDiagonalMat, $petsc_library),
        PetscErrorCode,
        (Mat, Mat),
        arg1,
        arg2,
    )
end

@for_petsc function MatSetValues(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:MatSetValues, $petsc_library),
        PetscErrorCode,
        (
            Mat,
            $PetscInt,
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscScalar},
            InsertMode,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function MatSetValuesBlocked(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:MatSetValuesBlocked, $petsc_library),
        PetscErrorCode,
        (
            Mat,
            $PetscInt,
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscScalar},
            InsertMode,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function MatSetValuesRow(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatSetValuesRow, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatSetValuesRowLocal(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatSetValuesRowLocal, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatSetValuesBatch(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatSetValuesBatch, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatSetRandom(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSetRandom, $petsc_library),
        PetscErrorCode,
        (Mat, PetscRandom),
        arg1,
        arg2,
    )
end

@for_petsc function MatSetStencil(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatSetStencil, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@enum MatAssemblyType::UInt32 begin
    MAT_FLUSH_ASSEMBLY = 1
    MAT_FINAL_ASSEMBLY = 0
end

@for_petsc function MatAssemblyBegin(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatAssemblyBegin, $petsc_library),
        PetscErrorCode,
        (Mat, MatAssemblyType),
        arg1,
        arg2,
    )
end

@for_petsc function MatAssemblyEnd(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatAssemblyEnd, $petsc_library),
        PetscErrorCode,
        (Mat, MatAssemblyType),
        arg1,
        arg2,
    )
end

@for_petsc function MatAssembled(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatAssembled, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{PetscBool}),
        arg1,
        arg2,
    )
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

@for_petsc function MatSetOption(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatSetOption, $petsc_library),
        PetscErrorCode,
        (Mat, MatOption, PetscBool),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatGetOption(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatGetOption, $petsc_library),
        PetscErrorCode,
        (Mat, MatOption, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatPropagateSymmetryOptions(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatPropagateSymmetryOptions, $petsc_library),
        PetscErrorCode,
        (Mat, Mat),
        arg1,
        arg2,
    )
end

@for_petsc function MatGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetType, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{MatType}),
        arg1,
        arg2,
    )
end

@for_petsc function MatGetValues(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatGetValues, $petsc_library),
        PetscErrorCode,
        (
            Mat,
            $PetscInt,
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscScalar},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MatGetRow(::$UnionPetscLib, arg1, arg2, arg3, arg4, arg5)
    @chk ccall(
        (:MatGetRow, $petsc_library),
        PetscErrorCode,
        (
            Mat,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscScalar}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatRestoreRow(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatRestoreRow, $petsc_library),
        PetscErrorCode,
        (
            Mat,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscScalar}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatGetRowUpperTriangular(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatGetRowUpperTriangular, $petsc_library),
        PetscErrorCode,
        (Mat,),
        arg1,
    )
end

@for_petsc function MatRestoreRowUpperTriangular(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatRestoreRowUpperTriangular, $petsc_library),
        PetscErrorCode,
        (Mat,),
        arg1,
    )
end

@for_petsc function MatGetColumnVector(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatGetColumnVector, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatSeqAIJGetArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSeqAIJGetArray, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatSeqAIJGetArrayRead(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSeqAIJGetArrayRead, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatSeqAIJRestoreArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSeqAIJRestoreArray, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatSeqAIJRestoreArrayRead(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSeqAIJRestoreArrayRead, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatSeqAIJGetMaxRowNonzeros(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSeqAIJGetMaxRowNonzeros, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function MatSeqAIJSetValuesLocalFast(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:MatSeqAIJSetValuesLocalFast, $petsc_library),
        PetscErrorCode,
        (
            Mat,
            $PetscInt,
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscScalar},
            InsertMode,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function MatSeqAIJSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSeqAIJSetType, $petsc_library),
        PetscErrorCode,
        (Mat, MatType),
        arg1,
        arg2,
    )
end

@for_petsc function MatSeqAIJRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSeqAIJRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function MatSeqBAIJGetArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSeqBAIJGetArray, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatSeqBAIJRestoreArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSeqBAIJRestoreArray, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatSeqSBAIJGetArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSeqSBAIJGetArray, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatSeqSBAIJRestoreArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSeqSBAIJRestoreArray, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatDenseGetArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatDenseGetArray, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatDenseRestoreArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatDenseRestoreArray, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatDensePlaceArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatDensePlaceArray, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscScalar}),
        arg1,
        arg2,
    )
end

@for_petsc function MatDenseReplaceArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatDenseReplaceArray, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscScalar}),
        arg1,
        arg2,
    )
end

@for_petsc function MatDenseResetArray(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatDenseResetArray, $petsc_library),
        PetscErrorCode,
        (Mat,),
        arg1,
    )
end

@for_petsc function MatDenseGetArrayRead(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatDenseGetArrayRead, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatDenseRestoreArrayRead(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatDenseRestoreArrayRead, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatDenseGetArrayWrite(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatDenseGetArrayWrite, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatDenseRestoreArrayWrite(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatDenseRestoreArrayWrite, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatGetBlockSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetBlockSize, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function MatSetBlockSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSetBlockSize, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function MatGetBlockSizes(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatGetBlockSizes, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatSetBlockSizes(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatSetBlockSizes, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatSetBlockSizesFromMats(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatSetBlockSizesFromMats, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, Mat),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatSetVariableBlockSizes(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatSetVariableBlockSizes, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatGetVariableBlockSizes(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatGetVariableBlockSizes, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatDenseGetColumn(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatDenseGetColumn, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatDenseRestoreColumn(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatDenseRestoreColumn, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatDenseGetColumnVec(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatDenseGetColumnVec, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatDenseRestoreColumnVec(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatDenseRestoreColumnVec, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatDenseGetColumnVecRead(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatDenseGetColumnVecRead, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatDenseRestoreColumnVecRead(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatDenseRestoreColumnVecRead, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatDenseGetColumnVecWrite(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatDenseGetColumnVecWrite, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatDenseRestoreColumnVecWrite(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatDenseRestoreColumnVecWrite, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatDenseGetSubMatrix(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatDenseGetSubMatrix, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatDenseRestoreSubMatrix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatDenseRestoreSubMatrix, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function MatMult(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMult, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMultDiagonalBlock(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMultDiagonalBlock, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMultAdd(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatMultAdd, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatMultTranspose(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMultTranspose, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMultHermitianTranspose(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatMultHermitianTranspose, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatIsTranspose(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatIsTranspose, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, $PetscReal, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatIsHermitianTranspose(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatIsHermitianTranspose, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, $PetscReal, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatMultTransposeAdd(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatMultTransposeAdd, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatMultHermitianTransposeAdd(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatMultHermitianTransposeAdd, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatMultConstrained(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMultConstrained, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMultTransposeConstrained(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatMultTransposeConstrained, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMatSolve(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMatSolve, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, Mat),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMatSolveTranspose(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMatSolveTranspose, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, Mat),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMatTransposeSolve(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMatTransposeSolve, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, Mat),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatResidual(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatResidual, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@enum MatDuplicateOption::UInt32 begin
    MAT_DO_NOT_COPY_VALUES = 0
    MAT_COPY_VALUES = 1
    MAT_SHARE_NONZERO_PATTERN = 2
end

@for_petsc function MatConvert(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatConvert, $petsc_library),
        PetscErrorCode,
        (Mat, MatType, MatReuse, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatDuplicate(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatDuplicate, $petsc_library),
        PetscErrorCode,
        (Mat, MatDuplicateOption, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatCopy(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatCopy, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, MatStructure),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatView, $petsc_library),
        PetscErrorCode,
        (Mat, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function MatIsSymmetric(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatIsSymmetric, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscReal, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatIsStructurallySymmetric(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatIsStructurallySymmetric, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function MatIsHermitian(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatIsHermitian, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscReal, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatIsSymmetricKnown(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatIsSymmetricKnown, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{PetscBool}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatIsHermitianKnown(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatIsHermitianKnown, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{PetscBool}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMissingDiagonal(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMissingDiagonal, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{PetscBool}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatLoad(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatLoad, $petsc_library),
        PetscErrorCode,
        (Mat, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function MatGetRowIJ(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:MatGetRowIJ, $petsc_library),
        PetscErrorCode,
        (
            Mat,
            $PetscInt,
            PetscBool,
            PetscBool,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscInt}},
            Ptr{PetscBool},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function MatRestoreRowIJ(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:MatRestoreRowIJ, $petsc_library),
        PetscErrorCode,
        (
            Mat,
            $PetscInt,
            PetscBool,
            PetscBool,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscInt}},
            Ptr{PetscBool},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function MatGetColumnIJ(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:MatGetColumnIJ, $petsc_library),
        PetscErrorCode,
        (
            Mat,
            $PetscInt,
            PetscBool,
            PetscBool,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscInt}},
            Ptr{PetscBool},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function MatRestoreColumnIJ(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:MatRestoreColumnIJ, $petsc_library),
        PetscErrorCode,
        (
            Mat,
            $PetscInt,
            PetscBool,
            PetscBool,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscInt}},
            Ptr{PetscBool},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
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

@for_petsc function MatGetInfo(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatGetInfo, $petsc_library),
        PetscErrorCode,
        (Mat, MatInfoType, Ptr{MatInfo}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatGetDiagonal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetDiagonal, $petsc_library),
        PetscErrorCode,
        (Mat, Vec),
        arg1,
        arg2,
    )
end

@for_petsc function MatGetRowMax(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatGetRowMax, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatGetRowMin(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatGetRowMin, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatGetRowMaxAbs(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatGetRowMaxAbs, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatGetRowMinAbs(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatGetRowMinAbs, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatGetRowSum(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetRowSum, $petsc_library),
        PetscErrorCode,
        (Mat, Vec),
        arg1,
        arg2,
    )
end

@for_petsc function MatTranspose(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatTranspose, $petsc_library),
        PetscErrorCode,
        (Mat, MatReuse, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatHermitianTranspose(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatHermitianTranspose, $petsc_library),
        PetscErrorCode,
        (Mat, MatReuse, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatPermute(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatPermute, $petsc_library),
        PetscErrorCode,
        (Mat, IS, IS, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatDiagonalScale(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatDiagonalScale, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatDiagonalSet(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatDiagonalSet, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, InsertMode),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatEqual(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatEqual, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMultEqual(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatMultEqual, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, $PetscInt, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatMultAddEqual(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatMultAddEqual, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, $PetscInt, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatMultTransposeEqual(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatMultTransposeEqual, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, $PetscInt, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatMultTransposeAddEqual(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatMultTransposeAddEqual, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, $PetscInt, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatMatMultEqual(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatMatMultEqual, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, Mat, $PetscInt, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatTransposeMatMultEqual(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatTransposeMatMultEqual, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, Mat, $PetscInt, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatMatTransposeMultEqual(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatMatTransposeMultEqual, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, Mat, $PetscInt, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatPtAPMultEqual(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatPtAPMultEqual, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, Mat, $PetscInt, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatRARtMultEqual(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatRARtMultEqual, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, Mat, $PetscInt, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatIsLinear(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatIsLinear, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatNorm(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatNorm, $petsc_library),
        PetscErrorCode,
        (Mat, NormType, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatGetColumnNorms(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatGetColumnNorms, $petsc_library),
        PetscErrorCode,
        (Mat, NormType, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatZeroEntries(::$UnionPetscLib, arg1)
    @chk ccall((:MatZeroEntries, $petsc_library), PetscErrorCode, (Mat,), arg1)
end

@for_petsc function MatSetInf(::$UnionPetscLib, arg1)
    @chk ccall((:MatSetInf, $petsc_library), PetscErrorCode, (Mat,), arg1)
end

@for_petsc function MatZeroRows(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatZeroRows, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}, $PetscScalar, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MatZeroRowsIS(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatZeroRowsIS, $petsc_library),
        PetscErrorCode,
        (Mat, IS, $PetscScalar, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatZeroRowsColumns(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatZeroRowsColumns, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}, $PetscScalar, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MatZeroRowsColumnsIS(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatZeroRowsColumnsIS, $petsc_library),
        PetscErrorCode,
        (Mat, IS, $PetscScalar, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatGetSize(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatGetSize, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatGetLocalSize(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatGetLocalSize, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatGetOwnershipRange(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatGetOwnershipRange, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatGetOwnershipRanges(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetOwnershipRanges, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatGetOwnershipRangeColumn(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatGetOwnershipRangeColumn, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatGetOwnershipRangesColumn(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetOwnershipRangesColumn, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatGetOwnershipIS(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatGetOwnershipIS, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{IS}, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatCreateSubMatrices(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatCreateSubMatrices, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{IS}, Ptr{IS}, MatReuse, Ptr{Ptr{Mat}}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MatGetSubMatrices(
    ::$UnionPetscLib,
    mat,
    n,
    irow,
    icol,
    scall,
    submat,
)
    @chk ccall(
        (:MatGetSubMatrices, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{IS}, Ptr{IS}, MatReuse, Ptr{Ptr{Mat}}),
        mat,
        n,
        irow,
        icol,
        scall,
        submat,
    )
end

@for_petsc function MatCreateSubMatricesMPI(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatCreateSubMatricesMPI, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{IS}, Ptr{IS}, MatReuse, Ptr{Ptr{Mat}}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MatGetSubMatricesMPI(
    ::$UnionPetscLib,
    mat,
    n,
    irow,
    icol,
    scall,
    submat,
)
    @chk ccall(
        (:MatGetSubMatricesMPI, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{IS}, Ptr{IS}, MatReuse, Ptr{Ptr{Mat}}),
        mat,
        n,
        irow,
        icol,
        scall,
        submat,
    )
end

@for_petsc function MatDestroyMatrices(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatDestroyMatrices, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{Ptr{Mat}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatDestroySubMatrices(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatDestroySubMatrices, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{Ptr{Mat}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatCreateSubMatrix(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatCreateSubMatrix, $petsc_library),
        PetscErrorCode,
        (Mat, IS, IS, MatReuse, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatGetSubMatrix(
    ::$UnionPetscLib,
    mat,
    isrow,
    iscol,
    cll,
    newmat,
)
    @chk ccall(
        (:MatGetSubMatrix, $petsc_library),
        PetscErrorCode,
        (Mat, IS, IS, MatReuse, Ptr{Mat}),
        mat,
        isrow,
        iscol,
        cll,
        newmat,
    )
end

@for_petsc function MatGetLocalSubMatrix(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatGetLocalSubMatrix, $petsc_library),
        PetscErrorCode,
        (Mat, IS, IS, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatRestoreLocalSubMatrix(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatRestoreLocalSubMatrix, $petsc_library),
        PetscErrorCode,
        (Mat, IS, IS, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatGetSeqNonzeroStructure(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetSeqNonzeroStructure, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function MatDestroySeqNonzeroStructure(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatDestroySeqNonzeroStructure, $petsc_library),
        PetscErrorCode,
        (Ptr{Mat},),
        arg1,
    )
end

@for_petsc function MatCreateMPIAIJSumSeqAIJ(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatCreateMPIAIJSumSeqAIJ, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Mat, $PetscInt, $PetscInt, MatReuse, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MatCreateMPIAIJSumSeqAIJSymbolic(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatCreateMPIAIJSumSeqAIJSymbolic, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Mat, $PetscInt, $PetscInt, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatCreateMPIAIJSumSeqAIJNumeric(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatCreateMPIAIJSumSeqAIJNumeric, $petsc_library),
        PetscErrorCode,
        (Mat, Mat),
        arg1,
        arg2,
    )
end

@for_petsc function MatMPIAIJGetLocalMat(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMPIAIJGetLocalMat, $petsc_library),
        PetscErrorCode,
        (Mat, MatReuse, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMPIAIJGetLocalMatCondensed(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatMPIAIJGetLocalMatCondensed, $petsc_library),
        PetscErrorCode,
        (Mat, MatReuse, Ptr{IS}, Ptr{IS}, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatMPIAIJGetLocalMatMerge(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatMPIAIJGetLocalMatMerge, $petsc_library),
        PetscErrorCode,
        (Mat, MatReuse, Ptr{IS}, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatGetBrowsOfAcols(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatGetBrowsOfAcols, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, MatReuse, Ptr{IS}, Ptr{IS}, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MatGetGhosts(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatGetGhosts, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatIncreaseOverlap(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatIncreaseOverlap, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{IS}, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatIncreaseOverlapSplit(::$UnionPetscLib, mat, n, is, ov)
    @chk ccall(
        (:MatIncreaseOverlapSplit, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{IS}, $PetscInt),
        mat,
        n,
        is,
        ov,
    )
end

@for_petsc function MatMPIAIJSetUseScalableIncreaseOverlap(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatMPIAIJSetUseScalableIncreaseOverlap, $petsc_library),
        PetscErrorCode,
        (Mat, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function MatMatMult(::$UnionPetscLib, arg1, arg2, arg3, arg4, arg5)
    @chk ccall(
        (:MatMatMult, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, MatReuse, $PetscReal, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatMatMatMult(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatMatMatMult, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, Mat, MatReuse, $PetscReal, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MatGalerkin(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatGalerkin, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, Mat, MatReuse, $PetscReal, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MatPtAP(::$UnionPetscLib, arg1, arg2, arg3, arg4, arg5)
    @chk ccall(
        (:MatPtAP, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, MatReuse, $PetscReal, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatRARt(::$UnionPetscLib, arg1, arg2, arg3, arg4, arg5)
    @chk ccall(
        (:MatRARt, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, MatReuse, $PetscReal, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatTransposeMatMult(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatTransposeMatMult, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, MatReuse, $PetscReal, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatMatTransposeMult(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatMatTransposeMult, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, MatReuse, $PetscReal, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatAXPY(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatAXPY, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscScalar, Mat, MatStructure),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatAYPX(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatAYPX, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscScalar, Mat, MatStructure),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatScale(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatScale, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscScalar),
        arg1,
        arg2,
    )
end

@for_petsc function MatShift(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatShift, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscScalar),
        arg1,
        arg2,
    )
end

@for_petsc function MatSetLocalToGlobalMapping(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatSetLocalToGlobalMapping, $petsc_library),
        PetscErrorCode,
        (Mat, ISLocalToGlobalMapping, ISLocalToGlobalMapping),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatGetLocalToGlobalMapping(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatGetLocalToGlobalMapping, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{ISLocalToGlobalMapping}, Ptr{ISLocalToGlobalMapping}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatZeroRowsLocal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatZeroRowsLocal, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}, $PetscScalar, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MatZeroRowsLocalIS(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatZeroRowsLocalIS, $petsc_library),
        PetscErrorCode,
        (Mat, IS, $PetscScalar, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatZeroRowsColumnsLocal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatZeroRowsColumnsLocal, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}, $PetscScalar, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MatZeroRowsColumnsLocalIS(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatZeroRowsColumnsLocalIS, $petsc_library),
        PetscErrorCode,
        (Mat, IS, $PetscScalar, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatGetValuesLocal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatGetValuesLocal, $petsc_library),
        PetscErrorCode,
        (
            Mat,
            $PetscInt,
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscScalar},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MatSetValuesLocal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:MatSetValuesLocal, $petsc_library),
        PetscErrorCode,
        (
            Mat,
            $PetscInt,
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscScalar},
            InsertMode,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function MatSetValuesBlockedLocal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:MatSetValuesBlockedLocal, $petsc_library),
        PetscErrorCode,
        (
            Mat,
            $PetscInt,
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscScalar},
            InsertMode,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function MatStashSetInitialSize(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatStashSetInitialSize, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatStashGetInfo(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatStashGetInfo, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatInterpolate(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatInterpolate, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatInterpolateAdd(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatInterpolateAdd, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatRestrict(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatRestrict, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMatInterpolate(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMatInterpolate, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMatInterpolateAdd(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatMatInterpolateAdd, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, Mat, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatMatRestrict(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMatRestrict, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatCreateVecs(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatCreateVecs, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Vec}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatGetVecs(::$UnionPetscLib, mat, x, y)
    @chk ccall(
        (:MatGetVecs, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Vec}, Ptr{Vec}),
        mat,
        x,
        y,
    )
end

@for_petsc function MatCreateRedundantMatrix(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatCreateRedundantMatrix, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, MPI_Comm, MatReuse, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatGetMultiProcBlock(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatGetMultiProcBlock, $petsc_library),
        PetscErrorCode,
        (Mat, MPI_Comm, MatReuse, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatFindZeroDiagonals(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatFindZeroDiagonals, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{IS}),
        arg1,
        arg2,
    )
end

@for_petsc function MatFindOffBlockDiagonalEntries(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatFindOffBlockDiagonalEntries, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{IS}),
        arg1,
        arg2,
    )
end

@for_petsc function MatCreateMPIMatConcatenateSeqMat(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatCreateMPIMatConcatenateSeqMat, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Mat, $PetscInt, MatReuse, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatSetValue(::$UnionPetscLib, v, i, j, va, mode)
    @chk ccall(
        (:MatSetValue, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt, $PetscScalar, InsertMode),
        v,
        i,
        j,
        va,
        mode,
    )
end

@for_petsc function MatGetValue(::$UnionPetscLib, v, i, j, va)
    @chk ccall(
        (:MatGetValue, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt, Ptr{$PetscScalar}),
        v,
        i,
        j,
        va,
    )
end

@for_petsc function MatSetValueLocal(::$UnionPetscLib, v, i, j, va, mode)
    @chk ccall(
        (:MatSetValueLocal, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt, $PetscScalar, InsertMode),
        v,
        i,
        j,
        va,
        mode,
    )
end

@for_petsc function MatShellGetContext(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatShellGetContext, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function MatInodeAdjustForInodes(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatInodeAdjustForInodes, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{IS}, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatInodeGetInodeSizes(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatInodeGetInodeSizes, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatSeqAIJSetColumnIndices(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSeqAIJSetColumnIndices, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function MatSeqBAIJSetColumnIndices(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSeqBAIJSetColumnIndices, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function MatCreateSeqAIJWithArrays(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:MatCreateSeqAIJWithArrays, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscScalar},
            Ptr{Mat},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function MatCreateSeqBAIJWithArrays(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:MatCreateSeqBAIJWithArrays, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscScalar},
            Ptr{Mat},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function MatCreateSeqSBAIJWithArrays(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:MatCreateSeqSBAIJWithArrays, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscScalar},
            Ptr{Mat},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function MatCreateSeqAIJFromTriple(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:MatCreateSeqAIJFromTriple, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscScalar},
            Ptr{Mat},
            $PetscInt,
            PetscBool,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function MatSeqBAIJSetPreallocation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatSeqBAIJSetPreallocation, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatSeqSBAIJSetPreallocation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatSeqSBAIJSetPreallocation, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatSeqAIJSetPreallocation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatSeqAIJSetPreallocation, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatSeqAIJSetTotalPreallocation(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSeqAIJSetTotalPreallocation, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function MatMPIBAIJSetPreallocation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatMPIBAIJSetPreallocation, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MatMPISBAIJSetPreallocation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatMPISBAIJSetPreallocation, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MatMPIAIJSetPreallocation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatMPIAIJSetPreallocation, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatSeqAIJSetPreallocationCSR(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatSeqAIJSetPreallocationCSR, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatSeqBAIJSetPreallocationCSR(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatSeqBAIJSetPreallocationCSR, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatMPIAIJSetPreallocationCSR(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatMPIAIJSetPreallocationCSR, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatMPIBAIJSetPreallocationCSR(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatMPIBAIJSetPreallocationCSR, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatMPIAdjSetPreallocation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatMPIAdjSetPreallocation, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatMPIAdjToSeq(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatMPIAdjToSeq, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function MatMPIDenseSetPreallocation(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatMPIDenseSetPreallocation, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscScalar}),
        arg1,
        arg2,
    )
end

@for_petsc function MatSeqDenseSetPreallocation(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSeqDenseSetPreallocation, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscScalar}),
        arg1,
        arg2,
    )
end

@for_petsc function MatMPIAIJGetSeqAIJ(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatMPIAIJGetSeqAIJ, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}, Ptr{Mat}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatMPIBAIJGetSeqBAIJ(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatMPIBAIJGetSeqBAIJ, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}, Ptr{Mat}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatMPIAdjCreateNonemptySubcommMat(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatMPIAdjCreateNonemptySubcommMat, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function MatDenseGetLDA(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatDenseGetLDA, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function MatDenseSetLDA(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatDenseSetLDA, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function MatSeqDenseSetLDA(::$UnionPetscLib, A, lda)
    @chk ccall(
        (:MatSeqDenseSetLDA, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt),
        A,
        lda,
    )
end

@for_petsc function MatDenseGetLocalMatrix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatDenseGetLocalMatrix, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function MatBlockMatSetPreallocation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatBlockMatSetPreallocation, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatStoreValues(::$UnionPetscLib, arg1)
    @chk ccall((:MatStoreValues, $petsc_library), PetscErrorCode, (Mat,), arg1)
end

@for_petsc function MatRetrieveValues(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatRetrieveValues, $petsc_library),
        PetscErrorCode,
        (Mat,),
        arg1,
    )
end

@for_petsc function MatFindNonzeroRows(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatFindNonzeroRows, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{IS}),
        arg1,
        arg2,
    )
end

@for_petsc function MatFindZeroRows(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatFindZeroRows, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{IS}),
        arg1,
        arg2,
    )
end

const MatOrderingType = Ptr{Cchar}

@for_petsc function MatGetOrdering(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatGetOrdering, $petsc_library),
        PetscErrorCode,
        (Mat, MatOrderingType, Ptr{IS}, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatGetOrderingList(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatGetOrderingList, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscFunctionList},),
        arg1,
    )
end

@for_petsc function MatOrderingRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatOrderingRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function MatReorderForNonzeroDiagonal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatReorderForNonzeroDiagonal, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscReal, IS, IS),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatCreateLaplacian(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatCreateLaplacian, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscReal, PetscBool, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
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

@for_petsc function MatFactorGetError(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatFactorGetError, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{MatFactorError}),
        arg1,
        arg2,
    )
end

@for_petsc function MatFactorClearError(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatFactorClearError, $petsc_library),
        PetscErrorCode,
        (Mat,),
        arg1,
    )
end

@for_petsc function MatFactorGetErrorZeroPivot(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatFactorGetErrorZeroPivot, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscReal}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatGetInertia(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatGetInertia, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatSolve(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatSolve, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatForwardSolve(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatForwardSolve, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatBackwardSolve(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatBackwardSolve, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatSolveAdd(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatSolveAdd, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatSolveTranspose(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatSolveTranspose, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatSolveTransposeAdd(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatSolveTransposeAdd, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatSetUnfactored(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatSetUnfactored, $petsc_library),
        PetscErrorCode,
        (Mat,),
        arg1,
    )
end

@enum MatFactorSchurStatus::UInt32 begin
    MAT_FACTOR_SCHUR_UNFACTORED = 0
    MAT_FACTOR_SCHUR_FACTORED = 1
    MAT_FACTOR_SCHUR_INVERTED = 2
end

@for_petsc function MatFactorSetSchurIS(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatFactorSetSchurIS, $petsc_library),
        PetscErrorCode,
        (Mat, IS),
        arg1,
        arg2,
    )
end

@for_petsc function MatFactorGetSchurComplement(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatFactorGetSchurComplement, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}, Ptr{MatFactorSchurStatus}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatFactorRestoreSchurComplement(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatFactorRestoreSchurComplement, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}, MatFactorSchurStatus),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatFactorInvertSchurComplement(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatFactorInvertSchurComplement, $petsc_library),
        PetscErrorCode,
        (Mat,),
        arg1,
    )
end

@for_petsc function MatFactorCreateSchurComplement(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatFactorCreateSchurComplement, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}, Ptr{MatFactorSchurStatus}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatFactorSolveSchurComplement(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatFactorSolveSchurComplement, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatFactorSolveSchurComplementTranspose(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatFactorSolveSchurComplementTranspose, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatFactorFactorizeSchurComplement(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatFactorFactorizeSchurComplement, $petsc_library),
        PetscErrorCode,
        (Mat,),
        arg1,
    )
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

@for_petsc function MatSOR(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:MatSOR, $petsc_library),
        PetscErrorCode,
        (
            Mat,
            Vec,
            $PetscReal,
            MatSORType,
            $PetscReal,
            $PetscInt,
            $PetscInt,
            Vec,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
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

@for_petsc function MatColoringCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatColoringCreate, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{MatColoring}),
        arg1,
        arg2,
    )
end

@for_petsc function MatColoringGetDegrees(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatColoringGetDegrees, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatColoringDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatColoringDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{MatColoring},),
        arg1,
    )
end

@for_petsc function MatColoringView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatColoringView, $petsc_library),
        PetscErrorCode,
        (MatColoring, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function MatColoringSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatColoringSetType, $petsc_library),
        PetscErrorCode,
        (MatColoring, MatColoringType),
        arg1,
        arg2,
    )
end

@for_petsc function MatColoringSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatColoringSetFromOptions, $petsc_library),
        PetscErrorCode,
        (MatColoring,),
        arg1,
    )
end

@for_petsc function MatColoringSetDistance(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatColoringSetDistance, $petsc_library),
        PetscErrorCode,
        (MatColoring, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function MatColoringGetDistance(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatColoringGetDistance, $petsc_library),
        PetscErrorCode,
        (MatColoring, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function MatColoringSetMaxColors(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatColoringSetMaxColors, $petsc_library),
        PetscErrorCode,
        (MatColoring, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function MatColoringGetMaxColors(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatColoringGetMaxColors, $petsc_library),
        PetscErrorCode,
        (MatColoring, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function MatColoringApply(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatColoringApply, $petsc_library),
        PetscErrorCode,
        (MatColoring, Ptr{ISColoring}),
        arg1,
        arg2,
    )
end

@for_petsc function MatColoringRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatColoringRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function MatColoringPatch(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatColoringPatch, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt, Ptr{ISColoringValue}, Ptr{ISColoring}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatColoringSetWeightType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatColoringSetWeightType, $petsc_library),
        PetscErrorCode,
        (MatColoring, MatColoringWeightType),
        arg1,
        arg2,
    )
end

@for_petsc function MatColoringSetWeights(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatColoringSetWeights, $petsc_library),
        PetscErrorCode,
        (MatColoring, Ptr{$PetscReal}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatColoringCreateWeights(
    ::$UnionPetscLib,
    arg1,
    arg2,
    lperm,
)
    @chk ccall(
        (:MatColoringCreateWeights, $petsc_library),
        PetscErrorCode,
        (MatColoring, Ptr{Ptr{$PetscReal}}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        lperm,
    )
end

@for_petsc function MatColoringTest(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatColoringTest, $petsc_library),
        PetscErrorCode,
        (MatColoring, ISColoring),
        arg1,
        arg2,
    )
end

@for_petsc function MatColoringTestValid(
    ::$UnionPetscLib,
    matcoloring,
    iscoloring,
)
    @chk ccall(
        (:MatColoringTestValid, $petsc_library),
        PetscErrorCode,
        (MatColoring, ISColoring),
        matcoloring,
        iscoloring,
    )
end

@for_petsc function MatISColoringTest(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatISColoringTest, $petsc_library),
        PetscErrorCode,
        (Mat, ISColoring),
        arg1,
        arg2,
    )
end

mutable struct _p_MatFDColoring end

const MatFDColoring = Ptr{_p_MatFDColoring}

@for_petsc function MatFDColoringCreate(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatFDColoringCreate, $petsc_library),
        PetscErrorCode,
        (Mat, ISColoring, Ptr{MatFDColoring}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatFDColoringDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatFDColoringDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{MatFDColoring},),
        arg1,
    )
end

@for_petsc function MatFDColoringView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatFDColoringView, $petsc_library),
        PetscErrorCode,
        (MatFDColoring, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function MatFDColoringSetFunction(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatFDColoringSetFunction, $petsc_library),
        PetscErrorCode,
        (MatFDColoring, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatFDColoringGetFunction(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatFDColoringGetFunction, $petsc_library),
        PetscErrorCode,
        (MatFDColoring, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatFDColoringSetParameters(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatFDColoringSetParameters, $petsc_library),
        PetscErrorCode,
        (MatFDColoring, $PetscReal, $PetscReal),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatFDColoringSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatFDColoringSetFromOptions, $petsc_library),
        PetscErrorCode,
        (MatFDColoring,),
        arg1,
    )
end

@for_petsc function MatFDColoringApply(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatFDColoringApply, $petsc_library),
        PetscErrorCode,
        (Mat, MatFDColoring, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatFDColoringSetF(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatFDColoringSetF, $petsc_library),
        PetscErrorCode,
        (MatFDColoring, Vec),
        arg1,
        arg2,
    )
end

@for_petsc function MatFDColoringGetPerturbedColumns(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatFDColoringGetPerturbedColumns, $petsc_library),
        PetscErrorCode,
        (MatFDColoring, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatFDColoringSetUp(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatFDColoringSetUp, $petsc_library),
        PetscErrorCode,
        (Mat, ISColoring, MatFDColoring),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatFDColoringSetBlockSize(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatFDColoringSetBlockSize, $petsc_library),
        PetscErrorCode,
        (MatFDColoring, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatFDColoringSetValues(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatFDColoringSetValues, $petsc_library),
        PetscErrorCode,
        (Mat, MatFDColoring, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
    )
end

mutable struct _p_MatTransposeColoring end

const MatTransposeColoring = Ptr{_p_MatTransposeColoring}

@for_petsc function MatTransposeColoringCreate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatTransposeColoringCreate, $petsc_library),
        PetscErrorCode,
        (Mat, ISColoring, Ptr{MatTransposeColoring}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatTransColoringApplySpToDen(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatTransColoringApplySpToDen, $petsc_library),
        PetscErrorCode,
        (MatTransposeColoring, Mat, Mat),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatTransColoringApplyDenToSp(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatTransColoringApplyDenToSp, $petsc_library),
        PetscErrorCode,
        (MatTransposeColoring, Mat, Mat),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatTransposeColoringDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatTransposeColoringDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{MatTransposeColoring},),
        arg1,
    )
end

mutable struct _p_MatPartitioning end

const MatPartitioning = Ptr{_p_MatPartitioning}

const MatPartitioningType = Ptr{Cchar}

@for_petsc function MatPartitioningCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatPartitioningCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{MatPartitioning}),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatPartitioningSetType, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, MatPartitioningType),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningSetNParts(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatPartitioningSetNParts, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningSetAdjacency(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatPartitioningSetAdjacency, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, Mat),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningSetVertexWeights(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatPartitioningSetVertexWeights, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningSetPartitionWeights(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatPartitioningSetPartitionWeights, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningSetUseEdgeWeights(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatPartitioningSetUseEdgeWeights, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningGetUseEdgeWeights(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatPartitioningGetUseEdgeWeights, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningApply(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatPartitioningApply, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, Ptr{IS}),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningImprove(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatPartitioningImprove, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, Ptr{IS}),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningViewImbalance(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatPartitioningViewImbalance, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, IS),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningApplyND(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatPartitioningApplyND, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, Ptr{IS}),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatPartitioningDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{MatPartitioning},),
        arg1,
    )
end

@for_petsc function MatPartitioningRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatPartitioningRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatPartitioningView, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningViewFromOptions(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatPartitioningViewFromOptions, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatPartitioningSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatPartitioningSetFromOptions, $petsc_library),
        PetscErrorCode,
        (MatPartitioning,),
        arg1,
    )
end

@for_petsc function MatPartitioningGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatPartitioningGetType, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, Ptr{MatPartitioningType}),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningParmetisSetRepartition(
    ::$UnionPetscLib,
    part,
)
    @chk ccall(
        (:MatPartitioningParmetisSetRepartition, $petsc_library),
        PetscErrorCode,
        (MatPartitioning,),
        part,
    )
end

@for_petsc function MatPartitioningParmetisSetCoarseSequential(
    ::$UnionPetscLib,
    arg1,
)
    @chk ccall(
        (:MatPartitioningParmetisSetCoarseSequential, $petsc_library),
        PetscErrorCode,
        (MatPartitioning,),
        arg1,
    )
end

@for_petsc function MatPartitioningParmetisGetEdgeCut(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatPartitioningParmetisGetEdgeCut, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
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

@for_petsc function MatPartitioningChacoSetGlobal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatPartitioningChacoSetGlobal, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, MPChacoGlobalType),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningChacoGetGlobal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatPartitioningChacoGetGlobal, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, Ptr{MPChacoGlobalType}),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningChacoSetLocal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatPartitioningChacoSetLocal, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, MPChacoLocalType),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningChacoGetLocal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatPartitioningChacoGetLocal, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, Ptr{MPChacoLocalType}),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningChacoSetCoarseLevel(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatPartitioningChacoSetCoarseLevel, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningChacoSetEigenSolver(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatPartitioningChacoSetEigenSolver, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, MPChacoEigenType),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningChacoGetEigenSolver(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatPartitioningChacoGetEigenSolver, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, Ptr{MPChacoEigenType}),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningChacoSetEigenTol(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatPartitioningChacoSetEigenTol, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningChacoGetEigenTol(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatPartitioningChacoGetEigenTol, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningChacoSetEigenNumber(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatPartitioningChacoSetEigenNumber, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningChacoGetEigenNumber(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatPartitioningChacoGetEigenNumber, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningPartySetGlobal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatPartitioningPartySetGlobal, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningPartySetLocal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatPartitioningPartySetLocal, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningPartySetCoarseLevel(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatPartitioningPartySetCoarseLevel, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningPartySetBipart(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatPartitioningPartySetBipart, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningPartySetMatchOptimization(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatPartitioningPartySetMatchOptimization, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, PetscBool),
        arg1,
        arg2,
    )
end

@enum MPPTScotchStrategyType::UInt32 begin
    MP_PTSCOTCH_DEFAULT = 0
    MP_PTSCOTCH_QUALITY = 1
    MP_PTSCOTCH_SPEED = 2
    MP_PTSCOTCH_BALANCE = 3
    MP_PTSCOTCH_SAFETY = 4
    MP_PTSCOTCH_SCALABILITY = 5
end

@for_petsc function MatPartitioningPTScotchSetImbalance(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatPartitioningPTScotchSetImbalance, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningPTScotchGetImbalance(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatPartitioningPTScotchGetImbalance, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningPTScotchSetStrategy(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatPartitioningPTScotchSetStrategy, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, MPPTScotchStrategyType),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningPTScotchGetStrategy(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatPartitioningPTScotchGetStrategy, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, Ptr{MPPTScotchStrategyType}),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningHierarchicalGetFineparts(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatPartitioningHierarchicalGetFineparts, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, Ptr{IS}),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningHierarchicalGetCoarseparts(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatPartitioningHierarchicalGetCoarseparts, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, Ptr{IS}),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningHierarchicalSetNcoarseparts(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatPartitioningHierarchicalSetNcoarseparts, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function MatPartitioningHierarchicalSetNfineparts(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatPartitioningHierarchicalSetNfineparts, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function MatMeshToVertexGraph(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMeshToVertexGraph, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMeshToCellGraph(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMeshToCellGraph, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
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

@for_petsc function MatSetOperation(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatSetOperation, $petsc_library),
        PetscErrorCode,
        (Mat, MatOperation, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatGetOperation(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatGetOperation, $petsc_library),
        PetscErrorCode,
        (Mat, MatOperation, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatHasOperation(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatHasOperation, $petsc_library),
        PetscErrorCode,
        (Mat, MatOperation, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatHasCongruentLayouts(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatHasCongruentLayouts, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function MatFreeIntermediateDataStructures(::$UnionPetscLib, A)
    @chk ccall(
        (:MatFreeIntermediateDataStructures, $petsc_library),
        PetscErrorCode,
        (Mat,),
        A,
    )
end

@for_petsc function MatShellSetOperation(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatShellSetOperation, $petsc_library),
        PetscErrorCode,
        (Mat, MatOperation, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatShellGetOperation(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatShellGetOperation, $petsc_library),
        PetscErrorCode,
        (Mat, MatOperation, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatShellSetContext(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatShellSetContext, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function MatShellSetVecType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatShellSetVecType, $petsc_library),
        PetscErrorCode,
        (Mat, VecType),
        arg1,
        arg2,
    )
end

@for_petsc function MatShellTestMult(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatShellTestMult, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Cvoid}, Vec, Ptr{Cvoid}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatShellTestMultTranspose(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatShellTestMultTranspose, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Cvoid}, Vec, Ptr{Cvoid}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatShellSetManageScalingShifts(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatShellSetManageScalingShifts, $petsc_library),
        PetscErrorCode,
        (Mat,),
        arg1,
    )
end

@for_petsc function MatShellSetMatProductOperation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:MatShellSetMatProductOperation, $petsc_library),
        PetscErrorCode,
        (
            Mat,
            MatProductType,
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
            MatType,
            MatType,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function MatIsShell(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatIsShell, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function MatMPIBAIJSetHashTableFactor(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatMPIBAIJSetHashTableFactor, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function MatISSetLocalMatType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatISSetLocalMatType, $petsc_library),
        PetscErrorCode,
        (Mat, MatType),
        arg1,
        arg2,
    )
end

@for_petsc function MatISSetPreallocation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatISSetPreallocation, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatISStoreL2L(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatISStoreL2L, $petsc_library),
        PetscErrorCode,
        (Mat, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function MatISFixLocalEmpty(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatISFixLocalEmpty, $petsc_library),
        PetscErrorCode,
        (Mat, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function MatISGetLocalMat(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatISGetLocalMat, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function MatISRestoreLocalMat(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatISRestoreLocalMat, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function MatISSetLocalMat(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatISSetLocalMat, $petsc_library),
        PetscErrorCode,
        (Mat, Mat),
        arg1,
        arg2,
    )
end

@for_petsc function MatISGetMPIXAIJ(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatISGetMPIXAIJ, $petsc_library),
        PetscErrorCode,
        (Mat, MatReuse, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

mutable struct _p_MatNullSpace end

const MatNullSpace = Ptr{_p_MatNullSpace}

@for_petsc function MatNullSpaceCreate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatNullSpaceCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, PetscBool, $PetscInt, Ptr{Vec}, Ptr{MatNullSpace}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatNullSpaceSetFunction(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatNullSpaceSetFunction, $petsc_library),
        PetscErrorCode,
        (MatNullSpace, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatNullSpaceDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatNullSpaceDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{MatNullSpace},),
        arg1,
    )
end

@for_petsc function MatNullSpaceRemove(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatNullSpaceRemove, $petsc_library),
        PetscErrorCode,
        (MatNullSpace, Vec),
        arg1,
        arg2,
    )
end

@for_petsc function MatGetNullSpace(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetNullSpace, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{MatNullSpace}),
        arg1,
        arg2,
    )
end

@for_petsc function MatGetTransposeNullSpace(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetTransposeNullSpace, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{MatNullSpace}),
        arg1,
        arg2,
    )
end

@for_petsc function MatSetTransposeNullSpace(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSetTransposeNullSpace, $petsc_library),
        PetscErrorCode,
        (Mat, MatNullSpace),
        arg1,
        arg2,
    )
end

@for_petsc function MatSetNullSpace(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSetNullSpace, $petsc_library),
        PetscErrorCode,
        (Mat, MatNullSpace),
        arg1,
        arg2,
    )
end

@for_petsc function MatSetNearNullSpace(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSetNearNullSpace, $petsc_library),
        PetscErrorCode,
        (Mat, MatNullSpace),
        arg1,
        arg2,
    )
end

@for_petsc function MatGetNearNullSpace(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetNearNullSpace, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{MatNullSpace}),
        arg1,
        arg2,
    )
end

@for_petsc function MatNullSpaceTest(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatNullSpaceTest, $petsc_library),
        PetscErrorCode,
        (MatNullSpace, Mat, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatNullSpaceView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatNullSpaceView, $petsc_library),
        PetscErrorCode,
        (MatNullSpace, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function MatNullSpaceGetVecs(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatNullSpaceGetVecs, $petsc_library),
        PetscErrorCode,
        (MatNullSpace, Ptr{PetscBool}, Ptr{$PetscInt}, Ptr{Ptr{Vec}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatNullSpaceCreateRigidBody(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatNullSpaceCreateRigidBody, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{MatNullSpace}),
        arg1,
        arg2,
    )
end

@for_petsc function MatReorderingSeqSBAIJ(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatReorderingSeqSBAIJ, $petsc_library),
        PetscErrorCode,
        (Mat, IS),
        arg1,
        arg2,
    )
end

@for_petsc function MatMPISBAIJSetHashTableFactor(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatMPISBAIJSetHashTableFactor, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function MatSeqSBAIJSetColumnIndices(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSeqSBAIJSetColumnIndices, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function MatCreateMAIJ(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatCreateMAIJ, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMAIJRedimension(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMAIJRedimension, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMAIJGetAIJ(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatMAIJGetAIJ, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function MatComputeOperator(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatComputeOperator, $petsc_library),
        PetscErrorCode,
        (Mat, MatType, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatComputeOperatorTranspose(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatComputeOperatorTranspose, $petsc_library),
        PetscErrorCode,
        (Mat, MatType, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatComputeExplicitOperator(::$UnionPetscLib, A, B)
    @chk ccall(
        (:MatComputeExplicitOperator, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}),
        A,
        B,
    )
end

@for_petsc function MatComputeExplicitOperatorTranspose(::$UnionPetscLib, A, B)
    @chk ccall(
        (:MatComputeExplicitOperatorTranspose, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}),
        A,
        B,
    )
end

@for_petsc function MatCreateKAIJ(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatCreateKAIJ, $petsc_library),
        PetscErrorCode,
        (
            Mat,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscScalar},
            Ptr{$PetscScalar},
            Ptr{Mat},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MatKAIJGetAIJ(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatKAIJGetAIJ, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function MatKAIJGetS(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatKAIJGetS, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatKAIJGetSRead(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatKAIJGetSRead, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatKAIJRestoreS(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatKAIJRestoreS, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatKAIJRestoreSRead(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatKAIJRestoreSRead, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatKAIJGetT(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatKAIJGetT, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatKAIJGetTRead(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatKAIJGetTRead, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatKAIJRestoreT(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatKAIJRestoreT, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatKAIJRestoreTRead(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatKAIJRestoreTRead, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatKAIJSetAIJ(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatKAIJSetAIJ, $petsc_library),
        PetscErrorCode,
        (Mat, Mat),
        arg1,
        arg2,
    )
end

@for_petsc function MatKAIJSetS(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatKAIJSetS, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatKAIJSetT(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatKAIJSetT, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatKAIJGetScaledIdentity(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatKAIJGetScaledIdentity, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function MatDiagonalScaleLocal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatDiagonalScaleLocal, $petsc_library),
        PetscErrorCode,
        (Mat, Vec),
        arg1,
        arg2,
    )
end

@for_petsc function MatCreateMFFD(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatCreateMFFD, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MatMFFDSetBase(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMFFDSetBase, $petsc_library),
        PetscErrorCode,
        (Mat, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMFFDSetFunction(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMFFDSetFunction, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMFFDSetFunctioni(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatMFFDSetFunctioni, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function MatMFFDSetFunctioniBase(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatMFFDSetFunctioniBase, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function MatMFFDSetHHistory(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMFFDSetHHistory, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscScalar}, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMFFDResetHHistory(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatMFFDResetHHistory, $petsc_library),
        PetscErrorCode,
        (Mat,),
        arg1,
    )
end

@for_petsc function MatMFFDSetFunctionError(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatMFFDSetFunctionError, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function MatMFFDSetPeriod(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatMFFDSetPeriod, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function MatMFFDGetH(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatMFFDGetH, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscScalar}),
        arg1,
        arg2,
    )
end

@for_petsc function MatMFFDSetOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatMFFDSetOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function MatMFFDCheckPositivity(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatMFFDCheckPositivity, $petsc_library),
        PetscErrorCode,
        (Ptr{Cvoid}, Vec, Vec, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatMFFDSetCheckh(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMFFDSetCheckh, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

mutable struct _p_MatMFFD end

const MatMFFD = Ptr{_p_MatMFFD}

const MatMFFDType = Ptr{Cchar}

@for_petsc function MatMFFDSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatMFFDSetType, $petsc_library),
        PetscErrorCode,
        (Mat, MatMFFDType),
        arg1,
        arg2,
    )
end

@for_petsc function MatMFFDRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatMFFDRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function MatMFFDDSSetUmin(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatMFFDDSSetUmin, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function MatMFFDWPSetComputeNormU(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatMFFDWPSetComputeNormU, $petsc_library),
        PetscErrorCode,
        (Mat, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function MatFDColoringSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatFDColoringSetType, $petsc_library),
        PetscErrorCode,
        (MatFDColoring, MatMFFDType),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerMathematicaPutMatrix(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscViewerMathematicaPutMatrix, $petsc_library),
        PetscErrorCode,
        (PetscViewer, $PetscInt, $PetscInt, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscViewerMathematicaPutCSRMatrix(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscViewerMathematicaPutCSRMatrix, $petsc_library),
        PetscErrorCode,
        (
            PetscViewer,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MatBindToCPU(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatBindToCPU, $petsc_library),
        PetscErrorCode,
        (Mat, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function MatPinToCPU(::$UnionPetscLib, A, flg)
    @chk ccall(
        (:MatPinToCPU, $petsc_library),
        PetscErrorCode,
        (Mat, PetscBool),
        A,
        flg,
    )
end

mutable struct _p_SplitCSRMat end

const PetscSplitCSRDataStructure = _p_SplitCSRMat

@for_petsc function MatCreateNest(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:MatCreateNest, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, Ptr{IS}, $PetscInt, Ptr{IS}, Ptr{Mat}, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function MatNestGetSize(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatNestGetSize, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatNestGetISs(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatNestGetISs, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{IS}, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatNestGetLocalISs(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatNestGetLocalISs, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{IS}, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatNestGetSubMats(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatNestGetSubMats, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{Ptr{Mat}}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatNestGetSubMat(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatNestGetSubMat, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatNestSetVecType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatNestSetVecType, $petsc_library),
        PetscErrorCode,
        (Mat, VecType),
        arg1,
        arg2,
    )
end

@for_petsc function MatNestSetSubMats(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatNestSetSubMats, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{IS}, $PetscInt, Ptr{IS}, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MatNestSetSubMat(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatNestSetSubMat, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt, Mat),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatChop(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatChop, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function MatComputeBandwidth(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatComputeBandwidth, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscReal, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatSubdomainsCreateCoalesce(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatSubdomainsCreateCoalesce, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}, Ptr{Ptr{IS}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatPreallocatorPreallocate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatPreallocatorPreallocate, $petsc_library),
        PetscErrorCode,
        (Mat, PetscBool, Mat),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatHeaderMerge(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatHeaderMerge, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function MatHeaderReplace(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatHeaderReplace, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}),
        arg1,
        arg2,
    )
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

@for_petsc function DMLabelCreate(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMLabelCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{DMLabel}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMLabelView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMLabelView, $petsc_library),
        PetscErrorCode,
        (DMLabel, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function DMLabelReset(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMLabelReset, $petsc_library),
        PetscErrorCode,
        (DMLabel,),
        arg1,
    )
end

@for_petsc function DMLabelDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMLabelDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{DMLabel},),
        arg1,
    )
end

@for_petsc function DMLabelGetDefaultValue(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMLabelGetDefaultValue, $petsc_library),
        PetscErrorCode,
        (DMLabel, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMLabelSetDefaultValue(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMLabelSetDefaultValue, $petsc_library),
        PetscErrorCode,
        (DMLabel, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function DMLabelDuplicate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMLabelDuplicate, $petsc_library),
        PetscErrorCode,
        (DMLabel, Ptr{DMLabel}),
        arg1,
        arg2,
    )
end

@for_petsc function DMLabelGetValue(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMLabelGetValue, $petsc_library),
        PetscErrorCode,
        (DMLabel, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMLabelSetValue(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMLabelSetValue, $petsc_library),
        PetscErrorCode,
        (DMLabel, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMLabelClearValue(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMLabelClearValue, $petsc_library),
        PetscErrorCode,
        (DMLabel, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMLabelAddStratum(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMLabelAddStratum, $petsc_library),
        PetscErrorCode,
        (DMLabel, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function DMLabelAddStrata(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMLabelAddStrata, $petsc_library),
        PetscErrorCode,
        (DMLabel, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMLabelAddStrataIS(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMLabelAddStrataIS, $petsc_library),
        PetscErrorCode,
        (DMLabel, IS),
        arg1,
        arg2,
    )
end

@for_petsc function DMLabelInsertIS(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMLabelInsertIS, $petsc_library),
        PetscErrorCode,
        (DMLabel, IS, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMLabelGetNumValues(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMLabelGetNumValues, $petsc_library),
        PetscErrorCode,
        (DMLabel, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMLabelGetStratumBounds(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMLabelGetStratumBounds, $petsc_library),
        PetscErrorCode,
        (DMLabel, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMLabelGetValueIS(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMLabelGetValueIS, $petsc_library),
        PetscErrorCode,
        (DMLabel, Ptr{IS}),
        arg1,
        arg2,
    )
end

@for_petsc function DMLabelStratumHasPoint(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMLabelStratumHasPoint, $petsc_library),
        PetscErrorCode,
        (DMLabel, $PetscInt, $PetscInt, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMLabelHasStratum(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMLabelHasStratum, $petsc_library),
        PetscErrorCode,
        (DMLabel, $PetscInt, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMLabelGetStratumSize(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMLabelGetStratumSize, $petsc_library),
        PetscErrorCode,
        (DMLabel, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMLabelGetStratumIS(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMLabelGetStratumIS, $petsc_library),
        PetscErrorCode,
        (DMLabel, $PetscInt, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMLabelSetStratumIS(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMLabelSetStratumIS, $petsc_library),
        PetscErrorCode,
        (DMLabel, $PetscInt, IS),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMLabelSetStratumBounds(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMLabelSetStratumBounds, $petsc_library),
        PetscErrorCode,
        (DMLabel, $PetscInt, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMLabelClearStratum(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMLabelClearStratum, $petsc_library),
        PetscErrorCode,
        (DMLabel, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function DMLabelComputeIndex(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMLabelComputeIndex, $petsc_library),
        PetscErrorCode,
        (DMLabel,),
        arg1,
    )
end

@for_petsc function DMLabelCreateIndex(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMLabelCreateIndex, $petsc_library),
        PetscErrorCode,
        (DMLabel, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMLabelDestroyIndex(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMLabelDestroyIndex, $petsc_library),
        PetscErrorCode,
        (DMLabel,),
        arg1,
    )
end

@for_petsc function DMLabelHasValue(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMLabelHasValue, $petsc_library),
        PetscErrorCode,
        (DMLabel, $PetscInt, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMLabelHasPoint(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMLabelHasPoint, $petsc_library),
        PetscErrorCode,
        (DMLabel, $PetscInt, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMLabelGetBounds(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMLabelGetBounds, $petsc_library),
        PetscErrorCode,
        (DMLabel, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMLabelFilter(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMLabelFilter, $petsc_library),
        PetscErrorCode,
        (DMLabel, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMLabelPermute(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMLabelPermute, $petsc_library),
        PetscErrorCode,
        (DMLabel, IS, Ptr{DMLabel}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMLabelDistribute(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMLabelDistribute, $petsc_library),
        PetscErrorCode,
        (DMLabel, PetscSF, Ptr{DMLabel}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMLabelGather(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMLabelGather, $petsc_library),
        PetscErrorCode,
        (DMLabel, PetscSF, Ptr{DMLabel}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMLabelConvertToSection(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMLabelConvertToSection, $petsc_library),
        PetscErrorCode,
        (DMLabel, Ptr{PetscSection}, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionCreateGlobalSectionLabel(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscSectionCreateGlobalSectionLabel, $petsc_library),
        PetscErrorCode,
        (
            PetscSection,
            PetscSF,
            PetscBool,
            DMLabel,
            $PetscInt,
            Ptr{PetscSection},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscSectionSymCreateLabel(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSectionSymCreateLabel, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, DMLabel, Ptr{PetscSectionSym}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionSymLabelSetLabel(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionSymLabelSetLabel, $petsc_library),
        PetscErrorCode,
        (PetscSectionSym, DMLabel),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionSymLabelSetStratum(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscSectionSymLabelSetStratum, $petsc_library),
        PetscErrorCode,
        (
            PetscSectionSym,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            PetscCopyMode,
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscScalar}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

mutable struct _p_PetscDS end

const PetscDS = Ptr{_p_PetscDS}

mutable struct _p_PetscWeakForm end

const PetscWeakForm = Ptr{_p_PetscWeakForm}

@for_petsc function DMInitializePackage(::$UnionPetscLib)
    @chk ccall((:DMInitializePackage, $petsc_library), PetscErrorCode, ())
end

const DMType = Ptr{Cchar}

@for_petsc function DMCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{DM}),
        arg1,
        arg2,
    )
end

@for_petsc function DMClone(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMClone, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DM}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetType, $petsc_library),
        PetscErrorCode,
        (DM, DMType),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetType, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DMType}),
        arg1,
        arg2,
    )
end

@for_petsc function DMRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMRegisterDestroy(::$UnionPetscLib)
    @chk ccall((:DMRegisterDestroy, $petsc_library), PetscErrorCode, ())
end

@for_petsc function DMView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMView, $petsc_library),
        PetscErrorCode,
        (DM, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function DMLoad(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMLoad, $petsc_library),
        PetscErrorCode,
        (DM, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function DMDestroy(::$UnionPetscLib, arg1)
    @chk ccall((:DMDestroy, $petsc_library), PetscErrorCode, (Ptr{DM},), arg1)
end

@for_petsc function DMCreateGlobalVector(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMCreateGlobalVector, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function DMCreateLocalVector(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMCreateLocalVector, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetLocalVector(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetLocalVector, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function DMRestoreLocalVector(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMRestoreLocalVector, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetGlobalVector(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetGlobalVector, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function DMRestoreGlobalVector(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMRestoreGlobalVector, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function DMClearGlobalVectors(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMClearGlobalVectors, $petsc_library),
        PetscErrorCode,
        (DM,),
        arg1,
    )
end

@for_petsc function DMClearLocalVectors(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMClearLocalVectors, $petsc_library),
        PetscErrorCode,
        (DM,),
        arg1,
    )
end

@for_petsc function DMHasNamedGlobalVector(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMHasNamedGlobalVector, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetNamedGlobalVector(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMGetNamedGlobalVector, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMRestoreNamedGlobalVector(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMRestoreNamedGlobalVector, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMHasNamedLocalVector(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMHasNamedLocalVector, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetNamedLocalVector(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMGetNamedLocalVector, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMRestoreNamedLocalVector(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMRestoreNamedLocalVector, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetLocalToGlobalMapping(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetLocalToGlobalMapping, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{ISLocalToGlobalMapping}),
        arg1,
        arg2,
    )
end

@for_petsc function DMCreateFieldIS(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMCreateFieldIS, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{Ptr{Ptr{Cchar}}}, Ptr{Ptr{IS}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMGetBlockSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetBlockSize, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMCreateColoring(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMCreateColoring, $petsc_library),
        PetscErrorCode,
        (DM, ISColoringType, Ptr{ISColoring}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMCreateMatrix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMCreateMatrix, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetMatrixPreallocateOnly(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetMatrixPreallocateOnly, $petsc_library),
        PetscErrorCode,
        (DM, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetMatrixStructureOnly(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetMatrixStructureOnly, $petsc_library),
        PetscErrorCode,
        (DM, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function DMCreateInterpolation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMCreateInterpolation, $petsc_library),
        PetscErrorCode,
        (DM, DM, Ptr{Mat}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMCreateRestriction(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMCreateRestriction, $petsc_library),
        PetscErrorCode,
        (DM, DM, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMCreateInjection(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMCreateInjection, $petsc_library),
        PetscErrorCode,
        (DM, DM, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMCreateMassMatrix(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMCreateMassMatrix, $petsc_library),
        PetscErrorCode,
        (DM, DM, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetWorkArray(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMGetWorkArray, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, MPI_Datatype, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMRestoreWorkArray(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMRestoreWorkArray, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, MPI_Datatype, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMRefine(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMRefine, $petsc_library),
        PetscErrorCode,
        (DM, MPI_Comm, Ptr{DM}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMCoarsen(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMCoarsen, $petsc_library),
        PetscErrorCode,
        (DM, MPI_Comm, Ptr{DM}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetCoarseDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetCoarseDM, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DM}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetCoarseDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetCoarseDM, $petsc_library),
        PetscErrorCode,
        (DM, DM),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetFineDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetFineDM, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DM}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetFineDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetFineDM, $petsc_library),
        PetscErrorCode,
        (DM, DM),
        arg1,
        arg2,
    )
end

@for_petsc function DMRefineHierarchy(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMRefineHierarchy, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{DM}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMCoarsenHierarchy(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMCoarsenHierarchy, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{DM}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMCoarsenHookAdd(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMCoarsenHookAdd, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMCoarsenHookRemove(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMCoarsenHookRemove, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMRefineHookAdd(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMRefineHookAdd, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMRefineHookRemove(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMRefineHookRemove, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMRestrict(::$UnionPetscLib, arg1, arg2, arg3, arg4, arg5)
    @chk ccall(
        (:DMRestrict, $petsc_library),
        PetscErrorCode,
        (DM, Mat, Vec, Mat, DM),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMInterpolate(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMInterpolate, $petsc_library),
        PetscErrorCode,
        (DM, Mat, DM),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMInterpolateSolution(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMInterpolateSolution, $petsc_library),
        PetscErrorCode,
        (DM, DM, Mat, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall((:DMSetFromOptions, $petsc_library), PetscErrorCode, (DM,), arg1)
end

@for_petsc function DMViewFromOptions(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMViewFromOptions, $petsc_library),
        PetscErrorCode,
        (DM, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMAdaptLabel(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMAdaptLabel, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel, Ptr{DM}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMAdaptMetric(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMAdaptMetric, $petsc_library),
        PetscErrorCode,
        (DM, Vec, DMLabel, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMSetUp(::$UnionPetscLib, arg1)
    @chk ccall((:DMSetUp, $petsc_library), PetscErrorCode, (DM,), arg1)
end

@for_petsc function DMCreateInterpolationScale(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMCreateInterpolationScale, $petsc_library),
        PetscErrorCode,
        (DM, DM, Mat, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMCreateAggregates(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMCreateAggregates, $petsc_library),
        PetscErrorCode,
        (DM, DM, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGlobalToLocalHookAdd(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMGlobalToLocalHookAdd, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMLocalToGlobalHookAdd(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMLocalToGlobalHookAdd, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMGlobalToLocal(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMGlobalToLocal, $petsc_library),
        PetscErrorCode,
        (DM, Vec, InsertMode, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMGlobalToLocalBegin(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMGlobalToLocalBegin, $petsc_library),
        PetscErrorCode,
        (DM, Vec, InsertMode, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMGlobalToLocalEnd(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMGlobalToLocalEnd, $petsc_library),
        PetscErrorCode,
        (DM, Vec, InsertMode, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMLocalToGlobal(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMLocalToGlobal, $petsc_library),
        PetscErrorCode,
        (DM, Vec, InsertMode, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMLocalToGlobalBegin(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMLocalToGlobalBegin, $petsc_library),
        PetscErrorCode,
        (DM, Vec, InsertMode, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMLocalToGlobalEnd(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMLocalToGlobalEnd, $petsc_library),
        PetscErrorCode,
        (DM, Vec, InsertMode, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMLocalToLocalBegin(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMLocalToLocalBegin, $petsc_library),
        PetscErrorCode,
        (DM, Vec, InsertMode, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMLocalToLocalEnd(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMLocalToLocalEnd, $petsc_library),
        PetscErrorCode,
        (DM, Vec, InsertMode, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMConvert(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMConvert, $petsc_library),
        PetscErrorCode,
        (DM, DMType, Ptr{DM}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetDimension(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetDimension, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetDimension(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetDimension, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetDimPoints(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMGetDimPoints, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMGetUseNatural(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetUseNatural, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetUseNatural(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetUseNatural, $petsc_library),
        PetscErrorCode,
        (DM, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetCoordinateDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetCoordinateDM, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DM}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetCoordinateDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetCoordinateDM, $petsc_library),
        PetscErrorCode,
        (DM, DM),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetCoordinateDim(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetCoordinateDim, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetCoordinateDim(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetCoordinateDim, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetCoordinateSection(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetCoordinateSection, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscSection}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetCoordinateSection(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMSetCoordinateSection, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, PetscSection),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetCoordinates(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetCoordinates, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetCoordinates(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetCoordinates, $petsc_library),
        PetscErrorCode,
        (DM, Vec),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetCoordinatesLocal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetCoordinatesLocal, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetCoordinatesLocalSetUp(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMGetCoordinatesLocalSetUp, $petsc_library),
        PetscErrorCode,
        (DM,),
        arg1,
    )
end

@for_petsc function DMGetCoordinatesLocalNoncollective(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:DMGetCoordinatesLocalNoncollective, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetCoordinatesLocalTuple(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMGetCoordinatesLocalTuple, $petsc_library),
        PetscErrorCode,
        (DM, IS, Ptr{PetscSection}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMSetCoordinatesLocal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetCoordinatesLocal, $petsc_library),
        PetscErrorCode,
        (DM, Vec),
        arg1,
        arg2,
    )
end

@for_petsc function DMLocatePoints(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMLocatePoints, $petsc_library),
        PetscErrorCode,
        (DM, Vec, DMPointLocationType, Ptr{PetscSF}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMGetPeriodicity(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMGetPeriodicity, $petsc_library),
        PetscErrorCode,
        (
            DM,
            Ptr{PetscBool},
            Ptr{Ptr{$PetscReal}},
            Ptr{Ptr{$PetscReal}},
            Ptr{Ptr{DMBoundaryType}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMSetPeriodicity(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMSetPeriodicity, $petsc_library),
        PetscErrorCode,
        (DM, PetscBool, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{DMBoundaryType}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMLocalizeCoordinate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMLocalizeCoordinate, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscScalar}, PetscBool, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMLocalizeCoordinates(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMLocalizeCoordinates, $petsc_library),
        PetscErrorCode,
        (DM,),
        arg1,
    )
end

@for_petsc function DMGetCoordinatesLocalized(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetCoordinatesLocalized, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetCoordinatesLocalizedLocal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetCoordinatesLocalizedLocal, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetNeighbors(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMGetNeighbors, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{Ptr{PetscMPIInt}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetCoordinateField(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetCoordinateField, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DMField}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetCoordinateField(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetCoordinateField, $petsc_library),
        PetscErrorCode,
        (DM, DMField),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetBoundingBox(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMGetBoundingBox, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetLocalBoundingBox(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMGetLocalBoundingBox, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMProjectCoordinates(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMProjectCoordinates, $petsc_library),
        PetscErrorCode,
        (DM, PetscFE),
        arg1,
        arg2,
    )
end

@for_petsc function DMSubDomainHookAdd(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMSubDomainHookAdd, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMSubDomainHookRemove(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMSubDomainHookRemove, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMSubDomainRestrict(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMSubDomainRestrict, $petsc_library),
        PetscErrorCode,
        (DM, VecScatter, VecScatter, DM),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMSetOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function DMAppendOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMAppendOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetVecType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetVecType, $petsc_library),
        PetscErrorCode,
        (DM, VecType),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetVecType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetVecType, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{VecType}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetMatType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetMatType, $petsc_library),
        PetscErrorCode,
        (DM, MatType),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetMatType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetMatType, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{MatType}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetISColoringType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetISColoringType, $petsc_library),
        PetscErrorCode,
        (DM, ISColoringType),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetISColoringType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetISColoringType, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{ISColoringType}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetApplicationContext(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetApplicationContext, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetApplicationContextDestroy(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetApplicationContextDestroy, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetApplicationContext(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetApplicationContext, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetVariableBounds(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetVariableBounds, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMHasVariableBounds(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMHasVariableBounds, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function DMHasColoring(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMHasColoring, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function DMHasCreateRestriction(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMHasCreateRestriction, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function DMHasCreateInjection(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMHasCreateInjection, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function DMComputeVariableBounds(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMComputeVariableBounds, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMCreateSubDM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMCreateSubDM, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}, Ptr{IS}, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMCreateSuperDM(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMCreateSuperDM, $petsc_library),
        PetscErrorCode,
        (Ptr{DM}, $PetscInt, Ptr{Ptr{IS}}, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMCreateSectionSubDM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMCreateSectionSubDM, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}, Ptr{IS}, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMCreateSectionSuperDM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMCreateSectionSuperDM, $petsc_library),
        PetscErrorCode,
        (Ptr{DM}, $PetscInt, Ptr{Ptr{IS}}, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMCreateFieldDecomposition(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMCreateFieldDecomposition, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{Ptr{Ptr{Cchar}}}, Ptr{Ptr{IS}}, Ptr{Ptr{DM}}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMCreateDomainDecomposition(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMCreateDomainDecomposition, $petsc_library),
        PetscErrorCode,
        (
            DM,
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cchar}}},
            Ptr{Ptr{IS}},
            Ptr{Ptr{IS}},
            Ptr{Ptr{DM}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMCreateDomainDecompositionScatters(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMCreateDomainDecompositionScatters, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            Ptr{DM},
            Ptr{Ptr{VecScatter}},
            Ptr{Ptr{VecScatter}},
            Ptr{Ptr{VecScatter}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMGetRefineLevel(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetRefineLevel, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetRefineLevel(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetRefineLevel, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetCoarsenLevel(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetCoarsenLevel, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetCoarsenLevel(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetCoarsenLevel, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function DMFinalizePackage(::$UnionPetscLib)
    @chk ccall((:DMFinalizePackage, $petsc_library), PetscErrorCode, ())
end

@for_petsc function VecGetDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecGetDM, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{DM}),
        arg1,
        arg2,
    )
end

@for_petsc function VecSetDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecSetDM, $petsc_library),
        PetscErrorCode,
        (Vec, DM),
        arg1,
        arg2,
    )
end

@for_petsc function MatGetDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetDM, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{DM}),
        arg1,
        arg2,
    )
end

@for_petsc function MatSetDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSetDM, $petsc_library),
        PetscErrorCode,
        (Mat, DM),
        arg1,
        arg2,
    )
end

@for_petsc function MatFDColoringUseDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatFDColoringUseDM, $petsc_library),
        PetscErrorCode,
        (Mat, MatFDColoring),
        arg1,
        arg2,
    )
end

mutable struct NLF_DAAD end

const NLF = Ptr{NLF_DAAD}

@for_petsc function DMPrintCellVector(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMPrintCellVector, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{Cchar}, $PetscInt, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPrintCellMatrix(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPrintCellMatrix, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{Cchar}, $PetscInt, $PetscInt, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPrintLocalVec(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMPrintLocalVec, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, $PetscReal, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMSetNullSpaceConstructor(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMSetNullSpaceConstructor, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetNullSpaceConstructor(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMGetNullSpaceConstructor, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMSetNearNullSpaceConstructor(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMSetNearNullSpaceConstructor, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetNearNullSpaceConstructor(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMGetNearNullSpaceConstructor, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetSection(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetSection, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscSection}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetSection(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetSection, $petsc_library),
        PetscErrorCode,
        (DM, PetscSection),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetLocalSection(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetLocalSection, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscSection}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetLocalSection(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetLocalSection, $petsc_library),
        PetscErrorCode,
        (DM, PetscSection),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetGlobalSection(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetGlobalSection, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscSection}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetGlobalSection(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetGlobalSection, $petsc_library),
        PetscErrorCode,
        (DM, PetscSection),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetDefaultSection(::$UnionPetscLib, dm, s)
    @chk ccall(
        (:DMGetDefaultSection, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscSection}),
        dm,
        s,
    )
end

@for_petsc function DMSetDefaultSection(::$UnionPetscLib, dm, s)
    @chk ccall(
        (:DMSetDefaultSection, $petsc_library),
        PetscErrorCode,
        (DM, PetscSection),
        dm,
        s,
    )
end

@for_petsc function DMGetDefaultGlobalSection(::$UnionPetscLib, dm, s)
    @chk ccall(
        (:DMGetDefaultGlobalSection, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscSection}),
        dm,
        s,
    )
end

@for_petsc function DMSetDefaultGlobalSection(::$UnionPetscLib, dm, s)
    @chk ccall(
        (:DMSetDefaultGlobalSection, $petsc_library),
        PetscErrorCode,
        (DM, PetscSection),
        dm,
        s,
    )
end

@for_petsc function DMGetSectionSF(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetSectionSF, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscSF}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetSectionSF(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetSectionSF, $petsc_library),
        PetscErrorCode,
        (DM, PetscSF),
        arg1,
        arg2,
    )
end

@for_petsc function DMCreateSectionSF(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMCreateSectionSF, $petsc_library),
        PetscErrorCode,
        (DM, PetscSection, PetscSection),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetDefaultSF(::$UnionPetscLib, dm, s)
    @chk ccall(
        (:DMGetDefaultSF, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscSF}),
        dm,
        s,
    )
end

@for_petsc function DMSetDefaultSF(::$UnionPetscLib, dm, s)
    @chk ccall(
        (:DMSetDefaultSF, $petsc_library),
        PetscErrorCode,
        (DM, PetscSF),
        dm,
        s,
    )
end

@for_petsc function DMCreateDefaultSF(::$UnionPetscLib, dm, l, g)
    @chk ccall(
        (:DMCreateDefaultSF, $petsc_library),
        PetscErrorCode,
        (DM, PetscSection, PetscSection),
        dm,
        l,
        g,
    )
end

@for_petsc function DMGetPointSF(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetPointSF, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscSF}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetPointSF(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetPointSF, $petsc_library),
        PetscErrorCode,
        (DM, PetscSF),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetDefaultConstraints(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMGetDefaultConstraints, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscSection}, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMSetDefaultConstraints(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMSetDefaultConstraints, $petsc_library),
        PetscErrorCode,
        (DM, PetscSection, Mat),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetOutputDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetOutputDM, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DM}),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetOutputSequenceNumber(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMGetOutputSequenceNumber, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMSetOutputSequenceNumber(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMSetOutputSequenceNumber, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscReal),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMOutputSequenceLoad(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMOutputSequenceLoad, $petsc_library),
        PetscErrorCode,
        (DM, PetscViewer, Ptr{Cchar}, $PetscInt, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMGetNumFields(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetNumFields, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetNumFields(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetNumFields, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetField(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMGetField, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{DMLabel}, Ptr{PetscObject}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMSetField(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMSetField, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, DMLabel, PetscObject),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMAddField(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMAddField, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel, PetscObject),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMSetFieldAvoidTensor(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMSetFieldAvoidTensor, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, PetscBool),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetFieldAvoidTensor(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMGetFieldAvoidTensor, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMClearFields(::$UnionPetscLib, arg1)
    @chk ccall((:DMClearFields, $petsc_library), PetscErrorCode, (DM,), arg1)
end

@for_petsc function DMCopyFields(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMCopyFields, $petsc_library),
        PetscErrorCode,
        (DM, DM),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetAdjacency(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMGetAdjacency, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{PetscBool}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMSetAdjacency(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMSetAdjacency, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, PetscBool, PetscBool),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMGetBasicAdjacency(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMGetBasicAdjacency, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscBool}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMSetBasicAdjacency(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMSetBasicAdjacency, $petsc_library),
        PetscErrorCode,
        (DM, PetscBool, PetscBool),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetNumDS(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetNumDS, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetDS(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetDS, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscDS}),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetCellDS(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMGetCellDS, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{PetscDS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetRegionDS(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMGetRegionDS, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel, Ptr{IS}, Ptr{PetscDS}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMSetRegionDS(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMSetRegionDS, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel, IS, PetscDS),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMGetRegionNumDS(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMGetRegionNumDS, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{DMLabel}, Ptr{IS}, Ptr{PetscDS}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMSetRegionNumDS(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMSetRegionNumDS, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, DMLabel, IS, PetscDS),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMFindRegionNum(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMFindRegionNum, $petsc_library),
        PetscErrorCode,
        (DM, PetscDS, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMCreateDS(::$UnionPetscLib, arg1)
    @chk ccall((:DMCreateDS, $petsc_library), PetscErrorCode, (DM,), arg1)
end

@for_petsc function DMClearDS(::$UnionPetscLib, arg1)
    @chk ccall((:DMClearDS, $petsc_library), PetscErrorCode, (DM,), arg1)
end

@for_petsc function DMCopyDS(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMCopyDS, $petsc_library),
        PetscErrorCode,
        (DM, DM),
        arg1,
        arg2,
    )
end

@for_petsc function DMCopyDisc(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMCopyDisc, $petsc_library),
        PetscErrorCode,
        (DM, DM),
        arg1,
        arg2,
    )
end

@for_petsc function DMComputeExactSolution(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMComputeExactSolution, $petsc_library),
        PetscErrorCode,
        (DM, $PetscReal, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMCreateLabel(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMCreateLabel, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetLabelValue(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMGetLabelValue, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMSetLabelValue(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMSetLabelValue, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMClearLabelValue(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMClearLabelValue, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMGetLabelSize(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMGetLabelSize, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetLabelIdIS(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMGetLabelIdIS, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetStratumSize(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMGetStratumSize, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMGetStratumIS(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMGetStratumIS, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, $PetscInt, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMSetStratumIS(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMSetStratumIS, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, $PetscInt, IS),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMClearLabelStratum(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMClearLabelStratum, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetLabelOutput(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMGetLabelOutput, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMSetLabelOutput(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMSetLabelOutput, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, PetscBool),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetNumLabels(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetNumLabels, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetLabelName(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMGetLabelName, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMHasLabel(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMHasLabel, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetLabel(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMGetLabel, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, Ptr{DMLabel}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetLabelByNum(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMGetLabelByNum, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{DMLabel}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMAddLabel(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMAddLabel, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel),
        arg1,
        arg2,
    )
end

@for_petsc function DMRemoveLabel(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMRemoveLabel, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, Ptr{DMLabel}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMRemoveLabelBySelf(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMRemoveLabelBySelf, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DMLabel}, PetscBool),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMCopyLabels(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMCopyLabels, $petsc_library),
        PetscErrorCode,
        (DM, DM, PetscCopyMode, PetscBool),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMAddBoundary(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
)
    @chk ccall(
        (:DMAddBoundary, $petsc_library),
        PetscErrorCode,
        (
            DM,
            DMBoundaryConditionType,
            Ptr{Cchar},
            Ptr{Cchar},
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Cvoid},
            Ptr{Cvoid},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
    )
end

@for_petsc function DMGetNumBoundary(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetNumBoundary, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetBoundary(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
)
    @chk ccall(
        (:DMGetBoundary, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            Ptr{DMBoundaryConditionType},
            Ptr{Ptr{Cchar}},
            Ptr{Ptr{Cchar}},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{Cvoid}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
    )
end

@for_petsc function DMIsBoundaryPoint(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMIsBoundaryPoint, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMCopyBoundary(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMCopyBoundary, $petsc_library),
        PetscErrorCode,
        (DM, DM),
        arg1,
        arg2,
    )
end

@for_petsc function DMProjectFunction(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMProjectFunction, $petsc_library),
        PetscErrorCode,
        (DM, $PetscReal, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}, InsertMode, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMProjectFunctionLocal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMProjectFunctionLocal, $petsc_library),
        PetscErrorCode,
        (DM, $PetscReal, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}, InsertMode, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMProjectFunctionLabel(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
)
    @chk ccall(
        (:DMProjectFunctionLabel, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscReal,
            DMLabel,
            $PetscInt,
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
            InsertMode,
            Vec,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
    )
end

@for_petsc function DMProjectFunctionLabelLocal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
)
    @chk ccall(
        (:DMProjectFunctionLabelLocal, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscReal,
            DMLabel,
            $PetscInt,
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
            InsertMode,
            Vec,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
    )
end

@for_petsc function DMProjectFieldLocal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMProjectFieldLocal, $petsc_library),
        PetscErrorCode,
        (DM, $PetscReal, Vec, Ptr{Ptr{Cvoid}}, InsertMode, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMProjectFieldLabelLocal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
)
    @chk ccall(
        (:DMProjectFieldLabelLocal, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscReal,
            DMLabel,
            $PetscInt,
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{$PetscInt},
            Vec,
            Ptr{Ptr{Cvoid}},
            InsertMode,
            Vec,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
    )
end

@for_petsc function DMProjectBdFieldLabelLocal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
)
    @chk ccall(
        (:DMProjectBdFieldLabelLocal, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscReal,
            DMLabel,
            $PetscInt,
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{$PetscInt},
            Vec,
            Ptr{Ptr{Cvoid}},
            InsertMode,
            Vec,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
    )
end

@for_petsc function DMComputeL2Diff(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMComputeL2Diff, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscReal,
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
            Vec,
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMComputeL2GradientDiff(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMComputeL2GradientDiff, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscReal,
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
            Vec,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMComputeL2FieldDiff(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMComputeL2FieldDiff, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscReal,
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
            Vec,
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMComputeError(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMComputeError, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{$PetscReal}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMHasBasisTransform(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMHasBasisTransform, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function DMCopyTransform(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMCopyTransform, $petsc_library),
        PetscErrorCode,
        (DM, DM),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetCompatibility(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMGetCompatibility, $petsc_library),
        PetscErrorCode,
        (DM, DM, Ptr{PetscBool}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMMonitorSet(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMMonitorSet, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMMonitorCancel(::$UnionPetscLib, arg1)
    @chk ccall((:DMMonitorCancel, $petsc_library), PetscErrorCode, (DM,), arg1)
end

@for_petsc function DMMonitorSetFromOptions(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMMonitorSetFromOptions, $petsc_library),
        PetscErrorCode,
        (
            DM,
            Ptr{Cchar},
            Ptr{Cchar},
            Ptr{Cchar},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{PetscBool},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMMonitor(::$UnionPetscLib, arg1)
    @chk ccall((:DMMonitor, $petsc_library), PetscErrorCode, (DM,), arg1)
end

@for_petsc function DMPolytopeTypeGetDim(::$UnionPetscLib, ct)
    ccall(
        (:DMPolytopeTypeGetDim, $petsc_library),
        $PetscInt,
        (DMPolytopeType,),
        ct,
    )
end

@for_petsc function DMPolytopeTypeGetConeSize(::$UnionPetscLib, ct)
    ccall(
        (:DMPolytopeTypeGetConeSize, $petsc_library),
        $PetscInt,
        (DMPolytopeType,),
        ct,
    )
end

@for_petsc function DMPolytopeTypeGetNumVertices(::$UnionPetscLib, ct)
    ccall(
        (:DMPolytopeTypeGetNumVertices, $petsc_library),
        $PetscInt,
        (DMPolytopeType,),
        ct,
    )
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

const PFType = Ptr{Cchar}

mutable struct _p_PF end

const PF = Ptr{_p_PF}

@for_petsc function PFCreate(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PFCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, $PetscInt, Ptr{PF}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PFSetType(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PFSetType, $petsc_library),
        PetscErrorCode,
        (PF, PFType, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PFSet(::$UnionPetscLib, arg1, arg2, arg3, arg4, arg5, arg6)
    @chk ccall(
        (:PFSet, $petsc_library),
        PetscErrorCode,
        (PF, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PFApply(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PFApply, $petsc_library),
        PetscErrorCode,
        (PF, $PetscInt, Ptr{$PetscScalar}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PFApplyVec(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PFApplyVec, $petsc_library),
        PetscErrorCode,
        (PF, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PFInitializePackage(::$UnionPetscLib)
    @chk ccall((:PFInitializePackage, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PFRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PFRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PFDestroy(::$UnionPetscLib, arg1)
    @chk ccall((:PFDestroy, $petsc_library), PetscErrorCode, (Ptr{PF},), arg1)
end

@for_petsc function PFSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall((:PFSetFromOptions, $petsc_library), PetscErrorCode, (PF,), arg1)
end

@for_petsc function PFGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PFGetType, $petsc_library),
        PetscErrorCode,
        (PF, Ptr{PFType}),
        arg1,
        arg2,
    )
end

@for_petsc function PFView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PFView, $petsc_library),
        PetscErrorCode,
        (PF, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PFViewFromOptions(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PFViewFromOptions, $petsc_library),
        PetscErrorCode,
        (PF, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

mutable struct _p_AO end

const AO = Ptr{_p_AO}

const AOType = Ptr{Cchar}

@for_petsc function AOInitializePackage(::$UnionPetscLib)
    @chk ccall((:AOInitializePackage, $petsc_library), PetscErrorCode, ())
end

@for_petsc function AOCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:AOCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{AO}),
        arg1,
        arg2,
    )
end

@for_petsc function AOSetIS(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:AOSetIS, $petsc_library),
        PetscErrorCode,
        (AO, IS, IS),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function AOSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall((:AOSetFromOptions, $petsc_library), PetscErrorCode, (AO,), arg1)
end

@for_petsc function AOCreateBasic(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:AOCreateBasic, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{AO}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function AOCreateBasicIS(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:AOCreateBasicIS, $petsc_library),
        PetscErrorCode,
        (IS, IS, Ptr{AO}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function AOCreateMemoryScalable(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:AOCreateMemoryScalable, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{AO}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function AOCreateMemoryScalableIS(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:AOCreateMemoryScalableIS, $petsc_library),
        PetscErrorCode,
        (IS, IS, Ptr{AO}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function AOCreateMapping(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:AOCreateMapping, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{AO}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function AOCreateMappingIS(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:AOCreateMappingIS, $petsc_library),
        PetscErrorCode,
        (IS, IS, Ptr{AO}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function AOView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:AOView, $petsc_library),
        PetscErrorCode,
        (AO, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function AOViewFromOptions(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:AOViewFromOptions, $petsc_library),
        PetscErrorCode,
        (AO, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function AODestroy(::$UnionPetscLib, arg1)
    @chk ccall((:AODestroy, $petsc_library), PetscErrorCode, (Ptr{AO},), arg1)
end

@for_petsc function AOSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:AOSetType, $petsc_library),
        PetscErrorCode,
        (AO, AOType),
        arg1,
        arg2,
    )
end

@for_petsc function AOGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:AOGetType, $petsc_library),
        PetscErrorCode,
        (AO, Ptr{AOType}),
        arg1,
        arg2,
    )
end

@for_petsc function AORegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:AORegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function AOPetscToApplication(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:AOPetscToApplication, $petsc_library),
        PetscErrorCode,
        (AO, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function AOApplicationToPetsc(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:AOApplicationToPetsc, $petsc_library),
        PetscErrorCode,
        (AO, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function AOPetscToApplicationIS(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:AOPetscToApplicationIS, $petsc_library),
        PetscErrorCode,
        (AO, IS),
        arg1,
        arg2,
    )
end

@for_petsc function AOApplicationToPetscIS(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:AOApplicationToPetscIS, $petsc_library),
        PetscErrorCode,
        (AO, IS),
        arg1,
        arg2,
    )
end

@for_petsc function AOPetscToApplicationPermuteInt(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:AOPetscToApplicationPermuteInt, $petsc_library),
        PetscErrorCode,
        (AO, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function AOApplicationToPetscPermuteInt(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:AOApplicationToPetscPermuteInt, $petsc_library),
        PetscErrorCode,
        (AO, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function AOPetscToApplicationPermuteReal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:AOPetscToApplicationPermuteReal, $petsc_library),
        PetscErrorCode,
        (AO, $PetscInt, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function AOApplicationToPetscPermuteReal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:AOApplicationToPetscPermuteReal, $petsc_library),
        PetscErrorCode,
        (AO, $PetscInt, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function AOMappingHasApplicationIndex(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:AOMappingHasApplicationIndex, $petsc_library),
        PetscErrorCode,
        (AO, $PetscInt, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function AOMappingHasPetscIndex(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:AOMappingHasPetscIndex, $petsc_library),
        PetscErrorCode,
        (AO, $PetscInt, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
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

@for_petsc function PetscQuadratureCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscQuadratureCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscQuadrature}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscQuadratureDuplicate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscQuadratureDuplicate, $petsc_library),
        PetscErrorCode,
        (PetscQuadrature, Ptr{PetscQuadrature}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscQuadratureGetOrder(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscQuadratureGetOrder, $petsc_library),
        PetscErrorCode,
        (PetscQuadrature, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscQuadratureSetOrder(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscQuadratureSetOrder, $petsc_library),
        PetscErrorCode,
        (PetscQuadrature, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscQuadratureGetNumComponents(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscQuadratureGetNumComponents, $petsc_library),
        PetscErrorCode,
        (PetscQuadrature, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscQuadratureSetNumComponents(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscQuadratureSetNumComponents, $petsc_library),
        PetscErrorCode,
        (PetscQuadrature, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscQuadratureGetData(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscQuadratureGetData, $petsc_library),
        PetscErrorCode,
        (
            PetscQuadrature,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscReal}},
            Ptr{Ptr{$PetscReal}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscQuadratureSetData(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscQuadratureSetData, $petsc_library),
        PetscErrorCode,
        (
            PetscQuadrature,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscQuadratureView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscQuadratureView, $petsc_library),
        PetscErrorCode,
        (PetscQuadrature, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscQuadratureDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscQuadratureDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscQuadrature},),
        arg1,
    )
end

@for_petsc function PetscQuadratureExpandComposite(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscQuadratureExpandComposite, $petsc_library),
        PetscErrorCode,
        (
            PetscQuadrature,
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{PetscQuadrature},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscQuadraturePushForward(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscQuadraturePushForward, $petsc_library),
        PetscErrorCode,
        (
            PetscQuadrature,
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            $PetscInt,
            Ptr{PetscQuadrature},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscDTLegendreEval(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscDTLegendreEval, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            Ptr{$PetscReal},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscDTJacobiNorm(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscDTJacobiNorm, $petsc_library),
        PetscErrorCode,
        ($PetscReal, $PetscReal, $PetscInt, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDTJacobiEval(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:PetscDTJacobiEval, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            $PetscReal,
            $PetscReal,
            Ptr{$PetscReal},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function PetscDTJacobiEvalJet(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscDTJacobiEvalJet, $petsc_library),
        PetscErrorCode,
        (
            $PetscReal,
            $PetscReal,
            $PetscInt,
            Ptr{$PetscReal},
            $PetscInt,
            $PetscInt,
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscDTPKDEvalJet(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscDTPKDEvalJet, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            $PetscInt,
            Ptr{$PetscReal},
            $PetscInt,
            $PetscInt,
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscDTGaussQuadrature(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDTGaussQuadrature, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscReal, $PetscReal, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDTGaussJacobiQuadrature(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscDTGaussJacobiQuadrature, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscDTGaussLobattoJacobiQuadrature(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscDTGaussLobattoJacobiQuadrature, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscDTGaussLobattoLegendreQuadrature(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDTGaussLobattoLegendreQuadrature, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            PetscGaussLobattoLegendreCreateType,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDTReconstructPoly(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscDTReconstructPoly, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            $PetscInt,
            Ptr{$PetscReal},
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscDTGaussTensorQuadrature(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscDTGaussTensorQuadrature, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscReal,
            $PetscReal,
            Ptr{PetscQuadrature},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscDTStroudConicalQuadrature(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscDTStroudConicalQuadrature, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscReal,
            $PetscReal,
            Ptr{PetscQuadrature},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscDTTanhSinhTensorQuadrature(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDTTanhSinhTensorQuadrature, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, $PetscReal, $PetscReal, Ptr{PetscQuadrature}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDTTanhSinhIntegrate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDTTanhSinhIntegrate, $petsc_library),
        PetscErrorCode,
        (Ptr{Cvoid}, $PetscReal, $PetscReal, $PetscInt, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDTTanhSinhIntegrateMPFR(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDTTanhSinhIntegrateMPFR, $petsc_library),
        PetscErrorCode,
        (Ptr{Cvoid}, $PetscReal, $PetscReal, $PetscInt, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscGaussLobattoLegendreIntegrate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscGaussLobattoLegendreIntegrate, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscGaussLobattoLegendreElementLaplacianCreate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscGaussLobattoLegendreElementLaplacianCreate, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{Ptr{Ptr{$PetscReal}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscGaussLobattoLegendreElementLaplacianDestroy(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscGaussLobattoLegendreElementLaplacianDestroy, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{Ptr{Ptr{$PetscReal}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscGaussLobattoLegendreElementGradientCreate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscGaussLobattoLegendreElementGradientCreate, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{Ptr{Ptr{$PetscReal}}},
            Ptr{Ptr{Ptr{$PetscReal}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscGaussLobattoLegendreElementGradientDestroy(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscGaussLobattoLegendreElementGradientDestroy, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{Ptr{Ptr{$PetscReal}}},
            Ptr{Ptr{Ptr{$PetscReal}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscGaussLobattoLegendreElementAdvectionCreate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscGaussLobattoLegendreElementAdvectionCreate, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{Ptr{Ptr{$PetscReal}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscGaussLobattoLegendreElementAdvectionDestroy(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscGaussLobattoLegendreElementAdvectionDestroy, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{Ptr{Ptr{$PetscReal}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscGaussLobattoLegendreElementMassCreate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscGaussLobattoLegendreElementMassCreate, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{Ptr{Ptr{$PetscReal}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscGaussLobattoLegendreElementMassDestroy(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscGaussLobattoLegendreElementMassDestroy, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{Ptr{Ptr{$PetscReal}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDTAltVApply(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDTAltVApply, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDTAltVWedge(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscDTAltVWedge, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscDTAltVWedgeMatrix(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDTAltVWedgeMatrix, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDTAltVPullback(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscDTAltVPullback, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            $PetscInt,
            Ptr{$PetscReal},
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscDTAltVPullbackMatrix(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDTAltVPullbackMatrix, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, Ptr{$PetscReal}, $PetscInt, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDTAltVInterior(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDTAltVInterior, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDTAltVInteriorMatrix(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDTAltVInteriorMatrix, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDTAltVInteriorPattern(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscDTAltVInteriorPattern, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, Ptr{NTuple{3, $PetscInt}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDTAltVStar(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscDTAltVStar, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscDTBaryToIndex(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscDTBaryToIndex, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDTIndexToBary(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscDTIndexToBary, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDTGradedOrderToIndex(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscDTGradedOrderToIndex, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDTIndexToGradedOrder(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscDTIndexToGradedOrder, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDTFactorial(::$UnionPetscLib, n, factorial)
    @chk ccall(
        (:PetscDTFactorial, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscReal}),
        n,
        factorial,
    )
end

@for_petsc function PetscDTFactorialInt(::$UnionPetscLib, n, factorial)
    @chk ccall(
        (:PetscDTFactorialInt, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscInt}),
        n,
        factorial,
    )
end

@for_petsc function PetscDTBinomial(::$UnionPetscLib, n, k, binomial)
    @chk ccall(
        (:PetscDTBinomial, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, Ptr{$PetscReal}),
        n,
        k,
        binomial,
    )
end

@for_petsc function PetscDTBinomialInt(::$UnionPetscLib, n, k, binomial)
    @chk ccall(
        (:PetscDTBinomialInt, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, Ptr{$PetscInt}),
        n,
        k,
        binomial,
    )
end

@for_petsc function PetscDTEnumPerm(::$UnionPetscLib, n, k, perm, isOdd)
    @chk ccall(
        (:PetscDTEnumPerm, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{PetscBool}),
        n,
        k,
        perm,
        isOdd,
    )
end

@for_petsc function PetscDTPermIndex(::$UnionPetscLib, n, perm, k, isOdd)
    @chk ccall(
        (:PetscDTPermIndex, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{PetscBool}),
        n,
        perm,
        k,
        isOdd,
    )
end

@for_petsc function PetscDTEnumSubset(::$UnionPetscLib, n, k, j, subset)
    @chk ccall(
        (:PetscDTEnumSubset, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}),
        n,
        k,
        j,
        subset,
    )
end

@for_petsc function PetscDTSubsetIndex(::$UnionPetscLib, n, k, subset, index)
    @chk ccall(
        (:PetscDTSubsetIndex, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        n,
        k,
        subset,
        index,
    )
end

@for_petsc function PetscDTEnumSplit(::$UnionPetscLib, n, k, j, perm, isOdd)
    @chk ccall(
        (:PetscDTEnumSplit, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{PetscBool}),
        n,
        k,
        j,
        perm,
        isOdd,
    )
end

@for_petsc function PetscFEInitializePackage(::$UnionPetscLib)
    @chk ccall((:PetscFEInitializePackage, $petsc_library), PetscErrorCode, ())
end

const PetscSpaceType = Ptr{Cchar}

@for_petsc function PetscSpaceCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSpaceCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscSpace}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSpaceDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscSpaceDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscSpace},),
        arg1,
    )
end

@for_petsc function PetscSpaceSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSpaceSetType, $petsc_library),
        PetscErrorCode,
        (PetscSpace, PetscSpaceType),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSpaceGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSpaceGetType, $petsc_library),
        PetscErrorCode,
        (PetscSpace, Ptr{PetscSpaceType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSpaceSetUp(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscSpaceSetUp, $petsc_library),
        PetscErrorCode,
        (PetscSpace,),
        arg1,
    )
end

@for_petsc function PetscSpaceSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscSpaceSetFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscSpace,),
        arg1,
    )
end

@for_petsc function PetscSpaceViewFromOptions(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSpaceViewFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscSpace, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSpaceView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSpaceView, $petsc_library),
        PetscErrorCode,
        (PetscSpace, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSpaceRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSpaceRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSpaceRegisterDestroy(::$UnionPetscLib)
    @chk ccall((:PetscSpaceRegisterDestroy, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscSpaceGetDimension(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSpaceGetDimension, $petsc_library),
        PetscErrorCode,
        (PetscSpace, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSpaceSetNumComponents(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSpaceSetNumComponents, $petsc_library),
        PetscErrorCode,
        (PetscSpace, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSpaceGetNumComponents(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSpaceGetNumComponents, $petsc_library),
        PetscErrorCode,
        (PetscSpace, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSpaceSetNumVariables(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSpaceSetNumVariables, $petsc_library),
        PetscErrorCode,
        (PetscSpace, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSpaceGetNumVariables(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSpaceGetNumVariables, $petsc_library),
        PetscErrorCode,
        (PetscSpace, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSpaceSetDegree(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSpaceSetDegree, $petsc_library),
        PetscErrorCode,
        (PetscSpace, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSpaceGetDegree(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSpaceGetDegree, $petsc_library),
        PetscErrorCode,
        (PetscSpace, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSpaceEvaluate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscSpaceEvaluate, $petsc_library),
        PetscErrorCode,
        (
            PetscSpace,
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscSpaceGetHeightSubspace(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSpaceGetHeightSubspace, $petsc_library),
        PetscErrorCode,
        (PetscSpace, $PetscInt, Ptr{PetscSpace}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSpacePolynomialSetSymmetric(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscSpacePolynomialSetSymmetric, $petsc_library),
        PetscErrorCode,
        (PetscSpace, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSpacePolynomialGetSymmetric(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscSpacePolynomialGetSymmetric, $petsc_library),
        PetscErrorCode,
        (PetscSpace, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSpacePolynomialSetTensor(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSpacePolynomialSetTensor, $petsc_library),
        PetscErrorCode,
        (PetscSpace, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSpacePolynomialGetTensor(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSpacePolynomialGetTensor, $petsc_library),
        PetscErrorCode,
        (PetscSpace, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSpaceTensorSetNumSubspaces(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscSpaceTensorSetNumSubspaces, $petsc_library),
        PetscErrorCode,
        (PetscSpace, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSpaceTensorGetNumSubspaces(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscSpaceTensorGetNumSubspaces, $petsc_library),
        PetscErrorCode,
        (PetscSpace, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSpaceTensorSetSubspace(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSpaceTensorSetSubspace, $petsc_library),
        PetscErrorCode,
        (PetscSpace, $PetscInt, PetscSpace),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSpaceTensorGetSubspace(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSpaceTensorGetSubspace, $petsc_library),
        PetscErrorCode,
        (PetscSpace, $PetscInt, Ptr{PetscSpace}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSpaceSumSetNumSubspaces(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSpaceSumSetNumSubspaces, $petsc_library),
        PetscErrorCode,
        (PetscSpace, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSpaceSumGetNumSubspaces(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSpaceSumGetNumSubspaces, $petsc_library),
        PetscErrorCode,
        (PetscSpace, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSpaceSumSetSubspace(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSpaceSumSetSubspace, $petsc_library),
        PetscErrorCode,
        (PetscSpace, $PetscInt, PetscSpace),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSpaceSumGetSubspace(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSpaceSumGetSubspace, $petsc_library),
        PetscErrorCode,
        (PetscSpace, $PetscInt, Ptr{PetscSpace}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSpaceSumSetConcatenate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSpaceSumSetConcatenate, $petsc_library),
        PetscErrorCode,
        (PetscSpace, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSpaceSumGetConcatenate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSpaceSumGetConcatenate, $petsc_library),
        PetscErrorCode,
        (PetscSpace, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSpaceCreateSum(
    ::$UnionPetscLib,
    numSubspaces,
    subspaces,
    concatenate,
    sumSpace,
)
    @chk ccall(
        (:PetscSpaceCreateSum, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{PetscSpace}, PetscBool, Ptr{PetscSpace}),
        numSubspaces,
        subspaces,
        concatenate,
        sumSpace,
    )
end

@for_petsc function PetscSpacePointGetPoints(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSpacePointGetPoints, $petsc_library),
        PetscErrorCode,
        (PetscSpace, Ptr{PetscQuadrature}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSpacePointSetPoints(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSpacePointSetPoints, $petsc_library),
        PetscErrorCode,
        (PetscSpace, PetscQuadrature),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSpaceCreateSubspace(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscSpaceCreateSubspace, $petsc_library),
        PetscErrorCode,
        (
            PetscSpace,
            PetscDualSpace,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            PetscCopyMode,
            Ptr{PetscSpace},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

const PetscDualSpaceType = Ptr{Cchar}

@for_petsc function PetscDualSpaceCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDualSpaceCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscDualSpace}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDualSpaceDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscDualSpace},),
        arg1,
    )
end

@for_petsc function PetscDualSpaceDuplicate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDualSpaceDuplicate, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{PetscDualSpace}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDualSpaceSetType, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, PetscDualSpaceType),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDualSpaceGetType, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{PetscDualSpaceType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceGetUniform(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDualSpaceGetUniform, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceGetNumDof(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDualSpaceGetNumDof, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceGetSection(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDualSpaceGetSection, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{PetscSection}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceSetUp(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDualSpaceSetUp, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace,),
        arg1,
    )
end

@for_petsc function PetscDualSpaceSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDualSpaceSetFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace,),
        arg1,
    )
end

@for_petsc function PetscDualSpaceViewFromOptions(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscDualSpaceViewFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDualSpaceView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDualSpaceView, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDualSpaceRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceRegisterDestroy(::$UnionPetscLib)
    @chk ccall(
        (:PetscDualSpaceRegisterDestroy, $petsc_library),
        PetscErrorCode,
        (),
    )
end

@for_petsc function PetscDualSpaceGetDimension(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDualSpaceGetDimension, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceGetInteriorDimension(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDualSpaceGetInteriorDimension, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceSetNumComponents(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDualSpaceSetNumComponents, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceGetNumComponents(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDualSpaceGetNumComponents, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceSetOrder(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDualSpaceSetOrder, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceGetOrder(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDualSpaceGetOrder, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceSetDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDualSpaceSetDM, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, DM),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceGetDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDualSpaceGetDM, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{DM}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceGetFunctional(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscDualSpaceGetFunctional, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, $PetscInt, Ptr{PetscQuadrature}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDualSpaceCreateReferenceCell(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDualSpaceCreateReferenceCell, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, $PetscInt, PetscBool, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDualSpaceGetSymmetries(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscDualSpaceGetSymmetries, $petsc_library),
        PetscErrorCode,
        (
            PetscDualSpace,
            Ptr{Ptr{Ptr{Ptr{$PetscInt}}}},
            Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}},
        ),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDualSpaceGetAllData(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDualSpaceGetAllData, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{PetscQuadrature}, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDualSpaceCreateAllDataDefault(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscDualSpaceCreateAllDataDefault, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{PetscQuadrature}, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDualSpaceGetInteriorData(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscDualSpaceGetInteriorData, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{PetscQuadrature}, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDualSpaceCreateInteriorDataDefault(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscDualSpaceCreateInteriorDataDefault, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{PetscQuadrature}, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDualSpaceApplyAll(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDualSpaceApplyAll, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{$PetscScalar}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDualSpaceApplyAllDefault(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscDualSpaceApplyAllDefault, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{$PetscScalar}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDualSpaceApplyInterior(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscDualSpaceApplyInterior, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{$PetscScalar}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDualSpaceApplyInteriorDefault(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscDualSpaceApplyInteriorDefault, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{$PetscScalar}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDualSpaceGetFormDegree(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDualSpaceGetFormDegree, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceSetFormDegree(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDualSpaceSetFormDegree, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceGetDeRahm(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDualSpaceGetDeRahm, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceLagrangeGetContinuity(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDualSpaceLagrangeGetContinuity, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceLagrangeSetContinuity(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDualSpaceLagrangeSetContinuity, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceLagrangeGetTensor(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDualSpaceLagrangeGetTensor, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceLagrangeSetTensor(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDualSpaceLagrangeSetTensor, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceLagrangeGetTrimmed(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDualSpaceLagrangeGetTrimmed, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceLagrangeSetTrimmed(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDualSpaceLagrangeSetTrimmed, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceLagrangeGetNodeType(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDualSpaceLagrangeGetNodeType, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{PetscDTNodeType}, Ptr{PetscBool}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDualSpaceLagrangeSetNodeType(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDualSpaceLagrangeSetNodeType, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, PetscDTNodeType, PetscBool, $PetscReal),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDualSpaceLagrangeGetUseMoments(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDualSpaceLagrangeGetUseMoments, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceLagrangeSetUseMoments(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDualSpaceLagrangeSetUseMoments, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceLagrangeGetMomentOrder(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDualSpaceLagrangeGetMomentOrder, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceLagrangeSetMomentOrder(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDualSpaceLagrangeSetMomentOrder, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceGetHeightSubspace(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscDualSpaceGetHeightSubspace, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, $PetscInt, Ptr{PetscDualSpace}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDualSpaceGetPointSubspace(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscDualSpaceGetPointSubspace, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, $PetscInt, Ptr{PetscDualSpace}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDualSpaceSimpleSetDimension(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDualSpaceSimpleSetDimension, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceSimpleSetFunctional(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscDualSpaceSimpleSetFunctional, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, $PetscInt, PetscQuadrature),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDualSpaceRefinedSetCellSpaces(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDualSpaceRefinedSetCellSpaces, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{PetscDualSpace}),
        arg1,
        arg2,
    )
end

const PetscFEType = Ptr{Cchar}

@for_petsc function PetscFECreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFECreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscFE}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFEDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscFEDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscFE},),
        arg1,
    )
end

@for_petsc function PetscFESetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFESetType, $petsc_library),
        PetscErrorCode,
        (PetscFE, PetscFEType),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFEGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFEGetType, $petsc_library),
        PetscErrorCode,
        (PetscFE, Ptr{PetscFEType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFESetUp(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscFESetUp, $petsc_library),
        PetscErrorCode,
        (PetscFE,),
        arg1,
    )
end

@for_petsc function PetscFESetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscFESetFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscFE,),
        arg1,
    )
end

@for_petsc function PetscFEViewFromOptions(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscFEViewFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscFE, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscFESetName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFESetName, $petsc_library),
        PetscErrorCode,
        (PetscFE, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFEView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFEView, $petsc_library),
        PetscErrorCode,
        (PetscFE, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFERegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFERegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFERegisterDestroy(::$UnionPetscLib)
    @chk ccall((:PetscFERegisterDestroy, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscFECreateDefault(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscFECreateDefault, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            PetscBool,
            Ptr{Cchar},
            $PetscInt,
            Ptr{PetscFE},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscFECreateLagrange(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscFECreateLagrange, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            PetscBool,
            $PetscInt,
            $PetscInt,
            Ptr{PetscFE},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscFEGetDimension(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFEGetDimension, $petsc_library),
        PetscErrorCode,
        (PetscFE, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFEGetSpatialDimension(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFEGetSpatialDimension, $petsc_library),
        PetscErrorCode,
        (PetscFE, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFESetNumComponents(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFESetNumComponents, $petsc_library),
        PetscErrorCode,
        (PetscFE, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFEGetNumComponents(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFEGetNumComponents, $petsc_library),
        PetscErrorCode,
        (PetscFE, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFEGetTileSizes(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscFEGetTileSizes, $petsc_library),
        PetscErrorCode,
        (
            PetscFE,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscFESetTileSizes(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscFESetTileSizes, $petsc_library),
        PetscErrorCode,
        (PetscFE, $PetscInt, $PetscInt, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscFESetBasisSpace(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFESetBasisSpace, $petsc_library),
        PetscErrorCode,
        (PetscFE, PetscSpace),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFEGetBasisSpace(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFEGetBasisSpace, $petsc_library),
        PetscErrorCode,
        (PetscFE, Ptr{PetscSpace}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFESetDualSpace(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFESetDualSpace, $petsc_library),
        PetscErrorCode,
        (PetscFE, PetscDualSpace),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFEGetDualSpace(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFEGetDualSpace, $petsc_library),
        PetscErrorCode,
        (PetscFE, Ptr{PetscDualSpace}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFESetQuadrature(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFESetQuadrature, $petsc_library),
        PetscErrorCode,
        (PetscFE, PetscQuadrature),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFEGetQuadrature(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFEGetQuadrature, $petsc_library),
        PetscErrorCode,
        (PetscFE, Ptr{PetscQuadrature}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFESetFaceQuadrature(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFESetFaceQuadrature, $petsc_library),
        PetscErrorCode,
        (PetscFE, PetscQuadrature),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFEGetFaceQuadrature(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFEGetFaceQuadrature, $petsc_library),
        PetscErrorCode,
        (PetscFE, Ptr{PetscQuadrature}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFECopyQuadrature(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFECopyQuadrature, $petsc_library),
        PetscErrorCode,
        (PetscFE, PetscFE),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFEGetNumDof(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFEGetNumDof, $petsc_library),
        PetscErrorCode,
        (PetscFE, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFERefine(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFERefine, $petsc_library),
        PetscErrorCode,
        (PetscFE, Ptr{PetscFE}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFEGetHeightSubspace(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscFEGetHeightSubspace, $petsc_library),
        PetscErrorCode,
        (PetscFE, $PetscInt, Ptr{PetscFE}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscFECompositeGetMapping(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscFECompositeGetMapping, $petsc_library),
        PetscErrorCode,
        (
            PetscFE,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscReal}},
            Ptr{Ptr{$PetscReal}},
            Ptr{Ptr{$PetscReal}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscFECreateHeightTrace(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscFECreateHeightTrace, $petsc_library),
        PetscErrorCode,
        (PetscFE, $PetscInt, Ptr{PetscFE}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscFECreatePointTrace(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscFECreatePointTrace, $petsc_library),
        PetscErrorCode,
        (PetscFE, $PetscInt, Ptr{PetscFE}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscFEOpenCLSetRealType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFEOpenCLSetRealType, $petsc_library),
        PetscErrorCode,
        (PetscFE, PetscDataType),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFEOpenCLGetRealType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFEOpenCLGetRealType, $petsc_library),
        PetscErrorCode,
        (PetscFE, Ptr{PetscDataType}),
        arg1,
        arg2,
    )
end

@for_petsc function DMDASetInterpolationType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDASetInterpolationType, $petsc_library),
        PetscErrorCode,
        (DM, DMDAInterpolationType),
        arg1,
        arg2,
    )
end

@for_petsc function DMDAGetInterpolationType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDAGetInterpolationType, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DMDAInterpolationType}),
        arg1,
        arg2,
    )
end

@for_petsc function DMDACreateAggregates(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDACreateAggregates, $petsc_library),
        PetscErrorCode,
        (DM, DM, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDASetElementType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDASetElementType, $petsc_library),
        PetscErrorCode,
        (DM, DMDAElementType),
        arg1,
        arg2,
    )
end

@for_petsc function DMDAGetElementType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDAGetElementType, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DMDAElementType}),
        arg1,
        arg2,
    )
end

@for_petsc function DMDAGetElements(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMDAGetElements, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMDARestoreElements(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMDARestoreElements, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMDAGetElementsSizes(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMDAGetElementsSizes, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMDAGetElementsCorners(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMDAGetElementsCorners, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMDAGetSubdomainCornersIS(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDAGetSubdomainCornersIS, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{IS}),
        arg1,
        arg2,
    )
end

@for_petsc function DMDARestoreSubdomainCornersIS(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDARestoreSubdomainCornersIS, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{IS}),
        arg1,
        arg2,
    )
end

@for_petsc function DMDACreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDACreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{DM}),
        arg1,
        arg2,
    )
end

@for_petsc function DMDASetSizes(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMDASetSizes, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMDACreate1d(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMDACreate1d, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            DMBoundaryType,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{DM},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMDACreate2d(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
)
    @chk ccall(
        (:DMDACreate2d, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            DMBoundaryType,
            DMBoundaryType,
            DMDAStencilType,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{DM},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
    )
end

@for_petsc function DMDACreate3d(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
    arg14,
    arg15,
    arg16,
    arg17,
)
    @chk ccall(
        (:DMDACreate3d, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            DMBoundaryType,
            DMBoundaryType,
            DMBoundaryType,
            DMDAStencilType,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{DM},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
        arg14,
        arg15,
        arg16,
        arg17,
    )
end

@for_petsc function DMDAGlobalToNaturalBegin(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMDAGlobalToNaturalBegin, $petsc_library),
        PetscErrorCode,
        (DM, Vec, InsertMode, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMDAGlobalToNaturalEnd(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMDAGlobalToNaturalEnd, $petsc_library),
        PetscErrorCode,
        (DM, Vec, InsertMode, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMDANaturalToGlobalBegin(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMDANaturalToGlobalBegin, $petsc_library),
        PetscErrorCode,
        (DM, Vec, InsertMode, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMDANaturalToGlobalEnd(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMDANaturalToGlobalEnd, $petsc_library),
        PetscErrorCode,
        (DM, Vec, InsertMode, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMDALocalToLocalBegin(::$UnionPetscLib, dm, g, mode, l)
    @chk ccall(
        (:DMDALocalToLocalBegin, $petsc_library),
        PetscErrorCode,
        (DM, Vec, InsertMode, Vec),
        dm,
        g,
        mode,
        l,
    )
end

@for_petsc function DMDALocalToLocalEnd(::$UnionPetscLib, dm, g, mode, l)
    @chk ccall(
        (:DMDALocalToLocalEnd, $petsc_library),
        PetscErrorCode,
        (DM, Vec, InsertMode, Vec),
        dm,
        g,
        mode,
        l,
    )
end

@for_petsc function DMDACreateNaturalVector(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDACreateNaturalVector, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function DMDAGetCorners(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMDAGetCorners, $petsc_library),
        PetscErrorCode,
        (
            DM,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMDAGetGhostCorners(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMDAGetGhostCorners, $petsc_library),
        PetscErrorCode,
        (
            DM,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMDAGetInfo(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
    arg14,
)
    @chk ccall(
        (:DMDAGetInfo, $petsc_library),
        PetscErrorCode,
        (
            DM,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{DMBoundaryType},
            Ptr{DMBoundaryType},
            Ptr{DMBoundaryType},
            Ptr{DMDAStencilType},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
        arg14,
    )
end

@for_petsc function DMDAGetProcessorSubset(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMDAGetProcessorSubset, $petsc_library),
        PetscErrorCode,
        (DM, DMDirection, $PetscInt, Ptr{MPI_Comm}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMDAGetProcessorSubsets(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDAGetProcessorSubsets, $petsc_library),
        PetscErrorCode,
        (DM, DMDirection, Ptr{MPI_Comm}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDAGetRay(::$UnionPetscLib, arg1, arg2, arg3, arg4, arg5)
    @chk ccall(
        (:DMDAGetRay, $petsc_library),
        PetscErrorCode,
        (DM, DMDirection, $PetscInt, Ptr{Vec}, Ptr{VecScatter}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMDAGlobalToNaturalAllCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDAGlobalToNaturalAllCreate, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{VecScatter}),
        arg1,
        arg2,
    )
end

@for_petsc function DMDANaturalAllToGlobalCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDANaturalAllToGlobalCreate, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{VecScatter}),
        arg1,
        arg2,
    )
end

@for_petsc function DMDAGetScatter(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDAGetScatter, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{VecScatter}, Ptr{VecScatter}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDAGetNeighbors(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDAGetNeighbors, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Ptr{PetscMPIInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function DMDASetAOType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDASetAOType, $petsc_library),
        PetscErrorCode,
        (DM, AOType),
        arg1,
        arg2,
    )
end

@for_petsc function DMDAGetAO(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDAGetAO, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{AO}),
        arg1,
        arg2,
    )
end

@for_petsc function DMDASetUniformCoordinates(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMDASetUniformCoordinates, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMDASetGLLCoordinates(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDASetGLLCoordinates, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDAGetCoordinateArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDAGetCoordinateArray, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMDARestoreCoordinateArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDARestoreCoordinateArray, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMDAGetLogicalCoordinate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
)
    @chk ccall(
        (:DMDAGetLogicalCoordinate, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscScalar,
            $PetscScalar,
            $PetscScalar,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscScalar},
            Ptr{$PetscScalar},
            Ptr{$PetscScalar},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
    )
end

@for_petsc function DMDAMapCoordsToPeriodicDomain(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMDAMapCoordsToPeriodicDomain, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscScalar}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDACreateCompatibleDMDA(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDACreateCompatibleDMDA, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{DM}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDAGetReducedDMDA(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDAGetReducedDMDA, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{DM}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDASetFieldName(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDASetFieldName, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDAGetFieldName(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDAGetFieldName, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDASetFieldNames(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDASetFieldNames, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function DMDAGetFieldNames(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDAGetFieldNames, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Ptr{Ptr{Cchar}}}),
        arg1,
        arg2,
    )
end

@for_petsc function DMDASetCoordinateName(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDASetCoordinateName, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDAGetCoordinateName(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDAGetCoordinateName, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDASetBoundaryType(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMDASetBoundaryType, $petsc_library),
        PetscErrorCode,
        (DM, DMBoundaryType, DMBoundaryType, DMBoundaryType),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMDASetDof(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDASetDof, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function DMDAGetDof(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDAGetDof, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMDASetOverlap(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMDASetOverlap, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMDAGetOverlap(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMDAGetOverlap, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMDASetNumLocalSubDomains(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDASetNumLocalSubDomains, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function DMDAGetNumLocalSubDomains(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDAGetNumLocalSubDomains, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMDAGetOffset(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMDAGetOffset, $petsc_library),
        PetscErrorCode,
        (
            DM,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMDASetOffset(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMDASetOffset, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMDAGetNonOverlappingRegion(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMDAGetNonOverlappingRegion, $petsc_library),
        PetscErrorCode,
        (
            DM,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMDASetNonOverlappingRegion(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMDASetNonOverlappingRegion, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMDASetStencilWidth(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDASetStencilWidth, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function DMDAGetStencilWidth(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDAGetStencilWidth, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMDASetOwnershipRanges(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMDASetOwnershipRanges, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMDAGetOwnershipRanges(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMDAGetOwnershipRanges, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMDASetNumProcs(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMDASetNumProcs, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMDASetStencilType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDASetStencilType, $petsc_library),
        PetscErrorCode,
        (DM, DMDAStencilType),
        arg1,
        arg2,
    )
end

@for_petsc function DMDAGetStencilType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDAGetStencilType, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DMDAStencilType}),
        arg1,
        arg2,
    )
end

@for_petsc function DMDAVecGetArray(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDAVecGetArray, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDAVecRestoreArray(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDAVecRestoreArray, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDAVecGetArrayWrite(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDAVecGetArrayWrite, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDAVecRestoreArrayWrite(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDAVecRestoreArrayWrite, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDAVecGetArrayDOF(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDAVecGetArrayDOF, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDAVecRestoreArrayDOF(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDAVecRestoreArrayDOF, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDAVecGetArrayRead(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDAVecGetArrayRead, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDAVecRestoreArrayRead(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDAVecRestoreArrayRead, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDAVecGetArrayDOFRead(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDAVecGetArrayDOFRead, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDAVecRestoreArrayDOFRead(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMDAVecRestoreArrayDOFRead, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatRegisterDAAD(::$UnionPetscLib)
    @chk ccall((:MatRegisterDAAD, $petsc_library), PetscErrorCode, ())
end

@for_petsc function MatCreateDAAD(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCreateDAAD, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function MatCreateSeqUSFFT(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatCreateSeqUSFFT, $petsc_library),
        PetscErrorCode,
        (Vec, DM, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDASetGetMatrix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDASetGetMatrix, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMDASetBlockFills(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDASetBlockFills, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDASetBlockFillsSparse(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDASetBlockFillsSparse, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDASetRefinementFactor(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMDASetRefinementFactor, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMDAGetRefinementFactor(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMDAGetRefinementFactor, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMDAGetArray(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDAGetArray, $petsc_library),
        PetscErrorCode,
        (DM, PetscBool, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDARestoreArray(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDARestoreArray, $petsc_library),
        PetscErrorCode,
        (DM, PetscBool, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDACreatePF(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDACreatePF, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PF}),
        arg1,
        arg2,
    )
end

@for_petsc function DMDAGetNumCells(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMDAGetNumCells, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMDAGetCellPoint(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMDAGetCellPoint, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMDAGetNumVertices(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMDAGetNumVertices, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMDAGetNumFaces(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMDAGetNumFaces, $petsc_library),
        PetscErrorCode,
        (
            DM,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMDAGetHeightStratum(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMDAGetHeightStratum, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMDAGetDepthStratum(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMDAGetDepthStratum, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMDAComputeCellGeometryFEM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMDAComputeCellGeometryFEM, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            PetscQuadrature,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMDAGetTransitiveClosure(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMDAGetTransitiveClosure, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, PetscBool, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMDARestoreTransitiveClosure(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMDARestoreTransitiveClosure, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, PetscBool, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMDASetVertexCoordinates(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMDASetVertexCoordinates, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMDASetPreallocationCenterDimension(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:DMDASetPreallocationCenterDimension, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function DMDAGetPreallocationCenterDimension(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:DMDAGetPreallocationCenterDimension, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMCompositeCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMCompositeCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{DM}),
        arg1,
        arg2,
    )
end

@for_petsc function DMCompositeAddDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMCompositeAddDM, $petsc_library),
        PetscErrorCode,
        (DM, DM),
        arg1,
        arg2,
    )
end

@for_petsc function DMCompositeSetCoupling(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMCompositeSetCoupling, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMCompositeAddVecScatter(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMCompositeAddVecScatter, $petsc_library),
        PetscErrorCode,
        (DM, VecScatter),
        arg1,
        arg2,
    )
end

@for_petsc function DMCompositeScatterArray(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMCompositeScatterArray, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMCompositeGatherArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMCompositeGatherArray, $petsc_library),
        PetscErrorCode,
        (DM, InsertMode, Vec, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMCompositeGetNumberDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMCompositeGetNumberDM, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMCompositeGetAccessArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMCompositeGetAccessArray, $petsc_library),
        PetscErrorCode,
        (DM, Vec, $PetscInt, Ptr{$PetscInt}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMCompositeRestoreAccessArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMCompositeRestoreAccessArray, $petsc_library),
        PetscErrorCode,
        (DM, Vec, $PetscInt, Ptr{$PetscInt}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMCompositeGetLocalAccessArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMCompositeGetLocalAccessArray, $petsc_library),
        PetscErrorCode,
        (DM, Vec, $PetscInt, Ptr{$PetscInt}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMCompositeRestoreLocalAccessArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMCompositeRestoreLocalAccessArray, $petsc_library),
        PetscErrorCode,
        (DM, Vec, $PetscInt, Ptr{$PetscInt}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMCompositeGetEntriesArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMCompositeGetEntriesArray, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DM}),
        arg1,
        arg2,
    )
end

@for_petsc function DMCompositeGetGlobalISs(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMCompositeGetGlobalISs, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Ptr{IS}}),
        arg1,
        arg2,
    )
end

@for_petsc function DMCompositeGetLocalISs(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMCompositeGetLocalISs, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Ptr{IS}}),
        arg1,
        arg2,
    )
end

@for_petsc function DMCompositeGetISLocalToGlobalMappings(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:DMCompositeGetISLocalToGlobalMappings, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Ptr{ISLocalToGlobalMapping}}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPatchCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPatchCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{DM}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPatchSolve(::$UnionPetscLib, arg1)
    @chk ccall((:DMPatchSolve, $petsc_library), PetscErrorCode, (DM,), arg1)
end

@for_petsc function DMPatchGetCoarse(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPatchGetCoarse, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DM}),
        arg1,
        arg2,
    )
end

mutable struct _p_PetscPartitioner end

const PetscPartitioner = Ptr{_p_PetscPartitioner}

@for_petsc function PetscPartitionerInitializePackage(::$UnionPetscLib)
    @chk ccall(
        (:PetscPartitionerInitializePackage, $petsc_library),
        PetscErrorCode,
        (),
    )
end

@for_petsc function PetscPartitionerFinalizePackage(::$UnionPetscLib)
    @chk ccall(
        (:PetscPartitionerFinalizePackage, $petsc_library),
        PetscErrorCode,
        (),
    )
end

const PetscPartitionerType = Ptr{Cchar}

@for_petsc function PetscPartitionerRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscPartitionerRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscPartitionerCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscPartitionerCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscPartitioner}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscPartitionerDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscPartitionerDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscPartitioner},),
        arg1,
    )
end

@for_petsc function PetscPartitionerSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscPartitionerSetType, $petsc_library),
        PetscErrorCode,
        (PetscPartitioner, PetscPartitionerType),
        arg1,
        arg2,
    )
end

@for_petsc function PetscPartitionerGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscPartitionerGetType, $petsc_library),
        PetscErrorCode,
        (PetscPartitioner, Ptr{PetscPartitionerType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscPartitionerSetUp(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscPartitionerSetUp, $petsc_library),
        PetscErrorCode,
        (PetscPartitioner,),
        arg1,
    )
end

@for_petsc function PetscPartitionerReset(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscPartitionerReset, $petsc_library),
        PetscErrorCode,
        (PetscPartitioner,),
        arg1,
    )
end

@for_petsc function PetscPartitionerSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscPartitionerSetFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscPartitioner,),
        arg1,
    )
end

@for_petsc function PetscPartitionerViewFromOptions(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscPartitionerViewFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscPartitioner, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscPartitionerView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscPartitionerView, $petsc_library),
        PetscErrorCode,
        (PetscPartitioner, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscPartitionerPartition(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:PetscPartitionerPartition, $petsc_library),
        PetscErrorCode,
        (
            PetscPartitioner,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            PetscSection,
            PetscSection,
            PetscSection,
            Ptr{IS},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function PetscPartitionerShellSetPartition(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscPartitionerShellSetPartition, $petsc_library),
        PetscErrorCode,
        (PetscPartitioner, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscPartitionerShellSetRandom(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscPartitionerShellSetRandom, $petsc_library),
        PetscErrorCode,
        (PetscPartitioner, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscPartitionerShellGetRandom(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscPartitionerShellGetRandom, $petsc_library),
        PetscErrorCode,
        (PetscPartitioner, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscPartitionerMatPartitioningGetMatPartitioning(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscPartitionerMatPartitioningGetMatPartitioning, $petsc_library),
        PetscErrorCode,
        (PetscPartitioner, Ptr{MatPartitioning}),
        arg1,
        arg2,
    )
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

const PetscLimiterType = Ptr{Cchar}

@for_petsc function PetscLimiterCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLimiterCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscLimiter}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLimiterDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLimiterDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscLimiter},),
        arg1,
    )
end

@for_petsc function PetscLimiterSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLimiterSetType, $petsc_library),
        PetscErrorCode,
        (PetscLimiter, PetscLimiterType),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLimiterGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLimiterGetType, $petsc_library),
        PetscErrorCode,
        (PetscLimiter, Ptr{PetscLimiterType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLimiterSetUp(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLimiterSetUp, $petsc_library),
        PetscErrorCode,
        (PetscLimiter,),
        arg1,
    )
end

@for_petsc function PetscLimiterSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLimiterSetFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscLimiter,),
        arg1,
    )
end

@for_petsc function PetscLimiterViewFromOptions(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscLimiterViewFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscLimiter, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLimiterView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLimiterView, $petsc_library),
        PetscErrorCode,
        (PetscLimiter, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLimiterRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLimiterRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLimiterRegisterDestroy(::$UnionPetscLib)
    @chk ccall(
        (:PetscLimiterRegisterDestroy, $petsc_library),
        PetscErrorCode,
        (),
    )
end

@for_petsc function PetscLimiterLimit(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscLimiterLimit, $petsc_library),
        PetscErrorCode,
        (PetscLimiter, $PetscReal, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscFVInitializePackage(::$UnionPetscLib)
    @chk ccall((:PetscFVInitializePackage, $petsc_library), PetscErrorCode, ())
end

const PetscFVType = Ptr{Cchar}

@for_petsc function PetscFVCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFVCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscFV}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFVDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscFVDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscFV},),
        arg1,
    )
end

@for_petsc function PetscFVSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFVSetType, $petsc_library),
        PetscErrorCode,
        (PetscFV, PetscFVType),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFVGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFVGetType, $petsc_library),
        PetscErrorCode,
        (PetscFV, Ptr{PetscFVType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFVSetUp(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscFVSetUp, $petsc_library),
        PetscErrorCode,
        (PetscFV,),
        arg1,
    )
end

@for_petsc function PetscFVSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscFVSetFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscFV,),
        arg1,
    )
end

@for_petsc function PetscFVViewFromOptions(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscFVViewFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscFV, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscFVView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFVView, $petsc_library),
        PetscErrorCode,
        (PetscFV, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFVRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFVRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFVRegisterDestroy(::$UnionPetscLib)
    @chk ccall((:PetscFVRegisterDestroy, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscFVSetComponentName(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscFVSetComponentName, $petsc_library),
        PetscErrorCode,
        (PetscFV, $PetscInt, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscFVGetComponentName(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscFVGetComponentName, $petsc_library),
        PetscErrorCode,
        (PetscFV, $PetscInt, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscFVSetLimiter(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFVSetLimiter, $petsc_library),
        PetscErrorCode,
        (PetscFV, PetscLimiter),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFVGetLimiter(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFVGetLimiter, $petsc_library),
        PetscErrorCode,
        (PetscFV, Ptr{PetscLimiter}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFVSetNumComponents(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFVSetNumComponents, $petsc_library),
        PetscErrorCode,
        (PetscFV, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFVGetNumComponents(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFVGetNumComponents, $petsc_library),
        PetscErrorCode,
        (PetscFV, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFVSetSpatialDimension(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFVSetSpatialDimension, $petsc_library),
        PetscErrorCode,
        (PetscFV, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFVGetSpatialDimension(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFVGetSpatialDimension, $petsc_library),
        PetscErrorCode,
        (PetscFV, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFVSetComputeGradients(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFVSetComputeGradients, $petsc_library),
        PetscErrorCode,
        (PetscFV, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFVGetComputeGradients(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFVGetComputeGradients, $petsc_library),
        PetscErrorCode,
        (PetscFV, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFVSetQuadrature(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFVSetQuadrature, $petsc_library),
        PetscErrorCode,
        (PetscFV, PetscQuadrature),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFVGetQuadrature(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFVGetQuadrature, $petsc_library),
        PetscErrorCode,
        (PetscFV, Ptr{PetscQuadrature}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFVSetDualSpace(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFVSetDualSpace, $petsc_library),
        PetscErrorCode,
        (PetscFV, PetscDualSpace),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFVGetDualSpace(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFVGetDualSpace, $petsc_library),
        PetscErrorCode,
        (PetscFV, Ptr{PetscDualSpace}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFVRefine(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFVRefine, $petsc_library),
        PetscErrorCode,
        (PetscFV, Ptr{PetscFV}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFVComputeGradient(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscFVComputeGradient, $petsc_library),
        PetscErrorCode,
        (PetscFV, $PetscInt, Ptr{$PetscScalar}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscFVLeastSquaresSetMaxFaces(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFVLeastSquaresSetMaxFaces, $petsc_library),
        PetscErrorCode,
        (PetscFV, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function DMFieldInitializePackage(::$UnionPetscLib)
    @chk ccall((:DMFieldInitializePackage, $petsc_library), PetscErrorCode, ())
end

@for_petsc function DMFieldFinalizePackage(::$UnionPetscLib)
    @chk ccall((:DMFieldFinalizePackage, $petsc_library), PetscErrorCode, ())
end

const DMFieldType = Ptr{Cchar}

@for_petsc function DMFieldSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMFieldSetType, $petsc_library),
        PetscErrorCode,
        (DMField, DMFieldType),
        arg1,
        arg2,
    )
end

@for_petsc function DMFieldGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMFieldGetType, $petsc_library),
        PetscErrorCode,
        (DMField, Ptr{DMFieldType}),
        arg1,
        arg2,
    )
end

@for_petsc function DMFieldRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMFieldRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@enum DMFieldContinuity::UInt32 begin
    DMFIELD_VERTEX = 0
    DMFIELD_EDGE = 1
    DMFIELD_FACET = 2
    DMFIELD_CELL = 3
end

@for_petsc function DMFieldDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMFieldDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{DMField},),
        arg1,
    )
end

@for_petsc function DMFieldView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMFieldView, $petsc_library),
        PetscErrorCode,
        (DMField, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function DMFieldGetDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMFieldGetDM, $petsc_library),
        PetscErrorCode,
        (DMField, Ptr{DM}),
        arg1,
        arg2,
    )
end

@for_petsc function DMFieldGetNumComponents(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMFieldGetNumComponents, $petsc_library),
        PetscErrorCode,
        (DMField, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMFieldGetContinuity(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMFieldGetContinuity, $petsc_library),
        PetscErrorCode,
        (DMField, Ptr{DMFieldContinuity}),
        arg1,
        arg2,
    )
end

@for_petsc function DMFieldEvaluate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMFieldEvaluate, $petsc_library),
        PetscErrorCode,
        (DMField, Vec, PetscDataType, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMFieldEvaluateFE(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMFieldEvaluateFE, $petsc_library),
        PetscErrorCode,
        (
            DMField,
            IS,
            PetscQuadrature,
            PetscDataType,
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMFieldEvaluateFV(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMFieldEvaluateFV, $petsc_library),
        PetscErrorCode,
        (DMField, IS, PetscDataType, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMFieldCreateDefaultQuadrature(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMFieldCreateDefaultQuadrature, $petsc_library),
        PetscErrorCode,
        (DMField, IS, Ptr{PetscQuadrature}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMFieldGetDegree(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMFieldGetDegree, $petsc_library),
        PetscErrorCode,
        (DMField, IS, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMFieldCreateDA(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMFieldCreateDA, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscScalar}, Ptr{DMField}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMFieldCreateDS(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMFieldCreateDS, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Vec, Ptr{DMField}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMFieldCreateShell(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMFieldCreateShell, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, DMFieldContinuity, Ptr{Cvoid}, Ptr{DMField}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMFieldShellSetDestroy(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMFieldShellSetDestroy, $petsc_library),
        PetscErrorCode,
        (DMField, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMFieldShellGetContext(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMFieldShellGetContext, $petsc_library),
        PetscErrorCode,
        (DMField, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMFieldShellSetEvaluate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMFieldShellSetEvaluate, $petsc_library),
        PetscErrorCode,
        (DMField, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMFieldShellSetEvaluateFE(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMFieldShellSetEvaluateFE, $petsc_library),
        PetscErrorCode,
        (DMField, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMFieldShellEvaluateFEDefault(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMFieldShellEvaluateFEDefault, $petsc_library),
        PetscErrorCode,
        (
            DMField,
            IS,
            PetscQuadrature,
            PetscDataType,
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMFieldShellSetEvaluateFV(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMFieldShellSetEvaluateFV, $petsc_library),
        PetscErrorCode,
        (DMField, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMFieldShellEvaluateFVDefault(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMFieldShellEvaluateFVDefault, $petsc_library),
        PetscErrorCode,
        (DMField, IS, PetscDataType, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMFieldShellSetGetDegree(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMFieldShellSetGetDegree, $petsc_library),
        PetscErrorCode,
        (DMField, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMFieldShellSetCreateDefaultQuadrature(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:DMFieldShellSetCreateDefaultQuadrature, $petsc_library),
        PetscErrorCode,
        (DMField, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscPartitionerDMPlexPartition(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscPartitionerDMPlexPartition, $petsc_library),
        PetscErrorCode,
        (PetscPartitioner, DM, PetscSection, PetscSection, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexBuildFromCellList(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexBuildFromCellList, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexBuildFromCellListParallel(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMPlexBuildFromCellListParallel, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{PetscSF},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMPlexBuildCoordinatesFromCellList(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMPlexBuildCoordinatesFromCellList, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexBuildCoordinatesFromCellListParallel(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexBuildCoordinatesFromCellListParallel, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, PetscSF, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{DM}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexCreateCohesiveSubmesh(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexCreateCohesiveSubmesh, $petsc_library),
        PetscErrorCode,
        (DM, PetscBool, Ptr{Cchar}, $PetscInt, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexCreateFromCellListPetsc(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
)
    @chk ccall(
        (:DMPlexCreateFromCellListPetsc, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            PetscBool,
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{DM},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
    )
end

@for_petsc function DMPlexCreateFromCellList(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
)
    @chk ccall(
        (:DMPlexCreateFromCellList, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            PetscBool,
            Ptr{Cint},
            $PetscInt,
            Ptr{Cdouble},
            Ptr{DM},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
    )
end

@for_petsc function DMPlexCreateFromCellListParallelPetsc(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
)
    @chk ccall(
        (:DMPlexCreateFromCellListParallelPetsc, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            PetscBool,
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{PetscSF},
            Ptr{DM},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
    )
end

@for_petsc function DMPlexCreateFromCellListParallel(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
)
    @chk ccall(
        (:DMPlexCreateFromCellListParallel, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            PetscBool,
            Ptr{Cint},
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{PetscSF},
            Ptr{DM},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
    )
end

@for_petsc function DMPlexCreateFromDAG(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMPlexCreateFromDAG, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscScalar},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMPlexCreateReferenceCell(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexCreateReferenceCell, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, PetscBool, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexCreateReferenceCellByType(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMPlexCreateReferenceCellByType, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, DMPolytopeType, Ptr{DM}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexSetOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexSetOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetChart(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexGetChart, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexSetChart(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexSetChart, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexGetConeSize(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexGetConeSize, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexSetConeSize(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexSetConeSize, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexAddConeSize(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexAddConeSize, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexGetCone(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexGetCone, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexGetConeTuple(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMPlexGetConeTuple, $petsc_library),
        PetscErrorCode,
        (DM, IS, Ptr{PetscSection}, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexGetConeRecursive(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexGetConeRecursive, $petsc_library),
        PetscErrorCode,
        (DM, IS, Ptr{$PetscInt}, Ptr{Ptr{IS}}, Ptr{Ptr{PetscSection}}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexRestoreConeRecursive(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexRestoreConeRecursive, $petsc_library),
        PetscErrorCode,
        (DM, IS, Ptr{$PetscInt}, Ptr{Ptr{IS}}, Ptr{Ptr{PetscSection}}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexGetConeRecursiveVertices(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMPlexGetConeRecursiveVertices, $petsc_library),
        PetscErrorCode,
        (DM, IS, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexSetCone(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexSetCone, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexInsertCone(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMPlexInsertCone, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexInsertConeOrientation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexInsertConeOrientation, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexGetConeOrientation(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexGetConeOrientation, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexSetConeOrientation(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexSetConeOrientation, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexGetSupportSize(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexGetSupportSize, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexSetSupportSize(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexSetSupportSize, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexGetSupport(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexGetSupport, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexSetSupport(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexSetSupport, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexInsertSupport(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexInsertSupport, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexGetConeSection(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetConeSection, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscSection}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetSupportSection(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetSupportSection, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscSection}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetCones(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetCones, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetConeOrientations(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetConeOrientations, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetMaxSizes(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexGetMaxSizes, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexSymmetrize(::$UnionPetscLib, arg1)
    @chk ccall((:DMPlexSymmetrize, $petsc_library), PetscErrorCode, (DM,), arg1)
end

@for_petsc function DMPlexStratify(::$UnionPetscLib, arg1)
    @chk ccall((:DMPlexStratify, $petsc_library), PetscErrorCode, (DM,), arg1)
end

@for_petsc function DMPlexEqual(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexEqual, $petsc_library),
        PetscErrorCode,
        (DM, DM, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexReverseCell(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexReverseCell, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexOrientCell(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMPlexOrientCell, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexCompareOrientations(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMPlexCompareOrientations, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{PetscBool},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMPlexOrient(::$UnionPetscLib, arg1)
    @chk ccall((:DMPlexOrient, $petsc_library), PetscErrorCode, (DM,), arg1)
end

@for_petsc function DMPlexPreallocateOperator(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:DMPlexPreallocateOperator, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Mat,
            PetscBool,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function DMPlexGetPointLocal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexGetPointLocal, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexPointLocalRead(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexPointLocalRead, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscScalar}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexPointLocalRef(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexPointLocalRef, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscScalar}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexGetPointLocalField(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexGetPointLocalField, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexPointLocalFieldRef(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexPointLocalFieldRef, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexPointLocalFieldRead(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexPointLocalFieldRead, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexGetPointGlobal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexGetPointGlobal, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexPointGlobalRead(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexPointGlobalRead, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscScalar}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexPointGlobalRef(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexPointGlobalRef, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscScalar}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexGetPointGlobalField(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexGetPointGlobalField, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexPointGlobalFieldRef(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexPointGlobalFieldRef, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexPointGlobalFieldRead(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexPointGlobalFieldRead, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@enum DMPlexInterpolatedFlag::Int32 begin
    DMPLEX_INTERPOLATED_INVALID = -1
    DMPLEX_INTERPOLATED_NONE = 0
    DMPLEX_INTERPOLATED_PARTIAL = 1
    DMPLEX_INTERPOLATED_MIXED = 2
    DMPLEX_INTERPOLATED_FULL = 3
end

@for_petsc function DMPlexInterpolate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexInterpolate, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DM}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexUninterpolate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexUninterpolate, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DM}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexInterpolatePointSF(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexInterpolatePointSF, $petsc_library),
        PetscErrorCode,
        (DM, PetscSF),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexIsInterpolated(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexIsInterpolated, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DMPlexInterpolatedFlag}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexIsInterpolatedCollective(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexIsInterpolatedCollective, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DMPlexInterpolatedFlag}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexFilter(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMPlexFilter, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel, $PetscInt, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexGetCellNumbering(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetCellNumbering, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{IS}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetVertexNumbering(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetVertexNumbering, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{IS}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexCreatePointNumbering(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexCreatePointNumbering, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{IS}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexCreateRankField(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexCreateRankField, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexCreateLabelField(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexCreateLabelField, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexGetDepth(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetDepth, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetDepthLabel(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetDepthLabel, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DMLabel}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetDepthStratum(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexGetDepthStratum, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexGetHeightStratum(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexGetHeightStratum, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexGetPointDepth(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexGetPointDepth, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexGetPointHeight(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexGetPointHeight, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexGetCellTypeLabel(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetCellTypeLabel, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DMLabel}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetCellType(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexGetCellType, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{DMPolytopeType}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexSetCellType(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexSetCellType, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, DMPolytopeType),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexComputeCellTypes(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMPlexComputeCellTypes, $petsc_library),
        PetscErrorCode,
        (DM,),
        arg1,
    )
end

@for_petsc function DMPlexInvertCell(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexInvertCell, $petsc_library),
        PetscErrorCode,
        (DMPolytopeType, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexReorderCell(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexReorderCell, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexGetMeet(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexGetMeet, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexGetFullMeet(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexGetFullMeet, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexRestoreMeet(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexRestoreMeet, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexGetJoin(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexGetJoin, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexGetFullJoin(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexGetFullJoin, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexRestoreJoin(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexRestoreJoin, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexGetTransitiveClosure(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexGetTransitiveClosure, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, PetscBool, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexRestoreTransitiveClosure(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexRestoreTransitiveClosure, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, PetscBool, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexGenerate(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMPlexGenerate, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, PetscBool, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexGenerateRegister(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexGenerateRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexGenerateRegisterAll(::$UnionPetscLib)
    @chk ccall((:DMPlexGenerateRegisterAll, $petsc_library), PetscErrorCode, ())
end

@for_petsc function DMPlexGenerateRegisterDestroy(::$UnionPetscLib)
    @chk ccall(
        (:DMPlexGenerateRegisterDestroy, $petsc_library),
        PetscErrorCode,
        (),
    )
end

@for_petsc function DMPlexCopyCoordinates(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexCopyCoordinates, $petsc_library),
        PetscErrorCode,
        (DM, DM),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexCreateDoublet(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMPlexCreateDoublet, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, PetscBool, PetscBool, $PetscReal, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMPlexCreateSquareBoundary(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexCreateSquareBoundary, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexCreateCubeBoundary(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexCreateCubeBoundary, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexCreateBoxMesh(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:DMPlexCreateBoxMesh, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            PetscBool,
            Ptr{$PetscInt},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{DMBoundaryType},
            PetscBool,
            Ptr{DM},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function DMPlexCreateSphereMesh(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexCreateSphereMesh, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, PetscBool, $PetscReal, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexCreateBallMesh(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexCreateBallMesh, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, $PetscReal, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexCreateHexCylinderMesh(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexCreateHexCylinderMesh, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, DMBoundaryType, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexCreateWedgeCylinderMesh(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexCreateWedgeCylinderMesh, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, PetscBool, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexCreateWedgeBoxMesh(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:DMPlexCreateWedgeBoxMesh, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            Ptr{$PetscInt},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{DMBoundaryType},
            PetscBool,
            PetscBool,
            Ptr{DM},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function DMPlexExtrude(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMPlexExtrude, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            $PetscReal,
            PetscBool,
            Ptr{$PetscReal},
            PetscBool,
            Ptr{DM},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMPlexCreateConeSection(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexCreateConeSection, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscSection}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexCheckSymmetry(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMPlexCheckSymmetry, $petsc_library),
        PetscErrorCode,
        (DM,),
        arg1,
    )
end

@for_petsc function DMPlexCheckSkeleton(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexCheckSkeleton, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexCheckFaces(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexCheckFaces, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexCheckGeometry(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMPlexCheckGeometry, $petsc_library),
        PetscErrorCode,
        (DM,),
        arg1,
    )
end

@for_petsc function DMPlexCheckPointSF(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMPlexCheckPointSF, $petsc_library),
        PetscErrorCode,
        (DM,),
        arg1,
    )
end

@for_petsc function DMPlexCheckInterfaceCones(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMPlexCheckInterfaceCones, $petsc_library),
        PetscErrorCode,
        (DM,),
        arg1,
    )
end

@for_petsc function DMPlexCheckCellShape(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexCheckCellShape, $petsc_library),
        PetscErrorCode,
        (DM, PetscBool, $PetscReal),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexComputeOrthogonalQuality(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexComputeOrthogonalQuality, $petsc_library),
        PetscErrorCode,
        (DM, PetscFV, $PetscReal, Ptr{Vec}, Ptr{DMLabel}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexTriangleSetOptions(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexTriangleSetOptions, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexTetgenSetOptions(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexTetgenSetOptions, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexCreateFromFile(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexCreateFromFile, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, PetscBool, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexCreateExodus(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMPlexCreateExodus, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, PetscBool, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexCreateExodusFromFile(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexCreateExodusFromFile, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, PetscBool, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexCreateCGNS(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMPlexCreateCGNS, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, PetscBool, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexCreateCGNSFromFile(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexCreateCGNSFromFile, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, PetscBool, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexCreateGmsh(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMPlexCreateGmsh, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, PetscViewer, PetscBool, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexCreateGmshFromFile(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexCreateGmshFromFile, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, PetscBool, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexCreateFluent(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMPlexCreateFluent, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, PetscViewer, PetscBool, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexCreateFluentFromFile(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexCreateFluentFromFile, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, PetscBool, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexCreateMedFromFile(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexCreateMedFromFile, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, PetscBool, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexCreatePLYFromFile(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexCreatePLYFromFile, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, PetscBool, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexCreateEGADSFromFile(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMPlexCreateEGADSFromFile, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{DM}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscViewerExodusIIGetId(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerExodusIIGetId, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{Cint}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexCreateNeighborCSR(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexCreateNeighborCSR, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscInt}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexGetPartitioner(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetPartitioner, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscPartitioner}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexSetPartitioner(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexSetPartitioner, $petsc_library),
        PetscErrorCode,
        (DM, PetscPartitioner),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexCreatePartitionerGraph(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMPlexCreatePartitionerGraph, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscInt}},
            Ptr{IS},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMPlexPartitionLabelInvert(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexPartitionLabelInvert, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel, PetscSF, DMLabel),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexPartitionLabelClosure(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexPartitionLabelClosure, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexPartitionLabelAdjacency(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexPartitionLabelAdjacency, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexPartitionLabelPropagate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexPartitionLabelPropagate, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexPartitionLabelCreateSF(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMPlexPartitionLabelCreateSF, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel, Ptr{PetscSF}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexSetPartitionBalance(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexSetPartitionBalance, $petsc_library),
        PetscErrorCode,
        (DM, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetPartitionBalance(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetPartitionBalance, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexIsDistributed(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexIsDistributed, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexDistribute(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMPlexDistribute, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{PetscSF}, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexDistributeOverlap(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexDistributeOverlap, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{PetscSF}, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexGetOverlap(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetOverlap, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexDistributeField(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMPlexDistributeField, $petsc_library),
        PetscErrorCode,
        (DM, PetscSF, PetscSection, Vec, PetscSection, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMPlexDistributeFieldIS(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMPlexDistributeFieldIS, $petsc_library),
        PetscErrorCode,
        (DM, PetscSF, PetscSection, IS, PetscSection, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMPlexDistributeData(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMPlexDistributeData, $petsc_library),
        PetscErrorCode,
        (
            DM,
            PetscSF,
            PetscSection,
            MPI_Datatype,
            Ptr{Cvoid},
            PetscSection,
            Ptr{Ptr{Cvoid}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMPlexRebalanceSharedPoints(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexRebalanceSharedPoints, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, PetscBool, PetscBool, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexMigrate(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexMigrate, $petsc_library),
        PetscErrorCode,
        (DM, PetscSF, DM),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexGetGatherDM(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexGetGatherDM, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscSF}, Ptr{DM}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexGetRedundantDM(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexGetRedundantDM, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscSF}, Ptr{DM}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexSetAdjacencyUser(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexSetAdjacencyUser, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexGetAdjacencyUser(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexGetAdjacencyUser, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexSetAdjacencyUseAnchors(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexSetAdjacencyUseAnchors, $petsc_library),
        PetscErrorCode,
        (DM, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetAdjacencyUseAnchors(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetAdjacencyUseAnchors, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetAdjacency(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMPlexGetAdjacency, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexSetMigrationSF(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexSetMigrationSF, $petsc_library),
        PetscErrorCode,
        (DM, PetscSF),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetMigrationSF(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetMigrationSF, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscSF}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetOrdering(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMPlexGetOrdering, $petsc_library),
        PetscErrorCode,
        (DM, MatOrderingType, DMLabel, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexPermute(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexPermute, $petsc_library),
        PetscErrorCode,
        (DM, IS, Ptr{DM}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexCreateProcessSF(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexCreateProcessSF, $petsc_library),
        PetscErrorCode,
        (DM, PetscSF, Ptr{IS}, Ptr{PetscSF}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexCreateTwoSidedProcessSF(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:DMPlexCreateTwoSidedProcessSF, $petsc_library),
        PetscErrorCode,
        (
            DM,
            PetscSF,
            PetscSection,
            IS,
            PetscSection,
            IS,
            Ptr{IS},
            Ptr{PetscSF},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function DMPlexDistributeOwnership(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexDistributeOwnership, $petsc_library),
        PetscErrorCode,
        (DM, PetscSection, Ptr{IS}, PetscSection, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexCreateOverlapLabel(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMPlexCreateOverlapLabel, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, PetscSection, IS, PetscSection, IS, Ptr{DMLabel}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMPlexCreateOverlapMigrationSF(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMPlexCreateOverlapMigrationSF, $petsc_library),
        PetscErrorCode,
        (DM, PetscSF, Ptr{PetscSF}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexStratifyMigrationSF(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMPlexStratifyMigrationSF, $petsc_library),
        PetscErrorCode,
        (DM, PetscSF, Ptr{PetscSF}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexCreateSubmesh(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexCreateSubmesh, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel, $PetscInt, PetscBool, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexCreateHybridMesh(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMPlexCreateHybridMesh, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel, DMLabel, Ptr{DMLabel}, Ptr{DMLabel}, Ptr{DM}, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMPlexGetSubpointMap(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetSubpointMap, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DMLabel}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexSetSubpointMap(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexSetSubpointMap, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetSubpointIS(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetSubpointIS, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{IS}),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetEnclosureRelation(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMGetEnclosureRelation, $petsc_library),
        PetscErrorCode,
        (DM, DM, Ptr{DMEnclosureType}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMGetEnclosurePoint(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMGetEnclosurePoint, $petsc_library),
        PetscErrorCode,
        (DM, DM, DMEnclosureType, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexLabelComplete(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexLabelComplete, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexLabelCohesiveComplete(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexLabelCohesiveComplete, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel, DMLabel, PetscBool, DM),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexLabelAddCells(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexLabelAddCells, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexLabelAddFaceCells(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexLabelAddFaceCells, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexLabelClearCells(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexLabelClearCells, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetRefinementLimit(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetRefinementLimit, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexSetRefinementLimit(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexSetRefinementLimit, $petsc_library),
        PetscErrorCode,
        (DM, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetRefinementUniform(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetRefinementUniform, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexSetRefinementUniform(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexSetRefinementUniform, $petsc_library),
        PetscErrorCode,
        (DM, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetRefinementFunction(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetRefinementFunction, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexSetRefinementFunction(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexSetRefinementFunction, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexCreateCoarsePointIS(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexCreateCoarsePointIS, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{IS}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetRegularRefinement(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetRegularRefinement, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexSetRegularRefinement(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexSetRegularRefinement, $petsc_library),
        PetscErrorCode,
        (DM, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetCellRefinerType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetCellRefinerType, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DMPlexCellRefinerType}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexSetCellRefinerType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexSetCellRefinerType, $petsc_library),
        PetscErrorCode,
        (DM, DMPlexCellRefinerType),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetNumFaceVertices(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexGetNumFaceVertices, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexGetOrientedFace(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:DMPlexGetOrientedFace, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{PetscBool},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function DMPlexGetMinRadius(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetMinRadius, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexSetMinRadius(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexSetMinRadius, $petsc_library),
        PetscErrorCode,
        (DM, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexComputeProjection2Dto1D(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexComputeProjection2Dto1D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscScalar}, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexComputeProjection3Dto1D(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexComputeProjection3Dto1D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscScalar}, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexComputeProjection3Dto2D(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMPlexComputeProjection3Dto2D, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscScalar}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

mutable struct _PetscGridHash end

const PetscGridHash = Ptr{_PetscGridHash}

@for_petsc function PetscGridHashCreate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscGridHashCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, Ptr{$PetscScalar}, Ptr{PetscGridHash}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscGridHashEnlarge(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscGridHashEnlarge, $petsc_library),
        PetscErrorCode,
        (PetscGridHash, Ptr{$PetscScalar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscGridHashSetGrid(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscGridHashSetGrid, $petsc_library),
        PetscErrorCode,
        (PetscGridHash, Ptr{$PetscInt}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscGridHashGetEnclosingBox(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscGridHashGetEnclosingBox, $petsc_library),
        PetscErrorCode,
        (
            PetscGridHash,
            $PetscInt,
            Ptr{$PetscScalar},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscGridHashDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscGridHashDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscGridHash},),
        arg1,
    )
end

@for_petsc function DMPlexFindVertices(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexFindVertices, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscReal}, $PetscReal, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexComputeCellGeometryFVM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexComputeCellGeometryFVM, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexComputeGeometryFVM(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexComputeGeometryFVM, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Vec}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexComputeGradientFVM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexComputeGradientFVM, $petsc_library),
        PetscErrorCode,
        (DM, PetscFV, Vec, Vec, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexGetDataFVM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexGetDataFVM, $petsc_library),
        PetscErrorCode,
        (DM, PetscFV, Ptr{Vec}, Ptr{Vec}, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexComputeGeometryFEM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexComputeGeometryFEM, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetGeometryFVM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexGetGeometryFVM, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Vec}, Ptr{Vec}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexGetGradientDM(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexGetGradientDM, $petsc_library),
        PetscErrorCode,
        (DM, PetscFV, Ptr{DM}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexInsertBoundaryValues(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMPlexInsertBoundaryValues, $petsc_library),
        PetscErrorCode,
        (DM, PetscBool, Vec, $PetscReal, Vec, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMPlexInsertTimeDerivativeBoundaryValues(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMPlexInsertTimeDerivativeBoundaryValues, $petsc_library),
        PetscErrorCode,
        (DM, PetscBool, Vec, $PetscReal, Vec, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMPlexInsertBoundaryValuesEssential(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
)
    @chk ccall(
        (:DMPlexInsertBoundaryValuesEssential, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscReal,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            DMLabel,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Vec,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
    )
end

@for_petsc function DMPlexInsertBoundaryValuesEssentialField(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
)
    @chk ccall(
        (:DMPlexInsertBoundaryValuesEssentialField, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscReal,
            Vec,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            DMLabel,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Vec,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
    )
end

@for_petsc function DMPlexInsertBoundaryValuesEssentialBdField(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
)
    @chk ccall(
        (:DMPlexInsertBoundaryValuesEssentialBdField, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscReal,
            Vec,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            DMLabel,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Vec,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
    )
end

@for_petsc function DMPlexInsertBoundaryValuesRiemann(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
    arg14,
)
    @chk ccall(
        (:DMPlexInsertBoundaryValuesRiemann, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscReal,
            Vec,
            Vec,
            Vec,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            DMLabel,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Vec,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
        arg14,
    )
end

@for_petsc function DMPlexMarkBoundaryFaces(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexMarkBoundaryFaces, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, DMLabel),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexCreateSection(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
)
    @chk ccall(
        (:DMPlexCreateSection, $petsc_library),
        PetscErrorCode,
        (
            DM,
            Ptr{DMLabel},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{IS},
            Ptr{IS},
            IS,
            Ptr{PetscSection},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
    )
end

@for_petsc function DMPlexGetSubdomainSection(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetSubdomainSection, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscSection}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexComputeCellGeometryAffineFEM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMPlexComputeCellGeometryAffineFEM, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMPlexComputeCellGeometryFEM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMPlexComputeCellGeometryFEM, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            PetscQuadrature,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMPlexCoordinatesToReference(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexCoordinatesToReference, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexReferenceToCoordinates(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexReferenceToCoordinates, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexShearGeometry(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexShearGeometry, $petsc_library),
        PetscErrorCode,
        (DM, DMDirection, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexRemapGeometry(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexRemapGeometry, $petsc_library),
        PetscErrorCode,
        (DM, $PetscReal, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexVecGetClosure(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMPlexVecGetClosure, $petsc_library),
        PetscErrorCode,
        (
            DM,
            PetscSection,
            Vec,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscScalar}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMPlexVecRestoreClosure(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMPlexVecRestoreClosure, $petsc_library),
        PetscErrorCode,
        (
            DM,
            PetscSection,
            Vec,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscScalar}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMPlexVecSetClosure(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMPlexVecSetClosure, $petsc_library),
        PetscErrorCode,
        (DM, PetscSection, Vec, $PetscInt, Ptr{$PetscScalar}, InsertMode),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMPlexMatSetClosure(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMPlexMatSetClosure, $petsc_library),
        PetscErrorCode,
        (
            DM,
            PetscSection,
            PetscSection,
            Mat,
            $PetscInt,
            Ptr{$PetscScalar},
            InsertMode,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMPlexMatSetClosureGeneral(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
)
    @chk ccall(
        (:DMPlexMatSetClosureGeneral, $petsc_library),
        PetscErrorCode,
        (
            DM,
            PetscSection,
            PetscSection,
            DM,
            PetscSection,
            PetscSection,
            Mat,
            $PetscInt,
            Ptr{$PetscScalar},
            InsertMode,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
    )
end

@for_petsc function DMPlexGetClosureIndices(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:DMPlexGetClosureIndices, $petsc_library),
        PetscErrorCode,
        (
            DM,
            PetscSection,
            PetscSection,
            $PetscInt,
            PetscBool,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscScalar}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function DMPlexRestoreClosureIndices(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:DMPlexRestoreClosureIndices, $petsc_library),
        PetscErrorCode,
        (
            DM,
            PetscSection,
            PetscSection,
            $PetscInt,
            PetscBool,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscScalar}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function DMPlexMatSetClosureRefined(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
)
    @chk ccall(
        (:DMPlexMatSetClosureRefined, $petsc_library),
        PetscErrorCode,
        (
            DM,
            PetscSection,
            PetscSection,
            DM,
            PetscSection,
            PetscSection,
            Mat,
            $PetscInt,
            Ptr{$PetscScalar},
            InsertMode,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
    )
end

@for_petsc function DMPlexMatGetClosureIndicesRefined(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:DMPlexMatGetClosureIndicesRefined, $petsc_library),
        PetscErrorCode,
        (
            DM,
            PetscSection,
            PetscSection,
            DM,
            PetscSection,
            PetscSection,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function DMPlexCreateClosureIndex(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexCreateClosureIndex, $petsc_library),
        PetscErrorCode,
        (DM, PetscSection),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexSetClosurePermutationTensor(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMPlexSetClosurePermutationTensor, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, PetscSection),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexConstructGhostCells(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexConstructGhostCells, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, Ptr{$PetscInt}, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexConstructCohesiveCells(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexConstructCohesiveCells, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel, DMLabel, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexGetVTKCellHeight(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetVTKCellHeight, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexSetVTKCellHeight(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexSetVTKCellHeight, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexVTKWriteAll(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexVTKWriteAll, $petsc_library),
        PetscErrorCode,
        (PetscObject, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetGhostCellStratum(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMPlexGetGhostCellStratum, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexGetSimplexOrBoxCells(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexGetSimplexOrBoxCells, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexGetCellFields(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:DMPlexGetCellFields, $petsc_library),
        PetscErrorCode,
        (
            DM,
            IS,
            Vec,
            Vec,
            Vec,
            Ptr{Ptr{$PetscScalar}},
            Ptr{Ptr{$PetscScalar}},
            Ptr{Ptr{$PetscScalar}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function DMPlexRestoreCellFields(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:DMPlexRestoreCellFields, $petsc_library),
        PetscErrorCode,
        (
            DM,
            IS,
            Vec,
            Vec,
            Vec,
            Ptr{Ptr{$PetscScalar}},
            Ptr{Ptr{$PetscScalar}},
            Ptr{Ptr{$PetscScalar}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function DMPlexGetFaceFields(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
)
    @chk ccall(
        (:DMPlexGetFaceFields, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            $PetscInt,
            Vec,
            Vec,
            Vec,
            Vec,
            Vec,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscScalar}},
            Ptr{Ptr{$PetscScalar}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
    )
end

@for_petsc function DMPlexRestoreFaceFields(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
)
    @chk ccall(
        (:DMPlexRestoreFaceFields, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            $PetscInt,
            Vec,
            Vec,
            Vec,
            Vec,
            Vec,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscScalar}},
            Ptr{Ptr{$PetscScalar}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
    )
end

@for_petsc function DMPlexGetScale(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexGetScale, $petsc_library),
        PetscErrorCode,
        (DM, PetscUnit, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexSetScale(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexSetScale, $petsc_library),
        PetscErrorCode,
        (DM, PetscUnit, $PetscReal),
        arg1,
        arg2,
        arg3,
    )
end

mutable struct JacActionCtx
    dm::DM
    u::Vec
    J::Mat
    user::Ptr{Cvoid}
    JacActionCtx() = new()
end

@for_petsc function DMPlexSetMaxProjectionHeight(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexSetMaxProjectionHeight, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetMaxProjectionHeight(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetMaxProjectionHeight, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetActivePoint(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetActivePoint, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexSetActivePoint(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexSetActivePoint, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexProjectFieldLocal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexProjectFieldLocal, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Ptr{Cvoid}}, InsertMode, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexComputeL2DiffLocal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMPlexComputeL2DiffLocal, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscReal,
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
            Vec,
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMPlexComputeL2FieldDiff(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMPlexComputeL2FieldDiff, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscReal,
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
            Vec,
            Ptr{$PetscReal},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMPlexComputeL2DiffVec(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMPlexComputeL2DiffVec, $petsc_library),
        PetscErrorCode,
        (DM, $PetscReal, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMPlexComputeCellwiseIntegralFEM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexComputeCellwiseIntegralFEM, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexComputeIntegralFEM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexComputeIntegralFEM, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{$PetscScalar}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexComputeBdIntegral(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:DMPlexComputeBdIntegral, $petsc_library),
        PetscErrorCode,
        (
            DM,
            Vec,
            DMLabel,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Cvoid},
            Ptr{$PetscScalar},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function DMPlexComputeInterpolatorNested(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexComputeInterpolatorNested, $petsc_library),
        PetscErrorCode,
        (DM, DM, PetscBool, Mat, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexComputeInterpolatorGeneral(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexComputeInterpolatorGeneral, $petsc_library),
        PetscErrorCode,
        (DM, DM, Mat, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexComputeGradientClementInterpolant(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMPlexComputeGradientClementInterpolant, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexComputeInjectorFEM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexComputeInjectorFEM, $petsc_library),
        PetscErrorCode,
        (DM, DM, Ptr{VecScatter}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexComputeMassMatrixNested(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexComputeMassMatrixNested, $petsc_library),
        PetscErrorCode,
        (DM, DM, Mat, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexComputeMassMatrixGeneral(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexComputeMassMatrixGeneral, $petsc_library),
        PetscErrorCode,
        (DM, DM, Mat, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexCreateRigidBody(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexCreateRigidBody, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{MatNullSpace}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexCreateRigidBodies(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMPlexCreateRigidBodies, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            DMLabel,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{MatNullSpace},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMPlexSetSNESLocalFEM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexSetSNESLocalFEM, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexSNESComputeBoundaryFEM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMPlexSNESComputeBoundaryFEM, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexSNESComputeResidualFEM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexSNESComputeResidualFEM, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexSNESComputeJacobianFEM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexSNESComputeJacobianFEM, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Mat, Mat, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexComputeJacobianAction(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:DMPlexComputeJacobianAction, $petsc_library),
        PetscErrorCode,
        (DM, IS, $PetscReal, $PetscReal, Vec, Vec, Vec, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function DMPlexComputeBdResidualSingle(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:DMPlexComputeBdResidualSingle, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscReal,
            DMLabel,
            $PetscInt,
            Ptr{$PetscInt},
            $PetscInt,
            Vec,
            Vec,
            Vec,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function DMPlexComputeBdJacobianSingle(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
)
    @chk ccall(
        (:DMPlexComputeBdJacobianSingle, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscReal,
            DMLabel,
            $PetscInt,
            Ptr{$PetscInt},
            $PetscInt,
            Vec,
            Vec,
            $PetscReal,
            Mat,
            Mat,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
    )
end

@for_petsc function DMPlexTSComputeBoundary(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexTSComputeBoundary, $petsc_library),
        PetscErrorCode,
        (DM, $PetscReal, Vec, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexTSComputeRHSFunctionFVM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexTSComputeRHSFunctionFVM, $petsc_library),
        PetscErrorCode,
        (DM, $PetscReal, Vec, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexTSComputeIFunctionFEM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMPlexTSComputeIFunctionFEM, $petsc_library),
        PetscErrorCode,
        (DM, $PetscReal, Vec, Vec, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMPlexTSComputeIJacobianFEM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:DMPlexTSComputeIJacobianFEM, $petsc_library),
        PetscErrorCode,
        (DM, $PetscReal, Vec, Vec, $PetscReal, Mat, Mat, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function DMPlexComputeRHSFunctionFVM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexComputeRHSFunctionFVM, $petsc_library),
        PetscErrorCode,
        (DM, $PetscReal, Vec, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMPlexReconstructGradientsFVM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMPlexReconstructGradientsFVM, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexGetAnchors(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexGetAnchors, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscSection}, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexSetAnchors(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexSetAnchors, $petsc_library),
        PetscErrorCode,
        (DM, PetscSection, IS),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexSetReferenceTree(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexSetReferenceTree, $petsc_library),
        PetscErrorCode,
        (DM, DM),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetReferenceTree(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetReferenceTree, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DM}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexReferenceTreeGetChildSymmetry(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:DMPlexReferenceTreeGetChildSymmetry, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function DMPlexCreateDefaultReferenceTree(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexCreateDefaultReferenceTree, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, PetscBool, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexSetTree(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMPlexSetTree, $petsc_library),
        PetscErrorCode,
        (DM, PetscSection, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexGetTree(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMPlexGetTree, $petsc_library),
        PetscErrorCode,
        (
            DM,
            Ptr{PetscSection},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscInt}},
            Ptr{PetscSection},
            Ptr{Ptr{$PetscInt}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMPlexGetTreeParent(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexGetTreeParent, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexGetTreeChildren(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexGetTreeChildren, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexTreeRefineCell(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexTreeRefineCell, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{DM}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexComputeInjectorReferenceTree(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:DMPlexComputeInjectorReferenceTree, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexTransferVecTree(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
)
    @chk ccall(
        (:DMPlexTransferVecTree, $petsc_library),
        PetscErrorCode,
        (
            DM,
            Vec,
            DM,
            Vec,
            PetscSF,
            PetscSF,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            PetscBool,
            $PetscReal,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
    )
end

@for_petsc function DMPlexMonitorThroughput(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexMonitorThroughput, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexCreateGlobalToNaturalSF(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexCreateGlobalToNaturalSF, $petsc_library),
        PetscErrorCode,
        (DM, PetscSection, PetscSF, Ptr{PetscSF}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexSetGlobalToNaturalSF(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexSetGlobalToNaturalSF, $petsc_library),
        PetscErrorCode,
        (DM, PetscSF),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetGlobalToNaturalSF(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetGlobalToNaturalSF, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscSF}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGlobalToNaturalBegin(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMPlexGlobalToNaturalBegin, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexGlobalToNaturalEnd(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexGlobalToNaturalEnd, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexNaturalToGlobalBegin(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMPlexNaturalToGlobalBegin, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexNaturalToGlobalEnd(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexNaturalToGlobalEnd, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexAdapt(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMPlexAdapt, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Cchar}, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexSnapToGeomModel(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexSnapToGeomModel, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscScalar}, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexGlobalToLocalBasis(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGlobalToLocalBasis, $petsc_library),
        PetscErrorCode,
        (DM, Vec),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexLocalToGlobalBasis(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexLocalToGlobalBasis, $petsc_library),
        PetscErrorCode,
        (DM, Vec),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexCreateBasisRotation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexCreateBasisRotation, $petsc_library),
        PetscErrorCode,
        (DM, $PetscReal, $PetscReal, $PetscReal),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexCellRefinerCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexCellRefinerCreate, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DMPlexCellRefiner}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexCellRefinerSetUp(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMPlexCellRefinerSetUp, $petsc_library),
        PetscErrorCode,
        (DMPlexCellRefiner,),
        arg1,
    )
end

@for_petsc function DMPlexCellRefinerDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMPlexCellRefinerDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{DMPlexCellRefiner},),
        arg1,
    )
end

@for_petsc function DMPlexCellRefinerRefine(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:DMPlexCellRefinerRefine, $petsc_library),
        PetscErrorCode,
        (
            DMPlexCellRefiner,
            DMPolytopeType,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{Ptr{DMPolytopeType}},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscInt}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function DMPlexCellRefinerGetAffineTransforms(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMPlexCellRefinerGetAffineTransforms, $petsc_library),
        PetscErrorCode,
        (
            DMPlexCellRefiner,
            DMPolytopeType,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscReal}},
            Ptr{Ptr{$PetscReal}},
            Ptr{Ptr{$PetscReal}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMPlexCellRefinerGetAffineFaceTransforms(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMPlexCellRefinerGetAffineFaceTransforms, $petsc_library),
        PetscErrorCode,
        (
            DMPlexCellRefiner,
            DMPolytopeType,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscReal}},
            Ptr{Ptr{$PetscReal}},
            Ptr{Ptr{$PetscReal}},
            Ptr{Ptr{$PetscReal}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMPlexRefineUniform(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexRefineUniform, $petsc_library),
        PetscErrorCode,
        (DM, DMPlexCellRefiner, Ptr{DM}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMRedundantCreate(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMRedundantCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, PetscMPIInt, $PetscInt, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMRedundantSetSize(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMRedundantSetSize, $petsc_library),
        PetscErrorCode,
        (DM, PetscMPIInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMRedundantGetSize(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMRedundantGetSize, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscMPIInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMShellCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMShellCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{DM}),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellSetContext(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMShellSetContext, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellGetContext(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMShellGetContext, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellSetMatrix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMShellSetMatrix, $petsc_library),
        PetscErrorCode,
        (DM, Mat),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellSetGlobalVector(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMShellSetGlobalVector, $petsc_library),
        PetscErrorCode,
        (DM, Vec),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellGetGlobalVector(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMShellGetGlobalVector, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellSetLocalVector(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMShellSetLocalVector, $petsc_library),
        PetscErrorCode,
        (DM, Vec),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellSetCreateGlobalVector(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMShellSetCreateGlobalVector, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellSetCreateLocalVector(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMShellSetCreateLocalVector, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellSetGlobalToLocal(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMShellSetGlobalToLocal, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMShellSetGlobalToLocalVecScatter(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:DMShellSetGlobalToLocalVecScatter, $petsc_library),
        PetscErrorCode,
        (DM, VecScatter),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellSetLocalToGlobal(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMShellSetLocalToGlobal, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMShellSetLocalToGlobalVecScatter(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:DMShellSetLocalToGlobalVecScatter, $petsc_library),
        PetscErrorCode,
        (DM, VecScatter),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellSetLocalToLocal(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMShellSetLocalToLocal, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMShellSetLocalToLocalVecScatter(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:DMShellSetLocalToLocalVecScatter, $petsc_library),
        PetscErrorCode,
        (DM, VecScatter),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellSetCreateMatrix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMShellSetCreateMatrix, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellSetCoarsen(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMShellSetCoarsen, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellGetCoarsen(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMShellGetCoarsen, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellSetRefine(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMShellSetRefine, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellGetRefine(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMShellGetRefine, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellSetCreateInterpolation(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMShellSetCreateInterpolation, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellGetCreateInterpolation(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMShellGetCreateInterpolation, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellSetCreateRestriction(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMShellSetCreateRestriction, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellGetCreateRestriction(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMShellGetCreateRestriction, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellSetCreateInjection(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMShellSetCreateInjection, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellGetCreateInjection(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMShellGetCreateInjection, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellSetCreateFieldDecomposition(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:DMShellSetCreateFieldDecomposition, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellSetCreateDomainDecomposition(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:DMShellSetCreateDomainDecomposition, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellSetCreateDomainDecompositionScatters(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:DMShellSetCreateDomainDecompositionScatters, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellSetCreateSubDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMShellSetCreateSubDM, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMShellGetCreateSubDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMShellGetCreateSubDM, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
    )
end

@for_petsc function DMGlobalToLocalBeginDefaultShell(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMGlobalToLocalBeginDefaultShell, $petsc_library),
        PetscErrorCode,
        (DM, Vec, InsertMode, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMGlobalToLocalEndDefaultShell(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMGlobalToLocalEndDefaultShell, $petsc_library),
        PetscErrorCode,
        (DM, Vec, InsertMode, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMLocalToGlobalBeginDefaultShell(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMLocalToGlobalBeginDefaultShell, $petsc_library),
        PetscErrorCode,
        (DM, Vec, InsertMode, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMLocalToGlobalEndDefaultShell(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMLocalToGlobalEndDefaultShell, $petsc_library),
        PetscErrorCode,
        (DM, Vec, InsertMode, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMLocalToLocalBeginDefaultShell(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMLocalToLocalBeginDefaultShell, $petsc_library),
        PetscErrorCode,
        (DM, Vec, InsertMode, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMLocalToLocalEndDefaultShell(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMLocalToLocalEndDefaultShell, $petsc_library),
        PetscErrorCode,
        (DM, Vec, InsertMode, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMSlicedCreate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:DMSlicedCreate, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{DM},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function DMSlicedSetPreallocation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMSlicedSetPreallocation, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMSlicedSetBlockFills(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMSlicedSetBlockFills, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMSlicedSetGhosts(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMSlicedSetGhosts, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
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

@for_petsc function DMSwarmCreateGlobalVectorFromField(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMSwarmCreateGlobalVectorFromField, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMSwarmDestroyGlobalVectorFromField(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMSwarmDestroyGlobalVectorFromField, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMSwarmCreateLocalVectorFromField(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMSwarmCreateLocalVectorFromField, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMSwarmDestroyLocalVectorFromField(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMSwarmDestroyLocalVectorFromField, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMSwarmInitializeFieldRegister(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMSwarmInitializeFieldRegister, $petsc_library),
        PetscErrorCode,
        (DM,),
        arg1,
    )
end

@for_petsc function DMSwarmFinalizeFieldRegister(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMSwarmFinalizeFieldRegister, $petsc_library),
        PetscErrorCode,
        (DM,),
        arg1,
    )
end

@for_petsc function DMSwarmSetLocalSizes(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMSwarmSetLocalSizes, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMSwarmRegisterPetscDatatypeField(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMSwarmRegisterPetscDatatypeField, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, $PetscInt, PetscDataType),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMSwarmRegisterUserStructField(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMSwarmRegisterUserStructField, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMSwarmRegisterUserDatatypeField(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMSwarmRegisterUserDatatypeField, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, Csize_t, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMSwarmGetField(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMSwarmGetField, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, Ptr{$PetscInt}, Ptr{PetscDataType}, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMSwarmRestoreField(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMSwarmRestoreField, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, Ptr{$PetscInt}, Ptr{PetscDataType}, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMSwarmVectorDefineField(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSwarmVectorDefineField, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSwarmAddPoint(::$UnionPetscLib, arg1)
    @chk ccall((:DMSwarmAddPoint, $petsc_library), PetscErrorCode, (DM,), arg1)
end

@for_petsc function DMSwarmAddNPoints(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSwarmAddNPoints, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function DMSwarmRemovePoint(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMSwarmRemovePoint, $petsc_library),
        PetscErrorCode,
        (DM,),
        arg1,
    )
end

@for_petsc function DMSwarmRemovePointAtIndex(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSwarmRemovePointAtIndex, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function DMSwarmCopyPoint(::$UnionPetscLib, dm, arg2, arg3)
    @chk ccall(
        (:DMSwarmCopyPoint, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt),
        dm,
        arg2,
        arg3,
    )
end

@for_petsc function DMSwarmGetLocalSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSwarmGetLocalSize, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSwarmGetSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSwarmGetSize, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSwarmMigrate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSwarmMigrate, $petsc_library),
        PetscErrorCode,
        (DM, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function DMSwarmCollectViewCreate(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMSwarmCollectViewCreate, $petsc_library),
        PetscErrorCode,
        (DM,),
        arg1,
    )
end

@for_petsc function DMSwarmCollectViewDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMSwarmCollectViewDestroy, $petsc_library),
        PetscErrorCode,
        (DM,),
        arg1,
    )
end

@for_petsc function DMSwarmSetCellDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSwarmSetCellDM, $petsc_library),
        PetscErrorCode,
        (DM, DM),
        arg1,
        arg2,
    )
end

@for_petsc function DMSwarmGetCellDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSwarmGetCellDM, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DM}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSwarmSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSwarmSetType, $petsc_library),
        PetscErrorCode,
        (DM, DMSwarmType),
        arg1,
        arg2,
    )
end

@for_petsc function DMSwarmSetPointsUniformCoordinates(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMSwarmSetPointsUniformCoordinates, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscInt}, InsertMode),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMSwarmSetPointCoordinates(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMSwarmSetPointCoordinates, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscReal}, PetscBool, InsertMode),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMSwarmInsertPointsUsingCellDM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMSwarmInsertPointsUsingCellDM, $petsc_library),
        PetscErrorCode,
        (DM, DMSwarmPICLayoutType, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMSwarmSetPointCoordinatesCellwise(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMSwarmSetPointCoordinatesCellwise, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMSwarmSetPointCoordinatesRandom(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:DMSwarmSetPointCoordinatesRandom, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function DMSwarmViewFieldsXDMF(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMSwarmViewFieldsXDMF, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, $PetscInt, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMSwarmViewXDMF(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSwarmViewXDMF, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSwarmSortGetAccess(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMSwarmSortGetAccess, $petsc_library),
        PetscErrorCode,
        (DM,),
        arg1,
    )
end

@for_petsc function DMSwarmSortRestoreAccess(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMSwarmSortRestoreAccess, $petsc_library),
        PetscErrorCode,
        (DM,),
        arg1,
    )
end

@for_petsc function DMSwarmSortGetPointsPerCell(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMSwarmSortGetPointsPerCell, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMSwarmSortGetNumberOfPointsPerCell(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMSwarmSortGetNumberOfPointsPerCell, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMSwarmSortGetIsValid(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSwarmSortGetIsValid, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSwarmSortGetSizes(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMSwarmSortGetSizes, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMSwarmProjectFields(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMSwarmProjectFields, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{Ptr{Cchar}}, Ptr{Ptr{Vec}}, PetscBool),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMSwarmCreateMassMatrixSquare(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMSwarmCreateMassMatrixSquare, $petsc_library),
        PetscErrorCode,
        (DM, DM, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMCreate_Product(::$UnionPetscLib, arg1)
    @chk ccall((:DMCreate_Product, $petsc_library), PetscErrorCode, (DM,), arg1)
end

@for_petsc function DMProductGetDM(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMProductGetDM, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{DM}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMProductSetDimensionIndex(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMProductSetDimensionIndex, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMProductSetDM(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMProductSetDM, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, DM),
        arg1,
        arg2,
        arg3,
    )
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

@enum DMStagStencilType::UInt32 begin
    DMSTAG_STENCIL_NONE = 0
    DMSTAG_STENCIL_STAR = 1
    DMSTAG_STENCIL_BOX = 2
end

@for_petsc function DMCreate_Stag(::$UnionPetscLib, arg1)
    @chk ccall((:DMCreate_Stag, $petsc_library), PetscErrorCode, (DM,), arg1)
end

@for_petsc function DMStagCreate1d(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:DMStagCreate1d, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            DMBoundaryType,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            DMStagStencilType,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{DM},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function DMStagCreate2d(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
    arg14,
    arg15,
)
    @chk ccall(
        (:DMStagCreate2d, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            DMBoundaryType,
            DMBoundaryType,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            DMStagStencilType,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{DM},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
        arg14,
        arg15,
    )
end

@for_petsc function DMStagCreate3d(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
    arg14,
    arg15,
    arg16,
    arg17,
    arg18,
    arg19,
    arg20,
)
    @chk ccall(
        (:DMStagCreate3d, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            DMBoundaryType,
            DMBoundaryType,
            DMBoundaryType,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            DMStagStencilType,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{DM},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
        arg14,
        arg15,
        arg16,
        arg17,
        arg18,
        arg19,
        arg20,
    )
end

@for_petsc function DMStagCreateCompatibleDMStag(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMStagCreateCompatibleDMStag, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMStagGetBoundaryTypes(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMStagGetBoundaryTypes, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DMBoundaryType}, Ptr{DMBoundaryType}, Ptr{DMBoundaryType}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMStagGetCorners(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
)
    @chk ccall(
        (:DMStagGetCorners, $petsc_library),
        PetscErrorCode,
        (
            DM,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
    )
end

@for_petsc function DMStagGetDOF(::$UnionPetscLib, arg1, arg2, arg3, arg4, arg5)
    @chk ccall(
        (:DMStagGetDOF, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMStagGetEntries(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMStagGetEntries, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMStagGetEntriesPerElement(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMStagGetEntriesPerElement, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMStagGetGhostCorners(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMStagGetGhostCorners, $petsc_library),
        PetscErrorCode,
        (
            DM,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMStagGetGlobalSizes(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMStagGetGlobalSizes, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMStagGetIsFirstRank(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMStagGetIsFirstRank, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscBool}, Ptr{PetscBool}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMStagGetIsLastRank(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMStagGetIsLastRank, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscBool}, Ptr{PetscBool}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMStagGetLocalSizes(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMStagGetLocalSizes, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMStagGetLocationDOF(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMStagGetLocationDOF, $petsc_library),
        PetscErrorCode,
        (DM, DMStagStencilLocation, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMStagGetLocationSlot(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMStagGetLocationSlot, $petsc_library),
        PetscErrorCode,
        (DM, DMStagStencilLocation, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMStagGetNumRanks(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMStagGetNumRanks, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMStagGetOwnershipRanges(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMStagGetOwnershipRanges, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMStagGetProductCoordinateArrays(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMStagGetProductCoordinateArrays, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMStagGetProductCoordinateArraysRead(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMStagGetProductCoordinateArraysRead, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMStagGetProductCoordinateLocationSlot(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMStagGetProductCoordinateLocationSlot, $petsc_library),
        PetscErrorCode,
        (DM, DMStagStencilLocation, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMStagGetStencilType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMStagGetStencilType, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DMStagStencilType}),
        arg1,
        arg2,
    )
end

@for_petsc function DMStagGetStencilWidth(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMStagGetStencilWidth, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMStagMigrateVec(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMStagMigrateVec, $petsc_library),
        PetscErrorCode,
        (DM, Vec, DM, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMStagPopulateLocalToGlobalInjective(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMStagPopulateLocalToGlobalInjective, $petsc_library),
        PetscErrorCode,
        (DM,),
        arg1,
    )
end

@for_petsc function DMStagRestoreProductCoordinateArrays(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMStagRestoreProductCoordinateArrays, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMStagRestoreProductCoordinateArraysRead(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMStagRestoreProductCoordinateArraysRead, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMStagSetBoundaryTypes(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMStagSetBoundaryTypes, $petsc_library),
        PetscErrorCode,
        (DM, DMBoundaryType, DMBoundaryType, DMBoundaryType),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMStagSetCoordinateDMType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMStagSetCoordinateDMType, $petsc_library),
        PetscErrorCode,
        (DM, DMType),
        arg1,
        arg2,
    )
end

@for_petsc function DMStagSetDOF(::$UnionPetscLib, arg1, arg2, arg3, arg4, arg5)
    @chk ccall(
        (:DMStagSetDOF, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMStagSetGlobalSizes(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMStagSetGlobalSizes, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMStagSetNumRanks(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMStagSetNumRanks, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMStagSetOwnershipRanges(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMStagSetOwnershipRanges, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMStagSetStencilType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMStagSetStencilType, $petsc_library),
        PetscErrorCode,
        (DM, DMStagStencilType),
        arg1,
        arg2,
    )
end

@for_petsc function DMStagSetStencilWidth(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMStagSetStencilWidth, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function DMStagSetUniformCoordinates(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMStagSetUniformCoordinates, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMStagSetUniformCoordinatesExplicit(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMStagSetUniformCoordinatesExplicit, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMStagSetUniformCoordinatesProduct(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:DMStagSetUniformCoordinatesProduct, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function DMStagVecGetArray(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMStagVecGetArray, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMStagVecGetArrayRead(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMStagVecGetArrayRead, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMStagVecRestoreArray(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMStagVecRestoreArray, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMStagVecRestoreArrayRead(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMStagVecRestoreArrayRead, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMStagVecSplitToDMDA(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMStagVecSplitToDMDA, $petsc_library),
        PetscErrorCode,
        (DM, Vec, DMStagStencilLocation, $PetscInt, Ptr{DM}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMStagGet1dCoordinateArraysDOFRead(
    ::$UnionPetscLib,
    dm,
    ax,
    ay,
    az,
)
    @chk ccall(
        (:DMStagGet1dCoordinateArraysDOFRead, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        dm,
        ax,
        ay,
        az,
    )
end

@for_petsc function DMStagGet1dCoordinateLocationSlot(
    ::$UnionPetscLib,
    dm,
    loc,
    s,
)
    @chk ccall(
        (:DMStagGet1dCoordinateLocationSlot, $petsc_library),
        PetscErrorCode,
        (DM, DMStagStencilLocation, Ptr{$PetscInt}),
        dm,
        loc,
        s,
    )
end

@for_petsc function DMStagGetGhostType(::$UnionPetscLib, dm, s)
    @chk ccall(
        (:DMStagGetGhostType, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DMStagStencilType}),
        dm,
        s,
    )
end

@for_petsc function DMStagRestore1dCoordinateArraysDOFRead(
    ::$UnionPetscLib,
    dm,
    ax,
    ay,
    az,
)
    @chk ccall(
        (:DMStagRestore1dCoordinateArraysDOFRead, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        dm,
        ax,
        ay,
        az,
    )
end

@for_petsc function DMStagSetGhostType(::$UnionPetscLib, dm, s)
    @chk ccall(
        (:DMStagSetGhostType, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DMStagStencilType}),
        dm,
        s,
    )
end

@for_petsc function DMStagVecGetArrayDOF(::$UnionPetscLib, dm, v, a)
    @chk ccall(
        (:DMStagVecGetArrayDOF, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Cvoid}),
        dm,
        v,
        a,
    )
end

@for_petsc function DMStagVecGetArrayDOFRead(::$UnionPetscLib, dm, v, a)
    @chk ccall(
        (:DMStagVecGetArrayDOFRead, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Cvoid}),
        dm,
        v,
        a,
    )
end

@for_petsc function DMStagVecRestoreArrayDOF(::$UnionPetscLib, dm, v, a)
    @chk ccall(
        (:DMStagVecRestoreArrayDOF, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Cvoid}),
        dm,
        v,
        a,
    )
end

@for_petsc function DMStagVecRestoreArrayDOFRead(::$UnionPetscLib, dm, v, a)
    @chk ccall(
        (:DMStagVecRestoreArrayDOFRead, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Cvoid}),
        dm,
        v,
        a,
    )
end

@for_petsc function PetscWeakFormCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscWeakFormCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscWeakForm}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscWeakFormDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscWeakFormDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscWeakForm},),
        arg1,
    )
end

@for_petsc function PetscWeakFormView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscWeakFormView, $petsc_library),
        PetscErrorCode,
        (PetscWeakForm, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscWeakFormGetNumFields(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscWeakFormGetNumFields, $petsc_library),
        PetscErrorCode,
        (PetscWeakForm, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscWeakFormSetNumFields(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscWeakFormSetNumFields, $petsc_library),
        PetscErrorCode,
        (PetscWeakForm, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscWeakFormGetObjective(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscWeakFormGetObjective, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscWeakFormAddObjective(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscWeakFormAddObjective, $petsc_library),
        PetscErrorCode,
        (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscWeakFormSetObjective(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscWeakFormSetObjective, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Cvoid}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscWeakFormGetIndexObjective(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscWeakFormGetIndexObjective, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Cvoid}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscWeakFormSetIndexObjective(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscWeakFormSetIndexObjective, $petsc_library),
        PetscErrorCode,
        (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscWeakFormGetResidual(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscWeakFormGetResidual, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function PetscWeakFormAddResidual(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscWeakFormAddResidual, $petsc_library),
        PetscErrorCode,
        (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscWeakFormSetResidual(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscWeakFormSetResidual, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Cvoid}},
            $PetscInt,
            Ptr{Ptr{Cvoid}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function PetscWeakFormSetIndexResidual(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscWeakFormSetIndexResidual, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Cvoid},
            $PetscInt,
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function PetscWeakFormHasJacobian(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscWeakFormHasJacobian, $petsc_library),
        PetscErrorCode,
        (PetscWeakForm, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscWeakFormGetJacobian(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
)
    @chk ccall(
        (:PetscWeakFormGetJacobian, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
    )
end

@for_petsc function PetscWeakFormAddJacobian(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:PetscWeakFormAddJacobian, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function PetscWeakFormSetJacobian(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
)
    @chk ccall(
        (:PetscWeakFormSetJacobian, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Cvoid}},
            $PetscInt,
            Ptr{Ptr{Cvoid}},
            $PetscInt,
            Ptr{Ptr{Cvoid}},
            $PetscInt,
            Ptr{Ptr{Cvoid}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
    )
end

@for_petsc function PetscWeakFormSetIndexJacobian(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
)
    @chk ccall(
        (:PetscWeakFormSetIndexJacobian, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Cvoid},
            $PetscInt,
            Ptr{Cvoid},
            $PetscInt,
            Ptr{Cvoid},
            $PetscInt,
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
    )
end

@for_petsc function PetscWeakFormHasJacobianPreconditioner(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscWeakFormHasJacobianPreconditioner, $petsc_library),
        PetscErrorCode,
        (PetscWeakForm, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscWeakFormGetJacobianPreconditioner(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
)
    @chk ccall(
        (:PetscWeakFormGetJacobianPreconditioner, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
    )
end

@for_petsc function PetscWeakFormAddJacobianPreconditioner(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:PetscWeakFormAddJacobianPreconditioner, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function PetscWeakFormSetJacobianPreconditioner(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
)
    @chk ccall(
        (:PetscWeakFormSetJacobianPreconditioner, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Cvoid}},
            $PetscInt,
            Ptr{Ptr{Cvoid}},
            $PetscInt,
            Ptr{Ptr{Cvoid}},
            $PetscInt,
            Ptr{Ptr{Cvoid}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
    )
end

@for_petsc function PetscWeakFormSetIndexJacobianPreconditioner(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
)
    @chk ccall(
        (:PetscWeakFormSetIndexJacobianPreconditioner, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Cvoid},
            $PetscInt,
            Ptr{Cvoid},
            $PetscInt,
            Ptr{Cvoid},
            $PetscInt,
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
    )
end

@for_petsc function PetscWeakFormHasDynamicJacobian(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscWeakFormHasDynamicJacobian, $petsc_library),
        PetscErrorCode,
        (PetscWeakForm, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscWeakFormGetDynamicJacobian(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
)
    @chk ccall(
        (:PetscWeakFormGetDynamicJacobian, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
    )
end

@for_petsc function PetscWeakFormAddDynamicJacobian(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:PetscWeakFormAddDynamicJacobian, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function PetscWeakFormSetDynamicJacobian(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
)
    @chk ccall(
        (:PetscWeakFormSetDynamicJacobian, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Cvoid}},
            $PetscInt,
            Ptr{Ptr{Cvoid}},
            $PetscInt,
            Ptr{Ptr{Cvoid}},
            $PetscInt,
            Ptr{Ptr{Cvoid}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
    )
end

@for_petsc function PetscWeakFormSetIndexDynamicJacobian(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
)
    @chk ccall(
        (:PetscWeakFormSetIndexDynamicJacobian, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Cvoid},
            $PetscInt,
            Ptr{Cvoid},
            $PetscInt,
            Ptr{Cvoid},
            $PetscInt,
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
    )
end

@for_petsc function PetscWeakFormGetBdResidual(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscWeakFormGetBdResidual, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function PetscWeakFormAddBdResidual(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscWeakFormAddBdResidual, $petsc_library),
        PetscErrorCode,
        (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscWeakFormSetBdResidual(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscWeakFormSetBdResidual, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Cvoid}},
            $PetscInt,
            Ptr{Ptr{Cvoid}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function PetscWeakFormSetIndexBdResidual(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:PetscWeakFormSetIndexBdResidual, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Cvoid},
            $PetscInt,
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function PetscWeakFormHasBdJacobian(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscWeakFormHasBdJacobian, $petsc_library),
        PetscErrorCode,
        (PetscWeakForm, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscWeakFormGetBdJacobian(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
)
    @chk ccall(
        (:PetscWeakFormGetBdJacobian, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
    )
end

@for_petsc function PetscWeakFormAddBdJacobian(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:PetscWeakFormAddBdJacobian, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function PetscWeakFormSetBdJacobian(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
)
    @chk ccall(
        (:PetscWeakFormSetBdJacobian, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Cvoid}},
            $PetscInt,
            Ptr{Ptr{Cvoid}},
            $PetscInt,
            Ptr{Ptr{Cvoid}},
            $PetscInt,
            Ptr{Ptr{Cvoid}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
    )
end

@for_petsc function PetscWeakFormSetIndexBdJacobian(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
)
    @chk ccall(
        (:PetscWeakFormSetIndexBdJacobian, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Cvoid},
            $PetscInt,
            Ptr{Cvoid},
            $PetscInt,
            Ptr{Cvoid},
            $PetscInt,
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
    )
end

@for_petsc function PetscWeakFormHasBdJacobianPreconditioner(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscWeakFormHasBdJacobianPreconditioner, $petsc_library),
        PetscErrorCode,
        (PetscWeakForm, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscWeakFormGetBdJacobianPreconditioner(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
)
    @chk ccall(
        (:PetscWeakFormGetBdJacobianPreconditioner, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
    )
end

@for_petsc function PetscWeakFormAddBdJacobianPreconditioner(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:PetscWeakFormAddBdJacobianPreconditioner, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function PetscWeakFormSetBdJacobianPreconditioner(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
)
    @chk ccall(
        (:PetscWeakFormSetBdJacobianPreconditioner, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Cvoid}},
            $PetscInt,
            Ptr{Ptr{Cvoid}},
            $PetscInt,
            Ptr{Ptr{Cvoid}},
            $PetscInt,
            Ptr{Ptr{Cvoid}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
    )
end

@for_petsc function PetscWeakFormSetIndexBdJacobianPreconditioner(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
)
    @chk ccall(
        (:PetscWeakFormSetIndexBdJacobianPreconditioner, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Cvoid},
            $PetscInt,
            Ptr{Cvoid},
            $PetscInt,
            Ptr{Cvoid},
            $PetscInt,
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
    )
end

@for_petsc function PetscWeakFormGetRiemannSolver(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscWeakFormGetRiemannSolver, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{Ptr{Cvoid}}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscWeakFormSetRiemannSolver(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscWeakFormSetRiemannSolver, $petsc_library),
        PetscErrorCode,
        (
            PetscWeakForm,
            DMLabel,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Cvoid}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscWeakFormSetIndexRiemannSolver(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscWeakFormSetIndexRiemannSolver, $petsc_library),
        PetscErrorCode,
        (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscDSInitializePackage(::$UnionPetscLib)
    @chk ccall((:PetscDSInitializePackage, $petsc_library), PetscErrorCode, ())
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

@for_petsc function PetscDSCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscDS}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDSDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscDS},),
        arg1,
    )
end

@for_petsc function PetscDSSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSSetType, $petsc_library),
        PetscErrorCode,
        (PetscDS, PetscDSType),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSGetType, $petsc_library),
        PetscErrorCode,
        (PetscDS, Ptr{PetscDSType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSSetUp(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDSSetUp, $petsc_library),
        PetscErrorCode,
        (PetscDS,),
        arg1,
    )
end

@for_petsc function PetscDSSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDSSetFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscDS,),
        arg1,
    )
end

@for_petsc function PetscDSViewFromOptions(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDSViewFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscDS, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDSView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSView, $petsc_library),
        PetscErrorCode,
        (PetscDS, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSRegisterDestroy(::$UnionPetscLib)
    @chk ccall((:PetscDSRegisterDestroy, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscDSGetHeightSubspace(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDSGetHeightSubspace, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{PetscDS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDSGetSpatialDimension(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSGetSpatialDimension, $petsc_library),
        PetscErrorCode,
        (PetscDS, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSGetCoordinateDimension(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSGetCoordinateDimension, $petsc_library),
        PetscErrorCode,
        (PetscDS, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSSetCoordinateDimension(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSSetCoordinateDimension, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSGetHybrid(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSGetHybrid, $petsc_library),
        PetscErrorCode,
        (PetscDS, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSSetHybrid(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSSetHybrid, $petsc_library),
        PetscErrorCode,
        (PetscDS, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSGetNumFields(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSGetNumFields, $petsc_library),
        PetscErrorCode,
        (PetscDS, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSGetTotalDimension(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSGetTotalDimension, $petsc_library),
        PetscErrorCode,
        (PetscDS, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSGetTotalComponents(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSGetTotalComponents, $petsc_library),
        PetscErrorCode,
        (PetscDS, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSGetFieldIndex(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDSGetFieldIndex, $petsc_library),
        PetscErrorCode,
        (PetscDS, PetscObject, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDSGetFieldSize(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDSGetFieldSize, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDSGetFieldOffset(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDSGetFieldOffset, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDSGetDimensions(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSGetDimensions, $petsc_library),
        PetscErrorCode,
        (PetscDS, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSGetComponents(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSGetComponents, $petsc_library),
        PetscErrorCode,
        (PetscDS, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSGetComponentOffset(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscDSGetComponentOffset, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDSGetComponentOffsets(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSGetComponentOffsets, $petsc_library),
        PetscErrorCode,
        (PetscDS, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSGetComponentDerivativeOffsets(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDSGetComponentDerivativeOffsets, $petsc_library),
        PetscErrorCode,
        (PetscDS, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSGetWeakForm(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSGetWeakForm, $petsc_library),
        PetscErrorCode,
        (PetscDS, Ptr{PetscWeakForm}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSSetWeakForm(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSSetWeakForm, $petsc_library),
        PetscErrorCode,
        (PetscDS, PetscWeakForm),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSGetDiscretization(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDSGetDiscretization, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{PetscObject}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDSSetDiscretization(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDSSetDiscretization, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, PetscObject),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDSAddDiscretization(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSAddDiscretization, $petsc_library),
        PetscErrorCode,
        (PetscDS, PetscObject),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSGetQuadrature(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSGetQuadrature, $petsc_library),
        PetscErrorCode,
        (PetscDS, Ptr{PetscQuadrature}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSGetImplicit(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDSGetImplicit, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDSSetImplicit(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDSSetImplicit, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, PetscBool),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDSGetJetDegree(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDSGetJetDegree, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDSSetJetDegree(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDSSetJetDegree, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDSGetConstants(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDSGetConstants, $petsc_library),
        PetscErrorCode,
        (PetscDS, Ptr{$PetscInt}, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDSSetConstants(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDSSetConstants, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDSGetObjective(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDSGetObjective, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDSSetObjective(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDSSetObjective, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDSGetResidual(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscDSGetResidual, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDSSetResidual(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscDSSetResidual, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDSHasJacobian(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSHasJacobian, $petsc_library),
        PetscErrorCode,
        (PetscDS, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSGetJacobian(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscDSGetJacobian, $petsc_library),
        PetscErrorCode,
        (
            PetscDS,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscDSSetJacobian(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscDSSetJacobian, $petsc_library),
        PetscErrorCode,
        (
            PetscDS,
            $PetscInt,
            $PetscInt,
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscDSUseJacobianPreconditioner(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDSUseJacobianPreconditioner, $petsc_library),
        PetscErrorCode,
        (PetscDS, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSHasJacobianPreconditioner(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDSHasJacobianPreconditioner, $petsc_library),
        PetscErrorCode,
        (PetscDS, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSGetJacobianPreconditioner(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscDSGetJacobianPreconditioner, $petsc_library),
        PetscErrorCode,
        (
            PetscDS,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscDSSetJacobianPreconditioner(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscDSSetJacobianPreconditioner, $petsc_library),
        PetscErrorCode,
        (
            PetscDS,
            $PetscInt,
            $PetscInt,
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscDSHasDynamicJacobian(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSHasDynamicJacobian, $petsc_library),
        PetscErrorCode,
        (PetscDS, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSGetDynamicJacobian(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscDSGetDynamicJacobian, $petsc_library),
        PetscErrorCode,
        (
            PetscDS,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscDSSetDynamicJacobian(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscDSSetDynamicJacobian, $petsc_library),
        PetscErrorCode,
        (
            PetscDS,
            $PetscInt,
            $PetscInt,
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscDSGetRiemannSolver(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDSGetRiemannSolver, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDSSetRiemannSolver(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDSSetRiemannSolver, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDSGetUpdate(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDSGetUpdate, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDSSetUpdate(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDSSetUpdate, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDSGetContext(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDSGetContext, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDSSetContext(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDSSetContext, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDSGetBdResidual(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDSGetBdResidual, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDSSetBdResidual(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDSSetBdResidual, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDSHasBdJacobian(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSHasBdJacobian, $petsc_library),
        PetscErrorCode,
        (PetscDS, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSGetBdJacobian(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscDSGetBdJacobian, $petsc_library),
        PetscErrorCode,
        (
            PetscDS,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscDSSetBdJacobian(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscDSSetBdJacobian, $petsc_library),
        PetscErrorCode,
        (
            PetscDS,
            $PetscInt,
            $PetscInt,
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscDSHasBdJacobianPreconditioner(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDSHasBdJacobianPreconditioner, $petsc_library),
        PetscErrorCode,
        (PetscDS, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSGetBdJacobianPreconditioner(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscDSGetBdJacobianPreconditioner, $petsc_library),
        PetscErrorCode,
        (
            PetscDS,
            $PetscInt,
            $PetscInt,
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscDSSetBdJacobianPreconditioner(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscDSSetBdJacobianPreconditioner, $petsc_library),
        PetscErrorCode,
        (
            PetscDS,
            $PetscInt,
            $PetscInt,
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscDSGetExactSolution(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDSGetExactSolution, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDSSetExactSolution(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDSSetExactSolution, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDSGetExactSolutionTimeDerivative(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDSGetExactSolutionTimeDerivative, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDSSetExactSolutionTimeDerivative(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDSSetExactSolutionTimeDerivative, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDSGetEvaluationArrays(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDSGetEvaluationArrays, $petsc_library),
        PetscErrorCode,
        (
            PetscDS,
            Ptr{Ptr{$PetscScalar}},
            Ptr{Ptr{$PetscScalar}},
            Ptr{Ptr{$PetscScalar}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDSGetWeakFormArrays(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:PetscDSGetWeakFormArrays, $petsc_library),
        PetscErrorCode,
        (
            PetscDS,
            Ptr{Ptr{$PetscScalar}},
            Ptr{Ptr{$PetscScalar}},
            Ptr{Ptr{$PetscScalar}},
            Ptr{Ptr{$PetscScalar}},
            Ptr{Ptr{$PetscScalar}},
            Ptr{Ptr{$PetscScalar}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function PetscDSGetWorkspace(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscDSGetWorkspace, $petsc_library),
        PetscErrorCode,
        (
            PetscDS,
            Ptr{Ptr{$PetscReal}},
            Ptr{Ptr{$PetscScalar}},
            Ptr{Ptr{$PetscScalar}},
            Ptr{Ptr{$PetscScalar}},
            Ptr{Ptr{$PetscScalar}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscDSCopyConstants(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSCopyConstants, $petsc_library),
        PetscErrorCode,
        (PetscDS, PetscDS),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSCopyEquations(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSCopyEquations, $petsc_library),
        PetscErrorCode,
        (PetscDS, PetscDS),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSSelectDiscretizations(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDSSelectDiscretizations, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{$PetscInt}, PetscDS),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDSSelectEquations(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDSSelectEquations, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{$PetscInt}, PetscDS),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDSAddBoundary(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
)
    @chk ccall(
        (:PetscDSAddBoundary, $petsc_library),
        PetscErrorCode,
        (
            PetscDS,
            DMBoundaryConditionType,
            Ptr{Cchar},
            Ptr{Cchar},
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Cvoid},
            Ptr{Cvoid},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
    )
end

@for_petsc function PetscDSUpdateBoundary(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
)
    @chk ccall(
        (:PetscDSUpdateBoundary, $petsc_library),
        PetscErrorCode,
        (
            PetscDS,
            $PetscInt,
            DMBoundaryConditionType,
            Ptr{Cchar},
            Ptr{Cchar},
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Cvoid},
            Ptr{Cvoid},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
    )
end

@for_petsc function PetscDSGetNumBoundary(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDSGetNumBoundary, $petsc_library),
        PetscErrorCode,
        (PetscDS, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDSGetBoundary(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
)
    @chk ccall(
        (:PetscDSGetBoundary, $petsc_library),
        PetscErrorCode,
        (
            PetscDS,
            $PetscInt,
            Ptr{DMBoundaryConditionType},
            Ptr{Ptr{Cchar}},
            Ptr{Ptr{Cchar}},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{Cvoid}},
            Ptr{Ptr{Cvoid}},
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{Cvoid}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
        arg12,
        arg13,
    )
end

@for_petsc function PetscDSCopyBoundary(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDSCopyBoundary, $petsc_library),
        PetscErrorCode,
        (PetscDS, $PetscInt, Ptr{$PetscInt}, PetscDS),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function CharacteristicInitializePackage(::$UnionPetscLib)
    @chk ccall(
        (:CharacteristicInitializePackage, $petsc_library),
        PetscErrorCode,
        (),
    )
end

mutable struct _p_Characteristic end

const Characteristic = Ptr{_p_Characteristic}

const CharacteristicType = Ptr{Cchar}

@for_petsc function CharacteristicCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:CharacteristicCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Characteristic}),
        arg1,
        arg2,
    )
end

@for_petsc function CharacteristicSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:CharacteristicSetType, $petsc_library),
        PetscErrorCode,
        (Characteristic, CharacteristicType),
        arg1,
        arg2,
    )
end

@for_petsc function CharacteristicSetUp(::$UnionPetscLib, arg1)
    @chk ccall(
        (:CharacteristicSetUp, $petsc_library),
        PetscErrorCode,
        (Characteristic,),
        arg1,
    )
end

@for_petsc function CharacteristicSetVelocityInterpolation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:CharacteristicSetVelocityInterpolation, $petsc_library),
        PetscErrorCode,
        (
            Characteristic,
            DM,
            Vec,
            Vec,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Cvoid},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function CharacteristicSetVelocityInterpolationLocal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
)
    @chk ccall(
        (:CharacteristicSetVelocityInterpolationLocal, $petsc_library),
        PetscErrorCode,
        (
            Characteristic,
            DM,
            Vec,
            Vec,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Cvoid},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
    )
end

@for_petsc function CharacteristicSetFieldInterpolation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:CharacteristicSetFieldInterpolation, $petsc_library),
        PetscErrorCode,
        (
            Characteristic,
            DM,
            Vec,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Cvoid},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function CharacteristicSetFieldInterpolationLocal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
)
    @chk ccall(
        (:CharacteristicSetFieldInterpolationLocal, $petsc_library),
        PetscErrorCode,
        (
            Characteristic,
            DM,
            Vec,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Cvoid},
            Ptr{Cvoid},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
    )
end

@for_petsc function CharacteristicSolve(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:CharacteristicSolve, $petsc_library),
        PetscErrorCode,
        (Characteristic, $PetscReal, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function CharacteristicDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:CharacteristicDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{Characteristic},),
        arg1,
    )
end

@for_petsc function CharacteristicRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:CharacteristicRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
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

@for_petsc function PCInitializePackage(::$UnionPetscLib)
    @chk ccall((:PCInitializePackage, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PCCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PC}),
        arg1,
        arg2,
    )
end

@for_petsc function PCSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCSetType, $petsc_library),
        PetscErrorCode,
        (PC, PCType),
        arg1,
        arg2,
    )
end

@for_petsc function PCGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCGetType, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{PCType}),
        arg1,
        arg2,
    )
end

@for_petsc function PCSetUp(::$UnionPetscLib, arg1)
    @chk ccall((:PCSetUp, $petsc_library), PetscErrorCode, (PC,), arg1)
end

@for_petsc function PCSetFailedReason(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCSetFailedReason, $petsc_library),
        PetscErrorCode,
        (PC, PCFailedReason),
        arg1,
        arg2,
    )
end

@for_petsc function PCGetFailedReason(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCGetFailedReason, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{PCFailedReason}),
        arg1,
        arg2,
    )
end

@for_petsc function PCGetSetUpFailedReason(::$UnionPetscLib, pc, reason)
    @chk ccall(
        (:PCGetSetUpFailedReason, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{PCFailedReason}),
        pc,
        reason,
    )
end

@for_petsc function PCGetFailedReasonRank(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCGetFailedReasonRank, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{PCFailedReason}),
        arg1,
        arg2,
    )
end

@for_petsc function PCSetUpOnBlocks(::$UnionPetscLib, arg1)
    @chk ccall((:PCSetUpOnBlocks, $petsc_library), PetscErrorCode, (PC,), arg1)
end

@for_petsc function PCApply(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCApply, $petsc_library),
        PetscErrorCode,
        (PC, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCMatApply(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCMatApply, $petsc_library),
        PetscErrorCode,
        (PC, Mat, Mat),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCApplySymmetricLeft(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCApplySymmetricLeft, $petsc_library),
        PetscErrorCode,
        (PC, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCApplySymmetricRight(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCApplySymmetricRight, $petsc_library),
        PetscErrorCode,
        (PC, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCApplyBAorAB(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PCApplyBAorAB, $petsc_library),
        PetscErrorCode,
        (PC, PCSide, Vec, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PCApplyTranspose(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCApplyTranspose, $petsc_library),
        PetscErrorCode,
        (PC, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCApplyTransposeExists(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCApplyTransposeExists, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PCApplyBAorABTranspose(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PCApplyBAorABTranspose, $petsc_library),
        PetscErrorCode,
        (PC, PCSide, Vec, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PCSetReusePreconditioner(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCSetReusePreconditioner, $petsc_library),
        PetscErrorCode,
        (PC, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PCGetReusePreconditioner(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCGetReusePreconditioner, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PCSetErrorIfFailure(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCSetErrorIfFailure, $petsc_library),
        PetscErrorCode,
        (PC, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PCApplyRichardson(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
)
    @chk ccall(
        (:PCApplyRichardson, $petsc_library),
        PetscErrorCode,
        (
            PC,
            Vec,
            Vec,
            Vec,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscInt,
            PetscBool,
            Ptr{$PetscInt},
            Ptr{PCRichardsonConvergedReason},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
        arg11,
    )
end

@for_petsc function PCApplyRichardsonExists(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCApplyRichardsonExists, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PCSetUseAmat(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCSetUseAmat, $petsc_library),
        PetscErrorCode,
        (PC, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PCGetUseAmat(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCGetUseAmat, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PCRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PCReset(::$UnionPetscLib, arg1)
    @chk ccall((:PCReset, $petsc_library), PetscErrorCode, (PC,), arg1)
end

@for_petsc function PCDestroy(::$UnionPetscLib, arg1)
    @chk ccall((:PCDestroy, $petsc_library), PetscErrorCode, (Ptr{PC},), arg1)
end

@for_petsc function PCSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall((:PCSetFromOptions, $petsc_library), PetscErrorCode, (PC,), arg1)
end

@for_petsc function PCFactorGetMatrix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCFactorGetMatrix, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function PCSetModifySubMatrices(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCSetModifySubMatrices, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCModifySubMatrices(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PCModifySubMatrices, $petsc_library),
        PetscErrorCode,
        (PC, $PetscInt, Ptr{IS}, Ptr{IS}, Ptr{Mat}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PCSetOperators(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCSetOperators, $petsc_library),
        PetscErrorCode,
        (PC, Mat, Mat),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCGetOperators(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCGetOperators, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{Mat}, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCGetOperatorsSet(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCGetOperatorsSet, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{PetscBool}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCView, $petsc_library),
        PetscErrorCode,
        (PC, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PCLoad(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCLoad, $petsc_library),
        PetscErrorCode,
        (PC, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PCViewFromOptions(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCViewFromOptions, $petsc_library),
        PetscErrorCode,
        (PC, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCSetOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCSetOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PCAppendOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCAppendOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PCGetOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCGetOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PCComputeOperator(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCComputeOperator, $petsc_library),
        PetscErrorCode,
        (PC, MatType, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCComputeExplicitOperator(::$UnionPetscLib, A, B)
    @chk ccall(
        (:PCComputeExplicitOperator, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{Mat}),
        A,
        B,
    )
end

@for_petsc function PCGetDiagonalScale(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCGetDiagonalScale, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PCDiagonalScaleLeft(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCDiagonalScaleLeft, $petsc_library),
        PetscErrorCode,
        (PC, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCDiagonalScaleRight(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCDiagonalScaleRight, $petsc_library),
        PetscErrorCode,
        (PC, Vec, Vec),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCSetDiagonalScale(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCSetDiagonalScale, $petsc_library),
        PetscErrorCode,
        (PC, Vec),
        arg1,
        arg2,
    )
end

@for_petsc function PCSetDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall((:PCSetDM, $petsc_library), PetscErrorCode, (PC, DM), arg1, arg2)
end

@for_petsc function PCGetDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCGetDM, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{DM}),
        arg1,
        arg2,
    )
end

@for_petsc function PCGetInterpolations(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCGetInterpolations, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{$PetscInt}, Ptr{Ptr{Mat}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCGetCoarseOperators(::$UnionPetscLib, pc, arg2, arg3)
    @chk ccall(
        (:PCGetCoarseOperators, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{$PetscInt}, Ptr{Ptr{Mat}}),
        pc,
        arg2,
        arg3,
    )
end

@for_petsc function PCSetCoordinates(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PCSetCoordinates, $petsc_library),
        PetscErrorCode,
        (PC, $PetscInt, $PetscInt, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PCSetApplicationContext(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCSetApplicationContext, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PCGetApplicationContext(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCGetApplicationContext, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PCJacobiSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCJacobiSetType, $petsc_library),
        PetscErrorCode,
        (PC, PCJacobiType),
        arg1,
        arg2,
    )
end

@for_petsc function PCJacobiGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCJacobiGetType, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{PCJacobiType}),
        arg1,
        arg2,
    )
end

@for_petsc function PCJacobiSetUseAbs(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCJacobiSetUseAbs, $petsc_library),
        PetscErrorCode,
        (PC, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PCJacobiGetUseAbs(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCJacobiGetUseAbs, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PCSORSetSymmetric(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCSORSetSymmetric, $petsc_library),
        PetscErrorCode,
        (PC, MatSORType),
        arg1,
        arg2,
    )
end

@for_petsc function PCSORGetSymmetric(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCSORGetSymmetric, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{MatSORType}),
        arg1,
        arg2,
    )
end

@for_petsc function PCSORSetOmega(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCSORSetOmega, $petsc_library),
        PetscErrorCode,
        (PC, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function PCSORGetOmega(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCSORGetOmega, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function PCSORSetIterations(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCSORSetIterations, $petsc_library),
        PetscErrorCode,
        (PC, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCSORGetIterations(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCSORGetIterations, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCEisenstatSetOmega(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCEisenstatSetOmega, $petsc_library),
        PetscErrorCode,
        (PC, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function PCEisenstatGetOmega(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCEisenstatGetOmega, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function PCEisenstatSetNoDiagonalScaling(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PCEisenstatSetNoDiagonalScaling, $petsc_library),
        PetscErrorCode,
        (PC, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PCEisenstatGetNoDiagonalScaling(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PCEisenstatGetNoDiagonalScaling, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PCBJacobiSetTotalBlocks(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCBJacobiSetTotalBlocks, $petsc_library),
        PetscErrorCode,
        (PC, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCBJacobiGetTotalBlocks(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCBJacobiGetTotalBlocks, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCBJacobiSetLocalBlocks(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCBJacobiSetLocalBlocks, $petsc_library),
        PetscErrorCode,
        (PC, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCBJacobiGetLocalBlocks(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCBJacobiGetLocalBlocks, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCShellSetApply(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCShellSetApply, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PCShellSetMatApply(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCShellSetMatApply, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PCShellSetApplySymmetricLeft(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCShellSetApplySymmetricLeft, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PCShellSetApplySymmetricRight(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCShellSetApplySymmetricRight, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PCShellSetApplyBA(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCShellSetApplyBA, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PCShellSetApplyTranspose(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCShellSetApplyTranspose, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PCShellSetSetUp(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCShellSetSetUp, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PCShellSetApplyRichardson(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCShellSetApplyRichardson, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PCShellSetView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCShellSetView, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PCShellSetDestroy(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCShellSetDestroy, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PCShellSetContext(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCShellSetContext, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PCShellGetContext(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCShellGetContext, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
    )
end

@for_petsc function PCShellSetName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCShellSetName, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PCShellGetName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCShellGetName, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PCFactorSetZeroPivot(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCFactorSetZeroPivot, $petsc_library),
        PetscErrorCode,
        (PC, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function PCFactorSetShiftType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCFactorSetShiftType, $petsc_library),
        PetscErrorCode,
        (PC, MatFactorShiftType),
        arg1,
        arg2,
    )
end

@for_petsc function PCFactorSetShiftAmount(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCFactorSetShiftAmount, $petsc_library),
        PetscErrorCode,
        (PC, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function PCFactorSetMatSolverType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCFactorSetMatSolverType, $petsc_library),
        PetscErrorCode,
        (PC, MatSolverType),
        arg1,
        arg2,
    )
end

@for_petsc function PCFactorGetMatSolverType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCFactorGetMatSolverType, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{MatSolverType}),
        arg1,
        arg2,
    )
end

@for_petsc function PCFactorSetUpMatSolverType(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PCFactorSetUpMatSolverType, $petsc_library),
        PetscErrorCode,
        (PC,),
        arg1,
    )
end

@for_petsc function PCFactorSetMatSolverPackage(::$UnionPetscLib, pc, stype)
    @chk ccall(
        (:PCFactorSetMatSolverPackage, $petsc_library),
        PetscErrorCode,
        (PC, MatSolverType),
        pc,
        stype,
    )
end

@for_petsc function PCFactorGetMatSolverPackage(::$UnionPetscLib, pc, stype)
    @chk ccall(
        (:PCFactorGetMatSolverPackage, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{MatSolverType}),
        pc,
        stype,
    )
end

@for_petsc function PCFactorSetUpMatSolverPackage(::$UnionPetscLib, pc)
    @chk ccall(
        (:PCFactorSetUpMatSolverPackage, $petsc_library),
        PetscErrorCode,
        (PC,),
        pc,
    )
end

@for_petsc function PCFactorSetFill(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCFactorSetFill, $petsc_library),
        PetscErrorCode,
        (PC, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function PCFactorSetColumnPivot(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCFactorSetColumnPivot, $petsc_library),
        PetscErrorCode,
        (PC, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function PCFactorReorderForNonzeroDiagonal(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PCFactorReorderForNonzeroDiagonal, $petsc_library),
        PetscErrorCode,
        (PC, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function PCFactorSetMatOrderingType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCFactorSetMatOrderingType, $petsc_library),
        PetscErrorCode,
        (PC, MatOrderingType),
        arg1,
        arg2,
    )
end

@for_petsc function PCFactorSetReuseOrdering(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCFactorSetReuseOrdering, $petsc_library),
        PetscErrorCode,
        (PC, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PCFactorSetReuseFill(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCFactorSetReuseFill, $petsc_library),
        PetscErrorCode,
        (PC, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PCFactorSetUseInPlace(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCFactorSetUseInPlace, $petsc_library),
        PetscErrorCode,
        (PC, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PCFactorGetUseInPlace(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCFactorGetUseInPlace, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PCFactorSetAllowDiagonalFill(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCFactorSetAllowDiagonalFill, $petsc_library),
        PetscErrorCode,
        (PC, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PCFactorGetAllowDiagonalFill(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCFactorGetAllowDiagonalFill, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PCFactorSetPivotInBlocks(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCFactorSetPivotInBlocks, $petsc_library),
        PetscErrorCode,
        (PC, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PCFactorSetLevels(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCFactorSetLevels, $petsc_library),
        PetscErrorCode,
        (PC, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PCFactorGetLevels(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCFactorGetLevels, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PCFactorSetDropTolerance(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PCFactorSetDropTolerance, $petsc_library),
        PetscErrorCode,
        (PC, $PetscReal, $PetscReal, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PCFactorGetZeroPivot(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCFactorGetZeroPivot, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function PCFactorGetShiftAmount(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCFactorGetShiftAmount, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function PCFactorGetShiftType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCFactorGetShiftType, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{MatFactorShiftType}),
        arg1,
        arg2,
    )
end

@for_petsc function PCASMSetLocalSubdomains(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PCASMSetLocalSubdomains, $petsc_library),
        PetscErrorCode,
        (PC, $PetscInt, Ptr{IS}, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PCASMSetTotalSubdomains(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PCASMSetTotalSubdomains, $petsc_library),
        PetscErrorCode,
        (PC, $PetscInt, Ptr{IS}, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PCASMSetOverlap(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCASMSetOverlap, $petsc_library),
        PetscErrorCode,
        (PC, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PCASMSetDMSubdomains(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCASMSetDMSubdomains, $petsc_library),
        PetscErrorCode,
        (PC, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PCASMGetDMSubdomains(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCASMGetDMSubdomains, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PCASMSetSortIndices(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCASMSetSortIndices, $petsc_library),
        PetscErrorCode,
        (PC, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PCASMSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCASMSetType, $petsc_library),
        PetscErrorCode,
        (PC, PCASMType),
        arg1,
        arg2,
    )
end

@for_petsc function PCASMGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCASMGetType, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{PCASMType}),
        arg1,
        arg2,
    )
end

@for_petsc function PCASMSetLocalType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCASMSetLocalType, $petsc_library),
        PetscErrorCode,
        (PC, PCCompositeType),
        arg1,
        arg2,
    )
end

@for_petsc function PCASMGetLocalType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCASMGetLocalType, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{PCCompositeType}),
        arg1,
        arg2,
    )
end

@for_petsc function PCASMCreateSubdomains(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCASMCreateSubdomains, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{Ptr{IS}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCASMDestroySubdomains(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCASMDestroySubdomains, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{IS}, Ptr{IS}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCASMCreateSubdomains2D(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
)
    @chk ccall(
        (:PCASMCreateSubdomains2D, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{IS}},
            Ptr{Ptr{IS}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
    )
end

@for_petsc function PCASMGetLocalSubdomains(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PCASMGetLocalSubdomains, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{$PetscInt}, Ptr{Ptr{IS}}, Ptr{Ptr{IS}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PCASMGetLocalSubmatrices(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCASMGetLocalSubmatrices, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{$PetscInt}, Ptr{Ptr{Mat}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCASMGetSubMatType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCASMGetSubMatType, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{MatType}),
        arg1,
        arg2,
    )
end

@for_petsc function PCASMSetSubMatType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCASMSetSubMatType, $petsc_library),
        PetscErrorCode,
        (PC, MatType),
        arg1,
        arg2,
    )
end

@for_petsc function PCGASMSetTotalSubdomains(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCGASMSetTotalSubdomains, $petsc_library),
        PetscErrorCode,
        (PC, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PCGASMSetSubdomains(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PCGASMSetSubdomains, $petsc_library),
        PetscErrorCode,
        (PC, $PetscInt, Ptr{IS}, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PCGASMSetOverlap(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCGASMSetOverlap, $petsc_library),
        PetscErrorCode,
        (PC, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PCGASMSetUseDMSubdomains(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCGASMSetUseDMSubdomains, $petsc_library),
        PetscErrorCode,
        (PC, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PCGASMGetUseDMSubdomains(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCGASMGetUseDMSubdomains, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PCGASMSetSortIndices(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCGASMSetSortIndices, $petsc_library),
        PetscErrorCode,
        (PC, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PCGASMSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCGASMSetType, $petsc_library),
        PetscErrorCode,
        (PC, PCGASMType),
        arg1,
        arg2,
    )
end

@for_petsc function PCGASMCreateSubdomains(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PCGASMCreateSubdomains, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}, Ptr{Ptr{IS}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PCGASMDestroySubdomains(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCGASMDestroySubdomains, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{Ptr{IS}}, Ptr{Ptr{IS}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCGASMCreateSubdomains2D(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
)
    @chk ccall(
        (:PCGASMCreateSubdomains2D, $petsc_library),
        PetscErrorCode,
        (
            PC,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{IS}},
            Ptr{Ptr{IS}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
        arg8,
        arg9,
        arg10,
    )
end

@for_petsc function PCGASMGetSubdomains(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PCGASMGetSubdomains, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{$PetscInt}, Ptr{Ptr{IS}}, Ptr{Ptr{IS}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PCGASMGetSubmatrices(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PCGASMGetSubmatrices, $petsc_library),
        PetscErrorCode,
        (PC, Ptr{$PetscInt}, Ptr{Ptr{Mat}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PCCompositeSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PCCompositeSetType, $petsc_library),
        PetscErrorCode,
    )
end