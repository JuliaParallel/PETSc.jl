#
# START OF PROLOGUE
#

using MPI
const MPI_Comm = MPI.MPI_Comm
const MPI_Datatype = MPI.MPI_Datatype
const MPI_File = MPI.MPI_File
const MPI_Aint = MPI.MPI_Aint
const MPI_Info = MPI.MPI_Info
const MPI_Win = MPI.MPI_Win
const MPI_Offset = MPI.MPI_Offset
const MPI_Op = MPI.MPI_Op
const MPI_UNSIGNED_SHORT = MPI.UNSIGNED_SHORT
const MPI_INT64_T = MPI.INT64_T
const MPI_INT32_T = MPI.INT32_T
const MPI_FLOAT = MPI.FLOAT
const MPI_COMM_SELF = MPI.COMM_SELF
const MPI_DOUBLE = MPI.DOUBLE
const MPI_SUM = MPI.SUM
const MPI_MAX = MPI.MAX
const MPI_MIN = MPI.MIN
const MPI_REPLACE = MPI.REPLACE
const MPIU_INT64 = MPI.UINT64_T
const MPIU_INT32 = MPI.UINT32_T

# We know these will be Cvoid, so just set them to be that
const PetscOptions = Ptr{Cvoid}
const PetscViewer = Ptr{Cvoid}
const PetscObject = Ptr{Cvoid}
const Vec = Ptr{Cvoid}
const VecType = Cstring
const Mat = Ptr{Cvoid}
const MatType = Cstring
const KSP = Ptr{Cvoid}
const KSPType = Cstring
const SNES = Ptr{Cvoid}
const SNESType = Cstring
const DM = Ptr{Cvoid}

#
# END OF PROLOGUE
#

const __darwin_off_t = Int64

mutable struct ADIOI_FileD end

const off_t = __darwin_off_t

const PetscInt64 = Int64

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

const PetscMPIInt = Cint

@for_petsc function PetscMPIErrorString(::$UnionPetscLib, arg1, arg2)
    ccall(
        (:PetscMPIErrorString, $petsc_library),
        Cvoid,
        (PetscMPIInt, Ptr{Cchar}),
        arg1,
        arg2,
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

@enum PetscFPTrap::UInt32 begin
    PETSC_FP_TRAP_OFF = 0
    PETSC_FP_TRAP_INDIV = 1
    PETSC_FP_TRAP_FLTOPERR = 2
    PETSC_FP_TRAP_FLTOVF = 4
    PETSC_FP_TRAP_FLTUND = 8
    PETSC_FP_TRAP_FLTDIV = 16
    PETSC_FP_TRAP_FLTINEX = 32
end

const PetscClassId = Cint

# typedef void ( PetscVoidFn ) ( void )
const PetscVoidFn = Cvoid

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

@for_petsc function PetscMemzero(::$UnionPetscLib, a, n)
    @chk ccall(
        (:PetscMemzero, $petsc_library),
        PetscErrorCode,
        (Ptr{Cvoid}, Csize_t),
        a,
        n,
    )
end

@enum PetscEnum::UInt32 begin
    ENUM_DUMMY = 0
end

mutable struct _n_PetscFunctionList end

const PetscFunctionList = Ptr{_n_PetscFunctionList}

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

const PetscBT = Ptr{Cchar}

@for_petsc function PetscBTLookup(::$UnionPetscLib, array, index)
    ccall(
        (:PetscBTLookup, $petsc_library),
        Cchar,
        (PetscBT, $PetscInt),
        array,
        index,
    )
end

@for_petsc function PetscLogObjectParent(::$UnionPetscLib, o, p)
    @chk ccall(
        (:PetscLogObjectParent, $petsc_library),
        PetscErrorCode,
        (PetscObject, PetscObject),
        o,
        p,
    )
end

const PetscLogEvent = Cint

@for_petsc function PetscLogEventBegin_Internal(
    ::$UnionPetscLib,
    e,
    o1,
    o2,
    o3,
    o4,
)
    @chk ccall(
        (:PetscLogEventBegin_Internal, $petsc_library),
        PetscErrorCode,
        (PetscLogEvent, PetscObject, PetscObject, PetscObject, PetscObject),
        e,
        o1,
        o2,
        o3,
        o4,
    )
end

@for_petsc function PetscLogEventEnd_Internal(
    ::$UnionPetscLib,
    e,
    o1,
    o2,
    o3,
    o4,
)
    @chk ccall(
        (:PetscLogEventEnd_Internal, $petsc_library),
        PetscErrorCode,
        (PetscLogEvent, PetscObject, PetscObject, PetscObject, PetscObject),
        e,
        o1,
        o2,
        o3,
        o4,
    )
end

const PetscLogDouble = Cdouble

@for_petsc function PetscMPITypeSize(
    ::$UnionPetscLib,
    count,
    type,
    length,
    length_th,
)
    @chk ccall(
        (:PetscMPITypeSize, $petsc_library),
        PetscErrorCode,
        ($PetscInt, MPI_Datatype, Ptr{PetscLogDouble}, Ptr{PetscLogDouble}),
        count,
        type,
        length,
        length_th,
    )
end

@for_petsc function PetscMPIParallelComm(::$UnionPetscLib, comm)
    ccall((:PetscMPIParallelComm, $petsc_library), Cint, (MPI_Comm,), comm)
end

@for_petsc function PetscMPITypeSizeComm(
    ::$UnionPetscLib,
    comm,
    counts,
    type,
    length,
    length_th,
)
    @chk ccall(
        (:PetscMPITypeSizeComm, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            Ptr{PetscMPIInt},
            MPI_Datatype,
            Ptr{PetscLogDouble},
            Ptr{PetscLogDouble},
        ),
        comm,
        counts,
        type,
        length,
        length_th,
    )
end

const PetscLogStage = Cint

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

@for_petsc function PetscLogStageGetId(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogStageGetId, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{PetscLogStage}),
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

@for_petsc function PetscLogStageRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogStageRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{PetscLogStage}),
        arg1,
        arg2,
    )
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

mutable struct _p_PetscMatlabEngine end

const PetscMatlabEngine = Ptr{_p_PetscMatlabEngine}

@for_petsc function PETSC_MATLAB_ENGINE_(::$UnionPetscLib, arg1)
    ccall(
        (:PETSC_MATLAB_ENGINE_, $petsc_library),
        PetscMatlabEngine,
        (MPI_Comm,),
        arg1,
    )
end

mutable struct _p_PetscDeviceContext end

const PetscDeviceContext = Ptr{_p_PetscDeviceContext}

@enum PetscMemType::UInt32 begin
    PETSC_MEMTYPE_HOST = 0
    PETSC_MEMTYPE_DEVICE = 1
    # PETSC_MEMTYPE_CUDA = 1
    PETSC_MEMTYPE_NVSHMEM = 17
    PETSC_MEMTYPE_HIP = 3
    PETSC_MEMTYPE_SYCL = 5
end

@for_petsc function PetscDeviceAllocate_Private(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscDeviceAllocate_Private, $petsc_library),
        PetscErrorCode,
        (
            PetscDeviceContext,
            PetscBool,
            PetscMemType,
            Csize_t,
            Csize_t,
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

@for_petsc function PetscDeviceDeallocate_Private(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDeviceDeallocate_Private, $petsc_library),
        PetscErrorCode,
        (PetscDeviceContext, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDeviceMemcpy(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscDeviceMemcpy, $petsc_library),
        PetscErrorCode,
        (PetscDeviceContext, Ptr{Cvoid}, Ptr{Cvoid}, Csize_t),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDeviceMemset(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:PetscDeviceMemset, $petsc_library),
        PetscErrorCode,
        (PetscDeviceContext, Ptr{Cvoid}, $PetscInt, Csize_t),
        arg1,
        arg2,
        arg3,
        arg4,
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
    PETSC_VIEWER_ALL = 39
end

@for_petsc function PETSC_VIEWER_STDERR_(::$UnionPetscLib, arg1)
    ccall(
        (:PETSC_VIEWER_STDERR_, $petsc_library),
        PetscViewer,
        (MPI_Comm,),
        arg1,
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

mutable struct _p_ISLocalToGlobalMapping end

const ISLocalToGlobalMapping = Ptr{_p_ISLocalToGlobalMapping}

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

@for_petsc function PetscSortRemoveDupsInt(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSortRemoveDupsInt, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
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

@for_petsc function PetscObjectSetOptionsPrefix(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscObjectSetOptionsPrefix, $petsc_library),
        PetscErrorCode,
        (PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@enum KSPConvergedReason::Int32 begin
    KSP_CONVERGED_RTOL_NORMAL = 1
    KSP_CONVERGED_ATOL_NORMAL = 9
    KSP_CONVERGED_RTOL = 2
    KSP_CONVERGED_ATOL = 3
    KSP_CONVERGED_ITS = 4
    KSP_CONVERGED_NEG_CURVE = 5
    # KSP_CONVERGED_CG_NEG_CURVE = 5
    KSP_CONVERGED_CG_CONSTRAINED = 6
    # KSP_CONVERGED_STEP_LENGTH = 6
    KSP_CONVERGED_HAPPY_BREAKDOWN = 7
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

@for_petsc function KSPConvergedDefault(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:KSPConvergedDefault, $petsc_library),
        PetscErrorCode,
        (KSP, $PetscInt, $PetscReal, Ptr{KSPConvergedReason}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function KSPConvergedDefaultDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:KSPConvergedDefaultDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{Cvoid},),
        arg1,
    )
end

@for_petsc function KSPConvergedDefaultCreate(::$UnionPetscLib, arg1)
    @chk ccall(
        (:KSPConvergedDefaultCreate, $petsc_library),
        PetscErrorCode,
        (Ptr{Ptr{Cvoid}},),
        arg1,
    )
end

@for_petsc function KSPConvergedDefaultSetUIRNorm(::$UnionPetscLib, arg1)
    @chk ccall(
        (:KSPConvergedDefaultSetUIRNorm, $petsc_library),
        PetscErrorCode,
        (KSP,),
        arg1,
    )
end

@for_petsc function KSPConvergedDefaultSetUMIRNorm(::$UnionPetscLib, arg1)
    @chk ccall(
        (:KSPConvergedDefaultSetUMIRNorm, $petsc_library),
        PetscErrorCode,
        (KSP,),
        arg1,
    )
end

@for_petsc function KSPConvergedSkip(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:KSPConvergedSkip, $petsc_library),
        PetscErrorCode,
        (KSP, $PetscInt, $PetscReal, Ptr{KSPConvergedReason}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@enum SNESConvergedReason::Int32 begin
    SNES_CONVERGED_FNORM_ABS = 2
    SNES_CONVERGED_FNORM_RELATIVE = 3
    SNES_CONVERGED_SNORM_RELATIVE = 4
    SNES_CONVERGED_ITS = 5
    SNES_BREAKOUT_INNER_ITER = 6
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

@for_petsc function SNESConvergedSkip(
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
        (:SNESConvergedSkip, $petsc_library),
        PetscErrorCode,
        (
            SNES,
            $PetscInt,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            Ptr{SNESConvergedReason},
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

@enum __JL_Ctag_2::Int32 begin
    PETSC_SUCCESS = 0
    PETSC_ERR_BOOLEAN_MACRO_FAILURE = 1
    PETSC_ERR_MIN_VALUE = 54
    PETSC_ERR_MEM = 55
    PETSC_ERR_SUP = 56
    PETSC_ERR_SUP_SYS = 57
    PETSC_ERR_ORDER = 58
    PETSC_ERR_SIG = 59
    PETSC_ERR_FP = 72
    PETSC_ERR_COR = 74
    PETSC_ERR_LIB = 76
    PETSC_ERR_PLIB = 77
    PETSC_ERR_MEMC = 78
    PETSC_ERR_CONV_FAILED = 82
    PETSC_ERR_USER = 83
    PETSC_ERR_SYS = 88
    PETSC_ERR_POINTER = 70
    PETSC_ERR_MPI_LIB_INCOMP = 87
    PETSC_ERR_ARG_SIZ = 60
    PETSC_ERR_ARG_IDN = 61
    PETSC_ERR_ARG_WRONG = 62
    PETSC_ERR_ARG_CORRUPT = 64
    PETSC_ERR_ARG_OUTOFRANGE = 63
    PETSC_ERR_ARG_BADPTR = 68
    PETSC_ERR_ARG_NOTSAMETYPE = 69
    PETSC_ERR_ARG_NOTSAMECOMM = 80
    PETSC_ERR_ARG_WRONGSTATE = 73
    PETSC_ERR_ARG_TYPENOTSET = 89
    PETSC_ERR_ARG_INCOMP = 75
    PETSC_ERR_ARG_NULL = 85
    PETSC_ERR_ARG_UNKNOWN_TYPE = 86
    PETSC_ERR_FILE_OPEN = 65
    PETSC_ERR_FILE_READ = 66
    PETSC_ERR_FILE_WRITE = 67
    PETSC_ERR_FILE_UNEXPECTED = 79
    PETSC_ERR_MAT_LU_ZRPVT = 71
    PETSC_ERR_MAT_CH_ZRPVT = 81
    PETSC_ERR_INT_OVERFLOW = 84
    PETSC_ERR_FLOP_COUNT = 90
    PETSC_ERR_NOT_CONVERGED = 91
    PETSC_ERR_MISSING_FACTOR = 92
    PETSC_ERR_OPT_OVERWRITE = 93
    PETSC_ERR_WRONG_MPI_SIZE = 94
    PETSC_ERR_USER_INPUT = 95
    PETSC_ERR_GPU_RESOURCE = 96
    PETSC_ERR_GPU = 97
    PETSC_ERR_MPI = 98
    PETSC_ERR_RETURN = 99
    PETSC_ERR_MEM_LEAK = 100
    PETSC_ERR_MAX_VALUE = 101
    PETSC_ERR_MIN_SIGNED_BOUND_DO_NOT_USE = -2147483648
    PETSC_ERR_MAX_SIGNED_BOUND_DO_NOT_USE = 2147483647
end

@enum __JL_Ctag_3::Int32 begin
    PETSC_MPI_INT_MIN = -2147483648
    PETSC_MPI_INT_MAX = 2147483647
end

const PetscSizeT = Csize_t

const PetscCount = Cptrdiff_t

const PetscShort = Cshort

const PetscChar = Cchar

const PetscInt32 = Int32

const PetscBLASInt = Cint

@enum __JL_Ctag_7::Int32 begin
    PETSC_BLAS_INT_MIN = -2147483648
    PETSC_BLAS_INT_MAX = 2147483647
end

const PetscCuBLASInt = Cint

@enum __JL_Ctag_8::Int32 begin
    PETSC_CUBLAS_INT_MIN = -2147483648
    PETSC_CUBLAS_INT_MAX = 2147483647
end

const PetscHipBLASInt = Cint

@enum __JL_Ctag_9::Int32 begin
    PETSC_HIPBLAS_INT_MIN = -2147483648
    PETSC_HIPBLAS_INT_MAX = 2147483647
end

@enum PetscBool3::Int32 begin
    PETSC_BOOL3_FALSE = 0
    PETSC_BOOL3_TRUE = 1
    PETSC_BOOL3_UNKNOWN = -1
end

@enum PetscCopyMode::UInt32 begin
    PETSC_COPY_VALUES = 0
    PETSC_OWN_POINTER = 1
    PETSC_USE_POINTER = 2
end

mutable struct _p_PetscToken end

const PetscToken = Ptr{_p_PetscToken}

const PetscObjectId = PetscInt64

const PetscObjectState = PetscInt64

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

@for_petsc function PetscIsCloseAtTolScalar(
    ::$UnionPetscLib,
    lhs,
    rhs,
    rtol,
    atol,
)
    ccall(
        (:PetscIsCloseAtTolScalar, $petsc_library),
        PetscBool,
        ($PetscScalar, $PetscScalar, $PetscReal, $PetscReal),
        lhs,
        rhs,
        rtol,
        atol,
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

@for_petsc function PetscCommGetComm(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscCommGetComm, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{MPI_Comm}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscCommRestoreComm(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscCommRestoreComm, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{MPI_Comm}),
        arg1,
        arg2,
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
        (PetscErrorCode, Ptr{Ptr{Cchar}}, Ptr{Ptr{Cchar}}),
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
    petscroutine::NTuple{64, Cint}
    currentsize::Cint
    hotdepth::Cint
    check::PetscBool
    PetscStack() = new()
end

@for_petsc function PetscCIFilename(::$UnionPetscLib, arg1)
    ccall((:PetscCIFilename, $petsc_library), Ptr{Cchar}, (Ptr{Cchar},), arg1)
end

@for_petsc function PetscCILinenumber(::$UnionPetscLib, arg1)
    ccall((:PetscCILinenumber, $petsc_library), Cint, (Cint,), arg1)
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

@for_petsc function PetscSysFinalizePackage(::$UnionPetscLib)
    @chk ccall((:PetscSysFinalizePackage, $petsc_library), PetscErrorCode, ())
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

# typedef PetscVoidFn * PetscVoidFunction
const PetscVoidFunction = Ptr{PetscVoidFn}

# typedef PetscVoidFn * * PetscVoidStarFunction
const PetscVoidStarFunction = Ptr{Ptr{PetscVoidFn}}

# typedef PetscErrorCode ( PetscErrorCodeFn ) ( void )
const PetscErrorCodeFn = Cvoid

# typedef PetscErrorCodeFn * PetscErrorCodeFunction
const PetscErrorCodeFunction = Ptr{PetscErrorCodeFn}

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

@enum PetscOptionSource::UInt32 begin
    PETSC_OPT_CODE = 0
    PETSC_OPT_COMMAND_LINE = 1
    PETSC_OPT_FILE = 2
    PETSC_OPT_ENVIRONMENT = 3
    NUM_PETSC_OPT_SOURCE = 4
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
    arg4,
)
    @chk ccall(
        (:PetscOptionsMonitorDefault, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, PetscOptionSource, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
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

@for_petsc function PetscOptionsLeftError(::$UnionPetscLib)
    @chk ccall((:PetscOptionsLeftError, $petsc_library), PetscErrorCode, ())
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

@for_petsc function PetscObjectObjectTypeCompare(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscObjectObjectTypeCompare, $petsc_library),
        PetscErrorCode,
        (PetscObject, PetscObject, Ptr{PetscBool}),
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

@for_petsc function PetscDemangleSymbol(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDemangleSymbol, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscMallocGetStack(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscMallocGetStack, $petsc_library),
        PetscErrorCode,
        (Ptr{Cvoid}, Ptr{Ptr{PetscStack}}),
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

@for_petsc function PetscObjectsView(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscObjectsView, $petsc_library),
        PetscErrorCode,
        (PetscViewer,),
        arg1,
    )
end

@for_petsc function PetscObjectsGetObject(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscObjectsGetObject, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{PetscObject}, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
        arg3,
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

@for_petsc function PetscFunctionListClear(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscFunctionListClear, $petsc_library),
        PetscErrorCode,
        (PetscFunctionList,),
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

@for_petsc function PetscFunctionListPrintNonEmpty(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscFunctionListPrintNonEmpty, $petsc_library),
        PetscErrorCode,
        (PetscFunctionList,),
        arg1,
    )
end

@for_petsc function PetscFunctionListPrintAll(::$UnionPetscLib)
    @chk ccall((:PetscFunctionListPrintAll, $petsc_library), PetscErrorCode, ())
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

@for_petsc function PetscBasename(::$UnionPetscLib, arg1)
    ccall((:PetscBasename, $petsc_library), Ptr{Cchar}, (Ptr{Cchar},), arg1)
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

@for_petsc function PetscStrcat(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscStrcat, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
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

@for_petsc function PetscStrtolower(::$UnionPetscLib, a)
    @chk ccall(
        (:PetscStrtolower, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar},),
        a,
    )
end

@for_petsc function PetscStrtoupper(::$UnionPetscLib, a)
    @chk ccall(
        (:PetscStrtoupper, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar},),
        a,
    )
end

@for_petsc function PetscStrlen(::$UnionPetscLib, s, len)
    @chk ccall(
        (:PetscStrlen, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Csize_t}),
        s,
        len,
    )
end

@for_petsc function PetscStrallocpy(::$UnionPetscLib, s, t)
    @chk ccall(
        (:PetscStrallocpy, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Ptr{Cchar}}),
        s,
        t,
    )
end

@for_petsc function PetscStrcmpNoError(::$UnionPetscLib, a, b, flg)
    ccall(
        (:PetscStrcmpNoError, $petsc_library),
        Cvoid,
        (Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}),
        a,
        b,
        flg,
    )
end

@for_petsc function PetscStrcmp(::$UnionPetscLib, a, b, flg)
    @chk ccall(
        (:PetscStrcmp, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}),
        a,
        b,
        flg,
    )
end

@for_petsc function PetscStrncpy(::$UnionPetscLib, s, t, n)
    @chk ccall(
        (:PetscStrncpy, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Csize_t),
        s,
        t,
        n,
    )
end

@for_petsc function PetscStrlcat(::$UnionPetscLib, s, t, n)
    @chk ccall(
        (:PetscStrlcat, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Csize_t),
        s,
        t,
        n,
    )
end

@for_petsc function PetscStrncmp(::$UnionPetscLib, a, b, n, t)
    @chk ccall(
        (:PetscStrncmp, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Csize_t, Ptr{PetscBool}),
        a,
        b,
        n,
        t,
    )
end

@for_petsc function PetscStrrstr(::$UnionPetscLib, a, b, tmp)
    @chk ccall(
        (:PetscStrrstr, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Ptr{Ptr{Cchar}}),
        a,
        b,
        tmp,
    )
end

@for_petsc function PetscStrstr(::$UnionPetscLib, haystack, needle, tmp)
    @chk ccall(
        (:PetscStrstr, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Ptr{Ptr{Cchar}}),
        haystack,
        needle,
        tmp,
    )
end

@for_petsc function PetscStrgrt(::$UnionPetscLib, a, b, t)
    @chk ccall(
        (:PetscStrgrt, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}),
        a,
        b,
        t,
    )
end

@for_petsc function PetscStrchr(::$UnionPetscLib, a, b, c)
    @chk ccall(
        (:PetscStrchr, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Cchar, Ptr{Ptr{Cchar}}),
        a,
        b,
        c,
    )
end

@for_petsc function PetscStrrchr(::$UnionPetscLib, a, b, c)
    @chk ccall(
        (:PetscStrrchr, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Cchar, Ptr{Ptr{Cchar}}),
        a,
        b,
        c,
    )
end

@for_petsc function PetscStrendswith(::$UnionPetscLib, a, b, flg)
    @chk ccall(
        (:PetscStrendswith, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}),
        a,
        b,
        flg,
    )
end

@for_petsc function PetscStrbeginswith(::$UnionPetscLib, a, b, flg)
    @chk ccall(
        (:PetscStrbeginswith, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}),
        a,
        b,
        flg,
    )
end

@for_petsc function PetscBTIndex_Internal(::$UnionPetscLib, index)
    ccall(
        (:PetscBTIndex_Internal, $petsc_library),
        Csize_t,
        ($PetscInt,),
        index,
    )
end

@for_petsc function PetscBTMask_Internal(::$UnionPetscLib, index)
    ccall((:PetscBTMask_Internal, $petsc_library), Cchar, ($PetscInt,), index)
end

@for_petsc function PetscBTLength(::$UnionPetscLib, m)
    ccall((:PetscBTLength, $petsc_library), Csize_t, ($PetscInt,), m)
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

@for_petsc function PetscBTCreate(::$UnionPetscLib, m, array)
    @chk ccall(
        (:PetscBTCreate, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{PetscBT}),
        m,
        array,
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

@for_petsc function PetscBTClear(::$UnionPetscLib, array, index)
    @chk ccall(
        (:PetscBTClear, $petsc_library),
        PetscErrorCode,
        (PetscBT, $PetscInt),
        array,
        index,
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

@for_petsc function PetscBTLookupClear(::$UnionPetscLib, array, index)
    ccall(
        (:PetscBTLookupClear, $petsc_library),
        Cchar,
        (PetscBT, $PetscInt),
        array,
        index,
    )
end

@for_petsc function PetscBTView(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscBTView, $petsc_library),
        PetscErrorCode,
        ($PetscInt, PetscBT, PetscViewer),
        arg1,
        arg2,
        arg3,
    )
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

mutable struct _n_PetscIntStack end

const PetscIntStack = Ptr{_n_PetscIntStack}

const PetscLogClass = Cint

mutable struct _p_PetscLogHandler end

const PetscLogHandler = Ptr{_p_PetscLogHandler}

const PetscLogHandlerType = Ptr{Cchar}

mutable struct _n_PetscLogRegistry end

const PetscLogRegistry = Ptr{_n_PetscLogRegistry}

mutable struct _n_PetscLogState
    registry::PetscLogRegistry
    active::PetscBT
    stage_stack::PetscIntStack
    current_stage::Cint
    bt_num_stages::Cint
    bt_num_events::Cint
    refct::Cint
    _n_PetscLogState() = new()
end

const PetscLogState = Ptr{_n_PetscLogState}

mutable struct PetscLogEventInfo
    name::Ptr{Cchar}
    classid::PetscClassId
    collective::PetscBool
    PetscLogEventInfo() = new()
end

mutable struct PetscLogClassInfo
    name::Ptr{Cchar}
    classid::PetscClassId
    PetscLogClassInfo() = new()
end

mutable struct _PetscLogStageInfo
    name::Ptr{Cchar}
    _PetscLogStageInfo() = new()
end

const PetscLogStageInfo = _PetscLogStageInfo

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

@for_petsc function PetscIntStackCreate(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscIntStackCreate, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscIntStack},),
        arg1,
    )
end

@for_petsc function PetscIntStackDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscIntStackDestroy, $petsc_library),
        PetscErrorCode,
        (PetscIntStack,),
        arg1,
    )
end

@for_petsc function PetscIntStackPush(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscIntStackPush, $petsc_library),
        PetscErrorCode,
        (PetscIntStack, Cint),
        arg1,
        arg2,
    )
end

@for_petsc function PetscIntStackPop(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscIntStackPop, $petsc_library),
        PetscErrorCode,
        (PetscIntStack, Ptr{Cint}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscIntStackTop(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscIntStackTop, $petsc_library),
        PetscErrorCode,
        (PetscIntStack, Ptr{Cint}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscIntStackEmpty(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscIntStackEmpty, $petsc_library),
        PetscErrorCode,
        (PetscIntStack, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogStateCreate(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogStateCreate, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscLogState},),
        arg1,
    )
end

@for_petsc function PetscLogStateDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogStateDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscLogState},),
        arg1,
    )
end

@for_petsc function PetscLogStateGetRegistry(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogStateGetRegistry, $petsc_library),
        PetscErrorCode,
        (PetscLogState, Ptr{PetscLogRegistry}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogStateClassRegister(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscLogStateClassRegister, $petsc_library),
        PetscErrorCode,
        (PetscLogState, Ptr{Cchar}, PetscClassId, Ptr{PetscLogStage}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscLogStateClassSetActive(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscLogStateClassSetActive, $petsc_library),
        PetscErrorCode,
        (PetscLogState, PetscLogStage, PetscClassId, PetscBool),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscLogStateClassSetActiveAll(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscLogStateClassSetActiveAll, $petsc_library),
        PetscErrorCode,
        (PetscLogState, PetscClassId, PetscBool),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogStateStageRegister(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscLogStateStageRegister, $petsc_library),
        PetscErrorCode,
        (PetscLogState, Ptr{Cchar}, Ptr{PetscLogStage}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogStateStagePush(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogStateStagePush, $petsc_library),
        PetscErrorCode,
        (PetscLogState, PetscLogStage),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogStateStagePop(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogStateStagePop, $petsc_library),
        PetscErrorCode,
        (PetscLogState,),
        arg1,
    )
end

@for_petsc function PetscLogStateStageSetActive(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscLogStateStageSetActive, $petsc_library),
        PetscErrorCode,
        (PetscLogState, PetscLogStage, PetscBool),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogStateStageGetActive(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscLogStateStageGetActive, $petsc_library),
        PetscErrorCode,
        (PetscLogState, PetscLogStage, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogStateGetCurrentStage(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogStateGetCurrentStage, $petsc_library),
        PetscErrorCode,
        (PetscLogState, Ptr{PetscLogStage}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogStateEventRegister(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscLogStateEventRegister, $petsc_library),
        PetscErrorCode,
        (PetscLogState, Ptr{Cchar}, PetscClassId, Ptr{PetscLogEvent}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscLogStateEventSetCollective(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscLogStateEventSetCollective, $petsc_library),
        PetscErrorCode,
        (PetscLogState, PetscLogEvent, PetscBool),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogStateEventSetActive(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscLogStateEventSetActive, $petsc_library),
        PetscErrorCode,
        (PetscLogState, PetscLogStage, PetscLogEvent, PetscBool),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscLogStateEventSetActiveAll(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscLogStateEventSetActiveAll, $petsc_library),
        PetscErrorCode,
        (PetscLogState, PetscLogEvent, PetscBool),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogStateEventGetActive(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscLogStateEventGetActive, $petsc_library),
        PetscErrorCode,
        (PetscLogState, PetscLogStage, PetscLogEvent, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscLogStateGetEventFromName(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscLogStateGetEventFromName, $petsc_library),
        PetscErrorCode,
        (PetscLogState, Ptr{Cchar}, Ptr{PetscLogEvent}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogStateGetStageFromName(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscLogStateGetStageFromName, $petsc_library),
        PetscErrorCode,
        (PetscLogState, Ptr{Cchar}, Ptr{PetscLogStage}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogStateGetClassFromName(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscLogStateGetClassFromName, $petsc_library),
        PetscErrorCode,
        (PetscLogState, Ptr{Cchar}, Ptr{PetscLogClass}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogStateGetClassFromClassId(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscLogStateGetClassFromClassId, $petsc_library),
        PetscErrorCode,
        (PetscLogState, PetscClassId, Ptr{PetscLogClass}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogStateGetNumEvents(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogStateGetNumEvents, $petsc_library),
        PetscErrorCode,
        (PetscLogState, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogStateGetNumStages(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogStateGetNumStages, $petsc_library),
        PetscErrorCode,
        (PetscLogState, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogStateGetNumClasses(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogStateGetNumClasses, $petsc_library),
        PetscErrorCode,
        (PetscLogState, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogStateEventGetInfo(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscLogStateEventGetInfo, $petsc_library),
        PetscErrorCode,
        (PetscLogState, PetscLogEvent, Ptr{PetscLogEventInfo}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogStateStageGetInfo(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscLogStateStageGetInfo, $petsc_library),
        PetscErrorCode,
        (PetscLogState, PetscLogStage, Ptr{PetscLogStageInfo}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogStateClassGetInfo(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscLogStateClassGetInfo, $petsc_library),
        PetscErrorCode,
        (PetscLogState, PetscLogClass, Ptr{PetscLogClassInfo}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogHandlerRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogHandlerRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogHandlerCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogHandlerCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscLogHandler}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogHandlerSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogHandlerSetType, $petsc_library),
        PetscErrorCode,
        (PetscLogHandler, PetscLogHandlerType),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogHandlerGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogHandlerGetType, $petsc_library),
        PetscErrorCode,
        (PetscLogHandler, Ptr{PetscLogHandlerType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogHandlerDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogHandlerDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscLogHandler},),
        arg1,
    )
end

@for_petsc function PetscLogHandlerSetState(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogHandlerSetState, $petsc_library),
        PetscErrorCode,
        (PetscLogHandler, PetscLogState),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogHandlerGetState(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogHandlerGetState, $petsc_library),
        PetscErrorCode,
        (PetscLogHandler, Ptr{PetscLogState}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogHandlerEventBegin(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscLogHandlerEventBegin, $petsc_library),
        PetscErrorCode,
        (
            PetscLogHandler,
            PetscLogEvent,
            PetscObject,
            PetscObject,
            PetscObject,
            PetscObject,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscLogHandlerEventEnd(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscLogHandlerEventEnd, $petsc_library),
        PetscErrorCode,
        (
            PetscLogHandler,
            PetscLogEvent,
            PetscObject,
            PetscObject,
            PetscObject,
            PetscObject,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscLogHandlerEventSync(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscLogHandlerEventSync, $petsc_library),
        PetscErrorCode,
        (PetscLogHandler, PetscLogEvent, MPI_Comm),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogHandlerObjectCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogHandlerObjectCreate, $petsc_library),
        PetscErrorCode,
        (PetscLogHandler, PetscObject),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogHandlerObjectDestroy(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogHandlerObjectDestroy, $petsc_library),
        PetscErrorCode,
        (PetscLogHandler, PetscObject),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogHandlerStagePush(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogHandlerStagePush, $petsc_library),
        PetscErrorCode,
        (PetscLogHandler, PetscLogStage),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogHandlerStagePop(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogHandlerStagePop, $petsc_library),
        PetscErrorCode,
        (PetscLogHandler, PetscLogStage),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogHandlerView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogHandlerView, $petsc_library),
        PetscErrorCode,
        (PetscLogHandler, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogHandlerGetEventPerfInfo(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscLogHandlerGetEventPerfInfo, $petsc_library),
        PetscErrorCode,
        (
            PetscLogHandler,
            PetscLogStage,
            PetscLogEvent,
            Ptr{Ptr{PetscEventPerfInfo}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscLogHandlerGetStagePerfInfo(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscLogHandlerGetStagePerfInfo, $petsc_library),
        PetscErrorCode,
        (PetscLogHandler, PetscLogStage, Ptr{Ptr{PetscEventPerfInfo}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogHandlerSetLogActions(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogHandlerSetLogActions, $petsc_library),
        PetscErrorCode,
        (PetscLogHandler, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogHandlerSetLogObjects(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogHandlerSetLogObjects, $petsc_library),
        PetscErrorCode,
        (PetscLogHandler, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogHandlerGetNumObjects(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogHandlerGetNumObjects, $petsc_library),
        PetscErrorCode,
        (PetscLogHandler, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogHandlerEventDeactivatePush(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscLogHandlerEventDeactivatePush, $petsc_library),
        PetscErrorCode,
        (PetscLogHandler, PetscLogStage, PetscLogEvent),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogHandlerEventDeactivatePop(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscLogHandlerEventDeactivatePop, $petsc_library),
        PetscErrorCode,
        (PetscLogHandler, PetscLogStage, PetscLogEvent),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogHandlerEventsPause(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogHandlerEventsPause, $petsc_library),
        PetscErrorCode,
        (PetscLogHandler,),
        arg1,
    )
end

@for_petsc function PetscLogHandlerEventsResume(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogHandlerEventsResume, $petsc_library),
        PetscErrorCode,
        (PetscLogHandler,),
        arg1,
    )
end

@for_petsc function PetscLogHandlerDump(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogHandlerDump, $petsc_library),
        PetscErrorCode,
        (PetscLogHandler, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogHandlerStageSetVisible(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscLogHandlerStageSetVisible, $petsc_library),
        PetscErrorCode,
        (PetscLogHandler, PetscLogStage, PetscBool),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogHandlerStageGetVisible(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscLogHandlerStageGetVisible, $petsc_library),
        PetscErrorCode,
        (PetscLogHandler, PetscLogStage, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogHandlerCreateTrace(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscLogHandlerCreateTrace, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Libc.FILE}, Ptr{PetscLogHandler}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscLogHandlerCreateLegacy(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscLogHandlerCreateLegacy, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{PetscLogHandler},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

mutable struct _n_PetscLogHandlerHot
    handler::PetscLogHandler
    eventBegin::Ptr{Cvoid}
    eventEnd::Ptr{Cvoid}
    eventSync::Ptr{Cvoid}
    objectCreate::Ptr{Cvoid}
    objectDestroy::Ptr{Cvoid}
    _n_PetscLogHandlerHot() = new()
end

const PetscLogHandlerHot = _n_PetscLogHandlerHot

@for_petsc function PetscLogObjectMemory(::$UnionPetscLib, o, m)
    @chk ccall(
        (:PetscLogObjectMemory, $petsc_library),
        PetscErrorCode,
        (PetscObject, PetscLogDouble),
        o,
        m,
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

@for_petsc function PetscLogMPEBegin(::$UnionPetscLib)
    @chk ccall((:PetscLogMPEBegin, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscLogPerfstubsBegin(::$UnionPetscLib)
    @chk ccall((:PetscLogPerfstubsBegin, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscLogLegacyCallbacksBegin(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscLogLegacyCallbacksBegin, $petsc_library),
        PetscErrorCode,
        (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
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

@for_petsc function PetscLogMPEDump(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogMPEDump, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar},),
        arg1,
    )
end

@for_petsc function PetscLogGetState(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogGetState, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscLogState},),
        arg1,
    )
end

@for_petsc function PetscLogGetDefaultHandler(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogGetDefaultHandler, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscLogHandler},),
        arg1,
    )
end

@for_petsc function PetscLogHandlerStart(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogHandlerStart, $petsc_library),
        PetscErrorCode,
        (PetscLogHandler,),
        arg1,
    )
end

@for_petsc function PetscLogHandlerStop(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogHandlerStop, $petsc_library),
        PetscErrorCode,
        (PetscLogHandler,),
        arg1,
    )
end

@for_petsc function PetscLogIsActive(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscLogIsActive, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscBool},),
        arg1,
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

@for_petsc function PetscLogStageGetName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogStageGetName, $petsc_library),
        PetscErrorCode,
        (PetscLogEvent, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogStageGetPerfInfo(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogStageGetPerfInfo, $petsc_library),
        PetscErrorCode,
        (PetscLogStage, Ptr{PetscEventPerfInfo}),
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

@for_petsc function PetscLogEventGetName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogEventGetName, $petsc_library),
        PetscErrorCode,
        (PetscLogEvent, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogEventGetPerfInfo(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscLogEventGetPerfInfo, $petsc_library),
        PetscErrorCode,
        (PetscLogStage, PetscLogEvent, Ptr{PetscEventPerfInfo}),
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

@for_petsc function PetscLogEventsPause(::$UnionPetscLib)
    @chk ccall((:PetscLogEventsPause, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscLogEventsResume(::$UnionPetscLib)
    @chk ccall((:PetscLogEventsResume, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscLogClassGetClassId(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogClassGetClassId, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{PetscClassId}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogClassIdGetName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscLogClassIdGetName, $petsc_library),
        PetscErrorCode,
        (PetscClassId, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscLogEventSync(::$UnionPetscLib, e, comm)
    @chk ccall(
        (:PetscLogEventSync, $petsc_library),
        PetscErrorCode,
        (PetscLogEvent, MPI_Comm),
        e,
        comm,
    )
end

@for_petsc function PetscLogObjectCreate(::$UnionPetscLib, o)
    @chk ccall(
        (:PetscLogObjectCreate, $petsc_library),
        PetscErrorCode,
        (PetscObject,),
        o,
    )
end

@for_petsc function PetscLogObjectDestroy(::$UnionPetscLib, o)
    @chk ccall(
        (:PetscLogObjectDestroy, $petsc_library),
        PetscErrorCode,
        (PetscObject,),
        o,
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

mutable struct PetscClassRegInfo
    name::Ptr{Cchar}
    classid::PetscClassId
    PetscClassRegInfo() = new()
end

struct _n_PetscClassRegLog
    numClasses::Cint
    maxClasses::Cint
    classInfo::Ptr{PetscClassRegInfo}
end

const PetscClassRegLog = Ptr{_n_PetscClassRegLog}

mutable struct PetscClassPerfInfo
    id::PetscClassId
    creations::Cint
    destructions::Cint
    mem::PetscLogDouble
    descMem::PetscLogDouble
    PetscClassPerfInfo() = new()
end

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
    timer::Ptr{Cvoid}
    PetscEventRegInfo() = new()
end

struct _n_PetscEventRegLog
    numEvents::Cint
    maxEvents::Cint
    eventInfo::Ptr{PetscEventRegInfo}
end

const PetscEventRegLog = Ptr{_n_PetscEventRegLog}

mutable struct _n_PetscEventPerfLog
    numEvents::Cint
    maxEvents::Cint
    eventInfo::Ptr{PetscEventPerfInfo}
    _n_PetscEventPerfLog() = new()
end

const PetscEventPerfLog = Ptr{_n_PetscEventPerfLog}

struct _PetscStageInfo
    name::Ptr{Cchar}
    used::PetscBool
    perfInfo::PetscEventPerfInfo
    classLog::PetscClassPerfLog
    timer::Ptr{Cvoid}
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

@for_petsc function PetscLogGetStageLog(::$UnionPetscLib, s)
    @chk ccall(
        (:PetscLogGetStageLog, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscStageLog},),
        s,
    )
end

@for_petsc function PetscStageLogGetCurrent(::$UnionPetscLib, a, b)
    @chk ccall(
        (:PetscStageLogGetCurrent, $petsc_library),
        PetscErrorCode,
        (PetscStageLog, Ptr{Cint}),
        a,
        b,
    )
end

@for_petsc function PetscStageLogGetEventPerfLog(::$UnionPetscLib, a, b, c)
    @chk ccall(
        (:PetscStageLogGetEventPerfLog, $petsc_library),
        PetscErrorCode,
        (PetscStageLog, Cint, Ptr{PetscEventPerfLog}),
        a,
        b,
        c,
    )
end

@for_petsc function PetscLogPushCurrentEvent_Internal(::$UnionPetscLib, e)
    @chk ccall(
        (:PetscLogPushCurrentEvent_Internal, $petsc_library),
        PetscErrorCode,
        (PetscLogEvent,),
        e,
    )
end

@for_petsc function PetscLogPopCurrentEvent_Internal(::$UnionPetscLib)
    @chk ccall(
        (:PetscLogPopCurrentEvent_Internal, $petsc_library),
        PetscErrorCode,
        (),
    )
end

@for_petsc function PetscLogAllBegin(::$UnionPetscLib)
    @chk ccall((:PetscLogAllBegin, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscLogSet(::$UnionPetscLib, a, b)
    @chk ccall(
        (:PetscLogSet, $petsc_library),
        PetscErrorCode,
        (Ptr{Cvoid}, Ptr{Cvoid}),
        a,
        b,
    )
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

@for_petsc function PetscFFlush(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscFFlush, $petsc_library),
        PetscErrorCode,
        (Ptr{Libc.FILE},),
        arg1,
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

@for_petsc function PetscIntCast(::$UnionPetscLib, a, b)
    @chk ccall(
        (:PetscIntCast, $petsc_library),
        PetscErrorCode,
        (PetscInt64, Ptr{$PetscInt}),
        a,
        b,
    )
end

@for_petsc function PetscCountCast(::$UnionPetscLib, a, b)
    @chk ccall(
        (:PetscCountCast, $petsc_library),
        PetscErrorCode,
        (PetscCount, Ptr{$PetscInt}),
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

@for_petsc function PetscCuBLASIntCast(::$UnionPetscLib, a, b)
    @chk ccall(
        (:PetscCuBLASIntCast, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{PetscCuBLASInt}),
        a,
        b,
    )
end

@for_petsc function PetscHipBLASIntCast(::$UnionPetscLib, a, b)
    @chk ccall(
        (:PetscHipBLASIntCast, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{PetscHipBLASInt}),
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

@for_petsc function PetscSortedInt64(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSortedInt64, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{PetscInt64}, Ptr{PetscBool}),
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

@for_petsc function PetscSortInt64(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSortInt64, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{PetscInt64}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSortCount(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSortCount, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{PetscCount}),
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

@for_petsc function PetscSortedCheckDupsInt(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSortedCheckDupsInt, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{$PetscInt}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
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

@for_petsc function PetscSortIntWithCountArray(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSortIntWithCountArray, $petsc_library),
        PetscErrorCode,
        (PetscCount, Ptr{$PetscInt}, Ptr{PetscCount}),
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

@for_petsc function PetscSortIntWithIntCountArrayPair(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSortIntWithIntCountArrayPair, $petsc_library),
        PetscErrorCode,
        (PetscCount, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{PetscCount}),
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

@for_petsc function PetscRandomFinalizePackage(::$UnionPetscLib)
    @chk ccall(
        (:PetscRandomFinalizePackage, $petsc_library),
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

@for_petsc function PetscBoxUpload(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscBoxUpload, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
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

@for_petsc function PetscGlobusAuthorize(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscGlobusAuthorize, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Csize_t),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscGlobusUpload(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscGlobusUpload, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
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

@for_petsc function PCMPIServerBegin(::$UnionPetscLib)
    @chk ccall((:PCMPIServerBegin, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PCMPIServerEnd(::$UnionPetscLib)
    @chk ccall((:PCMPIServerEnd, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PCMPICommsDestroy(::$UnionPetscLib)
    @chk ccall((:PCMPICommsDestroy, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscBLASSetNumThreads(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscBLASSetNumThreads, $petsc_library),
        PetscErrorCode,
        ($PetscInt,),
        arg1,
    )
end

@for_petsc function PetscBLASGetNumThreads(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscBLASGetNumThreads, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscInt},),
        arg1,
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

@for_petsc function PetscDrawSetVisible(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawSetVisible, $petsc_library),
        PetscErrorCode,
        (PetscDraw, PetscBool),
        arg1,
        arg2,
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

@for_petsc function PetscDrawSPGetDimension(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDrawSPGetDimension, $petsc_library),
        PetscErrorCode,
        (PetscDrawSP, Ptr{Cint}),
        arg1,
        arg2,
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

@for_petsc function PetscDrawSPAddPointColorized(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDrawSPAddPointColorized, $petsc_library),
        PetscErrorCode,
        (PetscDrawSP, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
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

@for_petsc function PetscMemTypeToString(::$UnionPetscLib, mtype)
    ccall(
        (:PetscMemTypeToString, $petsc_library),
        Ptr{Cchar},
        (PetscMemType,),
        mtype,
    )
end

@enum PetscOffloadMask::UInt32 begin
    PETSC_OFFLOAD_UNALLOCATED = 0
    PETSC_OFFLOAD_CPU = 1
    PETSC_OFFLOAD_GPU = 2
    PETSC_OFFLOAD_BOTH = 3
    PETSC_OFFLOAD_VECKOKKOS = 256
    # PETSC_OFFLOAD_KOKKOS = 256
end

@for_petsc function PetscOffloadMaskToString(::$UnionPetscLib, mask)
    ccall(
        (:PetscOffloadMaskToString, $petsc_library),
        Ptr{Cchar},
        (PetscOffloadMask,),
        mask,
    )
end

@for_petsc function PetscOffloadMaskToMemType(::$UnionPetscLib, mask)
    ccall(
        (:PetscOffloadMaskToMemType, $petsc_library),
        PetscMemType,
        (PetscOffloadMask,),
        mask,
    )
end

@enum PetscDeviceInitType::UInt32 begin
    PETSC_DEVICE_INIT_NONE = 0
    PETSC_DEVICE_INIT_LAZY = 1
    PETSC_DEVICE_INIT_EAGER = 2
end

@enum PetscDeviceType::UInt32 begin
    PETSC_DEVICE_HOST = 0
    PETSC_DEVICE_CUDA = 1
    PETSC_DEVICE_HIP = 2
    PETSC_DEVICE_SYCL = 3
    PETSC_DEVICE_MAX = 4
end

@enum PetscDeviceAttribute::UInt32 begin
    PETSC_DEVICE_ATTR_SIZE_T_SHARED_MEM_PER_BLOCK = 0
    PETSC_DEVICE_ATTR_MAX = 1
end

mutable struct _n_PetscDevice end

const PetscDevice = Ptr{_n_PetscDevice}

@enum PetscStreamType::UInt32 begin
    PETSC_STREAM_DEFAULT = 0
    PETSC_STREAM_NONBLOCKING = 1
    PETSC_STREAM_DEFAULT_WITH_BARRIER = 2
    PETSC_STREAM_NONBLOCKING_WITH_BARRIER = 3
    PETSC_STREAM_MAX = 4
end

@enum PetscDeviceContextJoinMode::UInt32 begin
    PETSC_DEVICE_CONTEXT_JOIN_DESTROY = 0
    PETSC_DEVICE_CONTEXT_JOIN_SYNC = 1
    PETSC_DEVICE_CONTEXT_JOIN_NO_SYNC = 2
end

@enum PetscDeviceCopyMode::UInt32 begin
    PETSC_DEVICE_COPY_HTOH = 0
    PETSC_DEVICE_COPY_DTOH = 1
    PETSC_DEVICE_COPY_HTOD = 2
    PETSC_DEVICE_COPY_DTOD = 3
    PETSC_DEVICE_COPY_AUTO = 4
end

@for_petsc function PetscOffloadMaskToDeviceCopyMode(
    ::$UnionPetscLib,
    dest,
    src,
)
    ccall(
        (:PetscOffloadMaskToDeviceCopyMode, $petsc_library),
        PetscDeviceCopyMode,
        (PetscOffloadMask, PetscOffloadMask),
        dest,
        src,
    )
end

@for_petsc function PetscMemTypeToDeviceCopyMode(::$UnionPetscLib, dest, src)
    ccall(
        (:PetscMemTypeToDeviceCopyMode, $petsc_library),
        PetscDeviceCopyMode,
        (PetscMemType, PetscMemType),
        dest,
        src,
    )
end

@enum PetscMemoryAccessMode::UInt32 begin
    PETSC_MEMORY_ACCESS_READ = 1
    PETSC_MEMORY_ACCESS_WRITE = 2
    PETSC_MEMORY_ACCESS_READ_WRITE = 3
end

@for_petsc function PetscMemoryAccessModeToString(::$UnionPetscLib, mode)
    ccall(
        (:PetscMemoryAccessModeToString, $petsc_library),
        Ptr{Cchar},
        (PetscMemoryAccessMode,),
        mode,
    )
end

@for_petsc function PetscDeviceInitializePackage(::$UnionPetscLib)
    @chk ccall(
        (:PetscDeviceInitializePackage, $petsc_library),
        PetscErrorCode,
        (),
    )
end

@for_petsc function PetscDeviceFinalizePackage(::$UnionPetscLib)
    @chk ccall(
        (:PetscDeviceFinalizePackage, $petsc_library),
        PetscErrorCode,
        (),
    )
end

@for_petsc function PetscGetMemType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscGetMemType, $petsc_library),
        PetscErrorCode,
        (Ptr{Cvoid}, Ptr{PetscMemType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDeviceCreate(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDeviceCreate, $petsc_library),
        PetscErrorCode,
        (PetscDeviceType, $PetscInt, Ptr{PetscDevice}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDeviceDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDeviceDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscDevice},),
        arg1,
    )
end

@for_petsc function PetscDeviceConfigure(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDeviceConfigure, $petsc_library),
        PetscErrorCode,
        (PetscDevice,),
        arg1,
    )
end

@for_petsc function PetscDeviceView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDeviceView, $petsc_library),
        PetscErrorCode,
        (PetscDevice, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDeviceGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDeviceGetType, $petsc_library),
        PetscErrorCode,
        (PetscDevice, Ptr{PetscDeviceType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDeviceGetDeviceId(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDeviceGetDeviceId, $petsc_library),
        PetscErrorCode,
        (PetscDevice, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PETSC_DEVICE_DEFAULT(::$UnionPetscLib)
    ccall((:PETSC_DEVICE_DEFAULT, $petsc_library), PetscDeviceType, ())
end

@for_petsc function PetscDeviceSetDefaultDeviceType(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDeviceSetDefaultDeviceType, $petsc_library),
        PetscErrorCode,
        (PetscDeviceType,),
        arg1,
    )
end

@for_petsc function PetscDeviceInitialize(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDeviceInitialize, $petsc_library),
        PetscErrorCode,
        (PetscDeviceType,),
        arg1,
    )
end

@for_petsc function PetscDeviceInitialized(::$UnionPetscLib, arg1)
    ccall(
        (:PetscDeviceInitialized, $petsc_library),
        PetscBool,
        (PetscDeviceType,),
        arg1,
    )
end

@for_petsc function PetscDeviceContextCreate(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDeviceContextCreate, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscDeviceContext},),
        arg1,
    )
end

@for_petsc function PetscDeviceContextDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDeviceContextDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscDeviceContext},),
        arg1,
    )
end

@for_petsc function PetscDeviceContextSetStreamType(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDeviceContextSetStreamType, $petsc_library),
        PetscErrorCode,
        (PetscDeviceContext, PetscStreamType),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDeviceContextGetStreamType(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDeviceContextGetStreamType, $petsc_library),
        PetscErrorCode,
        (PetscDeviceContext, Ptr{PetscStreamType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDeviceContextSetDevice(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDeviceContextSetDevice, $petsc_library),
        PetscErrorCode,
        (PetscDeviceContext, PetscDevice),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDeviceContextGetDevice(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDeviceContextGetDevice, $petsc_library),
        PetscErrorCode,
        (PetscDeviceContext, Ptr{PetscDevice}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDeviceContextGetDeviceType(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDeviceContextGetDeviceType, $petsc_library),
        PetscErrorCode,
        (PetscDeviceContext, Ptr{PetscDeviceType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDeviceContextSetUp(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDeviceContextSetUp, $petsc_library),
        PetscErrorCode,
        (PetscDeviceContext,),
        arg1,
    )
end

@for_petsc function PetscDeviceContextDuplicate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDeviceContextDuplicate, $petsc_library),
        PetscErrorCode,
        (PetscDeviceContext, Ptr{PetscDeviceContext}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDeviceContextQueryIdle(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDeviceContextQueryIdle, $petsc_library),
        PetscErrorCode,
        (PetscDeviceContext, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDeviceContextWaitForContext(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDeviceContextWaitForContext, $petsc_library),
        PetscErrorCode,
        (PetscDeviceContext, PetscDeviceContext),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDeviceContextForkWithStreamType(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDeviceContextForkWithStreamType, $petsc_library),
        PetscErrorCode,
        (
            PetscDeviceContext,
            PetscStreamType,
            $PetscInt,
            Ptr{Ptr{PetscDeviceContext}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDeviceContextFork(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDeviceContextFork, $petsc_library),
        PetscErrorCode,
        (PetscDeviceContext, $PetscInt, Ptr{Ptr{PetscDeviceContext}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDeviceContextJoin(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDeviceContextJoin, $petsc_library),
        PetscErrorCode,
        (
            PetscDeviceContext,
            $PetscInt,
            PetscDeviceContextJoinMode,
            Ptr{Ptr{PetscDeviceContext}},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDeviceContextSynchronize(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDeviceContextSynchronize, $petsc_library),
        PetscErrorCode,
        (PetscDeviceContext,),
        arg1,
    )
end

@for_petsc function PetscDeviceContextSetFromOptions(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDeviceContextSetFromOptions, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, PetscDeviceContext),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDeviceContextView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscDeviceContextView, $petsc_library),
        PetscErrorCode,
        (PetscDeviceContext, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDeviceContextViewFromOptions(
    ::$UnionPetscLib,
    arg1,
    arg2,
    name,
)
    @chk ccall(
        (:PetscDeviceContextViewFromOptions, $petsc_library),
        PetscErrorCode,
        (PetscDeviceContext, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        name,
    )
end

@for_petsc function PetscDeviceContextGetCurrentContext(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDeviceContextGetCurrentContext, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscDeviceContext},),
        arg1,
    )
end

@for_petsc function PetscDeviceContextSetCurrentContext(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscDeviceContextSetCurrentContext, $petsc_library),
        PetscErrorCode,
        (PetscDeviceContext,),
        arg1,
    )
end

@for_petsc function PetscDeviceContextGetStreamHandle(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDeviceContextGetStreamHandle, $petsc_library),
        PetscErrorCode,
        (PetscDeviceContext, Ptr{Ptr{Cvoid}}),
        arg1,
        arg2,
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

mutable struct _n_ISColoring end

const ISColoring = Ptr{_n_ISColoring}

@for_petsc function ISInitializePackage(::$UnionPetscLib)
    @chk ccall((:ISInitializePackage, $petsc_library), PetscErrorCode, ())
end

@for_petsc function ISFinalizePackage(::$UnionPetscLib)
    @chk ccall((:ISFinalizePackage, $petsc_library), PetscErrorCode, ())
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

@for_petsc function ISRegisterAll(::$UnionPetscLib)
    @chk ccall((:ISRegisterAll, $petsc_library), PetscErrorCode, ())
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

@for_petsc function ISGetNonlocalIS(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISGetNonlocalIS, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{IS}),
        arg1,
        arg2,
    )
end

@for_petsc function ISRestoreNonlocalIS(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISRestoreNonlocalIS, $petsc_library),
        PetscErrorCode,
        (IS, Ptr{IS}),
        arg1,
        arg2,
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

@for_petsc function ISShift(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:ISShift, $petsc_library),
        PetscErrorCode,
        (IS, $PetscInt, IS),
        arg1,
        arg2,
        arg3,
    )
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

@for_petsc function ISGeneralSetIndicesFromMask(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:ISGeneralSetIndicesFromMask, $petsc_library),
        PetscErrorCode,
        (IS, $PetscInt, $PetscInt, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
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

@for_petsc function ISLocalToGlobalMappingGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISLocalToGlobalMappingGetType, $petsc_library),
        PetscErrorCode,
        (ISLocalToGlobalMapping, Ptr{ISLocalToGlobalMappingType}),
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

@for_petsc function ISLocalToGlobalMappingLoad(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:ISLocalToGlobalMappingLoad, $petsc_library),
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

@for_petsc function ISLocalToGlobalMappingDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:ISLocalToGlobalMappingDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{ISLocalToGlobalMapping},),
        arg1,
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

@for_petsc function ISLocalToGlobalMappingGetBlockNodeInfo(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:ISLocalToGlobalMappingGetBlockNodeInfo, $petsc_library),
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

@for_petsc function ISLocalToGlobalMappingRestoreBlockNodeInfo(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:ISLocalToGlobalMappingRestoreBlockNodeInfo, $petsc_library),
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

@for_petsc function ISLocalToGlobalMappingGetBlockMultiLeavesSF(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:ISLocalToGlobalMappingGetBlockMultiLeavesSF, $petsc_library),
        PetscErrorCode,
        (ISLocalToGlobalMapping, Ptr{PetscSF}),
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
    n,
    bs,
    imax,
    is_in,
    is_out,
)
    @chk ccall(
        (:ISCompressIndicesSorted, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, $PetscInt, Ptr{IS}, Ptr{IS}),
        n,
        bs,
        imax,
        is_in,
        is_out,
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

const PetscViewerType = Ptr{Cchar}

@for_petsc function PetscViewerInitializePackage(::$UnionPetscLib)
    @chk ccall(
        (:PetscViewerInitializePackage, $petsc_library),
        PetscErrorCode,
        (),
    )
end

@for_petsc function PetscViewerFinalizePackage(::$UnionPetscLib)
    @chk ccall(
        (:PetscViewerFinalizePackage, $petsc_library),
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

@for_petsc function PetscOptionsGetViewers(
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
        (:PetscOptionsGetViewers, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            PetscOptions,
            Ptr{Cchar},
            Ptr{Cchar},
            Ptr{$PetscInt},
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
        arg8,
    )
end

@for_petsc function PetscOptionsRestoreViewer(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscOptionsRestoreViewer, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscViewer},),
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
            PetscInt64,
            PetscInt64,
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
            PetscInt64,
            PetscInt64,
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

@for_petsc function PetscViewerSiloClearName(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscViewerSiloClearName, $petsc_library),
        PetscErrorCode,
        (PetscViewer,),
        arg1,
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

@for_petsc function PetscViewerASCIIGetStdout(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerASCIIGetStdout, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{PetscViewer}),
        arg1,
        arg2,
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

@enum ScatterMode::UInt32 begin
    SCATTER_FORWARD = 0
    SCATTER_REVERSE = 1
    SCATTER_FORWARD_LOCAL = 2
    SCATTER_REVERSE_LOCAL = 3
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

@for_petsc function VecCreateFromOptions(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:VecCreateFromOptions, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, $PetscInt, $PetscInt, $PetscInt, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
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

@enum ReductionType::UInt32 begin
    REDUCTION_SUM_REALPART = 10
    REDUCTION_MEAN_REALPART = 11
    REDUCTION_SUM_IMAGINARYPART = 12
    REDUCTION_MEAN_IMAGINARYPART = 13
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

@for_petsc function VecMean(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecMean, $petsc_library),
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

@for_petsc function VecMAXPBY(::$UnionPetscLib, arg1, arg2, arg3, arg4, arg5)
    @chk ccall(
        (:VecMAXPBY, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Ptr{$PetscScalar}, $PetscScalar, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
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

@for_petsc function VecStrideSumAll(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecStrideSumAll, $petsc_library),
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

@for_petsc function VecStrideSum(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecStrideSum, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscInt, Ptr{$PetscScalar}),
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

@for_petsc function VecSetPreallocationCOO(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecSetPreallocationCOO, $petsc_library),
        PetscErrorCode,
        (Vec, PetscCount, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecSetPreallocationCOOLocal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:VecSetPreallocationCOOLocal, $petsc_library),
        PetscErrorCode,
        (Vec, PetscCount, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecSetValuesCOO(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecSetValuesCOO, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscScalar}, InsertMode),
        arg1,
        arg2,
        arg3,
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

@for_petsc function VecRegisterAll(::$UnionPetscLib)
    @chk ccall((:VecRegisterAll, $petsc_library), PetscErrorCode, ())
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

@for_petsc function VecViewNative(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecViewNative, $petsc_library),
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

@for_petsc function VecErrorWeightedNorms(
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
        (:VecErrorWeightedNorms, $petsc_library),
        PetscErrorCode,
        (
            Vec,
            Vec,
            Vec,
            NormType,
            $PetscReal,
            Vec,
            $PetscReal,
            Vec,
            $PetscReal,
            Ptr{$PetscReal},
            Ptr{$PetscInt},
            Ptr{$PetscReal},
            Ptr{$PetscInt},
            Ptr{$PetscReal},
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
        arg11,
        arg12,
        arg13,
        arg14,
        arg15,
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

@for_petsc function VecBoundToCPU(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecBoundToCPU, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function VecSetBindingPropagates(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecSetBindingPropagates, $petsc_library),
        PetscErrorCode,
        (Vec, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function VecGetBindingPropagates(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecGetBindingPropagates, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{PetscBool}),
        arg1,
        arg2,
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

@for_petsc function VecCreateLocalVector(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecCreateLocalVector, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Vec}),
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

@for_petsc function VecGetArrayWriteAndMemType(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:VecGetArrayWriteAndMemType, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}, Ptr{PetscMemType}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function VecRestoreArrayWriteAndMemType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecRestoreArrayWriteAndMemType, $petsc_library),
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

@enum VecOperation::UInt32 begin
    VECOP_DUPLICATE = 0
    VECOP_SET = 10
    VECOP_VIEW = 33
    VECOP_LOAD = 41
    VECOP_VIEWNATIVE = 69
    VECOP_LOADNATIVE = 70
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

@for_petsc function VecGhostGetGhostIS(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecGhostGetGhostIS, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{IS}),
        arg1,
        arg2,
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

@for_petsc function VecISShift(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:VecISShift, $petsc_library),
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

@for_petsc function VecFilter(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:VecFilter, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function VecChop(::$UnionPetscLib, v, tol)
    @chk ccall(
        (:VecChop, $petsc_library),
        PetscErrorCode,
        (Vec, $PetscReal),
        v,
        tol,
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

@for_petsc function VecTaggerRegisterAll(::$UnionPetscLib)
    @chk ccall((:VecTaggerRegisterAll, $petsc_library), PetscErrorCode, ())
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

@for_petsc function VecTaggerComputeIS(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:VecTaggerComputeIS, $petsc_library),
        PetscErrorCode,
        (VecTagger, Vec, Ptr{IS}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
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

@for_petsc function VecCreateSeqCUDA(::$UnionPetscLib, a, b, c)
    @chk ccall(
        (:VecCreateSeqCUDA, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, Ptr{Vec}),
        a,
        b,
        c,
    )
end

@for_petsc function VecCreateSeqHIP(::$UnionPetscLib, a, b, c)
    @chk ccall(
        (:VecCreateSeqHIP, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, Ptr{Vec}),
        a,
        b,
        c,
    )
end

@for_petsc function VecCreateSeqCUDAWithArray(::$UnionPetscLib, a, b, c, d, e)
    @chk ccall(
        (:VecCreateSeqCUDAWithArray, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{Vec}),
        a,
        b,
        c,
        d,
        e,
    )
end

@for_petsc function VecCreateSeqHIPWithArray(::$UnionPetscLib, a, b, c, d, e)
    @chk ccall(
        (:VecCreateSeqHIPWithArray, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{Vec}),
        a,
        b,
        c,
        d,
        e,
    )
end

@for_petsc function VecCreateSeqCUDAWithArrays(
    ::$UnionPetscLib,
    a,
    b,
    c,
    d,
    e,
    f,
)
    @chk ccall(
        (:VecCreateSeqCUDAWithArrays, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscScalar},
            Ptr{$PetscScalar},
            Ptr{Vec},
        ),
        a,
        b,
        c,
        d,
        e,
        f,
    )
end

@for_petsc function VecCreateSeqHIPWithArrays(
    ::$UnionPetscLib,
    a,
    b,
    c,
    d,
    e,
    f,
)
    @chk ccall(
        (:VecCreateSeqHIPWithArrays, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscScalar},
            Ptr{$PetscScalar},
            Ptr{Vec},
        ),
        a,
        b,
        c,
        d,
        e,
        f,
    )
end

@for_petsc function VecCreateMPICUDA(::$UnionPetscLib, a, b, c, d)
    @chk ccall(
        (:VecCreateMPICUDA, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, $PetscInt, Ptr{Vec}),
        a,
        b,
        c,
        d,
    )
end

@for_petsc function VecCreateMPIHIP(::$UnionPetscLib, a, b, c, d)
    @chk ccall(
        (:VecCreateMPIHIP, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, $PetscInt, Ptr{Vec}),
        a,
        b,
        c,
        d,
    )
end

@for_petsc function VecCreateMPICUDAWithArray(
    ::$UnionPetscLib,
    a,
    b,
    c,
    d,
    e,
    f,
)
    @chk ccall(
        (:VecCreateMPICUDAWithArray, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscScalar},
            Ptr{Vec},
        ),
        a,
        b,
        c,
        d,
        e,
        f,
    )
end

@for_petsc function VecCreateMPIHIPWithArray(::$UnionPetscLib, a, b, c, d, e, f)
    @chk ccall(
        (:VecCreateMPIHIPWithArray, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscScalar},
            Ptr{Vec},
        ),
        a,
        b,
        c,
        d,
        e,
        f,
    )
end

@for_petsc function VecCreateMPICUDAWithArrays(
    ::$UnionPetscLib,
    a,
    b,
    c,
    d,
    e,
    f,
    g,
)
    @chk ccall(
        (:VecCreateMPICUDAWithArrays, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscScalar},
            Ptr{$PetscScalar},
            Ptr{Vec},
        ),
        a,
        b,
        c,
        d,
        e,
        f,
        g,
    )
end

@for_petsc function VecCreateMPIHIPWithArrays(
    ::$UnionPetscLib,
    a,
    b,
    c,
    d,
    e,
    f,
    g,
)
    @chk ccall(
        (:VecCreateMPIHIPWithArrays, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscScalar},
            Ptr{$PetscScalar},
            Ptr{Vec},
        ),
        a,
        b,
        c,
        d,
        e,
        f,
        g,
    )
end

@for_petsc function VecCUDAGetArray(::$UnionPetscLib, a, b)
    @chk ccall(
        (:VecCUDAGetArray, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        a,
        b,
    )
end

@for_petsc function VecHIPGetArray(::$UnionPetscLib, a, b)
    @chk ccall(
        (:VecHIPGetArray, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        a,
        b,
    )
end

@for_petsc function VecCUDARestoreArray(::$UnionPetscLib, a, b)
    @chk ccall(
        (:VecCUDARestoreArray, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        a,
        b,
    )
end

@for_petsc function VecHIPRestoreArray(::$UnionPetscLib, a, b)
    @chk ccall(
        (:VecHIPRestoreArray, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        a,
        b,
    )
end

@for_petsc function VecCUDAGetArrayRead(::$UnionPetscLib, a, b)
    @chk ccall(
        (:VecCUDAGetArrayRead, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        a,
        b,
    )
end

@for_petsc function VecHIPGetArrayRead(::$UnionPetscLib, a, b)
    @chk ccall(
        (:VecHIPGetArrayRead, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        a,
        b,
    )
end

@for_petsc function VecCUDARestoreArrayRead(::$UnionPetscLib, a, b)
    @chk ccall(
        (:VecCUDARestoreArrayRead, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        a,
        b,
    )
end

@for_petsc function VecHIPRestoreArrayRead(::$UnionPetscLib, a, b)
    @chk ccall(
        (:VecHIPRestoreArrayRead, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        a,
        b,
    )
end

@for_petsc function VecCUDAGetArrayWrite(::$UnionPetscLib, a, b)
    @chk ccall(
        (:VecCUDAGetArrayWrite, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        a,
        b,
    )
end

@for_petsc function VecHIPGetArrayWrite(::$UnionPetscLib, a, b)
    @chk ccall(
        (:VecHIPGetArrayWrite, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        a,
        b,
    )
end

@for_petsc function VecCUDARestoreArrayWrite(::$UnionPetscLib, a, b)
    @chk ccall(
        (:VecCUDARestoreArrayWrite, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        a,
        b,
    )
end

@for_petsc function VecHIPRestoreArrayWrite(::$UnionPetscLib, a, b)
    @chk ccall(
        (:VecHIPRestoreArrayWrite, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Ptr{$PetscScalar}}),
        a,
        b,
    )
end

@for_petsc function VecCUDAPlaceArray(::$UnionPetscLib, a, b)
    @chk ccall(
        (:VecCUDAPlaceArray, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscScalar}),
        a,
        b,
    )
end

@for_petsc function VecHIPPlaceArray(::$UnionPetscLib, a, b)
    @chk ccall(
        (:VecHIPPlaceArray, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscScalar}),
        a,
        b,
    )
end

@for_petsc function VecCUDAReplaceArray(::$UnionPetscLib, a, b)
    @chk ccall(
        (:VecCUDAReplaceArray, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscScalar}),
        a,
        b,
    )
end

@for_petsc function VecHIPReplaceArray(::$UnionPetscLib, a, b)
    @chk ccall(
        (:VecHIPReplaceArray, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{$PetscScalar}),
        a,
        b,
    )
end

@for_petsc function VecCUDAResetArray(::$UnionPetscLib, a)
    @chk ccall((:VecCUDAResetArray, $petsc_library), PetscErrorCode, (Vec,), a)
end

@for_petsc function VecHIPResetArray(::$UnionPetscLib, a)
    @chk ccall((:VecHIPResetArray, $petsc_library), PetscErrorCode, (Vec,), a)
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

@enum PetscSFConcatenateRootMode::UInt32 begin
    PETSCSF_CONCATENATE_ROOTMODE_LOCAL = 0
    PETSCSF_CONCATENATE_ROOTMODE_SHARED = 1
    PETSCSF_CONCATENATE_ROOTMODE_GLOBAL = 2
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

@for_petsc function PetscSFGetRanksSF(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSFGetRanksSF, $petsc_library),
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

@for_petsc function PetscSFConcatenate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscSFConcatenate, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            Ptr{PetscSF},
            PetscSFConcatenateRootMode,
            Ptr{$PetscInt},
            Ptr{PetscSF},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscSFCreateStridedSF(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscSFCreateStridedSF, $petsc_library),
        PetscErrorCode,
        (PetscSF, $PetscInt, $PetscInt, $PetscInt, Ptr{PetscSF}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function PetscSFMerge(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscSFMerge, $petsc_library),
        PetscErrorCode,
        (PetscSF, PetscSF, Ptr{PetscSF}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSFSetGraphFromCoordinates(
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
        (:PetscSFSetGraphFromCoordinates, $petsc_library),
        PetscErrorCode,
        (
            PetscSF,
            $PetscInt,
            $PetscInt,
            $PetscInt,
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

@for_petsc function PetscSFFetchAndOpWithMemTypeBegin(
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
        (:PetscSFFetchAndOpWithMemTypeBegin, $petsc_library),
        PetscErrorCode,
        (
            PetscSF,
            MPI_Datatype,
            PetscMemType,
            Ptr{Cvoid},
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
        arg8,
        arg9,
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

@for_petsc function PetscSectionGetBlockStarts(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionGetBlockStarts, $petsc_library),
        PetscErrorCode,
        (PetscSection, Ptr{PetscBT}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionSetBlockStarts(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionSetBlockStarts, $petsc_library),
        PetscErrorCode,
        (PetscSection, PetscBT),
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

@for_petsc function PetscSectionGetIncludesConstraints(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscSectionGetIncludesConstraints, $petsc_library),
        PetscErrorCode,
        (PetscSection, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionSetIncludesConstraints(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscSectionSetIncludesConstraints, $petsc_library),
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

@for_petsc function PetscSectionLoad(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionLoad, $petsc_library),
        PetscErrorCode,
        (PetscSection, PetscViewer),
        arg1,
        arg2,
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
    arg6,
)
    @chk ccall(
        (:PetscSectionCreateGlobalSection, $petsc_library),
        PetscErrorCode,
        (
            PetscSection,
            PetscSF,
            PetscBool,
            PetscBool,
            PetscBool,
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

@for_petsc function PetscSectionCreateSubdomainSection(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSectionCreateSubdomainSection, $petsc_library),
        PetscErrorCode,
        (PetscSection, IS, Ptr{PetscSection}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSectionCreateComponentSubsection(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSectionCreateComponentSubsection, $petsc_library),
        PetscErrorCode,
        (PetscSection, $PetscInt, Ptr{$PetscInt}, Ptr{PetscSection}),
        arg1,
        arg2,
        arg3,
        arg4,
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

@for_petsc function PetscSectionResetClosurePermutation(::$UnionPetscLib, arg1)
    @chk ccall(
        (:PetscSectionResetClosurePermutation, $petsc_library),
        PetscErrorCode,
        (PetscSection,),
        arg1,
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

@for_petsc function PetscSectionSymCopy(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscSectionSymCopy, $petsc_library),
        PetscErrorCode,
        (PetscSectionSym, PetscSectionSym),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSectionSymDistribute(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSectionSymDistribute, $petsc_library),
        PetscErrorCode,
        (PetscSectionSym, PetscSF, Ptr{PetscSectionSym}),
        arg1,
        arg2,
        arg3,
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

@for_petsc function MatFactorGetCanUseOrdering(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatFactorGetCanUseOrdering, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function MatFactorGetUseOrdering(::$UnionPetscLib, A, b)
    @chk ccall(
        (:MatFactorGetUseOrdering, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{PetscBool}),
        A,
        b,
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

# typedef PetscErrorCode ( MatSolverFn ) ( Mat , MatFactorType , Mat * )
const MatSolverFn = Cvoid

# typedef MatSolverFn * MatSolverFunction
const MatSolverFunction = Ptr{MatSolverFn}

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
        (MatSolverType, MatType, MatFactorType, Ptr{MatSolverFn}),
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
            Ptr{Ptr{MatSolverFn}},
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
        (MatSolverType, MatType, MatFactorType, Ptr{MatSolverFn}),
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
            Ptr{Ptr{MatSolverFn}},
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

@for_petsc function MatProductGetAlgorithm(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatProductGetAlgorithm, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{MatProductAlgorithm}),
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

@for_petsc function MatProductGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatProductGetType, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{MatProductType}),
        arg1,
        arg2,
    )
end

@for_petsc function MatProductGetMats(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatProductGetMats, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}, Ptr{Mat}, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
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

@for_petsc function MatFinalizePackage(::$UnionPetscLib)
    @chk ccall((:MatFinalizePackage, $petsc_library), PetscErrorCode, ())
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

@for_petsc function MatCreateFromOptions(
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
        (:MatCreateFromOptions, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            Ptr{Cchar},
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
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

@for_petsc function MatSetOptionsPrefixFactor(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSetOptionsPrefixFactor, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function MatAppendOptionsPrefixFactor(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatAppendOptionsPrefixFactor, $petsc_library),
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

@for_petsc function MatCreateMPIAIJPERM(
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
        (:MatCreateMPIAIJPERM, $petsc_library),
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

@for_petsc function MatCreateSeqAIJPERM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatCreateSeqAIJPERM, $petsc_library),
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

@for_petsc function MatSeqSELLGetFillRatio(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSeqSELLGetFillRatio, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function MatSeqSELLGetMaxSliceWidth(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSeqSELLGetMaxSliceWidth, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function MatSeqSELLGetAvgSliceWidth(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSeqSELLGetAvgSliceWidth, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function MatSeqSELLSetSliceHeight(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSeqSELLSetSliceHeight, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function MatSeqSELLGetVarSliceSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSeqSELLGetVarSliceSize, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function MatCreateSeqAIJSELL(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatCreateSeqAIJSELL, $petsc_library),
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

@for_petsc function MatCreateMPIAIJSELL(
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
        (:MatCreateMPIAIJSELL, $petsc_library),
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

@for_petsc function MatMPISELLGetLocalMatCondensed(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatMPISELLGetLocalMatCondensed, $petsc_library),
        PetscErrorCode,
        (Mat, MatReuse, Ptr{IS}, Ptr{IS}, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatMPISELLGetSeqSELL(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatMPISELLGetSeqSELL, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}, Ptr{Mat}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
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

@for_petsc function MatUpdateMPIAIJWithArray(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatUpdateMPIAIJWithArray, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscScalar}),
        arg1,
        arg2,
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
        (Mat, PetscCount, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatSetPreallocationCOOLocal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatSetPreallocationCOOLocal, $petsc_library),
        PetscErrorCode,
        (Mat, PetscCount, Ptr{$PetscInt}, Ptr{$PetscInt}),
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

@for_petsc function MatCreateCentering(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatCreateCentering, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, $PetscInt, $PetscInt, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
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

@for_petsc function MatLRCSetMats(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatLRCSetMats, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, Mat, Vec, Mat),
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

@for_petsc function MatNormalGetMat(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatNormalGetMat, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function MatNormalHermitianGetMat(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatNormalHermitianGetMat, $petsc_library),
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

@for_petsc function MatCreateDiagonal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCreateDiagonal, $petsc_library),
        PetscErrorCode,
        (Vec, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function MatDiagonalGetDiagonal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatDiagonalGetDiagonal, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function MatDiagonalRestoreDiagonal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatDiagonalRestoreDiagonal, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function MatDiagonalGetInverseDiagonal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatDiagonalGetInverseDiagonal, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function MatDiagonalRestoreInverseDiagonal(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatDiagonalRestoreInverseDiagonal, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function MatHYPRESetPreallocation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatHYPRESetPreallocation, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
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

@for_petsc function MatPythonGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatPythonGetType, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{Cchar}}),
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

@for_petsc function MatInvertVariableBlockEnvelope(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatInvertVariableBlockEnvelope, $petsc_library),
        PetscErrorCode,
        (Mat, MatReuse, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatComputeVariableBlockEnvelope(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatComputeVariableBlockEnvelope, $petsc_library),
        PetscErrorCode,
        (Mat,),
        arg1,
    )
end

@for_petsc function MatSetValuesIS(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatSetValuesIS, $petsc_library),
        PetscErrorCode,
        (Mat, IS, IS, Ptr{$PetscScalar}, InsertMode),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
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

struct MatStencil{PetscInt}
    k::PetscInt
    j::PetscInt
    i::PetscInt
    c::PetscInt
end

@for_petsc function MatSetValuesStencil(
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
        (:MatSetValuesStencil, $petsc_library),
        PetscErrorCode,
        (
            Mat,
            $PetscInt,
            Ptr{MatStencil{$PetscInt}},
            $PetscInt,
            Ptr{MatStencil{$PetscInt}},
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

@for_petsc function MatSetValuesBlockedStencil(
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
        (:MatSetValuesBlockedStencil, $petsc_library),
        PetscErrorCode,
        (
            Mat,
            $PetscInt,
            Ptr{MatStencil{$PetscInt}},
            $PetscInt,
            Ptr{MatStencil{$PetscInt}},
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
    MAT_STRUCTURAL_SYMMETRY_ETERNAL = 25
    MAT_SPD_ETERNAL = 26
    MAT_OPTION_MAX = 27
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

@for_petsc function MatSeqAIJGetArrayWrite(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSeqAIJGetArrayWrite, $petsc_library),
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

@for_petsc function MatSeqAIJRestoreArrayWrite(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSeqAIJRestoreArrayWrite, $petsc_library),
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

@for_petsc function MatSeqAIJKron(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatSeqAIJKron, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, MatReuse, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
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

@for_petsc function MatDenseGetArrayAndMemType(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatDenseGetArrayAndMemType, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}, Ptr{PetscMemType}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatDenseRestoreArrayAndMemType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatDenseRestoreArrayAndMemType, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatDenseGetArrayReadAndMemType(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatDenseGetArrayReadAndMemType, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}, Ptr{PetscMemType}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatDenseRestoreArrayReadAndMemType(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatDenseRestoreArrayReadAndMemType, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatDenseGetArrayWriteAndMemType(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatDenseGetArrayWriteAndMemType, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Ptr{$PetscScalar}}, Ptr{PetscMemType}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatDenseRestoreArrayWriteAndMemType(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatDenseRestoreArrayWriteAndMemType, $petsc_library),
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
    arg5,
    arg6,
)
    @chk ccall(
        (:MatDenseGetSubMatrix, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{Mat}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
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

@for_petsc function MatIsStructurallySymmetricKnown(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatIsStructurallySymmetricKnown, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{PetscBool}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatIsSPDKnown(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatIsSPDKnown, $petsc_library),
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

@for_petsc function MatGetRowSumAbs(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetRowSumAbs, $petsc_library),
        PetscErrorCode,
        (Mat, Vec),
        arg1,
        arg2,
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

@for_petsc function MatTransposeSymbolic(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatTransposeSymbolic, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}),
        arg1,
        arg2,
    )
end

@for_petsc function MatTransposeSetPrecursor(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatTransposeSetPrecursor, $petsc_library),
        PetscErrorCode,
        (Mat, Mat),
        arg1,
        arg2,
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

@for_petsc function MatMultHermitianTransposeEqual(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatMultHermitianTransposeEqual, $petsc_library),
        PetscErrorCode,
        (Mat, Mat, $PetscInt, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatMultHermitianTransposeAddEqual(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:MatMultHermitianTransposeAddEqual, $petsc_library),
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

@for_petsc function MatGetColumnSums(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetColumnSums, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscScalar}),
        arg1,
        arg2,
    )
end

@for_petsc function MatGetColumnSumsRealPart(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetColumnSumsRealPart, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function MatGetColumnSumsImaginaryPart(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetColumnSumsImaginaryPart, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function MatGetColumnMeans(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetColumnMeans, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscScalar}),
        arg1,
        arg2,
    )
end

@for_petsc function MatGetColumnMeansRealPart(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetColumnMeansRealPart, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function MatGetColumnMeansImaginaryPart(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetColumnMeansImaginaryPart, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscReal}),
        arg1,
        arg2,
    )
end

@for_petsc function MatGetColumnReductions(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatGetColumnReductions, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscReal}),
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

@for_petsc function MatZeroRowsStencil(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatZeroRowsStencil, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{MatStencil{$PetscInt}}, $PetscScalar, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function MatZeroRowsColumnsStencil(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:MatZeroRowsColumnsStencil, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{MatStencil{$PetscInt}}, $PetscScalar, Vec, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
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

@for_petsc function MatAIJGetLocalMat(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatAIJGetLocalMat, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{Mat}),
        arg1,
        arg2,
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

@for_petsc function MatMPIAIJGetNumberNonzeros(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatMPIAIJGetNumberNonzeros, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{PetscCount}),
        arg1,
        arg2,
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

@for_petsc function MatSetValue(::$UnionPetscLib, mat, i, j, va, mode)
    @chk ccall(
        (:MatSetValue, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt, $PetscScalar, InsertMode),
        mat,
        i,
        j,
        va,
        mode,
    )
end

@for_petsc function MatGetValue(::$UnionPetscLib, mat, row, col, va)
    @chk ccall(
        (:MatGetValue, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt, Ptr{$PetscScalar}),
        mat,
        row,
        col,
        va,
    )
end

@for_petsc function MatSetValueLocal(::$UnionPetscLib, mat, i, j, va, mode)
    @chk ccall(
        (:MatSetValueLocal, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt, $PetscScalar, InsertMode),
        mat,
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

@for_petsc function MatMPIAdjToSeqRankZero(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatMPIAdjToSeqRankZero, $petsc_library),
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

mutable struct _p_MatCoarsen end

const MatCoarsen = Ptr{_p_MatCoarsen}

const MatCoarsenType = Ptr{Cchar}

@for_petsc function MatCoarsenCreate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCoarsenCreate, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{MatCoarsen}),
        arg1,
        arg2,
    )
end

@for_petsc function MatCoarsenSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCoarsenSetType, $petsc_library),
        PetscErrorCode,
        (MatCoarsen, MatCoarsenType),
        arg1,
        arg2,
    )
end

@for_petsc function MatCoarsenSetAdjacency(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCoarsenSetAdjacency, $petsc_library),
        PetscErrorCode,
        (MatCoarsen, Mat),
        arg1,
        arg2,
    )
end

@for_petsc function MatCoarsenSetGreedyOrdering(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCoarsenSetGreedyOrdering, $petsc_library),
        PetscErrorCode,
        (MatCoarsen, IS),
        arg1,
        arg2,
    )
end

@for_petsc function MatCoarsenSetStrictAggs(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCoarsenSetStrictAggs, $petsc_library),
        PetscErrorCode,
        (MatCoarsen, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function MatCoarsenApply(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatCoarsenApply, $petsc_library),
        PetscErrorCode,
        (MatCoarsen,),
        arg1,
    )
end

@for_petsc function MatCoarsenDestroy(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatCoarsenDestroy, $petsc_library),
        PetscErrorCode,
        (Ptr{MatCoarsen},),
        arg1,
    )
end

@for_petsc function MatCoarsenRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCoarsenRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function MatCoarsenView(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCoarsenView, $petsc_library),
        PetscErrorCode,
        (MatCoarsen, PetscViewer),
        arg1,
        arg2,
    )
end

@for_petsc function MatCoarsenSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatCoarsenSetFromOptions, $petsc_library),
        PetscErrorCode,
        (MatCoarsen,),
        arg1,
    )
end

@for_petsc function MatCoarsenGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCoarsenGetType, $petsc_library),
        PetscErrorCode,
        (MatCoarsen, Ptr{MatCoarsenType}),
        arg1,
        arg2,
    )
end

@for_petsc function MatCoarsenViewFromOptions(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatCoarsenViewFromOptions, $petsc_library),
        PetscErrorCode,
        (MatCoarsen, PetscObject, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatCoarsenMISKSetDistance(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCoarsenMISKSetDistance, $petsc_library),
        PetscErrorCode,
        (MatCoarsen, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function MatCoarsenMISKGetDistance(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCoarsenMISKGetDistance, $petsc_library),
        PetscErrorCode,
        (MatCoarsen, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function MatCoarsenSetMaximumIterations(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCoarsenSetMaximumIterations, $petsc_library),
        PetscErrorCode,
        (MatCoarsen, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function MatCoarsenSetThreshold(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatCoarsenSetThreshold, $petsc_library),
        PetscErrorCode,
        (MatCoarsen, $PetscReal),
        arg1,
        arg2,
    )
end

@for_petsc function MatCoarsenSetStrengthIndex(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatCoarsenSetStrengthIndex, $petsc_library),
        PetscErrorCode,
        (MatCoarsen, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
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

@for_petsc function MatFactorGetPreferredOrdering(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatFactorGetPreferredOrdering, $petsc_library),
        PetscErrorCode,
        (Mat, MatFactorType, Ptr{MatOrderingType}),
        arg1,
        arg2,
        arg3,
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

@for_petsc function MatSeqDenseInvert(::$UnionPetscLib, arg1)
    @chk ccall(
        (:MatSeqDenseInvert, $petsc_library),
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

@for_petsc function MatPartitioningSetNumberVertexWeights(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:MatPartitioningSetNumberVertexWeights, $petsc_library),
        PetscErrorCode,
        (MatPartitioning, $PetscInt),
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
    MATOP_SETUP = 29
    MATOP_ILUFACTOR_SYMBOLIC = 30
    MATOP_ICCFACTOR_SYMBOLIC = 31
    MATOP_GET_DIAGONAL_BLOCK = 32
    MATOP_SET_INF = 33
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
    MATOP_MATMAT_MULT_SYMBOLIC = 64
    MATOP_MATMAT_MULT_NUMERIC = 65
    MATOP_SET_LOCAL_TO_GLOBAL_MAP = 66
    MATOP_SET_VALUES_LOCAL = 67
    MATOP_ZERO_ROWS_LOCAL = 68
    MATOP_GET_ROW_MAX_ABS = 69
    MATOP_GET_ROW_MIN_ABS = 70
    MATOP_CONVERT = 71
    MATOP_HAS_OPERATION = 72
    MATOP_SET_VALUES_ADIFOR = 74
    MATOP_FD_COLORING_APPLY = 75
    MATOP_SET_FROM_OPTIONS = 76
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
    MATOP_MAT_MULT_SYMBOLIC = 90
    MATOP_MAT_MULT_NUMERIC = 91
    MATOP_PTAP_SYMBOLIC = 93
    MATOP_PTAP_NUMERIC = 94
    MATOP_MAT_TRANSPOSE_MULT_SYMBO = 96
    MATOP_MAT_TRANSPOSE_MULT_NUMER = 97
    MATOP_BIND_TO_CPU = 98
    MATOP_PRODUCTSETFROMOPTIONS = 99
    MATOP_PRODUCTSYMBOLIC = 100
    MATOP_PRODUCTNUMERIC = 101
    MATOP_CONJUGATE = 102
    MATOP_VIEW_NATIVE = 103
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
    MATOP_INVERT_VBLOCK_DIAGONAL = 127
    MATOP_CREATE_SUB_MATRICES_MPI = 128
    MATOP_SET_VALUES_BATCH = 129
    MATOP_TRANSPOSE_MAT_MULT_SYMBO = 131
    MATOP_TRANSPOSE_MAT_MULT_NUMER = 132
    MATOP_TRANSPOSE_COLORING_CREAT = 133
    MATOP_TRANS_COLORING_APPLY_SPT = 134
    MATOP_TRANS_COLORING_APPLY_DEN = 135
    MATOP_RART_SYMBOLIC = 137
    MATOP_RART_NUMERIC = 138
    MATOP_SET_BLOCK_SIZES = 139
    MATOP_AYPX = 140
    MATOP_RESIDUAL = 141
    MATOP_FDCOLORING_SETUP = 142
    MATOP_FIND_OFFBLOCK_ENTRIES = 143
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

@for_petsc function MatShellSetContextDestroy(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatShellSetContextDestroy, $petsc_library),
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

@for_petsc function MatISSetAllowRepeated(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatISSetAllowRepeated, $petsc_library),
        PetscErrorCode,
        (Mat, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function MatISGetAllowRepeated(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatISGetAllowRepeated, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{PetscBool}),
        arg1,
        arg2,
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

@for_petsc function MatISGetLocalToGlobalMapping(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatISGetLocalToGlobalMapping, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{ISLocalToGlobalMapping}, Ptr{ISLocalToGlobalMapping}),
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

@for_petsc function MatGetNullSpaces(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatGetNullSpaces, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{Mat}, Ptr{Ptr{MatNullSpace}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatRestoreNullSpaces(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatRestoreNullSpaces, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{Mat}, Ptr{Ptr{MatNullSpace}}),
        arg1,
        arg2,
        arg3,
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

@for_petsc function MatMFFDInitializePackage(::$UnionPetscLib)
    @chk ccall((:MatMFFDInitializePackage, $petsc_library), PetscErrorCode, ())
end

@for_petsc function MatMFFDFinalizePackage(::$UnionPetscLib)
    @chk ccall((:MatMFFDFinalizePackage, $petsc_library), PetscErrorCode, ())
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

@for_petsc function MatMumpsSetIcntl(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMumpsSetIcntl, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMumpsGetIcntl(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMumpsGetIcntl, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMumpsSetCntl(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMumpsSetCntl, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscReal),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMumpsGetCntl(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMumpsGetCntl, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMumpsGetInfo(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMumpsGetInfo, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMumpsGetInfog(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMumpsGetInfog, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMumpsGetRinfo(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMumpsGetRinfo, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMumpsGetRinfog(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMumpsGetRinfog, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMumpsGetNullPivots(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:MatMumpsGetNullPivots, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatMumpsGetInverse(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatMumpsGetInverse, $petsc_library),
        PetscErrorCode,
        (Mat, Mat),
        arg1,
        arg2,
    )
end

@for_petsc function MatMumpsGetInverseTranspose(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatMumpsGetInverseTranspose, $petsc_library),
        PetscErrorCode,
        (Mat, Mat),
        arg1,
        arg2,
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

@for_petsc function MatBoundToCPU(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatBoundToCPU, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{PetscBool}),
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

@for_petsc function MatSetBindingPropagates(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSetBindingPropagates, $petsc_library),
        PetscErrorCode,
        (Mat, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function MatGetBindingPropagates(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatGetBindingPropagates, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function MatCreateScaLAPACK(
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
        (:MatCreateScaLAPACK, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
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

@for_petsc function MatScaLAPACKSetBlockSizes(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatScaLAPACKSetBlockSizes, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function MatScaLAPACKGetBlockSizes(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:MatScaLAPACKGetBlockSizes, $petsc_library),
        PetscErrorCode,
        (Mat, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
    )
end

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

@for_petsc function MatFilter(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:MatFilter, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscReal, PetscBool, PetscBool),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function MatChop(::$UnionPetscLib, A, tol)
    @chk ccall(
        (:MatChop, $petsc_library),
        PetscErrorCode,
        (Mat, $PetscReal),
        A,
        tol,
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

@for_petsc function MatSeqAIJGetCSRAndMemType(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:MatSeqAIJGetCSRAndMemType, $petsc_library),
        PetscErrorCode,
        (
            Mat,
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscInt}},
            Ptr{Ptr{$PetscScalar}},
            Ptr{PetscMemType},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function MatCreateGraph(
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
        (:MatCreateGraph, $petsc_library),
        PetscErrorCode,
        (
            Mat,
            PetscBool,
            PetscBool,
            $PetscReal,
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

@for_petsc function MatEliminateZeros(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatEliminateZeros, $petsc_library),
        PetscErrorCode,
        (Mat, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function MatCreateDenseFromVecType(
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
        (:MatCreateDenseFromVecType, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            VecType,
            $PetscInt,
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
        arg8,
        arg9,
    )
end

@for_petsc function MatSetHPL(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:MatSetHPL, $petsc_library),
        PetscErrorCode,
        (Mat, Cint),
        arg1,
        arg2,
    )
end

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

@enum DMBlockingType::UInt32 begin
    DM_BLOCKING_TOPOLOGICAL_POINT = 0
    DM_BLOCKING_FIELD_NODE = 1
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
    DM_POLYTOPE_UNKNOWN_CELL = 15
    DM_POLYTOPE_UNKNOWN_FACE = 16
    DM_NUM_POLYTOPES = 17
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

@enum DMReorderDefaultFlag::Int32 begin
    DM_REORDER_DEFAULT_NOTSET = -1
    DM_REORDER_DEFAULT_FALSE = 0
    DM_REORDER_DEFAULT_TRUE = 1
end

mutable struct _p_DMField end

const DMField = Ptr{_p_DMField}

mutable struct _p_UniversalLabel end

const DMUniversalLabel = Ptr{_p_UniversalLabel}

mutable struct _PETSc_DMCEED end

const DMCeed = Ptr{_PETSc_DMCEED}

mutable struct _n_DMGeneratorFunctionList end

const DMGeneratorFunctionList = Ptr{_n_DMGeneratorFunctionList}

mutable struct _p_PetscFE end

const PetscFE = Ptr{_p_PetscFE}

@enum PetscFEJacobianType::UInt32 begin
    PETSCFE_JACOBIAN = 0
    PETSCFE_JACOBIAN_PRE = 1
    PETSCFE_JACOBIAN_DYN = 2
end

const DMLabelType = Ptr{Cchar}

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

@for_petsc function DMLabelSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMLabelSetType, $petsc_library),
        PetscErrorCode,
        (DMLabel, DMLabelType),
        arg1,
        arg2,
    )
end

@for_petsc function DMLabelGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMLabelGetType, $petsc_library),
        PetscErrorCode,
        (DMLabel, Ptr{DMLabelType}),
        arg1,
        arg2,
    )
end

@for_petsc function DMLabelRegister(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMLabelRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}),
        arg1,
        arg2,
    )
end

@for_petsc function DMLabelRegisterAll(::$UnionPetscLib)
    @chk ccall((:DMLabelRegisterAll, $petsc_library), PetscErrorCode, ())
end

@for_petsc function DMLabelRegisterDestroy(::$UnionPetscLib)
    @chk ccall((:DMLabelRegisterDestroy, $petsc_library), PetscErrorCode, ())
end

@for_petsc function DMLabelSetUp(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMLabelSetUp, $petsc_library),
        PetscErrorCode,
        (DMLabel,),
        arg1,
    )
end

@for_petsc function DMLabelSetFromOptions(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMLabelSetFromOptions, $petsc_library),
        PetscErrorCode,
        (DMLabel,),
        arg1,
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

@for_petsc function DMLabelDuplicate(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMLabelDuplicate, $petsc_library),
        PetscErrorCode,
        (DMLabel, Ptr{DMLabel}),
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

@for_petsc function DMLabelGetNonEmptyStratumValuesIS(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:DMLabelGetNonEmptyStratumValuesIS, $petsc_library),
        PetscErrorCode,
        (DMLabel, Ptr{IS}),
        arg1,
        arg2,
    )
end

@for_petsc function DMLabelGetValueIndex(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMLabelGetValueIndex, $petsc_library),
        PetscErrorCode,
        (DMLabel, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
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

@for_petsc function DMLabelGetStratumPointIndex(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMLabelGetStratumPointIndex, $petsc_library),
        PetscErrorCode,
        (DMLabel, $PetscInt, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMLabelCompare(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    message,
)
    @chk ccall(
        (:DMLabelCompare, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, DMLabel, DMLabel, Ptr{PetscBool}, Ptr{Ptr{Cchar}}),
        arg1,
        arg2,
        arg3,
        arg4,
        message,
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

@for_petsc function DMLabelPropagateBegin(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMLabelPropagateBegin, $petsc_library),
        PetscErrorCode,
        (DMLabel, PetscSF),
        arg1,
        arg2,
    )
end

@for_petsc function DMLabelPropagatePush(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMLabelPropagatePush, $petsc_library),
        PetscErrorCode,
        (DMLabel, PetscSF, Ptr{Cvoid}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMLabelPropagateEnd(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMLabelPropagateEnd, $petsc_library),
        PetscErrorCode,
        (DMLabel, PetscSF),
        arg1,
        arg2,
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

@for_petsc function PetscSectionSymLabelGetStratum(
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
        (:PetscSectionSymLabelGetStratum, $petsc_library),
        PetscErrorCode,
        (
            PetscSectionSym,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
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
        arg7,
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

@enum PetscWeakFormKind::UInt32 begin
    PETSC_WF_OBJECTIVE = 0
    PETSC_WF_F0 = 1
    PETSC_WF_F1 = 2
    PETSC_WF_G0 = 3
    PETSC_WF_G1 = 4
    PETSC_WF_G2 = 5
    PETSC_WF_G3 = 6
    PETSC_WF_GP0 = 7
    PETSC_WF_GP1 = 8
    PETSC_WF_GP2 = 9
    PETSC_WF_GP3 = 10
    PETSC_WF_GT0 = 11
    PETSC_WF_GT1 = 12
    PETSC_WF_GT2 = 13
    PETSC_WF_GT3 = 14
    PETSC_WF_BDF0 = 15
    PETSC_WF_BDF1 = 16
    PETSC_WF_BDG0 = 17
    PETSC_WF_BDG1 = 18
    PETSC_WF_BDG2 = 19
    PETSC_WF_BDG3 = 20
    PETSC_WF_BDGP0 = 21
    PETSC_WF_BDGP1 = 22
    PETSC_WF_BDGP2 = 23
    PETSC_WF_BDGP3 = 24
    PETSC_WF_R = 25
    PETSC_WF_CEED = 26
    PETSC_NUM_WF = 27
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

@enum PetscDTSimplexQuadratureType::Int32 begin
    PETSCDTSIMPLEXQUAD_DEFAULT = -1
    PETSCDTSIMPLEXQUAD_CONIC = 0
    PETSCDTSIMPLEXQUAD_MINSYM = 1
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

@for_petsc function PetscQuadratureGetCellType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscQuadratureGetCellType, $petsc_library),
        PetscErrorCode,
        (PetscQuadrature, Ptr{DMPolytopeType}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscQuadratureSetCellType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscQuadratureSetCellType, $petsc_library),
        PetscErrorCode,
        (PetscQuadrature, DMPolytopeType),
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

@for_petsc function PetscQuadratureEqual(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscQuadratureEqual, $petsc_library),
        PetscErrorCode,
        (PetscQuadrature, PetscQuadrature, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
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

@for_petsc function PetscDTTensorQuadratureCreate(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscDTTensorQuadratureCreate, $petsc_library),
        PetscErrorCode,
        (PetscQuadrature, PetscQuadrature, Ptr{PetscQuadrature}),
        arg1,
        arg2,
        arg3,
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

@for_petsc function PetscQuadratureComputePermutations(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscQuadratureComputePermutations, $petsc_library),
        PetscErrorCode,
        (PetscQuadrature, Ptr{$PetscInt}, Ptr{Ptr{IS}}),
        arg1,
        arg2,
        arg3,
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

@for_petsc function PetscDTPTrimmedSize(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDTPTrimmedSize, $petsc_library),
        PetscErrorCode,
        ($PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDTPTrimmedEvalJet(
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
        (:PetscDTPTrimmedEvalJet, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            $PetscInt,
            Ptr{$PetscReal},
            $PetscInt,
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

@for_petsc function PetscDTSimplexQuadrature(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDTSimplexQuadrature, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            $PetscInt,
            PetscDTSimplexQuadratureType,
            Ptr{PetscQuadrature},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function PetscDTCreateDefaultQuadrature(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDTCreateDefaultQuadrature, $petsc_library),
        PetscErrorCode,
        (DMPolytopeType, $PetscInt, Ptr{PetscQuadrature}, Ptr{PetscQuadrature}),
        arg1,
        arg2,
        arg3,
        arg4,
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
    arg6,
)
    @chk ccall(
        (:PetscDTTanhSinhIntegrate, $petsc_library),
        PetscErrorCode,
        (
            Ptr{Cvoid},
            $PetscReal,
            $PetscReal,
            $PetscInt,
            Ptr{Cvoid},
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

@for_petsc function PetscDTTanhSinhIntegrateMPFR(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscDTTanhSinhIntegrateMPFR, $petsc_library),
        PetscErrorCode,
        (
            Ptr{Cvoid},
            $PetscReal,
            $PetscReal,
            $PetscInt,
            Ptr{Cvoid},
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

# typedef PetscErrorCode ( * PetscProbFunc ) ( const PetscReal [ ] , const PetscReal [ ] , PetscReal [ ] )
const PetscProbFunc = Ptr{Cvoid}

@enum DTProbDensityType::UInt32 begin
    DTPROB_DENSITY_CONSTANT = 0
    DTPROB_DENSITY_GAUSSIAN = 1
    DTPROB_DENSITY_MAXWELL_BOLTZMANN = 2
    DTPROB_NUM_DENSITY = 3
end

@for_petsc function PetscPDFMaxwellBoltzmann1D(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscPDFMaxwellBoltzmann1D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscCDFMaxwellBoltzmann1D(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscCDFMaxwellBoltzmann1D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscPDFMaxwellBoltzmann2D(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscPDFMaxwellBoltzmann2D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscCDFMaxwellBoltzmann2D(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscCDFMaxwellBoltzmann2D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscPDFMaxwellBoltzmann3D(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscPDFMaxwellBoltzmann3D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscCDFMaxwellBoltzmann3D(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscCDFMaxwellBoltzmann3D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscPDFGaussian1D(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscPDFGaussian1D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscCDFGaussian1D(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscCDFGaussian1D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscPDFSampleGaussian1D(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscPDFSampleGaussian1D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscPDFGaussian2D(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscPDFGaussian2D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscPDFSampleGaussian2D(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscPDFSampleGaussian2D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscPDFGaussian3D(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscPDFGaussian3D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscPDFSampleGaussian3D(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscPDFSampleGaussian3D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscPDFConstant1D(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscPDFConstant1D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscCDFConstant1D(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscCDFConstant1D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscPDFSampleConstant1D(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscPDFSampleConstant1D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscPDFConstant2D(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscPDFConstant2D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscCDFConstant2D(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscCDFConstant2D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscPDFSampleConstant2D(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscPDFSampleConstant2D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscPDFConstant3D(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscPDFConstant3D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscCDFConstant3D(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscCDFConstant3D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscPDFSampleConstant3D(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscPDFSampleConstant3D, $petsc_library),
        PetscErrorCode,
        (Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscProbCreateFromOptions(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:PetscProbCreateFromOptions, $petsc_library),
        PetscErrorCode,
        (
            $PetscInt,
            Ptr{Cchar},
            Ptr{Cchar},
            Ptr{PetscProbFunc},
            Ptr{PetscProbFunc},
            Ptr{PetscProbFunc},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function PetscProbComputeKSStatistic(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscProbComputeKSStatistic, $petsc_library),
        PetscErrorCode,
        (Vec, PetscProbFunc, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
    )
end

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

@for_petsc function DMClearNamedGlobalVectors(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMClearNamedGlobalVectors, $petsc_library),
        PetscErrorCode,
        (DM,),
        arg1,
    )
end

@for_petsc function DMClearNamedLocalVectors(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMClearNamedLocalVectors, $petsc_library),
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

@for_petsc function DMSetMatrixPreallocateSkip(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetMatrixPreallocateSkip, $petsc_library),
        PetscErrorCode,
        (DM, PetscBool),
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

@for_petsc function DMSetBlockingType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetBlockingType, $petsc_library),
        PetscErrorCode,
        (DM, DMBlockingType),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetBlockingType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetBlockingType, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DMBlockingType}),
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

@for_petsc function DMCreateMassMatrixLumped(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMCreateMassMatrixLumped, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Vec}),
        arg1,
        arg2,
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

@for_petsc function DMExtrude(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMExtrude, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{DM}),
        arg1,
        arg2,
        arg3,
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

@for_petsc function DMGenerate(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMGenerate, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}, PetscBool, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMGenerateRegister(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMGenerateRegister, $petsc_library),
        PetscErrorCode,
        (Ptr{Cchar}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, $PetscInt),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMGenerateRegisterAll(::$UnionPetscLib)
    @chk ccall((:DMGenerateRegisterAll, $petsc_library), PetscErrorCode, ())
end

@for_petsc function DMGenerateRegisterDestroy(::$UnionPetscLib)
    @chk ccall((:DMGenerateRegisterDestroy, $petsc_library), PetscErrorCode, ())
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

@for_petsc function DMAdaptMetric(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMAdaptMetric, $petsc_library),
        PetscErrorCode,
        (DM, Vec, DMLabel, DMLabel, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
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

@for_petsc function DMGetCellCoordinateDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetCellCoordinateDM, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DM}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetCellCoordinateDM(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetCellCoordinateDM, $petsc_library),
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

@for_petsc function DMGetCellCoordinateSection(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetCellCoordinateSection, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscSection}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetCellCoordinateSection(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMSetCellCoordinateSection, $petsc_library),
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

@for_petsc function DMGetCellCoordinates(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetCellCoordinates, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetCellCoordinates(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetCellCoordinates, $petsc_library),
        PetscErrorCode,
        (DM, Vec),
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

@for_petsc function DMGetCoordinatesLocal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetCoordinatesLocal, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Vec}),
        arg1,
        arg2,
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

@for_petsc function DMGetCellCoordinatesLocalSetUp(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMGetCellCoordinatesLocalSetUp, $petsc_library),
        PetscErrorCode,
        (DM,),
        arg1,
    )
end

@for_petsc function DMGetCellCoordinatesLocal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetCellCoordinatesLocal, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetCellCoordinatesLocalNoncollective(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:DMGetCellCoordinatesLocalNoncollective, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Vec}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetCellCoordinatesLocal(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetCellCoordinatesLocal, $petsc_library),
        PetscErrorCode,
        (DM, Vec),
        arg1,
        arg2,
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

@for_petsc function DMSetCoordinateDisc(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMSetCoordinateDisc, $petsc_library),
        PetscErrorCode,
        (DM, PetscFE, PetscBool),
        arg1,
        arg2,
        arg3,
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

@for_petsc function DMGetPeriodicity(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMGetPeriodicity, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Ptr{$PetscReal}}, Ptr{Ptr{$PetscReal}}, Ptr{Ptr{$PetscReal}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMSetPeriodicity(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMSetPeriodicity, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
        arg1,
        arg2,
        arg3,
        arg4,
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
    arg6,
    arg7,
)
    @chk ccall(
        (:DMCreateSectionSubDM, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{$PetscInt},
            Ptr{IS},
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

@for_petsc function DMPrintCellIndices(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMPrintCellIndices, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{Cchar}, $PetscInt, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

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

@for_petsc function DMPrintCellVectorReal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPrintCellVectorReal, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{Cchar}, $PetscInt, Ptr{$PetscReal}),
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

@for_petsc function DMCreateSectionPermutation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMCreateSectionPermutation, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{IS}, Ptr{PetscBT}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMReorderSectionGetDefault(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMReorderSectionGetDefault, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DMReorderDefaultFlag}),
        arg1,
        arg2,
    )
end

@for_petsc function DMReorderSectionSetDefault(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMReorderSectionSetDefault, $petsc_library),
        PetscErrorCode,
        (DM, DMReorderDefaultFlag),
        arg1,
        arg2,
    )
end

@for_petsc function DMReorderSectionGetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMReorderSectionGetType, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{MatOrderingType}),
        arg1,
        arg2,
    )
end

@for_petsc function DMReorderSectionSetType(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMReorderSectionSetType, $petsc_library),
        PetscErrorCode,
        (DM, MatOrderingType),
        arg1,
        arg2,
    )
end

@for_petsc function DMUseTensorOrder(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMUseTensorOrder, $petsc_library),
        PetscErrorCode,
        (DM, PetscBool),
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

@for_petsc function DMGetNaturalSF(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetNaturalSF, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscSF}),
        arg1,
        arg2,
    )
end

@for_petsc function DMSetNaturalSF(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetNaturalSF, $petsc_library),
        PetscErrorCode,
        (DM, PetscSF),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetDefaultConstraints(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMGetDefaultConstraints, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscSection}, Ptr{Mat}, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMSetDefaultConstraints(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMSetDefaultConstraints, $petsc_library),
        PetscErrorCode,
        (DM, PetscSection, Mat, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
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

@for_petsc function DMGetCellDS(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMGetCellDS, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{PetscDS}, Ptr{PetscDS}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMGetRegionDS(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMGetRegionDS, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel, Ptr{IS}, Ptr{PetscDS}, Ptr{PetscDS}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMSetRegionDS(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMSetRegionDS, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel, IS, PetscDS, PetscDS),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMGetRegionNumDS(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMGetRegionNumDS, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{DMLabel}, Ptr{IS}, Ptr{PetscDS}, Ptr{PetscDS}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
    )
end

@for_petsc function DMSetRegionNumDS(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMSetRegionNumDS, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, DMLabel, IS, PetscDS, PetscDS),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
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

@for_petsc function DMCreateFEDefault(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMCreateFEDefault, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{Cchar}, $PetscInt, Ptr{PetscFE}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
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

@for_petsc function DMGetNumAuxiliaryVec(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMGetNumAuxiliaryVec, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMGetAuxiliaryVec(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMGetAuxiliaryVec, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel, $PetscInt, $PetscInt, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMSetAuxiliaryVec(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMSetAuxiliaryVec, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel, $PetscInt, $PetscInt, Vec),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMGetAuxiliaryLabels(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMGetAuxiliaryLabels, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DMLabel}, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMCopyAuxiliaryVec(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMCopyAuxiliaryVec, $petsc_library),
        PetscErrorCode,
        (DM, DM),
        arg1,
        arg2,
    )
end

@for_petsc function DMClearAuxiliaryVec(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMClearAuxiliaryVec, $petsc_library),
        PetscErrorCode,
        (DM,),
        arg1,
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

@for_petsc function DMCreateLabelAtIndex(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMCreateLabelAtIndex, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{Cchar}),
        arg1,
        arg2,
        arg3,
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

@for_petsc function DMGetFirstLabeledPoint(
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
        (:DMGetFirstLabeledPoint, $petsc_library),
        PetscErrorCode,
        (
            DM,
            DM,
            DMLabel,
            $PetscInt,
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{PetscDS},
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

@enum DMCopyLabelsMode::UInt32 begin
    DM_COPY_LABELS_REPLACE = 0
    DM_COPY_LABELS_KEEP = 1
    DM_COPY_LABELS_FAIL = 2
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

@for_petsc function DMSetLabel(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMSetLabel, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel),
        arg1,
        arg2,
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

@for_petsc function DMCopyLabels(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    emode,
)
    @chk ccall(
        (:DMCopyLabels, $petsc_library),
        PetscErrorCode,
        (DM, DM, PetscCopyMode, PetscBool, DMCopyLabelsMode),
        arg1,
        arg2,
        arg3,
        arg4,
        emode,
    )
end

@for_petsc function DMCompareLabels(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMCompareLabels, $petsc_library),
        PetscErrorCode,
        (DM, DM, Ptr{PetscBool}, Ptr{Ptr{Cchar}}),
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
    arg13,
)
    @chk ccall(
        (:DMAddBoundary, $petsc_library),
        PetscErrorCode,
        (
            DM,
            DMBoundaryConditionType,
            Ptr{Cchar},
            DMLabel,
            $PetscInt,
            Ptr{$PetscInt},
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
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

@for_petsc function DMProjectFieldLabel(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    funcs,
    arg10,
    arg11,
)
    @chk ccall(
        (:DMProjectFieldLabel, $petsc_library),
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
        funcs,
        arg10,
        arg11,
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

@for_petsc function DMPolytopeTypeIsHybrid(::$UnionPetscLib, ct)
    ccall(
        (:DMPolytopeTypeIsHybrid, $petsc_library),
        PetscBool,
        (DMPolytopeType,),
        ct,
    )
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

@for_petsc function DMPolytopeTypeSimpleShape(::$UnionPetscLib, dim, simplex)
    ccall(
        (:DMPolytopeTypeSimpleShape, $petsc_library),
        DMPolytopeType,
        ($PetscInt, PetscBool),
        dim,
        simplex,
    )
end

@for_petsc function DMPolytopeTypeGetNumArrangements(::$UnionPetscLib, ct)
    ccall(
        (:DMPolytopeTypeGetNumArrangements, $petsc_library),
        $PetscInt,
        (DMPolytopeType,),
        ct,
    )
end

@for_petsc function DMPolytopeTypeGetArrangement(::$UnionPetscLib, ct, o)
    ccall(
        (:DMPolytopeTypeGetArrangement, $petsc_library),
        Ptr{$PetscInt},
        (DMPolytopeType, $PetscInt),
        ct,
        o,
    )
end

@for_petsc function DMPolytopeTypeGetVertexArrangement(::$UnionPetscLib, ct, o)
    ccall(
        (:DMPolytopeTypeGetVertexArrangement, $petsc_library),
        Ptr{$PetscInt},
        (DMPolytopeType, $PetscInt),
        ct,
        o,
    )
end

@for_petsc function DMPolytopeTypeComposeOrientation(
    ::$UnionPetscLib,
    ct,
    o1,
    o2,
)
    ccall(
        (:DMPolytopeTypeComposeOrientation, $petsc_library),
        $PetscInt,
        (DMPolytopeType, $PetscInt, $PetscInt),
        ct,
        o1,
        o2,
    )
end

@for_petsc function DMPolytopeTypeComposeOrientationInv(
    ::$UnionPetscLib,
    ct,
    o1,
    o2,
)
    ccall(
        (:DMPolytopeTypeComposeOrientationInv, $petsc_library),
        $PetscInt,
        (DMPolytopeType, $PetscInt, $PetscInt),
        ct,
        o1,
        o2,
    )
end

@for_petsc function DMPolytopeMatchOrientation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPolytopeMatchOrientation, $petsc_library),
        PetscErrorCode,
        (
            DMPolytopeType,
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
    )
end

@for_petsc function DMPolytopeMatchVertexOrientation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPolytopeMatchVertexOrientation, $petsc_library),
        PetscErrorCode,
        (
            DMPolytopeType,
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
    )
end

@for_petsc function DMPolytopeGetOrientation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPolytopeGetOrientation, $petsc_library),
        PetscErrorCode,
        (DMPolytopeType, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPolytopeGetVertexOrientation(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPolytopeGetVertexOrientation, $petsc_library),
        PetscErrorCode,
        (DMPolytopeType, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPolytopeInCellTest(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPolytopeInCellTest, $petsc_library),
        PetscErrorCode,
        (DMPolytopeType, Ptr{$PetscReal}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
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

struct DMDALocalInfo{PetscInt}
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

@for_petsc function PFStringSetFunction(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PFStringSetFunction, $petsc_library),
        PetscErrorCode,
        (PF, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function PFInitializePackage(::$UnionPetscLib)
    @chk ccall((:PFInitializePackage, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PFFinalizePackage(::$UnionPetscLib)
    @chk ccall((:PFFinalizePackage, $petsc_library), PetscErrorCode, ())
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

@for_petsc function AOFinalizePackage(::$UnionPetscLib)
    @chk ccall((:AOFinalizePackage, $petsc_library), PetscErrorCode, ())
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

@for_petsc function AORegisterAll(::$UnionPetscLib)
    @chk ccall((:AORegisterAll, $petsc_library), PetscErrorCode, ())
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

mutable struct _p_PetscSpace end

const PetscSpace = Ptr{_p_PetscSpace}

@for_petsc function PetscFEInitializePackage(::$UnionPetscLib)
    @chk ccall((:PetscFEInitializePackage, $petsc_library), PetscErrorCode, ())
end

@for_petsc function PetscFEFinalizePackage(::$UnionPetscLib)
    @chk ccall((:PetscFEFinalizePackage, $petsc_library), PetscErrorCode, ())
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

@for_petsc function PetscSpacePolynomialSetSymmetric(::$UnionPetscLib, sp, s)
    @chk ccall(
        (:PetscSpacePolynomialSetSymmetric, $petsc_library),
        PetscErrorCode,
        (PetscSpace, PetscBool),
        sp,
        s,
    )
end

@for_petsc function PetscSpacePolynomialGetSymmetric(::$UnionPetscLib, sp, s)
    @chk ccall(
        (:PetscSpacePolynomialGetSymmetric, $petsc_library),
        PetscErrorCode,
        (PetscSpace, Ptr{PetscBool}),
        sp,
        s,
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

@for_petsc function PetscSpacePTrimmedSetFormDegree(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscSpacePTrimmedSetFormDegree, $petsc_library),
        PetscErrorCode,
        (PetscSpace, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscSpacePTrimmedGetFormDegree(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscSpacePTrimmedGetFormDegree, $petsc_library),
        PetscErrorCode,
        (PetscSpace, Ptr{$PetscInt}),
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

@for_petsc function PetscSpaceSumSetInterleave(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSpaceSumSetInterleave, $petsc_library),
        PetscErrorCode,
        (PetscSpace, PetscBool, PetscBool),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSpaceSumGetInterleave(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscSpaceSumGetInterleave, $petsc_library),
        PetscErrorCode,
        (PetscSpace, Ptr{PetscBool}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscSpaceCreateSum(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscSpaceCreateSum, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{PetscSpace}, PetscBool, Ptr{PetscSpace}),
        arg1,
        arg2,
        arg3,
        arg4,
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

@for_petsc function PetscDualSpaceGetInteriorSection(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDualSpaceGetInteriorSection, $petsc_library),
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

@for_petsc function PetscDualSpaceEqual(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:PetscDualSpaceEqual, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, PetscDualSpace, Ptr{PetscBool}),
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

@for_petsc function PetscDualSpaceSumSetNumSubspaces(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDualSpaceSumSetNumSubspaces, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceSumGetNumSubspaces(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDualSpaceSumGetNumSubspaces, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceSumSetSubspace(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscDualSpaceSumSetSubspace, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, $PetscInt, PetscDualSpace),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDualSpaceSumGetSubspace(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscDualSpaceSumGetSubspace, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, $PetscInt, Ptr{PetscDualSpace}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDualSpaceSumSetConcatenate(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDualSpaceSumSetConcatenate, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceSumGetConcatenate(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:PetscDualSpaceSumGetConcatenate, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceSumSetInterleave(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscDualSpaceSumSetInterleave, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, PetscBool, PetscBool),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDualSpaceSumGetInterleave(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:PetscDualSpaceSumGetInterleave, $petsc_library),
        PetscErrorCode,
        (PetscDualSpace, Ptr{PetscBool}, Ptr{PetscBool}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscDualSpaceCreateSum(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:PetscDualSpaceCreateSum, $petsc_library),
        PetscErrorCode,
        ($PetscInt, Ptr{PetscDualSpace}, PetscBool, Ptr{PetscDualSpace}),
        arg1,
        arg2,
        arg3,
        arg4,
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

@for_petsc function PetscFECreateVector(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscFECreateVector, $petsc_library),
        PetscErrorCode,
        (PetscFE, $PetscInt, PetscBool, PetscBool, Ptr{PetscFE}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
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

@for_petsc function PetscFECreateByCell(
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
        (:PetscFECreateByCell, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            DMPolytopeType,
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

@for_petsc function PetscFECreateLagrangeByCell(
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
        (:PetscFECreateLagrangeByCell, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            $PetscInt,
            DMPolytopeType,
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

@for_petsc function PetscFECreateFromSpaces(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:PetscFECreateFromSpaces, $petsc_library),
        PetscErrorCode,
        (
            PetscSpace,
            PetscDualSpace,
            PetscQuadrature,
            PetscQuadrature,
            Ptr{PetscFE},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
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

@for_petsc function DMDAGetBoundaryType(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMDAGetBoundaryType, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DMBoundaryType}, Ptr{DMBoundaryType}, Ptr{DMBoundaryType}),
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

@for_petsc function DMDAMapMatStencilToGlobal(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMDAMapMatStencilToGlobal, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{MatStencil{$PetscInt}}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
        arg4,
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

@for_petsc function DMDAVecGetArrayDOFWrite(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDAVecGetArrayDOFWrite, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDAVecRestoreArrayDOFWrite(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMDAVecRestoreArrayDOFWrite, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMDACreatePatchIS(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMDACreatePatchIS, $petsc_library),
        PetscErrorCode,
        (
            DM,
            Ptr{MatStencil{$PetscInt}},
            Ptr{MatStencil{$PetscInt}},
            Ptr{IS},
            PetscBool,
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
    )
end

@for_petsc function DMDAGetLocalInfo(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMDAGetLocalInfo, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DMDALocalInfo{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function MatRegisterDAAD(::$UnionPetscLib)
    @chk ccall((:MatRegisterDAAD, $petsc_library), PetscErrorCode, ())
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

@for_petsc function DMDAConvertToCell(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMDAConvertToCell, $petsc_library),
        PetscErrorCode,
        (DM, MatStencil{$PetscInt}, Ptr{$PetscInt}),
        arg1,
        arg2,
        arg3,
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

@for_petsc function DMPatchZoom(
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
        (:DMPatchZoom, $petsc_library),
        PetscErrorCode,
        (
            DM,
            MatStencil{$PetscInt},
            MatStencil{$PetscInt},
            MPI_Comm,
            Ptr{DM},
            Ptr{PetscSF},
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

@for_petsc function DMPatchSolve(::$UnionPetscLib, arg1)
    @chk ccall((:DMPatchSolve, $petsc_library), PetscErrorCode, (DM,), arg1)
end

@for_petsc function DMPatchGetPatchSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPatchGetPatchSize, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{MatStencil{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPatchSetPatchSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPatchSetPatchSize, $petsc_library),
        PetscErrorCode,
        (DM, MatStencil{$PetscInt}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPatchGetCommSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPatchGetCommSize, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{MatStencil{$PetscInt}}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPatchSetCommSize(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPatchSetCommSize, $petsc_library),
        PetscErrorCode,
        (DM, MatStencil{$PetscInt}),
        arg1,
        arg2,
    )
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

@for_petsc function DMPatchCreateGrid(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMPatchCreateGrid, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            MatStencil{$PetscInt},
            MatStencil{$PetscInt},
            MatStencil{$PetscInt},
            Ptr{DM},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
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
    arg10,
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
        arg10,
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

@enum DMPlexShape::UInt32 begin
    DM_SHAPE_BOX = 0
    DM_SHAPE_BOX_SURFACE = 1
    DM_SHAPE_BALL = 2
    DM_SHAPE_SPHERE = 3
    DM_SHAPE_CYLINDER = 4
    DM_SHAPE_SCHWARZ_P = 5
    DM_SHAPE_GYROID = 6
    DM_SHAPE_DOUBLET = 7
    DM_SHAPE_ANNULUS = 8
    DM_SHAPE_HYPERCUBIC = 9
    DM_SHAPE_ZBOX = 10
    DM_SHAPE_UNKNOWN = 11
end

@enum DMPlexCoordMap::UInt32 begin
    DM_COORD_MAP_NONE = 0
    DM_COORD_MAP_SHEAR = 1
    DM_COORD_MAP_FLARE = 2
    DM_COORD_MAP_ANNULUS = 3
    DM_COORD_MAP_SHELL = 4
    DM_COORD_MAP_UNKNOWN = 5
end

@enum DMPlexCSRAlgorithm::UInt32 begin
    DM_PLEX_CSR_MAT = 0
    DM_PLEX_CSR_GRAPH = 1
    DM_PLEX_CSR_OVERLAP = 2
end

struct _p_DMPlexPointQueue{PetscInt}
    size::PetscInt
    points::Ptr{PetscInt}
    front::PetscInt
    back::PetscInt
    num::PetscInt
end

const DMPlexPointQueue = Ptr{_p_DMPlexPointQueue}

mutable struct _p_PetscLimiter end

const PetscLimiter = Ptr{_p_PetscLimiter}

mutable struct _p_PetscFV end

const PetscFV = Ptr{_p_PetscFV}

struct PetscFVFaceGeom{PetscReal, PetscScalar}
    normal::NTuple{3, PetscReal}
    centroid::NTuple{3, PetscReal}
    grad::NTuple{2, NTuple{3, PetscScalar}}
end

struct PetscFVCellGeom{PetscReal}
    centroid::NTuple{3, PetscReal}
    volume::PetscReal
end

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

@for_petsc function PetscFVFinalizePackage(::$UnionPetscLib)
    @chk ccall((:PetscFVFinalizePackage, $petsc_library), PetscErrorCode, ())
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

@for_petsc function PetscFVCreateDualSpace(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFVCreateDualSpace, $petsc_library),
        PetscErrorCode,
        (PetscFV, DMPolytopeType),
        arg1,
        arg2,
    )
end

@for_petsc function PetscFVClone(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFVClone, $petsc_library),
        PetscErrorCode,
        (PetscFV, Ptr{PetscFV}),
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

@for_petsc function PetscFVIntegrateRHSFunction(
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
        (:PetscFVIntegrateRHSFunction, $petsc_library),
        PetscErrorCode,
        (
            PetscFV,
            PetscDS,
            $PetscInt,
            $PetscInt,
            Ptr{PetscFVFaceGeom{$PetscReal, $PetscScalar}},
            Ptr{$PetscReal},
            Ptr{$PetscScalar},
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

@for_petsc function PetscFVLeastSquaresSetMaxFaces(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscFVLeastSquaresSetMaxFaces, $petsc_library),
        PetscErrorCode,
        (PetscFV, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscDualSpaceApplyFVM(
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
        (:PetscDualSpaceApplyFVM, $petsc_library),
        PetscErrorCode,
        (
            PetscDualSpace,
            $PetscInt,
            $PetscReal,
            Ptr{PetscFVCellGeom{$PetscReal}},
            $PetscInt,
            Ptr{Cvoid},
            Ptr{Cvoid},
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

@for_petsc function DMFieldCreateDSWithDG(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMFieldCreateDSWithDG, $petsc_library),
        PetscErrorCode,
        (DM, DM, $PetscInt, Vec, Vec, Ptr{DMField}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
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
    arg8,
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
    arg13,
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
            Ptr{Ptr{$PetscInt}},
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
)
    @chk ccall(
        (:DMPlexCreateReferenceCell, $petsc_library),
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

@for_petsc function DMPlexGetOrientedCone(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexGetOrientedCone, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexRestoreOrientedCone(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexRestoreOrientedCone, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
        arg1,
        arg2,
        arg3,
        arg4,
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

@for_petsc function DMPlexOrientPoint(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexOrientPoint, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, $PetscInt),
        arg1,
        arg2,
        arg3,
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

@for_petsc function DMPlexFilter(
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
        (:DMPlexFilter, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel, $PetscInt, PetscBool, PetscBool, Ptr{PetscSF}, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
        arg7,
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

@for_petsc function DMPlexGetCellTypeStratum(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexGetCellTypeStratum, $petsc_library),
        PetscErrorCode,
        (DM, DMPolytopeType, Ptr{$PetscInt}, Ptr{$PetscInt}),
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

@for_petsc function DMPlexGetCompressedClosure(
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
        (:DMPlexGetCompressedClosure, $petsc_library),
        PetscErrorCode,
        (
            DM,
            PetscSection,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{PetscSection},
            Ptr{IS},
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

@for_petsc function DMPlexRestoreCompressedClosure(
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
        (:DMPlexRestoreCompressedClosure, $petsc_library),
        PetscErrorCode,
        (
            DM,
            PetscSection,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{Ptr{$PetscInt}},
            Ptr{PetscSection},
            Ptr{IS},
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

@enum DMPlexTPSType::UInt32 begin
    DMPLEX_TPS_SCHWARZ_P = 0
    DMPLEX_TPS_GYROID = 1
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

@for_petsc function DMPlexCopyCoordinates(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexCopyCoordinates, $petsc_library),
        PetscErrorCode,
        (DM, DM),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexCreateCoordinateSpace(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexCreateCoordinateSpace, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, PetscBool, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
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

@for_petsc function DMPlexCreateBoxSurfaceMesh(
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
        (:DMPlexCreateBoxSurfaceMesh, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscReal},
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
)
    @chk ccall(
        (:DMPlexCreateHexCylinderMesh, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, DMBoundaryType, Ptr{DM}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexCreateTPSMesh(
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
        (:DMPlexCreateTPSMesh, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            DMPlexTPSType,
            Ptr{$PetscInt},
            Ptr{DMBoundaryType},
            PetscBool,
            $PetscInt,
            $PetscInt,
            $PetscReal,
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

@for_petsc function DMPlexCreateHypercubicMesh(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMPlexCreateHypercubicMesh, $petsc_library),
        PetscErrorCode,
        (
            MPI_Comm,
            $PetscInt,
            Ptr{$PetscInt},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{DM},
        ),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
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
    arg8,
    arg9,
)
    @chk ccall(
        (:DMPlexExtrude, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            $PetscReal,
            PetscBool,
            PetscBool,
            PetscBool,
            Ptr{$PetscReal},
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
    )
end

@for_petsc function DMPlexInflateToGeomModel(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMPlexInflateToGeomModel, $petsc_library),
        PetscErrorCode,
        (DM,),
        arg1,
    )
end

@for_petsc function DMPlexSetIsoperiodicFaceSF(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMPlexSetIsoperiodicFaceSF, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{PetscSF}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexGetIsoperiodicFaceSF(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMPlexGetIsoperiodicFaceSF, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{$PetscInt}, Ptr{Ptr{PetscSF}}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexSetIsoperiodicFaceTransform(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMPlexSetIsoperiodicFaceTransform, $petsc_library),
        PetscErrorCode,
        (DM, $PetscInt, Ptr{$PetscScalar}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexCheck(::$UnionPetscLib, arg1)
    @chk ccall((:DMPlexCheck, $petsc_library), PetscErrorCode, (DM,), arg1)
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

@for_petsc function DMPlexCheckPointSF(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexCheckPointSF, $petsc_library),
        PetscErrorCode,
        (DM, PetscSF, PetscBool),
        arg1,
        arg2,
        arg3,
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
    arg5,
)
    @chk ccall(
        (:DMPlexCreateFromFile, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, PetscBool, Ptr{DM}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
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

@for_petsc function DMPlexCreateEGADSLiteFromFile(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMPlexCreateEGADSLiteFromFile, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, Ptr{DM}),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function PetscViewerExodusIIOpen(
    ::$UnionPetscLib,
    comm,
    name,
    type,
    exo,
)
    @chk ccall(
        (:PetscViewerExodusIIOpen, $petsc_library),
        PetscErrorCode,
        (MPI_Comm, Ptr{Cchar}, PetscFileMode, Ptr{PetscViewer}),
        comm,
        name,
        type,
        exo,
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

@for_petsc function PetscViewerExodusIISetOrder(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerExodusIISetOrder, $petsc_library),
        PetscErrorCode,
        (PetscViewer, $PetscInt),
        arg1,
        arg2,
    )
end

@for_petsc function PetscViewerExodusIIGetOrder(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:PetscViewerExodusIIGetOrder, $petsc_library),
        PetscErrorCode,
        (PetscViewer, Ptr{$PetscInt}),
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
    arg4,
)
    @chk ccall(
        (:DMPlexPartitionLabelCreateSF, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel, PetscBool, Ptr{PetscSF}),
        arg1,
        arg2,
        arg3,
        arg4,
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

@for_petsc function DMPlexRemapMigrationSF(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexRemapMigrationSF, $petsc_library),
        PetscErrorCode,
        (PetscSF, PetscSF, Ptr{PetscSF}),
        arg1,
        arg2,
        arg3,
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

@for_petsc function DMPlexSetOverlap(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexSetOverlap, $petsc_library),
        PetscErrorCode,
        (DM, DM, $PetscInt),
        arg1,
        arg2,
        arg3,
    )
end

@for_petsc function DMPlexDistributeGetDefault(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexDistributeGetDefault, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexDistributeSetDefault(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexDistributeSetDefault, $petsc_library),
        PetscErrorCode,
        (DM, PetscBool),
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

@for_petsc function DMPlexDistributionSetName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexDistributionSetName, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Cchar}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexDistributionGetName(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexDistributionGetName, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{Ptr{Cchar}}),
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

@for_petsc function DMPlexGetOrdering1D(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetOrdering1D, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{IS}),
        arg1,
        arg2,
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

@for_petsc function DMPlexReorderGetDefault(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexReorderGetDefault, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{DMReorderDefaultFlag}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexReorderSetDefault(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexReorderSetDefault, $petsc_library),
        PetscErrorCode,
        (DM, DMReorderDefaultFlag),
        arg1,
        arg2,
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

@for_petsc function DMPlexCreatePointSF(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexCreatePointSF, $petsc_library),
        PetscErrorCode,
        (DM, PetscSF, PetscBool, Ptr{PetscSF}),
        arg1,
        arg2,
        arg3,
        arg4,
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

@for_petsc function DMPlexCreateOverlapLabelFromLabels(
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
        (:DMPlexCreateOverlapLabelFromLabels, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            Ptr{DMLabel},
            Ptr{$PetscInt},
            $PetscInt,
            Ptr{DMLabel},
            Ptr{$PetscInt},
            PetscSection,
            IS,
            PetscSection,
            IS,
            Ptr{DMLabel},
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
    arg8,
)
    @chk ccall(
        (:DMPlexCreateHybridMesh, $petsc_library),
        PetscErrorCode,
        (
            DM,
            DMLabel,
            DMLabel,
            $PetscInt,
            Ptr{DMLabel},
            Ptr{DMLabel},
            Ptr{DM},
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
    arg6,
)
    @chk ccall(
        (:DMPlexLabelCohesiveComplete, $petsc_library),
        PetscErrorCode,
        (DM, DMLabel, DMLabel, $PetscInt, PetscBool, DM),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
        arg6,
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

@for_petsc function DMPlexFindVertices(::$UnionPetscLib, arg1, arg2, arg3, arg4)
    @chk ccall(
        (:DMPlexFindVertices, $petsc_library),
        PetscErrorCode,
        (DM, Vec, $PetscReal, Ptr{IS}),
        arg1,
        arg2,
        arg3,
        arg4,
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

@for_petsc function DMPlexInsertBoundaryValuesFVM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexInsertBoundaryValuesFVM, $petsc_library),
        PetscErrorCode,
        (DM, PetscFV, Vec, $PetscReal, Ptr{Vec}),
        arg1,
        arg2,
        arg3,
        arg4,
        arg5,
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

@for_petsc function DMPlexGetCellCoordinates(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMPlexGetCellCoordinates, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            Ptr{PetscBool},
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
    )
end

@for_petsc function DMPlexRestoreCellCoordinates(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
)
    @chk ccall(
        (:DMPlexRestoreCellCoordinates, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            Ptr{PetscBool},
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
    arg11,
    arg12,
)
    @chk ccall(
        (:DMPlexMatSetClosureGeneral, $petsc_library),
        PetscErrorCode,
        (
            DM,
            PetscSection,
            PetscSection,
            PetscBool,
            DM,
            PetscSection,
            PetscSection,
            PetscBool,
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
        arg11,
        arg12,
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

@for_petsc function DMPlexReorderCohesiveSupports(::$UnionPetscLib, arg1)
    @chk ccall(
        (:DMPlexReorderCohesiveSupports, $petsc_library),
        PetscErrorCode,
        (DM,),
        arg1,
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

@for_petsc function DMPlexIsSimplex(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexIsSimplex, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscBool}),
        arg1,
        arg2,
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

@for_petsc function DMPlexGetFaceGeometry(
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
        (:DMPlexGetFaceGeometry, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            $PetscInt,
            Vec,
            Vec,
            Ptr{$PetscInt},
            Ptr{Ptr{PetscFVFaceGeom{$PetscReal, $PetscScalar}}},
            Ptr{Ptr{$PetscReal}},
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

@for_petsc function DMPlexRestoreFaceGeometry(
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
        (:DMPlexRestoreFaceGeometry, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscInt,
            $PetscInt,
            Vec,
            Vec,
            Ptr{$PetscInt},
            Ptr{Ptr{PetscFVFaceGeom{$PetscReal, $PetscScalar}}},
            Ptr{Ptr{$PetscReal}},
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

@for_petsc function DMPlexComputeClementInterpolant(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
)
    @chk ccall(
        (:DMPlexComputeClementInterpolant, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Vec),
        arg1,
        arg2,
        arg3,
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

@for_petsc function DMPlexSetSNESLocalFEM(::$UnionPetscLib, arg1, arg2, arg3)
    @chk ccall(
        (:DMPlexSetSNESLocalFEM, $petsc_library),
        PetscErrorCode,
        (DM, PetscBool, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
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

@for_petsc function DMPlexSNESComputeObjectiveFEM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexSNESComputeObjectiveFEM, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Ptr{$PetscReal}, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
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

@for_petsc function DMPlexSNESComputeResidualCEED(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexSNESComputeResidualCEED, $petsc_library),
        PetscErrorCode,
        (DM, Vec, Vec, Ptr{Cvoid}),
        arg1,
        arg2,
        arg3,
        arg4,
    )
end

@for_petsc function DMPlexSNESComputeResidualDS(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
)
    @chk ccall(
        (:DMPlexSNESComputeResidualDS, $petsc_library),
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
    arg12,
)
    @chk ccall(
        (:DMPlexComputeBdJacobianSingle, $petsc_library),
        PetscErrorCode,
        (
            DM,
            $PetscReal,
            PetscWeakForm,
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
        arg12,
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

@for_petsc function DMPlexTSComputeRHSFunctionFVMCEED(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexTSComputeRHSFunctionFVMCEED, $petsc_library),
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

@for_petsc function DMPlexTSComputeRHSFunctionFEM(
    ::$UnionPetscLib,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
)
    @chk ccall(
        (:DMPlexTSComputeRHSFunctionFEM, $petsc_library),
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

@for_petsc function DMPlexGetUseCeed(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexGetUseCeed, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexSetUseCeed(::$UnionPetscLib, arg1, arg2)
    @chk ccall(
        (:DMPlexSetUseCeed, $petsc_library),
        PetscErrorCode,
        (DM, PetscBool),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexGetUseMatClosurePermutation(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:DMPlexGetUseMatClosurePermutation, $petsc_library),
        PetscErrorCode,
        (DM, Ptr{PetscBool}),
        arg1,
        arg2,
    )
end

@for_petsc function DMPlexSetUseMatClosurePermutation(
    ::$UnionPetscLib,
    arg1,
    arg2,
)
    @chk ccall(
        (:DMPlexSetUseMatClosurePermutation, $petsc_library),
        PetscErrorCode,
        (DM, PetscBool),
        arg1,
        arg2,
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

    arg1,
)
