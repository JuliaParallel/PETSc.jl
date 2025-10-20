#
# START OF PROLOGUE
#

using MPI
const MPI_Comm = MPI.Comm
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
#const VecType = Cstring
const Mat = Ptr{Cvoid}
#const MatType = Cstring
const KSP = Ptr{Cvoid}
#const KSPType = Cstring
const SNES = Ptr{Cvoid}
#const SNESType = Cstring
const DM = Ptr{Cvoid}
const PetscDLHandle = Ptr{Cvoid}


const PETSC_DECIDE = -1
const PETSC_DETERMINE = PETSC_DECIDE

PetscInt = Int64
PetscInt64 = Int64
PetscInt32 = Int32
PetscScalar = Float64
PetscReal = Float64
#PetscBool = Bool


# Stuff that I don't really wanty to define by hand, but seem to not be part of the petsc python interface?
mutable struct _p_PetscSF end
const PetscSF = Ptr{_p_PetscSF}

const PETSCSTACKSIZE = 64

const void = Cvoid
const char = Cchar

mutable struct PetscDraw end
mutable struct DMLabel end
mutable struct TSMonitorLGCtx end
mutable struct PetscCtxDestroyFn end

# stuff I need to define to get PETSc.jl to load with "using". We need to find a real solution
mutable struct KSPConvergedReasonViewFn end
mutable struct PC end
mutable struct KSPMonitorFn end
mutable struct KSPConvergenceTestFn end
mutable struct KSPComputeOperatorsFn end
mutable struct KSPComputeRHSFn end
mutable struct KSPComputeInitialGuessFn end
mutable struct PeCtx end
mutable struct KSPPSolveFn end

const PetscObject = Ptr{Cvoid}
const external = Ptr{Cvoid}
const KSPMonitorRegisterFn = Ptr{Cvoid}
const KSPMonitorRegisterCreateFn = Ptr{Cvoid}
const KSPMonitorRegisterDestroyFn = Ptr{Cvoid}
const KSPGuess = Ptr{Cvoid}
const KSPFlexibleModifyPCFn = Ptr{Cvoid}
const PetscVoidFn = Cvoid
const PetscProbFn = Ptr{Cvoid}

const PetscBT = Ptr{Cchar}

# required in Sys_wrappers
mutable struct _n_PetscLogRegistry end
const PetscLogRegistry = Ptr{_n_PetscLogRegistry}

mutable struct _n_PetscIntStack end
const PetscIntStack = Ptr{_n_PetscIntStack}
const PetscLogClass = Cint

mutable struct _p_PetscLogHandler end
const PetscLogHandler = Ptr{_p_PetscLogHandler}
const PetscLogHandlerType = Ptr{Cchar}
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

mutable struct _n_PetscLayout end
const PetscLayout = Ptr{_n_PetscLayout}

mutable struct _p_PetscQuadrature end
const PetscQuadrature = Ptr{_p_PetscQuadrature}

mutable struct _p_PetscRandom end
const PetscRandom = Ptr{_p_PetscRandom}

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

@enum PetscMemType::UInt32 begin
    PETSC_MEMTYPE_HOST = 0
    PETSC_MEMTYPE_DEVICE = 1
    # PETSC_MEMTYPE_CUDA = 1
    PETSC_MEMTYPE_NVSHMEM = 17
    PETSC_MEMTYPE_HIP = 3
    PETSC_MEMTYPE_SYCL = 5
end


# needed for Vec ---
mutable struct _p_IS end
const IS = Ptr{_p_IS}

mutable struct _n_ISColoring end
const ISColoring = Ptr{_n_ISColoring}

mutable struct _p_ISLocalToGlobalMapping end
const ISLocalToGlobalMapping = Ptr{_p_ISLocalToGlobalMapping}

const ISLocalToGlobalMappingType = Ptr{Cchar}

mutable struct _p_PetscSection end
const PetscSection = Ptr{_p_PetscSection}

mutable struct _p_PetscSectionSym end
const PetscSectionSym = Ptr{_p_PetscSectionSym}
const PetscSectionSymType = Ptr{Cchar}



#
# END OF PROLOGUE
#

# load all generated files
#include("../src/LibPETSc_lib.jl")
include("enums_wrappers.jl")
include("senums_wrappers.jl")
include("typedefs_wrappers.jl")
include("struct_wrappers.jl")

include("Sys_wrappers.jl")
include("KSP_wrappers.jl")
include("Vec_wrappers.jl")
include("Vecs_wrappers.jl")

