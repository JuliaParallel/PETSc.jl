# this summarizes local types generated in all packages
# Its for wrapping new packages; simplifies my life..


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
#const PetscOptions = Ptr{Cvoid}
#mutable struct _n_PetscOptions end
#const PetscOptions = Ptr{_n_PetscOptions}


const PetscViewer = Ptr{Cvoid}
const PetscObject = Ptr{Cvoid}


#const VecType = Cstring
#const Mat = Ptr{Cvoid}
#const MatType = Cstring
#const KSP = Ptr{Cvoid}
#const KSPType = Cstring
#const SNES = Ptr{Cvoid}
#const SNESType = Cstring
#const DM = Ptr{Cvoid}
const PetscDLHandle = Ptr{Cvoid}


const PETSC_DECIDE = -1
const PETSC_DETERMINE = PETSC_DECIDE
const PETSC_COMM_SELF = MPI.COMM_SELF

PetscInt = Int64
PetscInt64 = Int64
PetscInt32 = Int32
PetscScalar = Float64
PetscReal = Float64
#PetscBool = Bool


# ----- Custom Julia struct for PETSc Vec -----
const CVec = Ptr{Cvoid}
abstract type AbstractPetscVec{T} end
mutable struct PetscVec{PetscLib} <: AbstractPetscVec{PetscLib}
    ptr::CVec
    age::Int
    
    # Constructor from pointer and age
    PetscVec{PetscLib}(ptr::CVec, age::Int = 0) where {PetscLib} = new{PetscLib}(ptr, age)
    
    # Constructor for empty Vec (null pointer)
    PetscVec{PetscLib}() where {PetscLib} = new{PetscLib}(Ptr{Cvoid}(C_NULL), 0)
end

# Convenience constructor from petsclib instance
PetscVec(lib::PetscLib) where {PetscLib} = PetscVec{PetscLib}()
PetscVec(ptr::CVec, lib::PetscLib, age::Int = 0) where {PetscLib} = PetscVec{PetscLib}(ptr, age)
Base.convert(::Type{CVec}, v::AbstractPetscVec) = v.ptr
Base.unsafe_convert(::Type{CVec}, v::AbstractPetscVec) = v.ptr
# ------------------------------------------------------

# ----- Custom Julia struct for PETSc Mat -----
const CMat = Ptr{Cvoid}
abstract type AbstractPetscMat{T} end
mutable struct PetscMat{PetscLib} <: AbstractPetscMat{PetscLib}
    ptr::CMat
    age::Int
    
    # Constructor from pointer and age
    PetscMat{PetscLib}(ptr::CMat, age::Int = 0) where {PetscLib} = new{PetscLib}(ptr, age)
    
    # Constructor for empty Mat (null pointer)
    PetscMat{PetscLib}() where {PetscLib} = new{PetscLib}(Ptr{Cvoid}(C_NULL), 0)
end

# Convenience constructor from petsclib instance
PetscMat(lib::PetscLib) where {PetscLib} = PetscMat{PetscLib}()
PetscMat(ptr::CMat, lib::PetscLib, age::Int = 0) where {PetscLib} = PetscMat{PetscLib}(ptr, age)
Base.convert(::Type{CMat}, v::AbstractPetscMat) = v.ptr
Base.unsafe_convert(::Type{CMat}, v::AbstractPetscMat) = v.ptr
# ------------------------------------------------------

# ----- Custom Julia struct for PETSc KSP -----
const CKSP = Ptr{Cvoid}
abstract type AbstractPetscKSP{T} end
mutable struct PetscKSP{PetscLib} <: AbstractPetscKSP{PetscLib}
    ptr::CKSP
    age::Int
    
    # Constructor from pointer and age
    PetscKSP{PetscLib}(ptr::CKSP, age::Int = 0) where {PetscLib} = new{PetscLib}(ptr, age)
    
    # Constructor for empty KSP (null pointer)
    PetscKSP{PetscLib}() where {PetscLib} = new{PetscLib}(Ptr{Cvoid}(C_NULL), 0)
end

# Convenience constructor from petsclib instance
PetscKSP(lib::PetscLib) where {PetscLib} = PetscKSP{PetscLib}()
PetscKSP(ptr::CKSP, lib::PetscLib, age::Int = 0) where {PetscLib} = PetscKSP{PetscLib}(ptr, age)
Base.convert(::Type{CKSP}, v::AbstractPetscKSP) = v.ptr
Base.unsafe_convert(::Type{CKSP}, v::AbstractPetscKSP) = v.ptr


# ------------------------------------------------------

# ----- Custom Julia struct for PETSc SNES -----
const CSNES = Ptr{Cvoid}
abstract type AbstractPetscSNES{T} end
mutable struct PetscSNES{PetscLib} <: AbstractPetscSNES{PetscLib}
    ptr::CSNES
    age::Int
    
    # Constructor from pointer and age
    PetscSNES{PetscLib}(ptr::CSNES, age::Int = 0) where {PetscLib} = new{PetscLib}(ptr, age)
    
    # Constructor for empty SNES (null pointer)
    PetscSNES{PetscLib}() where {PetscLib} = new{PetscLib}(Ptr{Cvoid}(C_NULL), 0)
end

# Convenience constructor from petsclib instance
PetscSNES(lib::PetscLib) where {PetscLib} = PetscSNES{PetscLib}()
PetscSNES(ptr::CSNES, lib::PetscLib, age::Int = 0) where {PetscLib} = PetscSNES{PetscLib}(ptr, age)
Base.convert(::Type{CSNES}, v::AbstractPetscSNES) = v.ptr
Base.unsafe_convert(::Type{CSNES}, v::AbstractPetscSNES) = v.ptr


# ------------------------------------------------------

# ----- Custom Julia struct for PETSc DM -----
const CDM = Ptr{Cvoid}
abstract type AbstractPetscDM{T} end
mutable struct PetscDM{PetscLib} <: AbstractPetscDM{PetscLib}
    ptr::CDM
    age::Int
    
    # Constructor from pointer and age
    PetscDM{PetscLib}(ptr::CDM, age::Int = 0) where {PetscLib} = new{PetscLib}(ptr, age)
    
    # Constructor for empty DM (null pointer)
    PetscDM{PetscLib}() where {PetscLib} = new{PetscLib}(Ptr{Cvoid}(C_NULL), 0)
end

# Convenience constructor from petsclib instance
PetscDM(lib::PetscLib) where {PetscLib} = PetscDM{PetscLib}()
PetscDM(ptr::CDM, lib::PetscLib, age::Int = 0) where {PetscLib} = PetscDM{PetscLib}(ptr, age)
Base.convert(::Type{CDM}, v::AbstractPetscDM) = v.ptr
Base.unsafe_convert(::Type{CDM}, v::AbstractPetscDM) = v.ptr

# Custom display for REPL
function Base.show(io::IO, v::AbstractPetscDM{PetscLib}) where {PetscLib}
    if v.ptr == C_NULL
        print(io, "PETSc DM (null pointer)")
        return
    else
        print(io, "PETSc DM object")
    end
    return nothing
end
# ------------------------------------------------------

# ------------------------------------------------------
# PetscOptions
const COptions = Ptr{Cvoid} 
abstract type AbstractPetscOptions{T} end

mutable struct PetscOptions{PetscLib} <: AbstractPetscOptions{PetscLib}
    ptr::Ptr{Cvoid}
    
    PetscOptions{PetscLib}(ptr::Ptr{Cvoid} = C_NULL) where {PetscLib} = new{PetscLib}(ptr)
end

# Convenience constructors
PetscOptions(lib::PetscLib) where {PetscLib} = PetscOptions{PetscLib}()
PetscOptions(ptr::Ptr{Cvoid}, lib::PetscLib) where {PetscLib} = PetscOptions{PetscLib}(ptr)

# Conversion methods
Base.convert(::Type{Ptr{Cvoid}}, v::AbstractPetscOptions) = v.ptr
Base.unsafe_convert(::Type{Ptr{Cvoid}}, v::AbstractPetscOptions) = v.ptr
# ------------------------------------------------------



# Stuff that I don't really want to define by hand, but seem to not be part of the petsc python interface?
mutable struct _p_PetscSF end
const PetscSF = Ptr{_p_PetscSF}

const PETSCSTACKSIZE = 64

const void = Cvoid
const char = Cchar

mutable struct PetscDraw end
mutable struct DMLabel end
mutable struct TSMonitorLGCtx end
mutable struct PetscCtxDestroyFn end
mutable struct PetscErrorCodeFn end

# stuff I need to define to get PETSc.jl to load with "using". We need to find a real solution
#mutable struct KSPConvergedReasonViewFn end
#mutable struct PC end
mutable struct _n_PC end
const PC = Ptr{_n_PC}

#mutable struct KSPMonitorFn end
#mutable struct KSPConvergenceTestFn end
#mutable struct KSPComputeOperatorsFn end
#mutable struct KSPComputeRHSFn end
#mutable struct KSPComputeInitialGuessFn end
#mutable struct PeCtx end
#mutable struct KSPPSolveFn end
mutable struct MatHtoolKernelFn end


const PetscObject = Ptr{Cvoid}
const external = Ptr{Cvoid}
#const KSPMonitorRegisterFn = Ptr{Cvoid}
#const KSPMonitorRegisterCreateFn = Ptr{Cvoid}
#const KSPMonitorRegisterDestroyFn = Ptr{Cvoid}
#const KSPGuess = Ptr{Cvoid}
#const KSPFlexibleModifyPCFn = Ptr{Cvoid}
const PetscVoidFn = Cvoid
const PetscProbFn = Ptr{Cvoid}

const PetscBT = Ptr{Cchar}

# required in Sys_wrappers
mutable struct _n_PetscLogRegistry end
const PetscLogRegistry = Ptr{_n_PetscLogRegistry}

mutable struct _n_PetscIntStack end
const PetscIntStack = Ptr{_n_PetscIntStack}
#const PetscLogClass = Cint

#mutable struct _p_PetscLogHandler end
#const PetscLogHandler = Ptr{_p_PetscLogHandler}
#const PetscLogHandlerType = Ptr{Cchar}
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

#mutable struct _n_PetscLayout end
#const PetscLayout = Ptr{_n_PetscLayout}
#
#mutable struct _p_PetscQuadrature end
#const PetscQuadrature = Ptr{_p_PetscQuadrature}
#
#mutable struct _p_PetscRandom end
#const PetscRandom = Ptr{_p_PetscRandom}

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


include("../src/autowrapped/enums_wrappers.jl")
include("../src/autowrapped/senums_wrappers.jl")
include("../src/autowrapped/typedefs_wrappers.jl")
include("../src/autowrapped/struct_wrappers.jl")


# autodefined type arguments for class ------
mutable struct SNESFunctionFn end

mutable struct SNESNGSFn end

mutable struct SNESJacobianFn end

mutable struct SNESInitialGuessFn end

mutable struct SNESUpdateFn end

mutable struct _n_SNESLineSearch end
const SNESLineSearch = Ptr{_n_SNESLineSearch}

mutable struct SNESObjectiveFn end

# -------------------------------------------------------

# autodefined type arguments for class ------
mutable struct KSPConvergedReasonViewFn end

mutable struct KSPMonitorFn end

mutable struct KSPConvergenceTestFn end

mutable struct KSPComputeOperatorsFn end

mutable struct KSPComputeRHSFn end

mutable struct KSPComputeInitialGuessFn end

mutable struct _n_PeCtx end
const PeCtx = Ptr{_n_PeCtx}

mutable struct KSPPSolveFn end

mutable struct KSPMonitorRegisterFn end

mutable struct KSPMonitorRegisterCreateFn end

mutable struct KSPMonitorRegisterDestroyFn end

mutable struct _n_KSPGuess end
const KSPGuess = Ptr{_n_KSPGuess}

mutable struct KSPFlexibleModifyPCFn end

# -------------------------------------------------------

# autodefined type arguments for class Mat ------
mutable struct _n_MatNullSpace end
const MatNullSpace = Ptr{_n_MatNullSpace}

mutable struct _n_MatTransposeColoring end
const MatTransposeColoring = Ptr{_n_MatTransposeColoring}

mutable struct _n_hypre_ParCSRMatrix end
const hypre_ParCSRMatrix = Ptr{_n_hypre_ParCSRMatrix}

mutable struct _n_PetscFunctionList end
const PetscFunctionList = Ptr{_n_PetscFunctionList}

# autodefined type arguments for class ------
mutable struct _n_PetscOptionItems end
const PetscOptionItems = Ptr{_n_PetscOptionItems}

# -------------------------------------------------------

# autodefined type arguments for class ------
mutable struct _n_SSL_CTX end
const SSL_CTX = Ptr{_n_SSL_CTX}

mutable struct _n_SSL end
const SSL = Ptr{_n_SSL}

mutable struct _n_hid_t end
const hid_t = Ptr{_n_hid_t}

mutable struct _n_MPI_Request end
const MPI_Request = Ptr{_n_MPI_Request}

mutable struct _n_PetscLogHandler end
const PetscLogHandler = Ptr{_n_PetscLogHandler}

mutable struct _n_PetscLayout end
const PetscLayout = Ptr{_n_PetscLayout}

mutable struct _n_PetscQuadrature end
const PetscQuadrature = Ptr{_n_PetscQuadrature}

mutable struct _n_P4estVertexMaps end
const P4estVertexMaps = Ptr{_n_P4estVertexMaps}

mutable struct _n_pointInterpolationP4est end
const pointInterpolationP4est = Ptr{_n_pointInterpolationP4est}

mutable struct _n_hCsize_t end
const hCsize_t = Ptr{_n_hCsize_t}

mutable struct _n_MPIU_Count end
const MPIU_Count = Ptr{_n_MPIU_Count}

mutable struct _n_PetscBLASInt end
const PetscBLASInt = Ptr{_n_PetscBLASInt}
# -------------------------------------------------------

# autodefined type arguments for class Vec ------
mutable struct _n_ViennaCLVector end
const ViennaCLVector = Ptr{_n_ViennaCLVector}

# -------------------------------------------------------

# autodefined type arguments for class ------
mutable struct _n_Vecs end
const Vecs = Ptr{_n_Vecs}

mutable struct n_PetscRandom end
const PetscRandom = Ptr{n_PetscRandom}
# -------------------------------------------------------

# autodefined type arguments for class ------
mutable struct _n_PCRiCchardsonConvergedReason end
const PCRiCchardsonConvergedReason = Ptr{_n_PCRiCchardsonConvergedReason}

mutable struct PCModifySubMatricesFn end

mutable struct PCMGCoarseSpaceConstructorFn end

mutable struct PCShellPSolveFn end

# -------------------------------------------------------

# autodefined type arguments for class ------
mutable struct PetscPoCintFn end

mutable struct _n_DMInterpolationInfo end
const DMInterpolationInfo = Ptr{_n_DMInterpolationInfo}

mutable struct _n_PetscDS end
const PetscDS = Ptr{_n_PetscDS}

mutable struct _n_PetscFE end
const PetscFE = Ptr{_n_PetscFE}

mutable struct _n_DMField end
const DMField = Ptr{_n_DMField}

mutable struct _n_DMPoCintLocationType end
const DMPoCintLocationType = Ptr{_n_DMPoCintLocationType}

mutable struct PetscSimplePoCintFn end

mutable struct _n_DMSwarmCellDM end
const DMSwarmCellDM = Ptr{_n_DMSwarmCellDM}

mutable struct _n_AO end
const AO = Ptr{_n_AO}

mutable struct _n_PF end
const PF = Ptr{_n_PF}

mutable struct _n_moab_Tag end
const moab_Tag = Ptr{_n_moab_Tag}

mutable struct _n_moab_Range end
const moab_Range = Ptr{_n_moab_Range}

mutable struct _n_moab_EntityHandle end
const moab_EntityHandle = Ptr{_n_moab_EntityHandle}

mutable struct _n_moab_Interface end
const moab_Interface = Ptr{_n_moab_Interface}

mutable struct _n_moab_ParallelComm end
const moab_ParallelComm = Ptr{_n_moab_ParallelComm}

mutable struct _n_moab_EntityType end
const moab_EntityType = Ptr{_n_moab_EntityType}

mutable struct _n_DMPlexTransform end
const DMPlexTransform = Ptr{_n_DMPlexTransform}

mutable struct _n_PetscPartitioner end
const PetscPartitioner = Ptr{_n_PetscPartitioner}

mutable struct _n_PetscFV end
const PetscFV = Ptr{_n_PetscFV}

mutable struct _n_PetscWeakForm end
const PetscWeakForm = Ptr{_n_PetscWeakForm}

mutable struct _n_PetscGeom end
const PetscGeom = Ptr{_n_PetscGeom}

mutable struct _n_PetscHMapI end
const PetscHMapI = Ptr{_n_PetscHMapI}

mutable struct TSIFunctionFn end

mutable struct TSI2FunctionFn end

mutable struct TSI2JacobianFn end

mutable struct TSRHSFunctionFn end

mutable struct TSTransientVariableFn end

mutable struct TSSolutionFn end

mutable struct TSForcingFn end

mutable struct TSIJacobianFn end

mutable struct TSRHSJacobianFn end

mutable struct _n_TS end
const TS = Ptr{_n_TS}

mutable struct DMDATSRHSFunctionLocalFn end

mutable struct DMDATSRHSJacobianLocalFn end

mutable struct DMDATSIFunctionLocalFn end

mutable struct DMDATSIJacobianLocalFn end

# -------------------------------------------------------

# autodefined type arguments for class ------
mutable struct PetscPoCintFn end

mutable struct _n_DMInterpolationInfo end
const DMInterpolationInfo = Ptr{_n_DMInterpolationInfo}

mutable struct _n_PetscDS end
const PetscDS = Ptr{_n_PetscDS}

mutable struct _n_PetscFE end
const PetscFE = Ptr{_n_PetscFE}

mutable struct _n_DMField end
const DMField = Ptr{_n_DMField}

mutable struct _n_DMPoCintLocationType end
const DMPoCintLocationType = Ptr{_n_DMPoCintLocationType}

mutable struct PetscSimplePoCintFn end

mutable struct _n_DMSwarmCellDM end
const DMSwarmCellDM = Ptr{_n_DMSwarmCellDM}

mutable struct _n_AO end
const AO = Ptr{_n_AO}

mutable struct _n_PF end
const PF = Ptr{_n_PF}

mutable struct _n_moab_Tag end
const moab_Tag = Ptr{_n_moab_Tag}

mutable struct _n_moab_Range end
const moab_Range = Ptr{_n_moab_Range}

mutable struct _n_moab_EntityHandle end
const moab_EntityHandle = Ptr{_n_moab_EntityHandle}

mutable struct _n_moab_Interface end
const moab_Interface = Ptr{_n_moab_Interface}

mutable struct _n_moab_ParallelComm end
const moab_ParallelComm = Ptr{_n_moab_ParallelComm}

mutable struct _n_moab_EntityType end
const moab_EntityType = Ptr{_n_moab_EntityType}

mutable struct _n_DMPlexTransform end
const DMPlexTransform = Ptr{_n_DMPlexTransform}

mutable struct _n_PetscPartitioner end
const PetscPartitioner = Ptr{_n_PetscPartitioner}

mutable struct _n_PetscFV end
const PetscFV = Ptr{_n_PetscFV}

mutable struct _n_PetscWeakForm end
const PetscWeakForm = Ptr{_n_PetscWeakForm}

mutable struct _n_PetscGeom end
const PetscGeom = Ptr{_n_PetscGeom}

mutable struct _n_PetscHMapI end
const PetscHMapI = Ptr{_n_PetscHMapI}

mutable struct TSIFunctionFn end

mutable struct TSI2FunctionFn end

mutable struct TSI2JacobianFn end

mutable struct TSRHSFunctionFn end

mutable struct TSTransientVariableFn end

mutable struct TSSolutionFn end

mutable struct TSForcingFn end

mutable struct TSIJacobianFn end

mutable struct TSRHSJacobianFn end

mutable struct _n_TS end
const TS = Ptr{_n_TS}

mutable struct DMDATSRHSFunctionLocalFn end

mutable struct DMDATSRHSJacobianLocalFn end

mutable struct DMDATSIFunctionLocalFn end

mutable struct DMDATSIJacobianLocalFn end

# -------------------------------------------------------

# -------------------------------------------------------
# autodefined type arguments for class ------
mutable struct _n_PetscDrawAxis end
const PetscDrawAxis = Ptr{_n_PetscDrawAxis}

mutable struct _n_PetscDrawLG end
const PetscDrawLG = Ptr{_n_PetscDrawLG}

mutable struct _n_PetscDrawSP end
const PetscDrawSP = Ptr{_n_PetscDrawSP}

mutable struct _n_PetscDrawHG end
const PetscDrawHG = Ptr{_n_PetscDrawHG}

mutable struct _n_PetscDrawBar end
const PetscDrawBar = Ptr{_n_PetscDrawBar}

# -------------------------------------------------------

# -------------------------------------------------------
# autodefined type arguments for class ------
mutable struct _n_PetscDrawAxis end
const PetscDrawAxis = Ptr{_n_PetscDrawAxis}

mutable struct _n_PetscDrawLG end
const PetscDrawLG = Ptr{_n_PetscDrawLG}

mutable struct _n_PetscDrawSP end
const PetscDrawSP = Ptr{_n_PetscDrawSP}

mutable struct _n_PetscDrawHG end
const PetscDrawHG = Ptr{_n_PetscDrawHG}

mutable struct _n_PetscDrawBar end
const PetscDrawBar = Ptr{_n_PetscDrawBar}

# -------------------------------------------------------

# autodefined type arguments for class ------
mutable struct _n_PetscRegressor end
const PetscRegressor = Ptr{_n_PetscRegressor}
# -------------------------------------------------------

# autodefined type arguments for class ------
mutable struct _n_ISColoringValue end
const ISColoringValue = Ptr{_n_ISColoringValue}
# -------------------------------------------------------

# autodefined type arguments for class ------
mutable struct _n_TSTrajectory end
const TSTrajectory = Ptr{_n_TSTrajectory}

mutable struct _n_TSAdapt end
const TSAdapt = Ptr{_n_TSAdapt}

mutable struct _n_TSMonitorVTKCtx end
const TSMonitorVTKCtx = Ptr{_n_TSMonitorVTKCtx}

mutable struct TSRHSJacobianPFn end

mutable struct _n_TSGLLEAdapt end
const TSGLLEAdapt = Ptr{_n_TSGLLEAdapt}

mutable struct TSGLLEAcceptFn end

mutable struct TSAlpha2PredictorFn end

mutable struct TSRHSFunctionFn end
mutable struct TSSolutionFn end
mutable struct TSForcingFn end
mutable struct TSRHSJacobianFn end
mutable struct TSIFunctionFn end
mutable struct TSIJacobianFn end
mutable struct TSI2FunctionFn end
mutable struct TSI2JacobianFn end
mutable struct TSTransientVariableFn end
# -------------------------------------------------------

# autodefined type arguments for class ------
mutable struct _n_TaoLineSearch end
const TaoLineSearch = Ptr{_n_TaoLineSearch}

# -------------------------------------------------------