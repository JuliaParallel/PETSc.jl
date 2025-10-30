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

# Custom display for REPL
function Base.show(io::IO, v::AbstractPetscSNES{PetscLib}) where {PetscLib}
    if v.ptr == C_NULL
        print(io, "PETSc SNES (null pointer)")
        return
    else
        print(io, "PETSc SNES object")
    end
    return nothing
end
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

# needed for Mat ---
#=
mutable struct _p_MatNullSpace end
const MatNullSpace = Ptr{_p_MatNullSpace}

mutable struct _p_MatTransposeColoring end
const MatTransposeColoring = Ptr{_p_MatTransposeColoring}

mutable struct _p_hypre_ParCSRMatrix end
const hypre_ParCSRMatrix = Ptr{_p_hypre_ParCSRMatrix}

mutable struct _p_PetscFunctionList end
const PetscFunctionList = Ptr{_p_PetscFunctionList}
=#
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
include("Vec_wrappers.jl")
include("Vecs_wrappers.jl")
include("Mat_wrappers.jl")
include("PetscOptions_wrappers.jl")
include("KSP_wrappers.jl")
include("PetscObject_wrappers.jl")


