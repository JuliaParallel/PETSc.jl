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
const PetscViewer = Ptr{Cvoid}
const PetscObject = Ptr{Cvoid}


const PETSC_DECIDE = -1
const PETSC_DETERMINE = PETSC_DECIDE
const PETSC_COMM_SELF = MPI.COMM_SELF

PetscInt = Int64
PetscInt64 = Int64
PetscInt32 = Int32
PetscScalar = Float64
PetscReal = Float64
#PetscBool = Bool

mutable struct _n_ISColoring end
const ISColoring = Ptr{_n_ISColoring}

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
    f!::Function
    updateJ!::Function

    # Constructor from pointer and age
    #PetscSNES{PetscLib}(ptr::CSNES, age::Int = 0) where {PetscLib} = new{PetscLib}(ptr, age)
    
    # Constructor for empty SNES (null pointer)
    PetscSNES{PetscLib}(ptr, age) where {PetscLib} = new{PetscLib}(
                        ptr, 
                        age,
                        x -> error("function not defined"),
                        x -> error("function not defined"),
                        )                  
end

# Convenience constructor from petsclib instance
PetscSNES(lib::PetscLib) where {PetscLib} = PetscSNES{PetscLib}()
PetscSNES(ptr::Ptr, lib::PetscLib, f!::Function, updateJ!::Function, age::Int = 0) where {PetscLib} = PetscSNES{PetscLib}(ptr, age, f!, updateJ!)
PetscSNES(ptr::Ptr, lib::PetscLib, age::Int = 0) where {PetscLib} = PetscSNES{PetscLib}(ptr, age)
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

# ------------------------------------------------------
# IS
const CIS = Ptr{Cvoid} 
abstract type AbstractIS{T} end

mutable struct IS{PetscLib} <: AbstractIS{PetscLib}
    ptr::Ptr{Cvoid}
    
    IS{PetscLib}(ptr::Ptr{Cvoid} = C_NULL) where {PetscLib} = new{PetscLib}(ptr)
end

# Convenience constructors
IS(lib::PetscLib) where {PetscLib} = IS{PetscLib}()
IS(ptr::Ptr{Cvoid}, lib::PetscLib) where {PetscLib} = IS{PetscLib}(ptr)

# Conversion methods
Base.convert(::Type{Ptr{Cvoid}}, v::AbstractIS) = v.ptr
Base.unsafe_convert(::Type{Ptr{Cvoid}}, v::AbstractIS) = v.ptr
# ------------------------------------------------------

# ------------------------------------------------------
# PF
const CPF = Ptr{Cvoid} 
abstract type AbstractPF{T} end

mutable struct PF{PetscLib} <: AbstractPF{PetscLib}
    ptr::Ptr{Cvoid}
    
    PF{PetscLib}(ptr::Ptr{Cvoid} = C_NULL) where {PetscLib} = new{PetscLib}(ptr)
end

# Convenience constructors
PF(lib::PetscLib) where {PetscLib} = PF{PetscLib}()
PF(ptr::Ptr{Cvoid}, lib::PetscLib) where {PetscLib} = PF{PetscLib}(ptr)

# Conversion methods
Base.convert(::Type{Ptr{Cvoid}}, v::AbstractPF) = v.ptr
Base.unsafe_convert(::Type{Ptr{Cvoid}}, v::AbstractPF) = v.ptr
# ------------------------------------------------------

# ------------------------------------------------------
# TS
const CTS = Ptr{Cvoid} 
abstract type AbstractTS{T} end

mutable struct TS{PetscLib} <: AbstractTS{PetscLib}
    ptr::Ptr{Cvoid}
    
    TS{PetscLib}(ptr::Ptr{Cvoid} = C_NULL) where {PetscLib} = new{PetscLib}(ptr)
end

# Convenience constructors
TS(lib::PetscLib) where {PetscLib} = TS{PetscLib}()
TS(ptr::Ptr{Cvoid}, lib::PetscLib) where {PetscLib} = TS{PetscLib}(ptr)

# Conversion methods
Base.convert(::Type{Ptr{Cvoid}}, v::AbstractTS) = v.ptr
Base.unsafe_convert(::Type{Ptr{Cvoid}}, v::AbstractTS) = v.ptr
# ------------------------------------------------------

# ------------------------------------------------------
# Tao
const CTao = Ptr{Cvoid} 
abstract type AbstractTao{T} end

mutable struct Tao{PetscLib} <: AbstractTao{PetscLib}
    ptr::Ptr{Cvoid}
    
    Tao{PetscLib}(ptr::Ptr{Cvoid} = C_NULL) where {PetscLib} = new{PetscLib}(ptr)
end

# Convenience constructors
Tao(lib::PetscLib) where {PetscLib} = Tao{PetscLib}()
Tao(ptr::Ptr{Cvoid}, lib::PetscLib) where {PetscLib} = Tao{PetscLib}(ptr)

# Conversion methods
Base.convert(::Type{Ptr{Cvoid}}, v::AbstractTao) = v.ptr
Base.unsafe_convert(::Type{Ptr{Cvoid}}, v::AbstractTao) = v.ptr
# ------------------------------------------------------

# ------------------------------------------------------
# AO
const CAO = Ptr{Cvoid} 
abstract type AbstractAO{T} end

mutable struct AO{PetscLib} <: AbstractAO{PetscLib}
    ptr::Ptr{Cvoid}
    
    AO{PetscLib}(ptr::Ptr{Cvoid} = C_NULL) where {PetscLib} = new{PetscLib}(ptr)
end

# Convenience constructors
AO(lib::PetscLib) where {PetscLib} = AO{PetscLib}()
AO(ptr::Ptr{Cvoid}, lib::PetscLib) where {PetscLib} = AO{PetscLib}(ptr)

# Conversion methods
Base.convert(::Type{Ptr{Cvoid}}, v::AbstractAO) = v.ptr
Base.unsafe_convert(::Type{Ptr{Cvoid}}, v::AbstractAO) = v.ptr
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
mutable struct _n_PC end
const PC = Ptr{_n_PC}

const PetscObject = Ptr{Cvoid}
const external = Ptr{Cvoid}
const PetscVoidFn = Cvoid
const PetscProbFn = Ptr{Cvoid}
const PetscBT = Ptr{Cchar}

# required in Sys_wrappers
mutable struct _n_PetscLogRegistry end
const PetscLogRegistry = Ptr{_n_PetscLogRegistry}

mutable struct _n_PetscIntStack end
const PetscIntStack = Ptr{_n_PetscIntStack}
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

# needed for Mat ---
#
# END OF PROLOGUE
#

# load all generated files
#include("../src/LibPETSc_lib.jl")
include("petscarray.jl")
include("enums_wrappers.jl")
include("senums_wrappers.jl")
include("typedefs_wrappers.jl")
include("struct_wrappers.jl")
include("Sys_wrappers.jl")
include("Vec_wrappers.jl")
include("Vecs_wrappers.jl")
include("Mat_wrappers.jl")
include("KSP_wrappers.jl")
include("SNES_wrappers.jl")
include("DM_wrappers.jl")
include("PetscOptions_wrappers.jl")
include("PetscObject_wrappers.jl")
include("PetscDraw_wrappers.jl")
include("PetscRegressor_wrappers.jl")
include("PF_wrappers.jl")
include("IS_wrappers.jl")
include("TS_wrappers.jl")
include("AO_wrappers.jl")
include("Tao_wrappers.jl")
include("DMaddons_wrappers.jl")
include("VecTagger_wrappers.jl")
include("PetscDS_wrappers.jl")
include("Mataddons_wrappers.jl")
include("ISaddons_wrappers.jl")
include("SNESLineSearch_wrappers.jl")
include("PetscBag_wrappers.jl")
include("KSPGuess_wrappers.jl")
include("PetscKDTree_wrappers.jl")
include("PetscGridHash_wrappers.jl")
include("PetscSection_wrappers.jl")
include("TSaddons_wrappers.jl")
include("PetscSpace_wrappers.jl")
include("PetscDevice_wrappers.jl")
include("PetscLayout_wrappers.jl")
include("PetscMatlabEngine_wrappers.jl")
include("PetscPartitioner_wrappers.jl")
include("PetscConvEst_wrappers.jl")
include("PetscFE_wrappers.jl")
include("PetscBench_wrappers.jl")
include("PetscToken_wrappers.jl")
include("PetscFunctionList_wrappers.jl")
include("PetscDLLibrary_wrappers.jl")
include("PetscContainer_wrappers.jl")
include("PetscRandom_wrappers.jl")
include("Petsccomm_wrappers.jl")
include("PetscOmpCtrl_wrappers.jl")
include("PetscHeap_wrappers.jl")
include("PetscSegBuffer_wrappers.jl")
include("PetscLimiter_wrappers.jl")
include("PetscFV_wrappers.jl")
include("Tao_addons_wrappers.jl")
include("PetscViewer_wrappers.jl")
include("Characteristic_wrappers.jl")
include("PetscSF_wrappers.jl")
include("PetscDualSpace_wrappers.jl")
include("PetscOptions_addons_wrappers.jl")
include("PetscIntStack_wrappers.jl")
include("PetscLog_wrappers.jl")