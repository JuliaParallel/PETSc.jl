#
# START OF PROLOGUE
#

using MPI
using OffsetArrays   # DMStagVecGetArray and friends return OffsetArrays
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
const MPI_Count = MPI.MPI_Count

# We know these will be Cvoid, so just set them to be that
const PetscViewer = Ptr{Cvoid}
const PetscObject = Ptr{Cvoid}


const PETSC_DECIDE = -1
const PETSC_DETERMINE = PETSC_DECIDE
const PETSC_CURRENT = -2     # keep the value already set
const PETSC_UNLIMITED = -3   # no limit, where a function takes a count or a tolerance
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
    own::Bool

    PetscVec{PetscLib}(ptr::CVec = C_NULL, age::Int = 0; own::Bool = true) where {PetscLib} = new{PetscLib}(ptr, age, own)
end

# Convenience constructor from petsclib instance
PetscVec(lib::PetscLib) where {PetscLib} = PetscVec{PetscLib}(C_NULL, lib.age)
PetscVec(ptr::CVec, lib::PetscLib, age::Int = lib.age; own::Bool = true) where {PetscLib} = PetscVec{PetscLib}(ptr, age; own)
Base.convert(::Type{CVec}, v::AbstractPetscVec) = v.ptr
Base.unsafe_convert(::Type{CVec}, v::AbstractPetscVec) = v.ptr
# ------------------------------------------------------

# ----- Custom Julia struct for PETSc Mat -----
const CMat = Ptr{Cvoid}
abstract type AbstractPetscMat{T} end
mutable struct PetscMat{PetscLib} <: AbstractPetscMat{PetscLib}
    ptr::CMat
    age::Int
    own::Bool

    PetscMat{PetscLib}(ptr::CMat = C_NULL, age::Int = 0; own::Bool = true) where {PetscLib} = new{PetscLib}(ptr, age, own)
end

# Convenience constructor from petsclib instance
PetscMat(lib::PetscLib) where {PetscLib} = PetscMat{PetscLib}(C_NULL, lib.age)
PetscMat(ptr::CMat, lib::PetscLib, age::Int = lib.age; own::Bool = true) where {PetscLib} = PetscMat{PetscLib}(ptr, age; own)
Base.convert(::Type{CMat}, v::AbstractPetscMat) = v.ptr
Base.unsafe_convert(::Type{CMat}, v::AbstractPetscMat) = v.ptr
# ------------------------------------------------------

# ----- Custom Julia struct for PETSc KSP -----
# `ptr`, `age` and `own` only; the Julia state of the solver (callbacks, options)
# lives with the PETSc object
const CKSP = Ptr{Cvoid}
abstract type AbstractKSP{T} end
mutable struct KSP{PetscLib} <: AbstractKSP{PetscLib}
    ptr::CKSP
    age::Int
    own::Bool

    KSP{PetscLib}(ptr::CKSP = C_NULL, age::Int = 0; own::Bool = true) where {PetscLib} = new{PetscLib}(ptr, age, own)
end

# Convenience constructor from petsclib instance
KSP(lib::PetscLib) where {PetscLib} = KSP{PetscLib}(C_NULL, lib.age)
KSP(ptr::CKSP, lib::PetscLib, age::Int = lib.age; own::Bool = true) where {PetscLib} = KSP{PetscLib}(ptr, age; own)
Base.convert(::Type{CKSP}, v::AbstractKSP) = v.ptr
Base.unsafe_convert(::Type{CKSP}, v::AbstractKSP) = v.ptr
# ------------------------------------------------------

# ----- Custom Julia struct for PETSc PC -----
const CPC = Ptr{Cvoid}
abstract type AbstractPC{T} end
mutable struct PC{PetscLib} <: AbstractPC{PetscLib}
    ptr::CPC
    age::Int
    own::Bool

    PC{PetscLib}(ptr::CPC = C_NULL, age::Int = 0; own::Bool = true) where {PetscLib} = new{PetscLib}(ptr, age, own)
end

PC(lib::PetscLib) where {PetscLib} = PC{PetscLib}(C_NULL, lib.age)
PC(ptr::CPC, lib::PetscLib, age::Int = lib.age; own::Bool = true) where {PetscLib} = PC{PetscLib}(ptr, age; own)
Base.convert(::Type{CPC}, v::AbstractPC) = v.ptr
Base.unsafe_convert(::Type{CPC}, v::AbstractPC) = v.ptr
# ------------------------------------------------------

# ----- Custom Julia struct for PETSc SNES -----
const CSNES = Ptr{Cvoid}
abstract type AbstractSNES{T} end
mutable struct SNES{PetscLib} <: AbstractSNES{PetscLib}
    ptr::CSNES
    age::Int
    own::Bool

    SNES{PetscLib}(ptr::CSNES = C_NULL, age::Int = 0; own::Bool = true) where {PetscLib} = new{PetscLib}(ptr, age, own)
end

# Convenience constructor from petsclib instance
SNES(lib::PetscLib) where {PetscLib} = SNES{PetscLib}(C_NULL, lib.age)
SNES(ptr::Ptr, lib::PetscLib, age::Int = lib.age; own::Bool = true) where {PetscLib} = SNES{PetscLib}(ptr, age; own)
Base.convert(::Type{CSNES}, v::AbstractSNES) = v.ptr
Base.unsafe_convert(::Type{CSNES}, v::AbstractSNES) = v.ptr
# ------------------------------------------------------

# ----- Custom Julia struct for PETSc DM -----
const CDM = Ptr{Cvoid}
abstract type AbstractPetscDM{T} end
mutable struct PetscDM{PetscLib} <: AbstractPetscDM{PetscLib}
    ptr::CDM
    age::Int
    own::Bool

    PetscDM{PetscLib}(ptr::CDM = C_NULL, age::Int = 0; own::Bool = true) where {PetscLib} = new{PetscLib}(ptr, age, own)
end

# Convenience constructor from petsclib instance
PetscDM(lib::PetscLib) where {PetscLib} = PetscDM{PetscLib}(C_NULL, lib.age)
PetscDM(ptr::CDM, lib::PetscLib, age::Int = lib.age; own::Bool = true) where {PetscLib} = PetscDM{PetscLib}(ptr, age; own)
Base.convert(::Type{CDM}, v::AbstractPetscDM) = v.ptr
Base.unsafe_convert(::Type{CDM}, v::AbstractPetscDM) = v.ptr
# ------------------------------------------------------

# ------------------------------------------------------
# PetscOptions
const COptions = Ptr{Cvoid} 
abstract type AbstractPetscOptions{T} end

mutable struct PetscOptions{PetscLib} <: AbstractPetscOptions{PetscLib}
    ptr::Ptr{Cvoid}
    age::Int
    own::Bool

    PetscOptions{PetscLib}(ptr::Ptr{Cvoid} = C_NULL, age::Int = 0; own::Bool = true) where {PetscLib} = new{PetscLib}(ptr, age, own)
end

# Convenience constructors
PetscOptions(lib::PetscLib) where {PetscLib} = PetscOptions{PetscLib}(C_NULL, lib.age)
PetscOptions(ptr::Ptr{Cvoid}, lib::PetscLib, age::Int = lib.age; own::Bool = true) where {PetscLib} = PetscOptions{PetscLib}(ptr, age; own)

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
    age::Int
    own::Bool

    IS{PetscLib}(ptr::Ptr{Cvoid} = C_NULL, age::Int = 0; own::Bool = true) where {PetscLib} = new{PetscLib}(ptr, age, own)
end

# Convenience constructors
IS(lib::PetscLib) where {PetscLib} = IS{PetscLib}(C_NULL, lib.age)
IS(ptr::Ptr{Cvoid}, lib::PetscLib, age::Int = lib.age; own::Bool = true) where {PetscLib} = IS{PetscLib}(ptr, age; own)

# Conversion methods
Base.convert(::Type{Ptr{Cvoid}}, v::AbstractIS) = v.ptr
Base.unsafe_convert(::Type{Ptr{Cvoid}}, v::AbstractIS) = v.ptr
# Allows a mutable IS to be passed as Ptr{CIS} (= Ptr{Ptr{Cvoid}}) to C
# functions that write the IS handle into the pointed-to slot, e.g.
# DMGetStratumIS. Julia passes pointer_from_objref(v), which is the address of
# v.ptr (the first field), so PETSc writes directly into v.ptr.
Base.unsafe_convert(::Type{Ptr{CIS}}, v::AbstractIS) = Ptr{CIS}(Base.pointer_from_objref(v))
# ------------------------------------------------------

# ------------------------------------------------------
# PF
const CPF = Ptr{Cvoid} 
abstract type AbstractPF{T} end

mutable struct PF{PetscLib} <: AbstractPF{PetscLib}
    ptr::Ptr{Cvoid}
    age::Int
    own::Bool

    PF{PetscLib}(ptr::Ptr{Cvoid} = C_NULL, age::Int = 0; own::Bool = true) where {PetscLib} = new{PetscLib}(ptr, age, own)
end

# Convenience constructors
PF(lib::PetscLib) where {PetscLib} = PF{PetscLib}(C_NULL, lib.age)
PF(ptr::Ptr{Cvoid}, lib::PetscLib, age::Int = lib.age; own::Bool = true) where {PetscLib} = PF{PetscLib}(ptr, age; own)

# Conversion methods
Base.convert(::Type{Ptr{Cvoid}}, v::AbstractPF) = v.ptr
Base.unsafe_convert(::Type{Ptr{Cvoid}}, v::AbstractPF) = v.ptr
# ------------------------------------------------------

# ------------------------------------------------------
# TS
const CTS = Ptr{Cvoid} 
abstract type AbstractTS{T} end

mutable struct TS{PetscLib} <: AbstractTS{PetscLib}
    ptr::CTS
    age::Int
    own::Bool

    TS{PetscLib}(ptr::CTS = C_NULL, age::Int = 0; own::Bool = true) where {PetscLib} = new{PetscLib}(ptr, age, own)
end

# Convenience constructors
TS(lib::PetscLib) where {PetscLib} = TS{PetscLib}(C_NULL, lib.age)
TS(ptr::CTS, lib::PetscLib, age::Int = lib.age; own::Bool = true) where {PetscLib} =
    TS{PetscLib}(ptr, age; own)

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
    age::Int
    own::Bool

    Tao{PetscLib}(ptr::Ptr{Cvoid} = C_NULL, age::Int = 0; own::Bool = true) where {PetscLib} = new{PetscLib}(ptr, age, own)
end

# Convenience constructors
Tao(lib::PetscLib) where {PetscLib} = Tao{PetscLib}(C_NULL, lib.age)
Tao(ptr::Ptr{Cvoid}, lib::PetscLib, age::Int = lib.age; own::Bool = true) where {PetscLib} = Tao{PetscLib}(ptr, age; own)

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
    age::Int
    own::Bool

    AO{PetscLib}(ptr::Ptr{Cvoid} = C_NULL, age::Int = 0; own::Bool = true) where {PetscLib} = new{PetscLib}(ptr, age, own)
end

# Convenience constructors
AO(lib::PetscLib) where {PetscLib} = AO{PetscLib}(C_NULL, lib.age)
AO(ptr::Ptr{Cvoid}, lib::PetscLib, age::Int = lib.age; own::Bool = true) where {PetscLib} = AO{PetscLib}(ptr, age; own)

# Conversion methods
Base.convert(::Type{Ptr{Cvoid}}, v::AbstractAO) = v.ptr
Base.unsafe_convert(::Type{Ptr{Cvoid}}, v::AbstractAO) = v.ptr
# ------------------------------------------------------

# ------------------------------------------------------
# Constructors taking the library *type* (wrappers are called with either the petsclib instance
# or its type, see @for_petsc): look the instance up to get the current age.
for T in (:PetscVec, :PetscMat, :KSP, :PC, :SNES, :PetscDM, :TS, :PetscOptions, :IS, :PF, :Tao, :AO)
    @eval $T(ptr::Ptr{Cvoid}, ::Type{PetscLib}; own::Bool = true) where {PetscLib} = $T(ptr, getlib(PetscLib); own)
end
# ------------------------------------------------------

# ------------------------------------------------------
# Deprecated names (v0.4), kept as aliases for one minor cycle: remove in v0.6
const PetscKSP = KSP
const PetscSNES = SNES
const AbstractPetscKSP = AbstractKSP
const AbstractPetscSNES = AbstractSNES
# ------------------------------------------------------

# Stuff that I don't really want to define by hand, but seem to not be part of the petsc python interface?
mutable struct _p_PetscSF end
const PetscSF = Ptr{_p_PetscSF}

const PETSCSTACKSIZE = 64

const void = Cvoid
const char = Cchar

# PetscDraw and TSMonitorLGCtx are opaque pointer handles declared by the generator (opaque_types.jl);
# the former `mutable struct PetscDraw end` placeholders made `Ref{PetscDraw}()` an undefined reference.
const DMLabel = Ptr{Cvoid}  # C typedef struct _n_DMLabel *DMLabel (pointer type)
mutable struct PetscCtxDestroyFn end
mutable struct PetscErrorCodeFn end

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
include("petscarray.jl")
include("enums_wrappers.jl")
include("senums_wrappers.jl")
include("typedefs_wrappers.jl")
include("opaque_types.jl")
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
include("PetscDA_wrappers.jl")
include("PetscViewer_wrappers.jl")
include("Characteristic_wrappers.jl")
include("PetscSF_wrappers.jl")
include("PetscDualSpace_wrappers.jl")
include("PetscOptions_addons_wrappers.jl")
include("PetscIntStack_wrappers.jl")
include("PetscLog_wrappers.jl")
include("PC_wrappers.jl")
include("extra_wrappers.jl")
