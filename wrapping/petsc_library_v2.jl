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
#const PetscOptions = Ptr{Cvoid}
#const PetscViewer = Ptr{Cvoid}
#const PetscObject = Ptr{Cvoid}
#const Vec = Ptr{Cvoid}
#const VecType = Cstring
#const Mat = Ptr{Cvoid}
#const MatType = Cstring
#const KSP = Ptr{Cvoid}
#const KSPType = Cstring
#const SNES = Ptr{Cvoid}
#const SNESType = Cstring
#const DM = Ptr{Cvoid}


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
mutable struct DM end
mutable struct Vec end
mutable struct Mat end
mutable struct KSP end
mutable struct PetscViewer end
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

#
# END OF PROLOGUE
#

# load all generated files
#include("../src/LibPETSc_lib.jl")
include("enums_wrappers.jl")
include("senums_wrappers.jl")
include("typedefs_wrappers.jl")
include("struct_wrappers.jl")
include("KSP_wrappers.jl")


