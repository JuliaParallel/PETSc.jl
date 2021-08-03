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
