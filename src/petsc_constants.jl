
# export names
# do typealiases have to be exported?  I don't think so

export PETSC_NULL, PETSC_IGNORE, PETSC_DECIDE, PETSC_DETERMINE, PETSCDEFAULT, PETSC_COMM_SELF

export PetscInt, PetscScalar, PetscBool, PetscErrorCode, PetscDataType, PetscReal

export PETSC_INSERT_VALUES, PETSC_ADD_VALUES, PETSC_COPY_VALUESA

export NORM_1, NORM_2, NORM_FROBENIUS, NORM_INFINITY, NORM_MAX

export PETSC_MAT_FLUSH_ASSEMBLY, PETSC_MAT_FINAL_ASSEMBLY

export KSPRICHARDSON, KSPCHEBYSHEV, KSPCG, KSPGROPPCG, KSPIPECG, KSPCGNE, KSPNASH, KSPSTCG, KSPGLTR, KSPFCG, KSPGMRES, KSPFGMRES, KSPLGMRES, KSPDGMRES, KSPPGMRES, KSPTCQMR, SKPBCGS, KSPIBCGS, KSPFBCGS, KSPFBCGSR, KSPBCGSL, KSPCGS, KSPTFQMR, KSPCR, KSPPIPECR, KSPSLQR, SKPPREONLY, KSPQCQ, KSPBICG, KSPMINRES, KSPSYMMLQ, KSPLCD, KSPPYTHON, KSPCR

export KSP_NORM_DEFAULT, KSP_NORM_NONE, KSP_NORM_PRECONDITIONED, KSP_NORM_UNPRECONGDITIONED, KSP_NORM_NATURAL

export KSP_NORM_MAX



export VECSEQ
export VECMPI
export VECSTANDARD
export VECSHARED
export VECSEQCUSP
export VECMPICUSP
export VECCUSP
export VECSEQVIENNACL
export VECMPIVIENNACL
export VECVIENNACL
export VECNEST
export VECSEQPTHREAD
export VECMPIPTHREAD
export VECPTHREAD



# decide which version of Petsc to use
# if no PETSC_DIR or PETSC_ARCH defined, use the one (hopefully) built 
# during package installation
if !haskey(ENV, "PETSC_DIR") && !haskey(ENV, "PETSC_ARCH")
  file_path = joinpath(Pkg.dir("PETSc"), "deps/petsc_evars")
  args = open(readdlm, file_path)
  println("args = ", args)
  ENV["PETSC_DIR"] = args[1]
  ENV["PETSC_ARCH"] = args[2]
end

global const PETSC_DIR = ENV["PETSC_DIR"]
global const PETSC_ARCH = ENV["PETSC_ARCH"]


#PETSC_DIR = readall(`echo $PETSC_DIR`)
#PETSC_ARCH = readall(`echo $PETSC_ARCH`)

#PETSC_DIR = "/home/jared/build/petsc-3.6.0"
#PETSC_ARCH = "arch-linux2-c-debug"

#PETSC_DIR = getenv("PETSC_DIR");
#PETSC_ARCH = getenv("PETSC_ARCH");
#=
if (length(PETSC_DIR) == 0)
  disp("Must have environmental variable PETSC_DIR set")
end
if (length(PETSC_ARCH) == 0)
  disp("Must have environmental variable PETSC_ARCH set")
end
=#

global const libpetsclocation = string(PETSC_DIR, "/", PETSC_ARCH, "/lib/", "libpetsc")
global const petsc = libpetsclocation # for compatability with auto generated wrappers
global const libpetsc = Libdl.dlopen(libpetsclocation)




# definitions of Petsc types
# a way to automatically generate these would be preferable
# define all the data types that are known a priori
typealias Petsc64bitInt Int64
#typealias PetscBLASInt Int32
#typealias PetscScalar Cint
typealias PetscBool Cint
typealias PetscDataType Cint  # C enums are Int32

typealias PetscLogDouble Cdouble
typealias PetscErrorCode Cint
#=
const PETSC_PRECISION_SINGLE = (UInt32)(4)
const PETSC_PRECISION_DOUBLE = (UInt32)(8)
=#


const PETSC_INT = (Int32)(0)
const PETSC_DOUBLE = (Int32)(1)
const PETSC_COMPLEX = (Int32)(2)
const PETSC_LONG = (Int32)(3)
const PETSC_SHORT = (Int32)(4)
const PETSC_FLOAT = (Int32)(5)
const PETSC_CHAR = (Int32)(6)
const PETSC_BIT_LOGICAL = (Int32)(7)
const PETSC_ENUM = (Int32)(8)
const PETSC_BOOL = (Int32)(9)
const PETSC___FLOAT128 = (Int32)(10)
const PETSC_OBJECT = (Int32)(11)
const PETSC_FUNCTION = (Int32)(12)
const PETSC_STRING = (Int32)(12)



# these functions are used for figuring out the sizes of the datatypes
# use the ones in PETSc.jl for writing programs
function PetscDataTypeFromString_(name::AbstractString)
    ptype = Array(Cint, 1)
    found = Array(PetscBool, 1)
    ccall((:PetscDataTypeFromString,petsc),PetscErrorCode,(Cstring,Ptr{PetscDataType},Ptr{PetscBool}), name, ptype, found)

    return ptype[1], convert(Bool, found[1])
end


function PetscDataTypeGetSize_(dtype::PetscDataType)
    datasize = Array(Csize_t, 1)
    ccall((:PetscDataTypeGetSize,petsc),PetscErrorCode,(PetscDataType,Ptr{Csize_t}), dtype, datasize)

    return datasize[1]
end


# define types that depend on the options Petsc was compiled with
(real_type_enum, found_real) = PetscDataTypeFromString_("Real")
(scalar_type_enum, found_scalar) = PetscDataTypeFromString_("Scalar")
int_size = PetscDataTypeGetSize_(PETSC_INT)

# confirm a value was found
@assert(found_real)
@assert(found_scalar)

# figure out what the values mean
if real_type_enum == PETSC_DOUBLE
  precision = "double"
elseif real_type_enum == PETSC_FLOAT
  precision = "single"
else
  println("unknown type of Real")
  println("real_type_enum = ", real_type_enum)
end


if scalar_type_enum == real_type_enum
  scalar_type = "real"
elseif scalar_type_enum == PETSC_COMPLEX
  scalar_type = "complex"
else
  println("unknown type of Scalar")
  println("scalar_type_enum = ", scalar_type_enum)
end

# figure out types

if scalar_type == "real"
  # single or double precision real
  if precision == "single"
    scalar_dtype = Float32
    real_dtype = Float32
  else
    scalar_dtype = Float64
    real_dtype = Float64
  end
else  # scalar_type = complex
  if precision == "single"
    scalar_dtype = Complex64
    real_dtype = Float32
  else
    scalar_dtype = Complex128
    real_dtype = Float64
  end
end

if int_size == 4
   int_dtype = Int32
elseif int_size == 8
   int_dtype = Int64
else
  println("unknown Int size")
  println("int_size = ", int_size)
end

println("PetscScalar type = ", scalar_dtype)
println("PetscReal type = ", real_dtype)
println("PetscInt type = ", int_dtype)

typealias PetscScalar scalar_dtype
typealias PetscReal real_dtype
typealias PetscInt int_dtype


# some useful type unions
PetscInt_arr_or_null = Union(AbstractArray{PetscInt}, Ptr{Void})

#=
const PETSC_PI = pi
const PETSC_MAX_INT = 2147483647
const PETSC_MIN_INT = -PETSC_MAX_INT - 1
const PETSC_MAX_REAL = PETSC_MAX_INT  # made up
const PETSC_INFINITY = PETSC_MAX_REAL / 4.0
const PETSC_NINFINITY = -PETSC_INFINITY
#const PassiveReal = Cint  # what is this?
=#

# some useful constants
#const PassiveScalar = PetscScalar
#const MPIU_MATSCALAR = MPIU_SCALAR
#const MPIU_2INT = Cint
global const PETSC_NULL = C_NULL
global const PETSC_IGNORE = C_NULL
global const PETSC_DECIDE = convert(Int32, -1)
global const PETSC_DETERMINE = PETSC_DECIDE
global const PETSC_DEFAULT = convert(Int32, -2)
global const PETSC_COMM_SELF = MPI.COMM_SELF

typealias MPI_Comm MPI.Comm  # use MPI package communicator type
# typealias to size of MPI communicator value
# this is not defined by the MPI C standard, so it might be
# 32 bits (MPICH), or possibly 64 bits
typealias comm_type typeof(MPI.COMM_WORLD.val) 

global const PETSC_INSERT_VALUES = convert(Int32, 1);
global const PETSC_ADD_VALUES    = convert(Int32, 2);
global const PETSC_COPY_VALUES   = convert(Int32, 0);

typealias NormType Int32
global const NORM_1         = convert(Int32, 0);
global const NORM_2         = convert(Int32, 1);
global const NORM_FROBENIUS = convert(Int32 ,2);
global const NORM_INFINITY  = convert(Int32 ,3);
global const NORM_MAX       = NORM_INFINITY;




global const  PETSC_MAT_FLUSH_ASSEMBLY = convert(Int32, 1)
global const  PETSC_MAT_FINAL_ASSEMBLY = convert(Int32, 0)

# vector formats
global const VECSEQ = "seq"
global const VECMPI = "mpi"
global const VECSTANDARD = "standard"
global const VECSHARED = "shared"
global const VECSEQCUSP = "seqcusp"
global const VECMPICUSP = "mpicusp"
global const VECCUSP = "cusp"
global const VECSEQVIENNACL = "seqviennacl"
global const VECMPIVIENNACL = "mpiviennacl"
global const VECVIENNACL = "viennacl"
global const VECNEST = "nest"
global const VECSEQPTHREAD = "seqpthread"
global const VECMPIPTHREAD = "mpipthread"
global const VECPTHREAD = "pthread"

export PetscMatStructure
export DIFFERENT_NONZERO_PATTERN, SUBSET_NONZERO_PATTERN, SAME_NONZERO_PATTERN
# begin enum MatStructure
typealias PetscMatStructure Int32
global const DIFFERENT_NONZERO_PATTERN = (Int32)(0)
global const SUBSET_NONZERO_PATTERN = (Int32)(1)
global const SAME_NONZERO_PATTERN = (Int32)(2)
# end enum MatStructure

export MATSEQAIJ, MATMPIAIJ
# matrix types
global const MATSAME = "same"
global const MATMAIJ = "maij"
global const MATSEQMAIJ = "seqmaij"
global const MATMPIMAIJ = "mpimaij"
global const MATIS = "is"
global const MATAIJ = "aij"
global const MATSEQAIJ = "seqaij"
global const MATSEQAIJPTHREAD = "seqaijpthread"
global const MATAIJPTHREAD = "aijpthread"
global const MATMPIAIJ = "mpiaij"
global const MATAIJCRL = "aijcrl"
global const MATSEQAIJCRL = "seqaijcrl"
global const MATMPIAIJCRL = "mpiaijcrl"
global const MATAIJCUSP = "aijcusp"
global const MATSEQAIJCUSP = "seqaijcusp"
global const MATMPIAIJCUSP = "mpiaijcusp"
global const MATAIJCUSPARSE = "aijcusparse"
global const MATSEQAIJCUSPARSE = "seqaijcusparse"
global const MATMPIAIJCUSPARSE = "mpiaijcusparse"
global const MATAIJVIENNACL = "aijviennacl"
global const MATSEQAIJVIENNACL = "seqaijviennacl"
global const MATMPIAIJVIENNACL = "mpiaijviennacl"
global const MATAIJPERM = "aijperm"
global const MATSEQAIJPERM = "seqaijperm"
global const MATMPIAIJPERM = "mpiaijperm"
global const MATSHELL = "shell"
global const MATDENSE = "dense"
global const MATSEQDENSE = "seqdense"
global const MATMPIDENSE = "mpidense"
global const MATELEMENTAL = "elemental"
global const MATBAIJ = "baij"
global const MATSEQBAIJ = "seqbaij"
global const MATMPIBAIJ = "mpibaij"
global const MATMPIADJ = "mpiadj"
global const MATSBAIJ = "sbaij"
global const MATSEQSBAIJ = "seqsbaij"
global const MATMPISBAIJ = "mpisbaij"
global const MATSEQBSTRM = "seqbstrm"
global const MATMPIBSTRM = "mpibstrm"
global const MATBSTRM = "bstrm"
global const MATSEQSBSTRM = "seqsbstrm"
global const MATMPISBSTRM = "mpisbstrm"
global const MATSBSTRM = "sbstrm"
global const MATDAAD = "daad"
global const MATMFFD = "mffd"
global const MATNORMAL = "normal"
global const MATLRC = "lrc"
global const MATSCATTER = "scatter"
global const MATBLOCKMAT = "blockmat"
global const MATCOMPOSITE = "composite"
global const MATFFT = "fft"
global const MATFFTW = "fftw"
global const MATSEQCUFFT = "seqcufft"
global const MATTRANSPOSEMAT = "transpose"
global const MATSCHURCOMPLEMENT = "schurcomplement"
global const MATPYTHON = "python"
global const MATHYPRESTRUCT = "hyprestruct"
global const MATHYPRESSTRUCT = "hypresstruct"
global const MATSUBMATRIX = "submatrix"
global const MATLOCALREF = "localref"
global const MATNEST = "nest"

export MAT_LOCAL, MAT_GLOBAL_MAX, MAT_GLOBAL_SUM
typealias MatInfoType Int32
const MAT_LOCAL = (Int32)(1)
const MAT_GLOBAL_MAX = (Int32)(2)
const MAT_GLOBAL_SUM = (Int32)(3)
#

export MAT_INITIAL_MATRIX, MAT_REUSE_MATRIX, MAT_IGNORE_MATRIX
typealias MatReuse Int32
const MAT_INITIAL_MATRIX = (Int32)(0)
const MAT_REUSE_MATRIX = (Int32)(1)
const MAT_IGNORE_MATRIX = (Int32)(2)
#

export KSPType
# types of KSP solvers
typealias KSPType ASCIIString
global const KSPRICHARDSON = "richardson"
global const KSPCHEBYSHEV = "chebyshev"
global const KSPCG       =  "cg"
global const KSPGROPPCG  =  "groppcg"
global const KSPPIPECG   =  "pipecg"
global const   KSPCGNE   =    "cgne"
global const   KSPNASH   =    "nash"
global const   KSPSTCG   =    "stcg"
global const   KSPGLTR   =    "gltr"
global const KSPFCG      =  "fcg"
global const KSPGMRES    =  "gmres"
global const   KSPFGMRES  =   "fgmres"
global const   KSPLGMRES  =   "lgmres"
global const   KSPDGMRES  =   "dgmres"
global const   KSPPGMRES  =   "pgmres"
global const KSPTCQMR     = "tcqmr"
global const KSPBCGS      = "bcgs"
global const   KSPIBCGS   =   "ibcgs"
global const   KSPFBCGS   =   "fbcgs"
global const   KSPFBCGSR  =   "fbcgsr"
global const   KSPBCGSL   =   "bcgsl"
global const KSPCGS       = "cgs"
global const KSPTFQMR     = "tfqmr"
global const KSPCR        = "cr"
global const KSPPIPECR    = "pipecr"
global const KSPLSQR      = "lsqr"
global const KSPPREONLY   = "preonly"
global const KSPQCG       = "qcg"
global const KSPBICG      = "bicg"
global const KSPMINRES    = "minres"
global const KSPSYMMLQ    = "symmlq"
global const KSPLCD       = "lcd"
global const KSPPYTHON    = "python"
global const KSPGCR       = "gcr"

typealias KSPNormType Cint
global const KSP_NORM_DEFAULT = (Int32)(-1)
global const KSP_NORM_NONE = (Int32)(0)
global const KSP_NORM_PRECONDITIONED = (Int32)(1)
global const KSP_NORM_UNPRECONDITIONED = (Int32)(2)
global const KSP_NORM_NATURAL = (Int32)(3)
# end enum KSPNormType

export KSPConvergedReason

typealias KSPConvergedReason Cint
global const KSP_CONVERGED_RTOL_NORMAL = (Int32)(1)
global const KSP_CONVERGED_ATOL_NORMAL = (Int32)(9)
global const KSP_CONVERGED_RTOL = (Int32)(2)
global const KSP_CONVERGED_ATOL = (Int32)(3)
global const KSP_CONVERGED_ITS = (Int32)(4)
global const KSP_CONVERGED_CG_NEG_CURVE = (Int32)(5)
global const KSP_CONVERGED_CG_CONSTRAINED = (Int32)(6)
global const KSP_CONVERGED_STEP_LENGTH = (Int32)(7)
global const KSP_CONVERGED_HAPPY_BREAKDOWN = (Int32)(8)
global const KSP_DIVERGED_NULL = (Int32)(-2)
global const KSP_DIVERGED_ITS = (Int32)(-3)
global const KSP_DIVERGED_DTOL = (Int32)(-4)
global const KSP_DIVERGED_BREAKDOWN = (Int32)(-5)
global const KSP_DIVERGED_BREAKDOWN_BICG = (Int32)(-6)
global const KSP_DIVERGED_NONSYMMETRIC = (Int32)(-7)
global const KSP_DIVERGED_INDEFINITE_PC = (Int32)(-8)
global const KSP_DIVERGED_NANORINF = (Int32)(-9)
global const KSP_DIVERGED_INDEFINITE_MAT = (Int32)(-10)
global const KSP_DIVERGED_PCSETUP_FAILED = (Int32)(-11)
global const KSP_CONVERGED_ITERATING = (Int32)(0)
#
# map values to string for printing
export KSPConvergedReasonDict
global const KSPConvergedReasonDict = Dict{KSPConvergedReason, ASCIIString} (
 KSP_CONVERGED_RTOL_NORMAL => "Converged: RTol normal",
 KSP_CONVERGED_ATOL_NORMAL => "Converged: ATol normal",
 KSP_CONVERGED_RTOL => "Converged: RTol",
 KSP_CONVERGED_ATOL => "converged: ATol",
 KSP_CONVERGED_ITS => "Converged: ITS",
 KSP_CONVERGED_CG_NEG_CURVE => "Converged: CG Negative Curvature",
 KSP_CONVERGED_CG_CONSTRAINED => "Converged: CG Constrained",
 KSP_CONVERGED_STEP_LENGTH => "Converged: step length",
 KSP_CONVERGED_HAPPY_BREAKDOWN => "Converged: happy breakdown",
 KSP_DIVERGED_NULL => "Diverged: Null",
 KSP_DIVERGED_ITS => "Diverged: iteration max",
 KSP_DIVERGED_DTOL => "Diverged divergence tolerance",
 KSP_DIVERGED_BREAKDOWN => "Diverged: Krylov breakdown",
 KSP_DIVERGED_BREAKDOWN_BICG => "Diverged: BICG breakdown",
 KSP_DIVERGED_NONSYMMETRIC => "Diverged: non symmetric",
 KSP_DIVERGED_INDEFINITE_PC => "Diverged: indefinitate PC",
 KSP_DIVERGED_NANORINF => "Diverged nan or inf",
 KSP_DIVERGED_INDEFINITE_MAT => "Diverged: indefinate matrix",
 KSP_DIVERGED_PCSETUP_FAILED => "Diverged: PC setup failed",
 KSP_CONVERGED_ITERATING => "Still running: please be patient",
 )
#


global const KSP_NORM_MAX = KSP_NORM_NATURAL + 1


