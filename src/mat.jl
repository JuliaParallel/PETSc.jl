# AbstractMatrix wrapper around Petsc Mat
export Mat, petscview

type Mat{T, MType} <: AbstractSparseMatrix{T,PetscInt}
  p::C.Mat{T}
  assembling::Bool # whether we are in the middle of assemble(vec)
  insertmode::C.InsertMode # current mode for setindex!
  data::Any # keep a reference to anything needed for the Mat
            # -- needed if the Mat is a wrapper around a Julia object,
            #    to prevent the object from being garbage collected.
  function Mat(p::C.Mat{T}, data=nothing)
    A = new(p, false, C.INSERT_VALUES, data)
    chk(C.MatSetType(p, MType))
    finalizer(A, MatDestroy)
    return A
  end
end

comm{T}(a::Mat{T}) = MPI.Comm(C.PetscObjectComm(T, a.p.pobj))

function Mat{T}(::Type{T}, mtype=C.MatType=C.MATSEQ; comm::MPI.Comm=MPI.COMM_WORLD)
  p = Ref{C.Mat{T}}()
  chk(C.MatCreate(comm, p))
  Mat{T, mtype}(p[])
end

function Mat{T}(::Type{T}, m::Integer, n::Integer;
                mlocal::Integer=C.PETSC_DECIDE, nlocal::Integer=C.PETSC_DECIDE,
                nz::Integer=16, nnz::AbstractVector=PetscInt[],
                onz::Integer=0, onnz::AbstractVector=PetscInt[],
                comm::MPI.Comm=MPI.COMM_WORLD,
                mtype::Symbol=C.MATMPIAIJ)

  mat = Mat(T, mtype, comm=comm)
  resize!(mat, m, n, mlocal=mlocal, nlocal=nlocal)
  setpreallocation!(mat, nz=nz, nnz=nnz, onz=onz, onnz=onnz)
  setoption!(mat, C.MAT_ROW_ORIENTED, false)  # julia data is column major

  return mat
end

function MatDestroy{T}(mat::Mat{T})
  PetscFinalized(T) || C.MatDestroy(Ref(mat.p))
end

function petscview{T}(mat::Mat{T})
  viewer = C.PetscViewer{T}(C_NULL)
  chk(C.MatView(mat.p, viewer))
end

function setoption!(m::Mat, option::C.MatOption, val::Bool)
  chk(C.MatSetOption(m.p, option, PetscBool(val)))
  m
end

gettype{T,MT}(a::Mat{T,MT}) = MT

function Base.resize!(a::Mat, m::Integer=C.PETSC_DECIDE, n::Integer=C.PETSC_DECIDE;
  mlocal::Integer=C.PETSC_DECIDE, nlocal::Integer=C.PETSC_DECIDE)
  if m == mlocal == C.PETSC_DECIDE
    throw(ArgumentError("either the global (m) or local (mlocal) #rows must be specified"))
  end
  if n == nlocal == C.PETSC_DECIDE
    throw(ArgumentError("either the global (n) or local (nlocal) #cols must be specified"))
  end
  chk(C.MatSetSizes(a.p, mlocal, nlocal, m, n))
  a
end

function setpreallocation!{T, MType}(a::Mat{T, MType};
  nz::Integer=16, nnz::AbstractVector=PetscInt[],
  onz::Integer=0, onnz::AbstractVector=PetscInt[],
  bs::Integer=16)
  if MType == C.MATSEQAIJ
    pnnz = if isempty(nnz)
      Ptr{PetscInt}(0)
    else
      if length(nnz) != size(a,1)
        throw(ArgumentError("length(nnz) must be # rows"))
      end
      isa(nnz,Vector{PetscInt}) ? nnz : PetscInt[ i for i in nnz ]
    end
    chk(C.MatSeqAIJSetPreallocation(a.p, nz, pnnz))
  elseif MType == C.MATMPIAIJ
    mlocal = sizelocal(a,1)
    pnnz = if isempty(nnz)
      Ptr{PetscInt}(0)
    else
      if length(nnz) != mlocal
        throw(ArgumentError("length(nnz) must be # local rows"))
      end
      isa(nnz,Vector{PetscInt}) ? nnz : PetscInt[ i for i in nnz ]
    end
    ponnz = if isempty(onnz)
      Ptr{PetscInt}(0)
    else
      if length(onnz) != mlocal
        throw(ArgumentError("length(onnz) must be # local rows"))
      end
      isa(onnz,Vector{PetscInt}) ? onnz : PetscInt[ i for i in onnz ]
    end
    chk(C.MatMPIAIJSetPreallocation(a.p, nz, pnnz, onz, ponnz))
  elseif MType == C.MATMPIBAIJ
    mlocal = sizelocal(a,1)
    pnnz = if isempty(nnz)
      Ptr{PetscInt}(0)
    else
      if length(nnz) != mlocal
        throw(ArgumentError("length(nnz) must be # local rows"))
      end
      isa(nnz,Vector{PetscInt}) ? nnz : PetscInt[ i for i in nnz ]
    end
    ponnz = if isempty(onnz)
      Ptr{PetscInt}(0)
    else
      if length(onnz) != mlocal
        throw(ArgumentError("length(onnz) must be # local rows"))
      end
      isa(onnz,Vector{PetscInt}) ? onnz : PetscInt[ i for i in onnz ]
    end
    chk(C.MatMPIBAIJSetPreallocation(a.p, bs, nz, pnnz, onz, ponnz))
  elseif MType == C.MATMPISBAIJ
    mlocal = sizelocal(a,1)
    pnnz = if isempty(nnz)
      Ptr{PetscInt}(0)
    else
      if length(nnz) != mlocal
        throw(ArgumentError("length(nnz) must be # local rows"))
      end
      isa(nnz,Vector{PetscInt}) ? nnz : PetscInt[ i for i in nnz ]
    end
    ponnz = if isempty(onnz)
      Ptr{PetscInt}(0)
    else
      if length(onnz) != mlocal
        throw(ArgumentError("length(onnz) must be # local rows"))
      end
      isa(onnz,Vector{PetscInt}) ? onnz : PetscInt[ i for i in onnz ]
    end
    chk(C.MatMPISBAIJSetPreallocation(a.p, bs, nz, pnnz, onz, ponnz))
  elseif MType == C.MATBLOCKMAT
    pnnz = if isempty(nnz)
      Ptr{PetscInt}(0)
    else
      if length(nnz) != size(a,1)
        throw(ArgumentError("length(nnz) must be # rows"))
      end
      isa(nnz,Vector{PetscInt}) ? nnz : PetscInt[ i for i in nnz ]
    end
    chk(C.MatBlockMatSetPreallocation(a.p, bs, nz, pnnz))
  elseif MType == C.MATSEQBAIJ
    pnnz = if isempty(nnz)
      Ptr{PetscInt}(0)
    else
      if length(nnz) != size(a,1)
        throw(ArgumentError("length(nnz) must be # rows"))
      end
      isa(nnz,Vector{PetscInt}) ? nnz : PetscInt[ i for i in nnz ]
    end
    chk(C.MatSeqBAIJSetPreallocation(a.p, bs, nz, pnnz))
  elseif MType == C.MATSEQSBAIJ
    pnnz = if isempty(nnz)
      Ptr{PetscInt}(0)
    else
      if length(nnz) != size(a,1)
        throw(ArgumentError("length(nnz) must be # rows"))
      end
      isa(nnz,Vector{PetscInt}) ? nnz : PetscInt[ i for i in nnz ]
    end
    chk(C.MatSeqSBAIJSetPreallocation(a.p, bs, nz, pnnz))
  else # TODO
    throw(ArgumentError("unsupported matrix type $T"))
  end
  a
end

# construct Vec for multiplication by a::Mat or transpose(a::Mat)
const mat2vec = Dict{C.MatType, C.MatType}( :mpiaij => :aij, :seqaij => :seq )
Vec{T2}(a::Mat{T2}, transposed=false) =
  transposed ? Vec(T, size(a,1), comm=comm(a), T=mat2vec[gettype(a)],
  mlocal=sizelocal(a,1)) :
  Vec(T, size(a,2), comm=comm(a), T=mat2vec[gettype(a)],
  mlocal=sizelocal(a,2))

#############################################################################
Base.convert(::Type{C.Mat}, a::Mat) = a.p

function Base.size(a::Mat)
  m = Ref{PetscInt}()
  n = Ref{PetscInt}()
  chk(C.MatGetSize(a.p, m, n))
  (Int(m[]), Int(n[]))
end

function sizelocal(a::Mat)
  m = Ref{PetscInt}()
  n = Ref{PetscInt}()
  chk(C.MatGetLocalSize(a.p, m, n))
  (Int(m[]), Int(n[]))
end

lengthlocal(a::Mat) = prod(sizelocal(a))

# this causes the assembly state of the underlying petsc matrix to be copied
function Base.similar{T, MType}(a::Mat{T, MType})
  p = Ref{C.Mat{T}}()
  chk(C.MatDuplicate(a.p, C.MAT_DO_NOT_COPY_VALUES, p))
  Mat{T, MType}(p[])
end

Base.similar{T}(a::Mat{T}, ::Type{T}) = similar(a)
Base.similar{T,MType}(a::Mat{T,MType}, T2::Type) =
  Mat(T2, size(a)..., comm=comm(a), mtype=MType)
Base.similar{T,MType}(a::Mat{T,MType}, T2::Type, m::Integer, n::Integer) =
  (m,n) == size(a) && T2==T ? similar(a) : Mat(T2, m,n, comm=comm(a), mtype=MType)
Base.similar{T,MType}(a::Mat{T,MType}, m::Integer, n::Integer) = similar(a, T, m, n)
Base.similar(a::Mat, T::Type, d::Dims) = similar(a, T, d...)
Base.similar{T}(a::Mat{T}, d::Dims) = similar(a, T, d)

function Base.copy{T,MType}(a::Mat{T,MType})
  p = Ref{C.Mat{T}}()
  chk(C.MatDuplicate(a.p, C.MAT_COPY_VALUES, p))
  Mat{T,MType}(p[])
end

function getinfo{T}(m::Mat{T}, infotype::Integer=C.MAT_GLOBAL_SUM)
  info = Ref{C.MatInfo{T}}()
  println(typeof(info))
  chk(C.MatGetInfo(m.p, C.MatInfoType(infotype), info))
  info[]
end

Base.nnz(m::Mat) = Int(getinfo(m).nz_used)

#############################################################################

# for efficient matrix assembly, put all calls to A[...] = ... inside
# assemble(A) do ... end

function AssemblyBegin(x::Mat, t::C.MatAssemblyType=C.MAT_FLUSH_ASSEMBLY)
  chk(C.MatAssemblyBegin(x.p, t))
end

function AssemblyEnd(x::Mat, t::C.MatAssemblyType=C.MAT_FLUSH_ASSEMBLY)
  chk(C.MatAssemblyEnd(x.p, t))
end

function isassembled(p::C.Mat)
  b = Ref{PetscBool}()
  chk(C.MatAssembled(p, b))
  return b[] != 0
end

isassembled(x::Mat) = !x.assembling && isassembled(x.p)

function assemble(f::Function, x::Union{Vec,Mat},
  insertmode=x.insertmode,
  assemblytype::C.MatAssemblyType=C.MAT_FINAL_ASSEMBLY)
  if x.insertmode != insertmode && x.assembling
    error("nested assemble with different insertmodes not allowed")
  end
  old_assembling = x.assembling
  old_insertmode = x.insertmode
  try
    x.assembling = true
    x.insertmode = insertmode
    result = f()
    if !old_assembling # don't Assemble if we are in nested assemble
      AssemblyBegin(x, assemblytype)
      yield() # do async computations while messages are in transit
      AssemblyEnd(x, assemblytype)
    end
    return result
  finally
    x.assembling = old_assembling
    x.insertmode = old_insertmode
  end
end

# force assembly even if it might not be necessary
function assemble(x::Union{Vec,Mat}, t::C.MatAssemblyType=C.MAT_FINAL_ASSEMBLY)
  AssemblyBegin(x, t)
  AssemblyEnd(x,t)
end


# in ksp solve we need to finalize assembly from raw pointer:
function assemble(p::C.Mat, t::C.MatAssemblyType=C.MAT_FINAL_ASSEMBLY)
  if !isassembled(p)
    chk(C.MatAssemblyBegin(p, t))
    chk(C.MatAssemblyEnd(p, t))
  end
  return nothing
end

# intermediate assembly, before it is finally compressed for use
iassemble(x::Union{Vec,Mat}) = assemble(() -> nothing, x, x.insertmode, C.MAT_FLUSH_ASSEMBLY)
iassemble(f::Function, x::Mat, insertmode=x.insertmode) =
  assemble(f, x, insertmode, C.MAT_FLUSH_ASSEMBLY)

# like x[i,j] = v, but requires i,j to be 0-based indices for Petsc
function setindex0!{T}(x::Mat{T}, v::Array{T},
                       i::Array{PetscInt}, j::Array{PetscInt})
  ni = length(i)
  nj = length(j)
  if length(v) != ni*nj
    throw(ArgumentError("length(values) != length(indices)"))
  end
  chk(C.MatSetValues(x.p, ni, i, nj, j, v, x.insertmode))
  if !x.assembling
    AssemblyBegin(x,C.MAT_FLUSH_ASSEMBLY)
    AssemblyEnd(x,C.MAT_FLUSH_ASSEMBLY)
  end
  x
end

import Base: setindex!

function setindex!{T}(x::Mat{T}, v::Number, i::Integer, j::Integer)
  # can't call MatSetValue since that is a static inline function
  setindex0!(x, T[ v ],
  PetscInt[ i - 1 ],
  PetscInt[ j - 1 ])
  v
end
function setindex!{T3, T1<:Integer, T2<:Integer}(x::Mat{T3}, v::Array{T3},
  I::AbstractArray{T1},
  J::AbstractArray{T2})
  I0 = PetscInt[ i-1 for i in I ]
  J0 = PetscInt[ j-1 for j in J ]
  setindex0!(x, v, I0, J0)
end
function setindex!{T2, T<:Integer}(x::Mat{T2}, v::Array{T2},
  i::Integer, J::AbstractArray{T})
  I0 = PetscInt[ i-1 ]
  J0 = PetscInt[ j-1 for j in J ]
  setindex0!(x, v, I0, J0)
end
function setindex!{T2, T<:Real}(x::Mat{T2}, v::Array{T2},
  I::AbstractArray{T}, j::Integer)
  I0 = PetscInt[ i-1 for i in I ]
  J0 = PetscInt[ j-1 ]
  setindex0!(x, v, I0, J0)
end

setindex!{T1<:Integer, T2<:Integer}(x::Mat, v::Number,
I::AbstractArray{T1},
J::AbstractArray{T2}) = iassemble(x) do
  for i in I
    for j in J
      x[i,j] = v
    end
  end
  x
end
setindex!{T<:Integer}(x::Mat, v::Number,
i::Integer, J::AbstractArray{T}) = iassemble(x) do
  for j in J
    x[i,j] = v
  end
  x
end
setindex!{T<:Integer}(x::Mat, v::Number,
I::AbstractArray{T}, j::Integer) = iassemble(x) do
  for i in I
    x[i,j] = v
  end
  x
end

function setindex!{T0<:Number, T1<:Integer, T2<:Integer}(x::Mat, v::AbstractArray{T0},
  I::AbstractArray{T1},
  J::AbstractArray{T2})
  if length(v) != length(I)*length(J)
    throw(ArgumentError("length(values) != length(indices)"))
  end
  v0 = eltype(Mat)[ z for z in v ]
  setindex!(x, v0, I, J)
end

# fill! is not a very sensible function to call for sparse matrices,
# but we might as well have it, especially for the v=0 case
function Base.fill!(x::Mat, v::Number)
  if v == 0
    chk(C.MatZeroEntries(x.p))
  else
    # FIXME: only loop over local rows
    iassemble(x) do
      m,n = size(x)
      for i in 1:m
        for j in 1:n
          x[i,j] = v
        end
      end
    end
  end
  return x
end

#############################################################################
import Base.getindex

# like getindex but for 0-based indices i and j
function getindex0{T}(x::Mat{T}, i::Vector{PetscInt}, j::Vector{PetscInt})
  ni = length(i)
  nj = length(j)
  v = Array(T, nj, ni) # row-major!
  chk(C.MatGetValues(x.p, ni, i, nj, j, v))
  ni <= 1 || nj <= 1 ? reshape(v, ni, nj) : transpose(v)
end

getindex(a::Mat, i0::Integer, i1::Integer) =
  getindex0(a, PetscInt[ i0-1 ], PetscInt[ i1-1 ])[1]

getindex{T0<:Integer,T1<:Integer}(a::Mat, I0::AbstractArray{T0}, I1::AbstractArray{T1}) =
  getindex0(a, PetscInt[ i0-1 for i0 in I0 ], PetscInt[ i1-1 for i1 in I1 ])

getindex{T0<:Integer}(a::Mat, I0::AbstractArray{T0}, i1::Integer) =
  reshape(getindex0(a, PetscInt[ i0-1 for i0 in I0 ], PetscInt[ i1-1 ]), length(I0))

getindex{T1<:Integer}(a::Mat, i0::Integer, I1::AbstractArray{T1}) =
  getindex0(a, PetscInt[ i0-1 ], PetscInt[ i1-1 for i1 in I1 ])

#############################################################################
# transposition etc.

export MatTranspose

function Base.full(a::Mat)
  m,n = size(a)
  a[1:m, 1:n]
end

# create a new matrix wrapping around A for matrix-vector multiplication
# but which does not actually require new storage
for (f,pf) in ((:MatTranspose,:MatCreateTranspose), # acts like A.'
  (:MatNormal, :MatCreateNormal))      # acts like A'*A
  pfe = Expr(:quote, pf)
  @eval function $f{T, MType}(a::Mat{T, MType})
    p = Ref{C.Mat{T}}()
    chk(C.$pf(a.p, p))
    Mat{T, MType}(p[], a)
  end
end

for (f,pf) in ((:transpose,:MatTranspose),(:ctranspose,:MatHermitianTranspose))
  fb = symbol(string(f,"!"))
  pfe = Expr(:quote, pf)
  @eval begin
    function Base.$fb(a::Mat)
      pa = [a.p]
      chk(C.$pf(a.p, C.MAT_REUSE_MATRIX, pa))
      a
    end

    function Base.$f{T,MType}(a::Mat{T,MType})
      p = Ref{C.Mat{T}}()
      chk(C.$pf(a.p, C.MAT_INITIAL_MATRIX, p))
      Mat{T,MType}(p[])
    end
  end
end

function Base.conj!(a::Mat)
  chk(C.MatConjugate(a.p))
  a
end
Base.conj(a::Mat) = conj!(copy(a))

#############################################################################
# simple math operations

# skip for now
#=
function chop!(x::Mat, tol::Real)
chk(ccall((:MatChop, petsc), PetscErrorCode, (pMat, PetscReal), x, tol))
x
end
=#

import Base: .*, ./, .\, *, +, -, ==
import Base.LinAlg: At_mul_B, At_mul_B!, Ac_mul_B, Ac_mul_B!, A_mul_Bt, A_mul_Bt!

function Base.trace{T}(A::Mat{T})
  t = Ref{T}()
  chk(C.MatGetTrace(A.p,t))
  return t[]
end

function Base.real{T<:Complex}(A::Mat{T})
  N = copy(A)
  chk(C.MatRealPart(N.p))
  return N
end
Base.real{T<:Real}(A::Mat{T}) = A

function Base.imag{T<:Complex}(A::Mat{T})
  N = copy(A)
  chk(C.MatImaginaryPart(N.p))
  return N
end

function Base.LinAlg.ishermitian{T}(A::Mat{T}, tol::Real=eps(real(float(one(T)))))
  bool_arr = Ref{PetscBool}()
  chk(C.MatIsHermitian(A.p, tol, bool_arr))
  return bool_arr[] != 0
end

function Base.LinAlg.issym{T}(A::Mat{T}, tol::Real=eps(real(float(one(T)))))
  bool_arr = Ref{PetscBool}()
  chk(C.MatIsSymmetric(A.p, tol, bool_arr))
  return bool_arr[] != 0
end

#currently ONLY gets the main diagonal
function Base.diag{T}(A::Mat{T},vtype::C.VecType=C.VECSEQ)
  m = size(A, 1)
  b = Vec(T, m, vtype, comm=comm(A), mlocal=sizelocal(A,1))
  chk(C.MatGetDiagonal(A.p,b.p))
  return b
end

function (*){T, MType}(A::Mat{T}, x::Vec{T, MType})
  m = size(A, 1)
  b = Vec(T, m, MType, comm=comm(A), mlocal=sizelocal(A,1))
  chk(C.MatMult(A.p, x.p, b.p))
  return b
end

function (.*){T, MType}(A::Mat{T, MType}, x::Number)
  Y = copy(A)
  chk(C.MatScale(Y.p, T(x)))
  return Y
end
(.*){T, MType}(x::Number, A::Mat{T, MType}) = A .* x
(./){T, MType}(A::Mat{T, MType}, x::Number) = A .* inv(x)
(.\){T, MType}(x::Number, A::Mat{T, MType}) = A .* inv(x)

function (*){T, MType}(A::Mat{T,MType}, B::Mat{T})
  p = Ptr{Float64}(0)
  p_arr = Ref{C.Mat{T}}()
  chk(C.MatMatMult(A.p, B.p, C.MAT_INITIAL_MATRIX, real(T)(C.PETSC_DEFAULT), p_arr))
  new_mat = Mat{T, MType}(p_arr[])
  return new_mat
end

#these two only work for SEQAIJ
function At_mul_B{T}(A::Mat{T}, B::Mat{T})
  p = Ptr{Float64}(0)
  p_arr = Ref{C.Mat{T}}()
  chk(C.MatTransposeMatMult(A.p, B.p, C.MAT_INITIAL_MATRIX, real(T)(C.PETSC_DEFAULT), p_arr))
  new_mat = Mat{T, C.MATSEQAIJ}(p_arr[])
  return new_mat
end

function A_mul_Bt{T}(A::Mat{T}, B::Mat{T})
  p = Ptr{Float64}(0)
  p_arr = Ref{C.Mat{T}}()
  chk(C.MatMatTransposeMult(A.p, B.p, C.MAT_INITIAL_MATRIX, real(T)(C.PETSC_DEFAULT), p_arr))
  new_mat = Mat{T, C.MATSEQAIJ}(p_arr[])
  return new_mat
end

function At_mul_B!{T, MType}(A::Mat{T}, x::Vec{T,MType}, y::Vec{T,MType})
  chk(C.MatMultTranspose(A.p, x.p, y.p))
  return y
end

function At_mul_B{T, MType}(A::Mat{T}, x::Vec{T,MType})
  m = size(A, 1)
  b = Vec(T, m, MType, comm=comm(A), mlocal=sizelocal(A,1))
  chk(C.MatMultTranspose(A.p, x.p, b.p))
  return b
end

function Ac_mul_B!{T, MType}(A::Mat{T}, x::Vec{T,MType}, y::Vec{T,MType})
  chk(C.MatMultHermitianTranspose(A.p, x.p, y.p))
  return y
end

function Ac_mul_B{T, MType}(A::Mat{T}, x::Vec{T,MType})
  m = size(A, 1)
  b = Vec(T, m, MType, comm=comm(A), mlocal=sizelocal(A,1))
  chk(C.MatMultHermitianTranspose(A.p, x.p, b.p))
  return b
end

for (f,s) in ((:+,1), (:-,-1))
  @eval function ($f){T}(A::Mat{T}, B::Mat{T})
    Y = copy(A)
    chk(C.MatAXPY(Y.p, T($s), B.p, C.DIFFERENT_NONZERO_PATTERN))
    return Y
  end
end

(-){T, MType}(A::Mat{T, MType}) = A * (-1)

# there don't appear to be PETSc functions for pointwise
# operations on matrices

function (==){T}(A::Mat{T}, b::Mat{T})
  bool_arr = Ref{PetscBool}()
  chk(C.MatEqual(A.p, b.p, bool_arr))
  return bool_arr[] != 0
end
