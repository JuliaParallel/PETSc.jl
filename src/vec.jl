# AbstractVector wrapper around PETSc Vec types
export Vec, comm

type Vec{T,VType} <: AbstractVector{T}
  p::C.Vec{T}
  assembling::Bool # whether we are in the middle of assemble(vec)
  insertmode::C.InsertMode # current mode for setindex!
  data::Any # keep a reference to anything needed for the Mat
            # -- needed if the Mat is a wrapper around a Julia object,
            #    to prevent the object from being garbage collected.
  function Vec(p::C.Vec{T}, data=nothing)
    v = new(p, false, C.INSERT_VALUES, data)
    chk(C.VecSetType(p, VType))  # set the type here to ensure it matches VType
    finalizer(v, VecDestroy)
    return v
  end
end

comm{T}(v::Vec{T}) = MPI.Comm(C.PetscObjectComm(T, v.p.pobj))

function Vec{T}(::Type{T}, vtype::C.VecType=C.VECSEQ; comm::MPI.Comm=MPI.COMM_SELF)
  p = Ref{C.Vec{T}}()
  chk(C.VecCreate(comm, p))
  Vec{T, vtype}(p[])
end

function Vec{T<:Scalar}(::Type{T}, len::Integer, vtype::C.VecType=C.VECSEQ;
                        comm::MPI.Comm=MPI.COMM_SELF, mlocal::Integer=C.PETSC_DECIDE)
  vec = Vec(T, vtype; comm=comm)
  resize!(vec, len, mlocal=mlocal)
  vec
end

# make a Vec that is a wrapper around v, where v stores the local data
function Vec{T<:Scalar}(v::Vector{T}; comm::MPI.Comm=MPI.COMM_SELF)
  p = Ref{C.Vec{T}}()
  chk(C.VecCreateMPIWithArray(comm, 1, length(v), C.PETSC_DECIDE, v, p))
  pv = Vec{T, C.VECMPI}(p[], v)
  return pv
end

function VecDestroy{T}(vec::Vec{T})
  PetscFinalized(T) || C.VecDestroy(Ref(vec.p))
end

function petscview{T}(vec::Vec{T})
  viewer = C.PetscViewer{T}(C_NULL)
  chk(C.VecView(vec.p, viewer))
end

export gettype

gettype{T,VT}(a::Vec{T,VT}) = VT

function Base.resize!(x::Vec, m::Integer=C.PETSC_DECIDE; mlocal::Integer=C.PETSC_DECIDE)
  if m == mlocal == C.PETSC_DECIDE
    throw(ArgumentError("either the length (m) or local length (mlocal) must be specified"))
  end

  chk(C.VecSetSizes(x.p, mlocal, m))
  x
end

##########################################################################
export lengthlocal, sizelocal, localpart

Base.convert(::Type{C.Vec}, v::Vec) = v.p

function Base.length(x::Vec)
  sz = Ref{PetscInt}()
  chk(C.VecGetSize(x.p, sz))
  Int(sz[])
end
Base.size(x::Vec) = (length(x),)

function lengthlocal(x::Vec)
  sz = Ref{PetscInt}()
  chk(C.VecGetLocalSize(x.p, sz))
  sz[]
end
sizelocal(x::Vec) = (lengthlocal(x),)
sizelocal{T,n}(t::AbstractArray{T,n}, d) = (d>n ? 1 : sizelocal(t)[d])

function localpart(v::Vec)
  # this function returns a range from the first to the last indicies (1 based)
  # this is different than the Petsc VecGetOwnershipRange function where
  # the max value is one more than the number of entries
  low = Ref{PetscInt}()
  high = Ref{PetscInt}()
  chk(C.VecGetOwnershipRange(v.p, low, high))
  return (low[]+1):(high[])
end

function Base.similar{T,VType}(x::Vec{T,VType})
  p = Ref{C.Vec{T}}()
  chk(C.VecDuplicate(x.p, p))
  Vec{T,VType}(p[])
end

Base.similar{T}(x::Vec{T}, ::Type{T}) = similar(x)
Base.similar{T,VType}(x::Vec{T,VType}, T2::Type) =
  Vec(T2, length(x), VType; comm=comm(x), mlocal=lengthlocal(x))

function Base.similar{T,VType}(x::Vec{T,VType}, T2::Type, len::Union{Int,Dims})
  length(len) == 1 || throw(ArgumentError("expecting 1-dimensional size"))
  len[1]==length(x) && T2==T ? similar(x) : Vec(T2, len[1], VType; comm=comm(x))
end

function Base.similar{T,VType}(x::Vec{T,VType}, len::Union{Int,Dims})
  length(len) == 1 || throw(ArgumentError("expecting 1-dimensional size"))
  len[1]==length(x) ? similar(x) : Vec(T, len[1], VType; comm=comm(x))
end

function Base.copy(x::Vec)
  y = similar(x)
  chk(C.VecCopy(x.p, y.p))
  y
end

##########################################################################
import Base: setindex!
export assemble, isassembled

# for efficient vector assembly, put all calls to x[...] = ... inside
# assemble(x) do ... end
function AssemblyBegin(x::Vec, t::C.MatAssemblyType=C.MAT_FINAL_ASSEMBLY)
  # the t parameter is unused for vectors
  chk(C.VecAssemblyBegin(x.p))
end

function AssemblyEnd(x::Vec, t::C.MatAssemblyType=C.MAT_FINAL_ASSEMBLY)
  chk(C.VecAssemblyEnd(x.p))
end

isassembled(x::Vec) = !x.assembling
# assemble(f::Function, x::Vec) is defined in mat.jl

# like x[i] = v, but requires i to be 0-based indices for Petsc
function setindex0!{T}(x::Vec{T}, v::Array{T}, i::Array{PetscInt})
  n = length(v)
  if n != length(i)
    throw(ArgumentError("length(values) != length(indices)"))
  end
  #    println("  in setindex0, passed bounds check")
  chk(C.VecSetValues(x.p, n, i, v, x.insertmode))
  if !x.assembling
    AssemblyBegin(x)
    AssemblyEnd(x)
  end
  x
end

function setindex!{T}(x::Vec{T}, v::Number, i::Integer)
  # can't call VecSetValue since that is a static inline function
  setindex0!(x, T[ v ], PetscInt[ i - 1 ])
  v
end

# set multiple entries to a single value
setindex!{T<:Integer}(x::Vec, v::Number, I::AbstractArray{T}) = assemble(x) do
  for i in I
    x[i] = v
  end
  x
end

function Base.fill!{T}(x::Vec{T}, v::Number)
  chk(C.VecSet(x.p, T(v)))
  return x
end

function setindex!{T<:Integer}(x::Vec, v::Number, I::Range{T})
  if abs(step(I)) == 1 && minimum(I) == 1 && maximum(I) == length(x)
    fill!(x, v)
    return v
  else
    # use invoke here to avoid a recursion loop
    return invoke(setindex!, (Vec,typeof(v),AbstractVector{T}), x,v,I)
  end
end

setindex!{T<:Real}(x::Vec, V::AbstractArray, I::AbstractArray{T}) =
assemble(x) do
  if length(V) != length(I)
    throw(ArgumentError("length(values) != length(indices)"))
  end
  # possibly faster to make a PetscScalar array from V, and
  # a copy of the I array shifted by 1, to call setindex0! instead?
  c = 1
  for i in I
    x[i] = V[c]
    c += 1
  end
  x
end

# logical indexing
setindex!(A::Vec, x::Number, I::AbstractArray{Bool}) = assemble(A) do
  for i = 1:length(I)
    if I[i]
      A[i] = x
    end
  end
  A
end
for T in (:(Array{T2}),:(AbstractArray{T2})) # avoid method ambiguities
  @eval setindex!{T2<:Scalar}(A::Vec, X::$T, I::AbstractArray{Bool}) = assemble(A) do
    c = 1
    for i = 1:length(I)
      if I[i]
        A[i] = X[c]
        c += 1
      end
    end
    A
  end
end

##########################################################################
import Base.getindex

# like getindex but for 0-based indices i
function getindex0{T}(x::Vec{T}, i::Vector{PetscInt})
  v = similar(i, T)
  chk(C.VecGetValues(x.p, length(i), i, v))
  v
end

getindex(x::Vec, i::Integer) = getindex0(x, PetscInt[i-1])[1]

getindex(x::Vec, I::AbstractVector{PetscInt}) =
  getindex0(x, PetscInt[ (i-1) for i in I ])

##########################################################################

import Base: abs, exp, log, conj, conj!
export abs!, exp!, log!
for (f,pf) in ((:abs,:VecAbs), (:exp,:VecExp), (:log,:VecLog),
  (:conj,:VecConjugate))
  fb = symbol(string(f, "!"))
  @eval begin
    function $fb(x::Vec)
      chk(C.$pf(x.p))
      x
    end
    $f(x::Vec) = $fb(copy(x))
  end
end

# skip for now
#=
export chop!
function chop!(x::Vec, tol::Real)
VecChop(x.c, tol)
chk(ccall((:VecChop, petsc), PetscErrorCode, (pVec, PetscReal), x, tol))
x
end
=#

for (f, pf, sf) in ((:findmax, :VecMax, :maximum), (:findmin, :VecMin, :minimum))
  @eval begin
    function Base.$f{T<:Real}(x::Vec{T})
      i = Ref{PetscInt}()
      v = Ref{T}()
      chk(C.$pf(x.p, i, v))
      (v[], i[]+1)
    end
    Base.$sf{T<:Real}(x::Vec{T}) = $f(x)[1]
  end
end
# For complex numbers, VecMax and VecMin apparently return the max/min
# real parts, which doesn't match Julia's maximum/minimum semantics.

function Base.norm{T<:Real}(x::Union{Vec{T},Vec{Complex{T}}}, p::Number)
  v = Ref{T}()
  n = p == 1 ? C.NORM_1 : p == 2 ? C.NORM_2 : p == Inf ? C.NORM_INFINITY :
  throw(ArgumentError("unrecognized Petsc norm $p"))
  chk(C.VecNorm(x.p, n, v))
  v[]
end

if VERSION >= v"0.5.0-dev+8353" # JuliaLang/julia#13681
  import Base.normalize!
else
  export normalize!
end

# computes v = norm(x,2), divides x by v, and returns v
function normalize!{T<:Real}(x::Union{Vec{T},Vec{Complex{T}}})
  v = Ref{T}()
  chk(C.VecNormalize(x.p, v))
  v[]
end

function Base.dot{T}(x::Vec{T}, y::Vec{T})
  d = Ref{T}()
  chk(C.VecDot(y.p, x.p, d))
  return d[]
end

# unconjugated dot product (called for x'*y)
function Base.At_mul_B{T<:Complex}(x::Vec{T}, y::Vec{T})
  d = Array(T, 1)
  chk(C.VecTDot(x.p, y.p, d))
  return d
end

# pointwise operations on pairs of vectors (TODO: support in-place variants?)
import Base: max, min, .*, ./, .\
for (f,pf) in ((:max,:VecPointwiseMax), (:min,:VecPointwiseMin),
  (:.*,:VecPointwiseMult), (:./,:VecPointwiseDivide))
  @eval function ($f)(x::Vec, y::Vec)
    w = similar(x)
    chk(C.$pf(w.p, x.p, y.p))
    w
  end
end

import Base: +, -
function Base.scale!{T}(x::Vec{T}, s::Number)
  chk(C.VecScale(x.p, T(s)))
  x
end
Base.scale{T}(x::Vec{T},s::Number) = scale!(copy(x),s)
(.*)(x::Vec, a::Number...) = scale(x, prod(a))
(.*)(a::Number, x::Vec) = scale(x, a)
(./)(x::Vec, a::Number) = scale(x, inv(a))
(.\)(a::Number, x::Vec) = scale(x, inv(a))
function (./)(a::Number, x::Vec)
  y = copy(x)
  chk(C.VecReciprocal(y.p))
  if a != 1.0
    scale!(y, a)
  end
  y
end

function (+){T<:Scalar}(x::Vec{T}, a::Number...)
  y = copy(x)
  chk(C.VecShift(y.p, T(sum(a))))
  return y
end
(+){T<:Scalar}(a::Number, x::Vec{T}) = x + a
(-){T<:Scalar}(x::Vec{T}, a::Number) = x + (-a)
(-)(x::Vec) = scale(x, -1)
function (-){T<:Scalar}(a::Number, x::Vec{T})
  y = -x
  chk(C.VecShift(y.p, T(a)))
  return y
end

import Base: ==
function (==)(x::Vec, y::Vec)
  b = Ref{PetscBool}()
  chk(C.VecEqual(x.p, y.p, b))
  b[] != 0
end

function (==)(x::Vec, y::AbstractArray)
  flag = true
  for i=1:length(x)
    flag = flag && x[i] == y[i]
  end

  return flag
end

function Base.sum{T}(x::Vec{T})
  s = Ref{T}()
  chk(C.VecSum(x.p, s))
  s[]
end

##########################################################################
export axpy!, aypx!, axpby!, axpbypcz!
import Base.LinAlg.BLAS.axpy!

# y <- alpha*x + y
function axpy!{T}(alpha::Number, x::Vec{T}, y::Vec{T})
  chk(C.VecAXPY(y.p, T(alpha), x.p))
  y
end
# w <- alpha*x + y
function axpy!{T}(alpha::Number, x::Vec{T}, y::Vec{T}, w::Vec{T})
  chk(C.VecWAXPY(w.p, T(alpha), x.p, y.p))
  y
end
# y <- alpha*y + x
function aypx!{T}(x::Vec{T}, alpha::Number, y::Vec{T})
  chk(C.VecAYPX( y.p, T(alpha), x.p))
  y
end
# y <- alpha*x + beta*y
function axpby!{T}(alpha::Number, x::Vec{T}, beta::Number, y::Vec{T})
  chk(C.VecAXPBY(y.p, T(alpha), T(beta), x.p))
  y
end
# z <- alpha*x + beta*y + gamma*z
function axpbypcz!{T}(alpha::Number, x::Vec{T}, beta::Number, y::Vec{T},
  gamma::Number, z::Vec{T})
  chk(C.VecAXPBYPCZ(z.p, T(alpha), T(beta), T(gamma), x.p, y.p))
  z
end

# y <- y + \sum_i alpha[i] * x[i]
function axpy!{V<:Vec}(y::V, alpha::AbstractArray, x::AbstractArray{V})
  n = length(x)
  length(alpha) == n || throw(BoundsError())
  _x = [X.p for X in x]
  _alpha = eltype(y)[a for a in alpha]
  C.VecMAXPY(y.p, n, _alpha, _x)
  y
end

##########################################################################
# element-wise vector operations:
import Base: .*, ./, .^, +, -

for (f,pf) in ((:.*,:VecPointwiseMult), (:./,:VecPointwiseDivide), (:.^,:VecPow))
  @eval function ($f)(x::Vec, y::Vec)
    z = similar(x)
    chk(C.$pf(z.p, x.p, y.p))
    return z
  end
end

for (f,s) in ((:+,1), (:-,-1))
  @eval function ($f){T}(x::Vec{T}, y::Vec{T})
    z = similar(x)
    chk(C.VecWAXPY(z.p, T($s), y.p, x.p))
    return z
  end
end


##############################################################################
export VecArray, VecArrayRestore
type VecArray{T <: Scalar} <: AbstractArray{T, 1}
  a::Array{T, 1}  # the array object constructed around the pointer
  ref::Ref{Ptr{T}}  # reference to the pointer to the data
  pobj::C.Vec{T}
  isfinalized::Bool  # has this been finalized yet
  function VecArray(a::Array, ref::Ref, ptr)
    varr = new(a, ref, ptr, false)
    # backup finalizer, shouldn't ever be used because users must call
    # VecArrayRestore before their changes will take effect
    finalizer(varr, VecArrayRestore)
    return varr
  end

end

function VecArray{T}(vec::Vec{T})

  len = length(vec)

  ref = Ref{Ptr{T}}()
  chk(C.VecGetArray(vec.p, ref))
  a = pointer_to_array(ref[], len)
  return VecArray{T}(a, ref, vec.p)
end

#TODO: check if pobj was already finalized
function VecArrayRestore{T}(varr::VecArray{T})

  if !varr.isfinalized && !PetscFinalized(T) 
    ptr = [varr.ref[]]
    chk(C.VecRestoreArray(varr.pobj, ptr))
  end 
  varr.isfinalized = true
end

Base.linearindexing(varr::VecArray) = Base.linearindexing(varr.a)
Base.size(varr::VecArray) = size(varr.a)

# indexing
getindex(varr::VecArray, i) = getindex(varr.a, i)
setindex!(varr::VecArray, v, i) = setindex!(varr.a, v, i)

Base.copy(varr::VecArray) = deepcopy(varr)

function (==)(x::VecArray, y::AbstractArray)
  flag = true
  for i=1:length(x)
    flag = flag && x[i] == y[i]
  end

  return flag
end

# what to do about similar?  it shouldn't be used?
