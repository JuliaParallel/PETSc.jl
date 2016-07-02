# AbstractVector wrapper around PETSc Vec types
export Vec, comm, NullVec

@doc """
  Construct a high level Vec object from a low level C.Vec.
  The data field is used to protect things from GC.
  A finalizer is attached to deallocate the memory of the underlying C.Vec, unless 
  first_instance is set to true.
  If assembling is false, then every call to setindex! also calls AssembleBegin and 
  AssembleEnd.
"""
type Vec{T,VType} <: AbstractVector{T}
  p::C.Vec{T}
  assembling::Bool # whether we are in the middle of assemble(vec)
  insertmode::C.InsertMode # current mode for setindex!
  data::Any # keep a reference to anything needed for the Mat
            # -- needed if the Mat is a wrapper around a Julia object,
            #    to prevent the object from being garbage collected.
  function Vec(p::C.Vec{T}, data=nothing; first_instance::Bool=true)
    v = new(p, false, C.INSERT_VALUES, data)
    if first_instance
      chk(C.VecSetType(p, VType))  # set the type here to ensure it matches VType
      finalizer(v, PetscDestroy)
    end
    return v
  end
end


@doc """
  Null vectors, used in place of void pointers in the C
  API
"""
global const NullVec = Dict{DataType, Vec}()


if have_petsc[1]
  global const NullVec1 = Vec{Float64, C.VECSTANDARD}(C.Vec{Float64}(C_NULL), first_instance=false)
  NullVec[Float64] = NullVec1
end
if have_petsc[2]
  global const NullVec2 = Vec{Float32, C.VECSTANDARD}(C.Vec{Float32}(C_NULL), first_instance=false)
  NullVec[Float32] = NullVec2
end
if have_petsc[3]
  global const NullVec3 = Vec{Complex128, C.VECSTANDARD}(C.Vec{Complex128}(C_NULL), first_instance=false)
  NullVec[Complex128] = NullVec3

end
@doc """
  Gets the MPI communicator of a vector.
"""
comm{T}(v::Vec{T}) = MPI.Comm(C.PetscObjectComm(T, v.p.pobj))


export gettype

@doc """
  Get the symbol that is the format of the vector
"""
gettype{T,VT}(a::Vec{T,VT}) = VT


@doc """
  Create an empty, unsized vector.
"""
function Vec{T}(::Type{T}, vtype::C.VecType=C.VECMPI; 
                comm::MPI.Comm=MPI.COMM_WORLD)
  p = Ref{C.Vec{T}}()
  chk(C.VecCreate(comm, p))
  Vec{T, vtype}(p[])
end

@doc """
  Create a vector, specifying the (global) length len or the local length
  mlocal
"""
function Vec{T<:Scalar}(::Type{T}, len::Integer=C.PETSC_DECIDE;
                         vtype::C.VecType=C.VECMPI,
                         comm::MPI.Comm=MPI.COMM_WORLD, 
                         mlocal::Integer=C.PETSC_DECIDE)
  vec = Vec(T, vtype; comm=comm)
  resize!(vec, len, mlocal=mlocal)
  vec
end

@doc """
  Make a PETSc vector out of an array.  If used in parallel, the array becomes
  the local part of the PETSc vector
"""
# make a Vec that is a wrapper around v, where v stores the local data
function Vec{T<:Scalar}(v::Vector{T}; comm::MPI.Comm=MPI.COMM_WORLD)
  p = Ref{C.Vec{T}}()
  chk(C.VecCreateMPIWithArray(comm, 1, length(v), C.PETSC_DECIDE, v, p))
  pv = Vec{T, C.VECMPI}(p[], v)
  return pv
end

export VecGhost, VecLocal, restore


@doc """
  Make a PETSc vector with space for ghost values.  ghost_idx are the 
  global indices that will be copied into the ghost space.
"""
# making mlocal the position and mglobal the keyword argument is inconsistent
# with the other Vec constructors, but it makes more sense here
function VecGhost{T<:Scalar, I <: Integer}(::Type{T}, mlocal::Integer, 
                  ghost_idx::Array{I,1}; comm=MPI.COMM_WORLD, m=C.PETSC_DECIDE, bs=1)

    nghost = length(ghost_idx)
    ghost_idx2 = [ PetscInt(i -1) for i in ghost_idx]

    vref = Ref{C.Vec{T}}()
    if bs == 1
      chk(C.VecCreateGhost(comm, mlocal, m, nghost, ghost_idx2, vref))
    elseif bs > 1
      chk(C.VecCreateGhostBlock(comm, bs, mlocal, mlocal, m, nghost, ghost_idx2, vref))
    else
      println(STDERR, "WARNING: unsupported block size requested, bs = ", bs)
    end

    return Vec{T, C.VECMPI}(vref[])
end

@doc """
  Create a VECSEQ that contains both the local and the ghost values of the 
  original vector.  The underlying memory for the orignal and output vectors
  alias.
"""
function VecLocal{T <:Scalar}( v::Vec{T, C.VECMPI})

  vref = Ref{C.Vec{T}}()
  chk(C.VecGhostGetLocalForm(v.p, vref))
  # store v to use with Get/Restore LocalForm
  # Petsc reference counting solves the gc problem
  return Vec{T, C.VECSEQ}(vref[], v)
end

#TODO: use restore for all types of restoring a local view
@doc """
  Tell Petsc the VecLocal is no longer needed
"""
function restore{T}(v::Vec{T, C.VECSEQ})

  vp = v.data
  vref = Ref(v.p)
  chk(C.VecGhostRestoreLocalForm(vp.p, vref))
end


@doc """
  The Petsc function to deallocate Vec objects
"""
function PetscDestroy{T}(vec::Vec{T})
  if !PetscFinalized(T)  && !isfinalized(vec)
    C.VecDestroy(Ref(vec.p))
    vec.p = C.Vec{T}(C_NULL)  # indicate the vector is finalized
  end
end

@doc """
  Determine whether a vector has already been finalized
"""
function isfinalized(vec::Vec)
  return isfinalized(vec.p)
end

function isfinalized(vec::C.Vec)
  return vec.pobj == C_NULL
end

global const is_nullvec = isfinalized  # another name for doing the same check

@doc """
  Use the PETSc routine for printing a vector to stdout
"""
function petscview{T}(vec::Vec{T})
  viewer = C.PetscViewer{T}(C_NULL)
  chk(C.VecView(vec.p, viewer))
end

function Base.resize!(x::Vec, m::Integer=C.PETSC_DECIDE; mlocal::Integer=C.PETSC_DECIDE)
  if m == mlocal == C.PETSC_DECIDE
    throw(ArgumentError("either the length (m) or local length (mlocal) must be specified"))
  end

  chk(C.VecSetSizes(x.p, mlocal, m))
  x
end

###############################################################################
export ghost_begin!, ghost_end!, scatter!, ghost_update!
# ghost vectors: essential methods
@doc """
  Start communication to update the ghost values (on other processes) from the local
  values
"""
function ghost_begin!{T<:Scalar}(v::Vec{T, C.VECMPI}; imode=C.INSERT_VALUES,
                               smode=C.SCATTER_FORWARD)
    chk(C.VecGhostUpdateBegin(v.p, imode, smode))
    return v
end

@doc """
  Finish communication for updating ghost values
"""
function ghost_end!{T<:Scalar}(v::Vec{T, C.VECMPI}; imode=C.INSERT_VALUES,
                               smode=C.SCATTER_FORWARD)
    chk(C.VecGhostUpdateEnd(v.p, imode, smode))
    return v
end

# ghost vectors: helpful methods
@doc """
  Convenience method for calling both ghost_begin! and ghost_end!
"""
function scatter!{T<:Scalar}(v::Vec{T, C.VECMPI}; imode=C.INSERT_VALUES, smode=C.SCATTER_FORWARD)

  ghost_begin!(v, imode=imode, smode=smode)
  ghost_end!(v, imode=imode, smode=smode)
end

# is there a way to specify all varargs must be same type?
# this can't be named scatter! because of ambiguity with the index set scatter!
@doc """
  Convenience method for calling ghost_begin! and ghost_end! for multiple vectors
"""
function ghost_update!(v...; imode=C.INSERT_VALUES, smode=C.SCATTER_FORWARD)

  for i in v
    ghost_begin!(i, imode=imode, smode=smode)
  end

  for i in v
    ghost_end!(i, imode=imode, smode=smode)
  end

  return v
end



###############################################################################
export lengthlocal, sizelocal, localpart

Base.convert(::Type{C.Vec}, v::Vec) = v.p

import Base.length
@doc """
  Get the global length of the vector
"""
function length(x::Vec)
  sz = Ref{PetscInt}()
  chk(C.VecGetSize(x.p, sz))
  Int(sz[])
end

@doc """
  Get the global size of the vector
"""
Base.size(x::Vec) = (length(x),)

@doc """
  Get the length of the local portion of the vector
"""
function lengthlocal(x::Vec)
  sz = Ref{PetscInt}()
  chk(C.VecGetLocalSize(x.p, sz))
  sz[]
end

"""
  Get the local size of the vector
"""
sizelocal(x::Vec) = (lengthlocal(x),)

"""
  Get local size of the vector
"""
sizelocal{T,n}(t::AbstractArray{T,n}, d) = (d>n ? 1 : sizelocal(t)[d])

@doc """
  Get the range of global indices that define the local part of the vector.
  Internally, this calls the Petsc function VecGetOwnershipRange, and has
  the same limitations as that function, namely that some vector formats do 
  not have a well defined contiguous range.
"""
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
  len[1]==length(x) && T2==T ? similar(x) : Vec(T2, len[1], vtype=VType; comm=comm(x))
end

function Base.similar{T,VType}(x::Vec{T,VType}, len::Union{Int,Dims})
  length(len) == 1 || throw(ArgumentError("expecting 1-dimensional size"))
  len[1]==length(x) ? similar(x) : Vec(T, len[1], vtype=VType; comm=comm(x))
end

function Base.copy(x::Vec)
  y = similar(x)
  chk(C.VecCopy(x.p, y.p))
  y
end

##########################################################################
import Base: setindex!
export assemble, isassembled, AssemblyBegin, AssemblyEnd

# for efficient vector assembly, put all calls to x[...] = ... inside
# assemble(x) do ... end
@doc """
  Start communication to assemble stashed values into the vector
"""
function AssemblyBegin(x::Vec, t::C.MatAssemblyType=C.MAT_FINAL_ASSEMBLY)
  # the t parameter is unused for vectors
  chk(C.VecAssemblyBegin(x.p))
end

@doc """
  Finish communication for assembling the vector
"""
function AssemblyEnd(x::Vec, t::C.MatAssemblyType=C.MAT_FINAL_ASSEMBLY)
  chk(C.VecAssemblyEnd(x.p))
end

@doc """
  Check if a vector is assembled (ie. does not have stashed values)
"""
isassembled(x::Vec) = !x.assembling
# assemble(f::Function, x::Vec) is defined in mat.jl

@doc """
  Like setindex, but requires the indices be 0-base
"""
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

#TODO: make this a single call to VecSetValues
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

export chop!
function chop!(x::Vec, tol::Real)
  chk(C.VecChop(x.p, tol))
#  chk(ccall((:VecChop, petsc), PetscErrorCode, (pVec, PetscReal), x, tol))
  x
end

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

@doc """
  computes v = norm(x,2), divides x by v, and returns v
"""
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
  x_arr = LocalVector(x) 
  for i=1:length(x_arr)  # do localpart, then MPI reduce
    flag = flag && x_arr[i] == y[i]
  end
  restore(x_arr)

  buf = Int8[flag]
  # process 0 is root
  recbuf = MPI.Reduce(buf, 1, MPI.LAND, 0, comm(x))

  if  MPI.Comm_rank(comm(x)) == 0
    buf[1] = recbuf[1]
  end

  MPI.Bcast!(buf, 1, 0, comm(x))
 
  return convert(Bool, buf[1]) 
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
export LocalVector, LocalVector_readonly, restore

@doc """
  Object representing the local part of the array, accessing the memory directly.
  Supports all the same indexing as a regular Array
"""
type LocalVector{T <: Scalar, ReadOnly} <: DenseArray{T, 1}
  a::Array{T, 1}  # the array object constructed around the pointer
  ref::Ref{Ptr{T}}  # reference to the pointer to the data
  pobj::C.Vec{T}
  isfinalized::Bool  # has this been finalized yet
  function LocalVector(a::Array, ref::Ref, ptr)
    varr = new(a, ref, ptr, false)
    # backup finalizer, shouldn't ever be used because users must call
    # restore before their changes will take effect
    finalizer(varr, restore)
    return varr
  end

end

typealias LocalVectorRead{T} LocalVector{T, true}
typealias LocalVectorWrite{T} LocalVector{T, false}

@doc """
  Get the LocalArray of a vector.  Users must call LocalArrayRestore when
  finished updating the vector
"""
function LocalVector{T}(vec::Vec{T})

  len = lengthlocal(vec)

  ref = Ref{Ptr{T}}()
  chk(C.VecGetArray(vec.p, ref))
  a = pointer_to_array(ref[], len)
  return LocalVector{T, false}(a, ref, vec.p)
end

@doc """
  Tell Petsc the LocalArray is no longer being used
"""
function restore{T}(varr::LocalVectorWrite{T})

  if !varr.isfinalized && !PetscFinalized(T) && !isfinalized(varr.pobj)
    ptr = varr.ref
    chk(C.VecRestoreArray(varr.pobj, ptr))
  end 
  varr.isfinalized = true
end

@doc """
  Get a readonly version of LocalVector.  Users must call restore when 
  finished with the object.
"""
function LocalVector_readonly{T}(vec::Vec{T})

  len = lengthlocal(vec)

  ref = Ref{Ptr{T}}()
  chk(C.VecGetArrayRead(vec.p, ref))
  a = pointer_to_array(ref[], len)
  return LocalVector{T, true}(a, ref, vec.p)
end

function restore{T}(varr::LocalVectorRead{T})

  if !varr.isfinalized && !PetscFinalized(T) && !isfinalized(varr.pobj)
    ptr = [varr.ref[]]
    chk(C.VecRestoreArrayRead(varr.pobj, ptr))
  end 
  varr.isfinalized = true
end

Base.size(varr::LocalVector) = size(varr.a)
# indexing
getindex(varr::LocalVector, i) = getindex(varr.a, i)
setindex!(varr::LocalVectorWrite, v, i) = setindex!(varr.a, v, i)
Base.unsafe_convert{T}(::Type{Ptr{T}}, a::LocalVector{T}) = Base.unsafe_convert(Ptr{T}, a.a)
Base.stride(a::LocalVector, d::Integer) = stride(a.a, d)
Base.similar(a::LocalVector, T=eltype(a), dims=size(a)) = similar(a.a, T, dims)

function (==)(x::LocalVector, y::AbstractArray)
  return x.a == y
end
