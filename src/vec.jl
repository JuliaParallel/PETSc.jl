# AbstractVector wrapper around PETSc Vec types
export Vec
#typealias pVec Ptr{Void} # Vec arguments in Petsc are pointers
# T is data type
# VType is type of vector from C.VecType
type Vec{T, VType} <: AbstractVector{T}
    p::C.Vec{T}
    assembling::Bool # whether we are in the middle of assemble(vec) do ..
    insertmode::C.InsertMode # current mode for setindex!
    comm::MPI.Comm  # communicator this process is defined on
    vectype::C.VecType
    function Vec(p::C.Vec{T}; comm=MPI.COMM_SELF)  # default sequantial vector
        v = new(p, false, C.INSERT_VALUES, comm, VType)
        settype!(v, VType)  # set the type here to ensure it matches VType
#        finalizer(v, VecDestroy)
        return v
    end
end

function Vec{T}(::Type{T}, vtype::C.VecType=C.VECSEQ; comm=MPI.COMM_SELF)
    p = Array(C.Vec{T}, 1)
    C.VecCreate(comm, p)
    Vec{T, vtype}(p[1])
end

function Vec{T <: Scalar}(::Type{T}, len::Integer, vtype::C.VecType=C.VECSEQ ; comm=MPI.COMM_SELF)
  p = Array(C.Vec{T}, 1)
  C.VecCreate(comm, p)
  vec = Vec{T, vtype}(p[1])
#  settype!(vec, VType)
  setsizes!(vec, len)
  vec
end


#=
Vec(m::Integer;
    mlocal::Integer=C.PETSC_DECIDE,
    comm=PETSC_COMM_WORLD,
    T::C.VecType`="mpi") =
  setsizes!(settype!(Vec(comm=comm), T), m, mlocal=mlocal)
=#


function VecDestroy(vec::Vec)

  tmp = Array(PetscBool, 1)
  C.PetscFinalized(eltype(vec), tmp)
   
  if tmp[1] == 0  # if petsc has not been finalized yet
    C.VecDestroy([vec.p])
  end

   # if Petsc has been finalized, let the OS deallocate the memory
end

 

function settype!(a::Vec, T::C.VecType)
    chk(C.VecSetType(a.p, T))
    a
end

function gettype(a::Vec)
    p = Array(C.VecType, 1)  # was Uint8
#    p[1] = "a"  # need to initialize symbol array to something
    
    chk(C.VecGetType(a.p, p))
    p[1]
end

function setsizes!(x::Vec, m::Integer; mlocal::Integer=C.PETSC_DECIDE)
    chk(C.VecSetSizes(x.p, m, mlocal))
    x
end
##########################################################################
import Base: convert, length, size, similar, copy
export lengthlocal, sizelocal

convert(::Type{C.Vec}, v::Vec) = v.p

function length(x::Vec)
    sz = Array(PetscInt, 1)
    chk(C.VecGetSize(x.p, sz))
    Int(sz[1])
end
size(x::Vec) = (length(x),)

function lengthlocal(x::Vec)
    sz = Array(PetscInt, 1)
    chk(C.VecGetLocalSize(x.p, sz))
    sz[1]
end
sizelocal(x::Vec) = (lengthlocal(x),)
sizelocal{T,n}(t::AbstractArray{T,n}, d) = (d>n ? 1 : sizelocal(t)[d])

function similar{T, VType}(x::Vec{T, VType})
    p = Array(C.Vec{T}, 1)
    chk(C.VecDuplicate( x.p, p))
    Vec{T, VType}(p[1])
end


similar{T}(x::Vec{T}, ::Type{T}) = similar(x)
function similar{T, VType}(x::Vec{T, VType}, ::Type{T}, len::Int)
#    println("T = ", T)
#    println("VType = ", VType)
    len==length(x) ? similar(x) : Vec(T, len, VType; comm=x.comm)
end

similar{T, DType <: Integer}(x::Vec{T}, t::Type{T}, dims::Tuple{DType}) = similar(x, t, dims[1])

function copy(x::Vec)
    y = similar(x)
    chk(C.VecCopy(x.p, y.p))
    y
end

##########################################################################
import Base: setindex!, fill!, to_index
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

    chk(C.VecSetValues(x.p, n, i, v, x.insertmode))
    if !x.assembling
        AssemblyBegin(x)
        AssemblyEnd(x)
    end
    x
end

function setindex!{T}(x::Vec{T}, v::Number, i::Real)
    # can't call VecSetValue since that is a static inline function
    setindex0!(x, T[ v ], PetscInt[ to_index(i) - 1 ])
    v
end

# what is this function doing?
#=
function setindex!{T<: Scalar}(x::Vec{T}, v::Array{T}, I::AbstractArray{T})
    I0 = PetscInt[ to_index(i)-1 for i in I ]
    setindex0!(x, v, I0)
end
=#


# set multiple entries to a single value
setindex!{T<:Real}(x::Vec, v::Number, I::AbstractArray{T}) = assemble(x) do
    for i in I
        x[i] = v
    end
    x
end

function fill!{T}(x::Vec{T}, v::Number)
    chk(C.VecSet(x.p, T(v)))
    return x
end

function setindex!(x::Vec, v::Number, I::Range{Int})
    if abs(step(I)) == 1 && minimum(I) == 1 && maximum(I) == length(x)
        fill!(x, v)
        return v
    else
        # why use invoke here?
        return invoke(setindex!, (Vec,typeof(v),AbstractVector{Int}), x,v,I)
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
    @eval setindex!{T2 <: Scalar}(A::Vec, X::$T, I::AbstractArray{Bool}) = assemble(A) do
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

getindex(x::Vec, i::Real) = getindex0(x, PetscInt[to_index(i)-1])[1]

getindex{T<:Real}(x::Vec, I::AbstractVector{T}) =
  getindex0(x, PetscInt[ to_index(i)-1 for i in I ])

##########################################################################
import Base: abs, exp, log, conj, conj!, max, min, findmax, findmin, norm, maximum, minimum, .*, ./, +, -, scale!, dot, ==, !=, sum
export normalize!, abs!, exp!, log!, chop!

for (f,pf) in ((:abs,:VecAbs), (:exp,:VecExp), (:log,:VecLog), 
               (:conj,:VecConjugate))
    fb = symbol(string(f, "!"))
    @eval begin
        function $fb(x::Vec)
            chk(C.$pf(x.p))
#            chk(ccall(($(Expr(:quote, pf)),petsc), PetscErrorCode, 
 #                     (pVec,), x))
            x
        end
 
        $f(x::Vec) = $fb(copy(x))
    end
end

# skip for now
#=
function chop!(x::Vec, tol::Real)
    VecChop(x.c, tol)
    chk(ccall((:VecChop, petsc), PetscErrorCode, (pVec, PetscReal), x, tol))
    x
end
=#


# real versions
for (f,pf, sf) in ((:findmax, :VecMax, :maximum), (:findmin, :VecMin, :minimum))
#    ff = symbol(string("find", f))
    @eval begin
        function $f{T <: Scalar}(x::Vec{T})
            i = Array(PetscInt, 1)
            v = Array(T, 1)
            chk(C.$pf(x.p, i, v))
            (v[1], i[1]+1)
        end
        $sf(x::Vec) = $f(x)[1]
    end
end

# complex versions
for (f,pf) in ((:maximum, :VecMax), (:minimum, :VecMin))
    ff = symbol(string("find", f))
    @eval begin
        function $ff{T}(x::Vec{Complex{T}})
            i = Array(PetscInt, 1)
            v = Array(T, 1)
            chk(C.$pf(x.p, i, v))
            (v[1], i[1]+1)
        end
        $f(x::Vec) = $ff(x)[1]
    end
end



function norm{T <: Real}(x::Vec{T}, p::Number)
    v = Array(T, 1)
    n = p == 1 ? C.NORM_1 : p == 2 ? C.NORM_2 : p == Inf ? C.NORM_INFINITY :
       throw(ArgumentError("unrecognized Petsc norm $p"))
    chk(C.VecNorm(x.p, n, p))
    v[1]
end

function norm{T <: Real}(x::Vec{Complex{T}}, p::Number)
    v = Array(PetscReal, 1)
    n = p == 1 ? C.NORM_1 : p == 2 ? C.NORM_2 : p == Inf ? C.NORM_INFINITY :
       throw(ArgumentError("unrecognized Petsc norm $p"))
    chk(C.VecNorm(x.p, n, p))
    v[1]
end




# computes v = norm(x,2), divides x by v, and returns v
function normalize!{T <: Real}(x::Vec{T})
    v = Array(T, 1)
    chk(C.VecNormalize(x.p, v))
    v[1]
end

function normalize!{T <: Real}(x::Vec{Complex{T}})
    v = Array(T, 1)
    chk(C.VecNormalize(x.p, v))
    v[1]
end




function dot{T}(x::Vec{T}, y::Vec{T})
    d = Array(T, 1)
    chk(C.VecDot(x.p, y.p, d))
    d[1]
end

# TODO: also add dotu for unconjugated dot product?  should go in Base first

# pointwise operations on pairs of vectors (TODO: support in-place variants?)
for (f,pf) in ((:max,:VecPointwiseMax), (:min,:VecPointwiseMin),
               (:.*,:VecPointwiseMult), (:./,:VecPointwiseDivide))
    @eval function ($f)(x::Vec, y::Vec)
        w = similar(x)
        chk(C.$pf(w.p, x.p, y.p))
        w
    end
end

function scale!{T}(x::Vec{T}, s::Number)
    chk(C.VecScale(x.p, T(s)))
    x
end

function (+){T <: Scalar}(x::Vec{T}, a::Number...)
    C.VecShift(x.p, T(sum(a)))
    x
end


(+){T <: Scalar}(a::Number, x::Vec{T}) = x + a
(-){T <: Scalar}(x::Vec{T}, a::Number) = x + (-a)
(-)(x::Vec) = scale(x, -1)
(-){T <: Scalar}(a::Number, x::Vec{T}) = (-x) + a
(.*){T <: Scalar}(x::Vec{T}, a::Number...) = scale(x, prod(a))
(.*){T <: Scalar}(a::Number, x::Vec{T}) = scale(x, a)
(./){T <: Scalar}(x::Vec{T}, a::Number) = scale(x, 1/a)


function (./)(a::Number, x::Vec)

    chk(C.VecReciprocal(x.p))
    if a != 1.0
        scale!(x, a)
    end
    x
end

function (==)(x::Vec, y::Vec)
    b = Array(PetscBool,1)
    chk(C.VecEqual(x.p, y.b, b))
    b[1] != 0
end
(!=)(x::Vec, y::Vec) = !(x == y)

function sum{T}(x::Vec{T})
    s = Array(T,1)
    chk(C.VecSum(x.p, s))
    s[1]
end

##########################################################################
export axpy!, aypx!, axpby!, axpbypcz!
import Base.LinAlg.BLAS.axpy!

# y <- alpha*x + y
function axpy!{T}(alpha::Number, x::Vec{T}, y::Vec)
    chk(C.VecAXPY(y.p, T(alpha), x.p)) 
    y
end
# w <- alpha*x + y
function axpy!{T}(alpha::Number, x::Vec{T}, y::Vec, w::Vec)
    chk(C.VecWAXPY(w.p, T(alpha), x.p, y.p))
    y
end
# y <- alpha*y + x
function aypx!{T}(x::Vec{T}, alpha::Number, y::Vec)
    chk(C.VecAYPX( y.p, T(alpha), x.p))
    y
end
# y <- alpha*x + beta*y
function axpby!{T}(alpha::Number, x::Vec{T}, beta::Number, y::Vec)
    chk(C.VecAXPBY(y.p, T(alpha), T(beta), x.p))
    y
end
# z <- alpha*x + beta*y + gamma*z
function axpbypcz!{T}(alpha::Number, x::Vec{T}, beta::Number, y::Vec,
                   gamma::Number, z::Vec)
    chk(C.VecAXPBYPCZ(z.p, T(alpha), T(beta), T(gamma), x.p, y.p))
    z
end

#=
# y <- \sum_i alphax[i][1]*alphax[i][2] + y
function axpy!{T<:Number}(alphax::AbstractVector{(T,Vec)}, y::Vec)
    n = length(alphax)
    alpha = Array(T, n)
    x = Array(Vec, n)
    for i = 1:n
        ax = alphax[i]
        alpha[i] = ax[1]
        x[i] = ax[2]
    end

    maxpy!(y, alpha, x)
end
=#
function maxpy!{T <: Number, VType}(y::Vec, alpha::AbstractArray, x::Array{Vec{T, VType}})

  @assert length(alpha) == length(x)
  n = length(x)
  _x = Array(C.Vec{T}, n)
  _alpha = Array(T, n)
  for i=1:n
    _x[i] = x[i].p
    _alpha[i] = alpha[i]
  end

  C.VecMAXPY(y.p, n, alpha, _x)
  y
end

##########################################################################
# elementwise operations

function (.*)(x::Vec, y::Vec)
  z = similar(x)
#  println("before mult, z = ", z)
  chk(C.VecPointwiseMult(z.p, x.p, y.p))
  return z
end


function (./)(x::Vec, y::Vec)
  z = similar(x)
  chk(C.VecPointwiseDivide(z.p, x.p, y.p))
  return z
end



function (.^){T}(x::Vec{T}, y::Number)
  z = copy(x)
  chk(C.VecPow(z.p, T(y)))
  return z
end

function (+){T}(x::Vec{T}, y::Vec{T})
  z = similar(x)
  chk(C.VecWAXPY(z.p, T(1), x.p, y.p))
  return z
end

function (-){T}(x::Vec{T}, y::Vec{T})
  z = similar(x)
  chk(C.VecWAXPY(z.p, T(-1), y.p, x.p))
  return z
end




