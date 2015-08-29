# AbstractVector wrapper around PETSc Vec types
export Vec
typealias pVec Ptr{Void} # Vec arguments in Petsc are pointers
type Vec <: AbstractVector{PetscScalar}
    p::pVec
    assembling::Bool # whether we are in the middle of assemble(vec) do ..
    insertmode::InsertMode # current mode for setindex!
    function Vec(p::pVec)
        v = new(p, false, DEFAULT_INSERTMODE)
        finalizer(v, v -> ccall((:VecDestroy,petsc), PetscErrorCode,
                                (Ptr{pVec},), &v))
        v
    end
end

function Vec(; comm::Ptr{Void}=PETSC_COMM_WORLD)
    p = Array(pVec, 1)
    chk(ccall((:VecCreate,petsc), PetscErrorCode, (Ptr{Void}, Ptr{pVec}),
              PETSC_COMM_WORLD, p))
    Vec(p[1])
end

function settype!(a::Vec, T::String)
    chk(ccall((:VecSetType,petsc), PetscErrorCode, (pVec, Ptr{Uint8}),
              a, bytestring(T)))
    a
end
function gettype(a::Vec)
    p = Array(Ptr{Uint8},1)
    chk(ccall((:VecGetType,petsc), PetscErrorCode, (pVec, Ptr{Ptr{Uint8}}),
              a, p))
    bytestring(p[1])
end

function setsizes!(x::Vec, m::Integer; mlocal::Integer=PETSC_DECIDE)
    chk(ccall((:VecSetSizes,petsc), PetscErrorCode,
              (pVec, PetscInt, PetscInt), x, mlocal, m))
    x
end

Vec(m::Integer;
    mlocal::Integer=PETSC_DECIDE,
    comm::Ptr{Void}=PETSC_COMM_WORLD,
    T::String="mpi") =
  setsizes!(settype!(Vec(comm=comm), T), m, mlocal=mlocal)

##########################################################################
import Base: convert, length, size, similar, copy
export lengthlocal, sizelocal

convert(::Type{pVec}, v::Vec) = v.p

function length(x::Vec)
    sz = Array(PetscInt, 1)
    chk(ccall((:VecGetSize, petsc), PetscErrorCode, (pVec, Ptr{PetscInt}),
              x, sz))
    sz[1]
end
size(x::Vec) = (length(x),)

function lengthlocal(x::Vec)
    sz = Array(PetscInt, 1)
    chk(ccall((:VecGetLocalSize, petsc), 
              PetscErrorCode, (pVec, Ptr{PetscInt}), x, sz))
    sz[1]
end
sizelocal(x::Vec) = (lengthlocal(x),)
sizelocal{T,n}(t::AbstractArray{T,n}, d) = (d>n ? 1 : sizelocal(t)[d])

function similar(x::Vec)
    p = Array(pVec, 1)
    chk(ccall((:VecDuplicate,petsc), PetscErrorCode,(pVec,Ptr{pVec}),
              x, p))
    Vec(p[1])
end
similar(x::Vec, ::Type{PetscScalar}) = similar(x)
similar(x::Vec, ::Type{PetscScalar}, len::Int) =
    len==length(x) ? similar(x) : Vec(len; comm=comm(x), T=gettype(x))

function copy(x::Vec)
    y = similar(x)
    chk(ccall((:VecCopy,petsc), PetscErrorCode, (pVec, pVec),
              x, y))
    y
end

##########################################################################
import Base: setindex!, fill!, to_index
export assemble, isassembled

# for efficient vector assembly, put all calls to x[...] = ... inside
# assemble(x) do ... end
AssemblyBegin(x::Vec,t::MatAssemblyType=MAT_FINAL_ASSEMBLY) = chk(ccall((:VecAssemblyBegin,petsc), PetscErrorCode, (pVec,), x))
AssemblyEnd(x::Vec,t::MatAssemblyType=MAT_FINAL_ASSEMBLY) = chk(ccall((:VecAssemblyEnd,petsc), PetscErrorCode, (pVec,), x))
isassembled(x::Vec) = !x.assembling
# assemble(f::Function, x::Vec) is defined in mat.jl

# like x[i] = v, but requires i to be 0-based indices for Petsc
function setindex0!(x::Vec, v::Array{PetscScalar}, i::Array{PetscInt})
    n = length(v)
    if n != length(i)
        throw(ArgumentError("length(values) != length(indices)"))
    end
    chk(ccall((:VecSetValues,petsc), PetscErrorCode,
              (pVec, PetscInt,Ptr{PetscInt},Ptr{PetscScalar}, InsertMode),
              x, n, i, v, insertmode(x)))
    if !x.assembling
        AssemblyBegin(x)
        AssemblyEnd(x)
    end
    x
end

function setindex!(x::Vec, v::Number, i::Real)
    # can't call VecSetValue since that is a static inline function
    setindex0!(x, PetscScalar[ v ], PetscInt[ to_index(i) - 1 ])
    v
end
function setindex!{T<:Real}(x::Vec, v::Array{PetscScalar}, I::AbstractArray{T})
    I0 = PetscInt[ to_index(i)-1 for i in I ]
    setindex0!(x, v, I0)
end

setindex!{T<:Real}(x::Vec, v::Number, I::AbstractArray{T}) = assemble(x) do
    for i in I
        x[i] = v
    end
    x
end

function fill!(x::Vec, v::Number)
    chk(ccall((:VecSet,petsc), PetscErrorCode, (pVec, PetscScalar), x, v))
    return x
end

function setindex!(x::Vec, v::Number, I::Ranges{Int})
    if abs(step(I)) == 1 && min(I) == 1 && max(I) == length(x)
        fill!(x, v)
        return v
    else
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
setindex!(A::Vec, x::Number, I::AbstractArray{Bool}) = assemble(X) do
    for i = 1:length(I)
        if I[i]
            A[i] = x
        end
    end
    A
end
for T in (:(Array{PetscScalar}),:AbstractArray) # avoid method ambiguities
    @eval setindex!(A::Vec, X::$T, I::AbstractArray{Bool}) = assemble(X) do
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
function getindex0(x::Vec, i::Vector{PetscInt})
    v = similar(i, PetscScalar)
    chk(ccall((:VecGetValues,petsc), PetscErrorCode,
              (pVec, PetscInt, Ptr{PetscInt}, Ptr{PetscScalar}),
              x, length(i), i, v))
    v
end

getindex(x::Vec, i::Real) = getindex0(x, PetscInt[to_index(i)-1])[1]

getindex{T<:Real}(x::Vec, I::AbstractVector{T}) =
  getindex0(x, PetscInt[ to_index(i)-1 for i in I ])

##########################################################################
import Base: abs, exp, log, conj, conj!, max, min, findmax, findmin, norm, max, min, .*, ./, +, -, scale!, dot, ==, !=, sum
export normalize!, abs!, exp!, log!, chop!

for (f,pf) in ((:abs,:VecAbs), (:exp,:VecExp), (:log,:VecLog), 
               (:conj,:VecConjugate))
    fb = symbol(string(f, "!"))
    @eval begin
        function $fb(x::Vec)
            chk(ccall(($(Expr(:quote, pf)),petsc), PetscErrorCode, 
                      (pVec,), x))
            x
        end
        $f(x::Vec) = $fb(copy(x))
    end
end

function chop!(x::Vec, tol::Real)
    chk(ccall((:VecChop, petsc), PetscErrorCode, (pVec, PetscReal), x, tol))
    x
end

for (f,pf) in ((:max, :VecMax), (:min, :VecMin))
    ff = symbol(string("find", f))
    @eval begin
        function $ff(x::Vec)
            i = Array(PetscInt, 1)
            v = Array(PetscReal, 1)
            chk(ccall(($(Expr(:quote, pf)),petsc), PetscErrorCode,
                      (pVec, Ptr{PetscInt}, Ptr{PetscReal}),
                      x, i, v))
            (v[1], i[1]+1)
        end
        $f(x::Vec) = $ff(x)[1]
    end
end

function norm(x::Vec, p::Number)
    v = Array(PetscReal, 1)
    n = p == 1 ? NORM_1 : p == 2 ? NORM_2 : p == Inf ? NORM_INFINITY :
       throw(ArgumentError("unrecognized Petsc norm $p"))
    chk(ccall((:VecNorm,petsc), PetscErrorCode, 
              (pVec, NormType, Ptr{PetscReal}), x, n, v))
    v[1]
end

# computes v = norm(x,2), divides x by v, and returns v
function normalize!(x::Vec)
    v = Array(PetscReal, 1)
    chk(ccall((:VecNormalize,petsc), PetscErrorCode,
              (pVec, Ptr{PetscReal}), x, v))
    v[1]
end

function dot(x::Vec, y::Vec)
    d = Array(PetscScalar, 1)
    chk(ccall((:VecDot,petsc), PetscErrorCode, 
              (pVec,pVec,Ptr{PetscScalar}), x, y, d))
    d[1]
end

# TODO: also add dotu for unconjugated dot product?  should go in Base first

# pointwise operations on pairs of vectors (TODO: support in-place variants?)
for (f,pf) in ((:max,:VecPointwiseMax), (:min,:VecPointwiseMin),
               (:.*,:VecPointwiseMult), (:./,:VecPointwiseDivide))
    @eval function ($f)(x::Vec, y::Vec)
        w = similar(x)
        chk(ccall(($(Expr(:quote, pf)),petsc), PetscErrorCode,
                  (pVec, pVec, pVec), w, x, y))
        w
    end
end

function scale!(x::Vec, s::Number)
    chk(ccall((:VecScale,petsc), PetscErrorCode, (pVec,PetscScalar), x,s))
    x
end

function (+)(x::Vec, a::Number...)
    chk(ccall((:VecShift,petsc), PetscErrorCode, (pVec,PetscScalar), x,
              sum(a)))
    x
end
(+)(a::Number, x::Vec) = x + a
(-)(x::Vec, a::Number) = x + (-a)
(-)(x::Vec) = scale(x, -1)
(-)(a::Number, x::Vec) = (-x) + a
(.*)(x::Vec, a::Number...) = scale(x, prod(a))
(.*)(a::Number, x::Vec) = scale(x, a)
(./)(x::Vec, a::Number) = scale(x, 1/a)
function (./)(a::Number, x::Vec)
    chk(ccall((:VecReciprocal,petsc), PetscErrorCode, (pVec,), x))
    if a != 1.0
        scale!(x, a)
    end
    x
end

function (==)(x::Vec, y::Vec)
    b = Array(PetscBool,1)
    chk(ccall((:VecEqual,petsc), PetscErrorCode, (pVec,pVec,Ptr{PetscBool}), 
              x, y, b))
    bool(b[1])
end
(!=)(x::Vec, y::Vec) = !(x == y)

function sum(x::Vec)
    s = Array(PetscScalar,1)
    chk(ccall((:VecSum,petsc), PetscErrorCode, (pVec,Ptr{PetscScalar}), x,s))
    s[1]
end

##########################################################################
export axpy!, aypx!, axpby!, axpbypcz!
import Base.LinAlg.BLAS.axpy!

# y <- alpha*x + y
function axpy!(alpha::Number, x::Vec, y::Vec)
    chk(ccall((:VecAXPY,petsc), PetscErrorCode, (pVec,PetscScalar,pVec),
              y, alpha, x))
    y
end
# w <- alpha*x + y
function axpy!(alpha::Number, x::Vec, y::Vec, w::Vec)
    chk(ccall((:VecWAXPY,petsc), PetscErrorCode, (pVec,PetscScalar,pVec),
              w, alpha, x, y))
    y
end
# y <- alpha*y + x
function aypx!(x::Vec, alpha::Number, y::Vec)
    chk(ccall((:VecAYPX,petsc), PetscErrorCode, (pVec,PetscScalar,pVec),
              y, alpha, x))
    y
end
# y <- alpha*x + beta*y
function axpby!(alpha::Number, x::Vec, beta::Number, y::Vec)
    chk(ccall((:VecAXPBY,petsc), PetscErrorCode,
              (pVec,PetscScalar,PetscScalar,pVec), y, alpha, beta, x))
    y
end
# z <- alpha*x + beta*y + gamma*z
function axpbypcz!(alpha::Number, x::Vec, beta::Number, y::Vec,
                   gamma::Number, z::Vec)
    chk(ccall((:VecAXPBYPCZ,petsc), PetscErrorCode,
              (pVec,PetscScalar,PetscScalar,PetscScalar,pVec,pVec),
              z, alpha, beta, gamma, x, y))
    z
end
# y <- \sum_i alphax[i][1]*alphax[i][2] + y
function axpy!{T<:Number}(alphax::AbstractVector{(T,Vec)}, y::Vec)
    n = length(alphax)
    alpha = Array(PetscScalar, n)
    x = Array(pVec, n)
    for i = 1:n
        ax = alphax[i]
        alpha[i] = ax[1]
        x[i] = ax[2].p
    end
    chk(ccall((:VecMAXPY,petsc), PetscErrorCode, 
              (pVec,PetscInt,Ptr{PetscScalar},Ptr{pVec}), y,n,alpha,x))
    y
end

##########################################################################
