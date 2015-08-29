# AbstractMatrix wrapper around Petsc Mat
export Mat
typealias pMat Ptr{Void} # Mat arguments in Petsc are pointers
type Mat <: AbstractSparseMatrix{PetscScalar,PetscInt}
    p::pMat
    assembling::Bool # whether we are in the middle of assemble(vec) do ..
    insertmode::InsertMode # current mode for setindex!
    data::Any # keep a reference to anything needed for the Mat
    function Mat(p::pMat, data=nothing)
        A = new(p, false, INSERT_VALUES, data)
        finalizer(A, A -> ccall((:MatDestroy,petsc), PetscErrorCode,
                                (Ptr{pMat},), &A))
        A
    end
end

function setoption!(m::Mat, option::Integer, val::Bool)
    chk(ccall((:MatSetOption,petsc), PetscErrorCode, 
              (pMat, MatOption, PetscBool), m, option, val))
    m
end

function Mat(; comm::Ptr{Void}=PETSC_COMM_WORLD)
    p = Array(pMat, 1)
    chk(ccall((:MatCreate,petsc), PetscErrorCode, (Ptr{Void}, Ptr{pMat}),
              comm, p))
    Mat(p[1])
end

function settype!(a::Mat, T::String)
    chk(ccall((:MatSetType,petsc), PetscErrorCode, (pMat, Ptr{Uint8}),
              a, bytestring(T)))
    a
end
function gettype(a::Mat)
    p = Array(Ptr{Uint8},1)
    chk(ccall((:MatGetType,petsc), PetscErrorCode, (pMat, Ptr{Ptr{Uint8}}),
              a, p))
    bytestring(p[1])
end

function setsizes!(a::Mat, m::Integer, n::Integer;
                   mlocal::Integer=PETSC_DECIDE, nlocal::Integer=PETSC_DECIDE)
    chk(ccall((:MatSetSizes,petsc), PetscErrorCode,
              (pMat, PetscInt, PetscInt, PetscInt, PetscInt),
              a, mlocal, nlocal, m, n))
    a
end

function setpreallocation!(a::Mat, T::String=gettype(a);
                           nz::Integer=16, nnz::AbstractVector=PetscInt[],
                           onz::Integer=0, onnz::AbstractVector=PetscInt[])
    if T == "seqaij"
        pnnz = if isempty(nnz)
            C_NULL
        else
            if length(nnz) != size(a,1)
                throw(ArgumentError("length(nnz) must be # rows"))
            end
            isa(nnz,Vector{PetscInt}) ? nnz : PetscInt[ i for i in nnz ]
        end

        chk(ccall((:MatSeqAIJSetPreallocation,petsc), PetscErrorCode,
                  (pMat, PetscInt, Ptr{PetscInt}), a, nz, pnnz))
    elseif T == "mpiaij"
        mlocal = sizelocal(a,1)
        pnnz = if isempty(nnz)
            C_NULL
        else
            if length(nnz) != mlocal
                throw(ArgumentError("length(nnz) must be # local rows"))
            end
            isa(nnz,Vector{PetscInt}) ? nnz : PetscInt[ i for i in nnz ]
        end
        ponnz = if isempty(onnz)
            C_NULL
        else
            if length(onnz) != mlocal
                throw(ArgumentError("length(onnz) must be # local rows"))
            end
            isa(onnz,Vector{PetscInt}) ? onnz : PetscInt[ i for i in onnz ]
        end

        chk(ccall((:MatMPIAIJSetPreallocation,petsc), PetscErrorCode,
                  (pMat, PetscInt, Ptr{PetscInt}, PetscInt, Ptr{PetscInt}), 
                  a, nz, pnnz, onz, ponnz))
    else # TODO
        throw(ArgumentError("unsupported matrix type $T"))
    end
    a
end

Mat(m::Integer, n::Integer;
    mlocal::Integer=PETSC_DECIDE, nlocal::Integer=PETSC_DECIDE,
    nz::Integer=16, nnz::AbstractVector=PetscInt[],
    onz::Integer=0, onnz::AbstractVector=PetscInt[],
    comm::Ptr{Void}=PETSC_COMM_WORLD,
    T::String="mpiaij") =
  setoption!(setpreallocation!(setsizes!(settype!(Mat(comm=comm), T),
                                         m, n, mlocal=mlocal, nlocal=nlocal),
                               T, nz=nz, nnz=nnz, onz=onz, onnz=onnz),
             MAT_ROW_ORIENTED, false) # julia data is column-major

# construct Vec for multiplication by a::Mat or transpose(a::Mat)
const mat2vec = [ "mpiaij" => "aij", "seqaij" => "seq" ]
Vec(a::Mat, transposed=false) = 
  transposed ? Vec(size(a,1), comm=comm(a), T=mat2vec[gettype(a)],
                   mlocal=sizelocal(a,1)) :
               Vec(size(a,2), comm=comm(a), T=mat2vec[gettype(a)],
                   mlocal=sizelocal(a,2))

#############################################################################
import Base: nnz

convert(::Type{pMat}, a::Mat) = a.p

function size(a::Mat)
    mn = Array(PetscInt, 2)
    chk(ccall((:MatGetSize,petsc), PetscErrorCode,
              (pMat,Ptr{PetscInt},Ptr{PetscInt}),
              a, mn, pointer(mn)+sizeof(PetscInt)))
    (mn[1], mn[2])
end

function sizelocal(a::Mat)
    mn = Array(PetscInt, 2)
    chk(ccall((:MatGetLocalSize,petsc), PetscErrorCode,
              (pMat,Ptr{PetscInt},Ptr{PetscInt}),
              a, mn, pointer(mn)+sizeof(PetscInt)))
    (mn[1], mn[2])
end
lengthlocal(a::Mat) = prod(sizelocal(a))

function similar(a::Mat)
    p = Array(pMat, 1)
    chk(ccall((:MatDuplicate,petsc), PetscErrorCode,
              (pMat, MatDuplicateOption, Ptr{pMat}),
              a, MAT_DO_NOT_COPY_VALUES, p))
    Mat(p[1])
end
similar(a::Mat, ::Type{PetscScalar}) = similar(a)
similar(a::Mat, ::Type{PetscScalar}, ::Type{PetscInt}) = similar(a)
similar(a::Mat, ::Type{PetscScalar}, m::Int, n::Int) = 
    (m,n) == size(a) ? similar(a) : Mat(m,n)
similar(a::Mat, ::Type{PetscScalar}, d::(Int,Int)) = 
    similar(a, PetscScalar, d...)

function copy(a::Mat)
    p = Array(pMat, 1)
    chk(ccall((:MatDuplicate,petsc), PetscErrorCode,
              (pMat, MatDuplicateOption, Ptr{pMat}),
              a, MAT_COPY_VALUES, p))
    Mat(p[1])
end

function getinfo(m::Mat, infotype::Integer=MAT_GLOBAL_SUM)
    info = Array(MatInfo,1)
    chk(ccall((:MatGetInfo,petsc), PetscErrorCode,
              (pMat,MatInfoType,Ptr{MatInfo}), m, infotype, info))
    info[1]
end

nnz(m::Mat) = int(getinfo(m).nz_used)

#############################################################################

# for efficient matrix assembly, put all calls to A[...] = ... inside
# assemble(A) do ... end
AssemblyBegin(x::Mat,t::MatAssemblyType=MAT_FINAL_ASSEMBLY) = chk(ccall((:MatAssemblyBegin,petsc), PetscErrorCode, (pMat,MatAssemblyType), x,t))
AssemblyEnd(x::Mat,t::MatAssemblyType=MAT_FINAL_ASSEMBLY) = chk(ccall((:MatAssemblyEnd,petsc), PetscErrorCode, (pMat,MatAssemblyType), x,t))

function isassembled(x::Mat)
    if x.assembling
        return false
    else
        b = Array(PetscBool, 1)
        chk(ccall((:MatAssembled,petsc), PetscErrorCode,
                  (pMat, Ptr{PetscBool}), x, b))
        return bool(b[1])
    end
end

function assemble(f::Function, x::Union(Vec,Mat), 
                  insertmode=DEFAULT_INSERTMODE,
                  assemblytype::MatAssemblyType=MAT_FINAL_ASSEMBLY)
    old_assembling = x.assembling
    old_insertmode = x.insertmode
    try
        x.assembling = true
        if insertmode != DEFAULT_INSERTMODE
            if x.insertmode != DEFAULT_INSERTMODE && x.insertmode != insertmode
                error("nested assemble with different insertmodes not allowed")
            end
            x.insertmode = insertmode
        end
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

# finalize matrix assembly
assemble(x::Union(Vec,Mat)) = assemble(() -> nothing, x)

# intermediate assembly, before it is finally compressed for use
iassemble(f::Function, x::Mat) =
    assemble(f, x, DEFAULT_INSERTMODE, MAT_FLUSH_ASSEMBLY)

# like x[i,j] = v, but requires i,j to be 0-based indices for Petsc
function setindex0!(x::Mat, v::Array{PetscScalar}, 
                    i::Array{PetscInt}, j::Array{PetscInt})
    ni = length(i)
    nj = length(j)
    if length(v) != ni*nj
        throw(ArgumentError("length(values) != length(indices)"))
    end
    chk(ccall((:MatSetValues,petsc), PetscErrorCode,
              (pMat, PetscInt,Ptr{PetscInt}, PetscInt,Ptr{PetscInt},
               Ptr{PetscScalar}, InsertMode),
              x, ni,i, nj,j, v, insertmode(x)))
    if !x.assembling
        AssemblyBegin(x,MAT_FLUSH_ASSEMBLY)
        AssemblyEnd(x,MAT_FLUSH_ASSEMBLY)
    end
    x
end

function setindex!(x::Mat, v::Number, i::Real, j::Real)
    # can't call MatSetValue since that is a static inline function
    setindex0!(x, PetscScalar[ v ], 
               PetscInt[ to_index(i) - 1 ], 
               PetscInt[ to_index(j) - 1 ])
    v
end
function setindex!{T1<:Real, T2<:Real}(x::Mat, v::Array{PetscScalar}, 
                                       I::AbstractArray{T1},
                                       J::AbstractArray{T2})
    I0 = PetscInt[ to_index(i)-1 for i in I ]
    J0 = PetscInt[ to_index(j)-1 for j in J ]
    setindex0!(x, v, I0, J0)
end
function setindex!{T<:Real}(x::Mat, v::Array{PetscScalar}, 
                            i::Real, J::AbstractArray{T})
    I0 = PetscInt[ to_index(i)-1 ]
    J0 = PetscInt[ to_index(j)-1 for j in J ]
    setindex0!(x, v, I0, J0)
end
function setindex!{T<:Real}(x::Mat, v::Array{PetscScalar}, 
                            I::AbstractArray{T}, j::Real)
    I0 = PetscInt[ to_index(i)-1 for i in I ]
    J0 = PetscInt[ to_index(j)-1 ]
    setindex0!(x, v, I0, J0)
end

setindex!{T1<:Real, T2<:Real}(x::Mat, v::Number, 
                              I::AbstractArray{T1},
                              J::AbstractArray{T2}) = iassemble(x) do
    for i in I
        for j in J
            x[i,j] = v
        end
    end
    x
end
setindex!{T<:Real}(x::Mat, v::Number, 
                   i::Real, J::AbstractArray{T}) = iassemble(x) do
    for j in J
        x[i,j] = v
    end
    x
end
setindex!{T<:Real}(x::Mat, v::Number, 
                   I::AbstractArray{T}, j::Real) = iassemble(x) do
    for i in I
        x[i,j] = v
    end
    x
end

function setindex!{T0<:Real, T1<:Real, T2<:Real}(x::Mat, v::AbstractArray{T0}, 
                                                 I::AbstractArray{T1},
                                                 J::AbstractArray{T2})
    if length(v) != length(I)*length(J)
        throw(ArgumentError("length(values) != length(indices)"))
    end
    v0 = PetscScalar[ z for z in v ]
    setindex!(x, v0, I, J)
end

function fill!(x::Mat, v::Number)
    if v == 0
        chk(ccall((:MatZeroEntries,petsc), PetscErrorCode, (pMat,), x))
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
    x
end

#############################################################################

# like getindex but for 0-based indices i and j
function getindex0(x::Mat, i::Vector{PetscInt}, j::Vector{PetscInt})
    ni = length(i)
    nj = length(j)
    v = Array(PetscScalar, nj, ni) # row-major!
    chk(ccall((:MatGetValues,petsc), PetscErrorCode,
              (pMat, PetscInt, Ptr{PetscInt}, PetscInt, Ptr{PetscInt},
               Ptr{PetscScalar}), x, ni, i, nj, j, v))
    ni <= 1 || nj <= 1 ? reshape(v, ni, nj) : transpose(v)
end

getindex(a::Mat, i0::Real, i1::Real) =
  getindex0(a, PetscInt[ to_index(i0)-1 ], PetscInt[ to_index(i1)-1 ])[1]

getindex{T0<:Real,T1<:Real}(a::Mat, 
                            I0::AbstractArray{T0},
                            I1::AbstractArray{T1}) =
  getindex0(a, PetscInt[ to_index(i0)-1 for i0 in I0 ],
            PetscInt[ to_index(i1)-1 for i1 in I1 ])

getindex{T0<:Real}(a::Mat, I0::AbstractArray{T0}, i1::Real) =
    reshape(getindex0(a, PetscInt[ to_index(i0)-1 for i0 in I0 ],
                      PetscInt[ to_index(i1)-1 ]), length(I0))


getindex{T1<:Real}(a::Mat, i0::Real, I1::AbstractArray{T1}) =
  getindex0(a, PetscInt[ to_index(i0)-1 ],
            PetscInt[ to_index(i1)-1 for i1 in I1 ])

#############################################################################
# transposition etc.

import Base: dense, transpose, ctranspose
export MatTranspose, transpose!, ctranspose!

function dense(a::Mat)
    m,n = size(a)
    a[1:m, 1:n]
end

# create a new matrix wrapping around A for matrix-vector multiplication
# but which does not actually require new storage
for (f,pf) in ((:MatTranspose,:MatCreateTranspose), # acts like A.'
               (:MatNormal, :MatCreateNormal))      # acts like A'*A
    pfe = Expr(:quote, pf)
    @eval function $f(a::Mat)
        p = Array(pMat, 1)
        chk(ccall(($pfe,petsc), PetscErrorCode, (pMat,Ptr{pMat}),
                  a, p))
        Mat(p[1], a)
    end
end

for (f,pf) in ((:transpose,:MatTranspose),(:ctranspose,:MatHermitianTranspose))
    fb = symbol(string(f,"!"))
    pfe = Expr(:quote, pf)
    @eval begin
        function $fb(a::Mat)
            chk(ccall(($pfe,petsc), PetscErrorCode, 
                      (pMat,MatReuse,Ptr{pMat}),
                      a, MAT_REUSE_MATRIX, &a.p))
            a
        end
        
        function $f(a::Mat)
            p = Array(pMat, 1)
            chk(ccall(($pfe,petsc), PetscErrorCode, 
                      (pMat,MatReuse,Ptr{pMat}),
                      a, MAT_INITIAL_MATRIX, p))
            Mat(p[1])
        end
    end
end

function conj!(a::Mat)
    chk(ccall((:MatConjugate,petsc), PetscErrorCode, (pMat,), a))
    a
end
conj(a::Mat) = conj!(copy(a))

#############################################################################
# simple math operations

function chop!(x::Mat, tol::Real)
    chk(ccall((:MatChop, petsc), PetscErrorCode, (pMat, PetscReal), x, tol))
    x
end
