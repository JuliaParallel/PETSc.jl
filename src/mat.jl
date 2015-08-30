# AbstractMatrix wrapper around Petsc Mat
export Mat
#typealias pMat Ptr{Void} # Mat arguments in Petsc are pointers
type Mat{T, MType} <: AbstractSparseMatrix{T,PetscInt}
    p::C.Mat{T}
    assembling::Bool # whether we are in the middle of assemble(vec) do ..
    insertmode::C.InsertMode # current mode for setindex!
    data::Any # keep a reference to anything needed for the Mat
    comm::MPI.Comm
    mattype::C.MatType
    function Mat(p::C.Mat{T}, data=nothing; comm=MPI.COMM_SELF)  # default sequantial matrix
        A = new(p, false, C.INSERT_VALUES, data, comm, MType)
        settype!(A, MType)
        finalizer(A, MatDestroy)
        A
    end
end

function Mat{T}(::Type{T}, mtype=C.MatType=C.MATSEQ; comm=MPI.COMM_WORLD)
    p = Array(C.Mat{T}, 1)
    chk(C.MatCreate(comm, p))
    Mat{T, mtype}(p[1])
end


function Mat{T}(::Type{T}, m::Integer, n::Integer;
    mlocal::Integer=C.PETSC_DECIDE, nlocal::Integer=C.PETSC_DECIDE,
    nz::Integer=16, nnz::AbstractVector=PetscInt[],
    onz::Integer=0, onnz::AbstractVector=PetscInt[],
    comm=MPI.COMM_WORLD,
    mtype::Symbol=C.MATMPIAIJ)

    mat = Mat(T, mtype, comm=comm)
    setsizes!(mat, m, n, mlocal=mlocal, nlocal=nlocal)
    setpreallocation!(mat, nz=nz, nnz=nnz, onz=onz, onnz=onnz)
    setoption!(mat, C.MAT_ROW_ORIENTED, false)  # julia data is column major

    return mat
end


function MatDestroy(mat::Mat)

  tmp = Array(PetscBool, 1)
  C.PetscFinalized(eltype(mat), tmp)
   
  if tmp[1] == 0  # if petsc has not been finalized yet
    C.MatDestroy([mat.p])
  end

   # if Petsc has been finalized, let the OS deallocate the memory
end



function setoption!(m::Mat, option::C.MatOption, val::Bool)
    chk(C.MatSetOption(m.p, option, PetscBool(val)))
    m
end

function settype!(a::Mat, T::C.MatType)
    chk(C.MatSetType(a.p, T))
    a
end

function gettype(a::Mat)
    p = Array(C.VecType,1)
    chk(C.MatGetType(a.p, p))
    p[1]
end

function setsizes!(a::Mat, m::Integer, n::Integer;
                   mlocal::Integer=C.PETSC_DECIDE, nlocal::Integer=C.PETSC_DECIDE)
   chk(C.MatSetSizes(a.p, mlocal, nlocal, m, n))
    a
end

function setpreallocation!{T, MType}(a::Mat{T, MType};
                           nz::Integer=16, nnz::AbstractVector=PetscInt[],
                           onz::Integer=0, onnz::AbstractVector=PetscInt[])
    if MType == C.MATSEQAIJ
        pnnz = if isempty(nnz)
            C_NULL
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
        chk(C.MatMPIAIJSetPreallocation(a.p, nz, pnnz, onz, ponnz))
    else # TODO
        throw(ArgumentError("unsupported matrix type $T"))
    end
    a
end

# construct Vec for multiplication by a::Mat or transpose(a::Mat)
#const mat2vec = [ "mpiaij" => "aij", "seqaij" => "seq" ]
const mat2vec = Dict{C.MatType, C.MatType} ( :mpiaij => :aij, :seqaij => :seq ) 
Vec{T2}(a::Mat{T2}, transposed=false) = 
  transposed ? Vec(T, size(a,1), comm=a.comm, T=mat2vec[gettype(a)],
                   mlocal=sizelocal(a,1)) :
               Vec(T, size(a,2), comm=a.comm, T=mat2vec[gettype(a)],
                   mlocal=sizelocal(a,2))

#############################################################################
import Base: nnz

convert(::Type{C.Vec}, a::Mat) = a.p

function size(a::Mat)
    m = Array(PetscInt, 1)
    n = Array(PetscInt, 1)
    chk(C.MatGetSize(a.p, m, n))
    (Int(m[1]), Int(n[1]))
end

function sizelocal(a::Mat)
    m = Array(PetscInt, 1)
    n = Array(PetscInt, 1)
    chk(C.MatGetLocalSize(a.p, m, n))
    (Int(m[1]), Int(n[1]))
end

lengthlocal(a::Mat) = prod(sizelocal(a))

function similar{T, MType}(a::Mat{T, MType})
    p = Array(C.Mat{T}, 1)
    chk(C.MatDuplicate(A.p, C.MAT_DO_NOT_COPY_VALUES, p))
    Mat{T, MType}(p[1], comm=a.comm)
end

similar{T}(a::Mat{T}, ::Type{T}) = similar(a)
similar{T}(a::Mat{T}, ::Type{T}, ::Type{PetscInt}) = similar(a)
similar{T, MType}(a::Mat{T, MType}, ::Type{T}, m::Int, n::Int) = 
    (m,n) == size(a) ? similar(a) : Mat(m,n, comm=a.comm, mtype=MType)
similar{T}(a::Mat, ::Type{T}, d::Tuple{Int,Int}) = 
    similar(a, T, d...)

function copy{T, MType}(a::Mat{T, MType})
    p = Array(C.Mat{T}, 1)
    chk(C.MatDuplicate(a.p, C.MAT_COPY_VALUES, p))
    Mat{T, MType}(p[1])
end

function getinfo(m::Mat, infotype::Integer=C.MAT_GLOBAL_SUM)
    info = Array(C.MatInfo,1)
    chk(C.MatGetInfo(m.p, MatInfoType(infotype), info))
    info[1]
end

nnz(m::Mat) = int(getinfo(m).nz_used)

#############################################################################

# for efficient matrix assembly, put all calls to A[...] = ... inside
# assemble(A) do ... end
function AssemblyBegin(x::Mat, t::C.MatAssemblyType=C.MAT_FINAL_ASSEMBLY)
  chk(C.MatAssemblyBegin(x.p, t))
end

function AssemblyEnd(x::Mat, t::C.MatAssemblyType=C.MAT_FINAL_ASSEMBLY)
  chk(C.MatAssemblyEnd(x.p, t))
end

#AssemblyBegin(x::Mat,t::C.MatAssemblyType=C.MAT_FINAL_ASSEMBLY) = chk(ccall((:MatAssemblyBegin,petsc), PetscErrorCode, (pMat,MatAssemblyType), x,t))
#AssemblyEnd(x::Mat,t::MatAssemblyType=C.MAT_FINAL_ASSEMBLY) = chk(ccall((:MatAssemblyEnd,petsc), PetscErrorCode, (pMat,MatAssemblyType), x,t))

function isassembled(x::Mat)
    if x.assembling
        return false
    else
        b = Array(PetscBool, 1)
        chk(C.MatAssembled( x.p, b))
        return b[1] != 0
    end
end

function assemble(f::Function, x::Union(Vec,Mat), 
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

# finalize matrix assembly
assemble(x::Union(Vec,Mat)) = assemble(() -> nothing, x)

# intermediate assembly, before it is finally compressed for use
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
        AssemblyBegin(x,C.MAT_FINAL_ASSEMBLY)
        AssemblyEnd(x,C.MAT_FINAL_ASSEMBLY)
    end
    x
end

function setindex!{T}(x::Mat{T}, v::Number, i::Real, j::Real)
    # can't call MatSetValue since that is a static inline function
    setindex0!(x, T[ v ], 
               PetscInt[ to_index(i) - 1 ], 
               PetscInt[ to_index(j) - 1 ])
    v
end
function setindex!{T3, T1<:Real, T2<:Real}(x::Mat{T3}, v::Array{T3}, 
                                       I::AbstractArray{T1},
                                       J::AbstractArray{T2})
    I0 = PetscInt[ to_index(i)-1 for i in I ]
    J0 = PetscInt[ to_index(j)-1 for j in J ]
    setindex0!(x, v, I0, J0)
end
function setindex!{T2, T<:Real}(x::Mat{T2}, v::Array{T2}, 
                            i::Real, J::AbstractArray{T})
    I0 = PetscInt[ to_index(i)-1 ]
    J0 = PetscInt[ to_index(j)-1 for j in J ]
    setindex0!(x, v, I0, J0)
end
function setindex!{T2, T<:Real}(x::Mat{T2}, v::Array{T2}, 
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
    v0 = eltype(Mat)[ z for z in v ]
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
function getindex0{T}(x::Mat{T}, i::Vector{PetscInt}, j::Vector{PetscInt})
    ni = length(i)
    nj = length(j)
    v = Array(T, nj, ni) # row-major!
    chk(C.MatGetValues(x.p, ni, i, nj, j, v))

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

import Base: full, transpose, ctranspose
export MatTranspose, transpose!, ctranspose!

function full(a::Mat)
    m,n = size(a)
    a[1:m, 1:n]
end

# create a new matrix wrapping around A for matrix-vector multiplication
# but which does not actually require new storage
for (f,pf) in ((:MatTranspose,:MatCreateTranspose), # acts like A.'
               (:MatNormal, :MatCreateNormal))      # acts like A'*A
    pfe = Expr(:quote, pf)
    @eval function $f{T}(a::Mat{T})
        p = Array(C.Mat{T}, 1)
        chk($pf(a.p, p))
        Mat(p[1], a)
    end
end

for (f,pf) in ((:transpose,:MatTranspose),(:ctranspose,:MatHermitianTranspose))
    fb = symbol(string(f,"!"))
    pfe = Expr(:quote, pf)
    @eval begin
        function $fb(a::Mat)
            pa = [a.p]
            chk(C.$pf(a.p, C.MAT_REUSE_MATRIX, pa))
            a
        end
        
        function $f{T}(a::Mat{T})
            p = Array(C.Mat{T}, 1)
            chk(C.$pf(a.p, C.MAT_INITIAL_MATRIX, p))
            Mat(p[1])
        end
    end
end

function conj!(a::Mat)
    chk(C.MatConjugate(a.p))
    a
end
conj(a::Mat) = conj!(copy(a))

#############################################################################
# simple math operations
# skip for now
#=
function chop!(x::Mat, tol::Real)
    chk(ccall((:MatChop, petsc), PetscErrorCode, (pMat, PetscReal), x, tol))
    x
end
=#
