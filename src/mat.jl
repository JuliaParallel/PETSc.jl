# AbstractMatrix wrapper around Petsc Mat
export Mat, petscview, SubMat

abstract PetscMat{T, MType} <: AbstractSparseMatrix{T, PetscInt}
type Mat{T, MType} <: PetscMat{T, MType}
  p::C.Mat{T}
  assembling::Bool # whether we are in the middle of assemble(vec)
  insertmode::C.InsertMode # current mode for setindex!
  data::Any # keep a reference to anything needed for the Mat
            # -- needed if the Mat is a wrapper around a Julia object,
            #    to prevent the object from being garbage collected.
  function Mat(p::C.Mat{T}, data=nothing; first_instance::Bool=true)
    A = new(p, false, C.INSERT_VALUES, data)
    if first_instance  # if the pointer p has not been put into a Mat before
      chk(C.MatSetType(p, MType))
      finalizer(A, PetscDestroy)
    end
    
    return A
  end
end

type SubMat{T, MType} <: PetscMat{T}
  p::C.Mat{T}
  assembling::Bool # whether we are in the middle of assemble(vec)
  insertmode::C.InsertMode # current mode for setindex!
  data::Any # keep a reference to anything needed for the Mat
            # -- needed if the Mat is a wrapper around a Julia object,
            #    to prevent the object from being garbage collected.
            # in general, we *must* keep a copy of the parent matrix, because
            # creating a submatrix does not increase the reference count
            # in Petsc (for some unknown reason)
  function SubMat(p::C.Mat{T}, data=nothing)
    A = new(p, false, C.INSERT_VALUES, data)
    finalizer(A, SubMatRestore)
    return A
  end
end


"""
  Get the communicator for the object
"""
comm{T}(a::PetscMat{T}) = MPI.Comm(C.PetscObjectComm(T, a.p.pobj))

"""
  Create an empty, unsized matrix
"""
function Mat{T}(::Type{T}, mtype::C.MatType=C.MATSEQAIJ; comm::MPI.Comm=MPI.COMM_WORLD)
  p = Ref{C.Mat{T}}()
  chk(C.MatCreate(comm, p))
  Mat{T, mtype}(p[])
end

"""
  Create a matrix of a particular size, optionally specifying the pre-allocation.
  If pre-allocation is not specified, no preallocation is done
"""
function Mat{T}(::Type{T}, m::Integer=C.PETSC_DECIDE, n::Integer=C.PETSC_DECIDE;
                mlocal::Integer=C.PETSC_DECIDE, nlocal::Integer=C.PETSC_DECIDE,
                bs=1, nz::Integer=0, nnz::AbstractVector=PetscInt[],
                onz::Integer=0, onnz::AbstractVector=PetscInt[],
                comm::MPI.Comm=MPI.COMM_WORLD,
                mtype::Symbol=C.MATMPIAIJ)

  mat = Mat(T, mtype, comm=comm)
  resize!(mat, m, n, mlocal=mlocal, nlocal=nlocal)

  # don't preallocate unless the user specified to do so, because
  # after pre-allocation it is an error to try to change the sparsity pattern
  # Before preallocation Petsc will dynamically allocate memory as needed
  # slow but flexible.
  if nz==0 && onz == 0  && nnz == PetscInt[] && onnz == PetscInt[]
    if bs != 1
      chk(C.MatSetBlockSize(mat.p, PetscInt(bs)))
    end
    chk(C.MatSetUp(mat.p)) 
  else  # preallocate
    setpreallocation!(mat, nz=nz, nnz=nnz, onz=onz, onnz=onnz, bs=bs)
  end
  setoption!(mat, C.MAT_ROW_ORIENTED, false)  # julia data is column major

  return mat
end


"""
  Gets the a submatrix that references the entries in the original matrix.
  isrow and iscol contain the *local* indicies of the rows and columns to get.
  The matrix must have a LocalToGlobalMapping for this to work, therefore a 
  default one is created if the matrix does not already have one registered.
  The default mapping assumes the matrix is divided up into contiguous block
  of rows.  This is true of AIJ matrices but may not be for other matrix types.
"""
function SubMat{T, MType}(mat::Mat{T, MType}, isrow::IS{T}, iscol::IS{T})

  # create the local to global mapping for mat first
  # this only needs to be done the first time
  if !has_local_to_global_mapping(mat)
    rmap, cmap = local_to_global_mapping(mat)
    set_local_to_global_mapping(mat, rmap, cmap)  
  end

  # now we can actually create the submatrix
  submat = Ref{C.Mat{T}}()
  chk(C.MatGetLocalSubMatrix(mat.p, isrow.p, iscol.p, submat))
  # keep the data needed for the finalizer
  return SubMat{T, MType}(submat[], (mat, isrow, iscol))  
end

export SubMatRestore
function SubMatRestore{T}(smat::SubMat{T})

  C.MatRestoreLocalSubMatrix(smat.data[1].p, smat.data[2].p, smat.data[3].p, Ref(smat.p))
end

##### MatShell functions #####
export MatShell, setop!, getcontext
"""
  Create a high level matrix from an already created matrix pointer
"""
function Mat{T}(ptr::C.Mat{T})
# this is not type stable. Grr
  sym_arr = Array(Symbol, 1)
  chk(C.MatGetType(ptr, sym_arr))
  mtype = sym_arr[1]
  return Mat{T, mtype}(ptr, nothing, first_instance=false)
end


"""
  Create a shell matrix with specified size.  The ctx tuple contains can be
  accessed by any callback function.
"""
# must rename this because of (silent) method ambiguity
function MatShell{T}(::Type{T}, mlocal::Integer, nlocal::Integer, ctx::Tuple=();  m::Integer=C.PETSC_DECIDE, n::Integer=C.PETSC_DECIDE, comm=MPI.COMM_WORLD)

  mat_ptr = Ref{C.Mat{T}}()
  ctx_ptr = pointer_from_objref(ctx)
  chk(C.MatCreateShell(comm, mlocal, nlocal, m, n, ctx_ptr, mat_ptr))
  return Mat{T, C.MATSHELL}(mat_ptr[], ctx)  # protect ctx from gc
end

"""
  Provide a callback function for a particular matrix operation.  op is a 
  Petsc enum value inidcating the operation, and func is a void pointer 
  (obtained from cfunction() ) that performs the operation.

  The function should take the low level Petsc objects (defined in the C module)
  rather than the high level ones defined in this file.  There are constructors
  to create a high level object from a low level one

"""
function setop!{T}(mat::Mat{T, C.MATSHELL}, op::C.MatOperation, func::Ptr{Void})

  chk(C.MatShellSetOperation(mat.p, op, func))
end

"""
  Get the tuple of user provided data passed in when the shell matrix was 
  created.
"""
function getcontext{T}(mat::Mat{T, C.MATSHELL})

  ctx_ptr = C.MatShellGetContext(mat.p)
  return unsafe_pointer_to_objref(ctx_ptr)
end

"""
  Destroy a Mat object and the underlying data structure, if the object
  has not already been finalized
"""
function PetscDestroy{T}(mat::Mat{T})
  if !PetscFinalized(T)
    C.MatDestroy(Ref(mat.p))
    mat.p = C.Mat{T}(C_NULL)  # indicate the vector is finalized
  end
end

"""
  Check if PetscDestroy has been called on this object already
"""
function isfinalized(mat::PetscMat)
  return isfinalized(mat.p)
end

function isfinalized(mat::C.Mat)
  return mat.pobj == C_NULL
end

"""
  Print a Petsc matrix to STDOUT
"""
function petscview{T}(mat::PetscMat{T})
  viewer = C.PetscViewer{T}(C_NULL)
  chk(C.MatView(mat.p, viewer))
end

"""
  Print a Petsc matrix to a named file, in text format
"""
function petscwrite{T}(mat::PetscMat{T}, fname)
  viewer_ref = Ref{C.PetscViewer{T}}()
  chk(C.PetscViewerASCIIOpen(comm(mat), fname, viewer_ref))
  chk(C.MatView(mat.p, viewer_ref[]))
  chk(C.PetscViewerDestroy(viewer_ref))
end


export setoption!, gettype

"""
  Pass values to the Petsc function MatSetOption.  Note that the handful of 
  options that can be passed here should not be confused with those for the
  global options database
"""
function setoption!(m::Mat, option::C.MatOption, val::Bool)
  chk(C.MatSetOption(m.p, option, PetscBool(val)))
  m
end

"""
  Get the format of the matrix.
"""
gettype{T,MT}(a::PetscMat{T,MT}) = MT

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

"""
  Preallocates the sparsity pattern for (B)AIJ matrices.
"""
function setpreallocation!{T, MType}(a::Mat{T, MType};
  nz::Integer=16, nnz::AbstractVector=PetscInt[],
  onz::Integer=0, onnz::AbstractVector=PetscInt[],
  bs::Integer=1)
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

"""
  Maps Matrix formats to the corresponding vector format
"""
# construct Vec for multiplication by a::Mat or transpose(a::Mat)
const mat2vec = Dict{C.MatType, C.MatType}( :mpiaij => :aij, :seqaij => :seq )

"""
  Construct a vector suitable for multiplying by the given matrix
"""
Vec{T2}(a::PetscMat{T2}, transposed=false) =
  transposed ? Vec(T2, size(a,1), comm=comm(a), T=mat2vec[gettype(a)],
  mlocal=sizelocal(a,1)) :
  Vec(T2, size(a,2), comm=comm(a), T=mat2vec[gettype(a)],
  mlocal=sizelocal(a,2))

#############################################################################
Base.convert(::Type{C.Mat}, a::PetscMat) = a.p

export sizelocal, localranges, lengthlocal

"""
  Returns the global dimensions of the matrix
"""
function Base.size(a::PetscMat)
  m = Ref{PetscInt}()
  n = Ref{PetscInt}()
  chk(C.MatGetSize(a.p, m, n))
  (Int(m[]), Int(n[]))
end

"""
  Returns the local dimensions of the matrix
"""
function sizelocal(a::PetscMat)
  m = Ref{PetscInt}()
  n = Ref{PetscInt}()
  chk(C.MatGetLocalSize(a.p, m, n))
  (Int(m[]), Int(n[]))
end

"""
  This function returns two Range object corresponding the global indices
  of the rows and columns of the matrix A.

  The function has the same limiations as Petsc's MatGetOwnershipRange, 
  in that it assumes the rows of the matrix are divided up contigously.
"""
function localranges(a::PetscMat)
  start_ref = Ref{PetscInt}()
  end_1_ref = Ref{PetscInt}()
  chk(C.MatGetOwnershipRange(a.p, start_ref, end_1_ref))
  start_idx = start_ref[] + 1  # convert to 1 based index
  end_idx = end_1_ref[]  # subtract 1 because petsc supplies 1 past the end,
                         # then add 1 to convert to 1 based index

  return start_idx:end_idx, 1:size(a, 2)
end

"""
  Constructs 2 index sets that map from the local row and columns to the
  global rows and columns
"""
function localIS{T}(A::PetscMat{T})

  rows, cols = localranges(A)
  rowis = IS(T, rows, comm=comm(A))
  colis = IS(T, cols, comm=comm(A))
  return rowis, colis
end

"""
  Construct ISLocalToGlobalMappings for the the rows and columns of the matrix
"""
function local_to_global_mapping(A::PetscMat)

  # localIS creates strided index sets, which require only constant
  # memory 
  rowis, colis = localIS(A)
  row_ltog = ISLocalToGlobalMapping(rowis)
  col_ltog = ISLocalToGlobalMapping(colis)
  return row_ltog, col_ltog
end

# need a better name
"""
  Registers the ISLocalToGlobalMappings with the matrix
"""
function set_local_to_global_mapping{T}(A::PetscMat{T}, rmap::ISLocalToGlobalMapping{T}, cmap::ISLocalToGlobalMapping{T})

  chk(C.MatSetLocalToGlobalMapping(A.p, rmap.p, cmap.p))
end

"""
  Check if the local to global mappings have been registered
"""
function has_local_to_global_mapping{T}(A::PetscMat{T})

  rmap_ref = Ref{C.ISLocalToGlobalMapping{T}}()
  cmap_ref = Ref{C.ISLocalToGlobalMapping{T}}()
  chk(C.MatGetLocalToGlobalMapping(A.p, rmap_ref, cmap_ref))

  rmap = rmap_ref[]
  cmap = cmap_ref[]
  
  return rmap.pobj != C_NULL && cmap.pobj != C_NULL
end

"""
  prod(sizelocal))
"""
lengthlocal(a::PetscMat) = prod(sizelocal(a))

# this causes the assembly state of the underlying petsc matrix to be copied
function Base.similar{T, MType}(a::PetscMat{T, MType})
  p = Ref{C.Mat{T}}()
  chk(C.MatDuplicate(a.p, C.MAT_DO_NOT_COPY_VALUES, p))
  Mat{T, MType}(p[])
end

Base.similar{T}(a::PetscMat{T}, ::Type{T}) = similar(a)
Base.similar{T,MType}(a::PetscMat{T,MType}, T2::Type) =
  Mat(T2, size(a)..., comm=comm(a), mtype=MType)
Base.similar{T,MType}(a::PetscMat{T,MType}, T2::Type, m::Integer, n::Integer) =
  (m,n) == size(a) && T2==T ? similar(a) : Mat(T2, m,n, comm=comm(a), mtype=MType)
Base.similar{T,MType}(a::PetscMat{T,MType}, m::Integer, n::Integer) = similar(a, T, m, n)
Base.similar(a::PetscMat, T::Type, d::Dims) = similar(a, T, d...)
Base.similar{T}(a::PetscMat{T}, d::Dims) = similar(a, T, d)

function Base.copy{T,MType}(a::PetscMat{T,MType})
  p = Ref{C.Mat{T}}()
  chk(C.MatDuplicate(a.p, C.MAT_COPY_VALUES, p))
  Mat{T,MType}(p[])
end

"""
  Get the MatInfo struct for the matrix
"""
function getinfo(m::Mat, infotype::Integer=C.MAT_GLOBAL_SUM)
  info = Ref{C.MatInfo}()
  chk(C.MatGetInfo(m.p, C.MatInfoType(infotype), info))
  info[]
end

"""
  Number of non-zero entries that have been assigned to
"""
Base.nnz(m::Mat) = Int(getinfo(m).nz_used)

#############################################################################

# for efficient matrix assembly, put all calls to A[...] = ... inside
# assemble(A) do ... end
"""
  Start assembling the matrix (the implmentations probably post 
  non-blocking sends and received)
"""
function AssemblyBegin(x::PetscMat, t::C.MatAssemblyType=C.MAT_FLUSH_ASSEMBLY)
  chk(C.MatAssemblyBegin(x.p, t))
end

"""
  Finish assembling the matrix
"""
function AssemblyEnd(x::PetscMat, t::C.MatAssemblyType=C.MAT_FLUSH_ASSEMBLY)
  chk(C.MatAssemblyEnd(x.p, t))
end

"""
  Check if the matrix is assembled or not
"""
function isassembled(p::C.Mat)
  b = Ref{PetscBool}()
  chk(C.MatAssembled(p, b))
  return b[] != 0
end

"""
  Check if the matrix is assembled
"""
# why does this check x.assembling?
isassembled(x::PetscMat) = !x.assembling && isassembled(x.p)

"""
  This function provides a mechanism for efficiently inserting values into
  and then assembling Petsc matrices and vectors.  The function f must be a 
  zero argument function.

  This function can be used with the do block syntax.
"""
function assemble(f::Function, x::Union{Vec,PetscMat},
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
"""
  Assemble the Petsc object
"""
function assemble(x::Union{Vec,PetscMat}, t::C.MatAssemblyType=C.MAT_FINAL_ASSEMBLY)
  AssemblyBegin(x, t)
  AssemblyEnd(x,t)
end


# in ksp solve we need to finalize assembly from raw pointer:
"""
  Low level assemble function
"""
function assemble(p::C.Mat, t::C.MatAssemblyType=C.MAT_FINAL_ASSEMBLY)
  if !isassembled(p)
    chk(C.MatAssemblyBegin(p, t))
    chk(C.MatAssemblyEnd(p, t))
  end
  return nothing
end

# intermediate assembly, before it is finally compressed for use
"""
  Perform a flush assembly (take stashed values and put them into the matrix,
  but don't squeeze out any preallocated space that has not been used yet
"""
iassemble(x::Union{Vec,PetscMat}) = assemble(() -> nothing, x, x.insertmode, C.MAT_FLUSH_ASSEMBLY)
iassemble(f::Function, x::PetscMat, insertmode=x.insertmode) =
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

# like setindex0 above, but for submatrices.  the indices i and j must 
# be *local* indices
function setindex0!{T}(x::SubMat{T}, v::Array{T}, i::Array{PetscInt}, j::Array{PetscInt})
  ni = length(i)
  nj = length(j)
  if length(v) != ni*nj
    throw(ArgumentError("length(values) != length(indices)"))
  end
  chk(C.MatSetValuesLocal(x.p, ni, i, nj, j, v, x.insertmode))
  if !x.assembling
    AssemblyBegin(x,C.MAT_FLUSH_ASSEMBLY)
    AssemblyEnd(x,C.MAT_FLUSH_ASSEMBLY)
  end
  x
end

import Base: setindex!

function setindex!{T}(x::PetscMat{T}, v::Number, i::Integer, j::Integer)
  # can't call MatSetValue since that is a static inline function
  setindex0!(x, T[ v ],
  PetscInt[ i - 1 ],
  PetscInt[ j - 1 ])
  v
end
function setindex!{T3, T1<:Integer, T2<:Integer}(x::PetscMat{T3}, v::Array{T3},
  I::AbstractArray{T1},
  J::AbstractArray{T2})
  I0 = PetscInt[ i-1 for i in I ]
  J0 = PetscInt[ j-1 for j in J ]
  setindex0!(x, v, I0, J0)
end
function setindex!{T2, T<:Integer}(x::PetscMat{T2}, v::Array{T2},
  i::Integer, J::AbstractArray{T})
  I0 = PetscInt[ i-1 ]
  J0 = PetscInt[ j-1 for j in J ]
  setindex0!(x, v, I0, J0)
end
function setindex!{T2, T<:Real}(x::PetscMat{T2}, v::Array{T2},
  I::AbstractArray{T}, j::Integer)
  I0 = PetscInt[ i-1 for i in I ]
  J0 = PetscInt[ j-1 ]
  setindex0!(x, v, I0, J0)
end

setindex!{T1<:Integer, T2<:Integer}(x::PetscMat, v::Number,
I::AbstractArray{T1},
J::AbstractArray{T2}) = iassemble(x) do
  for i in I
    for j in J
      x[i,j] = v
    end
  end
  x
end
setindex!{T<:Integer}(x::PetscMat, v::Number,
i::Integer, J::AbstractArray{T}) = iassemble(x) do
  for j in J
    x[i,j] = v
  end
  x
end
setindex!{T<:Integer}(x::PetscMat, v::Number,
I::AbstractArray{T}, j::Integer) = iassemble(x) do
  for i in I
    x[i,j] = v
  end
  x
end

function setindex!{T0<:Number, T1<:Integer, T2<:Integer}(x::PetscMat, v::AbstractArray{T0},
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
"""
  Fill the matrix with the specified values.

  Currently, this function either destroys the sparsity pattern or 
  gives an error, unless v = 0, in which case it zeros out the non-zero 
  entries without changing the sparsity pattern
"""
function Base.fill!(x::PetscMat, v::Number)
  if v == 0
    chk(C.MatZeroEntries(x.p))
  else
    # FIXME: don't write to non-allocated entries
    iassemble(x) do
      rows, cols = localranges(x)
      for i in rows
        for j in cols
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

getindex(a::PetscMat, i0::Integer, i1::Integer) =
  getindex0(a, PetscInt[ i0-1 ], PetscInt[ i1-1 ])[1]

getindex{T0<:Integer,T1<:Integer}(a::PetscMat, I0::AbstractArray{T0}, I1::AbstractArray{T1}) =
  getindex0(a, PetscInt[ i0-1 for i0 in I0 ], PetscInt[ i1-1 for i1 in I1 ])

getindex{T0<:Integer}(a::PetscMat, I0::AbstractArray{T0}, i1::Integer) =
  reshape(getindex0(a, PetscInt[ i0-1 for i0 in I0 ], PetscInt[ i1-1 ]), length(I0))

getindex{T1<:Integer}(a::PetscMat, i0::Integer, I1::AbstractArray{T1}) =
  getindex0(a, PetscInt[ i0-1 ], PetscInt[ i1-1 for i1 in I1 ])

#############################################################################
# transposition etc.

export MatTranspose

"""
  Create a dense Julia matrix for a Petsc sparse matrix.  This only works
  for SEQ matrices
"""
function Base.full(a::PetscMat)
  m,n = size(a)
  a[1:m, 1:n]
end

"""
 create a new matrix wrapping around A for matrix-vector multiplication
 but which does not actually require new storage
"""
for (f,pf) in ((:MatTranspose,:MatCreateTranspose), # acts like A.'
  (:MatNormal, :MatCreateNormal))      # acts like A'*A
  pfe = Expr(:quote, pf)
  @eval function $f{T, MType}(a::PetscMat{T, MType})
    p = Ref{C.Mat{T}}()
    chk(C.$pf(a.p, p))
    Mat{T, MType}(p[], a)
  end
end

for (f,pf) in ((:transpose,:MatTranspose),(:ctranspose,:MatHermitianTranspose))
  fb = symbol(string(f,"!"))
  pfe = Expr(:quote, pf)
  @eval begin
    function Base.$fb(a::PetscMat)
      pa = [a.p]
      chk(C.$pf(a.p, C.MAT_REUSE_MATRIX, pa))
      a
    end

    function Base.$f{T,MType}(a::PetscMat{T,MType})
      p = Ref{C.Mat{T}}()
      chk(C.$pf(a.p, C.MAT_INITIAL_MATRIX, p))
      Mat{T,MType}(p[])
    end
  end
end

function Base.conj!(a::PetscMat)
  chk(C.MatConjugate(a.p))
  a
end
Base.conj(a::PetscMat) = conj!(copy(a))

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

function Base.trace{T}(A::PetscMat{T})
  t = Ref{T}()
  chk(C.MatGetTrace(A.p,t))
  return t[]
end

function Base.real{T<:Complex}(A::PetscMat{T})
  N = copy(A)
  chk(C.MatRealPart(N.p))
  return N
end
Base.real{T<:Real}(A::PetscMat{T}) = A

function Base.imag{T<:Complex}(A::PetscMat{T})
  N = copy(A)
  chk(C.MatImaginaryPart(N.p))
  return N
end

function Base.LinAlg.ishermitian{T}(A::PetscMat{T}, tol::Real=eps(real(float(one(T)))))
  bool_arr = Ref{PetscBool}()
  chk(C.MatIsHermitian(A.p, tol, bool_arr))
  return bool_arr[] != 0
end

function Base.LinAlg.issym{T}(A::PetscMat{T}, tol::Real=eps(real(float(one(T)))))
  bool_arr = Ref{PetscBool}()
  chk(C.MatIsSymmetric(A.p, tol, bool_arr))
  return bool_arr[] != 0
end

#currently ONLY gets the main diagonal
function Base.diag{T}(A::PetscMat{T},vtype::C.VecType=C.VECSEQ)
  m = size(A, 1)
  b = Vec(T, m, vtype=vtype, comm=comm(A), mlocal=sizelocal(A,1))
  chk(C.MatGetDiagonal(A.p,b.p))
  return b
end

function (*){T, MType}(A::PetscMat{T}, x::Vec{T, MType})
  m = size(A, 1)
  b = Vec(T, m, vtype=MType, comm=comm(A), mlocal=sizelocal(A,1))
  chk(C.MatMult(A.p, x.p, b.p))
  return b
end

function (*){T, MType}(A::PetscMat{T}, x::Vec{T, MType}, b::Vec{T})
  chk(C.MatMult(A.p, x.p, b.p))
  return b
end



function (.*){T, MType}(A::PetscMat{T, MType}, x::Number)
  Y = copy(A)
  chk(C.MatScale(Y.p, T(x)))
  return Y
end
(.*){T, MType}(x::Number, A::PetscMat{T, MType}) = A .* x
(./){T, MType}(A::PetscMat{T, MType}, x::Number) = A .* inv(x)
(.\){T, MType}(x::Number, A::PetscMat{T, MType}) = A .* inv(x)

function (*){T, MType}(A::PetscMat{T,MType}, B::PetscMat{T})
  p = Ptr{Float64}(0)
  p_arr = Ref{C.Mat{T}}()
  chk(C.MatMatMult(A.p, B.p, C.MAT_INITIAL_MATRIX, real(T)(C.PETSC_DEFAULT), p_arr))
  new_mat = Mat{T, MType}(p_arr[])
  return new_mat
end

#these two only work for SEQAIJ
function At_mul_B{T}(A::PetscMat{T}, B::PetscMat{T})
  p = Ptr{Float64}(0)
  p_arr = Ref{C.Mat{T}}()
  chk(C.MatTransposeMatMult(A.p, B.p, C.MAT_INITIAL_MATRIX, real(T)(C.PETSC_DEFAULT), p_arr))
  new_mat = Mat{T, C.MATSEQAIJ}(p_arr[])
  return new_mat
end

function A_mul_Bt{T}(A::PetscMat{T}, B::PetscMat{T})
  p = Ptr{Float64}(0)
  p_arr = Ref{C.Mat{T}}()
  chk(C.MatMatTransposeMult(A.p, B.p, C.MAT_INITIAL_MATRIX, real(T)(C.PETSC_DEFAULT), p_arr))
  new_mat = Mat{T, C.MATSEQAIJ}(p_arr[])
  return new_mat
end

function At_mul_B!{T, MType}(A::PetscMat{T}, x::Vec{T,MType}, y::Vec{T,MType})
  chk(C.MatMultTranspose(A.p, x.p, y.p))
  return y
end

function At_mul_B{T, MType}(A::PetscMat{T}, x::Vec{T,MType})
  m = size(A, 1)
  b = Vec(T, m, vtype=MType, comm=comm(A), mlocal=sizelocal(A,1))
  chk(C.MatMultTranspose(A.p, x.p, b.p))
  return b
end

function Ac_mul_B!{T, MType}(A::PetscMat{T}, x::Vec{T,MType}, y::Vec{T,MType})
  chk(C.MatMultHermitianTranspose(A.p, x.p, y.p))
  return y
end

function Ac_mul_B{T, MType}(A::PetscMat{T}, x::Vec{T,MType})
  m = size(A, 1)
  b = Vec(T, m, vtype=MType, comm=comm(A), mlocal=sizelocal(A,1))
  chk(C.MatMultHermitianTranspose(A.p, x.p, b.p))
  return b
end

for (f,s) in ((:+,1), (:-,-1))
  @eval function ($f){T}(A::PetscMat{T}, B::PetscMat{T})
    Y = copy(A)
    chk(C.MatAXPY(Y.p, T($s), B.p, C.DIFFERENT_NONZERO_PATTERN))
    return Y
  end
end

(-){T, MType}(A::PetscMat{T, MType}) = A * (-1)

# there don't appear to be PETSc functions for pointwise
# operations on matrices

function (==){T}(A::PetscMat{T}, b::PetscMat{T})
  bool_arr = Ref{PetscBool}()
  chk(C.MatEqual(A.p, b.p, bool_arr))
  return bool_arr[] != 0
end
