# AbstractMatrix wrapper around Petsc Mat
export Mat, petscview, SubMat

"""
  A Petsc matrix.

  Unlike Vecs, the Petsc implementation keeps track of the local assembly 
  state, so the Julia type does not have to.
  `verify_assembled`: if true, verify all processes are assembled, if false,
                      only local process
  `insertmode`: C.InsertMode used by `setindex!`
"""  
abstract PetscMat{T} <: AbstractSparseMatrix{T, PetscInt}
type Mat{T} <: PetscMat{T}
  p::C.Mat{T}
  verify_assembled::Bool # check all processes assembled state or just current
  insertmode::C.InsertMode # current mode for setindex!
  data::Any # keep a reference to anything needed for the Mat
            # -- needed if the Mat is a wrapper around a Julia object,
            #    to prevent the object from being garbage collected.
  function Mat(p::C.Mat{T}, data=nothing; first_instance::Bool=true, verify_assembled=true)
    A = new(p, verify_assembled, C.INSERT_VALUES, data)
    if first_instance  # if the pointer p has not been put into a Mat before
      finalizer(A, PetscDestroy)
    end
    
    return A
  end
end

type SubMat{T} <: PetscMat{T}
  p::C.Mat{T}
  verify_assembled::Bool
  insertmode::C.InsertMode # current mode for setindex!
  data::Any # keep a reference to anything needed for the Mat
            # -- needed if the Mat is a wrapper around a Julia object,
            #    to prevent the object from being garbage collected.
            # in general, we *must* keep a copy of the parent matrix, because
            # creating a submatrix does not increase the reference count
            # in Petsc (for some unknown reason)
  function SubMat(p::C.Mat{T}, data=nothing; verify_assembled=true)
    A = new(p, verify_assembled, C.INSERT_VALUES, data)
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
  chk(C.MatSetType(p[], mtype))
  Mat{T}(p[])
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
      set_block_size(mat, bs)
    end
    chk(C.MatSetUp(mat.p)) 
  else  # preallocate
    setpreallocation!(mat, nz=nz, nnz=nnz, onz=onz, onnz=onnz, bs=bs)
  end
  setoption!(mat, C.MAT_ROW_ORIENTED, false)  # julia data is column major

  return mat
end

"""
  Make a MATSEQ Petsc matrix for a SparseMatrixCSC.  This preserve the 
  sparsity pattern of the matrix
"""
function Mat{T}(A::SparseMatrixCSC{T})
  m, n = size(A)

  # count number of non-zeros in each row
  nz = zeros(PetscInt, m)
  for i=1:n  # loop over columns
    idx_start = A.colptr[i]
    idx_end = A.colptr[i+1]-1
    for j=idx_start:idx_end
      rownum = A.rowval[j]
      nz[rownum] += 1
    end
  end

  # create the matrix
  PA = Mat(T, m, n, nnz=nz, mtype=C.MATSEQAIJ)

  # copy values
  # create array large enough to hold all values in a column
  maxcol = 0
  for i=1:n
    nvals = A.colptr[i+1] - A.colptr[i]
    if nvals > maxcol
      maxcol = nvals
    end
  end
  idx = Array(PetscInt, maxcol)
  idy = PetscInt[0]
  
  for i=1:n
    idy[1] = i-1
    idx_start = A.colptr[i]
    idx_end = A.colptr[i+1]-1
    # convert indices to PetscInts, zero-based
    pos=1
    for j=idx_start:idx_end
      idx[pos] = A.rowval[j] - 1
      pos += 1
    end

    # get subarray of only the needed entries
    idx_view = sub(idx, 1:(pos-1))
    vals = sub(A.nzval, idx_start:idx_end)
    set_values!(PA, idx_view, idy, vals)
  end

  return PA
end

"""
  Construct at MATSEQAIJ from an AbstractArray.  The argument droptol is
  used to determine what size entry is considered non-zero
"""
function Mat{T}(A::AbstractArray{T}; droptol=0.0)

  m, n = size(A)

  # count non-zeros
  nz = zeros(PetscInt, m)
  for i=1:n
    for j=1:m
      val = A[j, i]
      if abs(val) > droptol
        nz[j] += 1
      end
    end
  end

  # create matrix
  PA = Mat(T, m, n, nnz=nz, mtype=C.MATSEQAIJ)

  # copy values
  maxrow = maximum(nz)
  idx = PetscInt[0]
  idy = zeros(PetscInt, maxrow)
  vals = zeros(T, maxrow)
  for i=1:m
    idx[1] = i-1
    pos = 1
    for j=1:n
      val = A[i, j]
      if abs(val) > droptol
        idy[pos] = j - 1
        vals[pos] = val
        pos += 1
      end
    end
    idy_view = sub(idy, 1:(pos-1))
    set_values!(PA, idx, idy_view, vals)
  end

  return PA
end




"""
  Gets the a submatrix that references the entries in the original matrix.
  isrow and iscol contain the *local* indicies of the rows and columns to get.
  The matrix must have a LocalToGlobalMapping for this to work, therefore a 
  default one is created if the matrix does not already have one registered.
  The default mapping assumes the matrix is divided up into contiguous block
  of rows.  This is true of AIJ matrices but may not be for other matrix types.
"""
function SubMat{T}(mat::Mat{T}, isrow::IS{T}, iscol::IS{T})

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
  return SubMat{T}(submat[], (mat, isrow, iscol), verify_assembled=mat.verify_assembled)  
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
#  sym_arr = Array(Symbol, 1)
#  chk(C.MatGetType(ptr, sym_arr))
#  mtype = sym_arr[1]
  return Mat{T}(ptr, nothing, first_instance=false)
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
  #TODO: remove unneeded MatSetType
  chk(C.MatSetType(mat_ptr[], C.MATSHELL))
  return Mat{T}(mat_ptr[], ctx)  # protect ctx from gc
end

"""
  Provide a callback function for a particular matrix operation.  op is a 
  Petsc enum value inidcating the operation, and func is a void pointer 
  (obtained from cfunction() ) that performs the operation.

  The function should take the low level Petsc objects (defined in the C module)
  rather than the high level ones defined in this file.  There are constructors
  to create a high level object from a low level one

"""
function setop!{T}(mat::Mat{T}, op::C.MatOperation, func::Ptr{Void})
  if gettype(mat) != C.MATSHELL
    throw(ArgumentError("Mat must be a MatShell"))
  end

  chk(C.MatShellSetOperation(mat.p, op, func))
end

"""
  Get the tuple of user provided data passed in when the shell matrix was 
  created.
"""
function getcontext(mat::Mat)
  if gettype(mat) != C.MATSHELL
    throw(ArgumentError("Mat must be a MatShell"))
  end

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


function set_block_size{T<:Scalar}(A::Mat{T}, bs::Integer)
  chk(C.MatSetBlockSize(A.p, bs))
end

function get_blocksize{T<:Scalar}(A::Mat{T})
  bs = Ref{PetscInt}()
  chk(C.MatGetBlockSize(A.p, bs))
  return Int(bs[])
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
function gettype(a::PetscMat)
  sym_arr = Array(C.MatType, 1)
  chk(C.MatGetType(a.p, sym_arr))
  return sym_arr[1]
end

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
function setpreallocation!{T}(a::Mat{T};
  nz::Integer=16, nnz::AbstractVector=PetscInt[],
  onz::Integer=0, onnz::AbstractVector=PetscInt[],
  bs::Integer=1)
  MType = gettype(a)
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
  Similar to localpart, but returns the range of block indices
"""
function localpart_block(A::Mat)
  low = Ref{PetscInt}()
  high = Ref{PetscInt}()
  chk(C.MatGetOwnershipRange(A.p, low, high))
  bs = get_blocksize(A)
  low_b = div(low[], bs); high_b = div(high[]-1, bs)
  rows = (low_b+1):(high_b+1)
  cols = 1:div(size(A, 2), bs)

  return rows, cols
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
  Like localIS, but returns a block index IS
"""
function localIS_block{T}(A::Mat{T})
  rows, cols = localpart_block(A)
  bs = get_blocksize(A)
  rowis = ISBlock(T, bs, rows, comm=comm(A))
  colis = ISBlock(T, bs, cols, comm=comm(A))
#  set_blocksize(rowis, get_blocksize(A))
  return rowis, colis
end


"""
  Construct ISLocalToGlobalMappings for the the rows and columns of the matrix
"""
function local_to_global_mapping(A::PetscMat)

  # localIS creates strided index sets, which require only constant
  # memory 
  if get_blocksize(A) == 1
    rowis, colis = localIS(A)
  else
    rowis, colis = localIS_block(A)
  end
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
function Base.similar{T}(a::PetscMat{T})
  p = Ref{C.Mat{T}}()
  chk(C.MatDuplicate(a.p, C.MAT_DO_NOT_COPY_VALUES, p))
  Mat{T}(p[])
end

Base.similar{T}(a::PetscMat{T}, ::Type{T}) = similar(a)

# TODO: make T2 -> Type{T2}
function Base.similar{T}(a::PetscMat{T}, T2::Type)
  MType = gettype(a)
  return Mat(T2, size(a)..., comm=comm(a), mtype=MType)
end

function Base.similar{T}(a::PetscMat{T}, T2::Type, m::Integer, n::Integer)
  MType = gettype(a)
  (m,n) == size(a) && T2==T ? similar(a) : Mat(T2, m,n, comm=comm(a), mtype=MType)
end

Base.similar{T}(a::PetscMat{T}, m::Integer, n::Integer) = similar(a, T, m, n)
Base.similar(a::PetscMat, T::Type, d::Dims) = similar(a, T, d...)
Base.similar{T}(a::PetscMat{T}, d::Dims) = similar(a, T, d)

function Base.copy{T}(a::PetscMat{T})
  p = Ref{C.Mat{T}}()
  chk(C.MatDuplicate(a.p, C.MAT_COPY_VALUES, p))
  Mat{T}(p[])
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
  Check if the matrix is assembled.  Whether all processes assembly state 
  is checked or only the local process is determined by `x.verify_assembled`.

  `local_only` forces only the local process to be checked, regardless of 
  `x.verify_assembled`.
"""
function isassembled(x::PetscMat; local_only=false)
  val = isassembled(x.p)
  if !local_only
    if x.verify_assembled
      val = MPI.Allreduce(Int8(val), MPI.LAND, comm(x))
    end
  end

  return Bool(val)
end

"""
  This function provides a mechanism for efficiently inserting values into
  and then assembling Petsc matrices and vectors.  The function f must be a 
  zero argument function.

  This function can be used with the do block syntax.
"""
function assemble(f::Function, x::Union{Vec,PetscMat},
  insertmode=x.insertmode,
  assemblytype::C.MatAssemblyType=C.MAT_FINAL_ASSEMBLY)
  if x.insertmode != insertmode && !isassembled(x, local_only=true)
    error("nested assemble with different insertmodes not allowed")
  end
  old_insertmode = x.insertmode
  try  # what is the purpose of the try - finally?
    x.insertmode = insertmode
    result = f()
    return result
  finally
    AssemblyBegin(x, assemblytype)
    yield() # do async computations while messages are in transit
    AssemblyEnd(x, assemblytype)
    x.insertmode = old_insertmode
  end
end

# force assembly even if it might not be necessary
"""
  Assemble the Petsc object
"""
function assemble(x::AbstractArray, t::C.MatAssemblyType=C.MAT_FINAL_ASSEMBLY)
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

iassemble(x::AbstractArray) = assemble(x, C.MAT_FLUSH_ASSEMBLY)

# like x[i,j] = v, but requires i,j to be 0-based indices for Petsc
function setindex0!{T}(x::Mat{T}, v::Array{T},
                       i::Array{PetscInt}, j::Array{PetscInt})
  ni = length(i)
  nj = length(j)
  if length(v) != ni*nj
    throw(ArgumentError("length(values) != length(indices)"))
  end
  chk(C.MatSetValues(x.p, ni, i, nj, j, v, x.insertmode))
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
# zero based indexing
#TODO: in 0.5, use boundscheck macro to verify stride=1

# global, non-block
function set_values!{T <: Scalar}(x::Mat{T}, idxm::StridedVecOrMat{PetscInt}, idxn::StridedVecOrMat{PetscInt}, v::StridedVecOrMat{T}, o::C.InsertMode=x.insertmode)
  # v should be m x n

  chk(C.MatSetValues(x.p, length(idxm), idxm, length(idxn), idxn, v, o))
end

function set_values!{T <: Scalar, I1 <: Integer, I2 <: Integer}(x::Mat{T}, idxm::StridedVecOrMat{I1}, idxn::StridedVecOrMat{I2}, v::StridedVecOrMat{T}, o::C.InsertMode=x.insertmode)

  _idxm = PetscInt[ i for i in idxm]
  _idxn = PetscInt[ i for i in idxn]
  set_values!(x, _idxm, _idxn, v, o)
end

function set_values!(x::Matrix, idxm::AbstractArray, idxn::AbstractArray, v::AbstractArray, o::C.InsertMode=C.INSERT_VALUES)
  if o == C.INSERT_VALUES
    for col = 1:length(idxn)
      colidx = idxn[col]
      for row = 1:length(idxm)
        rowidx = idxm[row]
        x[rowidx+1, colidx+1] = v[row, col]
      end
    end
  else
    for col = 1:length(idxn)
      colidx = idxn[col]
      for row = 1:length(idxm)
        rowidx = idxm[row]
        x[rowidx+1, colidx+1] += v[row, col]
      end
    end
  end

end


# global, block
function set_values_blocked!{T <: Scalar}(x::Mat{T}, idxm::StridedVecOrMat{PetscInt}, idxn::StridedVecOrMat{PetscInt}, v::StridedVecOrMat{T}, o::C.InsertMode=x.insertmode)

  # vals should be m*bs x n*bs
  chk(C.MatSetValuesBlocked(x.p, length(idxm), idxm, length(idxn), idxn, v, o))
end

function set_values_blocked!{T <: Scalar, I1 <: Integer, I2 <: Integer}(x::Mat{T}, idxm::StridedVecOrMat{I1}, idxn::StridedVecOrMat{I2}, v::StridedVecOrMat{T}, o::C.InsertMode=x.insertmode)

  _idxm = PetscInt[ i for i in idxm]
  _idxn = PetscInt[ i for i in idxn]
  set_values_blocked!(x, _idxm, _idxn, v, o)
end

# local, non-block
function set_values_local!{T <: Scalar}(x::Mat{T}, idxm::StridedVecOrMat{PetscInt}, idxn::StridedVecOrMat{PetscInt}, v::StridedVecOrMat{T}, o::C.InsertMode=x.insertmode)

  chk(C.MatSetValuesLocal(x.p, length(idxm), idxm, length(idxn), idxn, v, o))
end

function set_values_local!{T <: Scalar, I1 <: Integer, I2 <: Integer}(x::Mat{T}, idxm::StridedVecOrMat{I1}, idxn::StridedVecOrMat{I2}, v::StridedVecOrMat{T}, o::C.InsertMode=x.insertmode)

  _idxm = PetscInt[ i for i in idxm]
  _idxn = PetscInt[ i for i in idxn]
  set_values_local!(x, idxm, idxn, v, o)
end

function set_values_local!(x::Matrix, idxm::AbstractArray, idxn::AbstractArray, v::AbstractArray, o::C.InsertMode=C.INSERT_VALUES)

  set_values!(x, idxm, idxn, v, o)
end

# local, block
function set_values_blocked_local!{T <: Scalar}(x::Mat{T}, idxm::StridedVecOrMat{PetscInt}, idxn::StridedVecOrMat{PetscInt}, v::StridedVecOrMat{T}, o::C.InsertMode=x.insertmode)

  chk(C.MatSetValuesBlockedLocal(x.p, length(idxm), idxm, length(idxn), idxn, v, o))
end

function set_values_blocked_local!{T <: Scalar, I1 <: Integer, I2 <: Integer}(x::Mat{T}, idxm::StridedVecOrMat{I1}, idxn::StridedVecOrMat{I2}, v::StridedVecOrMat{T}, o::C.InsertMode=x.insertmode)

  _idxm = PetscInt[ i for i in idxm]
  _idxn = PetscInt[ i for i in idxn]
  set_values_blocked_local!(x, _idxm, _idxn, v, o)
end


###############################################################################
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

#=
"""
 create a new matrix wrapping around A for matrix-vector multiplication
 but which does not actually require new storage
"""
=#
for (f,pf) in ((:MatTranspose,:MatCreateTranspose), # acts like A.'
  (:MatNormal, :MatCreateNormal))      # acts like A'*A
  pfe = Expr(:quote, pf)
  @eval function $f{T}(a::PetscMat{T})
    p = Ref{C.Mat{T}}()
    chk(C.$pf(a.p, p))
    Mat{T}(p[], a)
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

    function Base.$f{T}(a::PetscMat{T})
      p = Ref{C.Mat{T}}()
      chk(C.$pf(a.p, C.MAT_INITIAL_MATRIX, p))
      Mat{T}(p[])
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



function (.*){T}(A::PetscMat{T}, x::Number)
  Y = copy(A)
  chk(C.MatScale(Y.p, T(x)))
  return Y
end
(.*){T}(x::Number, A::PetscMat{T}) = A .* x
(./){T}(A::PetscMat{T}, x::Number) = A .* inv(x)
(.\){T}(x::Number, A::PetscMat{T}) = A .* inv(x)

function (*){T}(A::PetscMat{T}, B::PetscMat{T})
  p = Ptr{Float64}(0)
  p_arr = Ref{C.Mat{T}}()
  chk(C.MatMatMult(A.p, B.p, C.MAT_INITIAL_MATRIX, real(T)(C.PETSC_DEFAULT), p_arr))
  new_mat = Mat{T}(p_arr[])
  return new_mat
end

#these two only work for SEQAIJ
function At_mul_B{T}(A::PetscMat{T}, B::PetscMat{T})
  p = Ptr{Float64}(0)
  p_arr = Ref{C.Mat{T}}()
  chk(C.MatTransposeMatMult(A.p, B.p, C.MAT_INITIAL_MATRIX, real(T)(C.PETSC_DEFAULT), p_arr))
  new_mat = Mat{T}(p_arr[])
  return new_mat
end

function A_mul_Bt{T}(A::PetscMat{T}, B::PetscMat{T})
  p = Ptr{Float64}(0)
  p_arr = Ref{C.Mat{T}}()
  chk(C.MatMatTransposeMult(A.p, B.p, C.MAT_INITIAL_MATRIX, real(T)(C.PETSC_DEFAULT), p_arr))
  new_mat = Mat{T}(p_arr[])
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

(-){T}(A::PetscMat{T}) = A * (-1)

# there don't appear to be PETSc functions for pointwise
# operations on matrices

function (==){T}(A::PetscMat{T}, b::PetscMat{T})
  bool_arr = Ref{PetscBool}()
  chk(C.MatEqual(A.p, b.p, bool_arr))
  return bool_arr[] != 0
end
#=
# needed for disambiguation
function (==){T}(A::PetscMat{T}, B::PetscMat{T})
  bool_arr = Ref{PetscBool}()
  chk(C.MatEqual(A.p, b.p, bool_arr))
  return bool_arr[] != 0
end
=#
"""
  Equality test for AbstractArray and PetscMat, SEQ only
"""
function (==){T}(A::PetscMat{T}, B::AbstractArray)
  if MPI.Comm_size(comm(A)) != 1
    throw(ArgumentError("Mat must reside on a single MPI rank"))
  end
  if size(A) != size(B)
    throw(ArgumentError("Matrices must be same size"))
  end

  m, n = size(B)

  isame = true
  for i=1:n
    for j=1:m
      isame = isame && A[i, j] == B[i,j]
    end
  end

  return isame
end

###############################################################################
# MatRow: accessing the structure of each row of a matrix
"""
  Object that enables access to a single row of a sparse matrix.

  Users *must* call restore when done with a MatRow, before attempting to 
  create another one.
"""
immutable MatRow{T}
  mat::C.Mat{T}  # the matrix to which the rows belong
  row::Int
  ncols::Int
  cols_ptr::Ptr{PetscInt}
  vals_ptr::Ptr{T}
#=
  function MatRow(A::Mat{T}, row::Integer, ref_ncols::Ref{PetscInt}, ref_cols::Ref{Ptr{PetscInt}}, ref_vals::Ref{Ptr{T}})
    ncols = ref_ncols[]
    cols = pointer_to_array(ref_cols[], ncols)
    vals = pointer_to_array(ref_vals[], ncols)

    obj = new(A, row, ref_ncols, ref_cols, ref_vals, ncols, cols, vals)
    finalizer(obj, restore)

    return obj
  end
=#
end

"""
  Preferred constructor for a MatRow
"""
function MatRow{T}(A::Mat{T}, row::Integer)

  ref_ncols = Ref{PetscInt}()
  ref_cols = Ref{Ptr{PetscInt}}()
  ref_vals = Ref{Ptr{T}}()
  chk(C.MatGetRow(A.p, row-1, ref_ncols, ref_cols, ref_vals))
  return MatRow{T}(A.p, row, ref_ncols[], ref_cols[], ref_vals[])
end

"""
  Count the number of non-zeros in a row of the matrix.  This temporarily
  creates a MatRow object, so it cannot be used if one already exists
"""
function count_row_nz{T}(A::Mat{T}, row::Integer)
  ref_ncols = Ref{PetscInt}()
  ref_cols = Ref{Ptr{PetscInt}}(C_NULL)
  ref_vals = Ref{Ptr{T}}(C_NULL)
  chk(C.MatGetRow(A.p, row-1, ref_ncols, ref_cols, ref_vals))
  ncols = ref_ncols[]
  chk(C.MatRestoreRow(A.p, row-1, ref_ncols, ref_cols, ref_vals))
  
  return ncols
end
function restore{T}(row::MatRow{T})
#  if !PetscFinalized(T) && !isfinalized(row.mat)
    C.MatRestoreRow(row.mat, row.row-1, Ref(PetscInt(row.ncols)), Ref{Ptr{PetscInt}}(C_NULL), Ref{Ptr{T}}(C_NULL))

#    return nothing  # return type stability
#  end
end

### indexing on a MatRow ###
import Base: length, size
length(A::MatRow) = A.ncols
size(A::MatRow) = (A.ncols,)
getcol(A::MatRow, i) = unsafe_load(A.cols_ptr, i) + 1
getval(A::MatRow, i) = unsafe_load(A.vals_ptr, i)


import Base.kron
"""
  Kronecker product for SEQ matrices only.  The output is a non-block matrix
  even if the inputs are block
"""
function kron{T}(A::Mat{T}, B::Mat{T})
  if (A.p == B.p)
    throw(ArgumentError("A and B cannot be same matrix"))
  end

  if MPI.Comm_size(comm(A)) != 1 || MPI.Comm_size(comm(B)) != 1
    throw(ArgumentError("A and B must reside on a single MPI rank"))
  end

  Am, An = size(A)
  Bm, Bn = size(B)
  # step 1: figure out size, sparsity pattern of result
  A_nz = zeros(Int, Am)
  B_nz = zeros(Int, Bm)
  for i=1:Am
    A_nz[i] = count_row_nz(A, i)
  end
  for i=1:Bm
    B_nz[i] = count_row_nz(B, i)
  end

  Dm = Am*Bm
  Dn = An*Bn
  D_nz = getDnzs(A_nz, B_nz)
  # create matrix
  # can't use C becaue that is the module name
  D = Mat(T, Dm, Dn, nnz=D_nz, mtype=C.MATSEQAIJ)

  # step 2: populate C, one row at a time
  max_entries = maximum(D_nz)
  D_colidx = zeros(PetscInt, max_entries)
  D_vals = zeros(T, max_entries)
  D_rowidx = zeros(PetscInt, 1)

  for i=1:Am
    rowA = MatRow(A, i)
    for j=1:Bm
      rowB = MatRow(B, j)
      rowidx = (i-1)*Bm + j

      pos = 1  # position in C_vals, C_colidx
      for k=1:length(rowA)
        Aval = getval(rowA, k)
        Aidx = getcol(rowA, k)
        for p=1:length(rowB)
          Bval = getval(rowB, p)
          Bidx = getcol(rowB, p)
          D_colidx[pos] = (Aidx - 1)*Bn + Bidx - 1
          D_vals[pos] = Aval*Bval
          pos += 1
        end
      end

      D_rowidx[1] = rowidx - 1
      idx_extract = sub(D_colidx, 1:(pos-1))
      vals_extract = sub(D_vals, 1:(pos-1))
      set_values!(D, D_rowidx, idx_extract, vals_extract)
      restore(rowB)
    end
    restore(rowA)
  end

  return D
end

function getDnzs(A_nz, B_nz)
# this function computes the number of non-zeros in D = kron(A, B), where A_nz and
# B_nz are the number of non-zeros in each row of A and B, respectively

  Am = length(A_nz)
  Bm = length(B_nz)

  Dm = length(A_nz)*length(B_nz)
  D_nz = zeros(PetscInt, Dm)
  for i=1:Am
    A_nz_i = A_nz[i]
    for j=1:Bm
      pos = (i-1)*Bm + j
      D_nz[pos] = A_nz_i*B_nz[j]
    end
  end

  return D_nz
end

function getmax_nz_col(A::SparseMatrixCSC)
  Am, An = size(A)
  A_maxnz = 0
  for i=1:An
    val = A.colptr[i+1] - A.colptr[i]
    if val > A_maxnz
      A_maxnz = val
    end
  end

  return A_maxnz
end

"""
  Kronecker product of A and B where the result is a Petsc Mat
"""
function PetscKron{T <: Scalar}(A::SparseMatrixCSC{T}, B::SparseMatrixCSC{T})

  Am = size(A, 1); An = size(A, 2)
  Bm = size(B, 1); Bn = size(B, 2)

  A_nz = zeros(Int, Am)
  B_nz = zeros(Int, Bm)

  for i in A.rowval
    A_nz[i] += 1
  end
  for i in B.rowval
    B_nz[i] += 1
  end

  Dm = Am*Bm
  Dn = An*Bn
  D_nz = getDnzs(A_nz, B_nz)
  # create matrix
  # can't use C becaue that is the module name
  D = Mat(T, Dm, Dn, nnz=D_nz, mtype=C.MATSEQAIJ)
  # because the Mat constructor is type unstable, use a function barrier
  PetscKron(A, B, D)

  return D
end

"""
  Kronecker product of A and B, storing the result in D.  D should already be 
  pre-allocated with the right sparsity pattern
"""
@noinline function PetscKron{T <: Scalar}(A::SparseMatrixCSC{T}, B::SparseMatrixCSC{T}, D::Mat{T})

  Am = size(A, 1); An = size(A, 2)
  Bm = size(B, 1); Bn = size(B, 2)

  # now figure out the maximum number of non-zeros in any column of D
  A_maxnz = getmax_nz_col(A)
  B_maxnz = getmax_nz_col(B)
  max_entries = A_maxnz * B_maxnz

  D_rowidx = zeros(PetscInt, max_entries)
  D_vals = zeros(T, max_entries)
  D_colidx = zeros(PetscInt, 1)

  # loop over columns of A and B
  for i=1:An
    col_start = A.colptr[i]
    col_end = A.colptr[i+1] - 1
    A_rowidx = unsafe_view(A.rowval, col_start:col_end)
    A_vals = unsafe_view(A.nzval, col_start:col_end)
    for j=1:Bn
      col_start = B.colptr[j]
      col_end = B.colptr[j+1] - 1
      B_rowidx = unsafe_view(B.rowval, col_start:col_end)
      B_vals = unsafe_view(B.nzval, col_start:col_end)

      pos = 1
      for k=1:length(A_vals)
        Aval = A_vals[k]
        Aidx = A_rowidx[k]
        for p=1:length(B_vals)
          Bval = B_vals[p]
          Bidx = B_rowidx[p]

          D_vals[pos] = Aval*Bval
          D_rowidx[pos] = (Aidx-1)*Bm + Bidx - 1
          pos += 1
        end
      end

      D_colidx[1] = (i-1)*Bn + j - 1
      rowidx_extract = unsafe_view(D_rowidx, 1:(pos-1))
      vals_extract = unsafe_view(D_vals, 1:(pos-1))

#      chk(C.MatSetValues(D.p, length(rowidx_extract), rowidx_extract, length(D_colidx), D_colidx, vals_extract, C.INSERT_VALUES))
      set_values!(D, rowidx_extract, D_colidx, vals_extract)
    end
#    assemble(D, C.MAT_FLUSH_ASSEMBLY)
  end

  return D
end

 
