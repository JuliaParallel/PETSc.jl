export PetscMat, PetscMatSetType, PetscSetUp, PetscMatSetValues, PetscMatAssemblyBegin, PetscMatAssemblyEnd, PetscMatSetSizes, PetscMatGetSize, PetscMatGetValues, PetscMatGetOwnershipRange, PetscMatXAIJSetPreallocation, PetscMatMPIAIJSetPreallocation, PetscMatSetFromOptions, PetscMatGetInfo


type PetscMat <: PetscObject
  pobj::Ptr{Void}
  function PetscMat(comm::MPI_Comm)
#    comm = PETSC_COMM_SELF();
    vec = Array(Ptr{Void},1)
    err = ccall( (:MatCreate,  libpetsclocation), PetscErrorCode, (comm_type, Ptr{Ptr{Void}}),comm.val,vec);
    vec = new(vec[1])
#    finalizer(vec,PetscDestroy)
    # does not seem to be called immediately when vec is no longer visible, is it called later during garbage collection?
    return vec
  end
end

  function PetscDestroy(vec::PetscMat)
    if (vec.pobj != 0)
      err = ccall( (:MatDestroy,  libpetsclocation), PetscErrorCode, (Ptr{Ptr{Void}},), &vec.pobj);
    end
#    vec.pobj = 0

#    println("Petsc Mat Destroy called")
#    sleep(5)
  end


immutable PetscMatInfo
    block_size::PetscLogDouble
    nz_allocated::PetscLogDouble
    nz_used::PetscLogDouble
    nz_unneeded::PetscLogDouble
    memory::PetscLogDouble
    assemblies::PetscLogDouble
    mallocs::PetscLogDouble
    fill_ratio_given::PetscLogDouble
    fill_ratio_needed::PetscLogDouble
    factor_mallocs::PetscLogDouble

    function PetscMatInfo()  # incomplete initialization
      return new()
    end
end


function show(io::IO, obj::PetscMatInfo)
# print the fields of PetscMatInfo
#  println("PetscMatInfo:")
  println(io, "  block_size : ", obj.block_size)
  println(io, "  nz_allocated : ", obj.nz_allocated)
  println(io, "  nz_used : ", obj.nz_used)
  println(io, "  nz_unneeded : ", obj.nz_unneeded)
  println(io, "  memory : ", obj.memory)
  println(io, "  assemblies : ", obj.assemblies)
  println(io, "  mallocs : ", obj.mallocs)
  println(io, "  fill_ratio_given : ", obj.fill_ratio_given)
  println(io, "  fill_ratio_needed : ", obj.fill_ratio_needed)
  println(io, "  factor_mallocs : ", obj.factor_mallocs)
end

  function PetscMatSetFromOptions(mat::PetscMat)
    ccall((:MatSetFromOptions,petsc),PetscErrorCode,(Ptr{Void},), mat.pobj)
  end


  function PetscMatSetType(vec::PetscMat,name)
    err = ccall( (:MatSetType,  libpetsclocation), PetscErrorCode,(Ptr{Void}, Cstring), vec.pobj,name);
  end

  function PetscSetUp(vec::PetscMat)
    err = ccall( ( :MatSetUp,  libpetsclocation), PetscErrorCode, (Ptr{Void},), vec.pobj);
  end

#=
  PETSC_MAT_FLUSH_ASSEMBLY = 1;
  PETSC_MAT_FINAL_ASSEMBLY = 0
=#

  function PetscMatSetValues(vec::PetscMat,idi::Array{PetscInt},idj::Array{PetscInt},array::Array{PetscScalar},flag::Integer)
    # remember, only matrices can be inserted into a Petsc matrix
    # if array is a 3 by 3, then idi and idj are vectors of length 3
    idi = idi
    idj = idj

    @assert length(idi)*length(idj) == length(array)

    # do check here to ensure array is the right shape (remember tranpose)
    err = ccall( ( :MatSetValues,  libpetsclocation), PetscErrorCode, (Ptr{Void}, PetscInt, Ptr{PetscInt}, PetscInt, Ptr{PetscInt}, Ptr{PetscScalar},Int32), vec.pobj,length(idi), idi, length(idi), idj,array,flag);
    idi = idi
    idj = idj
    return err
  end

  function PetscMatAssemblyBegin(obj::PetscMat,flg::Integer)
    err = ccall( ( :MatAssemblyBegin,  libpetsclocation), PetscErrorCode,(Ptr{Void},Int32), obj.pobj,flg);
  end

  function PetscMatAssemblyBegin(obj::PetscMat)
    return PetscMatAssemblyBegin(obj,PETSC_MAT_FINAL_ASSEMBLY);
  end

  function PetscMatAssemblyEnd(obj::PetscMat,flg::Integer)
    err = ccall( ( :MatAssemblyEnd,  libpetsclocation), PetscErrorCode,(Ptr{Void},Int32), obj.pobj,flg);
  end

  function PetscMatAssemblyEnd(obj::PetscMat)
    return PetscMatAssemblyEnd(obj,PETSC_MAT_FINAL_ASSEMBLY);
  end

  function PetscMatSetSizes(vec::PetscMat,m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt)
    err = ccall( ( :MatSetSizes,  libpetsclocation), PetscErrorCode, (Ptr{Void},PetscInt, PetscInt, PetscInt, PetscInt), vec.pobj,m,n,M,N);
  end

  function PetscView(obj::PetscMat,viewer)
    err = ccall( (:MatView,  libpetsclocation), PetscErrorCode, (Ptr{Void}, Int64),obj.pobj,0);
  end

  function PetscMatGetSize(obj::PetscMat)
    m = Array(PetscInt, 1)
    n = Array(PetscInt, 1)
    err = ccall(Libdl.dlsym(libpetsc, :MatGetSize), PetscErrorCode,(Ptr{Void}, Ptr{PetscInt},Ptr{PetscInt}), obj.pobj,m,n);
    return m[1],n[1]
  end


export PetscMatGetLocalSize
function PetscMatGetLocalSize(mat::PetscMat)
    m = Array(PetscInt, 1)
    n = Array(PetscInt, 1)
    ccall((:MatGetLocalSize,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscInt},Ptr{PetscInt}),mat.pobj, m, n)
    return m[1], n[1]
end




# Petsc populates v row major, so Julia will see the returned array as being transposed
function PetscMatGetValues(obj::PetscMat, idxm::Array{PetscInt, 1}, idxn::Array{PetscInt, 1}, v::Array{PetscScalar, 2})
    # do check here to ensure v is the right shape

    ccall((:MatGetValues,petsc),PetscErrorCode,(Ptr{Void},PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{PetscScalar}), obj.pobj, length(idxm), idxm, length(idxn), idxn, v)

end

function PetscMatGetOwnershipRange(mat::PetscMat)
    low = Array(PetscInt,1)
    high = Array(PetscInt,1)
    ccall((:MatGetOwnershipRange,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscInt},Ptr{PetscInt}),mat.pobj, low, high)

    return low[1], high[1]
end

### new function ###
export PetscMatAXPY, PetscMatAYPX, PetscMatScale, PetscMatShift
function PetscMatAXPY(Y::PetscMat, a::PetscScalar, X::PetscMat, str::PetscMatStructure)
    ccall((:MatAXPY,petsc),PetscErrorCode,(Ptr{Void},PetscScalar,Ptr{Void},PetscMatStructure), Y.pobj, a, X.pobj, str)
end

function PetscMatAYPX( Y::PetscMat, a::PetscScalar, X::PetscMat, str::PetscMatStructure)
    ccall((:MatAYPX,petsc),PetscErrorCode,(Ptr{Void},PetscScalar,Ptr{Void}, PetscMatStructure), Y.pobj, a, X.pobj, str)
end

function PetscMatScale(mat::PetscMat, a::PetscScalar)
    ccall((:MatScale,petsc),PetscErrorCode,(Ptr{Void},PetscScalar), mat.pobj, a)
end


function PetscMatShift(mat::PetscMat, a::PetscScalar)
    ccall((:MatShift,petsc),PetscErrorCode,(Ptr{Void},PetscScalar), mat.pobj, a)
end

export PetscMatMult, PetscMatMultAdd, PetscMatMultTranspose, PetscMatMultHermitianTranspose

function PetscMatMult(mat::PetscMat, x::PetscVec, y::PetscVec)
    ccall((:MatMult,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void}), mat.pobj, x.pobj, y.pobj)
end


function PetscMatMultAdd(A::PetscMat, x::PetscVec, y::PetscVec, z::PetscVec)
    ccall((:MatMultAdd,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void},Ptr{Void}), A.pobj, x.pobj, y.pobj, z.pobj)
end

function PetscMatMultTranspose(A::PetscMat, x::PetscVec, y::PetscVec)
    ccall((:MatMultTranspose,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void}), A.pobj, x.pobj, y.pobj)
end

function PetscMatMultHermitianTranspose(mat::PetscMat, x::PetscVec, y::PetscVec)
    ccall((:MatMultHermitianTranspose,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void}),mat.pobj, x.pobj, y.pobj)
end

function PetscMatXAIJSetPreallocation(mat::PetscMat, bs::PetscInt, dnnz::AbstractArray{PetscInt, 1}, onnz::AbstractArray{PetscInt,1}, dnnzu::AbstractArray{PetscInt, 1}, onnzu::AbstractArray{PetscInt, 1})
# this is a unified interface for matrix preallocation for the Petsc built in
# matrix types: Aij, Bij, and their respective symmetric forms SAij, SBij
# for a non symmetric format matrix, dnnzu and onnzu are not required
# and vice versa for a non symmetric format matrix

    ccall((:MatXAIJSetPreallocation,petsc),PetscErrorCode,(Ptr{Void},PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}), mat.pobj, bs, dnnz, onnz, dnnzu, onnzu)
end


arr_or_null = Union(AbstractArray{PetscInt}, Ptr{Void})

function PetscMatMPIAIJSetPreallocation(mat::PetscMat, d_nz::PetscInt, d_nnz::arr_or_null, o_nz::PetscInt, o_nnz::arr_or_null)

    ccall((:MatMPIAIJSetPreallocation,petsc),PetscErrorCode,(Ptr{Void}, PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt}), mat.pobj, d_nz, d_nnz, o_nz, o_nnz)
end

function PetscMatGetInfo(mat::PetscMat, info_type::Int32)
    matinfo = PetscMatInfo()  # create uninitialized struct
    ref_matinfo = Ref{PetscMatInfo}(matinfo)
    ccall((:MatGetInfo,petsc),PetscErrorCode,(Ptr{Void}, Int32,Ref{PetscMatInfo}), mat.pobj ,info_type, ref_matinfo)
    return ref_matinfo[]
end


