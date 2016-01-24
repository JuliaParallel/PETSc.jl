# these functions work directly with the C interface, offering better
# performance than the high level interface, particularly for indexing

##### Vec #####

  function VecCreate{T}(::Type{T}; comm=MPI.COMM_WORLD)
    vref = Ref{Vec{T}}()
    chk(VecCreate(comm, vref))
    return vref[]
  end 


  function SetValues{T}(vec::Vec{T},idx::AbstractVector{PetscInt},
                                vals::AbstractVector{T},
                                flag::Integer=INSERT_VALUES)

    chk(VecSetValues(vec, length(idx), idx, vals, InsertMode(flag)))
  end

  
  function AssemblyBegin(obj::Vec)

    chk(VecAssemblyBegin(obj))
  end

  function AssemblyEnd(obj::Vec)
    chk(VecAssemblyEnd(obj))
  end

  function GetValues{T}(vec::Vec{T}, idx::AbstractArray{PetscInt,1}, 
                             y::AbstractArray{T,1})

    chk(VecGetValues(vec, length(idx), idx, y))

  end


##### Mat #####

function MatCreateShell{T}(arg1::MPI.Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer, arg6::Ptr{Void}, dtype::Type{T}=Float64)
  # arg6 is the user provided context
    arg7 = Ref{Mat{dtype}}()
    chk(MatCreateShell(arg1, arg2, arg3, arg5, arg6, arg6, arg7))

    return Mat(arg7[])
end

#= # this function signature is not distinct from the auto generated one
function MatShellSetOperation(arg1::Mat,arg2::MatOperation,arg3::Ptr{Void})
# arg3 is a function pointer, and must have the signature:
# void fname(Mat, vec, vec) for MATOP_MULT
    chk(MatShellSetOperation(arg1, arg2, arg3))
end
=#

# TODO: make this work for non Float64
function MatShellGetContext(arg1::Mat{Float64})
# get the user provided context for the matrix shell
    arg2 = Ref{Ptr{Void}}()
    chk(ccall((:MatShellGetContext,petscRealDouble),PetscErrorCode,(Mat{Float64},Ref{Ptr{Void}}),arg1.pobj,arg2))
    return arg2[]  # turn it into a julia object here?
end

  function SetValues{ST}(vec::Mat,idi::AbstractArray{PetscInt},idj::AbstractArray{PetscInt},array::AbstractArray{ST},flag::Integer=INSERT_VALUES)
    # remember, only matrices can be inserted into a Petsc matrix
    # if array is a 3 by 3, then idi and idj are vectors of length 3

#    @assert length(idi)*length(idj) == length(array)

    # do check here to ensure array is the right shape (remember tranpose)
    chk(MatSetValues(vec, length(idi), idi, length(idj), idj, array, InsertMode(flag)))

  end

  function SetValuesBlocked{ST}(mat::Mat, idi::AbstractArray{PetscInt}, idj::AbstractArray{PetscInt}, v::AbstractArray{ST}, flag::Integer=INSERT_VALUES)

    chk(MatSetValuesBlocked(mat, length(idi), idi, length(idj), idj, v, InsertMode(flag)))
  end

  function MatSetOption(mat::Mat,arg2::MatOption,arg3::Bool)
    chk(MatSetOption(mat, arg2, PetscBool(arg3)))
  end

  function AssemblyBegin(obj::Mat,flg=MAT_FINAL_ASSEMBLY)
    chk(MatAssemblyBegin(obj, MatAssemblyType(flg)))
  end

  function AssemblyEnd(obj::Mat,flg=MAT_FINAL_ASSEMBLY)
    chk(MatAssemblyEnd(obj, MatAssemblyType(flg)))
  end


  function GetValues{ST}(obj::Mat, idxm::AbstractArray{PetscInt, 1}, idxn::AbstractArray{PetscInt, 1}, v::AbstractArray{ST})
    # do check here to ensure v is the right shape
    chk(MatGetValues(obj, length(idxm), idxm, length(idxn), idxn, v))
end






