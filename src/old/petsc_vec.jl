export PetscVec, PetscVecSetType, PetscVecSetValues, PetscVecAssemblyBegin, PetscVecAssemblyEnd, PetscVecSetSizes, PetscVecGetSize, PetscVecNorm, PetscVecGetValues, PetscVecGetOwnershipRange, PetscVecGetArray, PetscVecRestoreArray, PetscVecGetArrayRead, PetscVecRestoreArrayRead, PetscVecSet, PetscVecSqrtAbs, PetscVecLog, PetscVecExp, PetscVecAbs, PetscVecMax, PetscVecMin, PetscVecCopy, PetscVecDuplicate, PetscVecAXPY, PetscVecAXPBY, PetscVecAYPX, PetscVecWAXPY, PetscVecMAXPY, PetscVecAXPBYPCZ, PetscVecScale, PetscVecDot, PetscVecTDot, PetscVecSum, PetscVecSwap, PetscVecReciprocal, PetscVecShift, PetscVecPointwiseMult, PetscVecPointwiseDivide


type PetscVec
  pobj::Ptr{Void}
  function PetscVec(comm::MPI_Comm)
#    comm = PETSC_COMM_SELF();
    vec = Array(Ptr{Void},1)
    err = ccall(( :VecCreate, libpetsclocation ), PetscErrorCode,(comm_type,Ptr{Void}),comm.val,vec);
    vec = new(vec[1])
#    finalizer(vec,PetscDestroy)
    # does not seem to be called immediately when vec is no longer visible, is it called later during garbage collection? - yes
    return vec

  end

  function PetscVec(pobj::Ptr{Void})  # default constructor
    return new(pobj)
  end
 
end

  function PetscDestroy(vec::PetscVec)
    if (vec.pobj != 0)
      err = ccall(( :VecDestroy, libpetsclocation), PetscErrorCode, (Ptr{Ptr{Void}},), &vec.pobj);
    end
    vec.pobj = 0  # unnecessary? vec no longer has any references to it
#    println("VecDestroy called")
  end

  function PetscVecSetType(vec::PetscVec,name)
    err = ccall((:VecSetType,  libpetsclocation), PetscErrorCode, (comm_type, Cstring), vec.pobj,name);
  end

  function PetscVec(array::Array{PetscScalar})
    vec = PetscVec()
    err = ccall(( :VecSetType,  libpetsclocation), PetscErrorCode,(Ptr{Void},Cstring), vec.pobj,"seq");
    err = ccall( (:VecSetSizes,  libpetsclocation), PetscErrorCode, (Ptr{Void}, PetscInt, PetscInt), vec.pobj,length(array),length(array));
    # want a PetscInt array so build it ourselves
    idx = Array(PetscInt, length(array));
    for i=1:length(array);  idx[i] = i-1;  end
    err = ccall( ( :VecSetValues,  libpetsclocation), PetscErrorCode,(Ptr{Void},PetscInt, Ptr{PetscInt},Ptr{PetscScalar},Int32), vec.pobj,length(idx),idx,array,PETSC_INSERT_VALUES);
    err = ccall( ( :VecAssemblyBegin,  libpetsclocation), PetscErrorCode,(Ptr{Void},), vec.pobj);
    err = ccall( ( :VecAssemblyEnd,  libpetsclocation), PetscErrorCode, (Ptr{Void},), vec.pobj);
    return vec
  end

  function PetscVecSetValues(vec::PetscVec,idx::Array{PetscInt},array::Array{PetscScalar},flag::Integer)

    err = ccall( ( :VecSetValues,  libpetsclocation), PetscErrorCode, (Ptr{Void},PetscInt ,Ptr{PetscInt},Ptr{PetscScalar},Int32), vec.pobj,length(idx), idx,array,flag);
    return err
  end

  function PetscVecSetValues(vec::PetscVec,idx::Array{PetscInt},array::Array{PetscScalar})
    PetscVecSetValues(vec,idx,array,PETSC_INSERT_VALUES)
  end

  function PetscVecSetValues(vec::PetscVec,array::Array{PetscScalar})
    idx = Array(PetscInt,length(array))
    for i=1:length(array);  idx[i] = i-1;  end
    PetscVecSetValues(vec,idx,array,PETSC_INSERT_VALUES)
  end

  function PetscVecAssemblyBegin(obj::PetscVec)
    err = ccall( ( :VecAssemblyBegin,  libpetsclocation), PetscErrorCode, (Ptr{Void},), obj.pobj);
  end

  function PetscVecAssemblyEnd(obj::PetscVec)
    err = ccall( ( :VecAssemblyEnd,  libpetsclocation), PetscErrorCode,(Ptr{Void},), obj.pobj);
  end

  function PetscVecSetSizes(vec::PetscVec,n::PetscInt, N::PetscInt)
    err = ccall( ( :VecSetSizes,  libpetsclocation), PetscErrorCode, (Ptr{Void}, PetscInt, PetscInt), vec.pobj,n,N);
  end

  function PetscView(obj::PetscVec,viewer)
    err = ccall( ( :VecView,  libpetsclocation), PetscErrorCode, (Ptr{Void},Int64),obj.pobj,0);
  end

  function PetscVecGetSize(obj::PetscVec)
    n = Array(PetscInt, 1)
    err = ccall( ( :VecGetSize,  libpetsclocation), PetscErrorCode, (Ptr{Void},Ptr{PetscInt}), obj.pobj,n);
    return n[1]
  end

  function PetscVecNorm(obj::PetscVec,normtype::Integer)
    n = Array(PetscReal,1)
    err = ccall( ( :VecNorm,  libpetsclocation), PetscScalar, (Ptr{Void},Int32,Ptr{PetscReal}), obj.pobj,normtype, n);
    return n[1]
  end

  function PetscVecNorm(obj::PetscVec)
    return PetscVecNorm(obj,PETSC_NORM_2)
  end

function PetscVecGetValues(vec::PetscVec, ni::Integer, ix::AbstractArray{PetscInt,1}, y::AbstractArray{PetscScalar,1})

     # need indices to be PetscInt
#     ix_local = Array(PetscInt, ni)
#     for i=1:ni
#       ix_local[i] = ix[i]
#     end

    err = ccall((:VecGetValues,petsc),PetscErrorCode,(Ptr{Void},PetscInt,Ptr{PetscInt},Ptr{PetscScalar}),vec.pobj, ni, ix, y)

    return nothing
end

function PetscVecGetOwnershipRange(vec::PetscVec)
    low = Array(PetscInt, 1)
    high = Array(PetscInt, 1)

    ccall((:VecGetOwnershipRange,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscInt},Ptr{PetscInt}),vec.pobj, low, high)

  return low[1], high[1]
end


# new functions

function PetscVecGetArray(vec::PetscVec)
# gets a pointer to the data underlying a Petsc vec, turns it into a Julia
# array
# ptr_arr must be passed into PetscVecRestoreArray
    ptr_arr = Array(Ptr{PetscScalar}, 1)
    ccall((:VecGetArray,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Ptr{PetscScalar}}),vec.pobj, ptr_arr)

    first, last = PetscVecGetOwnershipRange(vec)
    len = last - first
    arr = pointer_to_array(ptr_arr[1], len)
    return arr, ptr_arr
end


function PetscVecRestoreArray(vec::PetscVec, ptr_arr::Array{Ptr{PetscScalar}, 1})
    ccall((:VecRestoreArray,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Ptr{PetscScalar}}),vec.pobj, ptr_arr)
end


function PetscVecGetArrayRead(vec::PetscVec)
    ptr_arr = Array(Ptr{PetscScalar}, 1)
    ccall((:VecGetArrayRead,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Ptr{PetscScalar}}),vec.pobj, ptr_arr)

    first, last = PetscVecGetOwnershipRange(vec)
    len = last - first
    arr = pointer_to_array(ptr_arr[1], len)
    return arr, ptr_arr
end

function PetscVecRestoreArrayRead(vec::PetscVec, ptr_arr::Array{Ptr{PetscScalar}, 1})

    ccall((:VecRestoreArrayRead,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Ptr{PetscScalar}}),vec.pobj, ptr_arr)
#    ccall((:VecRestoreArrayRead,petsc),PetscErrorCode,(Vec,Ptr{Ptr{PetscScalar}}),arg1,arg2)
end



function PetscVecSet(vec::PetscVec, val::PetscScalar)
    ccall((:VecSet,petsc),PetscErrorCode,(Ptr{Void},PetscScalar), vec.pobj, val)
end


function PetscVecSqrtAbs(vec::PetscVec)
    ccall((:VecSqrtAbs,petsc),PetscErrorCode,(Ptr{Void},), vec.pobj)
end

function PetscVecLog(vec::PetscVec)
    ccall((:VecLog,petsc),PetscErrorCode,(Ptr{Void},),vec.pobj)
end

function PetscVecExp(vec::PetscVec)
    ccall((:VecExp,petsc),PetscErrorCode,(Ptr{Void},),vec.pobj)
end

function PetscVecAbs(vec::PetscVec)
    ccall((:VecAbs,petsc),PetscErrorCode,(Ptr{Void},),vec.pobj)
end


function PetscVecMax(vec::PetscVec)
    r = Array(PetscReal, 1) # max value
    idx = Array(PetscInt, 1)  # index of max value
    ccall((:VecMax,petsc),PetscErrorCode,(Ptr{Void}, Ptr{PetscInt},Ptr{PetscReal}), vec.pobj, idx, r)

    return r[1], idx[1]
end

function PetscVecMin(vec::PetscVec)
    r = Array(PetscReal, 1) # min value
    idx = Array(PetscInt, 1)  # index of min value
 
    ccall((:VecMin,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscInt},Ptr{PetscReal}), vec.pobj, idx, r)

    return r[1], idx[1]
end

function PetscVecReciprocal(vec::PetscVec)
    ccall((:VecReciprocal,petsc),PetscErrorCode,(Ptr{Void},),vec.pobj)
end

function PetscVecShift(vec::PetscVec, a::PetscScalar)
    ccall((:VecShift,petsc),PetscErrorCode,(Ptr{Void},PetscScalar), vec.pobj, a)
end


function PetscVecPointwiseMult(w::PetscVec, x::PetscVec,y::PetscVec)
    ccall((:VecPointwiseMult,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Void}, Ptr{Void}), w.pobj, x.pobj, y.pobj)
end

function PetscVecPointwiseDivide(w::PetscVec, x::PetscVec, y::PetscVec)
    ccall((:VecPointwiseDivide,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Void}, Ptr{Void}), w.pobj, x.pobj, y.pobj)
end






function PetscVecCopy(vec::PetscVec , vec2::PetscVec)
    ccall((:VecCopy,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Void}), vec.pobj, vec2.pobj)
end

function PetscVecDuplicate( vec::PetscVec)
    ptr_arr = Array(Ptr{Void}, 1)

    ccall((:VecDuplicate,petsc),PetscErrorCode,( Ptr{Void}, Ptr{Ptr{Void}}), vec.pobj, ptr_arr)

    return PetscVec(ptr_arr[1])
end



# Some vector linear algebra
function PetscVecAXPY( vec1::PetscVec, a::PetscScalar, vec2::PetscVec)
    ccall((:VecAXPY,petsc),PetscErrorCode,(Ptr{Void},PetscScalar,Ptr{Void}), vec1.pobj, a, vec2.pobj)
end

function PetscVecAXPBY(vec1::PetscVec, a::PetscScalar, b::PetscScalar, vec2::PetscVec)
    ccall((:VecAXPBY,petsc),PetscErrorCode,(Ptr{Void},PetscScalar,PetscScalar,Ptr{Void}),vec1.pobj, a, b, vec2.pobj)
end

function PetscVecMAXPY(vec1::PetscVec, n::PetscInt, a::AbstractArray{PetscScalar, 1}, x::AbstractArray{Ptr{Void}, 1})
# the vector x must contains the pointers from the PetscVec objects, not the PetscVec objects themselves

    ccall((:VecMAXPY,petsc),PetscErrorCode,(Ptr{Void},PetscInt,Ptr{PetscScalar},Ptr{Ptr{Void}}),vec1.pobj, n, a, x)
end

function PetscVecAYPX(vec1::PetscVec, a::PetscScalar, vec2::PetscVec)
    ccall((:VecAYPX,petsc),PetscErrorCode,(Ptr{Void},PetscScalar,Ptr{Void}),vec1.pobj, a, vec2.pobj)
end

function PetscVecWAXPY(w::PetscVec, a::PetscScalar, x::PetscVec, y::PetscVec)
    ccall((:VecWAXPY,petsc),PetscErrorCode,(Ptr{Void},PetscScalar,Ptr{Void},Ptr{Void}), w.pobj, a, x.pobj, y.pobj)
end

function PetscVecAXPBYPCZ(z::PetscVec, alpha::PetscScalar, beta::PetscScalar, gamma::PetscScalar, x::PetscVec, y::PetscVec)
    ccall((:VecAXPBYPCZ,petsc),PetscErrorCode,(Ptr{Void},PetscScalar,PetscScalar,PetscScalar,Ptr{Void},Ptr{Void}), z.pobj, alpha,beta, gamma, x.pobj, y.pobj)
end


function PetscVecScale(vec::PetscVec, a::PetscScalar)
    ccall((:VecScale,petsc),PetscErrorCode,(Ptr{Void},PetscScalar), vec.pobj, a)
end

function PetscVecDot(x::PetscVec, y::PetscVec)
    r = Array(PetscScalar, 1)
    err = ccall((:VecDot,petsc),PetscErrorCode,( Ptr{Void}, Ptr{Void}, Ptr{PetscScalar}), x.pobj, y.pobj, r)

    println("PetscVecDot error code = ", err)
    return r[1]
end

function PetscVecTDot(x::PetscVec, y::PetscVec)
    r = Array(PetscScalar, 1)
    ccall((:VecTDot,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Void},Ptr{PetscScalar}), x.pobj, y.pobj, r)

    return r[1]
end


function PetscVecSum(vec::PetscVec)
    r = Array(PetscScalar, 1)
    ccall((:VecSum,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscScalar}),vec.pobj, r)
    return r[1]
end


function PetscVecSwap(x::PetscVec, y::PetscVec)
    ccall((:VecSwap,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Void}), x.pobj, y.pobj)

end


