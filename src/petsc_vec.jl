export PetscVec, PetscVecSetType, PetscVecSetValues, PetscVecAssemblyBegin, PetscVecAssemblyEnd, PetscVecSetSizes, PetscVecGetSize, PetscVecNorm, PetscVecgetValues, PetscVecGetValues


type PetscVec <: PetscObject
  pobj::Ptr{Void}
  function PetscVec(comm::MPI_Comm)
#    comm = PETSC_COMM_SELF();
    vec = Array(Ptr{Void},1)
    err = ccall(( :VecCreate, libpetsclocation ),Int32,(Int64,Ptr{Void}),comm.val,vec);
    vec = new(vec[1])
#    finalizer(vec,PetscDestroy)
    # does not seem to be called immediately when vec is no longer visible, is it called later during garbage collection?
    return vec
  end
end

  function PetscDestroy(vec::PetscVec)
    if (vec.pobj != 0)
      err = ccall(( :VecDestroy, libpetsclocation),Int32,(Ptr{Ptr{Void}},), &vec.pobj);
    end
#    vec.pobj = 0  # unnecessary? vec no longer has any references to it
    println("VecDestroy called")
#    sleep(5)
#    println("finished sleeping")
  end

  function PetscVecSetType(vec::PetscVec,name)
    err = ccall((:VecSetType,  libpetsclocation),Int32,(Int64,Cstring), vec.pobj,name);
  end

  function PetscVec(array::Array{Float64})
    vec = PetscVec()
    err = ccall(( :VecSetType,  libpetsclocation),Int32,(Ptr{Void},Cstring), vec.pobj,"seq");
    err = ccall( (:VecSetSizes,  libpetsclocation),Int32,(Ptr{Void},Int32,Int32), vec.pobj,length(array),length(array));
    # want a 32 bit int array so build it ourselves
    idx = Array(Int32,length(array));
    for i=1:length(array);  idx[i] = i-1;  end
    err = ccall( ( :VecSetValues,  libpetsclocation), Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Float64},Int32), vec.pobj,length(idx),idx,array,PETSC_INSERT_VALUES);
    err = ccall( ( :VecAssemblyBegin,  libpetsclocation),Int32,(Ptr{Void},), vec.pobj);
    err = ccall( ( :VecAssemblyEnd,  libpetsclocation),Int32,(Ptr{Void},), vec.pobj);
    return vec
  end

  function PetscVecSetValues(vec::PetscVec,idx::Array{Int64},array::Array{Float64},flag::Integer)
    idx = idx - 1
    err = ccall( ( :VecSetValues,  libpetsclocation), Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Float64},Int32), vec.pobj,length(idx),convert(Array{Int32},idx),array,flag);
    idx = idx + 1
    return err
  end

  function PetscVecSetValues(vec::PetscVec,idx::Array{Int64},array::Array{Float64})
    PetscVecSetValues(vec,idx,array,PETSC_INSERT_VALUES)
  end

  function PetscVecSetValues(vec::PetscVec,array::Array{Float64})
    idx = Array(Int64,length(array))
    for i=1:length(array);  idx[i] = i-1;  end
    PetscVecSetValues(vec,idx,array,PETSC_INSERT_VALUES)
  end

  function PetscVecAssemblyBegin(obj::PetscVec)
    err = ccall( ( :VecAssemblyBegin,  libpetsclocation),Int32,(Ptr{Void},), obj.pobj);
  end

  function PetscVecAssemblyEnd(obj::PetscVec)
    err = ccall( ( :VecAssemblyEnd,  libpetsclocation),Int32,(Ptr{Void},), obj.pobj);
  end

  function PetscVecSetSizes(vec::PetscVec,n::Int,N::Int)
    err = ccall( ( :VecSetSizes,  libpetsclocation),Int32,(Ptr{Void},Int32,Int32), vec.pobj,n,N);
  end

  function PetscView(obj::PetscVec,viewer)
    err = ccall( ( :VecView,  libpetsclocation),Int32,(Ptr{Void},Int64),obj.pobj,0);
  end

  function PetscVecGetSize(obj::PetscVec)
    n = Array(Int32,1)
    err = ccall( ( :VecGetSize,  libpetsclocation),Int32,(Ptr{Void},Ptr{Int32}), obj.pobj,n);
    return n[1]
  end

  function PetscVecNorm(obj::PetscVec,normtype::Int)
    n = Array(Float64,1)
    err = ccall( ( :VecNorm,  libpetsclocation),Int32,(Ptr{Void},Int32,Ptr{Int32}), obj.pobj,normtype, n);
    return n[1]
  end
  function PetscVecNorm(obj::PetscVec)
    return PetscVecNorm(obj,PETSC_NORM_2)
  end

function PetscVecGetValues(vec::PetscVec, ni::Integer, ix::AbstractArray{Int,1}, y::AbstractArray{PetscScalar,1})

     # need indices to be PetscInt
     ix_local = Array(PetscInt, ni)
     for i=1:ni
       ix_local[i] = ix[i]
     end

    err = ccall((:VecGetValues,petsc),PetscErrorCode,(Ptr{Void},PetscInt,Ptr{PetscInt},Ptr{PetscScalar}),vec.pobj, ni, ix_local, y)

    return nothing
end


