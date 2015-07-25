export PetscVec, PetscVecSetType, PetscVecSetValues, PetscVecAssemblyBegin, PetscVecAssemblyEnd, PetscVecSetSizes, PetscVecGetSize, PetscVecNorm, PetscVecGetValues


type PetscVec <: PetscObject
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


