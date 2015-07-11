export PetscMat, PetscMatSetType, PetscSetUp, PetscMatSetValues, PetscMatAssemblyBegin, PetscMatAssemblyEnd, PetscMatSetSizes, PetscMatGetSize


type PetscMat <: PetscObject
  pobj::Ptr{Void}
  function PetscMat(comm::MPI_Comm)
#    comm = PETSC_COMM_SELF();
    vec = Array(Ptr{Void},1)
    err = ccall( (:MatCreate,  libpetsclocation),Int32,(Int64,Ptr{Ptr{Void}}),comm.val,vec);
    vec = new(vec[1])
#    finalizer(vec,PetscDestroy)
    # does not seem to be called immediately when vec is no longer visible, is it called later during garbage collection?
    return vec
  end
end

  function PetscDestroy(vec::PetscMat)
    if (vec.pobj != 0)
      err = ccall( (:MatDestroy,  libpetsclocation),Int32,(Ptr{Ptr{Void}},), &vec.pobj);
    end
#    vec.pobj = 0

    println("Petsc Mat Destroy called")
#    sleep(5)
  end

  function PetscMatSetType(vec::PetscMat,name)
    err = ccall( (:MatSetType,  libpetsclocation),Int32,(Ptr{Void}, Cstring), vec.pobj,name);
  end

  function PetscSetUp(vec::PetscMat)
    err = ccall( ( :MatSetUp,  libpetsclocation),Int32,(Ptr{Void},), vec.pobj);
  end

#=
  PETSC_MAT_FLUSH_ASSEMBLY = 1;
  PETSC_MAT_FINAL_ASSEMBLY = 0
=#

  function PetscMatSetValues(vec::PetscMat,idi::Array{Int64},idj::Array{Int64},array::Array{Float64},flag::Integer)
    idi = idi - 1
    idj = idj - 1
    err = ccall( ( :MatSetValues,  libpetsclocation), Int32,(Ptr{Void},Int32,Ptr{Int32},Int32,Ptr{Int32},Ptr{Float64},Int32), vec.pobj,length(idi),convert(Array{Int32},idi),length(idi),convert(Array{Int32},idj),array,flag);
    idi = idi + 1
    idj = idj + 1
    return err
  end

  function PetscMatAssemblyBegin(obj::PetscMat,flg::Integer)
    err = ccall( ( :MatAssemblyBegin,  libpetsclocation),Int32,(Int64,Int32), obj.pobj,flg);
  end
  function PetscMatAssemblyBegin(obj::PetscMat)
    return PetscMatAssemblyBegin(obj,PETSC_MAT_FINAL_ASSEMBLY);
  end

  function PetscMatAssemblyEnd(obj::PetscMat,flg::Integer)
    err = ccall( ( :MatAssemblyEnd,  libpetsclocation),Int32,(Int64,Int32), obj.pobj,flg);
  end
  function PetscMatAssemblyEnd(obj::PetscMat)
    return PetscMatAssemblyEnd(obj,PETSC_MAT_FINAL_ASSEMBLY);
  end

  function PetscMatSetSizes(vec::PetscMat,m::Int,n::Int,M::Int,N::Int)
    err = ccall( ( :MatSetSizes,  libpetsclocation),Int32,(Ptr{Void},Int32,Int32,Int32,Int32), vec.pobj,m,n,M,N);
  end

  function PetscView(obj::PetscMat,viewer)
    err = ccall( (:MatView,  libpetsclocation),Int32,(Int64,Int64),obj.pobj,0);
  end

  function PetscMatGetSize(obj::PetscMat)
    m = Array(Int32,1)
    n = Array(Int32,1)
    err = ccall(Libdl.dlsym(libpetsc, :MatGetSize),Int32,(Int64,Ptr{Int32},Ptr{Int32}), obj.pobj,m,n);
    return (m[1],n[1])
  end
