export PC, KSPGetPC, PCSetType, PCGetType, PCFactorSetUseInPlace, PCFactorGetUseInPlace, PCSetReusePreconditioner, PCGetReusePreconditioner, PCFactorSetAllowDiagonalFill, PCFactorGetAllowDiagonalFill, PCFactorSetLevels, PCFactorGetLevels, PCSetReusePreconditioner, PCGetReusePreconditioner, PCBJacobiGetSubKSP, PCFactorSetFill, PCJacobiSetType, PCJacobiGetType

# preconditioner contex
# the KSP object creates the PC contex, so we don't provide a constructor
type PC
  pobj::Ptr{Void}
end

function KSPGetPC(ksp::KSP)
    arr = Array(Ptr{Void}, 1)
    ccall((:KSPGetPC,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void}),ksp.pobj, arr)
    return PC(arr[1])
end

function PCSetType(pc::PC, pctype::PCType)
    ccall((:PCSetType,petsc),PetscErrorCode,(Ptr{Void},Cstring), pc.pobj, pctype)
end

function PCGetType(pc::PC)
    arr = Array(Ptr{UInt8}, 1)
    ccall((:PCGetType,petsc),PetscErrorCode,(Ptr{Void},Ptr{Ptr{UInt8}}), pc.pobj, arr)
    return bytestring(arr[1])
end

function PCFactorSetUseInPlace(pc::PC, arg2::PetscBool)
    ccall((:PCFactorSetUseInPlace,petsc),PetscErrorCode,(Ptr{Void},PetscBool),pc.pobj, arg2)
end

function PCFactorGetUseInPlace(pc::PC)
    arr = Array(PetscBool, 1)
    ccall((:PCFactorGetUseInPlace,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscBool}), pc.pobj, arr)
    return arr[1]
end


function PCSetReusePreconditioner(pc::PC,arg2::PetscBool)
    ccall((:PCSetReusePreconditioner,petsc),PetscErrorCode,(Ptr{Void},PetscBool), pc.pobj, arg2)
end


function PCGetReusePreconditioner(pc::PC)
    arr = Array(PetscBool, 1)
    ccall((:PCGetReusePreconditioner,petsc),PetscErrorCode,(Ptr{Void}, Ptr{PetscBool}), pc.pobj, arr)
    return arr[1]
end

function PCFactorSetAllowDiagonalFill(pc::PC,arg2::PetscBool)
    ccall((:PCFactorSetAllowDiagonalFill,petsc),PetscErrorCode,(Ptr{Void},PetscBool), pc.pobj, arg2)
end

function PCFactorGetAllowDiagonalFill(pc::PC)
   arr = Array(PetscBool, 1)
    ccall((:PCFactorGetAllowDiagonalFill,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscBool}), pc.pobj, arr)
    return arr[1]
end

function PCFactorSetLevels(pc::PC,arg2::PetscInt)
    ccall((:PCFactorSetLevels,petsc),PetscErrorCode,(Ptr{Void}, PetscInt), pc.pobj, arg2)
end

function PCFactorGetLevels(pc::PC)
    arr = Array(PetscInt, 1)
    ccall((:PCFactorGetLevels,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscInt}),pc.pobj, arr)
    return arr[1]
end


function PCSetReusePreconditioner(pc::PC, arg2::PetscBool)
    ccall((:PCSetReusePreconditioner,petsc),PetscErrorCode,(Ptr{Void}, PetscBool), pc.pobj, arg2)
end

function PCGetReusePreconditioner(pc::PC)
    arr = Array(PetscBool, 1)
    ccall((:PCGetReusePreconditioner,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscBool}), pc.pobj, arr)
    return arr[1]
end

function PCBJacobiGetSubKSP(pc::PC)
    n_local_arr = Array(PetscInt, 1)
    first_local = Array(PetscInt, 1)
    ksp_ptrarr = Array(Ptr{Ptr{Void}}, 1)
    ccall((:PCBJacobiGetSubKSP,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{Ptr{Void}}}), pc.pobj, n_local_arr, first_local, ksp_ptrarr)

    n_local = n_local_arr[1]
    ksp_ptrarr2 = pointer_to_array(ksp_ptrarr[1], n_local)

    ksp_arr = Array(KSP, n_local)
    for i=1:n_local
      ksp_arr[i] = KSP(ksp_ptrarr2[i])
    end

    return n_local, first_local[1], ksp_arr
end

function PCFactorSetFill(pc::PC, fill::PetscReal)
    ccall((:PCFactorSetFill,petsc),PetscErrorCode,(Ptr{Void}, PetscReal), pc.pobj, fill)
end

function PCJacobiSetType(pc::PC, jacobitype::PCJacobiType)
    ccall((:PCJacobiSetType,petsc),PetscErrorCode,(Ptr{Void}, PCJacobiType), pc.pobj, jacobitype)
end

function PCJacobiGetType(pc::PC)
    arr = Array(PCJacobiType, 1)
    ccall((:PCJacobiGetType,petsc),PetscErrorCode,(Ptr{Void},Ptr{PCJacobiType}),pc.pobj, arr)
    return arr[1]
end


