# Julia wrapper for header: /home/jared/build/petsc-3.6.0/include/petscksp.h
# Automatically generated using Clang.jl wrap_c, version 0.0.0


function PetscIsInfOrNanReal(arg0::Type{Float64})
    err = ccall((:PetscIsInfOrNanReal,petsc1),PetscErrorCode,())
    return err
end

function PetscIsNormalReal(arg0::Type{Float64})
    err = ccall((:PetscIsNormalReal,petsc1),PetscBool,())
    return err
end

function PetscSetHelpVersionFunctions(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscSetHelpVersionFunctions,petsc1),PetscErrorCode,(Ptr{Void},Ptr{Void}),arg1,arg2)
    return err
end

function PetscCommDuplicate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{MPI_Comm},StridedArray{MPI_Comm},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:PetscCommDuplicate,petsc1),PetscErrorCode,(comm_type,Ptr{comm_type},Ptr{Cint}),arg1.val,arg2,arg3)
    return err
end

function PetscCommDestroy(arg0::Type{Float64},arg1::Union(Ptr{MPI_Comm},StridedArray{MPI_Comm},Ptr{Void}))
    err = ccall((:PetscCommDestroy,petsc1),PetscErrorCode,(Ptr{comm_type},),arg1)
    return err
end

function PetscMallocSet(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscMallocSet,petsc1),PetscErrorCode,(Ptr{Void},Ptr{Void}),arg1,arg2)
    return err
end

function PetscMallocClear(arg0::Type{Float64})
    err = ccall((:PetscMallocClear,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function PetscMallocDump(arg0::Type{Float64},arg1::Union(Ptr{FILE},StridedArray{FILE},Ptr{Void}))
    ccall((:PetscMallocDump,petsc1),PetscErrorCode,(Ptr{FILE},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscMallocDumpLog(arg0::Type{Float64},arg1::Union(Ptr{FILE},StridedArray{FILE},Ptr{Void}))
    ccall((:PetscMallocDumpLog,petsc1),PetscErrorCode,(Ptr{FILE},),arg1)
end 
=#
function PetscMallocGetCurrentUsage(arg0::Type{Float64},arg1::Union(Ptr{PetscLogDouble},StridedArray{PetscLogDouble},Ptr{Void}))
    err = ccall((:PetscMallocGetCurrentUsage,petsc1),PetscErrorCode,(Ptr{PetscLogDouble},),arg1)
    return err
end

function PetscMallocGetMaximumUsage(arg0::Type{Float64},arg1::Union(Ptr{PetscLogDouble},StridedArray{PetscLogDouble},Ptr{Void}))
    err = ccall((:PetscMallocGetMaximumUsage,petsc1),PetscErrorCode,(Ptr{PetscLogDouble},),arg1)
    return err
end

function PetscMallocDebug(arg0::Type{Float64},arg1::PetscBool)
    err = ccall((:PetscMallocDebug,petsc1),PetscErrorCode,(PetscBool,),arg1)
    return err
end

function PetscMallocGetDebug(arg0::Type{Float64},arg1::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscMallocGetDebug,petsc1),PetscErrorCode,(Ptr{PetscBool},),arg1)
    return err
end

function PetscMallocValidate(arg0::Type{Float64},arg1::Integer,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol))
    err = ccall((:PetscMallocValidate,petsc1),PetscErrorCode,(Cint,Cstring,Cstring),arg1,arg2,arg3)
    return err
end

function PetscMallocSetDumpLog(arg0::Type{Float64})
    err = ccall((:PetscMallocSetDumpLog,petsc1),PetscErrorCode,())
    return err
end

function PetscMallocSetDumpLogThreshold(arg0::Type{Float64},arg1::PetscLogDouble)
    err = ccall((:PetscMallocSetDumpLogThreshold,petsc1),PetscErrorCode,(PetscLogDouble,),arg1)
    return err
end

function PetscMallocGetDumpLog(arg0::Type{Float64},arg1::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscMallocGetDumpLog,petsc1),PetscErrorCode,(Ptr{PetscBool},),arg1)
    return err
end

#= skipping function with undefined symbols: 
 function PetscDataTypeToMPIDataType(arg0::Type{Float64},arg1::PetscDataType,arg2::Union(Ptr{MPI_Datatype},StridedArray{MPI_Datatype},Ptr{Void}))
    ccall((:PetscDataTypeToMPIDataType,petsc1),PetscErrorCode,(PetscDataType,Ptr{MPI_Datatype}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscMPIDataTypeToPetscDataType(arg0::Type{Float64},arg1::MPI_Datatype,arg2::Union(Ptr{PetscDataType},StridedArray{PetscDataType},Ptr{Void}))
    ccall((:PetscMPIDataTypeToPetscDataType,petsc1),PetscErrorCode,(MPI_Datatype,Ptr{PetscDataType}),arg1,arg2)
end 
=#
function PetscDataTypeGetSize(arg0::Type{Float64},arg1::PetscDataType,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:PetscDataTypeGetSize,petsc1),PetscErrorCode,(PetscDataType,Ptr{Cint}),arg1,arg2)
    return err
end

function PetscDataTypeFromString(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{PetscDataType},StridedArray{PetscDataType},Ptr{Void}),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscDataTypeFromString,petsc1),PetscErrorCode,(Cstring,Ptr{PetscDataType},Ptr{PetscBool}),arg1,arg2,arg3)
    return err
end

function PetscBitMemcpy(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Integer,arg5::Integer,arg6::PetscDataType)
    err = ccall((:PetscBitMemcpy,petsc1),PetscErrorCode,(Ptr{Void},PetscInt,Ptr{Void},PetscInt,PetscInt,PetscDataType),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function PetscMemmove(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),size_t::Integer)
    err = ccall((:PetscMemmove,petsc1),PetscErrorCode,(Ptr{Void},Ptr{Void},Cint),arg1,arg2,size_t)
    return err
end

function PetscMemcmp(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),size_t::Integer,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscMemcmp,petsc1),PetscErrorCode,(Ptr{Void},Ptr{Void},Cint,Ptr{PetscBool}),arg1,arg2,size_t,arg3)
    return err
end

function PetscStrlen(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:PetscStrlen,petsc1),PetscErrorCode,(Cstring,Ptr{Cint}),arg1,arg2)
    return err
end

function PetscStrToArray(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Uint8,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Ptr{Ptr{Uint8}}},StridedArray{Ptr{Ptr{Uint8}}},Ptr{Void}))
    err = ccall((:PetscStrToArray,petsc1),PetscErrorCode,(Cstring,Uint8,Ptr{Cint},Ptr{Ptr{Ptr{Uint8}}}),arg1,arg2,arg3,arg4)
    return err
end

function PetscStrToArrayDestroy(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:PetscStrToArrayDestroy,petsc1),PetscErrorCode,(Cint,Ptr{Ptr{Uint8}}),arg1,arg2)
    return err
end

function PetscStrcmp(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscStrcmp,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{PetscBool}),arg1,arg2,arg3)
    return err
end

function PetscStrgrt(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscStrgrt,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{PetscBool}),arg1,arg2,arg3)
    return err
end

function PetscStrcasecmp(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscStrcasecmp,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{PetscBool}),arg1,arg2,arg3)
    return err
end

function PetscStrncmp(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),size_t::Integer,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscStrncmp,petsc1),PetscErrorCode,(Cstring,Cstring,Cint,Ptr{PetscBool}),arg1,arg2,size_t,arg3)
    return err
end

function PetscStrcpy(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol))
    err = ccall((:PetscStrcpy,petsc1),PetscErrorCode,(Cstring,Cstring),arg1,arg2)
    return err
end

function PetscStrcat(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol))
    err = ccall((:PetscStrcat,petsc1),PetscErrorCode,(Cstring,Cstring),arg1,arg2)
    return err
end

function PetscStrncat(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),size_t::Integer)
    err = ccall((:PetscStrncat,petsc1),PetscErrorCode,(Cstring,Cstring,Cint),arg1,arg2,size_t)
    return err
end

function PetscStrncpy(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),size_t::Integer)
    err = ccall((:PetscStrncpy,petsc1),PetscErrorCode,(Cstring,Cstring,Cint),arg1,arg2,size_t)
    return err
end

function PetscStrchr(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Uint8,arg3::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:PetscStrchr,petsc1),PetscErrorCode,(Cstring,Uint8,Ptr{Ptr{Uint8}}),arg1,arg2,arg3)
    return err
end

function PetscStrtolower(arg0::Type{Float64},arg1::Union(ByteString,Symbol))
    err = ccall((:PetscStrtolower,petsc1),PetscErrorCode,(Cstring,),arg1)
    return err
end

function PetscStrtoupper(arg0::Type{Float64},arg1::Union(ByteString,Symbol))
    err = ccall((:PetscStrtoupper,petsc1),PetscErrorCode,(Cstring,),arg1)
    return err
end

function PetscStrrchr(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Uint8,arg3::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:PetscStrrchr,petsc1),PetscErrorCode,(Cstring,Uint8,Ptr{Ptr{Uint8}}),arg1,arg2,arg3)
    return err
end

function PetscStrstr(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:PetscStrstr,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Ptr{Uint8}}),arg1,arg2,arg3)
    return err
end

function PetscStrrstr(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:PetscStrrstr,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Ptr{Uint8}}),arg1,arg2,arg3)
    return err
end

function PetscStrendswith(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscStrendswith,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{PetscBool}),arg1,arg2,arg3)
    return err
end

function PetscStrbeginswith(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscStrbeginswith,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{PetscBool}),arg1,arg2,arg3)
    return err
end

function PetscStrendswithwhich(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PetscStrendswithwhich,petsc1),PetscErrorCode,(Cstring,Ptr{Ptr{Uint8}},Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function PetscStrallocpy(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:PetscStrallocpy,petsc1),PetscErrorCode,(Cstring,Ptr{Ptr{Uint8}}),arg1,arg2)
    return err
end

function PetscStrArrayallocpy(arg0::Type{Float64},arg1::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}),arg2::Union(Ptr{Ptr{Ptr{Uint8}}},StridedArray{Ptr{Ptr{Uint8}}},Ptr{Void}))
    err = ccall((:PetscStrArrayallocpy,petsc1),PetscErrorCode,(Ptr{Ptr{Uint8}},Ptr{Ptr{Ptr{Uint8}}}),arg1,arg2)
    return err
end

function PetscStrArrayDestroy(arg0::Type{Float64},arg1::Union(Ptr{Ptr{Ptr{Uint8}}},StridedArray{Ptr{Ptr{Uint8}}},Ptr{Void}))
    err = ccall((:PetscStrArrayDestroy,petsc1),PetscErrorCode,(Ptr{Ptr{Ptr{Uint8}}},),arg1)
    return err
end

function PetscStrNArrayallocpy(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}),arg3::Union(Ptr{Ptr{Ptr{Uint8}}},StridedArray{Ptr{Ptr{Uint8}}},Ptr{Void}))
    err = ccall((:PetscStrNArrayallocpy,petsc1),PetscErrorCode,(PetscInt,Ptr{Ptr{Uint8}},Ptr{Ptr{Ptr{Uint8}}}),arg1,arg2,arg3)
    return err
end

function PetscStrNArrayDestroy(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{Ptr{Ptr{Uint8}}},StridedArray{Ptr{Ptr{Uint8}}},Ptr{Void}))
    err = ccall((:PetscStrNArrayDestroy,petsc1),PetscErrorCode,(PetscInt,Ptr{Ptr{Ptr{Uint8}}}),arg1,arg2)
    return err
end

function PetscStrreplace(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),size_t::Integer)
    err = ccall((:PetscStrreplace,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint),arg1.val,arg2,arg3,size_t)
    return err
end

function PetscStrcmpNoError(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscStrcmpNoError,petsc1),Void,(Cstring,Cstring,Ptr{PetscBool}),arg1,arg2,arg3)
    return err
end

#= skipping function with undefined symbols: 
 function PetscTokenCreate(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Uint8,arg3::Union(Ptr{PetscToken},StridedArray{PetscToken},Ptr{Void}))
    ccall((:PetscTokenCreate,petsc1),PetscErrorCode,(Cstring,Uint8,Ptr{PetscToken}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscTokenFind(arg0::Type{Float64},arg1::PetscToken,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:PetscTokenFind,petsc1),PetscErrorCode,(PetscToken,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscTokenDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscToken},StridedArray{PetscToken},Ptr{Void}))
    ccall((:PetscTokenDestroy,petsc1),PetscErrorCode,(Ptr{PetscToken},),arg1)
end 
=#
function PetscEListFind(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}),arg3::Union(ByteString,Symbol),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscEListFind,petsc1),PetscErrorCode,(PetscInt,Ptr{Ptr{Uint8}},Cstring,Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function PetscEnumFind(arg0::Type{Float64},arg1::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}),arg2::Union(ByteString,Symbol),arg3::Union(Ptr{PetscEnum},StridedArray{PetscEnum},Ptr{Void}),arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscEnumFind,petsc1),PetscErrorCode,(Ptr{Ptr{Uint8}},Cstring,Ptr{PetscEnum},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
    return err
end

function PetscMaxSum(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PetscMaxSum,petsc1),PetscErrorCode,(comm_type,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1.val,arg2,arg3,arg4)
    return err
end

#= skipping function with undefined symbols: 
 function MPIULong_Send(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::Integer,arg3::MPI_Datatype,arg4::PetscMPIInt,arg5::PetscMPIInt,arg6::MPI_Comm)
    ccall((:MPIULong_Send,petsc1),PetscErrorCode,(Ptr{Void},PetscInt,MPI_Datatype,PetscMPIInt,PetscMPIInt,comm_type),arg1,arg2,arg3,arg4,arg5,arg6.val)
end 
=#
#= skipping function with undefined symbols: 
 function MPIULong_Recv(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::Integer,arg3::MPI_Datatype,arg4::PetscMPIInt,arg5::PetscMPIInt,arg6::MPI_Comm)
    ccall((:MPIULong_Recv,petsc1),PetscErrorCode,(Ptr{Void},PetscInt,MPI_Datatype,PetscMPIInt,PetscMPIInt,comm_type),arg1,arg2,arg3,arg4,arg5,arg6.val)
end 
=#
function PetscErrorPrintfInitialize(arg0::Type{Float64})
    err = ccall((:PetscErrorPrintfInitialize,petsc1),PetscErrorCode,())
    return err
end

function PetscErrorMessage(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}),arg3::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:PetscErrorMessage,petsc1),PetscErrorCode,(Cint,Ptr{Ptr{Uint8}},Ptr{Ptr{Uint8}}),arg1,arg2,arg3)
    return err
end

function PetscTraceBackErrorHandler(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::PetscErrorCode,arg6::PetscErrorType,arg7::Union(ByteString,Symbol),arg8::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscTraceBackErrorHandler,petsc1),PetscErrorCode,(comm_type,Cint,Cstring,Cstring,PetscErrorCode,PetscErrorType,Cstring,Ptr{Void}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
    return err
end

function PetscIgnoreErrorHandler(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::PetscErrorCode,arg6::PetscErrorType,arg7::Union(ByteString,Symbol),arg8::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscIgnoreErrorHandler,petsc1),PetscErrorCode,(comm_type,Cint,Cstring,Cstring,PetscErrorCode,PetscErrorType,Cstring,Ptr{Void}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
    return err
end

function PetscEmacsClientErrorHandler(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::PetscErrorCode,arg6::PetscErrorType,arg7::Union(ByteString,Symbol),arg8::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscEmacsClientErrorHandler,petsc1),PetscErrorCode,(comm_type,Cint,Cstring,Cstring,PetscErrorCode,PetscErrorType,Cstring,Ptr{Void}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
    return err
end

function PetscMPIAbortErrorHandler(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::PetscErrorCode,arg6::PetscErrorType,arg7::Union(ByteString,Symbol),arg8::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscMPIAbortErrorHandler,petsc1),PetscErrorCode,(comm_type,Cint,Cstring,Cstring,PetscErrorCode,PetscErrorType,Cstring,Ptr{Void}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
    return err
end

function PetscAbortErrorHandler(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::PetscErrorCode,arg6::PetscErrorType,arg7::Union(ByteString,Symbol),arg8::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscAbortErrorHandler,petsc1),PetscErrorCode,(comm_type,Cint,Cstring,Cstring,PetscErrorCode,PetscErrorType,Cstring,Ptr{Void}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
    return err
end

function PetscAttachDebuggerErrorHandler(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::PetscErrorCode,arg6::PetscErrorType,arg7::Union(ByteString,Symbol),arg8::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscAttachDebuggerErrorHandler,petsc1),PetscErrorCode,(comm_type,Cint,Cstring,Cstring,PetscErrorCode,PetscErrorType,Cstring,Ptr{Void}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
    return err
end

function PetscReturnErrorHandler(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::PetscErrorCode,arg6::PetscErrorType,arg7::Union(ByteString,Symbol),arg8::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscReturnErrorHandler,petsc1),PetscErrorCode,(comm_type,Cint,Cstring,Cstring,PetscErrorCode,PetscErrorType,Cstring,Ptr{Void}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
    return err
end

function PetscPushErrorHandler(arg0::Type{Float64},handler::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscPushErrorHandler,petsc1),PetscErrorCode,(Ptr{Void},Ptr{Void}),handler,arg1)
    return err
end

function PetscPopErrorHandler(arg0::Type{Float64})
    err = ccall((:PetscPopErrorHandler,petsc1),PetscErrorCode,())
    return err
end

function PetscSignalHandlerDefault(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscSignalHandlerDefault,petsc1),PetscErrorCode,(Cint,Ptr{Void}),arg1,arg2)
    return err
end

function PetscPushSignalHandler(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscPushSignalHandler,petsc1),PetscErrorCode,(Ptr{Void},Ptr{Void}),arg1,arg2)
    return err
end

function PetscPopSignalHandler(arg0::Type{Float64})
    err = ccall((:PetscPopSignalHandler,petsc1),PetscErrorCode,())
    return err
end

function PetscCheckPointerSetIntensity(arg0::Type{Float64},arg1::Integer)
    err = ccall((:PetscCheckPointerSetIntensity,petsc1),PetscErrorCode,(PetscInt,),arg1)
    return err
end

function PetscSetFPTrap(arg0::Type{Float64},arg1::PetscFPTrap)
    err = ccall((:PetscSetFPTrap,petsc1),PetscErrorCode,(PetscFPTrap,),arg1)
    return err
end

function PetscFPTrapPush(arg0::Type{Float64},arg1::PetscFPTrap)
    err = ccall((:PetscFPTrapPush,petsc1),PetscErrorCode,(PetscFPTrap,),arg1)
    return err
end

function PetscFPTrapPop(arg0::Type{Float64})
    err = ccall((:PetscFPTrapPop,petsc1),PetscErrorCode,())
    return err
end

function PetscStackCopy(arg0::Type{Float64},arg1::Union(Ptr{PetscStack},StridedArray{PetscStack},Ptr{Void}),arg2::Union(Ptr{PetscStack},StridedArray{PetscStack},Ptr{Void}))
    err = ccall((:PetscStackCopy,petsc1),PetscErrorCode,(Ptr{PetscStack},Ptr{PetscStack}),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function PetscStackPrint(arg0::Type{Float64},arg1::Union(Ptr{PetscStack},StridedArray{PetscStack},Ptr{Void}),arg2::Union(Ptr{FILE},StridedArray{FILE},Ptr{Void}))
    ccall((:PetscStackPrint,petsc1),PetscErrorCode,(Ptr{PetscStack},Ptr{FILE}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscStackView(arg0::Type{Float64},arg1::Union(Ptr{FILE},StridedArray{FILE},Ptr{Void}))
    ccall((:PetscStackView,petsc1),PetscErrorCode,(Ptr{FILE},),arg1)
end 
=#
function PetscStackDestroy(arg0::Type{Float64})
    err = ccall((:PetscStackDestroy,petsc1),PetscErrorCode,())
    return err
end

function PetscClassIdRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{PetscClassId},StridedArray{PetscClassId},Ptr{Void}))
    err = ccall((:PetscClassIdRegister,petsc1),PetscErrorCode,(Cstring,Ptr{PetscClassId}),arg1,arg2)
    return err
end

function PetscMemoryGetCurrentUsage(arg0::Type{Float64},arg1::Union(Ptr{PetscLogDouble},StridedArray{PetscLogDouble},Ptr{Void}))
    err = ccall((:PetscMemoryGetCurrentUsage,petsc1),PetscErrorCode,(Ptr{PetscLogDouble},),arg1)
    return err
end

function PetscMemoryGetMaximumUsage(arg0::Type{Float64},arg1::Union(Ptr{PetscLogDouble},StridedArray{PetscLogDouble},Ptr{Void}))
    err = ccall((:PetscMemoryGetMaximumUsage,petsc1),PetscErrorCode,(Ptr{PetscLogDouble},),arg1)
    return err
end

function PetscMemorySetGetMaximumUsage(arg0::Type{Float64})
    err = ccall((:PetscMemorySetGetMaximumUsage,petsc1),PetscErrorCode,())
    return err
end

function PetscMemoryTrace(arg0::Type{Float64},arg1::Union(ByteString,Symbol))
    err = ccall((:PetscMemoryTrace,petsc1),PetscErrorCode,(Cstring,),arg1)
    return err
end

function PetscInfoAllow(arg0::Type{Float64},arg1::PetscBool,arg2::Union(ByteString,Symbol))
    err = ccall((:PetscInfoAllow,petsc1),PetscErrorCode,(PetscBool,Cstring),arg1,arg2)
    return err
end

function PetscSleep(arg0::Type{Float64})
    err = ccall((:PetscSleep,petsc1),PetscErrorCode,())
    return err
end

function PetscInitialize(arg0::Type{Float64},arg1::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg2::Union(Ptr{Ptr{Ptr{Uint8}}},StridedArray{Ptr{Ptr{Uint8}}},Ptr{Void}),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol))
    err = ccall((:PetscInitialize,petsc1),PetscErrorCode,(Ptr{Cint},Ptr{Ptr{Ptr{Uint8}}},Cstring,Cstring),arg1,arg2,arg3,arg4)
    return err
end

function PetscInitializeNoPointers(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol))
    err = ccall((:PetscInitializeNoPointers,petsc1),PetscErrorCode,(Cint,Ptr{Ptr{Uint8}},Cstring,Cstring),arg1,arg2,arg3,arg4)
    return err
end

function PetscInitializeNoArguments(arg0::Type{Float64})
    err = ccall((:PetscInitializeNoArguments,petsc1),PetscErrorCode,())
    return err
end

function PetscInitialized(arg0::Type{Float64},arg1::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscInitialized,petsc1),PetscErrorCode,(Ptr{PetscBool},),arg1)
    return err
end

function PetscFinalized(arg0::Type{Float64},arg1::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscFinalized,petsc1),PetscErrorCode,(Ptr{PetscBool},),arg1)
    return err
end

function PetscFinalize(arg0::Type{Float64})
    err = ccall((:PetscFinalize,petsc1),PetscErrorCode,())
    return err
end

function PetscInitializeFortran(arg0::Type{Float64})
    err = ccall((:PetscInitializeFortran,petsc1),PetscErrorCode,())
    return err
end

function PetscGetArgs(arg0::Type{Float64},arg1::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg2::Union(Ptr{Ptr{Ptr{Uint8}}},StridedArray{Ptr{Ptr{Uint8}}},Ptr{Void}))
    err = ccall((:PetscGetArgs,petsc1),PetscErrorCode,(Ptr{Cint},Ptr{Ptr{Ptr{Uint8}}}),arg1,arg2)
    return err
end

function PetscGetArguments(arg0::Type{Float64},arg1::Union(Ptr{Ptr{Ptr{Uint8}}},StridedArray{Ptr{Ptr{Uint8}}},Ptr{Void}))
    err = ccall((:PetscGetArguments,petsc1),PetscErrorCode,(Ptr{Ptr{Ptr{Uint8}}},),arg1)
    return err
end

function PetscFreeArguments(arg0::Type{Float64},arg1::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:PetscFreeArguments,petsc1),PetscErrorCode,(Ptr{Ptr{Uint8}},),arg1)
    return err
end

function PetscEnd(arg0::Type{Float64})
    err = ccall((:PetscEnd,petsc1),PetscErrorCode,())
    return err
end

function PetscSysInitializePackage(arg0::Type{Float64})
    err = ccall((:PetscSysInitializePackage,petsc1),PetscErrorCode,())
    return err
end

function PetscPythonInitialize(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol))
    err = ccall((:PetscPythonInitialize,petsc1),PetscErrorCode,(Cstring,Cstring),arg1,arg2)
    return err
end

function PetscPythonFinalize(arg0::Type{Float64})
    err = ccall((:PetscPythonFinalize,petsc1),PetscErrorCode,())
    return err
end

function PetscPythonPrintError(arg0::Type{Float64})
    err = ccall((:PetscPythonPrintError,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function PetscPythonMonitorSet(arg0::Type{Float64},arg1::PetscObject,arg2::Union(ByteString,Symbol))
    ccall((:PetscPythonMonitorSet,petsc1),PetscErrorCode,(PetscObject,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscObject},StridedArray{PetscObject},Ptr{Void}))
    ccall((:PetscObjectDestroy,petsc1),PetscErrorCode,(Ptr{PetscObject},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectGetComm(arg0::Type{Float64},arg1::PetscObject,arg2::Union(Ptr{MPI_Comm},StridedArray{MPI_Comm},Ptr{Void}))
    ccall((:PetscObjectGetComm,petsc1),PetscErrorCode,(PetscObject,Ptr{comm_type}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectGetClassId(arg0::Type{Float64},arg1::PetscObject,arg2::Union(Ptr{PetscClassId},StridedArray{PetscClassId},Ptr{Void}))
    ccall((:PetscObjectGetClassId,petsc1),PetscErrorCode,(PetscObject,Ptr{PetscClassId}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectGetClassName(arg0::Type{Float64},arg1::PetscObject,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:PetscObjectGetClassName,petsc1),PetscErrorCode,(PetscObject,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectSetType(arg0::Type{Float64},arg1::PetscObject,arg2::Union(ByteString,Symbol))
    ccall((:PetscObjectSetType,petsc1),PetscErrorCode,(PetscObject,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectSetPrecision(arg0::Type{Float64},arg1::PetscObject,arg2::PetscPrecision)
    ccall((:PetscObjectSetPrecision,petsc1),PetscErrorCode,(PetscObject,PetscPrecision),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectGetType(arg0::Type{Float64},arg1::PetscObject,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:PetscObjectGetType,petsc1),PetscErrorCode,(PetscObject,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectSetName(arg0::Type{Float64},arg1::PetscObject,arg2::Union(ByteString,Symbol))
    ccall((:PetscObjectSetName,petsc1),PetscErrorCode,(PetscObject,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectGetName(arg0::Type{Float64},arg1::PetscObject,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:PetscObjectGetName,petsc1),PetscErrorCode,(PetscObject,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectSetTabLevel(arg0::Type{Float64},arg1::PetscObject,arg2::Integer)
    ccall((:PetscObjectSetTabLevel,petsc1),PetscErrorCode,(PetscObject,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectGetTabLevel(arg0::Type{Float64},arg1::PetscObject,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscObjectGetTabLevel,petsc1),PetscErrorCode,(PetscObject,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectIncrementTabLevel(arg0::Type{Float64},arg1::PetscObject,arg2::PetscObject,arg3::Integer)
    ccall((:PetscObjectIncrementTabLevel,petsc1),PetscErrorCode,(PetscObject,PetscObject,PetscInt),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectReference(arg0::Type{Float64},arg1::PetscObject)
    ccall((:PetscObjectReference,petsc1),PetscErrorCode,(PetscObject,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectGetReference(arg0::Type{Float64},arg1::PetscObject,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscObjectGetReference,petsc1),PetscErrorCode,(PetscObject,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectDereference(arg0::Type{Float64},arg1::PetscObject)
    ccall((:PetscObjectDereference,petsc1),PetscErrorCode,(PetscObject,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectGetNewTag(arg0::Type{Float64},arg1::PetscObject,arg2::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}))
    ccall((:PetscObjectGetNewTag,petsc1),PetscErrorCode,(PetscObject,Ptr{PetscMPIInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectCompose(arg0::Type{Float64},arg1::PetscObject,arg2::Union(ByteString,Symbol),arg3::PetscObject)
    ccall((:PetscObjectCompose,petsc1),PetscErrorCode,(PetscObject,Cstring,PetscObject),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectRemoveReference(arg0::Type{Float64},arg1::PetscObject,arg2::Union(ByteString,Symbol))
    ccall((:PetscObjectRemoveReference,petsc1),PetscErrorCode,(PetscObject,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectQuery(arg0::Type{Float64},arg1::PetscObject,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{PetscObject},StridedArray{PetscObject},Ptr{Void}))
    ccall((:PetscObjectQuery,petsc1),PetscErrorCode,(PetscObject,Cstring,Ptr{PetscObject}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectComposeFunction_Private(arg0::Type{Float64},arg1::PetscObject,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscObjectComposeFunction_Private,petsc1),PetscErrorCode,(PetscObject,Cstring,Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectSetFromOptions(arg0::Type{Float64},arg1::PetscObject)
    ccall((:PetscObjectSetFromOptions,petsc1),PetscErrorCode,(PetscObject,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectSetUp(arg0::Type{Float64},arg1::PetscObject)
    ccall((:PetscObjectSetUp,petsc1),PetscErrorCode,(PetscObject,),arg1)
end 
=#
function PetscCommGetNewTag(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}))
    err = ccall((:PetscCommGetNewTag,petsc1),PetscErrorCode,(comm_type,Ptr{PetscMPIInt}),arg1.val,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function PetscObjectAddOptionsHandler(arg0::Type{Float64},arg1::PetscObject,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscObjectAddOptionsHandler,petsc1),PetscErrorCode,(PetscObject,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectProcessOptionsHandlers(arg0::Type{Float64},arg1::PetscObject)
    ccall((:PetscObjectProcessOptionsHandlers,petsc1),PetscErrorCode,(PetscObject,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectDestroyOptionsHandlers(arg0::Type{Float64},arg1::PetscObject)
    ccall((:PetscObjectDestroyOptionsHandlers,petsc1),PetscErrorCode,(PetscObject,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectsListGetGlobalNumbering(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{PetscObject},StridedArray{PetscObject},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscObjectsListGetGlobalNumbering,petsc1),PetscErrorCode,(comm_type,PetscInt,Ptr{PetscObject},Ptr{PetscInt},Ptr{PetscInt}),arg1.val,arg2,arg3,arg4,arg5)
end 
=#
function PetscOptionsHasName(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscOptionsHasName,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{PetscBool}),arg1,arg2,arg3)
    return err
end

function PetscOptionsGetInt(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscOptionsGetInt,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
    return err
end

function PetscOptionsGetBool(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscOptionsGetBool,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
    return err
end

function PetscOptionsGetReal(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscOptionsGetReal,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Cint},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
    return err
end

function PetscOptionsGetScalar(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscOptionsGetScalar,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Float64},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
    return err
end

function PetscOptionsGetIntArray(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscOptionsGetIntArray,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function PetscOptionsGetRealArray(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscOptionsGetRealArray,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Cint},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,PetscReal,arg3,arg4)
    return err
end

function PetscOptionsGetScalarArray(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscOptionsGetScalarArray,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Float64},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function PetscOptionsGetBoolArray(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscOptionsGetBoolArray,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{PetscBool},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function PetscOptionsGetString(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),size_t::Integer,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscOptionsGetString,petsc1),PetscErrorCode,(Cstring,Cstring,Cstring,Cint,Ptr{PetscBool}),arg1,arg2,arg3,size_t,arg4)
    return err
end

function PetscOptionsGetStringArray(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscOptionsGetStringArray,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Ptr{Uint8}},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function PetscOptionsGetEList(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscOptionsGetEList,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Ptr{Uint8}},PetscInt,Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function PetscOptionsGetEnum(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}),arg4::Union(Ptr{PetscEnum},StridedArray{PetscEnum},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscOptionsGetEnum,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Ptr{Uint8}},Ptr{PetscEnum},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function PetscOptionsGetEnumArray(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}),arg4::Union(Ptr{PetscEnum},StridedArray{PetscEnum},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscOptionsGetEnumArray,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Ptr{Uint8}},Ptr{PetscEnum},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function PetscOptionsValidKey(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscOptionsValidKey,petsc1),PetscErrorCode,(Cstring,Ptr{PetscBool}),arg1,arg2)
    return err
end

function PetscOptionsSetAlias(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol))
    err = ccall((:PetscOptionsSetAlias,petsc1),PetscErrorCode,(Cstring,Cstring),arg1,arg2)
    return err
end

function PetscOptionsSetValue(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol))
    err = ccall((:PetscOptionsSetValue,petsc1),PetscErrorCode,(Cstring,Cstring),arg1,arg2)
    return err
end

function PetscOptionsClearValue(arg0::Type{Float64},arg1::Union(ByteString,Symbol))
    err = ccall((:PetscOptionsClearValue,petsc1),PetscErrorCode,(Cstring,),arg1)
    return err
end

function PetscOptionsAllUsed(arg0::Type{Float64},arg1::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PetscOptionsAllUsed,petsc1),PetscErrorCode,(Ptr{PetscInt},),arg1)
    return err
end

function PetscOptionsUsed(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscOptionsUsed,petsc1),PetscErrorCode,(Cstring,Ptr{PetscBool}),arg1,arg2)
    return err
end

function PetscOptionsLeft(arg0::Type{Float64})
    err = ccall((:PetscOptionsLeft,petsc1),PetscErrorCode,())
    return err
end

function PetscOptionsView(arg1::PetscViewer{Float64})
    err = ccall((:PetscOptionsView,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
    return err
end

function PetscOptionsCreate(arg0::Type{Float64})
    err = ccall((:PetscOptionsCreate,petsc1),PetscErrorCode,())
    return err
end

function PetscOptionsInsert(arg0::Type{Float64},arg1::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg2::Union(Ptr{Ptr{Ptr{Uint8}}},StridedArray{Ptr{Ptr{Uint8}}},Ptr{Void}),arg3::Union(ByteString,Symbol))
    err = ccall((:PetscOptionsInsert,petsc1),PetscErrorCode,(Ptr{Cint},Ptr{Ptr{Ptr{Uint8}}},Cstring),arg1,arg2,arg3)
    return err
end

function PetscOptionsInsertFile(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::PetscBool)
    err = ccall((:PetscOptionsInsertFile,petsc1),PetscErrorCode,(comm_type,Cstring,PetscBool),arg1.val,arg2,arg3)
    return err
end

function PetscOptionsInsertString(arg0::Type{Float64},arg1::Union(ByteString,Symbol))
    err = ccall((:PetscOptionsInsertString,petsc1),PetscErrorCode,(Cstring,),arg1)
    return err
end

function PetscOptionsDestroy(arg0::Type{Float64})
    err = ccall((:PetscOptionsDestroy,petsc1),PetscErrorCode,())
    return err
end

function PetscOptionsClear(arg0::Type{Float64})
    err = ccall((:PetscOptionsClear,petsc1),PetscErrorCode,())
    return err
end

function PetscOptionsPrefixPush(arg0::Type{Float64},arg1::Union(ByteString,Symbol))
    err = ccall((:PetscOptionsPrefixPush,petsc1),PetscErrorCode,(Cstring,),arg1)
    return err
end

function PetscOptionsPrefixPop(arg0::Type{Float64})
    err = ccall((:PetscOptionsPrefixPop,petsc1),PetscErrorCode,())
    return err
end

function PetscOptionsReject(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol))
    err = ccall((:PetscOptionsReject,petsc1),PetscErrorCode,(Cstring,Cstring),arg1,arg2)
    return err
end

function PetscOptionsGetAll(arg0::Type{Float64},arg1::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:PetscOptionsGetAll,petsc1),PetscErrorCode,(Ptr{Ptr{Uint8}},),arg1)
    return err
end

function PetscOptionsGetenv(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),size_t::Integer,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscOptionsGetenv,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint,Ptr{PetscBool}),arg1.val,arg2,arg3,size_t,arg4)
    return err
end

function PetscOptionsStringToInt(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PetscOptionsStringToInt,petsc1),PetscErrorCode,(Cstring,Ptr{PetscInt}),arg1,arg2)
    return err
end

function PetscOptionsStringToReal(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:PetscOptionsStringToReal,petsc1),PetscErrorCode,(Cstring,Ptr{Cint}),arg1,arg2)
    return err
end

function PetscOptionsStringToBool(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscOptionsStringToBool,petsc1),PetscErrorCode,(Cstring,Ptr{PetscBool}),arg1,arg2)
    return err
end

function PetscOptionsMonitorSet(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscOptionsMonitorSet,petsc1),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
    return err
end

function PetscOptionsMonitorCancel(arg0::Type{Float64})
    err = ccall((:PetscOptionsMonitorCancel,petsc1),PetscErrorCode,())
    return err
end

function PetscOptionsMonitorDefault(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscOptionsMonitorDefault,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Void}),arg1,arg2,arg3)
    return err
end

#= skipping function with undefined symbols: 
 function PetscOptionsBegin_Private(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::MPI_Comm,arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::Union(ByteString,Symbol))
    ccall((:PetscOptionsBegin_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},comm_type,Cstring,Cstring,Cstring),arg1,arg2.val,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectOptionsBegin_Private(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::PetscObject)
    ccall((:PetscObjectOptionsBegin_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},PetscObject),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsEnd_Private(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}))
    ccall((:PetscOptionsEnd_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsHead(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ByteString,Symbol))
    ccall((:PetscOptionsHead,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsEnum_Private(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}),arg6::PetscEnum,arg7::Union(Ptr{PetscEnum},StridedArray{PetscEnum},Ptr{Void}),arg8::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsEnum_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{Ptr{Uint8}},PetscEnum,Ptr{PetscEnum},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsInt_Private(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::Integer,arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsInt_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,PetscInt,Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsReal_Private(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),PetscReal::Integer,arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg6::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsReal_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Cint,Ptr{Cint},Ptr{PetscBool}),arg1,arg2,arg3,arg4,PetscReal,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsScalar_Private(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::Float64,arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsScalar_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Float64,Ptr{Float64},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsName_Private(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsName_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsString_Private(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::Union(ByteString,Symbol),arg6::Union(ByteString,Symbol),size_t::Integer,arg7::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsString_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Cstring,Cstring,Cint,Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,size_t,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsBool_Private(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::PetscBool,arg6::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg7::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsBool_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,PetscBool,Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsBoolGroupBegin_Private(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsBoolGroupBegin_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsBoolGroup_Private(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsBoolGroup_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsBoolGroupEnd_Private(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsBoolGroupEnd_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsFList_Private(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::PetscFunctionList,arg6::Union(ByteString,Symbol),arg7::Union(ByteString,Symbol),size_t::Integer,arg8::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsFList_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,PetscFunctionList,Cstring,Cstring,Cint,Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,size_t,arg8)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsEList_Private(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}),arg6::Integer,arg7::Union(ByteString,Symbol),arg8::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg9::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsEList_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{Ptr{Uint8}},PetscInt,Cstring,Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsRealArray_Private(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsRealArray_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{Cint},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,PetscReal,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsScalarArray_Private(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsScalarArray_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{Float64},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsIntArray_Private(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsIntArray_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsStringArray_Private(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsStringArray_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{Ptr{Uint8}},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsBoolArray_Private(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsBoolArray_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{PetscBool},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsEnumArray_Private(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}),arg6::Union(Ptr{PetscEnum},StridedArray{PetscEnum},Ptr{Void}),arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg8::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsEnumArray_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{Ptr{Uint8}},Ptr{PetscEnum},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end 
=#
function PetscOptionsSetFromOptions(arg0::Type{Float64})
    err = ccall((:PetscOptionsSetFromOptions,petsc1),PetscErrorCode,())
    return err
end

function PetscOptionsSAWsDestroy(arg0::Type{Float64})
    err = ccall((:PetscOptionsSAWsDestroy,petsc1),PetscErrorCode,())
    return err
end

function PetscMemoryShowUsage(arg1::PetscViewer{Float64},arg2::Union(ByteString,Symbol))
    err = ccall((:PetscMemoryShowUsage,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function PetscObjectPrintClassNamePrefixType(arg1::PetscObject,arg2::PetscViewer{Float64})
    ccall((:PetscObjectPrintClassNamePrefixType,petsc1),PetscErrorCode,(PetscObject,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectView(arg1::PetscObject,arg2::PetscViewer{Float64})
    ccall((:PetscObjectView,petsc1),PetscErrorCode,(PetscObject,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectQueryFunction_Private(arg0::Type{Float64},arg1::PetscObject,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:PetscObjectQueryFunction_Private,petsc1),PetscErrorCode,(PetscObject,Cstring,Ptr{Ptr{Void}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectSetOptionsPrefix(arg0::Type{Float64},arg1::PetscObject,arg2::Union(ByteString,Symbol))
    ccall((:PetscObjectSetOptionsPrefix,petsc1),PetscErrorCode,(PetscObject,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectAppendOptionsPrefix(arg0::Type{Float64},arg1::PetscObject,arg2::Union(ByteString,Symbol))
    ccall((:PetscObjectAppendOptionsPrefix,petsc1),PetscErrorCode,(PetscObject,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectPrependOptionsPrefix(arg0::Type{Float64},arg1::PetscObject,arg2::Union(ByteString,Symbol))
    ccall((:PetscObjectPrependOptionsPrefix,petsc1),PetscErrorCode,(PetscObject,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectGetOptionsPrefix(arg0::Type{Float64},arg1::PetscObject,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:PetscObjectGetOptionsPrefix,petsc1),PetscErrorCode,(PetscObject,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectChangeTypeName(arg0::Type{Float64},arg1::PetscObject,arg2::Union(ByteString,Symbol))
    ccall((:PetscObjectChangeTypeName,petsc1),PetscErrorCode,(PetscObject,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectRegisterDestroy(arg0::Type{Float64},arg1::PetscObject)
    ccall((:PetscObjectRegisterDestroy,petsc1),PetscErrorCode,(PetscObject,),arg1)
end 
=#
function PetscObjectRegisterDestroyAll(arg0::Type{Float64})
    err = ccall((:PetscObjectRegisterDestroyAll,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function PetscObjectViewFromOptions(arg0::Type{Float64},arg1::PetscObject,arg2::PetscObject,arg3::Union(ByteString,Symbol))
    ccall((:PetscObjectViewFromOptions,petsc1),PetscErrorCode,(PetscObject,PetscObject,Cstring),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectName(arg0::Type{Float64},arg1::PetscObject)
    ccall((:PetscObjectName,petsc1),PetscErrorCode,(PetscObject,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectTypeCompare(arg0::Type{Float64},arg1::PetscObject,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscObjectTypeCompare,petsc1),PetscErrorCode,(PetscObject,Cstring,Ptr{PetscBool}),arg1,arg2,arg3)
end 
=#
function PetscRegisterFinalize(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscRegisterFinalize,petsc1),PetscErrorCode,(Ptr{Void},),arg1)
    return err
end

function PetscRegisterFinalizeAll(arg0::Type{Float64})
    err = ccall((:PetscRegisterFinalizeAll,petsc1),PetscErrorCode,())
    return err
end

function PetscDLOpen(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::PetscDLMode,arg3::Union(Ptr{PetscDLHandle},StridedArray{PetscDLHandle},Ptr{Void}))
    err = ccall((:PetscDLOpen,petsc1),PetscErrorCode,(Cstring,PetscDLMode,Ptr{PetscDLHandle}),arg1,arg2,arg3)
    return err
end

function PetscDLClose(arg0::Type{Float64},arg1::Union(Ptr{PetscDLHandle},StridedArray{PetscDLHandle},Ptr{Void}))
    err = ccall((:PetscDLClose,petsc1),PetscErrorCode,(Ptr{PetscDLHandle},),arg1)
    return err
end

function PetscDLSym(arg0::Type{Float64},arg1::PetscDLHandle,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    err = ccall((:PetscDLSym,petsc1),PetscErrorCode,(PetscDLHandle,Cstring,Ptr{Ptr{Void}}),arg1,arg2,arg3)
    return err
end

#= skipping function with undefined symbols: 
 function PetscObjectsDump(arg0::Type{Float64},arg1::Union(Ptr{FILE},StridedArray{FILE},Ptr{Void}),arg2::PetscBool)
    ccall((:PetscObjectsDump,petsc1),PetscErrorCode,(Ptr{FILE},PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectListDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscObjectList},StridedArray{PetscObjectList},Ptr{Void}))
    ccall((:PetscObjectListDestroy,petsc1),PetscErrorCode,(Ptr{PetscObjectList},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectListFind(arg0::Type{Float64},arg1::PetscObjectList,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{PetscObject},StridedArray{PetscObject},Ptr{Void}))
    ccall((:PetscObjectListFind,petsc1),PetscErrorCode,(PetscObjectList,Cstring,Ptr{PetscObject}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectListReverseFind(arg0::Type{Float64},arg1::PetscObjectList,arg2::PetscObject,arg3::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}),arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscObjectListReverseFind,petsc1),PetscErrorCode,(PetscObjectList,PetscObject,Ptr{Ptr{Uint8}},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectListAdd(arg0::Type{Float64},arg1::Union(Ptr{PetscObjectList},StridedArray{PetscObjectList},Ptr{Void}),arg2::Union(ByteString,Symbol),arg3::PetscObject)
    ccall((:PetscObjectListAdd,petsc1),PetscErrorCode,(Ptr{PetscObjectList},Cstring,PetscObject),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectListRemoveReference(arg0::Type{Float64},arg1::Union(Ptr{PetscObjectList},StridedArray{PetscObjectList},Ptr{Void}),arg2::Union(ByteString,Symbol))
    ccall((:PetscObjectListRemoveReference,petsc1),PetscErrorCode,(Ptr{PetscObjectList},Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectListDuplicate(arg0::Type{Float64},arg1::PetscObjectList,arg2::Union(Ptr{PetscObjectList},StridedArray{PetscObjectList},Ptr{Void}))
    ccall((:PetscObjectListDuplicate,petsc1),PetscErrorCode,(PetscObjectList,Ptr{PetscObjectList}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFunctionListAdd_Private(arg0::Type{Float64},arg1::Union(Ptr{PetscFunctionList},StridedArray{PetscFunctionList},Ptr{Void}),arg2::Union(ByteString,Symbol),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscFunctionListAdd_Private,petsc1),PetscErrorCode,(Ptr{PetscFunctionList},Cstring,Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFunctionListDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscFunctionList},StridedArray{PetscFunctionList},Ptr{Void}))
    ccall((:PetscFunctionListDestroy,petsc1),PetscErrorCode,(Ptr{PetscFunctionList},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFunctionListFind_Private(arg0::Type{Float64},arg1::PetscFunctionList,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:PetscFunctionListFind_Private,petsc1),PetscErrorCode,(PetscFunctionList,Cstring,Ptr{Ptr{Void}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFunctionListPrintTypes(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{FILE},StridedArray{FILE},Ptr{Void}),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::Union(ByteString,Symbol),arg6::Union(ByteString,Symbol),arg7::PetscFunctionList,arg8::Union(ByteString,Symbol))
    ccall((:PetscFunctionListPrintTypes,petsc1),PetscErrorCode,(comm_type,Ptr{FILE},Cstring,Cstring,Cstring,Cstring,PetscFunctionList,Cstring),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFunctionListDuplicate(arg0::Type{Float64},arg1::PetscFunctionList,arg2::Union(Ptr{PetscFunctionList},StridedArray{PetscFunctionList},Ptr{Void}))
    ccall((:PetscFunctionListDuplicate,petsc1),PetscErrorCode,(PetscFunctionList,Ptr{PetscFunctionList}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFunctionListView(arg1::PetscFunctionList,arg2::PetscViewer{Float64})
    ccall((:PetscFunctionListView,petsc1),PetscErrorCode,(PetscFunctionList,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFunctionListGet(arg0::Type{Float64},arg1::PetscFunctionList,arg2::Union(Ptr{Ptr{Ptr{Uint8}}},StridedArray{Ptr{Ptr{Uint8}}},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscFunctionListGet,petsc1),PetscErrorCode,(PetscFunctionList,Ptr{Ptr{Ptr{Uint8}}},Ptr{Cint}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDLLibraryAppend(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscDLLibrary},StridedArray{PetscDLLibrary},Ptr{Void}),arg3::Union(ByteString,Symbol))
    ccall((:PetscDLLibraryAppend,petsc1),PetscErrorCode,(comm_type,Ptr{PetscDLLibrary},Cstring),arg1.val,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDLLibraryPrepend(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscDLLibrary},StridedArray{PetscDLLibrary},Ptr{Void}),arg3::Union(ByteString,Symbol))
    ccall((:PetscDLLibraryPrepend,petsc1),PetscErrorCode,(comm_type,Ptr{PetscDLLibrary},Cstring),arg1.val,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDLLibrarySym(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscDLLibrary},StridedArray{PetscDLLibrary},Ptr{Void}),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:PetscDLLibrarySym,petsc1),PetscErrorCode,(comm_type,Ptr{PetscDLLibrary},Cstring,Cstring,Ptr{Ptr{Void}}),arg1.val,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDLLibraryPrintPath(arg0::Type{Float64},arg1::PetscDLLibrary)
    ccall((:PetscDLLibraryPrintPath,petsc1),PetscErrorCode,(PetscDLLibrary,),arg1)
end 
=#
function PetscDLLibraryRetrieve(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),size_t::Integer,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscDLLibraryRetrieve,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint,Ptr{PetscBool}),arg1.val,arg2,arg3,size_t,arg4)
    return err
end

#= skipping function with undefined symbols: 
 function PetscDLLibraryOpen(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{PetscDLLibrary},StridedArray{PetscDLLibrary},Ptr{Void}))
    ccall((:PetscDLLibraryOpen,petsc1),PetscErrorCode,(comm_type,Cstring,Ptr{PetscDLLibrary}),arg1.val,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDLLibraryClose(arg0::Type{Float64},arg1::PetscDLLibrary)
    ccall((:PetscDLLibraryClose,petsc1),PetscErrorCode,(PetscDLLibrary,),arg1)
end 
=#
function PetscSplitOwnership(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PetscSplitOwnership,petsc1),PetscErrorCode,(comm_type,Ptr{PetscInt},Ptr{PetscInt}),arg1.val,arg2,arg3)
    return err
end

function PetscSplitOwnershipBlock(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PetscSplitOwnershipBlock,petsc1),PetscErrorCode,(comm_type,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1.val,arg2,arg3,arg4)
    return err
end

function PetscSequentialPhaseBegin(arg0::Type{Float64},arg1::MPI_Comm,arg2::PetscMPIInt)
    err = ccall((:PetscSequentialPhaseBegin,petsc1),PetscErrorCode,(comm_type,PetscMPIInt),arg1.val,arg2)
    return err
end

function PetscSequentialPhaseEnd(arg0::Type{Float64},arg1::MPI_Comm,arg2::PetscMPIInt)
    err = ccall((:PetscSequentialPhaseEnd,petsc1),PetscErrorCode,(comm_type,PetscMPIInt),arg1.val,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function PetscBarrier(arg0::Type{Float64},arg1::PetscObject)
    ccall((:PetscBarrier,petsc1),PetscErrorCode,(PetscObject,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscMPIDump(arg0::Type{Float64},arg1::Union(Ptr{FILE},StridedArray{FILE},Ptr{Void}))
    ccall((:PetscMPIDump,petsc1),PetscErrorCode,(Ptr{FILE},),arg1)
end 
=#
function PetscInfoDeactivateClass(arg0::Type{Float64},arg1::PetscClassId)
    err = ccall((:PetscInfoDeactivateClass,petsc1),PetscErrorCode,(PetscClassId,),arg1)
    return err
end

function PetscInfoActivateClass(arg0::Type{Float64},arg1::PetscClassId)
    err = ccall((:PetscInfoActivateClass,petsc1),PetscErrorCode,(PetscClassId,),arg1)
    return err
end

#= skipping function with undefined symbols: 
 function PetscLogGetStageLog(arg0::Type{Float64},arg1::Union(Ptr{PetscStageLog},StridedArray{PetscStageLog},Ptr{Void}))
    ccall((:PetscLogGetStageLog,petsc1),PetscErrorCode,(Ptr{PetscStageLog},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscStageLogGetCurrent(arg0::Type{Float64},arg1::PetscStageLog,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscStageLogGetCurrent,petsc1),PetscErrorCode,(PetscStageLog,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscStageLogGetEventPerfLog(arg0::Type{Float64},arg1::PetscStageLog,arg2::Integer,arg3::Union(Ptr{PetscEventPerfLog},StridedArray{PetscEventPerfLog},Ptr{Void}))
    ccall((:PetscStageLogGetEventPerfLog,petsc1),PetscErrorCode,(PetscStageLog,Cint,Ptr{PetscEventPerfLog}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscLogObjectParent(arg0::Type{Float64},arg1::PetscObject,arg2::PetscObject)
    ccall((:PetscLogObjectParent,petsc1),PetscErrorCode,(PetscObject,PetscObject),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscLogObjectMemory(arg0::Type{Float64},arg1::PetscObject,arg2::PetscLogDouble)
    ccall((:PetscLogObjectMemory,petsc1),PetscErrorCode,(PetscObject,PetscLogDouble),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscIntStackCreate(arg0::Type{Float64},arg1::Union(Ptr{PetscIntStack},StridedArray{PetscIntStack},Ptr{Void}))
    ccall((:PetscIntStackCreate,petsc1),PetscErrorCode,(Ptr{PetscIntStack},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscIntStackDestroy(arg0::Type{Float64},arg1::PetscIntStack)
    ccall((:PetscIntStackDestroy,petsc1),PetscErrorCode,(PetscIntStack,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscIntStackPush(arg0::Type{Float64},arg1::PetscIntStack,arg2::Integer)
    ccall((:PetscIntStackPush,petsc1),PetscErrorCode,(PetscIntStack,Cint),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscIntStackPop(arg0::Type{Float64},arg1::PetscIntStack,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscIntStackPop,petsc1),PetscErrorCode,(PetscIntStack,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscIntStackTop(arg0::Type{Float64},arg1::PetscIntStack,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscIntStackTop,petsc1),PetscErrorCode,(PetscIntStack,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscIntStackEmpty(arg0::Type{Float64},arg1::PetscIntStack,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscIntStackEmpty,petsc1),PetscErrorCode,(PetscIntStack,Ptr{PetscBool}),arg1,arg2)
end 
=#
function PetscFixFilename(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol))
    err = ccall((:PetscFixFilename,petsc1),PetscErrorCode,(Cstring,Cstring),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function PetscFOpen(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(Ptr{Ptr{FILE}},StridedArray{Ptr{FILE}},Ptr{Void}))
    ccall((:PetscFOpen,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Ptr{Ptr{FILE}}),arg1.val,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFClose(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{FILE},StridedArray{FILE},Ptr{Void}))
    ccall((:PetscFClose,petsc1),PetscErrorCode,(comm_type,Ptr{FILE}),arg1.val,arg2)
end 
=#
function PetscVSNPrintf(arg0::Type{Float64},arg1::Union(ByteString,Symbol),size_t::Integer,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),va_list::Integer)
    err = ccall((:PetscVSNPrintf,petsc1),PetscErrorCode,(Cstring,Cint,Cstring,Ptr{Cint},Cint),arg1,size_t,arg2,arg3,va_list)
    return err
end

#= skipping function with undefined symbols: 
 function PetscVFPrintfDefault(arg0::Type{Float64},arg1::Union(Ptr{FILE},StridedArray{FILE},Ptr{Void}),arg2::Union(ByteString,Symbol),va_list::Integer)
    ccall((:PetscVFPrintfDefault,petsc1),PetscErrorCode,(Ptr{FILE},Cstring,Cint),arg1,arg2,va_list)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSynchronizedFlush(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{FILE},StridedArray{FILE},Ptr{Void}))
    ccall((:PetscSynchronizedFlush,petsc1),PetscErrorCode,(comm_type,Ptr{FILE}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSynchronizedFGets(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{FILE},StridedArray{FILE},Ptr{Void}),size_t::Integer,arg3::Union(ByteString,Symbol))
    ccall((:PetscSynchronizedFGets,petsc1),PetscErrorCode,(comm_type,Ptr{FILE},Cint,Cstring),arg1.val,arg2,size_t,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscStartMatlab(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(Ptr{Ptr{FILE}},StridedArray{Ptr{FILE}},Ptr{Void}))
    ccall((:PetscStartMatlab,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Ptr{Ptr{FILE}}),arg1.val,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscStartJava(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(Ptr{Ptr{FILE}},StridedArray{Ptr{FILE}},Ptr{Void}))
    ccall((:PetscStartJava,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Ptr{Ptr{FILE}}),arg1.val,arg2,arg3,arg4)
end 
=#
function PetscGetPetscDir(arg0::Type{Float64},arg1::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:PetscGetPetscDir,petsc1),PetscErrorCode,(Ptr{Ptr{Uint8}},),arg1)
    return err
end

function PetscPopUpSelect(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Integer,arg5::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}),arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:PetscPopUpSelect,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint,Ptr{Ptr{Uint8}},Ptr{Cint}),arg1.val,arg2,arg3,arg4,arg5,arg6)
    return err
end

#= skipping function with undefined symbols: 
 function PetscContainerGetPointer(arg0::Type{Float64},arg1::PetscContainer,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:PetscContainerGetPointer,petsc1),PetscErrorCode,(PetscContainer,Ptr{Ptr{Void}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscContainerSetPointer(arg0::Type{Float64},arg1::PetscContainer,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscContainerSetPointer,petsc1),PetscErrorCode,(PetscContainer,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscContainerDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscContainer},StridedArray{PetscContainer},Ptr{Void}))
    ccall((:PetscContainerDestroy,petsc1),PetscErrorCode,(Ptr{PetscContainer},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscContainerCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscContainer},StridedArray{PetscContainer},Ptr{Void}))
    ccall((:PetscContainerCreate,petsc1),PetscErrorCode,(comm_type,Ptr{PetscContainer}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscContainerSetUserDestroy(arg0::Type{Float64},arg1::PetscContainer,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscContainerSetUserDestroy,petsc1),PetscErrorCode,(PetscContainer,Ptr{Void}),arg1,arg2)
end 
=#
function PetscIntView(arg1::Integer,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::PetscViewer{Float64})
    err = ccall((:PetscIntView,petsc1),PetscErrorCode,(PetscInt,Ptr{PetscInt},PetscViewer{Float64}),arg1,arg2,arg3)
    return err
end

function PetscRealView(arg1::Integer,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg2::PetscViewer{Float64})
    err = ccall((:PetscRealView,petsc1),PetscErrorCode,(PetscInt,Ptr{Cint},PetscViewer{Float64}),arg1,PetscReal,arg2)
    return err
end

function PetscScalarView(arg1::Integer,arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg3::PetscViewer{Float64})
    err = ccall((:PetscScalarView,petsc1),PetscErrorCode,(PetscInt,Ptr{Float64},PetscViewer{Float64}),arg1,arg2,arg3)
    return err
end

function PetscGetHostName(arg0::Type{Float64},arg1::Union(ByteString,Symbol),size_t::Integer)
    err = ccall((:PetscGetHostName,petsc1),PetscErrorCode,(Cstring,Cint),arg1,size_t)
    return err
end

function PetscGetUserName(arg0::Type{Float64},arg1::Union(ByteString,Symbol),size_t::Integer)
    err = ccall((:PetscGetUserName,petsc1),PetscErrorCode,(Cstring,Cint),arg1,size_t)
    return err
end

function PetscGetProgramName(arg0::Type{Float64},arg1::Union(ByteString,Symbol),size_t::Integer)
    err = ccall((:PetscGetProgramName,petsc1),PetscErrorCode,(Cstring,Cint),arg1,size_t)
    return err
end

function PetscSetProgramName(arg0::Type{Float64},arg1::Union(ByteString,Symbol))
    err = ccall((:PetscSetProgramName,petsc1),PetscErrorCode,(Cstring,),arg1)
    return err
end

function PetscGetDate(arg0::Type{Float64},arg1::Union(ByteString,Symbol),size_t::Integer)
    err = ccall((:PetscGetDate,petsc1),PetscErrorCode,(Cstring,Cint),arg1,size_t)
    return err
end

function PetscGetVersion(arg0::Type{Float64},arg1::Union(ByteString,Symbol),size_t::Integer)
    err = ccall((:PetscGetVersion,petsc1),PetscErrorCode,(Cstring,Cint),arg1,size_t)
    return err
end

function PetscSortInt(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PetscSortInt,petsc1),PetscErrorCode,(PetscInt,Ptr{PetscInt}),arg1,arg2)
    return err
end

function PetscSortRemoveDupsInt(arg0::Type{Float64},arg1::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PetscSortRemoveDupsInt,petsc1),PetscErrorCode,(Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2)
    return err
end

function PetscFindInt(arg0::Type{Float64},arg1::Integer,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PetscFindInt,petsc1),PetscErrorCode,(PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
    return err
end

function PetscSortIntWithPermutation(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PetscSortIntWithPermutation,petsc1),PetscErrorCode,(PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function PetscSortStrWithPermutation(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PetscSortStrWithPermutation,petsc1),PetscErrorCode,(PetscInt,Ptr{Ptr{Uint8}},Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function PetscSortIntWithArray(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PetscSortIntWithArray,petsc1),PetscErrorCode,(PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function PetscSortIntWithArrayPair(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PetscSortIntWithArrayPair,petsc1),PetscErrorCode,(PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
    return err
end

function PetscSortMPIInt(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}))
    err = ccall((:PetscSortMPIInt,petsc1),PetscErrorCode,(PetscInt,Ptr{PetscMPIInt}),arg1,arg2)
    return err
end

function PetscSortRemoveDupsMPIInt(arg0::Type{Float64},arg1::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg2::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}))
    err = ccall((:PetscSortRemoveDupsMPIInt,petsc1),PetscErrorCode,(Ptr{PetscInt},Ptr{PetscMPIInt}),arg1,arg2)
    return err
end

function PetscSortMPIIntWithArray(arg0::Type{Float64},arg1::PetscMPIInt,arg2::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}),arg3::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}))
    err = ccall((:PetscSortMPIIntWithArray,petsc1),PetscErrorCode,(PetscMPIInt,Ptr{PetscMPIInt},Ptr{PetscMPIInt}),arg1,arg2,arg3)
    return err
end

function PetscSortIntWithScalarArray(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:PetscSortIntWithScalarArray,petsc1),PetscErrorCode,(PetscInt,Ptr{PetscInt},Ptr{Float64}),arg1,arg2,arg3)
    return err
end

function PetscSortIntWithDataArray(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),size_t::Integer,arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscSortIntWithDataArray,petsc1),PetscErrorCode,(PetscInt,Ptr{PetscInt},Ptr{Void},Cint,Ptr{Void}),arg1,arg2,arg3,size_t,arg4)
    return err
end

function PetscSortReal(arg0::Type{Float64},arg1::Integer,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:PetscSortReal,petsc1),PetscErrorCode,(PetscInt,Ptr{Cint}),arg1,PetscReal)
    return err
end

function PetscSortRealWithPermutation(arg0::Type{Float64},arg1::Integer,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PetscSortRealWithPermutation,petsc1),PetscErrorCode,(PetscInt,Ptr{Cint},Ptr{PetscInt}),arg1,PetscReal,arg2)
    return err
end

function PetscSortRemoveDupsReal(arg0::Type{Float64},arg1::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:PetscSortRemoveDupsReal,petsc1),PetscErrorCode,(Ptr{PetscInt},Ptr{Cint}),arg1,PetscReal)
    return err
end

function PetscSortSplit(arg0::Type{Float64},arg1::Integer,arg2::Integer,arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PetscSortSplit,petsc1),PetscErrorCode,(PetscInt,PetscInt,Ptr{Float64},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
    return err
end

function PetscSortSplitReal(arg0::Type{Float64},arg1::Integer,arg2::Integer,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PetscSortSplitReal,petsc1),PetscErrorCode,(PetscInt,PetscInt,Ptr{Cint},Ptr{PetscInt}),arg1,arg2,PetscReal,arg3)
    return err
end

function PetscProcessTree(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg6::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg7::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg8::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:PetscProcessTree,petsc1),PetscErrorCode,(PetscInt,Ptr{PetscBool},Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
    return err
end

function PetscMergeIntArrayPair(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg8::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg9::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:PetscMergeIntArrayPair,petsc1),PetscErrorCode,(PetscInt,Ptr{PetscInt},Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
    return err
end

function PetscMergeIntArray(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:PetscMergeIntArray,petsc1),PetscErrorCode,(PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function PetscSetDisplay(arg0::Type{Float64})
    err = ccall((:PetscSetDisplay,petsc1),PetscErrorCode,())
    return err
end

function PetscGetDisplay(arg0::Type{Float64},arg1::Union(ByteString,Symbol),size_t::Integer)
    err = ccall((:PetscGetDisplay,petsc1),PetscErrorCode,(Cstring,Cint),arg1,size_t)
    return err
end

function PetscRandomInitializePackage(arg0::Type{Float64})
    err = ccall((:PetscRandomInitializePackage,petsc1),PetscErrorCode,())
    return err
end

function PetscRandomRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscRandomRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function PetscRandomSetType(arg0::Type{Float64},arg1::PetscRandom,arg2::PetscRandomType)
    ccall((:PetscRandomSetType,petsc1),PetscErrorCode,(PetscRandom,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscRandomSetFromOptions(arg0::Type{Float64},arg1::PetscRandom)
    ccall((:PetscRandomSetFromOptions,petsc1),PetscErrorCode,(PetscRandom,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscRandomGetType(arg0::Type{Float64},arg1::PetscRandom,arg2::Union(Ptr{PetscRandomType},StridedArray{PetscRandomType},Ptr{Void}))
    ccall((:PetscRandomGetType,petsc1),PetscErrorCode,(PetscRandom,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscRandomCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscRandom},StridedArray{PetscRandom},Ptr{Void}))
    ccall((:PetscRandomCreate,petsc1),PetscErrorCode,(comm_type,Ptr{PetscRandom}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscRandomGetValue(arg0::Type{Float64},arg1::PetscRandom,arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:PetscRandomGetValue,petsc1),PetscErrorCode,(PetscRandom,Ptr{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscRandomGetValueReal(arg0::Type{Float64},arg1::PetscRandom,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscRandomGetValueReal,petsc1),PetscErrorCode,(PetscRandom,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscRandomGetInterval(arg0::Type{Float64},arg1::PetscRandom,arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:PetscRandomGetInterval,petsc1),PetscErrorCode,(PetscRandom,Ptr{Float64},Ptr{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscRandomSetInterval(arg0::Type{Float64},arg1::PetscRandom,arg2::Float64,arg3::Float64)
    ccall((:PetscRandomSetInterval,petsc1),PetscErrorCode,(PetscRandom,Float64,Float64),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscRandomSetSeed(arg0::Type{Float64},arg1::PetscRandom,arg2::Culong)
    ccall((:PetscRandomSetSeed,petsc1),PetscErrorCode,(PetscRandom,Culong),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscRandomGetSeed(arg0::Type{Float64},arg1::PetscRandom,arg2::Union(Ptr{Culong},StridedArray{Culong},Ptr{Void}))
    ccall((:PetscRandomGetSeed,petsc1),PetscErrorCode,(PetscRandom,Ptr{Culong}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscRandomSeed(arg0::Type{Float64},arg1::PetscRandom)
    ccall((:PetscRandomSeed,petsc1),PetscErrorCode,(PetscRandom,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscRandomDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscRandom},StridedArray{PetscRandom},Ptr{Void}))
    ccall((:PetscRandomDestroy,petsc1),PetscErrorCode,(Ptr{PetscRandom},),arg1)
end 
=#
function PetscGetFullPath(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),size_t::Integer)
    err = ccall((:PetscGetFullPath,petsc1),PetscErrorCode,(Cstring,Cstring,Cint),arg1,arg2,size_t)
    return err
end

function PetscGetRelativePath(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),size_t::Integer)
    err = ccall((:PetscGetRelativePath,petsc1),PetscErrorCode,(Cstring,Cstring,Cint),arg1,arg2,size_t)
    return err
end

function PetscGetWorkingDirectory(arg0::Type{Float64},arg1::Union(ByteString,Symbol),size_t::Integer)
    err = ccall((:PetscGetWorkingDirectory,petsc1),PetscErrorCode,(Cstring,Cint),arg1,size_t)
    return err
end

function PetscGetRealPath(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol))
    err = ccall((:PetscGetRealPath,petsc1),PetscErrorCode,(Cstring,Cstring),arg1,arg2)
    return err
end

function PetscGetHomeDirectory(arg0::Type{Float64},arg1::Union(ByteString,Symbol),size_t::Integer)
    err = ccall((:PetscGetHomeDirectory,petsc1),PetscErrorCode,(Cstring,Cint),arg1,size_t)
    return err
end

function PetscTestFile(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Uint8,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscTestFile,petsc1),PetscErrorCode,(Cstring,Uint8,Ptr{PetscBool}),arg1,arg2,arg3)
    return err
end

function PetscTestDirectory(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Uint8,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscTestDirectory,petsc1),PetscErrorCode,(Cstring,Uint8,Ptr{PetscBool}),arg1,arg2,arg3)
    return err
end

function PetscBinaryRead(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Integer,arg4::PetscDataType)
    err = ccall((:PetscBinaryRead,petsc1),PetscErrorCode,(Cint,Ptr{Void},PetscInt,PetscDataType),arg1,arg2,arg3,arg4)
    return err
end

function PetscBinarySynchronizedRead(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Integer,arg5::PetscDataType)
    err = ccall((:PetscBinarySynchronizedRead,petsc1),PetscErrorCode,(comm_type,Cint,Ptr{Void},PetscInt,PetscDataType),arg1.val,arg2,arg3,arg4,arg5)
    return err
end

function PetscBinarySynchronizedWrite(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Integer,arg5::PetscDataType,arg6::PetscBool)
    err = ccall((:PetscBinarySynchronizedWrite,petsc1),PetscErrorCode,(comm_type,Cint,Ptr{Void},PetscInt,PetscDataType,PetscBool),arg1.val,arg2,arg3,arg4,arg5,arg6)
    return err
end

function PetscBinaryWrite(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Integer,arg4::PetscDataType,arg5::PetscBool)
    err = ccall((:PetscBinaryWrite,petsc1),PetscErrorCode,(Cint,Ptr{Void},PetscInt,PetscDataType,PetscBool),arg1,arg2,arg3,arg4,arg5)
    return err
end

function PetscBinaryOpen(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::PetscFileMode,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:PetscBinaryOpen,petsc1),PetscErrorCode,(Cstring,PetscFileMode,Ptr{Cint}),arg1,arg2,arg3)
    return err
end

function PetscBinaryClose(arg0::Type{Float64},arg1::Integer)
    err = ccall((:PetscBinaryClose,petsc1),PetscErrorCode,(Cint,),arg1)
    return err
end

function PetscSharedTmp(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscSharedTmp,petsc1),PetscErrorCode,(comm_type,Ptr{PetscBool}),arg1.val,arg2)
    return err
end

function PetscSharedWorkingDirectory(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscSharedWorkingDirectory,petsc1),PetscErrorCode,(comm_type,Ptr{PetscBool}),arg1.val,arg2)
    return err
end

function PetscGetTmp(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),size_t::Integer)
    err = ccall((:PetscGetTmp,petsc1),PetscErrorCode,(comm_type,Cstring,Cint),arg1.val,arg2,size_t)
    return err
end

function PetscFileRetrieve(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),size_t::Integer,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscFileRetrieve,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint,Ptr{PetscBool}),arg1.val,arg2,arg3,size_t,arg4)
    return err
end

function PetscLs(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),size_t::Integer,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscLs,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint,Ptr{PetscBool}),arg1.val,arg2,arg3,size_t,arg4)
    return err
end

function PetscOpenSocket(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Integer,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:PetscOpenSocket,petsc1),PetscErrorCode,(Cstring,Cint,Ptr{Cint}),arg1,arg2,arg3)
    return err
end

function PetscBinarySeek(arg0::Type{Float64},arg1::Integer,off_t::Integer,arg2::PetscBinarySeekType,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:PetscBinarySeek,petsc1),PetscErrorCode,(Cint,Cint,PetscBinarySeekType,Ptr{Cint}),arg1,off_t,arg2,arg3)
    return err
end

function PetscBinarySynchronizedSeek(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,off_t::Integer,arg3::PetscBinarySeekType,arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:PetscBinarySynchronizedSeek,petsc1),PetscErrorCode,(comm_type,Cint,Cint,PetscBinarySeekType,Ptr{Cint}),arg1.val,arg2,off_t,arg3,arg4)
    return err
end

function PetscByteSwap(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::PetscDataType,arg3::Integer)
    err = ccall((:PetscByteSwap,petsc1),PetscErrorCode,(Ptr{Void},PetscDataType,PetscInt),arg1,arg2,arg3)
    return err
end

function PetscSetDebugTerminal(arg0::Type{Float64},arg1::Union(ByteString,Symbol))
    err = ccall((:PetscSetDebugTerminal,petsc1),PetscErrorCode,(Cstring,),arg1)
    return err
end

function PetscSetDebugger(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::PetscBool)
    err = ccall((:PetscSetDebugger,petsc1),PetscErrorCode,(Cstring,PetscBool),arg1,arg2)
    return err
end

function PetscSetDefaultDebugger(arg0::Type{Float64})
    err = ccall((:PetscSetDefaultDebugger,petsc1),PetscErrorCode,())
    return err
end

function PetscSetDebuggerFromString(arg0::Type{Float64},arg1::Union(ByteString,Symbol))
    err = ccall((:PetscSetDebuggerFromString,petsc1),PetscErrorCode,(Cstring,),arg1)
    return err
end

function PetscAttachDebugger(arg0::Type{Float64})
    err = ccall((:PetscAttachDebugger,petsc1),PetscErrorCode,())
    return err
end

function PetscStopForDebugger(arg0::Type{Float64})
    err = ccall((:PetscStopForDebugger,petsc1),PetscErrorCode,())
    return err
end

function PetscGatherNumberOfMessages(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}),arg3::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}),arg4::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}))
    err = ccall((:PetscGatherNumberOfMessages,petsc1),PetscErrorCode,(comm_type,Ptr{PetscMPIInt},Ptr{PetscMPIInt},Ptr{PetscMPIInt}),arg1.val,arg2,arg3,arg4)
    return err
end

function PetscGatherMessageLengths(arg0::Type{Float64},arg1::MPI_Comm,arg2::PetscMPIInt,arg3::PetscMPIInt,arg4::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}),arg5::Union(Ptr{Ptr{PetscMPIInt}},StridedArray{Ptr{PetscMPIInt}},Ptr{Void}),arg6::Union(Ptr{Ptr{PetscMPIInt}},StridedArray{Ptr{PetscMPIInt}},Ptr{Void}))
    err = ccall((:PetscGatherMessageLengths,petsc1),PetscErrorCode,(comm_type,PetscMPIInt,PetscMPIInt,Ptr{PetscMPIInt},Ptr{Ptr{PetscMPIInt}},Ptr{Ptr{PetscMPIInt}}),arg1.val,arg2,arg3,arg4,arg5,arg6)
    return err
end

function PetscGatherMessageLengths2(arg0::Type{Float64},arg1::MPI_Comm,arg2::PetscMPIInt,arg3::PetscMPIInt,arg4::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}),arg5::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}),arg6::Union(Ptr{Ptr{PetscMPIInt}},StridedArray{Ptr{PetscMPIInt}},Ptr{Void}),arg7::Union(Ptr{Ptr{PetscMPIInt}},StridedArray{Ptr{PetscMPIInt}},Ptr{Void}),arg8::Union(Ptr{Ptr{PetscMPIInt}},StridedArray{Ptr{PetscMPIInt}},Ptr{Void}))
    err = ccall((:PetscGatherMessageLengths2,petsc1),PetscErrorCode,(comm_type,PetscMPIInt,PetscMPIInt,Ptr{PetscMPIInt},Ptr{PetscMPIInt},Ptr{Ptr{PetscMPIInt}},Ptr{Ptr{PetscMPIInt}},Ptr{Ptr{PetscMPIInt}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
    return err
end

#= skipping function with undefined symbols: 
 function PetscPostIrecvInt(arg0::Type{Float64},arg1::MPI_Comm,arg2::PetscMPIInt,arg3::PetscMPIInt,arg4::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}),arg5::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}),arg6::Union(Ptr{Ptr{Ptr{PetscInt}}},StridedArray{Ptr{Ptr{PetscInt}}},Ptr{Void}),arg7::Union(Ptr{Ptr{MPI_Request}},StridedArray{Ptr{MPI_Request}},Ptr{Void}))
    ccall((:PetscPostIrecvInt,petsc1),PetscErrorCode,(comm_type,PetscMPIInt,PetscMPIInt,Ptr{PetscMPIInt},Ptr{PetscMPIInt},Ptr{Ptr{Ptr{PetscInt}}},Ptr{Ptr{MPI_Request}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscPostIrecvScalar(arg0::Type{Float64},arg1::MPI_Comm,arg2::PetscMPIInt,arg3::PetscMPIInt,arg4::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}),arg5::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}),arg6::Union(Ptr{Ptr{Ptr{Float64}}},StridedArray{Ptr{Ptr{Float64}}},Ptr{Void}),arg7::Union(Ptr{Ptr{MPI_Request}},StridedArray{Ptr{MPI_Request}},Ptr{Void}))
    ccall((:PetscPostIrecvScalar,petsc1),PetscErrorCode,(comm_type,PetscMPIInt,PetscMPIInt,Ptr{PetscMPIInt},Ptr{PetscMPIInt},Ptr{Ptr{Ptr{Float64}}},Ptr{Ptr{MPI_Request}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscCommBuildTwoSided(arg0::Type{Float64},arg1::MPI_Comm,arg2::PetscMPIInt,arg3::MPI_Datatype,arg4::Integer,arg5::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}),arg6::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg8::Union(Ptr{Ptr{PetscMPIInt}},StridedArray{Ptr{PetscMPIInt}},Ptr{Void}),arg9::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscCommBuildTwoSided,petsc1),PetscErrorCode,(comm_type,PetscMPIInt,MPI_Datatype,PetscInt,Ptr{PetscMPIInt},Ptr{Void},Ptr{PetscInt},Ptr{Ptr{PetscMPIInt}},Ptr{Void}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end 
=#
function PetscCommBuildTwoSidedSetType(arg0::Type{Float64},arg1::MPI_Comm,arg2::PetscBuildTwoSidedType)
    err = ccall((:PetscCommBuildTwoSidedSetType,petsc1),PetscErrorCode,(comm_type,PetscBuildTwoSidedType),arg1.val,arg2)
    return err
end

function PetscCommBuildTwoSidedGetType(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscBuildTwoSidedType},StridedArray{PetscBuildTwoSidedType},Ptr{Void}))
    err = ccall((:PetscCommBuildTwoSidedGetType,petsc1),PetscErrorCode,(comm_type,Ptr{PetscBuildTwoSidedType}),arg1.val,arg2)
    return err
end

function PetscSSEIsEnabled(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscSSEIsEnabled,petsc1),PetscErrorCode,(comm_type,Ptr{PetscBool},Ptr{PetscBool}),arg1.val,arg2,arg3)
    return err
end

#= skipping function with undefined symbols: 
 function PetscObjectComm(arg0::Type{Float64},arg1::PetscObject)
    ccall((:PetscObjectComm,petsc1),MPI_Comm,(PetscObject,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSubcommDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscSubcomm},StridedArray{PetscSubcomm},Ptr{Void}))
    ccall((:PetscSubcommDestroy,petsc1),PetscErrorCode,(Ptr{PetscSubcomm},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSubcommSetNumber(arg0::Type{Float64},arg1::PetscSubcomm,arg2::Integer)
    ccall((:PetscSubcommSetNumber,petsc1),PetscErrorCode,(PetscSubcomm,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSubcommSetType(arg0::Type{Float64},arg1::PetscSubcomm,arg2::PetscSubcommType)
    ccall((:PetscSubcommSetType,petsc1),PetscErrorCode,(PetscSubcomm,PetscSubcommType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSubcommSetTypeGeneral(arg0::Type{Float64},arg1::PetscSubcomm,arg2::PetscMPIInt,arg3::PetscMPIInt)
    ccall((:PetscSubcommSetTypeGeneral,petsc1),PetscErrorCode,(PetscSubcomm,PetscMPIInt,PetscMPIInt),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSubcommView(arg1::PetscSubcomm,arg2::PetscViewer{Float64})
    ccall((:PetscSubcommView,petsc1),PetscErrorCode,(PetscSubcomm,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSubcommSetFromOptions(arg0::Type{Float64},arg1::PetscSubcomm)
    ccall((:PetscSubcommSetFromOptions,petsc1),PetscErrorCode,(PetscSubcomm,),arg1)
end 
=#
function PetscSegBufferCreate(arg0::Type{Float64})
    err = ccall((:PetscSegBufferCreate,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function PetscSegBufferDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscSegBuffer},StridedArray{PetscSegBuffer},Ptr{Void}))
    ccall((:PetscSegBufferDestroy,petsc1),PetscErrorCode,(Ptr{PetscSegBuffer},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSegBufferGet(arg0::Type{Float64},arg1::PetscSegBuffer,size_t::Integer,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscSegBufferGet,petsc1),PetscErrorCode,(PetscSegBuffer,Cint,Ptr{Void}),arg1,size_t,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSegBufferExtractAlloc(arg0::Type{Float64},arg1::PetscSegBuffer,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscSegBufferExtractAlloc,petsc1),PetscErrorCode,(PetscSegBuffer,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSegBufferExtractTo(arg0::Type{Float64},arg1::PetscSegBuffer,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscSegBufferExtractTo,petsc1),PetscErrorCode,(PetscSegBuffer,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSegBufferExtractInPlace(arg0::Type{Float64},arg1::PetscSegBuffer,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscSegBufferExtractInPlace,petsc1),PetscErrorCode,(PetscSegBuffer,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSegBufferGetSize(arg0::Type{Float64},arg1::PetscSegBuffer,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscSegBufferGetSize,petsc1),PetscErrorCode,(PetscSegBuffer,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSegBufferUnuse(arg0::Type{Float64},arg1::PetscSegBuffer,size_t::Integer)
    ccall((:PetscSegBufferUnuse,petsc1),PetscErrorCode,(PetscSegBuffer,Cint),arg1,size_t)
end 
=#
function PetscGoogleDriveAuthorize(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),size_t::Integer)
    err = ccall((:PetscGoogleDriveAuthorize,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint),arg1.val,arg2,arg3,size_t)
    return err
end

function PetscGoogleDriveRefresh(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),size_t::Integer)
    err = ccall((:PetscGoogleDriveRefresh,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint),arg1.val,arg2,arg3,size_t)
    return err
end

function PetscGoogleDriveUpload(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol))
    err = ccall((:PetscGoogleDriveUpload,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring),arg1.val,arg2,arg3)
    return err
end

function PetscBoxAuthorize(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),size_t::Integer)
    err = ccall((:PetscBoxAuthorize,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint),arg1.val,arg2,arg3,size_t)
    return err
end

function PetscBoxRefresh(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),size_t::Integer)
    err = ccall((:PetscBoxRefresh,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cstring,Cint),arg1.val,arg2,arg3,arg4,size_t)
    return err
end

function PetscTextBelt(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscTextBelt,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Ptr{PetscBool}),arg1.val,arg2,arg3,arg4)
    return err
end

function PetscPullJSONValue(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),size_t::Integer,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscPullJSONValue,petsc1),PetscErrorCode,(Cstring,Cstring,Cstring,Cint,Ptr{PetscBool}),arg1,arg2,arg3,size_t,arg4)
    return err
end

function PetscPushJSONValue(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),size_t::Integer)
    err = ccall((:PetscPushJSONValue,petsc1),PetscErrorCode,(Cstring,Cstring,Cstring,Cint),arg1,arg2,arg3,size_t)
    return err
end

function ISInitializePackage(arg0::Type{Float64})
    err = ccall((:ISInitializePackage,petsc1),PetscErrorCode,())
    return err
end

function ISSetType(arg1::IS{Float64},arg2::ISType)
    err = ccall((:ISSetType,petsc1),PetscErrorCode,(IS{Float64},Cstring),arg1,arg2)
    return err
end

function ISGetType(arg1::IS{Float64},arg2::Union(Ptr{ISType},StridedArray{ISType},Ptr{Void}))
    (arg2_,tmp) = symbol_get_before(arg2)
    err = ccall((:ISGetType,petsc1),PetscErrorCode,(IS{Float64},Ptr{Ptr{Uint8}}),arg1,arg2_)
    symbol_get_after(arg2_,arg2)
    return err
end

function ISRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:ISRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

function ISCreate(arg1::MPI_Comm,arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISCreate,petsc1),PetscErrorCode,(comm_type,Ptr{IS{Float64}}),arg1.val,arg2)
    return err
end

function ISCreateGeneral(arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),PetscCopyMode::Integer,arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISCreateGeneral,petsc1),PetscErrorCode,(comm_type,PetscInt,Ptr{PetscInt},Cint,Ptr{IS{Float64}}),arg1.val,arg2,arg3,PetscCopyMode,arg4)
    return err
end

function ISGeneralSetIndices(arg1::IS{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),PetscCopyMode::Integer)
    err = ccall((:ISGeneralSetIndices,petsc1),PetscErrorCode,(IS{Float64},PetscInt,Ptr{PetscInt},Cint),arg1,arg2,arg3,PetscCopyMode)
    return err
end

function ISCreateBlock(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),PetscCopyMode::Integer,arg5::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISCreateBlock,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,Ptr{PetscInt},Cint,Ptr{IS{Float64}}),arg1.val,arg2,arg3,arg4,PetscCopyMode,arg5)
    return err
end

function ISBlockSetIndices(arg1::IS{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),PetscCopyMode::Integer)
    err = ccall((:ISBlockSetIndices,petsc1),PetscErrorCode,(IS{Float64},PetscInt,PetscInt,Ptr{PetscInt},Cint),arg1,arg2,arg3,arg4,PetscCopyMode)
    return err
end

function ISCreateStride(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISCreateStride,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,Ptr{IS{Float64}}),arg1.val,arg2,arg3,arg4,arg5)
    return err
end

function ISStrideSetStride(arg1::IS{Float64},arg2::Integer,arg3::Integer,arg4::Integer)
    err = ccall((:ISStrideSetStride,petsc1),PetscErrorCode,(IS{Float64},PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
    return err
end

function ISDestroy(arg1::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISDestroy,petsc1),PetscErrorCode,(Ptr{IS{Float64}},),arg1)
    return err
end

function ISSetPermutation(arg1::IS{Float64})
    err = ccall((:ISSetPermutation,petsc1),PetscErrorCode,(IS{Float64},),arg1)
    return err
end

function ISPermutation(arg1::IS{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:ISPermutation,petsc1),PetscErrorCode,(IS{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function ISSetIdentity(arg1::IS{Float64})
    err = ccall((:ISSetIdentity,petsc1),PetscErrorCode,(IS{Float64},),arg1)
    return err
end

function ISIdentity(arg1::IS{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:ISIdentity,petsc1),PetscErrorCode,(IS{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function ISContiguousLocal(arg1::IS{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:ISContiguousLocal,petsc1),PetscErrorCode,(IS{Float64},PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function ISGetIndices(arg1::IS{Float64},arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:ISGetIndices,petsc1),PetscErrorCode,(IS{Float64},Ptr{Ptr{PetscInt}}),arg1,arg2)
    return err
end

function ISRestoreIndices(arg1::IS{Float64},arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:ISRestoreIndices,petsc1),PetscErrorCode,(IS{Float64},Ptr{Ptr{PetscInt}}),arg1,arg2)
    return err
end

function ISGetTotalIndices(arg1::IS{Float64},arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:ISGetTotalIndices,petsc1),PetscErrorCode,(IS{Float64},Ptr{Ptr{PetscInt}}),arg1,arg2)
    return err
end

function ISRestoreTotalIndices(arg1::IS{Float64},arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:ISRestoreTotalIndices,petsc1),PetscErrorCode,(IS{Float64},Ptr{Ptr{PetscInt}}),arg1,arg2)
    return err
end

function ISGetNonlocalIndices(arg1::IS{Float64},arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:ISGetNonlocalIndices,petsc1),PetscErrorCode,(IS{Float64},Ptr{Ptr{PetscInt}}),arg1,arg2)
    return err
end

function ISRestoreNonlocalIndices(arg1::IS{Float64},arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:ISRestoreNonlocalIndices,petsc1),PetscErrorCode,(IS{Float64},Ptr{Ptr{PetscInt}}),arg1,arg2)
    return err
end

function ISGetNonlocalIS(arg1::IS{Float64},is::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISGetNonlocalIS,petsc1),PetscErrorCode,(IS{Float64},Ptr{IS{Float64}}),arg1,is)
    return err
end

function ISRestoreNonlocalIS(arg1::IS{Float64},is::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISRestoreNonlocalIS,petsc1),PetscErrorCode,(IS{Float64},Ptr{IS{Float64}}),arg1,is)
    return err
end

function ISGetSize(arg1::IS{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:ISGetSize,petsc1),PetscErrorCode,(IS{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function ISGetLocalSize(arg1::IS{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:ISGetLocalSize,petsc1),PetscErrorCode,(IS{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function ISInvertPermutation(arg1::IS{Float64},arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISInvertPermutation,petsc1),PetscErrorCode,(IS{Float64},PetscInt,Ptr{IS{Float64}}),arg1,arg2,arg3)
    return err
end

function ISView(arg1::IS{Float64},arg2::PetscViewer{Float64})
    err = ccall((:ISView,petsc1),PetscErrorCode,(IS{Float64},PetscViewer{Float64}),arg1,arg2)
    return err
end

function ISEqual(arg1::IS{Float64},arg2::IS{Float64},arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:ISEqual,petsc1),PetscErrorCode,(IS{Float64},IS{Float64},Ptr{PetscBool}),arg1,arg2,arg3)
    return err
end

function ISSort(arg1::IS{Float64})
    err = ccall((:ISSort,petsc1),PetscErrorCode,(IS{Float64},),arg1)
    return err
end

function ISSortRemoveDups(arg1::IS{Float64})
    err = ccall((:ISSortRemoveDups,petsc1),PetscErrorCode,(IS{Float64},),arg1)
    return err
end

function ISSorted(arg1::IS{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:ISSorted,petsc1),PetscErrorCode,(IS{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function ISDifference(arg1::IS{Float64},arg2::IS{Float64},arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISDifference,petsc1),PetscErrorCode,(IS{Float64},IS{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3)
    return err
end

function ISSum(arg1::IS{Float64},arg2::IS{Float64},arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISSum,petsc1),PetscErrorCode,(IS{Float64},IS{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3)
    return err
end

function ISExpand(arg1::IS{Float64},arg2::IS{Float64},arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISExpand,petsc1),PetscErrorCode,(IS{Float64},IS{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3)
    return err
end

function ISGetMinMax(arg1::IS{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:ISGetMinMax,petsc1),PetscErrorCode,(IS{Float64},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function ISBlockGetIndices(arg1::IS{Float64},arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:ISBlockGetIndices,petsc1),PetscErrorCode,(IS{Float64},Ptr{Ptr{PetscInt}}),arg1,arg2)
    return err
end

function ISBlockRestoreIndices(arg1::IS{Float64},arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:ISBlockRestoreIndices,petsc1),PetscErrorCode,(IS{Float64},Ptr{Ptr{PetscInt}}),arg1,arg2)
    return err
end

function ISBlockGetLocalSize(arg1::IS{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:ISBlockGetLocalSize,petsc1),PetscErrorCode,(IS{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function ISBlockGetSize(arg1::IS{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:ISBlockGetSize,petsc1),PetscErrorCode,(IS{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function ISGetBlockSize(arg1::IS{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:ISGetBlockSize,petsc1),PetscErrorCode,(IS{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function ISSetBlockSize(arg1::IS{Float64},arg2::Integer)
    err = ccall((:ISSetBlockSize,petsc1),PetscErrorCode,(IS{Float64},PetscInt),arg1,arg2)
    return err
end

function ISStrideGetInfo(arg1::IS{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:ISStrideGetInfo,petsc1),PetscErrorCode,(IS{Float64},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function ISToGeneral(arg1::IS{Float64})
    err = ccall((:ISToGeneral,petsc1),PetscErrorCode,(IS{Float64},),arg1)
    return err
end

function ISDuplicate(arg1::IS{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISDuplicate,petsc1),PetscErrorCode,(IS{Float64},Ptr{IS{Float64}}),arg1,arg2)
    return err
end

function ISCopy(arg1::IS{Float64},arg2::IS{Float64})
    err = ccall((:ISCopy,petsc1),PetscErrorCode,(IS{Float64},IS{Float64}),arg1,arg2)
    return err
end

function ISAllGather(arg1::IS{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISAllGather,petsc1),PetscErrorCode,(IS{Float64},Ptr{IS{Float64}}),arg1,arg2)
    return err
end

function ISComplement(arg1::IS{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISComplement,petsc1),PetscErrorCode,(IS{Float64},PetscInt,PetscInt,Ptr{IS{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function ISConcatenate(arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISConcatenate,petsc1),PetscErrorCode,(comm_type,PetscInt,Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1.val,arg2,arg3,arg4)
    return err
end

function ISListToPair(arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg5::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISListToPair,petsc1),PetscErrorCode,(comm_type,PetscInt,Ptr{IS{Float64}},Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1.val,arg2,arg3,arg4,arg5)
    return err
end

function ISPairToList(arg1::IS{Float64},arg2::IS{Float64},arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    err = ccall((:ISPairToList,petsc1),PetscErrorCode,(IS{Float64},IS{Float64},Ptr{PetscInt},Ptr{Ptr{IS{Float64}}}),arg1,arg2,arg3,arg4)
    return err
end

function ISEmbed(arg1::IS{Float64},arg2::IS{Float64},arg3::PetscBool,arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISEmbed,petsc1),PetscErrorCode,(IS{Float64},IS{Float64},PetscBool,Ptr{IS{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function ISSortPermutation(arg1::IS{Float64},arg2::PetscBool,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISSortPermutation,petsc1),PetscErrorCode,(IS{Float64},PetscBool,Ptr{IS{Float64}}),arg1,arg2,arg3)
    return err
end

function ISOnComm(arg1::IS{Float64},arg2::MPI_Comm,PetscCopyMode::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISOnComm,petsc1),PetscErrorCode,(IS{Float64},comm_type,Cint,Ptr{IS{Float64}}),arg1,arg2.val,PetscCopyMode,arg3)
    return err
end

function ISLocalToGlobalMappingCreate(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),PetscCopyMode::Integer,arg5::Union(Ptr{ISLocalToGlobalMapping{Float64}},StridedArray{ISLocalToGlobalMapping{Float64}},Ptr{Void}))
    err = ccall((:ISLocalToGlobalMappingCreate,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,Ptr{PetscInt},Cint,Ptr{ISLocalToGlobalMapping{Float64}}),arg1.val,arg2,arg3,arg4,PetscCopyMode,arg5)
    return err
end

function ISLocalToGlobalMappingCreateIS(arg1::IS{Float64},arg2::Union(Ptr{ISLocalToGlobalMapping{Float64}},StridedArray{ISLocalToGlobalMapping{Float64}},Ptr{Void}))
    err = ccall((:ISLocalToGlobalMappingCreateIS,petsc1),PetscErrorCode,(IS{Float64},Ptr{ISLocalToGlobalMapping{Float64}}),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function ISLocalToGlobalMappingCreateSF(arg1::PetscSF,arg2::Integer,arg3::Union(Ptr{ISLocalToGlobalMapping{Float64}},StridedArray{ISLocalToGlobalMapping{Float64}},Ptr{Void}))
    ccall((:ISLocalToGlobalMappingCreateSF,petsc1),PetscErrorCode,(PetscSF,PetscInt,Ptr{ISLocalToGlobalMapping{Float64}}),arg1,arg2,arg3)
end 
=#
function ISLocalToGlobalMappingView(arg1::ISLocalToGlobalMapping{Float64},arg2::PetscViewer{Float64})
    err = ccall((:ISLocalToGlobalMappingView,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},PetscViewer{Float64}),arg1,arg2)
    return err
end

function ISLocalToGlobalMappingDestroy(arg1::Union(Ptr{ISLocalToGlobalMapping{Float64}},StridedArray{ISLocalToGlobalMapping{Float64}},Ptr{Void}))
    err = ccall((:ISLocalToGlobalMappingDestroy,petsc1),PetscErrorCode,(Ptr{ISLocalToGlobalMapping{Float64}},),arg1)
    return err
end

function ISLocalToGlobalMappingApply(arg1::ISLocalToGlobalMapping{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:ISLocalToGlobalMappingApply,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
    return err
end

function ISLocalToGlobalMappingApplyBlock(arg1::ISLocalToGlobalMapping{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:ISLocalToGlobalMappingApplyBlock,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
    return err
end

function ISLocalToGlobalMappingApplyIS(arg1::ISLocalToGlobalMapping{Float64},arg2::IS{Float64},arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISLocalToGlobalMappingApplyIS,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},IS{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3)
    return err
end

function ISGlobalToLocalMappingApply(arg1::ISLocalToGlobalMapping{Float64},arg2::ISGlobalToLocalMappingType,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:ISGlobalToLocalMappingApply,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},ISGlobalToLocalMappingType,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function ISGlobalToLocalMappingApplyBlock(arg1::ISLocalToGlobalMapping{Float64},arg2::ISGlobalToLocalMappingType,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:ISGlobalToLocalMappingApplyBlock,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},ISGlobalToLocalMappingType,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function ISGlobalToLocalMappingApplyIS(arg1::ISLocalToGlobalMapping{Float64},arg2::ISGlobalToLocalMappingType,arg3::IS{Float64},arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISGlobalToLocalMappingApplyIS,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},ISGlobalToLocalMappingType,IS{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function ISLocalToGlobalMappingGetSize(arg1::ISLocalToGlobalMapping{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:ISLocalToGlobalMappingGetSize,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function ISLocalToGlobalMappingGetInfo(arg1::ISLocalToGlobalMapping{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg4::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg5::Union(Ptr{Ptr{Ptr{PetscInt}}},StridedArray{Ptr{Ptr{PetscInt}}},Ptr{Void}))
    err = ccall((:ISLocalToGlobalMappingGetInfo,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{Ptr{Ptr{PetscInt}}}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function ISLocalToGlobalMappingRestoreInfo(arg1::ISLocalToGlobalMapping{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg4::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg5::Union(Ptr{Ptr{Ptr{PetscInt}}},StridedArray{Ptr{Ptr{PetscInt}}},Ptr{Void}))
    err = ccall((:ISLocalToGlobalMappingRestoreInfo,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{Ptr{Ptr{PetscInt}}}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function ISLocalToGlobalMappingGetBlockInfo(arg1::ISLocalToGlobalMapping{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg4::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg5::Union(Ptr{Ptr{Ptr{PetscInt}}},StridedArray{Ptr{Ptr{PetscInt}}},Ptr{Void}))
    err = ccall((:ISLocalToGlobalMappingGetBlockInfo,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{Ptr{Ptr{PetscInt}}}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function ISLocalToGlobalMappingRestoreBlockInfo(arg1::ISLocalToGlobalMapping{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg4::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg5::Union(Ptr{Ptr{Ptr{PetscInt}}},StridedArray{Ptr{Ptr{PetscInt}}},Ptr{Void}))
    err = ccall((:ISLocalToGlobalMappingRestoreBlockInfo,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{Ptr{Ptr{PetscInt}}}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function ISLocalToGlobalMappingGetIndices(arg1::ISLocalToGlobalMapping{Float64},arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:ISLocalToGlobalMappingGetIndices,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},Ptr{Ptr{PetscInt}}),arg1,arg2)
    return err
end

function ISLocalToGlobalMappingRestoreIndices(arg1::ISLocalToGlobalMapping{Float64},arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:ISLocalToGlobalMappingRestoreIndices,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},Ptr{Ptr{PetscInt}}),arg1,arg2)
    return err
end

function ISLocalToGlobalMappingGetBlockIndices(arg1::ISLocalToGlobalMapping{Float64},arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:ISLocalToGlobalMappingGetBlockIndices,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},Ptr{Ptr{PetscInt}}),arg1,arg2)
    return err
end

function ISLocalToGlobalMappingRestoreBlockIndices(arg1::ISLocalToGlobalMapping{Float64},arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:ISLocalToGlobalMappingRestoreBlockIndices,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},Ptr{Ptr{PetscInt}}),arg1,arg2)
    return err
end

function ISLocalToGlobalMappingConcatenate(arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{ISLocalToGlobalMapping{Float64}},StridedArray{ISLocalToGlobalMapping{Float64}},Ptr{Void}),arg4::Union(Ptr{ISLocalToGlobalMapping{Float64}},StridedArray{ISLocalToGlobalMapping{Float64}},Ptr{Void}))
    err = ccall((:ISLocalToGlobalMappingConcatenate,petsc1),PetscErrorCode,(comm_type,PetscInt,Ptr{ISLocalToGlobalMapping{Float64}},Ptr{ISLocalToGlobalMapping{Float64}}),arg1.val,arg2,arg3,arg4)
    return err
end

function ISG2LMapApply(arg1::ISLocalToGlobalMapping{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:ISG2LMapApply,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
    return err
end

function ISLocalToGlobalMappingGetBlockSize(arg1::ISLocalToGlobalMapping{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:ISLocalToGlobalMappingGetBlockSize,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function ISAllGatherColors(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}))
    err = ccall((:ISAllGatherColors,petsc1),PetscErrorCode,(comm_type,PetscInt,Ptr{Cint},Ptr{PetscInt},Ptr{Ptr{Cint}}),arg1.val,arg2,arg3,arg4,arg5)
    return err
end

function ISColoringCreate(arg1::MPI_Comm,arg2::Integer,arg3::Integer,ISColoringValue::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),PetscCopyMode::Integer,arg4::Union(Ptr{ISColoring{Float64}},StridedArray{ISColoring{Float64}},Ptr{Void}))
    err = ccall((:ISColoringCreate,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,Ptr{Cint},Cint,Ptr{ISColoring{Float64}}),arg1.val,arg2,arg3,ISColoringValue,PetscCopyMode,arg4)
    return err
end

function ISColoringDestroy(arg1::Union(Ptr{ISColoring{Float64}},StridedArray{ISColoring{Float64}},Ptr{Void}))
    err = ccall((:ISColoringDestroy,petsc1),PetscErrorCode,(Ptr{ISColoring{Float64}},),arg1)
    return err
end

function ISColoringView(arg1::ISColoring{Float64},arg2::PetscViewer{Float64})
    err = ccall((:ISColoringView,petsc1),PetscErrorCode,(ISColoring{Float64},PetscViewer{Float64}),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function ISColoringViewFromOptions(arg1::ISColoring{Float64},arg2::PetscObject,arg3::Union(ByteString,Symbol))
    ccall((:ISColoringViewFromOptions,petsc1),PetscErrorCode,(ISColoring{Float64},PetscObject,Cstring),arg1,arg2,arg3)
end 
=#
function ISColoringGetIS(arg1::ISColoring{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    err = ccall((:ISColoringGetIS,petsc1),PetscErrorCode,(ISColoring{Float64},Ptr{PetscInt},Ptr{Ptr{IS{Float64}}}),arg1,arg2,arg3)
    return err
end

function ISColoringRestoreIS(arg1::ISColoring{Float64},arg2::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    err = ccall((:ISColoringRestoreIS,petsc1),PetscErrorCode,(ISColoring{Float64},Ptr{Ptr{IS{Float64}}}),arg1,arg2)
    return err
end

function ISColoringReference(arg1::ISColoring{Float64})
    err = ccall((:ISColoringReference,petsc1),PetscErrorCode,(ISColoring{Float64},),arg1)
    return err
end

function ISColoringSetType(arg1::ISColoring{Float64},arg2::ISColoringType)
    err = ccall((:ISColoringSetType,petsc1),PetscErrorCode,(ISColoring{Float64},ISColoringType),arg1,arg2)
    return err
end

function ISPartitioningToNumbering(arg1::IS{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISPartitioningToNumbering,petsc1),PetscErrorCode,(IS{Float64},Ptr{IS{Float64}}),arg1,arg2)
    return err
end

function ISPartitioningCount(arg1::IS{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:ISPartitioningCount,petsc1),PetscErrorCode,(IS{Float64},PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function ISCompressIndicesGeneral(arg1::Integer,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg6::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISCompressIndicesGeneral,petsc1),PetscErrorCode,(PetscInt,PetscInt,PetscInt,PetscInt,Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function ISCompressIndicesSorted(arg1::Integer,arg2::Integer,arg3::Integer,arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg5::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISCompressIndicesSorted,petsc1),PetscErrorCode,(PetscInt,PetscInt,PetscInt,Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function ISExpandIndicesGeneral(arg1::Integer,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg6::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISExpandIndicesGeneral,petsc1),PetscErrorCode,(PetscInt,PetscInt,PetscInt,PetscInt,Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function PetscLayoutSetUp(arg1::PetscLayout{Float64})
    err = ccall((:PetscLayoutSetUp,petsc1),PetscErrorCode,(PetscLayout{Float64},),arg1)
    return err
end

function PetscLayoutDestroy(arg1::Union(Ptr{PetscLayout{Float64}},StridedArray{PetscLayout{Float64}},Ptr{Void}))
    err = ccall((:PetscLayoutDestroy,petsc1),PetscErrorCode,(Ptr{PetscLayout{Float64}{Float64}},),arg1)
    return err
end

function PetscLayoutDuplicate(arg1::PetscLayout{Float64},arg2::Union(Ptr{PetscLayout{Float64}},StridedArray{PetscLayout{Float64}},Ptr{Void}))
    err = ccall((:PetscLayoutDuplicate,petsc1),PetscErrorCode,(PetscLayout{Float64},Ptr{PetscLayout{Float64}{Float64}}),arg1,arg2)
    return err
end

function PetscLayoutReference(arg1::PetscLayout{Float64},arg2::Union(Ptr{PetscLayout{Float64}},StridedArray{PetscLayout{Float64}},Ptr{Void}))
    err = ccall((:PetscLayoutReference,petsc1),PetscErrorCode,(PetscLayout{Float64},Ptr{PetscLayout{Float64}{Float64}}),arg1,arg2)
    return err
end

function PetscLayoutSetLocalSize(arg1::PetscLayout{Float64},arg2::Integer)
    err = ccall((:PetscLayoutSetLocalSize,petsc1),PetscErrorCode,(PetscLayout{Float64},PetscInt),arg1,arg2)
    return err
end

function PetscLayoutGetLocalSize(arg1::PetscLayout{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PetscLayoutGetLocalSize,petsc1),PetscErrorCode,(PetscLayout{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function PetscLayoutSetSize(arg1::PetscLayout{Float64},arg2::Integer)
    err = ccall((:PetscLayoutSetSize,petsc1),PetscErrorCode,(PetscLayout{Float64},PetscInt),arg1,arg2)
    return err
end

function PetscLayoutGetSize(arg1::PetscLayout{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PetscLayoutGetSize,petsc1),PetscErrorCode,(PetscLayout{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function PetscLayoutSetBlockSize(arg1::PetscLayout{Float64},arg2::Integer)
    err = ccall((:PetscLayoutSetBlockSize,petsc1),PetscErrorCode,(PetscLayout{Float64},PetscInt),arg1,arg2)
    return err
end

function PetscLayoutGetBlockSize(arg1::PetscLayout{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PetscLayoutGetBlockSize,petsc1),PetscErrorCode,(PetscLayout{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function PetscLayoutGetRange(arg1::PetscLayout{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PetscLayoutGetRange,petsc1),PetscErrorCode,(PetscLayout{Float64},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function PetscLayoutGetRanges(arg1::PetscLayout{Float64},arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:PetscLayoutGetRanges,petsc1),PetscErrorCode,(PetscLayout{Float64},Ptr{Ptr{PetscInt}}),arg1,arg2)
    return err
end

function PetscLayoutSetISLocalToGlobalMapping(arg1::PetscLayout{Float64},arg2::ISLocalToGlobalMapping{Float64})
    err = ccall((:PetscLayoutSetISLocalToGlobalMapping,petsc1),PetscErrorCode,(PetscLayout{Float64},ISLocalToGlobalMapping{Float64}),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function PetscSFSetGraphLayout(arg1::PetscSF,arg2::PetscLayout{Float64},arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),PetscCopyMode::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscSFSetGraphLayout,petsc1),PetscErrorCode,(PetscSF,PetscLayout{Float64},PetscInt,Ptr{PetscInt},Cint,Ptr{PetscInt}),arg1,arg2,arg3,arg4,PetscCopyMode,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}))
    ccall((:PetscSectionCreate,petsc1),PetscErrorCode,(comm_type,Ptr{PetscSection}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionClone(arg0::Type{Float64},arg1::PetscSection,arg2::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}))
    ccall((:PetscSectionClone,petsc1),PetscErrorCode,(PetscSection,Ptr{PetscSection}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionCopy(arg0::Type{Float64},arg1::PetscSection,arg2::PetscSection)
    ccall((:PetscSectionCopy,petsc1),PetscErrorCode,(PetscSection,PetscSection),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetNumFields(arg0::Type{Float64},arg1::PetscSection,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscSectionGetNumFields,petsc1),PetscErrorCode,(PetscSection,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetNumFields(arg0::Type{Float64},arg1::PetscSection,arg2::Integer)
    ccall((:PetscSectionSetNumFields,petsc1),PetscErrorCode,(PetscSection,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetFieldName(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:PetscSectionGetFieldName,petsc1),PetscErrorCode,(PetscSection,PetscInt,Ptr{Ptr{Uint8}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetFieldName(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Union(ByteString,Symbol))
    ccall((:PetscSectionSetFieldName,petsc1),PetscErrorCode,(PetscSection,PetscInt,Cstring),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetFieldComponents(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscSectionGetFieldComponents,petsc1),PetscErrorCode,(PetscSection,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetFieldComponents(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer)
    ccall((:PetscSectionSetFieldComponents,petsc1),PetscErrorCode,(PetscSection,PetscInt,PetscInt),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetChart(arg0::Type{Float64},arg1::PetscSection,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscSectionGetChart,petsc1),PetscErrorCode,(PetscSection,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetChart(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer)
    ccall((:PetscSectionSetChart,petsc1),PetscErrorCode,(PetscSection,PetscInt,PetscInt),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetPermutation(arg1::PetscSection,arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:PetscSectionGetPermutation,petsc1),PetscErrorCode,(PetscSection,Ptr{IS{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetPermutation(arg1::PetscSection,arg2::IS{Float64})
    ccall((:PetscSectionSetPermutation,petsc1),PetscErrorCode,(PetscSection,IS{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscSectionGetDof,petsc1),PetscErrorCode,(PetscSection,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer)
    ccall((:PetscSectionSetDof,petsc1),PetscErrorCode,(PetscSection,PetscInt,PetscInt),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionAddDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer)
    ccall((:PetscSectionAddDof,petsc1),PetscErrorCode,(PetscSection,PetscInt,PetscInt),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetFieldDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscSectionGetFieldDof,petsc1),PetscErrorCode,(PetscSection,PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetFieldDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:PetscSectionSetFieldDof,petsc1),PetscErrorCode,(PetscSection,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionAddFieldDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:PetscSectionAddFieldDof,petsc1),PetscErrorCode,(PetscSection,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionHasConstraints(arg0::Type{Float64},arg1::PetscSection,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscSectionHasConstraints,petsc1),PetscErrorCode,(PetscSection,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetConstraintDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscSectionGetConstraintDof,petsc1),PetscErrorCode,(PetscSection,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetConstraintDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer)
    ccall((:PetscSectionSetConstraintDof,petsc1),PetscErrorCode,(PetscSection,PetscInt,PetscInt),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionAddConstraintDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer)
    ccall((:PetscSectionAddConstraintDof,petsc1),PetscErrorCode,(PetscSection,PetscInt,PetscInt),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetFieldConstraintDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscSectionGetFieldConstraintDof,petsc1),PetscErrorCode,(PetscSection,PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetFieldConstraintDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:PetscSectionSetFieldConstraintDof,petsc1),PetscErrorCode,(PetscSection,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionAddFieldConstraintDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:PetscSectionAddFieldConstraintDof,petsc1),PetscErrorCode,(PetscSection,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetConstraintIndices(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:PetscSectionGetConstraintIndices,petsc1),PetscErrorCode,(PetscSection,PetscInt,Ptr{Ptr{PetscInt}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetConstraintIndices(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscSectionSetConstraintIndices,petsc1),PetscErrorCode,(PetscSection,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetFieldConstraintIndices(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:PetscSectionGetFieldConstraintIndices,petsc1),PetscErrorCode,(PetscSection,PetscInt,PetscInt,Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetFieldConstraintIndices(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscSectionSetFieldConstraintIndices,petsc1),PetscErrorCode,(PetscSection,PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetUpBC(arg0::Type{Float64},arg1::PetscSection)
    ccall((:PetscSectionSetUpBC,petsc1),PetscErrorCode,(PetscSection,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetUp(arg0::Type{Float64},arg1::PetscSection)
    ccall((:PetscSectionSetUp,petsc1),PetscErrorCode,(PetscSection,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetMaxDof(arg0::Type{Float64},arg1::PetscSection,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscSectionGetMaxDof,petsc1),PetscErrorCode,(PetscSection,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetStorageSize(arg0::Type{Float64},arg1::PetscSection,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscSectionGetStorageSize,petsc1),PetscErrorCode,(PetscSection,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetConstrainedStorageSize(arg0::Type{Float64},arg1::PetscSection,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscSectionGetConstrainedStorageSize,petsc1),PetscErrorCode,(PetscSection,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetOffset(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscSectionGetOffset,petsc1),PetscErrorCode,(PetscSection,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetOffset(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer)
    ccall((:PetscSectionSetOffset,petsc1),PetscErrorCode,(PetscSection,PetscInt,PetscInt),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetFieldOffset(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscSectionGetFieldOffset,petsc1),PetscErrorCode,(PetscSection,PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetFieldOffset(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:PetscSectionSetFieldOffset,petsc1),PetscErrorCode,(PetscSection,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetOffsetRange(arg0::Type{Float64},arg1::PetscSection,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscSectionGetOffsetRange,petsc1),PetscErrorCode,(PetscSection,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionView(arg1::PetscSection,arg2::PetscViewer{Float64})
    ccall((:PetscSectionView,petsc1),PetscErrorCode,(PetscSection,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}))
    ccall((:PetscSectionDestroy,petsc1),PetscErrorCode,(Ptr{PetscSection},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionCreateGlobalSection(arg0::Type{Float64},arg1::PetscSection,arg2::PetscSF,arg3::PetscBool,arg4::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}))
    ccall((:PetscSectionCreateGlobalSection,petsc1),PetscErrorCode,(PetscSection,PetscSF,PetscBool,Ptr{PetscSection}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionCreateGlobalSectionCensored(arg0::Type{Float64},arg1::PetscSection,arg2::PetscSF,arg3::PetscBool,arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}))
    ccall((:PetscSectionCreateGlobalSectionCensored,petsc1),PetscErrorCode,(PetscSection,PetscSF,PetscBool,PetscInt,Ptr{PetscInt},Ptr{PetscSection}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionCreateSubsection(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}))
    ccall((:PetscSectionCreateSubsection,petsc1),PetscErrorCode,(PetscSection,PetscInt,Ptr{PetscInt},Ptr{PetscSection}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionCreateSubmeshSection(arg1::PetscSection,arg2::IS{Float64},arg3::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}))
    ccall((:PetscSectionCreateSubmeshSection,petsc1),PetscErrorCode,(PetscSection,IS{Float64},Ptr{PetscSection}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetPointLayout(arg1::MPI_Comm,arg2::PetscSection,arg3::Union(Ptr{PetscLayout{Float64}},StridedArray{PetscLayout{Float64}},Ptr{Void}))
    ccall((:PetscSectionGetPointLayout,petsc1),PetscErrorCode,(comm_type,PetscSection,Ptr{PetscLayout{Float64}{Float64}}),arg1.val,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetValueLayout(arg1::MPI_Comm,arg2::PetscSection,arg3::Union(Ptr{PetscLayout{Float64}},StridedArray{PetscLayout{Float64}},Ptr{Void}))
    ccall((:PetscSectionGetValueLayout,petsc1),PetscErrorCode,(comm_type,PetscSection,Ptr{PetscLayout{Float64}{Float64}}),arg1.val,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionPermute(arg1::PetscSection,arg2::IS{Float64},arg3::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}))
    ccall((:PetscSectionPermute,petsc1),PetscErrorCode,(PetscSection,IS{Float64},Ptr{PetscSection}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetField(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}))
    ccall((:PetscSectionGetField,petsc1),PetscErrorCode,(PetscSection,PetscInt,Ptr{PetscSection}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetClosureIndex(arg1::PetscSection,arg2::PetscObject,arg3::PetscSection,arg4::IS{Float64})
    ccall((:PetscSectionSetClosureIndex,petsc1),PetscErrorCode,(PetscSection,PetscObject,PetscSection,IS{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetClosureIndex(arg1::PetscSection,arg2::PetscObject,arg3::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}),arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:PetscSectionGetClosureIndex,petsc1),PetscErrorCode,(PetscSection,PetscObject,Ptr{PetscSection},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSFConvertPartition(arg1::PetscSF,arg2::PetscSection,arg3::IS{Float64},arg4::Union(Ptr{ISLocalToGlobalMapping{Float64}},StridedArray{ISLocalToGlobalMapping{Float64}},Ptr{Void}),arg5::Union(Ptr{PetscSF},StridedArray{PetscSF},Ptr{Void}))
    ccall((:PetscSFConvertPartition,petsc1),PetscErrorCode,(PetscSF,PetscSection,IS{Float64},Ptr{ISLocalToGlobalMapping{Float64}},Ptr{PetscSF}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSFCreateRemoteOffsets(arg0::Type{Float64},arg1::PetscSF,arg2::PetscSection,arg3::PetscSection,arg4::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:PetscSFCreateRemoteOffsets,petsc1),PetscErrorCode,(PetscSF,PetscSection,PetscSection,Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSFDistributeSection(arg0::Type{Float64},arg1::PetscSF,arg2::PetscSection,arg3::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg4::PetscSection)
    ccall((:PetscSFDistributeSection,petsc1),PetscErrorCode,(PetscSF,PetscSection,Ptr{Ptr{PetscInt}},PetscSection),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSFCreateSectionSF(arg0::Type{Float64},arg1::PetscSF,arg2::PetscSection,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::PetscSection,arg5::Union(Ptr{PetscSF},StridedArray{PetscSF},Ptr{Void}))
    ccall((:PetscSFCreateSectionSF,petsc1),PetscErrorCode,(PetscSF,PetscSection,Ptr{PetscInt},PetscSection,Ptr{PetscSF}),arg1,arg2,arg3,arg4,arg5)
end 
=#
function PetscViewerInitializePackage(arg0::Type{Float64})
    err = ccall((:PetscViewerInitializePackage,petsc1),PetscErrorCode,())
    return err
end

function PetscViewerRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscViewerRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

function PetscViewerCreate(arg1::MPI_Comm,arg2::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    err = ccall((:PetscViewerCreate,petsc1),PetscErrorCode,(comm_type,Ptr{PetscViewer{Float64}}),arg1.val,arg2)
    return err
end

function PetscViewerSetFromOptions(arg1::PetscViewer{Float64})
    err = ccall((:PetscViewerSetFromOptions,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
    return err
end

#= skipping function with undefined symbols: 
 function PetscViewerASCIIOpenWithFILE(arg1::MPI_Comm,arg2::Union(Ptr{FILE},StridedArray{FILE},Ptr{Void}),arg3::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    ccall((:PetscViewerASCIIOpenWithFILE,petsc1),PetscErrorCode,(comm_type,Ptr{FILE},Ptr{PetscViewer{Float64}}),arg1.val,arg2,arg3)
end 
=#
function PetscViewerASCIIOpen(arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    err = ccall((:PetscViewerASCIIOpen,petsc1),PetscErrorCode,(comm_type,Cstring,Ptr{PetscViewer{Float64}}),arg1.val,arg2,arg3)
    return err
end

#= skipping function with undefined symbols: 
 function PetscViewerASCIISetFILE(arg1::PetscViewer{Float64},arg2::Union(Ptr{FILE},StridedArray{FILE},Ptr{Void}))
    ccall((:PetscViewerASCIISetFILE,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{FILE}),arg1,arg2)
end 
=#
function PetscViewerBinaryOpen(arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::PetscFileMode,arg4::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    err = ccall((:PetscViewerBinaryOpen,petsc1),PetscErrorCode,(comm_type,Cstring,PetscFileMode,Ptr{PetscViewer{Float64}}),arg1.val,arg2,arg3,arg4)
    return err
end

function PetscViewerBinaryGetFlowControl(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PetscViewerBinaryGetFlowControl,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function PetscViewerBinarySetFlowControl(arg1::PetscViewer{Float64},arg2::Integer)
    err = ccall((:PetscViewerBinarySetFlowControl,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscInt),arg1,arg2)
    return err
end

function PetscViewerBinarySetUseMPIIO(arg1::PetscViewer{Float64},arg2::PetscBool)
    err = ccall((:PetscViewerBinarySetUseMPIIO,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscBool),arg1,arg2)
    return err
end

function PetscViewerBinaryGetUseMPIIO(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscViewerBinaryGetUseMPIIO,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function PetscViewerSocketOpen(arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Integer,arg4::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    err = ccall((:PetscViewerSocketOpen,petsc1),PetscErrorCode,(comm_type,Cstring,Cint,Ptr{PetscViewer{Float64}}),arg1.val,arg2,arg3,arg4)
    return err
end

function PetscViewerStringOpen(arg1::MPI_Comm,arg2::Union(ByteString,Symbol),size_t::Integer,arg3::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    err = ccall((:PetscViewerStringOpen,petsc1),PetscErrorCode,(comm_type,Cstring,Cint,Ptr{PetscViewer{Float64}}),arg1.val,arg2,size_t,arg3)
    return err
end

function PetscViewerDrawOpen(arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    err = ccall((:PetscViewerDrawOpen,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint,Cint,Cint,Cint,Ptr{PetscViewer{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
    return err
end

function PetscViewerDrawSetDrawType(arg1::PetscViewer{Float64},arg2::PetscDrawType)
    err = ccall((:PetscViewerDrawSetDrawType,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring),arg1,arg2)
    return err
end

function PetscViewerMathematicaOpen(arg1::MPI_Comm,arg2::Integer,arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    err = ccall((:PetscViewerMathematicaOpen,petsc1),PetscErrorCode,(comm_type,Cint,Cstring,Cstring,Ptr{PetscViewer{Float64}}),arg1.val,arg2,arg3,arg4,arg5)
    return err
end

function PetscViewerSiloOpen(arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    err = ccall((:PetscViewerSiloOpen,petsc1),PetscErrorCode,(comm_type,Cstring,Ptr{PetscViewer{Float64}}),arg1.val,arg2,arg3)
    return err
end

function PetscViewerMatlabOpen(arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::PetscFileMode,arg4::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    err = ccall((:PetscViewerMatlabOpen,petsc1),PetscErrorCode,(comm_type,Cstring,PetscFileMode,Ptr{PetscViewer{Float64}}),arg1.val,arg2,arg3,arg4)
    return err
end

function PetscViewerGetType(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscViewerType},StridedArray{PetscViewerType},Ptr{Void}))
    (arg2_,tmp) = symbol_get_before(arg2)
    err = ccall((:PetscViewerGetType,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Ptr{Uint8}}),arg1,arg2_)
    symbol_get_after(arg2_,arg2)
    return err
end

function PetscViewerSetType(arg1::PetscViewer{Float64},arg2::PetscViewerType)
    err = ccall((:PetscViewerSetType,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring),arg1,arg2)
    return err
end

function PetscViewerDestroy(arg1::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    err = ccall((:PetscViewerDestroy,petsc1),PetscErrorCode,(Ptr{PetscViewer{Float64}},),arg1)
    return err
end

function PetscViewerGetSingleton(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    err = ccall((:PetscViewerGetSingleton,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{PetscViewer{Float64}}),arg1,arg2)
    return err
end

function PetscViewerRestoreSingleton(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    err = ccall((:PetscViewerRestoreSingleton,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{PetscViewer{Float64}}),arg1,arg2)
    return err
end

function PetscViewerGetSubcomm(arg1::PetscViewer{Float64},arg2::MPI_Comm,arg3::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    err = ccall((:PetscViewerGetSubcomm,petsc1),PetscErrorCode,(PetscViewer{Float64},comm_type,Ptr{PetscViewer{Float64}}),arg1,arg2.val,arg3)
    return err
end

function PetscViewerRestoreSubcomm(arg1::PetscViewer{Float64},arg2::MPI_Comm,arg3::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    err = ccall((:PetscViewerRestoreSubcomm,petsc1),PetscErrorCode,(PetscViewer{Float64},comm_type,Ptr{PetscViewer{Float64}}),arg1,arg2.val,arg3)
    return err
end

function PetscViewerSetUp(arg1::PetscViewer{Float64})
    err = ccall((:PetscViewerSetUp,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
    return err
end

function PetscViewerView(arg1::PetscViewer{Float64},arg2::PetscViewer{Float64})
    err = ccall((:PetscViewerView,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscViewer{Float64}),arg1,arg2)
    return err
end

function PetscViewerAppendOptionsPrefix(arg1::PetscViewer{Float64},arg2::Union(ByteString,Symbol))
    err = ccall((:PetscViewerAppendOptionsPrefix,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring),arg1,arg2)
    return err
end

function PetscViewerGetOptionsPrefix(arg1::PetscViewer{Float64},arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:PetscViewerGetOptionsPrefix,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Ptr{Uint8}}),arg1,arg2)
    return err
end

function PetscViewerSetFormat(arg1::PetscViewer{Float64},arg2::PetscViewerFormat)
    err = ccall((:PetscViewerSetFormat,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscViewerFormat),arg1,arg2)
    return err
end

function PetscViewerPushFormat(arg1::PetscViewer{Float64},arg2::PetscViewerFormat)
    err = ccall((:PetscViewerPushFormat,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscViewerFormat),arg1,arg2)
    return err
end

function PetscViewerPopFormat(arg1::PetscViewer{Float64})
    err = ccall((:PetscViewerPopFormat,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
    return err
end

function PetscViewerGetFormat(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscViewerFormat},StridedArray{PetscViewerFormat},Ptr{Void}))
    err = ccall((:PetscViewerGetFormat,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{PetscViewerFormat}),arg1,arg2)
    return err
end

function PetscViewerFlush(arg1::PetscViewer{Float64})
    err = ccall((:PetscViewerFlush,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
    return err
end

function PetscOptionsGetViewer(arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}),arg5::Union(Ptr{PetscViewerFormat},StridedArray{PetscViewerFormat},Ptr{Void}),arg6::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscOptionsGetViewer,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Ptr{PetscViewer{Float64}},Ptr{PetscViewerFormat},Ptr{PetscBool}),arg1.val,arg2,arg3,arg4,arg5,arg6)
    return err
end

#= skipping function with undefined symbols: 
 function PetscOptionsViewer_Private(arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}),arg6::Union(Ptr{PetscViewerFormat},StridedArray{PetscViewerFormat},Ptr{Void}),arg7::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsViewer_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{PetscViewer{Float64}},Ptr{PetscViewerFormat},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscViewerASCIIGetPointer(arg1::PetscViewer{Float64},arg2::Union(Ptr{Ptr{FILE}},StridedArray{Ptr{FILE}},Ptr{Void}))
    ccall((:PetscViewerASCIIGetPointer,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Ptr{FILE}}),arg1,arg2)
end 
=#
function PetscViewerFileGetMode(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscFileMode},StridedArray{PetscFileMode},Ptr{Void}))
    err = ccall((:PetscViewerFileGetMode,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{PetscFileMode}),arg1,arg2)
    return err
end

function PetscViewerFileSetMode(arg1::PetscViewer{Float64},arg2::PetscFileMode)
    err = ccall((:PetscViewerFileSetMode,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscFileMode),arg1,arg2)
    return err
end

function PetscViewerRead(arg1::PetscViewer{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::PetscDataType)
    err = ccall((:PetscViewerRead,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Void},PetscInt,Ptr{PetscInt},PetscDataType),arg1,arg2,arg3,arg4,arg5)
    return err
end

function PetscViewerASCIISynchronizedAllow(arg1::PetscViewer{Float64},arg2::PetscBool)
    err = ccall((:PetscViewerASCIISynchronizedAllow,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscBool),arg1,arg2)
    return err
end

function PetscViewerASCIIPushTab(arg1::PetscViewer{Float64})
    err = ccall((:PetscViewerASCIIPushTab,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
    return err
end

function PetscViewerASCIIPopTab(arg1::PetscViewer{Float64})
    err = ccall((:PetscViewerASCIIPopTab,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
    return err
end

function PetscViewerASCIIUseTabs(arg1::PetscViewer{Float64},arg2::PetscBool)
    err = ccall((:PetscViewerASCIIUseTabs,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscBool),arg1,arg2)
    return err
end

function PetscViewerASCIISetTab(arg1::PetscViewer{Float64},arg2::Integer)
    err = ccall((:PetscViewerASCIISetTab,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscInt),arg1,arg2)
    return err
end

function PetscViewerASCIIGetTab(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PetscViewerASCIIGetTab,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function PetscViewerASCIIAddTab(arg1::PetscViewer{Float64},arg2::Integer)
    err = ccall((:PetscViewerASCIIAddTab,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscInt),arg1,arg2)
    return err
end

function PetscViewerASCIISubtractTab(arg1::PetscViewer{Float64},arg2::Integer)
    err = ccall((:PetscViewerASCIISubtractTab,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscInt),arg1,arg2)
    return err
end

function PetscViewerASCIIRead(arg1::PetscViewer{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::PetscDataType)
    err = ccall((:PetscViewerASCIIRead,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Void},PetscInt,Ptr{PetscInt},PetscDataType),arg1,arg2,arg3,arg4,arg5)
    return err
end

function PetscViewerBinaryGetDescriptor(arg1::PetscViewer{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:PetscViewerBinaryGetDescriptor,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Cint}),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function PetscViewerBinaryGetInfoPointer(arg1::PetscViewer{Float64},arg2::Union(Ptr{Ptr{FILE}},StridedArray{Ptr{FILE}},Ptr{Void}))
    ccall((:PetscViewerBinaryGetInfoPointer,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Ptr{FILE}}),arg1,arg2)
end 
=#
function PetscViewerBinaryRead(arg1::PetscViewer{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::PetscDataType)
    err = ccall((:PetscViewerBinaryRead,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Void},PetscInt,Ptr{PetscInt},PetscDataType),arg1,arg2,arg3,arg4,arg5)
    return err
end

function PetscViewerBinaryWrite(arg1::PetscViewer{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Integer,arg4::PetscDataType,arg5::PetscBool)
    err = ccall((:PetscViewerBinaryWrite,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Void},PetscInt,PetscDataType,PetscBool),arg1,arg2,arg3,arg4,arg5)
    return err
end

function PetscViewerStringSetString(arg1::PetscViewer{Float64},arg2::Union(ByteString,Symbol),arg3::Integer)
    err = ccall((:PetscViewerStringSetString,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring,PetscInt),arg1,arg2,arg3)
    return err
end

function PetscViewerDrawClear(arg1::PetscViewer{Float64})
    err = ccall((:PetscViewerDrawClear,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
    return err
end

function PetscViewerDrawSetHold(arg1::PetscViewer{Float64},arg2::PetscBool)
    err = ccall((:PetscViewerDrawSetHold,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscBool),arg1,arg2)
    return err
end

function PetscViewerDrawGetHold(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscViewerDrawGetHold,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function PetscViewerDrawSetPause(arg1::PetscViewer{Float64},PetscReal::Integer)
    err = ccall((:PetscViewerDrawSetPause,petsc1),PetscErrorCode,(PetscViewer{Float64},Cint),arg1,PetscReal)
    return err
end

function PetscViewerDrawGetPause(arg1::PetscViewer{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:PetscViewerDrawGetPause,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Cint}),arg1,arg2)
    return err
end

function PetscViewerDrawSetInfo(arg1::PetscViewer{Float64},arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer)
    err = ccall((:PetscViewerDrawSetInfo,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring,Cstring,Cint,Cint,Cint,Cint),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
    return err
end

function PetscViewerDrawResize(arg1::PetscViewer{Float64},arg2::Integer,arg3::Integer)
    err = ccall((:PetscViewerDrawResize,petsc1),PetscErrorCode,(PetscViewer{Float64},Cint,Cint),arg1,arg2,arg3)
    return err
end

function PetscViewerDrawSetBounds(arg1::PetscViewer{Float64},arg2::Integer,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:PetscViewerDrawSetBounds,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscInt,Ptr{Cint}),arg1,arg2,arg3)
    return err
end

function PetscViewerDrawGetBounds(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}))
    err = ccall((:PetscViewerDrawGetBounds,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{PetscInt},Ptr{Ptr{Cint}}),arg1,arg2,arg3)
    return err
end

function PetscViewerSocketSetConnection(arg1::PetscViewer{Float64},arg2::Union(ByteString,Symbol),arg3::Integer)
    err = ccall((:PetscViewerSocketSetConnection,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring,Cint),arg1,arg2,arg3)
    return err
end

function PetscViewerBinarySkipInfo(arg1::PetscViewer{Float64})
    err = ccall((:PetscViewerBinarySkipInfo,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
    return err
end

function PetscViewerBinarySetSkipInfo(arg1::PetscViewer{Float64},arg2::PetscBool)
    err = ccall((:PetscViewerBinarySetSkipInfo,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscBool),arg1,arg2)
    return err
end

function PetscViewerBinaryGetSkipInfo(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscViewerBinaryGetSkipInfo,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function PetscViewerBinarySetSkipOptions(arg1::PetscViewer{Float64},arg2::PetscBool)
    err = ccall((:PetscViewerBinarySetSkipOptions,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscBool),arg1,arg2)
    return err
end

function PetscViewerBinaryGetSkipOptions(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscViewerBinaryGetSkipOptions,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function PetscViewerBinarySetSkipHeader(arg1::PetscViewer{Float64},arg2::PetscBool)
    err = ccall((:PetscViewerBinarySetSkipHeader,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscBool),arg1,arg2)
    return err
end

function PetscViewerBinaryGetSkipHeader(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscViewerBinaryGetSkipHeader,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function PetscViewerBinaryReadStringArray(arg1::PetscViewer{Float64},arg2::Union(Ptr{Ptr{Ptr{Uint8}}},StridedArray{Ptr{Ptr{Uint8}}},Ptr{Void}))
    err = ccall((:PetscViewerBinaryReadStringArray,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Ptr{Ptr{Uint8}}}),arg1,arg2)
    return err
end

function PetscViewerBinaryWriteStringArray(arg1::PetscViewer{Float64},arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:PetscViewerBinaryWriteStringArray,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Ptr{Uint8}}),arg1,arg2)
    return err
end

function PetscViewerFileSetName(arg1::PetscViewer{Float64},arg2::Union(ByteString,Symbol))
    err = ccall((:PetscViewerFileSetName,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring),arg1,arg2)
    return err
end

function PetscViewerFileGetName(arg1::PetscViewer{Float64},arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:PetscViewerFileGetName,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Ptr{Uint8}}),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function PetscViewerVUGetPointer(arg1::PetscViewer{Float64},arg2::Union(Ptr{Ptr{FILE}},StridedArray{Ptr{FILE}},Ptr{Void}))
    ccall((:PetscViewerVUGetPointer,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Ptr{FILE}}),arg1,arg2)
end 
=#
function PetscViewerVUSetVecSeen(arg1::PetscViewer{Float64},arg2::PetscBool)
    err = ccall((:PetscViewerVUSetVecSeen,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscBool),arg1,arg2)
    return err
end

function PetscViewerVUGetVecSeen(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscViewerVUGetVecSeen,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function PetscViewerVUFlushDeferred(arg1::PetscViewer{Float64})
    err = ccall((:PetscViewerVUFlushDeferred,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
    return err
end

function PetscViewerMathematicaInitializePackage(arg0::Type{Float64})
    err = ccall((:PetscViewerMathematicaInitializePackage,petsc1),PetscErrorCode,())
    return err
end

function PetscViewerMathematicaFinalizePackage(arg0::Type{Float64})
    err = ccall((:PetscViewerMathematicaFinalizePackage,petsc1),PetscErrorCode,())
    return err
end

function PetscViewerMathematicaGetName(arg1::PetscViewer{Float64},arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:PetscViewerMathematicaGetName,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Ptr{Uint8}}),arg1,arg2)
    return err
end

function PetscViewerMathematicaSetName(arg1::PetscViewer{Float64},arg2::Union(ByteString,Symbol))
    err = ccall((:PetscViewerMathematicaSetName,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring),arg1,arg2)
    return err
end

function PetscViewerMathematicaClearName(arg1::PetscViewer{Float64})
    err = ccall((:PetscViewerMathematicaClearName,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
    return err
end

function PetscViewerMathematicaSkipPackets(arg1::PetscViewer{Float64},arg2::Integer)
    err = ccall((:PetscViewerMathematicaSkipPackets,petsc1),PetscErrorCode,(PetscViewer{Float64},Cint),arg1,arg2)
    return err
end

function PetscViewerSiloGetName(arg1::PetscViewer{Float64},arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:PetscViewerSiloGetName,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Ptr{Uint8}}),arg1,arg2)
    return err
end

function PetscViewerSiloSetName(arg1::PetscViewer{Float64},arg2::Union(ByteString,Symbol))
    err = ccall((:PetscViewerSiloSetName,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring),arg1,arg2)
    return err
end

function PetscViewerSiloClearName(arg1::PetscViewer{Float64})
    err = ccall((:PetscViewerSiloClearName,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
    return err
end

function PetscViewerSiloGetMeshName(arg1::PetscViewer{Float64},arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:PetscViewerSiloGetMeshName,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Ptr{Uint8}}),arg1,arg2)
    return err
end

function PetscViewerSiloSetMeshName(arg1::PetscViewer{Float64},arg2::Union(ByteString,Symbol))
    err = ccall((:PetscViewerSiloSetMeshName,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring),arg1,arg2)
    return err
end

function PetscViewerSiloClearMeshName(arg1::PetscViewer{Float64})
    err = ccall((:PetscViewerSiloClearMeshName,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
    return err
end

function PetscViewerNetcdfOpen(arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::PetscFileMode,arg4::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    err = ccall((:PetscViewerNetcdfOpen,petsc1),PetscErrorCode,(comm_type,Cstring,PetscFileMode,Ptr{PetscViewer{Float64}}),arg1.val,arg2,arg3,arg4)
    return err
end

function PetscViewerNetcdfGetID(arg1::PetscViewer{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:PetscViewerNetcdfGetID,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Cint}),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function PetscViewerVTKAddField(arg1::PetscViewer{Float64},arg2::PetscObject,PetscViewerVTKWriteFunction::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::PetscViewerVTKFieldType,arg4::PetscObject)
    ccall((:PetscViewerVTKAddField,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscObject,Ptr{Void},PetscViewerVTKFieldType,PetscObject),arg1,arg2,PetscViewerVTKWriteFunction,arg3,arg4)
end 
=#
function PetscViewerVTKOpen(arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::PetscFileMode,arg4::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    err = ccall((:PetscViewerVTKOpen,petsc1),PetscErrorCode,(comm_type,Cstring,PetscFileMode,Ptr{PetscViewer{Float64}}),arg1.val,arg2,arg3,arg4)
    return err
end

function PETSC_VIEWER_STDOUT_(arg0::Type{Float64},arg1::MPI_Comm)
    err = ccall((:PETSC_VIEWER_STDOUT_,petsc1),PetscViewer,(comm_type,),arg1.val)
    return err
end

function PetscViewerASCIIGetStdout(arg1::MPI_Comm,arg2::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    err = ccall((:PetscViewerASCIIGetStdout,petsc1),PetscErrorCode,(comm_type,Ptr{PetscViewer{Float64}}),arg1.val,arg2)
    return err
end

function PETSC_VIEWER_STDERR_(arg0::Type{Float64},arg1::MPI_Comm)
    err = ccall((:PETSC_VIEWER_STDERR_,petsc1),PetscViewer,(comm_type,),arg1.val)
    return err
end

function PetscViewerASCIIGetStderr(arg1::MPI_Comm,arg2::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    err = ccall((:PetscViewerASCIIGetStderr,petsc1),PetscErrorCode,(comm_type,Ptr{PetscViewer{Float64}}),arg1.val,arg2)
    return err
end

function PETSC_VIEWER_DRAW_(arg0::Type{Float64},arg1::MPI_Comm)
    err = ccall((:PETSC_VIEWER_DRAW_,petsc1),PetscViewer,(comm_type,),arg1.val)
    return err
end

function PETSC_VIEWER_SOCKET_(arg0::Type{Float64},arg1::MPI_Comm)
    err = ccall((:PETSC_VIEWER_SOCKET_,petsc1),PetscViewer,(comm_type,),arg1.val)
    return err
end

function PETSC_VIEWER_BINARY_(arg0::Type{Float64},arg1::MPI_Comm)
    err = ccall((:PETSC_VIEWER_BINARY_,petsc1),PetscViewer,(comm_type,),arg1.val)
    return err
end

function PETSC_VIEWER_MATLAB_(arg0::Type{Float64},arg1::MPI_Comm)
    err = ccall((:PETSC_VIEWER_MATLAB_,petsc1),PetscViewer,(comm_type,),arg1.val)
    return err
end

function PETSC_VIEWER_HDF5_(arg0::Type{Float64},arg1::MPI_Comm)
    err = ccall((:PETSC_VIEWER_HDF5_,petsc1),PetscViewer,(comm_type,),arg1.val)
    return err
end

function PetscViewerMatlabGetArray(arg1::PetscViewer{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::Union(ByteString,Symbol))
    err = ccall((:PetscViewerMatlabGetArray,petsc1),PetscErrorCode,(PetscViewer{Float64},Cint,Cint,Ptr{Float64},Cstring),arg1,arg2,arg3,arg4,arg5)
    return err
end

function PetscViewerMatlabPutVariable(arg1::PetscViewer{Float64},arg2::Union(ByteString,Symbol),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscViewerMatlabPutVariable,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring,Ptr{Void}),arg1,arg2,arg3)
    return err
end

#= skipping function with undefined symbols: 
 function PetscViewersCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscViewers},StridedArray{PetscViewers},Ptr{Void}))
    ccall((:PetscViewersCreate,petsc1),PetscErrorCode,(comm_type,Ptr{PetscViewers}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscViewersDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscViewers},StridedArray{PetscViewers},Ptr{Void}))
    ccall((:PetscViewersDestroy,petsc1),PetscErrorCode,(Ptr{PetscViewers},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscViewersGetViewer(arg1::PetscViewers,arg2::Integer,arg3::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    ccall((:PetscViewersGetViewer,petsc1),PetscErrorCode,(PetscViewers,PetscInt,Ptr{PetscViewer{Float64}}),arg1,arg2,arg3)
end 
=#
function VecInitializePackage(arg0::Type{Float64})
    err = ccall((:VecInitializePackage,petsc1),PetscErrorCode,())
    return err
end

function VecFinalizePackage(arg0::Type{Float64})
    err = ccall((:VecFinalizePackage,petsc1),PetscErrorCode,())
    return err
end

function VecCreate(arg1::MPI_Comm,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:VecCreate,petsc1),PetscErrorCode,(comm_type,Ptr{Vec{Float64}}),arg1.val,arg2)
    return err
end

function VecCreateSeq(arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:VecCreateSeq,petsc1),PetscErrorCode,(comm_type,PetscInt,Ptr{Vec{Float64}}),arg1.val,arg2,arg3)
    return err
end

function VecCreateMPI(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:VecCreateMPI,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,Ptr{Vec{Float64}}),arg1.val,arg2,arg3,arg4)
    return err
end

function VecCreateSeqWithArray(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:VecCreateSeqWithArray,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,Ptr{Float64},Ptr{Vec{Float64}}),arg1.val,arg2,arg3,arg4,arg5)
    return err
end

function VecCreateMPIWithArray(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg6::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:VecCreateMPIWithArray,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,Ptr{Float64},Ptr{Vec{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6)
    return err
end

function VecCreateShared(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:VecCreateShared,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,Ptr{Vec{Float64}}),arg1.val,arg2,arg3,arg4)
    return err
end

function VecSetFromOptions(arg1::Vec{Float64})
    err = ccall((:VecSetFromOptions,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
    return err
end

function VecDestroy(arg1::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:VecDestroy,petsc1),PetscErrorCode,(Ptr{Vec{Float64}},),arg1)
    return err
end

function VecZeroEntries(arg1::Vec{Float64})
    err = ccall((:VecZeroEntries,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
    return err
end

function VecSetOptionsPrefix(arg1::Vec{Float64},arg2::Union(ByteString,Symbol))
    err = ccall((:VecSetOptionsPrefix,petsc1),PetscErrorCode,(Vec{Float64},Cstring),arg1,arg2)
    return err
end

function VecAppendOptionsPrefix(arg1::Vec{Float64},arg2::Union(ByteString,Symbol))
    err = ccall((:VecAppendOptionsPrefix,petsc1),PetscErrorCode,(Vec{Float64},Cstring),arg1,arg2)
    return err
end

function VecGetOptionsPrefix(arg1::Vec{Float64},arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:VecGetOptionsPrefix,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Ptr{Uint8}}),arg1,arg2)
    return err
end

function VecSetSizes(arg1::Vec{Float64},arg2::Integer,arg3::Integer)
    err = ccall((:VecSetSizes,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,PetscInt),arg1,arg2,arg3)
    return err
end

function VecDotNorm2(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:VecDotNorm2,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{Float64},Ptr{Cint}),arg1,arg2,arg3,arg4)
    return err
end

function VecDot(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:VecDot,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{Float64}),arg1,arg2,arg3)
    return err
end

function VecDotRealPart(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:VecDotRealPart,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{Cint}),arg1,arg2,arg3)
    return err
end

function VecTDot(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:VecTDot,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{Float64}),arg1,arg2,arg3)
    return err
end

function VecMDot(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:VecMDot,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Ptr{Vec{Float64}},Ptr{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function VecMTDot(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:VecMTDot,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Ptr{Vec{Float64}},Ptr{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function VecGetSubVector(arg1::Vec{Float64},arg2::IS{Float64},arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:VecGetSubVector,petsc1),PetscErrorCode,(Vec{Float64},IS{Float64},Ptr{Vec{Float64}}),arg1,arg2,arg3)
    return err
end

function VecRestoreSubVector(arg1::Vec{Float64},arg2::IS{Float64},arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:VecRestoreSubVector,petsc1),PetscErrorCode,(Vec{Float64},IS{Float64},Ptr{Vec{Float64}}),arg1,arg2,arg3)
    return err
end

function VecNorm(arg1::Vec{Float64},arg2::NormType,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:VecNorm,petsc1),PetscErrorCode,(Vec{Float64},NormType,Ptr{Cint}),arg1,arg2,arg3)
    return err
end

function VecNormAvailable(arg1::Vec{Float64},arg2::NormType,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:VecNormAvailable,petsc1),PetscErrorCode,(Vec{Float64},NormType,Ptr{PetscBool},Ptr{Cint}),arg1,arg2,arg3,arg4)
    return err
end

function VecNormalize(arg1::Vec{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:VecNormalize,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Cint}),arg1,arg2)
    return err
end

function VecSum(arg1::Vec{Float64},arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:VecSum,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Float64}),arg1,arg2)
    return err
end

function VecMax(arg1::Vec{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:VecMax,petsc1),PetscErrorCode,(Vec{Float64},Ptr{PetscInt},Ptr{Cint}),arg1,arg2,arg3)
    return err
end

function VecMin(arg1::Vec{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:VecMin,petsc1),PetscErrorCode,(Vec{Float64},Ptr{PetscInt},Ptr{Cint}),arg1,arg2,arg3)
    return err
end

function VecScale(arg1::Vec{Float64},arg2::Float64)
    err = ccall((:VecScale,petsc1),PetscErrorCode,(Vec{Float64},Float64),arg1,arg2)
    return err
end

function VecCopy(arg1::Vec{Float64},arg2::Vec{Float64})
    err = ccall((:VecCopy,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64}),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function VecSetRandom(arg1::Vec{Float64},arg2::PetscRandom)
    ccall((:VecSetRandom,petsc1),PetscErrorCode,(Vec{Float64},PetscRandom),arg1,arg2)
end 
=#
function VecSet(arg1::Vec{Float64},arg2::Float64)
    err = ccall((:VecSet,petsc1),PetscErrorCode,(Vec{Float64},Float64),arg1,arg2)
    return err
end

function VecSetInf(arg1::Vec{Float64})
    err = ccall((:VecSetInf,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
    return err
end

function VecSwap(arg1::Vec{Float64},arg2::Vec{Float64})
    err = ccall((:VecSwap,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64}),arg1,arg2)
    return err
end

function VecAXPY(arg1::Vec{Float64},arg2::Float64,arg3::Vec{Float64})
    err = ccall((:VecAXPY,petsc1),PetscErrorCode,(Vec{Float64},Float64,Vec{Float64}),arg1,arg2,arg3)
    return err
end

function VecAXPBY(arg1::Vec{Float64},arg2::Float64,arg3::Float64,arg4::Vec{Float64})
    err = ccall((:VecAXPBY,petsc1),PetscErrorCode,(Vec{Float64},Float64,Float64,Vec{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function VecMAXPY(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg4::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:VecMAXPY,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Ptr{Float64},Ptr{Vec{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function VecAYPX(arg1::Vec{Float64},arg2::Float64,arg3::Vec{Float64})
    err = ccall((:VecAYPX,petsc1),PetscErrorCode,(Vec{Float64},Float64,Vec{Float64}),arg1,arg2,arg3)
    return err
end

function VecWAXPY(arg1::Vec{Float64},arg2::Float64,arg3::Vec{Float64},arg4::Vec{Float64})
    err = ccall((:VecWAXPY,petsc1),PetscErrorCode,(Vec{Float64},Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function VecAXPBYPCZ(arg1::Vec{Float64},arg2::Float64,arg3::Float64,arg4::Float64,arg5::Vec{Float64},arg6::Vec{Float64})
    err = ccall((:VecAXPBYPCZ,petsc1),PetscErrorCode,(Vec{Float64},Float64,Float64,Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function VecPointwiseMax(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:VecPointwiseMax,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function VecPointwiseMaxAbs(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:VecPointwiseMaxAbs,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function VecPointwiseMin(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:VecPointwiseMin,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function VecPointwiseMult(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:VecPointwiseMult,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function VecPointwiseDivide(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:VecPointwiseDivide,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function VecMaxPointwiseDivide(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:VecMaxPointwiseDivide,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{Cint}),arg1,arg2,arg3)
    return err
end

function VecShift(arg1::Vec{Float64},arg2::Float64)
    err = ccall((:VecShift,petsc1),PetscErrorCode,(Vec{Float64},Float64),arg1,arg2)
    return err
end

function VecReciprocal(arg1::Vec{Float64})
    err = ccall((:VecReciprocal,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
    return err
end

function VecPermute(arg1::Vec{Float64},arg2::IS{Float64},arg3::PetscBool)
    err = ccall((:VecPermute,petsc1),PetscErrorCode,(Vec{Float64},IS{Float64},PetscBool),arg1,arg2,arg3)
    return err
end

function VecSqrtAbs(arg1::Vec{Float64})
    err = ccall((:VecSqrtAbs,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
    return err
end

function VecLog(arg1::Vec{Float64})
    err = ccall((:VecLog,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
    return err
end

function VecExp(arg1::Vec{Float64})
    err = ccall((:VecExp,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
    return err
end

function VecAbs(arg1::Vec{Float64})
    err = ccall((:VecAbs,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
    return err
end

function VecDuplicate(arg1::Vec{Float64},arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:VecDuplicate,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Vec{Float64}}),arg1,arg2)
    return err
end

function VecDuplicateVecs(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Ptr{Vec{Float64}}},StridedArray{Ptr{Vec{Float64}}},Ptr{Void}))
    err = ccall((:VecDuplicateVecs,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Ptr{Ptr{Vec{Float64}}}),arg1,arg2,arg3)
    return err
end

function VecDestroyVecs(arg1::Integer,arg2::Union(Ptr{Ptr{Vec{Float64}}},StridedArray{Ptr{Vec{Float64}}},Ptr{Void}))
    err = ccall((:VecDestroyVecs,petsc1),PetscErrorCode,(PetscInt,Ptr{Ptr{Vec{Float64}}}),arg1,arg2)
    return err
end

function VecStrideNormAll(arg1::Vec{Float64},arg2::NormType,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:VecStrideNormAll,petsc1),PetscErrorCode,(Vec{Float64},NormType,Ptr{Cint}),arg1,arg2,PetscReal)
    return err
end

function VecStrideMaxAll(arg1::Vec{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:VecStrideMaxAll,petsc1),PetscErrorCode,(Vec{Float64},Ptr{PetscInt},Ptr{Cint}),arg1,arg2,PetscReal)
    return err
end

function VecStrideMinAll(arg1::Vec{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:VecStrideMinAll,petsc1),PetscErrorCode,(Vec{Float64},Ptr{PetscInt},Ptr{Cint}),arg1,arg2,PetscReal)
    return err
end

function VecStrideScaleAll(arg1::Vec{Float64},arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:VecStrideScaleAll,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Float64}),arg1,arg2)
    return err
end

function VecUniqueEntries(arg1::Vec{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    err = ccall((:VecUniqueEntries,petsc1),PetscErrorCode,(Vec{Float64},Ptr{PetscInt},Ptr{Ptr{Float64}}),arg1,arg2,arg3)
    return err
end

function VecStrideNorm(arg1::Vec{Float64},arg2::Integer,arg3::NormType,arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:VecStrideNorm,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,NormType,Ptr{Cint}),arg1,arg2,arg3,arg4)
    return err
end

function VecStrideMax(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:VecStrideMax,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Ptr{PetscInt},Ptr{Cint}),arg1,arg2,arg3,arg4)
    return err
end

function VecStrideMin(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:VecStrideMin,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Ptr{PetscInt},Ptr{Cint}),arg1,arg2,arg3,arg4)
    return err
end

function VecStrideScale(arg1::Vec{Float64},arg2::Integer,arg3::Float64)
    err = ccall((:VecStrideScale,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Float64),arg1,arg2,arg3)
    return err
end

function VecStrideSet(arg1::Vec{Float64},arg2::Integer,arg3::Float64)
    err = ccall((:VecStrideSet,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Float64),arg1,arg2,arg3)
    return err
end

function VecStrideGather(arg1::Vec{Float64},arg2::Integer,arg3::Vec{Float64},arg4::InsertMode)
    err = ccall((:VecStrideGather,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Vec{Float64},InsertMode),arg1,arg2,arg3,arg4)
    return err
end

function VecStrideScatter(arg1::Vec{Float64},arg2::Integer,arg3::Vec{Float64},arg4::InsertMode)
    err = ccall((:VecStrideScatter,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Vec{Float64},InsertMode),arg1,arg2,arg3,arg4)
    return err
end

function VecStrideGatherAll(arg1::Vec{Float64},arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg3::InsertMode)
    err = ccall((:VecStrideGatherAll,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Vec{Float64}},InsertMode),arg1,arg2,arg3)
    return err
end

function VecStrideScatterAll(arg1::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg2::Vec{Float64},arg3::InsertMode)
    err = ccall((:VecStrideScatterAll,petsc1),PetscErrorCode,(Ptr{Vec{Float64}},Vec{Float64},InsertMode),arg1,arg2,arg3)
    return err
end

function VecStrideSubSetScatter(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Vec{Float64},arg6::InsertMode)
    err = ccall((:VecStrideSubSetScatter,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Ptr{PetscInt},Ptr{PetscInt},Vec{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function VecStrideSubSetGather(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Vec{Float64},arg6::InsertMode)
    err = ccall((:VecStrideSubSetGather,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Ptr{PetscInt},Ptr{PetscInt},Vec{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function VecSetValues(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::InsertMode)
    err = ccall((:VecSetValues,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Ptr{PetscInt},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5)
    return err
end

function VecGetValues(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:VecGetValues,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Ptr{PetscInt},Ptr{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function VecAssemblyBegin(arg1::Vec{Float64})
    err = ccall((:VecAssemblyBegin,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
    return err
end

function VecAssemblyEnd(arg1::Vec{Float64})
    err = ccall((:VecAssemblyEnd,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
    return err
end

function VecStashSetInitialSize(arg1::Vec{Float64},arg2::Integer,arg3::Integer)
    err = ccall((:VecStashSetInitialSize,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,PetscInt),arg1,arg2,arg3)
    return err
end

function VecStashView(arg1::Vec{Float64},arg2::PetscViewer{Float64})
    err = ccall((:VecStashView,petsc1),PetscErrorCode,(Vec{Float64},PetscViewer{Float64}),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function VecStashViewFromOptions(arg1::Vec{Float64},arg2::PetscObject,arg3::Union(ByteString,Symbol))
    ccall((:VecStashViewFromOptions,petsc1),PetscErrorCode,(Vec{Float64},PetscObject,Cstring),arg1,arg2,arg3)
end 
=#
function VecStashGetInfo(arg1::Vec{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:VecStashGetInfo,petsc1),PetscErrorCode,(Vec{Float64},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function VecGetBlockSize(arg1::Vec{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:VecGetBlockSize,petsc1),PetscErrorCode,(Vec{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function VecSetValuesBlocked(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::InsertMode)
    err = ccall((:VecSetValuesBlocked,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Ptr{PetscInt},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5)
    return err
end

function VecSetType(arg1::Vec{Float64},arg2::VecType)
    err = ccall((:VecSetType,petsc1),PetscErrorCode,(Vec{Float64},Cstring),arg1,arg2)
    return err
end

function VecGetType(arg1::Vec{Float64},arg2::Union(Ptr{VecType},StridedArray{VecType},Ptr{Void}))
    (arg2_,tmp) = symbol_get_before(arg2)
    err = ccall((:VecGetType,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Ptr{Uint8}}),arg1,arg2_)
    symbol_get_after(arg2_,arg2)
    return err
end

function VecRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:VecRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

function VecScatterCreate(arg1::Vec{Float64},arg2::IS{Float64},arg3::Vec{Float64},arg4::IS{Float64},arg5::Union(Ptr{VecScatter{Float64}},StridedArray{VecScatter{Float64}},Ptr{Void}))
    err = ccall((:VecScatterCreate,petsc1),PetscErrorCode,(Vec{Float64},IS{Float64},Vec{Float64},IS{Float64},Ptr{VecScatter{Float64}}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function VecScatterCreateEmpty(arg1::MPI_Comm,arg2::Union(Ptr{VecScatter{Float64}},StridedArray{VecScatter{Float64}},Ptr{Void}))
    err = ccall((:VecScatterCreateEmpty,petsc1),PetscErrorCode,(comm_type,Ptr{VecScatter{Float64}}),arg1.val,arg2)
    return err
end

function VecScatterCreateLocal(arg1::VecScatter{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Integer,arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg8::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg9::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg10::Integer)
    err = ccall((:VecScatterCreateLocal,petsc1),PetscErrorCode,(VecScatter{Float64},PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},PetscInt),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
    return err
end

function VecScatterBegin(arg1::VecScatter{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::InsertMode,arg5::ScatterMode)
    err = ccall((:VecScatterBegin,petsc1),PetscErrorCode,(VecScatter{Float64},Vec{Float64},Vec{Float64},InsertMode,ScatterMode),arg1,arg2,arg3,arg4,arg5)
    return err
end

function VecScatterEnd(arg1::VecScatter{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::InsertMode,arg5::ScatterMode)
    err = ccall((:VecScatterEnd,petsc1),PetscErrorCode,(VecScatter{Float64},Vec{Float64},Vec{Float64},InsertMode,ScatterMode),arg1,arg2,arg3,arg4,arg5)
    return err
end

function VecScatterDestroy(arg1::Union(Ptr{VecScatter{Float64}},StridedArray{VecScatter{Float64}},Ptr{Void}))
    err = ccall((:VecScatterDestroy,petsc1),PetscErrorCode,(Ptr{VecScatter{Float64}},),arg1)
    return err
end

function VecScatterCopy(arg1::VecScatter{Float64},arg2::Union(Ptr{VecScatter{Float64}},StridedArray{VecScatter{Float64}},Ptr{Void}))
    err = ccall((:VecScatterCopy,petsc1),PetscErrorCode,(VecScatter{Float64},Ptr{VecScatter{Float64}}),arg1,arg2)
    return err
end

function VecScatterView(arg1::VecScatter{Float64},arg2::PetscViewer{Float64})
    err = ccall((:VecScatterView,petsc1),PetscErrorCode,(VecScatter{Float64},PetscViewer{Float64}),arg1,arg2)
    return err
end

function VecScatterGetMerged(arg1::VecScatter{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:VecScatterGetMerged,petsc1),PetscErrorCode,(VecScatter{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function VecGetArray4d(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Integer,arg9::Integer,arg10::Union(Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}},StridedArray{Ptr{Ptr{Ptr{Ptr{Float64}}}}},Ptr{Void}))
    err = ccall((:VecGetArray4d,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
    return err
end

function VecRestoreArray4d(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Integer,arg9::Integer,arg10::Union(Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}},StridedArray{Ptr{Ptr{Ptr{Ptr{Float64}}}}},Ptr{Void}))
    err = ccall((:VecRestoreArray4d,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
    return err
end

function VecGetArray3d(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Union(Ptr{Ptr{Ptr{Ptr{Float64}}}},StridedArray{Ptr{Ptr{Ptr{Float64}}}},Ptr{Void}))
    err = ccall((:VecGetArray3d,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{Float64}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
    return err
end

function VecRestoreArray3d(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Union(Ptr{Ptr{Ptr{Ptr{Float64}}}},StridedArray{Ptr{Ptr{Ptr{Float64}}}},Ptr{Void}))
    err = ccall((:VecRestoreArray3d,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{Float64}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
    return err
end

function VecGetArray2d(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Ptr{Ptr{Float64}}},StridedArray{Ptr{Ptr{Float64}}},Ptr{Void}))
    err = ccall((:VecGetArray2d,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Float64}}}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function VecRestoreArray2d(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Ptr{Ptr{Float64}}},StridedArray{Ptr{Ptr{Float64}}},Ptr{Void}))
    err = ccall((:VecRestoreArray2d,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Float64}}}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function VecGetArray1d(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    err = ccall((:VecGetArray1d,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,PetscInt,Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function VecRestoreArray1d(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    err = ccall((:VecRestoreArray1d,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,PetscInt,Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function VecGetArray4dRead(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Integer,arg9::Integer,arg10::Union(Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}},StridedArray{Ptr{Ptr{Ptr{Ptr{Float64}}}}},Ptr{Void}))
    err = ccall((:VecGetArray4dRead,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
    return err
end

function VecRestoreArray4dRead(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Integer,arg9::Integer,arg10::Union(Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}},StridedArray{Ptr{Ptr{Ptr{Ptr{Float64}}}}},Ptr{Void}))
    err = ccall((:VecRestoreArray4dRead,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
    return err
end

function VecGetArray3dRead(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Union(Ptr{Ptr{Ptr{Ptr{Float64}}}},StridedArray{Ptr{Ptr{Ptr{Float64}}}},Ptr{Void}))
    err = ccall((:VecGetArray3dRead,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{Float64}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
    return err
end

function VecRestoreArray3dRead(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Union(Ptr{Ptr{Ptr{Ptr{Float64}}}},StridedArray{Ptr{Ptr{Ptr{Float64}}}},Ptr{Void}))
    err = ccall((:VecRestoreArray3dRead,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{Float64}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
    return err
end

function VecGetArray2dRead(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Ptr{Ptr{Float64}}},StridedArray{Ptr{Ptr{Float64}}},Ptr{Void}))
    err = ccall((:VecGetArray2dRead,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Float64}}}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function VecRestoreArray2dRead(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Ptr{Ptr{Float64}}},StridedArray{Ptr{Ptr{Float64}}},Ptr{Void}))
    err = ccall((:VecRestoreArray2dRead,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Float64}}}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function VecGetArray1dRead(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    err = ccall((:VecGetArray1dRead,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,PetscInt,Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function VecRestoreArray1dRead(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    err = ccall((:VecRestoreArray1dRead,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,PetscInt,Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function VecPlaceArray(arg1::Vec{Float64},arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:VecPlaceArray,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Float64}),arg1,arg2)
    return err
end

function VecResetArray(arg1::Vec{Float64})
    err = ccall((:VecResetArray,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
    return err
end

function VecReplaceArray(arg1::Vec{Float64},arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:VecReplaceArray,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Float64}),arg1,arg2)
    return err
end

function VecGetArrays(arg1::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg2::Integer,arg3::Union(Ptr{Ptr{Ptr{Float64}}},StridedArray{Ptr{Ptr{Float64}}},Ptr{Void}))
    err = ccall((:VecGetArrays,petsc1),PetscErrorCode,(Ptr{Vec{Float64}},PetscInt,Ptr{Ptr{Ptr{Float64}}}),arg1,arg2,arg3)
    return err
end

function VecRestoreArrays(arg1::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg2::Integer,arg3::Union(Ptr{Ptr{Ptr{Float64}}},StridedArray{Ptr{Ptr{Float64}}},Ptr{Void}))
    err = ccall((:VecRestoreArrays,petsc1),PetscErrorCode,(Ptr{Vec{Float64}},PetscInt,Ptr{Ptr{Ptr{Float64}}}),arg1,arg2,arg3)
    return err
end

function VecView(arg1::Vec{Float64},arg2::PetscViewer{Float64})
    err = ccall((:VecView,petsc1),PetscErrorCode,(Vec{Float64},PetscViewer{Float64}),arg1,arg2)
    return err
end

function VecEqual(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:VecEqual,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{PetscBool}),arg1,arg2,arg3)
    return err
end

function VecLoad(arg1::Vec{Float64},arg2::PetscViewer{Float64})
    err = ccall((:VecLoad,petsc1),PetscErrorCode,(Vec{Float64},PetscViewer{Float64}),arg1,arg2)
    return err
end

function VecGetSize(arg1::Vec{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:VecGetSize,petsc1),PetscErrorCode,(Vec{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function VecGetLocalSize(arg1::Vec{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:VecGetLocalSize,petsc1),PetscErrorCode,(Vec{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function VecGetOwnershipRange(arg1::Vec{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:VecGetOwnershipRange,petsc1),PetscErrorCode,(Vec{Float64},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function VecGetOwnershipRanges(arg1::Vec{Float64},arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:VecGetOwnershipRanges,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Ptr{PetscInt}}),arg1,arg2)
    return err
end

function VecSetLocalToGlobalMapping(arg1::Vec{Float64},arg2::ISLocalToGlobalMapping{Float64})
    err = ccall((:VecSetLocalToGlobalMapping,petsc1),PetscErrorCode,(Vec{Float64},ISLocalToGlobalMapping{Float64}),arg1,arg2)
    return err
end

function VecSetValuesLocal(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::InsertMode)
    err = ccall((:VecSetValuesLocal,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Ptr{PetscInt},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5)
    return err
end

function VecGetLocalToGlobalMapping(arg1::Vec{Float64},arg2::Union(Ptr{ISLocalToGlobalMapping{Float64}},StridedArray{ISLocalToGlobalMapping{Float64}},Ptr{Void}))
    err = ccall((:VecGetLocalToGlobalMapping,petsc1),PetscErrorCode,(Vec{Float64},Ptr{ISLocalToGlobalMapping{Float64}}),arg1,arg2)
    return err
end

function VecDotBegin(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:VecDotBegin,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{Float64}),arg1,arg2,arg3)
    return err
end

function VecDotEnd(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:VecDotEnd,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{Float64}),arg1,arg2,arg3)
    return err
end

function VecTDotBegin(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:VecTDotBegin,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{Float64}),arg1,arg2,arg3)
    return err
end

function VecTDotEnd(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:VecTDotEnd,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{Float64}),arg1,arg2,arg3)
    return err
end

function VecNormBegin(arg1::Vec{Float64},arg2::NormType,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:VecNormBegin,petsc1),PetscErrorCode,(Vec{Float64},NormType,Ptr{Cint}),arg1,arg2,arg3)
    return err
end

function VecNormEnd(arg1::Vec{Float64},arg2::NormType,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:VecNormEnd,petsc1),PetscErrorCode,(Vec{Float64},NormType,Ptr{Cint}),arg1,arg2,arg3)
    return err
end

function VecMDotBegin(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:VecMDotBegin,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Ptr{Vec{Float64}},Ptr{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function VecMDotEnd(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:VecMDotEnd,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Ptr{Vec{Float64}},Ptr{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function VecMTDotBegin(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:VecMTDotBegin,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Ptr{Vec{Float64}},Ptr{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function VecMTDotEnd(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:VecMTDotEnd,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Ptr{Vec{Float64}},Ptr{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function PetscCommSplitReductionBegin(arg0::Type{Float64},arg1::MPI_Comm)
    err = ccall((:PetscCommSplitReductionBegin,petsc1),PetscErrorCode,(comm_type,),arg1.val)
    return err
end

function VecSetOption(arg1::Vec{Float64},arg2::VecOption,arg3::PetscBool)
    err = ccall((:VecSetOption,petsc1),PetscErrorCode,(Vec{Float64},VecOption,PetscBool),arg1,arg2,arg3)
    return err
end

function VecGetArray(arg1::Vec{Float64},arg2::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    err = ccall((:VecGetArray,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Ptr{Float64}}),arg1,arg2)
    return err
end

function VecGetArrayRead(arg1::Vec{Float64},arg2::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    err = ccall((:VecGetArrayRead,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Ptr{Float64}}),arg1,arg2)
    return err
end

function VecRestoreArray(arg1::Vec{Float64},arg2::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    err = ccall((:VecRestoreArray,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Ptr{Float64}}),arg1,arg2)
    return err
end

function VecRestoreArrayRead(arg1::Vec{Float64},arg2::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    err = ccall((:VecRestoreArrayRead,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Ptr{Float64}}),arg1,arg2)
    return err
end

function VecGetLocalVector(arg1::Vec{Float64},arg2::Vec{Float64})
    err = ccall((:VecGetLocalVector,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64}),arg1,arg2)
    return err
end

function VecRestoreLocalVector(arg1::Vec{Float64},arg2::Vec{Float64})
    err = ccall((:VecRestoreLocalVector,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64}),arg1,arg2)
    return err
end

function VecGetLocalVectorRead(arg1::Vec{Float64},arg2::Vec{Float64})
    err = ccall((:VecGetLocalVectorRead,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64}),arg1,arg2)
    return err
end

function VecRestoreLocalVectorRead(arg1::Vec{Float64},arg2::Vec{Float64})
    err = ccall((:VecRestoreLocalVectorRead,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64}),arg1,arg2)
    return err
end

function VecContourScale(arg1::Vec{Float64},PetscReal::Integer,arg2::Integer)
    err = ccall((:VecContourScale,petsc1),PetscErrorCode,(Vec{Float64},Cint,Cint),arg1,PetscReal,arg2)
    return err
end

function VecSetOperation(arg1::Vec{Float64},arg2::VecOperation,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:VecSetOperation,petsc1),PetscErrorCode,(Vec{Float64},VecOperation,Ptr{Void}),arg1,arg2,arg3)
    return err
end

function VecMPISetGhost(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:VecMPISetGhost,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function VecCreateGhost(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:VecCreateGhost,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Vec{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6)
    return err
end

function VecCreateGhostWithArray(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:VecCreateGhostWithArray,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Float64},Ptr{Vec{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
    return err
end

function VecCreateGhostBlock(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:VecCreateGhostBlock,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Vec{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
    return err
end

function VecCreateGhostBlockWithArray(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg8::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:VecCreateGhostBlockWithArray,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Float64},Ptr{Vec{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
    return err
end

function VecGhostGetLocalForm(arg1::Vec{Float64},arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:VecGhostGetLocalForm,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Vec{Float64}}),arg1,arg2)
    return err
end

function VecGhostRestoreLocalForm(arg1::Vec{Float64},arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:VecGhostRestoreLocalForm,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Vec{Float64}}),arg1,arg2)
    return err
end

function VecGhostIsLocalForm(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:VecGhostIsLocalForm,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{PetscBool}),arg1,arg2,arg3)
    return err
end

function VecGhostUpdateBegin(arg1::Vec{Float64},arg2::InsertMode,arg3::ScatterMode)
    err = ccall((:VecGhostUpdateBegin,petsc1),PetscErrorCode,(Vec{Float64},InsertMode,ScatterMode),arg1,arg2,arg3)
    return err
end

function VecGhostUpdateEnd(arg1::Vec{Float64},arg2::InsertMode,arg3::ScatterMode)
    err = ccall((:VecGhostUpdateEnd,petsc1),PetscErrorCode,(Vec{Float64},InsertMode,ScatterMode),arg1,arg2,arg3)
    return err
end

function VecConjugate(arg1::Vec{Float64})
    err = ccall((:VecConjugate,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
    return err
end

function VecScatterCreateToAll(arg1::Vec{Float64},arg2::Union(Ptr{VecScatter{Float64}},StridedArray{VecScatter{Float64}},Ptr{Void}),arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:VecScatterCreateToAll,petsc1),PetscErrorCode,(Vec{Float64},Ptr{VecScatter{Float64}},Ptr{Vec{Float64}}),arg1,arg2,arg3)
    return err
end

function VecScatterCreateToZero(arg1::Vec{Float64},arg2::Union(Ptr{VecScatter{Float64}},StridedArray{VecScatter{Float64}},Ptr{Void}),arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:VecScatterCreateToZero,petsc1),PetscErrorCode,(Vec{Float64},Ptr{VecScatter{Float64}},Ptr{Vec{Float64}}),arg1,arg2,arg3)
    return err
end

function ISComplementVec(arg1::IS{Float64},arg2::Vec{Float64},arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:ISComplementVec,petsc1),PetscErrorCode,(IS{Float64},Vec{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3)
    return err
end

function VecPow(arg1::Vec{Float64},arg2::Float64)
    err = ccall((:VecPow,petsc1),PetscErrorCode,(Vec{Float64},Float64),arg1,arg2)
    return err
end

function VecMedian(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64})
    err = ccall((:VecMedian,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function VecWhichBetween(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:VecWhichBetween,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function VecWhichBetweenOrEqual(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:VecWhichBetweenOrEqual,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function VecWhichGreaterThan(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:VecWhichGreaterThan,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3)
    return err
end

function VecWhichLessThan(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:VecWhichLessThan,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3)
    return err
end

function VecWhichEqual(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:VecWhichEqual,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3)
    return err
end

function VecISAXPY(arg1::Vec{Float64},arg2::IS{Float64},arg3::Float64,arg4::Vec{Float64})
    err = ccall((:VecISAXPY,petsc1),PetscErrorCode,(Vec{Float64},IS{Float64},Float64,Vec{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function VecISSet(arg1::Vec{Float64},arg2::IS{Float64},arg3::Float64)
    err = ccall((:VecISSet,petsc1),PetscErrorCode,(Vec{Float64},IS{Float64},Float64),arg1,arg2,arg3)
    return err
end

function VecBoundGradientProjection(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64},arg5::Vec{Float64})
    err = ccall((:VecBoundGradientProjection,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function VecStepBoundInfo(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64},arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg7::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:VecStepBoundInfo,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
    return err
end

function VecStepMax(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:VecStepMax,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{Cint}),arg1,arg2,arg3)
    return err
end

function VecStepMaxBounded(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64},arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:VecStepMaxBounded,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function PetscViewerMathematicaGetVector(arg1::PetscViewer{Float64},arg2::Vec{Float64})
    err = ccall((:PetscViewerMathematicaGetVector,petsc1),PetscErrorCode,(PetscViewer{Float64},Vec{Float64}),arg1,arg2)
    return err
end

function PetscViewerMathematicaPutVector(arg1::PetscViewer{Float64},arg2::Vec{Float64})
    err = ccall((:PetscViewerMathematicaPutVector,petsc1),PetscErrorCode,(PetscViewer{Float64},Vec{Float64}),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function VecsDestroy(arg0::Type{Float64},arg1::Vecs)
    ccall((:VecsDestroy,petsc1),PetscErrorCode,(Vecs,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function VecsCreateSeq(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Vecs},StridedArray{Vecs},Ptr{Void}))
    ccall((:VecsCreateSeq,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,Ptr{Vecs}),arg1.val,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function VecsCreateSeqWithArray(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::Union(Ptr{Vecs},StridedArray{Vecs},Ptr{Void}))
    ccall((:VecsCreateSeqWithArray,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,Ptr{Float64},Ptr{Vecs}),arg1.val,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function VecsDuplicate(arg0::Type{Float64},arg1::Vecs,arg2::Union(Ptr{Vecs},StridedArray{Vecs},Ptr{Void}))
    ccall((:VecsDuplicate,petsc1),PetscErrorCode,(Vecs,Ptr{Vecs}),arg1,arg2)
end 
=#
function VecNestGetSubVecs(arg1::Vec{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{Vec{Float64}}},StridedArray{Ptr{Vec{Float64}}},Ptr{Void}))
    err = ccall((:VecNestGetSubVecs,petsc1),PetscErrorCode,(Vec{Float64},Ptr{PetscInt},Ptr{Ptr{Vec{Float64}}}),arg1,arg2,arg3)
    return err
end

function VecNestGetSubVec(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:VecNestGetSubVec,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Ptr{Vec{Float64}}),arg1,arg2,arg3)
    return err
end

function VecNestSetSubVecs(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:VecNestSetSubVecs,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Ptr{PetscInt},Ptr{Vec{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function VecNestSetSubVec(arg1::Vec{Float64},arg2::Integer,arg3::Vec{Float64})
    err = ccall((:VecNestSetSubVec,petsc1),PetscErrorCode,(Vec{Float64},PetscInt,Vec{Float64}),arg1,arg2,arg3)
    return err
end

function VecCreateNest(arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg5::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:VecCreateNest,petsc1),PetscErrorCode,(comm_type,PetscInt,Ptr{IS{Float64}},Ptr{Vec{Float64}},Ptr{Vec{Float64}}),arg1.val,arg2,arg3,arg4,arg5)
    return err
end

function VecNestGetSize(arg1::Vec{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:VecNestGetSize,petsc1),PetscErrorCode,(Vec{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function PetscOptionsGetVec(arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Vec{Float64},arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PetscOptionsGetVec,petsc1),PetscErrorCode,(Cstring,Cstring,Vec{Float64},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
    return err
end

function VecChop(arg1::Vec{Float64},PetscReal::Integer)
    err = ccall((:VecChop,petsc1),PetscErrorCode,(Vec{Float64},Cint),arg1,PetscReal)
    return err
end

function VecGetLayout(arg1::Vec{Float64},arg2::Union(Ptr{PetscLayout{Float64}},StridedArray{PetscLayout{Float64}},Ptr{Void}))
    err = ccall((:VecGetLayout,petsc1),PetscErrorCode,(Vec{Float64},Ptr{PetscLayout{Float64}{Float64}}),arg1,arg2)
    return err
end

function VecSetLayout(arg1::Vec{Float64},arg2::PetscLayout{Float64})
    err = ccall((:VecSetLayout,petsc1),PetscErrorCode,(Vec{Float64},PetscLayout{Float64}),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function PetscSectionVecView(arg1::PetscSection,arg2::Vec{Float64},arg3::PetscViewer{Float64})
    ccall((:PetscSectionVecView,petsc1),PetscErrorCode,(PetscSection,Vec{Float64},PetscViewer{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function VecGetValuesSection(arg1::Vec{Float64},arg2::PetscSection,arg3::Integer,arg4::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:VecGetValuesSection,petsc1),PetscErrorCode,(Vec{Float64},PetscSection,PetscInt,Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function VecSetValuesSection(arg1::Vec{Float64},arg2::PetscSection,arg3::Integer,arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::InsertMode)
    ccall((:VecSetValuesSection,petsc1),PetscErrorCode,(Vec{Float64},PetscSection,PetscInt,Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionVecNorm(arg1::PetscSection,arg2::PetscSection,arg3::Vec{Float64},arg4::NormType,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscSectionVecNorm,petsc1),PetscErrorCode,(PetscSection,PetscSection,Vec{Float64},NormType,Ptr{Cint}),arg1,arg2,arg3,arg4,PetscReal)
end 
=#
function MatGetFactor(arg1::Mat{Float64},arg2::Union(ByteString,Symbol),arg3::MatFactorType,arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatGetFactor,petsc1),PetscErrorCode,(Mat{Float64},Cstring,MatFactorType,Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function MatGetFactorAvailable(arg1::Mat{Float64},arg2::Union(ByteString,Symbol),arg3::MatFactorType,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:MatGetFactorAvailable,petsc1),PetscErrorCode,(Mat{Float64},Cstring,MatFactorType,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
    return err
end

function MatFactorGetSolverPackage(arg1::Mat{Float64},arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:MatFactorGetSolverPackage,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Ptr{Uint8}}),arg1,arg2)
    return err
end

function MatGetFactorType(arg1::Mat{Float64},arg2::Union(Ptr{MatFactorType},StridedArray{MatFactorType},Ptr{Void}))
    err = ccall((:MatGetFactorType,petsc1),PetscErrorCode,(Mat{Float64},Ptr{MatFactorType}),arg1,arg2)
    return err
end

function MatSolverPackageRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::MatType,arg3::MatFactorType,arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:MatSolverPackageRegister,petsc1),PetscErrorCode,(Cstring,Cstring,MatFactorType,Ptr{Void}),arg1,arg2,arg3,arg4)
    return err
end

function MatSolverPackageGet(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::MatType,arg3::MatFactorType,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg6::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    err = ccall((:MatSolverPackageGet,petsc1),PetscErrorCode,(Cstring,Cstring,MatFactorType,Ptr{PetscBool},Ptr{PetscBool},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatInitializePackage(arg0::Type{Float64})
    err = ccall((:MatInitializePackage,petsc1),PetscErrorCode,())
    return err
end

function MatCreate(arg1::MPI_Comm,arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreate,petsc1),PetscErrorCode,(comm_type,Ptr{Mat{Float64}}),arg1.val,arg2)
    return err
end

function MatSetSizes(arg1::Mat{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer)
    err = ccall((:MatSetSizes,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4,arg5)
    return err
end

function MatSetType(arg1::Mat{Float64},arg2::MatType)
    err = ccall((:MatSetType,petsc1),PetscErrorCode,(Mat{Float64},Cstring),arg1,arg2)
    return err
end

function MatSetFromOptions(arg1::Mat{Float64})
    err = ccall((:MatSetFromOptions,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
    return err
end

function MatRegisterBaseName(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol))
    err = ccall((:MatRegisterBaseName,petsc1),PetscErrorCode,(Cstring,Cstring,Cstring),arg1,arg2,arg3)
    return err
end

function MatSetOptionsPrefix(arg1::Mat{Float64},arg2::Union(ByteString,Symbol))
    err = ccall((:MatSetOptionsPrefix,petsc1),PetscErrorCode,(Mat{Float64},Cstring),arg1,arg2)
    return err
end

function MatAppendOptionsPrefix(arg1::Mat{Float64},arg2::Union(ByteString,Symbol))
    err = ccall((:MatAppendOptionsPrefix,petsc1),PetscErrorCode,(Mat{Float64},Cstring),arg1,arg2)
    return err
end

function MatGetOptionsPrefix(arg1::Mat{Float64},arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:MatGetOptionsPrefix,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Ptr{Uint8}}),arg1,arg2)
    return err
end

function MatSetErrorIfFPE(arg1::Mat{Float64},arg2::PetscBool)
    err = ccall((:MatSetErrorIfFPE,petsc1),PetscErrorCode,(Mat{Float64},PetscBool),arg1,arg2)
    return err
end

function MatCreateSeqDense(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateSeqDense,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,Ptr{Float64},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5)
    return err
end

function MatCreateDense(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateDense,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Float64},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
    return err
end

function MatCreateSeqAIJ(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateSeqAIJ,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatCreateAIJ(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg8::Integer,arg9::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg10::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateAIJ,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
    return err
end

function MatCreateMPIAIJWithArrays(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg8::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg9::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateMPIAIJWithArrays,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
    return err
end

function MatCreateMPIAIJWithSplitArrays(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg8::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg9::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg10::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg11::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg12::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateMPIAIJWithSplitArrays,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64},Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12)
    return err
end

function MatCreateSeqBAIJ(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateSeqBAIJ,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
    return err
end

function MatCreateBAIJ(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg9::Integer,arg10::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg11::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateBAIJ,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
    return err
end

function MatCreateMPIBAIJWithArrays(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg8::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg9::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg10::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateMPIBAIJWithArrays,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
    return err
end

function MatCreateMPIAdj(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateMPIAdj,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
    return err
end

function MatCreateSeqSBAIJ(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateSeqSBAIJ,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
    return err
end

function MatCreateSBAIJ(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg9::Integer,arg10::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg11::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateSBAIJ,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
    return err
end

function MatCreateMPISBAIJWithArrays(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg8::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg9::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg10::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateMPISBAIJWithArrays,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
    return err
end

function MatSeqSBAIJSetPreallocationCSR(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:MatSeqSBAIJSetPreallocationCSR,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function MatMPISBAIJSetPreallocationCSR(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:MatMPISBAIJSetPreallocationCSR,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function MatXAIJSetPreallocation(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatXAIJSetPreallocation,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatCreateShell(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateShell,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Void},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
    return err
end

function MatCreateNormal(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateNormal,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2)
    return err
end

function MatCreateLRC(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64},arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateLRC,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function MatCreateIS(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::ISLocalToGlobalMapping{Float64},arg8::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateIS,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,ISLocalToGlobalMapping{Float64},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
    return err
end

function MatCreateSeqAIJCRL(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateSeqAIJCRL,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatCreateMPIAIJCRL(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Integer,arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg8::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateMPIAIJCRL,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
    return err
end

function MatCreateSeqBSTRM(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateSeqBSTRM,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
    return err
end

function MatCreateMPIBSTRM(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg9::Integer,arg10::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg11::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateMPIBSTRM,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
    return err
end

function MatCreateSeqSBSTRM(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateSeqSBSTRM,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
    return err
end

function MatCreateMPISBSTRM(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg9::Integer,arg10::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg11::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateMPISBSTRM,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
    return err
end

function MatCreateScatter(arg1::MPI_Comm,arg2::VecScatter{Float64},arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateScatter,petsc1),PetscErrorCode,(comm_type,VecScatter{Float64},Ptr{Mat{Float64}}),arg1.val,arg2,arg3)
    return err
end

function MatScatterSetVecScatter(arg1::Mat{Float64},arg2::VecScatter{Float64})
    err = ccall((:MatScatterSetVecScatter,petsc1),PetscErrorCode,(Mat{Float64},VecScatter{Float64}),arg1,arg2)
    return err
end

function MatScatterGetVecScatter(arg1::Mat{Float64},arg2::Union(Ptr{VecScatter{Float64}},StridedArray{VecScatter{Float64}},Ptr{Void}))
    err = ccall((:MatScatterGetVecScatter,petsc1),PetscErrorCode,(Mat{Float64},Ptr{VecScatter{Float64}}),arg1,arg2)
    return err
end

function MatCreateBlockMat(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateBlockMat,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
    return err
end

function MatCompositeAddMat(arg1::Mat{Float64},arg2::Mat{Float64})
    err = ccall((:MatCompositeAddMat,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64}),arg1,arg2)
    return err
end

function MatCompositeMerge(arg1::Mat{Float64})
    err = ccall((:MatCompositeMerge,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
    return err
end

function MatCreateComposite(arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateComposite,petsc1),PetscErrorCode,(comm_type,PetscInt,Ptr{Mat{Float64}},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4)
    return err
end

function MatCompositeSetType(arg1::Mat{Float64},arg2::MatCompositeType)
    err = ccall((:MatCompositeSetType,petsc1),PetscErrorCode,(Mat{Float64},MatCompositeType),arg1,arg2)
    return err
end

function MatCreateFFT(arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::MatType,arg5::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateFFT,petsc1),PetscErrorCode,(comm_type,PetscInt,Ptr{PetscInt},Cstring,Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5)
    return err
end

function MatCreateSeqCUFFT(arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateSeqCUFFT,petsc1),PetscErrorCode,(comm_type,PetscInt,Ptr{PetscInt},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4)
    return err
end

function MatCreateTranspose(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateTranspose,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2)
    return err
end

function MatCreateHermitianTranspose(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateHermitianTranspose,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2)
    return err
end

function MatCreateSubMatrix(arg1::Mat{Float64},arg2::IS{Float64},arg3::IS{Float64},arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateSubMatrix,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},IS{Float64},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function MatSubMatrixUpdate(arg1::Mat{Float64},arg2::Mat{Float64},arg3::IS{Float64},arg4::IS{Float64})
    err = ccall((:MatSubMatrixUpdate,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},IS{Float64},IS{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function MatCreateLocalRef(arg1::Mat{Float64},arg2::IS{Float64},arg3::IS{Float64},arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateLocalRef,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},IS{Float64},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function MatPythonSetType(arg1::Mat{Float64},arg2::Union(ByteString,Symbol))
    err = ccall((:MatPythonSetType,petsc1),PetscErrorCode,(Mat{Float64},Cstring),arg1,arg2)
    return err
end

function MatSetUp(arg1::Mat{Float64})
    err = ccall((:MatSetUp,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
    return err
end

function MatDestroy(arg1::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatDestroy,petsc1),PetscErrorCode,(Ptr{Mat{Float64}},),arg1)
    return err
end

function MatGetNonzeroState(arg1::Mat{Float64},arg2::Union(Ptr{PetscObjectState},StridedArray{PetscObjectState},Ptr{Void}))
    err = ccall((:MatGetNonzeroState,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscObjectState}),arg1,arg2)
    return err
end

function MatConjugate(arg1::Mat{Float64})
    err = ccall((:MatConjugate,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
    return err
end

function MatRealPart(arg1::Mat{Float64})
    err = ccall((:MatRealPart,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
    return err
end

function MatImaginaryPart(arg1::Mat{Float64})
    err = ccall((:MatImaginaryPart,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
    return err
end

function MatGetDiagonalBlock(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatGetDiagonalBlock,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2)
    return err
end

function MatGetTrace(arg1::Mat{Float64},arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:MatGetTrace,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Float64}),arg1,arg2)
    return err
end

function MatInvertBlockDiagonal(arg1::Mat{Float64},arg2::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    err = ccall((:MatInvertBlockDiagonal,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Ptr{Float64}}),arg1,arg2)
    return err
end

function MatSetValues(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::InsertMode)
    err = ccall((:MatSetValues,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
    return err
end

function MatSetValuesBlocked(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::InsertMode)
    err = ccall((:MatSetValuesBlocked,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
    return err
end

function MatSetValuesRow(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:MatSetValuesRow,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{Float64}),arg1,arg2,arg3)
    return err
end

function MatSetValuesRowLocal(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:MatSetValuesRowLocal,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{Float64}),arg1,arg2,arg3)
    return err
end

function MatSetValuesBatch(arg1::Mat{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:MatSetValuesBatch,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,PetscInt,Ptr{PetscInt},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5)
    return err
end

#= skipping function with undefined symbols: 
 function MatSetRandom(arg1::Mat{Float64},arg2::PetscRandom)
    ccall((:MatSetRandom,petsc1),PetscErrorCode,(Mat{Float64},PetscRandom),arg1,arg2)
end 
=#
function MatSetValuesStencil(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{MatStencil},StridedArray{MatStencil},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{MatStencil},StridedArray{MatStencil},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::InsertMode)
    err = ccall((:MatSetValuesStencil,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{MatStencil},PetscInt,Ptr{MatStencil},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
    return err
end

function MatSetValuesBlockedStencil(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{MatStencil},StridedArray{MatStencil},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{MatStencil},StridedArray{MatStencil},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::InsertMode)
    err = ccall((:MatSetValuesBlockedStencil,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{MatStencil},PetscInt,Ptr{MatStencil},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
    return err
end

function MatSetStencil(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Integer)
    err = ccall((:MatSetStencil,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},Ptr{PetscInt},PetscInt),arg1,arg2,arg3,arg4,arg5)
    return err
end

function MatSetColoring(arg1::Mat{Float64},arg2::ISColoring{Float64})
    err = ccall((:MatSetColoring,petsc1),PetscErrorCode,(Mat{Float64},ISColoring{Float64}),arg1,arg2)
    return err
end

function MatSetValuesAdifor(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:MatSetValuesAdifor,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{Void}),arg1,arg2,arg3)
    return err
end

function MatAssemblyBegin(arg1::Mat{Float64},arg2::MatAssemblyType)
    err = ccall((:MatAssemblyBegin,petsc1),PetscErrorCode,(Mat{Float64},MatAssemblyType),arg1,arg2)
    return err
end

function MatAssemblyEnd(arg1::Mat{Float64},arg2::MatAssemblyType)
    err = ccall((:MatAssemblyEnd,petsc1),PetscErrorCode,(Mat{Float64},MatAssemblyType),arg1,arg2)
    return err
end

function MatAssembled(arg1::Mat{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:MatAssembled,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function MatSetOption(arg1::Mat{Float64},arg2::MatOption,arg3::PetscBool)
    err = ccall((:MatSetOption,petsc1),PetscErrorCode,(Mat{Float64},MatOption,PetscBool),arg1,arg2,arg3)
    return err
end

function MatGetOption(arg1::Mat{Float64},arg2::MatOption,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:MatGetOption,petsc1),PetscErrorCode,(Mat{Float64},MatOption,Ptr{PetscBool}),arg1,arg2,arg3)
    return err
end

function MatGetType(arg1::Mat{Float64},arg2::Union(Ptr{MatType},StridedArray{MatType},Ptr{Void}))
    (arg2_,tmp) = symbol_get_before(arg2)
    err = ccall((:MatGetType,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Ptr{Uint8}}),arg1,arg2_)
    symbol_get_after(arg2_,arg2)
    return err
end

function MatGetValues(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:MatGetValues,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatGetRow(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg5::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    err = ccall((:MatGetRow,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function MatRestoreRow(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg5::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    err = ccall((:MatRestoreRow,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function MatGetRowUpperTriangular(arg1::Mat{Float64})
    err = ccall((:MatGetRowUpperTriangular,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
    return err
end

function MatRestoreRowUpperTriangular(arg1::Mat{Float64})
    err = ccall((:MatRestoreRowUpperTriangular,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
    return err
end

function MatGetColumn(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg5::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    err = ccall((:MatGetColumn,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function MatRestoreColumn(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg5::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    err = ccall((:MatRestoreColumn,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function MatGetColumnVector(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Integer)
    err = ccall((:MatGetColumnVector,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},PetscInt),arg1,arg2,arg3)
    return err
end

function MatSeqAIJGetArray(arg1::Mat{Float64},arg2::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    err = ccall((:MatSeqAIJGetArray,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Ptr{Float64}}),arg1,arg2)
    return err
end

function MatSeqAIJRestoreArray(arg1::Mat{Float64},arg2::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    err = ccall((:MatSeqAIJRestoreArray,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Ptr{Float64}}),arg1,arg2)
    return err
end

function MatSeqAIJGetMaxRowNonzeros(arg1::Mat{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatSeqAIJGetMaxRowNonzeros,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function MatSeqAIJSetValuesLocalFast(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::InsertMode)
    err = ccall((:MatSeqAIJSetValuesLocalFast,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
    return err
end

function MatDenseGetArray(arg1::Mat{Float64},arg2::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    err = ccall((:MatDenseGetArray,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Ptr{Float64}}),arg1,arg2)
    return err
end

function MatDenseRestoreArray(arg1::Mat{Float64},arg2::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    err = ccall((:MatDenseRestoreArray,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Ptr{Float64}}),arg1,arg2)
    return err
end

function MatGetBlockSize(arg1::Mat{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatGetBlockSize,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function MatSetBlockSize(arg1::Mat{Float64},arg2::Integer)
    err = ccall((:MatSetBlockSize,petsc1),PetscErrorCode,(Mat{Float64},PetscInt),arg1,arg2)
    return err
end

function MatGetBlockSizes(arg1::Mat{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatGetBlockSizes,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function MatSetBlockSizes(arg1::Mat{Float64},arg2::Integer,arg3::Integer)
    err = ccall((:MatSetBlockSizes,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,PetscInt),arg1,arg2,arg3)
    return err
end

function MatSetBlockSizesFromMats(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64})
    err = ccall((:MatSetBlockSizesFromMats,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3)
    return err
end

function MatSetNThreads(arg1::Mat{Float64},arg2::Integer)
    err = ccall((:MatSetNThreads,petsc1),PetscErrorCode,(Mat{Float64},PetscInt),arg1,arg2)
    return err
end

function MatGetNThreads(arg1::Mat{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatGetNThreads,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function MatMult(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:MatMult,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function MatMultDiagonalBlock(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:MatMultDiagonalBlock,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function MatMultAdd(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64})
    err = ccall((:MatMultAdd,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function MatMultTranspose(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:MatMultTranspose,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function MatMultHermitianTranspose(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:MatMultHermitianTranspose,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function MatIsTranspose(arg1::Mat{Float64},arg2::Mat{Float64},PetscReal::Integer,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:MatIsTranspose,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Cint,Ptr{PetscBool}),arg1,arg2,PetscReal,arg3)
    return err
end

function MatIsHermitianTranspose(arg1::Mat{Float64},arg2::Mat{Float64},PetscReal::Integer,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:MatIsHermitianTranspose,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Cint,Ptr{PetscBool}),arg1,arg2,PetscReal,arg3)
    return err
end

function MatMultTransposeAdd(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64})
    err = ccall((:MatMultTransposeAdd,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function MatMultHermitianTransposeAdd(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64})
    err = ccall((:MatMultHermitianTransposeAdd,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function MatMultConstrained(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:MatMultConstrained,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function MatMultTransposeConstrained(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:MatMultTransposeConstrained,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function MatMatSolve(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64})
    err = ccall((:MatMatSolve,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3)
    return err
end

function MatResidual(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64})
    err = ccall((:MatResidual,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function MatConvert(arg1::Mat{Float64},arg2::MatType,arg3::MatReuse,arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatConvert,petsc1),PetscErrorCode,(Mat{Float64},Cstring,MatReuse,Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function MatDuplicate(arg1::Mat{Float64},arg2::MatDuplicateOption,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatDuplicate,petsc1),PetscErrorCode,(Mat{Float64},MatDuplicateOption,Ptr{Mat{Float64}}),arg1,arg2,arg3)
    return err
end

function MatCopy(arg1::Mat{Float64},arg2::Mat{Float64},arg3::MatStructure)
    err = ccall((:MatCopy,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},MatStructure),arg1,arg2,arg3)
    return err
end

function MatView(arg1::Mat{Float64},arg2::PetscViewer{Float64})
    err = ccall((:MatView,petsc1),PetscErrorCode,(Mat{Float64},PetscViewer{Float64}),arg1,arg2)
    return err
end

function MatIsSymmetric(arg1::Mat{Float64},PetscReal::Integer,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:MatIsSymmetric,petsc1),PetscErrorCode,(Mat{Float64},Cint,Ptr{PetscBool}),arg1,PetscReal,arg2)
    return err
end

function MatIsStructurallySymmetric(arg1::Mat{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:MatIsStructurallySymmetric,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function MatIsHermitian(arg1::Mat{Float64},PetscReal::Integer,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:MatIsHermitian,petsc1),PetscErrorCode,(Mat{Float64},Cint,Ptr{PetscBool}),arg1,PetscReal,arg2)
    return err
end

function MatIsSymmetricKnown(arg1::Mat{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:MatIsSymmetricKnown,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3)
    return err
end

function MatIsHermitianKnown(arg1::Mat{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:MatIsHermitianKnown,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3)
    return err
end

function MatMissingDiagonal(arg1::Mat{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatMissingDiagonal,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscBool},Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function MatLoad(arg1::Mat{Float64},arg2::PetscViewer{Float64})
    err = ccall((:MatLoad,petsc1),PetscErrorCode,(Mat{Float64},PetscViewer{Float64}),arg1,arg2)
    return err
end

function MatGetRowIJ(arg1::Mat{Float64},arg2::Integer,arg3::PetscBool,arg4::PetscBool,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg7::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg8::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:MatGetRowIJ,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,PetscBool,PetscBool,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
    return err
end

function MatRestoreRowIJ(arg1::Mat{Float64},arg2::Integer,arg3::PetscBool,arg4::PetscBool,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg7::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg8::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:MatRestoreRowIJ,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,PetscBool,PetscBool,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
    return err
end

function MatGetColumnIJ(arg1::Mat{Float64},arg2::Integer,arg3::PetscBool,arg4::PetscBool,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg7::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg8::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:MatGetColumnIJ,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,PetscBool,PetscBool,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
    return err
end

function MatRestoreColumnIJ(arg1::Mat{Float64},arg2::Integer,arg3::PetscBool,arg4::PetscBool,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg7::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg8::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:MatRestoreColumnIJ,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,PetscBool,PetscBool,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
    return err
end

function MatGetInfo(arg1::Mat{Float64},arg2::MatInfoType,arg3::Union(Ptr{MatInfo},StridedArray{MatInfo},Ptr{Void}))
    err = ccall((:MatGetInfo,petsc1),PetscErrorCode,(Mat{Float64},MatInfoType,Ptr{MatInfo}),arg1,arg2,arg3)
    return err
end

function MatGetDiagonal(arg1::Mat{Float64},arg2::Vec{Float64})
    err = ccall((:MatGetDiagonal,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64}),arg1,arg2)
    return err
end

function MatGetRowMax(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatGetRowMax,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function MatGetRowMin(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatGetRowMin,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function MatGetRowMaxAbs(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatGetRowMaxAbs,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function MatGetRowMinAbs(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatGetRowMinAbs,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function MatGetRowSum(arg1::Mat{Float64},arg2::Vec{Float64})
    err = ccall((:MatGetRowSum,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64}),arg1,arg2)
    return err
end

function MatTranspose(arg1::Mat{Float64},arg2::MatReuse,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatTranspose,petsc1),PetscErrorCode,(Mat{Float64},MatReuse,Ptr{Mat{Float64}}),arg1,arg2,arg3)
    return err
end

function MatHermitianTranspose(arg1::Mat{Float64},arg2::MatReuse,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatHermitianTranspose,petsc1),PetscErrorCode,(Mat{Float64},MatReuse,Ptr{Mat{Float64}}),arg1,arg2,arg3)
    return err
end

function MatPermute(arg1::Mat{Float64},arg2::IS{Float64},arg3::IS{Float64},arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatPermute,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},IS{Float64},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function MatDiagonalScale(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:MatDiagonalScale,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function MatDiagonalSet(arg1::Mat{Float64},arg2::Vec{Float64},arg3::InsertMode)
    err = ccall((:MatDiagonalSet,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},InsertMode),arg1,arg2,arg3)
    return err
end

function MatEqual(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:MatEqual,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Ptr{PetscBool}),arg1,arg2,arg3)
    return err
end

function MatMultEqual(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Integer,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:MatMultEqual,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},PetscInt,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
    return err
end

function MatMultAddEqual(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Integer,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:MatMultAddEqual,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},PetscInt,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
    return err
end

function MatMultTransposeEqual(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Integer,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:MatMultTransposeEqual,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},PetscInt,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
    return err
end

function MatMultTransposeAddEqual(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Integer,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:MatMultTransposeAddEqual,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},PetscInt,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
    return err
end

function MatNorm(arg1::Mat{Float64},arg2::NormType,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:MatNorm,petsc1),PetscErrorCode,(Mat{Float64},NormType,Ptr{Cint}),arg1,arg2,arg3)
    return err
end

function MatGetColumnNorms(arg1::Mat{Float64},arg2::NormType,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:MatGetColumnNorms,petsc1),PetscErrorCode,(Mat{Float64},NormType,Ptr{Cint}),arg1,arg2,arg3)
    return err
end

function MatZeroEntries(arg1::Mat{Float64})
    err = ccall((:MatZeroEntries,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
    return err
end

function MatZeroRows(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Float64,arg5::Vec{Float64},arg6::Vec{Float64})
    err = ccall((:MatZeroRows,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatZeroRowsIS(arg1::Mat{Float64},arg2::IS{Float64},arg3::Float64,arg4::Vec{Float64},arg5::Vec{Float64})
    err = ccall((:MatZeroRowsIS,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function MatZeroRowsStencil(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{MatStencil},StridedArray{MatStencil},Ptr{Void}),arg4::Float64,arg5::Vec{Float64},arg6::Vec{Float64})
    err = ccall((:MatZeroRowsStencil,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{MatStencil},Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatZeroRowsColumnsStencil(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{MatStencil},StridedArray{MatStencil},Ptr{Void}),arg4::Float64,arg5::Vec{Float64},arg6::Vec{Float64})
    err = ccall((:MatZeroRowsColumnsStencil,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{MatStencil},Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatZeroRowsColumns(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Float64,arg5::Vec{Float64},arg6::Vec{Float64})
    err = ccall((:MatZeroRowsColumns,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatZeroRowsColumnsIS(arg1::Mat{Float64},arg2::IS{Float64},arg3::Float64,arg4::Vec{Float64},arg5::Vec{Float64})
    err = ccall((:MatZeroRowsColumnsIS,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function MatGetSize(arg1::Mat{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatGetSize,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function MatGetLocalSize(arg1::Mat{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatGetLocalSize,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function MatGetOwnershipRange(arg1::Mat{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatGetOwnershipRange,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function MatGetOwnershipRanges(arg1::Mat{Float64},arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:MatGetOwnershipRanges,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Ptr{PetscInt}}),arg1,arg2)
    return err
end

function MatGetOwnershipRangeColumn(arg1::Mat{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatGetOwnershipRangeColumn,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function MatGetOwnershipRangesColumn(arg1::Mat{Float64},arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:MatGetOwnershipRangesColumn,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Ptr{PetscInt}}),arg1,arg2)
    return err
end

function MatGetOwnershipIS(arg1::Mat{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:MatGetOwnershipIS,petsc1),PetscErrorCode,(Mat{Float64},Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3)
    return err
end

function MatGetSubMatrices(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg5::MatReuse,arg6::Union(Ptr{Ptr{Mat{Float64}}},StridedArray{Ptr{Mat{Float64}}},Ptr{Void}))
    err = ccall((:MatGetSubMatrices,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{IS{Float64}},Ptr{IS{Float64}},MatReuse,Ptr{Ptr{Mat{Float64}}}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatGetSubMatricesMPI(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg5::MatReuse,arg6::Union(Ptr{Ptr{Mat{Float64}}},StridedArray{Ptr{Mat{Float64}}},Ptr{Void}))
    err = ccall((:MatGetSubMatricesMPI,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{IS{Float64}},Ptr{IS{Float64}},MatReuse,Ptr{Ptr{Mat{Float64}}}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatDestroyMatrices(arg1::Integer,arg2::Union(Ptr{Ptr{Mat{Float64}}},StridedArray{Ptr{Mat{Float64}}},Ptr{Void}))
    err = ccall((:MatDestroyMatrices,petsc1),PetscErrorCode,(PetscInt,Ptr{Ptr{Mat{Float64}}}),arg1,arg2)
    return err
end

function MatGetSubMatrix(arg1::Mat{Float64},arg2::IS{Float64},arg3::IS{Float64},arg4::MatReuse,arg5::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatGetSubMatrix,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},IS{Float64},MatReuse,Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function MatGetLocalSubMatrix(arg1::Mat{Float64},arg2::IS{Float64},arg3::IS{Float64},arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatGetLocalSubMatrix,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},IS{Float64},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function MatRestoreLocalSubMatrix(arg1::Mat{Float64},arg2::IS{Float64},arg3::IS{Float64},arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatRestoreLocalSubMatrix,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},IS{Float64},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function MatGetSeqNonzeroStructure(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatGetSeqNonzeroStructure,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2)
    return err
end

function MatDestroySeqNonzeroStructure(arg1::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatDestroySeqNonzeroStructure,petsc1),PetscErrorCode,(Ptr{Mat{Float64}},),arg1)
    return err
end

function MatCreateMPIAIJSumSeqAIJ(arg1::MPI_Comm,arg2::Mat{Float64},arg3::Integer,arg4::Integer,arg5::MatReuse,arg6::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateMPIAIJSumSeqAIJ,petsc1),PetscErrorCode,(comm_type,Mat{Float64},PetscInt,PetscInt,MatReuse,Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatCreateMPIAIJSumSeqAIJSymbolic(arg1::MPI_Comm,arg2::Mat{Float64},arg3::Integer,arg4::Integer,arg5::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateMPIAIJSumSeqAIJSymbolic,petsc1),PetscErrorCode,(comm_type,Mat{Float64},PetscInt,PetscInt,Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5)
    return err
end

function MatCreateMPIAIJSumSeqAIJNumeric(arg1::Mat{Float64},arg2::Mat{Float64})
    err = ccall((:MatCreateMPIAIJSumSeqAIJNumeric,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64}),arg1,arg2)
    return err
end

function MatMPIAIJGetLocalMat(arg1::Mat{Float64},arg2::MatReuse,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatMPIAIJGetLocalMat,petsc1),PetscErrorCode,(Mat{Float64},MatReuse,Ptr{Mat{Float64}}),arg1,arg2,arg3)
    return err
end

function MatMPIAIJGetLocalMatCondensed(arg1::Mat{Float64},arg2::MatReuse,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg5::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatMPIAIJGetLocalMatCondensed,petsc1),PetscErrorCode,(Mat{Float64},MatReuse,Ptr{IS{Float64}},Ptr{IS{Float64}},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function MatGetBrowsOfAcols(arg1::Mat{Float64},arg2::Mat{Float64},arg3::MatReuse,arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg5::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg6::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatGetBrowsOfAcols,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},MatReuse,Ptr{IS{Float64}},Ptr{IS{Float64}},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatGetGhosts(arg1::Mat{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:MatGetGhosts,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3)
    return err
end

function MatIncreaseOverlap(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Integer)
    err = ccall((:MatIncreaseOverlap,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{IS{Float64}},PetscInt),arg1,arg2,arg3,arg4)
    return err
end

function MatMatMult(arg1::Mat{Float64},arg2::Mat{Float64},arg3::MatReuse,PetscReal::Integer,arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatMatMult,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},MatReuse,Cint,Ptr{Mat{Float64}}),arg1,arg2,arg3,PetscReal,arg4)
    return err
end

function MatMatMultSymbolic(arg1::Mat{Float64},arg2::Mat{Float64},PetscReal::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatMatMultSymbolic,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Cint,Ptr{Mat{Float64}}),arg1,arg2,PetscReal,arg3)
    return err
end

function MatMatMultNumeric(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64})
    err = ccall((:MatMatMultNumeric,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3)
    return err
end

function MatMatMatMult(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64},arg4::MatReuse,PetscReal::Integer,arg5::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatMatMatMult,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64},MatReuse,Cint,Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4,PetscReal,arg5)
    return err
end

function MatMatMatMultSymbolic(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64},PetscReal::Integer,arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatMatMatMultSymbolic,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64},Cint,Ptr{Mat{Float64}}),arg1,arg2,arg3,PetscReal,arg4)
    return err
end

function MatMatMatMultNumeric(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64},arg4::Mat{Float64})
    err = ccall((:MatMatMatMultNumeric,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function MatPtAP(arg1::Mat{Float64},arg2::Mat{Float64},arg3::MatReuse,PetscReal::Integer,arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatPtAP,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},MatReuse,Cint,Ptr{Mat{Float64}}),arg1,arg2,arg3,PetscReal,arg4)
    return err
end

function MatPtAPSymbolic(arg1::Mat{Float64},arg2::Mat{Float64},PetscReal::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatPtAPSymbolic,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Cint,Ptr{Mat{Float64}}),arg1,arg2,PetscReal,arg3)
    return err
end

function MatPtAPNumeric(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64})
    err = ccall((:MatPtAPNumeric,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3)
    return err
end

function MatRARt(arg1::Mat{Float64},arg2::Mat{Float64},arg3::MatReuse,PetscReal::Integer,arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatRARt,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},MatReuse,Cint,Ptr{Mat{Float64}}),arg1,arg2,arg3,PetscReal,arg4)
    return err
end

function MatRARtSymbolic(arg1::Mat{Float64},arg2::Mat{Float64},PetscReal::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatRARtSymbolic,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Cint,Ptr{Mat{Float64}}),arg1,arg2,PetscReal,arg3)
    return err
end

function MatRARtNumeric(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64})
    err = ccall((:MatRARtNumeric,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3)
    return err
end

function MatTransposeMatMult(arg1::Mat{Float64},arg2::Mat{Float64},arg3::MatReuse,PetscReal::Integer,arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatTransposeMatMult,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},MatReuse,Cint,Ptr{Mat{Float64}}),arg1,arg2,arg3,PetscReal,arg4)
    return err
end

function MatTransposetMatMultSymbolic(arg1::Mat{Float64},arg2::Mat{Float64},PetscReal::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatTransposetMatMultSymbolic,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Cint,Ptr{Mat{Float64}}),arg1,arg2,PetscReal,arg3)
    return err
end

function MatTransposetMatMultNumeric(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64})
    err = ccall((:MatTransposetMatMultNumeric,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3)
    return err
end

function MatMatTransposeMult(arg1::Mat{Float64},arg2::Mat{Float64},arg3::MatReuse,PetscReal::Integer,arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatMatTransposeMult,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},MatReuse,Cint,Ptr{Mat{Float64}}),arg1,arg2,arg3,PetscReal,arg4)
    return err
end

function MatMatTransposeMultSymbolic(arg1::Mat{Float64},arg2::Mat{Float64},PetscReal::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatMatTransposeMultSymbolic,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Cint,Ptr{Mat{Float64}}),arg1,arg2,PetscReal,arg3)
    return err
end

function MatMatTransposeMultNumeric(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64})
    err = ccall((:MatMatTransposeMultNumeric,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3)
    return err
end

function MatAXPY(arg1::Mat{Float64},arg2::Float64,arg3::Mat{Float64},arg4::MatStructure)
    err = ccall((:MatAXPY,petsc1),PetscErrorCode,(Mat{Float64},Float64,Mat{Float64},MatStructure),arg1,arg2,arg3,arg4)
    return err
end

function MatAYPX(arg1::Mat{Float64},arg2::Float64,arg3::Mat{Float64},arg4::MatStructure)
    err = ccall((:MatAYPX,petsc1),PetscErrorCode,(Mat{Float64},Float64,Mat{Float64},MatStructure),arg1,arg2,arg3,arg4)
    return err
end

function MatScale(arg1::Mat{Float64},arg2::Float64)
    err = ccall((:MatScale,petsc1),PetscErrorCode,(Mat{Float64},Float64),arg1,arg2)
    return err
end

function MatShift(arg1::Mat{Float64},arg2::Float64)
    err = ccall((:MatShift,petsc1),PetscErrorCode,(Mat{Float64},Float64),arg1,arg2)
    return err
end

function MatSetLocalToGlobalMapping(arg1::Mat{Float64},arg2::ISLocalToGlobalMapping{Float64},arg3::ISLocalToGlobalMapping{Float64})
    err = ccall((:MatSetLocalToGlobalMapping,petsc1),PetscErrorCode,(Mat{Float64},ISLocalToGlobalMapping{Float64},ISLocalToGlobalMapping{Float64}),arg1,arg2,arg3)
    return err
end

function MatGetLocalToGlobalMapping(arg1::Mat{Float64},arg2::Union(Ptr{ISLocalToGlobalMapping{Float64}},StridedArray{ISLocalToGlobalMapping{Float64}},Ptr{Void}),arg3::Union(Ptr{ISLocalToGlobalMapping{Float64}},StridedArray{ISLocalToGlobalMapping{Float64}},Ptr{Void}))
    err = ccall((:MatGetLocalToGlobalMapping,petsc1),PetscErrorCode,(Mat{Float64},Ptr{ISLocalToGlobalMapping{Float64}},Ptr{ISLocalToGlobalMapping{Float64}}),arg1,arg2,arg3)
    return err
end

function MatGetLayouts(arg1::Mat{Float64},arg2::Union(Ptr{PetscLayout{Float64}},StridedArray{PetscLayout{Float64}},Ptr{Void}),arg3::Union(Ptr{PetscLayout{Float64}},StridedArray{PetscLayout{Float64}},Ptr{Void}))
    err = ccall((:MatGetLayouts,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscLayout{Float64}{Float64}},Ptr{PetscLayout{Float64}{Float64}}),arg1,arg2,arg3)
    return err
end

function MatZeroRowsLocal(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Float64,arg5::Vec{Float64},arg6::Vec{Float64})
    err = ccall((:MatZeroRowsLocal,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatZeroRowsLocalIS(arg1::Mat{Float64},arg2::IS{Float64},arg3::Float64,arg4::Vec{Float64},arg5::Vec{Float64})
    err = ccall((:MatZeroRowsLocalIS,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function MatZeroRowsColumnsLocal(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Float64,arg5::Vec{Float64},arg6::Vec{Float64})
    err = ccall((:MatZeroRowsColumnsLocal,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatZeroRowsColumnsLocalIS(arg1::Mat{Float64},arg2::IS{Float64},arg3::Float64,arg4::Vec{Float64},arg5::Vec{Float64})
    err = ccall((:MatZeroRowsColumnsLocalIS,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function MatSetValuesLocal(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::InsertMode)
    err = ccall((:MatSetValuesLocal,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
    return err
end

function MatSetValuesBlockedLocal(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::InsertMode)
    err = ccall((:MatSetValuesBlockedLocal,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
    return err
end

function MatStashSetInitialSize(arg1::Mat{Float64},arg2::Integer,arg3::Integer)
    err = ccall((:MatStashSetInitialSize,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,PetscInt),arg1,arg2,arg3)
    return err
end

function MatStashGetInfo(arg1::Mat{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatStashGetInfo,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function MatInterpolate(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:MatInterpolate,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function MatInterpolateAdd(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64})
    err = ccall((:MatInterpolateAdd,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function MatRestrict(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:MatRestrict,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function MatCreateVecs(arg1::Mat{Float64},arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:MatCreateVecs,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Vec{Float64}},Ptr{Vec{Float64}}),arg1,arg2,arg3)
    return err
end

function MatGetMultiProcBlock(arg1::Mat{Float64},arg2::MPI_Comm,arg3::MatReuse,arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatGetMultiProcBlock,petsc1),PetscErrorCode,(Mat{Float64},comm_type,MatReuse,Ptr{Mat{Float64}}),arg1,arg2.val,arg3,arg4)
    return err
end

function MatFindZeroDiagonals(arg1::Mat{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:MatFindZeroDiagonals,petsc1),PetscErrorCode,(Mat{Float64},Ptr{IS{Float64}}),arg1,arg2)
    return err
end

function MatFindOffBlockDiagonalEntries(arg1::Mat{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:MatFindOffBlockDiagonalEntries,petsc1),PetscErrorCode,(Mat{Float64},Ptr{IS{Float64}}),arg1,arg2)
    return err
end

function MatCreateMPIMatConcatenateSeqMat(arg1::MPI_Comm,arg2::Mat{Float64},arg3::Integer,arg4::MatReuse,arg5::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateMPIMatConcatenateSeqMat,petsc1),PetscErrorCode,(comm_type,Mat{Float64},PetscInt,MatReuse,Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5)
    return err
end

function MatInodeAdjustForInodes(arg1::Mat{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:MatInodeAdjustForInodes,petsc1),PetscErrorCode,(Mat{Float64},Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3)
    return err
end

function MatInodeGetInodeSizes(arg1::Mat{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatInodeGetInodeSizes,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
    return err
end

function MatSeqAIJSetColumnIndices(arg1::Mat{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatSeqAIJSetColumnIndices,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function MatSeqBAIJSetColumnIndices(arg1::Mat{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatSeqBAIJSetColumnIndices,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function MatCreateSeqAIJWithArrays(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateSeqAIJWithArrays,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
    return err
end

function MatCreateSeqBAIJWithArrays(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg8::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateSeqBAIJWithArrays,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
    return err
end

function MatCreateSeqSBAIJWithArrays(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg8::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateSeqSBAIJWithArrays,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
    return err
end

function MatCreateSeqAIJFromTriple(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg8::Integer,arg9::PetscBool)
    err = ccall((:MatCreateSeqAIJFromTriple,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64},Ptr{Mat{Float64}},PetscInt,PetscBool),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
    return err
end

function MatSeqBAIJSetPreallocation(arg1::Mat{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatSeqBAIJSetPreallocation,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
    return err
end

function MatSeqSBAIJSetPreallocation(arg1::Mat{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatSeqSBAIJSetPreallocation,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
    return err
end

function MatSeqAIJSetPreallocation(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatSeqAIJSetPreallocation,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function MatMPIBAIJSetPreallocation(arg1::Mat{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Integer,arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatMPIBAIJSetPreallocation,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatMPISBAIJSetPreallocation(arg1::Mat{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Integer,arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatMPISBAIJSetPreallocation,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatMPIAIJSetPreallocation(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatMPIAIJSetPreallocation,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function MatSeqAIJSetPreallocationCSR(arg1::Mat{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:MatSeqAIJSetPreallocationCSR,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function MatSeqBAIJSetPreallocationCSR(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:MatSeqBAIJSetPreallocationCSR,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function MatMPIAIJSetPreallocationCSR(arg1::Mat{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:MatMPIAIJSetPreallocationCSR,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function MatMPIBAIJSetPreallocationCSR(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:MatMPIBAIJSetPreallocationCSR,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function MatMPIAdjSetPreallocation(arg1::Mat{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatMPIAdjSetPreallocation,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
    return err
end

function MatMPIDenseSetPreallocation(arg1::Mat{Float64},arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:MatMPIDenseSetPreallocation,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Float64}),arg1,arg2)
    return err
end

function MatSeqDenseSetPreallocation(arg1::Mat{Float64},arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:MatSeqDenseSetPreallocation,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Float64}),arg1,arg2)
    return err
end

function MatMPIAIJGetSeqAIJ(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg4::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:MatMPIAIJGetSeqAIJ,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}},Ptr{Mat{Float64}},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4)
    return err
end

function MatMPIBAIJGetSeqBAIJ(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg4::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:MatMPIBAIJGetSeqBAIJ,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}},Ptr{Mat{Float64}},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4)
    return err
end

function MatMPIAdjCreateNonemptySubcommMat(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatMPIAdjCreateNonemptySubcommMat,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2)
    return err
end

function MatISSetPreallocation(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatISSetPreallocation,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function MatSeqDenseSetLDA(arg1::Mat{Float64},arg2::Integer)
    err = ccall((:MatSeqDenseSetLDA,petsc1),PetscErrorCode,(Mat{Float64},PetscInt),arg1,arg2)
    return err
end

function MatDenseGetLocalMatrix(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatDenseGetLocalMatrix,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2)
    return err
end

function MatStoreValues(arg1::Mat{Float64})
    err = ccall((:MatStoreValues,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
    return err
end

function MatRetrieveValues(arg1::Mat{Float64})
    err = ccall((:MatRetrieveValues,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
    return err
end

function MatDAADSetCtx(arg1::Mat{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:MatDAADSetCtx,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Void}),arg1,arg2)
    return err
end

function MatFindNonzeroRows(arg1::Mat{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:MatFindNonzeroRows,petsc1),PetscErrorCode,(Mat{Float64},Ptr{IS{Float64}}),arg1,arg2)
    return err
end

function MatGetOrdering(arg1::Mat{Float64},arg2::MatOrderingType,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:MatGetOrdering,petsc1),PetscErrorCode,(Mat{Float64},Cstring,Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

#= skipping function with undefined symbols: 
 function MatGetOrderingList(arg0::Type{Float64},arg1::Union(Ptr{PetscFunctionList},StridedArray{PetscFunctionList},Ptr{Void}))
    ccall((:MatGetOrderingList,petsc1),PetscErrorCode,(Ptr{PetscFunctionList},),arg1)
end 
=#
function MatOrderingRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:MatOrderingRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

function MatReorderForNonzeroDiagonal(arg1::Mat{Float64},PetscReal::Integer,arg2::IS{Float64},arg3::IS{Float64})
    err = ccall((:MatReorderForNonzeroDiagonal,petsc1),PetscErrorCode,(Mat{Float64},Cint,IS{Float64},IS{Float64}),arg1,PetscReal,arg2,arg3)
    return err
end

function MatCreateLaplacian(arg1::Mat{Float64},PetscReal::Integer,arg2::PetscBool,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateLaplacian,petsc1),PetscErrorCode,(Mat{Float64},Cint,PetscBool,Ptr{Mat{Float64}}),arg1,PetscReal,arg2,arg3)
    return err
end

function MatFactorInfoInitialize(arg0::Type{Float64},arg1::Union(Ptr{MatFactorInfo},StridedArray{MatFactorInfo},Ptr{Void}))
    err = ccall((:MatFactorInfoInitialize,petsc1),PetscErrorCode,(Ptr{MatFactorInfo},),arg1)
    return err
end

function MatCholeskyFactor(arg1::Mat{Float64},arg2::IS{Float64},arg3::Union(Ptr{MatFactorInfo},StridedArray{MatFactorInfo},Ptr{Void}))
    err = ccall((:MatCholeskyFactor,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},Ptr{MatFactorInfo}),arg1,arg2,arg3)
    return err
end

function MatCholeskyFactorSymbolic(arg1::Mat{Float64},arg2::Mat{Float64},arg3::IS{Float64},arg4::Union(Ptr{MatFactorInfo},StridedArray{MatFactorInfo},Ptr{Void}))
    err = ccall((:MatCholeskyFactorSymbolic,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},IS{Float64},Ptr{MatFactorInfo}),arg1,arg2,arg3,arg4)
    return err
end

function MatCholeskyFactorNumeric(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Union(Ptr{MatFactorInfo},StridedArray{MatFactorInfo},Ptr{Void}))
    err = ccall((:MatCholeskyFactorNumeric,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Ptr{MatFactorInfo}),arg1,arg2,arg3)
    return err
end

function MatLUFactor(arg1::Mat{Float64},arg2::IS{Float64},arg3::IS{Float64},arg4::Union(Ptr{MatFactorInfo},StridedArray{MatFactorInfo},Ptr{Void}))
    err = ccall((:MatLUFactor,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},IS{Float64},Ptr{MatFactorInfo}),arg1,arg2,arg3,arg4)
    return err
end

function MatILUFactor(arg1::Mat{Float64},arg2::IS{Float64},arg3::IS{Float64},arg4::Union(Ptr{MatFactorInfo},StridedArray{MatFactorInfo},Ptr{Void}))
    err = ccall((:MatILUFactor,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},IS{Float64},Ptr{MatFactorInfo}),arg1,arg2,arg3,arg4)
    return err
end

function MatLUFactorSymbolic(arg1::Mat{Float64},arg2::Mat{Float64},arg3::IS{Float64},arg4::IS{Float64},arg5::Union(Ptr{MatFactorInfo},StridedArray{MatFactorInfo},Ptr{Void}))
    err = ccall((:MatLUFactorSymbolic,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},IS{Float64},IS{Float64},Ptr{MatFactorInfo}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function MatILUFactorSymbolic(arg1::Mat{Float64},arg2::Mat{Float64},arg3::IS{Float64},arg4::IS{Float64},arg5::Union(Ptr{MatFactorInfo},StridedArray{MatFactorInfo},Ptr{Void}))
    err = ccall((:MatILUFactorSymbolic,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},IS{Float64},IS{Float64},Ptr{MatFactorInfo}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function MatICCFactorSymbolic(arg1::Mat{Float64},arg2::Mat{Float64},arg3::IS{Float64},arg4::Union(Ptr{MatFactorInfo},StridedArray{MatFactorInfo},Ptr{Void}))
    err = ccall((:MatICCFactorSymbolic,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},IS{Float64},Ptr{MatFactorInfo}),arg1,arg2,arg3,arg4)
    return err
end

function MatICCFactor(arg1::Mat{Float64},arg2::IS{Float64},arg3::Union(Ptr{MatFactorInfo},StridedArray{MatFactorInfo},Ptr{Void}))
    err = ccall((:MatICCFactor,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},Ptr{MatFactorInfo}),arg1,arg2,arg3)
    return err
end

function MatLUFactorNumeric(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Union(Ptr{MatFactorInfo},StridedArray{MatFactorInfo},Ptr{Void}))
    err = ccall((:MatLUFactorNumeric,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Ptr{MatFactorInfo}),arg1,arg2,arg3)
    return err
end

function MatGetInertia(arg1::Mat{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatGetInertia,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
    return err
end

function MatSolve(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:MatSolve,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function MatForwardSolve(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:MatForwardSolve,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function MatBackwardSolve(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:MatBackwardSolve,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function MatSolveAdd(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64})
    err = ccall((:MatSolveAdd,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function MatSolveTranspose(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:MatSolveTranspose,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function MatSolveTransposeAdd(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64})
    err = ccall((:MatSolveTransposeAdd,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
    return err
end

#= skipping function with undefined symbols: 
 function MatSolves(arg1::Mat{Float64},arg2::Vecs,arg3::Vecs)
    ccall((:MatSolves,petsc1),PetscErrorCode,(Mat{Float64},Vecs,Vecs),arg1,arg2,arg3)
end 
=#
function MatSetUnfactored(arg1::Mat{Float64})
    err = ccall((:MatSetUnfactored,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
    return err
end

function MatSOR(arg1::Mat{Float64},arg2::Vec{Float64},PetscReal::Integer,arg3::MatSORType,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Vec{Float64})
    err = ccall((:MatSOR,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Cint,MatSORType,Cint,PetscInt,PetscInt,Vec{Float64}),arg1,arg2,PetscReal,arg3,arg4,arg5,arg6,arg7)
    return err
end

#= skipping function with undefined symbols: 
 function MatColoringCreate(arg1::Mat{Float64},arg2::Union(Ptr{MatColoring},StridedArray{MatColoring},Ptr{Void}))
    ccall((:MatColoringCreate,petsc1),PetscErrorCode,(Mat{Float64},Ptr{MatColoring}),arg1,arg2)
end 
=#
function MatColoringGetDegrees(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatColoringGetDegrees,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

#= skipping function with undefined symbols: 
 function MatColoringDestroy(arg0::Type{Float64},arg1::Union(Ptr{MatColoring},StridedArray{MatColoring},Ptr{Void}))
    ccall((:MatColoringDestroy,petsc1),PetscErrorCode,(Ptr{MatColoring},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function MatColoringView(arg1::MatColoring,arg2::PetscViewer{Float64})
    ccall((:MatColoringView,petsc1),PetscErrorCode,(MatColoring,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatColoringSetType(arg0::Type{Float64},arg1::MatColoring,arg2::MatColoringType)
    ccall((:MatColoringSetType,petsc1),PetscErrorCode,(MatColoring,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatColoringSetFromOptions(arg0::Type{Float64},arg1::MatColoring)
    ccall((:MatColoringSetFromOptions,petsc1),PetscErrorCode,(MatColoring,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function MatColoringSetDistance(arg0::Type{Float64},arg1::MatColoring,arg2::Integer)
    ccall((:MatColoringSetDistance,petsc1),PetscErrorCode,(MatColoring,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatColoringGetDistance(arg0::Type{Float64},arg1::MatColoring,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:MatColoringGetDistance,petsc1),PetscErrorCode,(MatColoring,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatColoringSetMaxColors(arg0::Type{Float64},arg1::MatColoring,arg2::Integer)
    ccall((:MatColoringSetMaxColors,petsc1),PetscErrorCode,(MatColoring,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatColoringGetMaxColors(arg0::Type{Float64},arg1::MatColoring,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:MatColoringGetMaxColors,petsc1),PetscErrorCode,(MatColoring,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatColoringApply(arg1::MatColoring,arg2::Union(Ptr{ISColoring{Float64}},StridedArray{ISColoring{Float64}},Ptr{Void}))
    ccall((:MatColoringApply,petsc1),PetscErrorCode,(MatColoring,Ptr{ISColoring{Float64}}),arg1,arg2)
end 
=#
function MatColoringRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:MatColoringRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

function MatColoringPatch(arg1::Mat{Float64},arg2::Integer,arg3::Integer,ISColoringValue::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{ISColoring{Float64}},StridedArray{ISColoring{Float64}},Ptr{Void}))
    err = ccall((:MatColoringPatch,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,PetscInt,Ptr{Cint},Ptr{ISColoring{Float64}}),arg1,arg2,arg3,ISColoringValue,arg4)
    return err
end

#= skipping function with undefined symbols: 
 function MatColoringSetWeightType(arg0::Type{Float64},arg1::MatColoring,arg2::MatColoringWeightType)
    ccall((:MatColoringSetWeightType,petsc1),PetscErrorCode,(MatColoring,MatColoringWeightType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatColoringSetWeights(arg0::Type{Float64},arg1::MatColoring,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:MatColoringSetWeights,petsc1),PetscErrorCode,(MatColoring,Ptr{Cint},Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function MatColoringCreateWeights(arg0::Type{Float64},arg1::MatColoring,arg2::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),lperm::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:MatColoringCreateWeights,petsc1),PetscErrorCode,(MatColoring,Ptr{Ptr{Cint}},Ptr{Ptr{PetscInt}}),arg1,arg2,lperm)
end 
=#
#= skipping function with undefined symbols: 
 function MatFDColoringCreate(arg1::Mat{Float64},arg2::ISColoring{Float64},arg3::Union(Ptr{MatFDColoring},StridedArray{MatFDColoring},Ptr{Void}))
    ccall((:MatFDColoringCreate,petsc1),PetscErrorCode,(Mat{Float64},ISColoring{Float64},Ptr{MatFDColoring}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function MatFDColoringDestroy(arg0::Type{Float64},arg1::Union(Ptr{MatFDColoring},StridedArray{MatFDColoring},Ptr{Void}))
    ccall((:MatFDColoringDestroy,petsc1),PetscErrorCode,(Ptr{MatFDColoring},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function MatFDColoringView(arg1::MatFDColoring,arg2::PetscViewer{Float64})
    ccall((:MatFDColoringView,petsc1),PetscErrorCode,(MatFDColoring,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatFDColoringSetFunction(arg0::Type{Float64},arg1::MatFDColoring,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:MatFDColoringSetFunction,petsc1),PetscErrorCode,(MatFDColoring,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function MatFDColoringGetFunction(arg0::Type{Float64},arg1::MatFDColoring,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:MatFDColoringGetFunction,petsc1),PetscErrorCode,(MatFDColoring,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function MatFDColoringSetParameters(arg0::Type{Float64},arg1::MatFDColoring,PetscReal::Integer,arg2::Integer)
    ccall((:MatFDColoringSetParameters,petsc1),PetscErrorCode,(MatFDColoring,Cint,Cint),arg1,PetscReal,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatFDColoringSetFromOptions(arg0::Type{Float64},arg1::MatFDColoring)
    ccall((:MatFDColoringSetFromOptions,petsc1),PetscErrorCode,(MatFDColoring,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function MatFDColoringApply(arg1::Mat{Float64},arg2::MatFDColoring,arg3::Vec{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:MatFDColoringApply,petsc1),PetscErrorCode,(Mat{Float64},MatFDColoring,Vec{Float64},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function MatFDColoringSetF(arg1::MatFDColoring,arg2::Vec{Float64})
    ccall((:MatFDColoringSetF,petsc1),PetscErrorCode,(MatFDColoring,Vec{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatFDColoringGetPerturbedColumns(arg0::Type{Float64},arg1::MatFDColoring,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:MatFDColoringGetPerturbedColumns,petsc1),PetscErrorCode,(MatFDColoring,Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function MatFDColoringSetUp(arg1::Mat{Float64},arg2::ISColoring{Float64},arg3::MatFDColoring)
    ccall((:MatFDColoringSetUp,petsc1),PetscErrorCode,(Mat{Float64},ISColoring{Float64},MatFDColoring),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function MatFDColoringSetBlockSize(arg0::Type{Float64},arg1::MatFDColoring,arg2::Integer,arg3::Integer)
    ccall((:MatFDColoringSetBlockSize,petsc1),PetscErrorCode,(MatFDColoring,PetscInt,PetscInt),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function MatTransposeColoringCreate(arg1::Mat{Float64},arg2::ISColoring{Float64},arg3::Union(Ptr{MatTransposeColoring},StridedArray{MatTransposeColoring},Ptr{Void}))
    ccall((:MatTransposeColoringCreate,petsc1),PetscErrorCode,(Mat{Float64},ISColoring{Float64},Ptr{MatTransposeColoring}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function MatTransColoringApplySpToDen(arg1::MatTransposeColoring,arg2::Mat{Float64},arg3::Mat{Float64})
    ccall((:MatTransColoringApplySpToDen,petsc1),PetscErrorCode,(MatTransposeColoring,Mat{Float64},Mat{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function MatTransColoringApplyDenToSp(arg1::MatTransposeColoring,arg2::Mat{Float64},arg3::Mat{Float64})
    ccall((:MatTransColoringApplyDenToSp,petsc1),PetscErrorCode,(MatTransposeColoring,Mat{Float64},Mat{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function MatTransposeColoringDestroy(arg0::Type{Float64},arg1::Union(Ptr{MatTransposeColoring},StridedArray{MatTransposeColoring},Ptr{Void}))
    ccall((:MatTransposeColoringDestroy,petsc1),PetscErrorCode,(Ptr{MatTransposeColoring},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{MatPartitioning},StridedArray{MatPartitioning},Ptr{Void}))
    ccall((:MatPartitioningCreate,petsc1),PetscErrorCode,(comm_type,Ptr{MatPartitioning}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningSetType(arg0::Type{Float64},arg1::MatPartitioning,arg2::MatPartitioningType)
    ccall((:MatPartitioningSetType,petsc1),PetscErrorCode,(MatPartitioning,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningSetNParts(arg0::Type{Float64},arg1::MatPartitioning,arg2::Integer)
    ccall((:MatPartitioningSetNParts,petsc1),PetscErrorCode,(MatPartitioning,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningSetAdjacency(arg1::MatPartitioning,arg2::Mat{Float64})
    ccall((:MatPartitioningSetAdjacency,petsc1),PetscErrorCode,(MatPartitioning,Mat{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningSetVertexWeights(arg0::Type{Float64},arg1::MatPartitioning,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:MatPartitioningSetVertexWeights,petsc1),PetscErrorCode,(MatPartitioning,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningSetPartitionWeights(arg0::Type{Float64},arg1::MatPartitioning,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:MatPartitioningSetPartitionWeights,petsc1),PetscErrorCode,(MatPartitioning,Ptr{Cint}),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningApply(arg1::MatPartitioning,arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:MatPartitioningApply,petsc1),PetscErrorCode,(MatPartitioning,Ptr{IS{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningDestroy(arg0::Type{Float64},arg1::Union(Ptr{MatPartitioning},StridedArray{MatPartitioning},Ptr{Void}))
    ccall((:MatPartitioningDestroy,petsc1),PetscErrorCode,(Ptr{MatPartitioning},),arg1)
end 
=#
function MatPartitioningRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:MatPartitioningRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function MatPartitioningView(arg1::MatPartitioning,arg2::PetscViewer{Float64})
    ccall((:MatPartitioningView,petsc1),PetscErrorCode,(MatPartitioning,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningSetFromOptions(arg0::Type{Float64},arg1::MatPartitioning)
    ccall((:MatPartitioningSetFromOptions,petsc1),PetscErrorCode,(MatPartitioning,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningGetType(arg0::Type{Float64},arg1::MatPartitioning,arg2::Union(Ptr{MatPartitioningType},StridedArray{MatPartitioningType},Ptr{Void}))
    ccall((:MatPartitioningGetType,petsc1),PetscErrorCode,(MatPartitioning,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningParmetisSetCoarseSequential(arg0::Type{Float64},arg1::MatPartitioning)
    ccall((:MatPartitioningParmetisSetCoarseSequential,petsc1),PetscErrorCode,(MatPartitioning,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningParmetisGetEdgeCut(arg0::Type{Float64},arg1::MatPartitioning,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:MatPartitioningParmetisGetEdgeCut,petsc1),PetscErrorCode,(MatPartitioning,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningChacoSetGlobal(arg0::Type{Float64},arg1::MatPartitioning,arg2::MPChacoGlobalType)
    ccall((:MatPartitioningChacoSetGlobal,petsc1),PetscErrorCode,(MatPartitioning,MPChacoGlobalType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningChacoGetGlobal(arg0::Type{Float64},arg1::MatPartitioning,arg2::Union(Ptr{MPChacoGlobalType},StridedArray{MPChacoGlobalType},Ptr{Void}))
    ccall((:MatPartitioningChacoGetGlobal,petsc1),PetscErrorCode,(MatPartitioning,Ptr{MPChacoGlobalType}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningChacoSetLocal(arg0::Type{Float64},arg1::MatPartitioning,arg2::MPChacoLocalType)
    ccall((:MatPartitioningChacoSetLocal,petsc1),PetscErrorCode,(MatPartitioning,MPChacoLocalType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningChacoGetLocal(arg0::Type{Float64},arg1::MatPartitioning,arg2::Union(Ptr{MPChacoLocalType},StridedArray{MPChacoLocalType},Ptr{Void}))
    ccall((:MatPartitioningChacoGetLocal,petsc1),PetscErrorCode,(MatPartitioning,Ptr{MPChacoLocalType}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningChacoSetCoarseLevel(arg0::Type{Float64},arg1::MatPartitioning,PetscReal::Integer)
    ccall((:MatPartitioningChacoSetCoarseLevel,petsc1),PetscErrorCode,(MatPartitioning,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningChacoSetEigenSolver(arg0::Type{Float64},arg1::MatPartitioning,arg2::MPChacoEigenType)
    ccall((:MatPartitioningChacoSetEigenSolver,petsc1),PetscErrorCode,(MatPartitioning,MPChacoEigenType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningChacoGetEigenSolver(arg0::Type{Float64},arg1::MatPartitioning,arg2::Union(Ptr{MPChacoEigenType},StridedArray{MPChacoEigenType},Ptr{Void}))
    ccall((:MatPartitioningChacoGetEigenSolver,petsc1),PetscErrorCode,(MatPartitioning,Ptr{MPChacoEigenType}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningChacoSetEigenTol(arg0::Type{Float64},arg1::MatPartitioning,PetscReal::Integer)
    ccall((:MatPartitioningChacoSetEigenTol,petsc1),PetscErrorCode,(MatPartitioning,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningChacoGetEigenTol(arg0::Type{Float64},arg1::MatPartitioning,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:MatPartitioningChacoGetEigenTol,petsc1),PetscErrorCode,(MatPartitioning,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningChacoSetEigenNumber(arg0::Type{Float64},arg1::MatPartitioning,arg2::Integer)
    ccall((:MatPartitioningChacoSetEigenNumber,petsc1),PetscErrorCode,(MatPartitioning,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningChacoGetEigenNumber(arg0::Type{Float64},arg1::MatPartitioning,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:MatPartitioningChacoGetEigenNumber,petsc1),PetscErrorCode,(MatPartitioning,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningPartySetGlobal(arg0::Type{Float64},arg1::MatPartitioning,arg2::Union(ByteString,Symbol))
    ccall((:MatPartitioningPartySetGlobal,petsc1),PetscErrorCode,(MatPartitioning,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningPartySetLocal(arg0::Type{Float64},arg1::MatPartitioning,arg2::Union(ByteString,Symbol))
    ccall((:MatPartitioningPartySetLocal,petsc1),PetscErrorCode,(MatPartitioning,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningPartySetCoarseLevel(arg0::Type{Float64},arg1::MatPartitioning,PetscReal::Integer)
    ccall((:MatPartitioningPartySetCoarseLevel,petsc1),PetscErrorCode,(MatPartitioning,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningPartySetBipart(arg0::Type{Float64},arg1::MatPartitioning,arg2::PetscBool)
    ccall((:MatPartitioningPartySetBipart,petsc1),PetscErrorCode,(MatPartitioning,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningPartySetMatchOptimization(arg0::Type{Float64},arg1::MatPartitioning,arg2::PetscBool)
    ccall((:MatPartitioningPartySetMatchOptimization,petsc1),PetscErrorCode,(MatPartitioning,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningPTScotchSetImbalance(arg0::Type{Float64},arg1::MatPartitioning,PetscReal::Integer)
    ccall((:MatPartitioningPTScotchSetImbalance,petsc1),PetscErrorCode,(MatPartitioning,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningPTScotchGetImbalance(arg0::Type{Float64},arg1::MatPartitioning,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:MatPartitioningPTScotchGetImbalance,petsc1),PetscErrorCode,(MatPartitioning,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningPTScotchSetStrategy(arg0::Type{Float64},arg1::MatPartitioning,arg2::MPPTScotchStrategyType)
    ccall((:MatPartitioningPTScotchSetStrategy,petsc1),PetscErrorCode,(MatPartitioning,MPPTScotchStrategyType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningPTScotchGetStrategy(arg0::Type{Float64},arg1::MatPartitioning,arg2::Union(Ptr{MPPTScotchStrategyType},StridedArray{MPPTScotchStrategyType},Ptr{Void}))
    ccall((:MatPartitioningPTScotchGetStrategy,petsc1),PetscErrorCode,(MatPartitioning,Ptr{MPPTScotchStrategyType}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatCoarsenCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{MatCoarsen},StridedArray{MatCoarsen},Ptr{Void}))
    ccall((:MatCoarsenCreate,petsc1),PetscErrorCode,(comm_type,Ptr{MatCoarsen}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatCoarsenSetType(arg0::Type{Float64},arg1::MatCoarsen,arg2::MatCoarsenType)
    ccall((:MatCoarsenSetType,petsc1),PetscErrorCode,(MatCoarsen,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatCoarsenSetAdjacency(arg1::MatCoarsen,arg2::Mat{Float64})
    ccall((:MatCoarsenSetAdjacency,petsc1),PetscErrorCode,(MatCoarsen,Mat{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatCoarsenSetGreedyOrdering(arg1::MatCoarsen,arg2::IS{Float64})
    ccall((:MatCoarsenSetGreedyOrdering,petsc1),PetscErrorCode,(MatCoarsen,IS{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatCoarsenSetStrictAggs(arg0::Type{Float64},arg1::MatCoarsen,arg2::PetscBool)
    ccall((:MatCoarsenSetStrictAggs,petsc1),PetscErrorCode,(MatCoarsen,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatCoarsenGetData(arg0::Type{Float64},arg1::MatCoarsen,arg2::Union(Ptr{Ptr{PetscCoarsenData}},StridedArray{Ptr{PetscCoarsenData}},Ptr{Void}))
    ccall((:MatCoarsenGetData,petsc1),PetscErrorCode,(MatCoarsen,Ptr{Ptr{PetscCoarsenData}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatCoarsenApply(arg0::Type{Float64},arg1::MatCoarsen)
    ccall((:MatCoarsenApply,petsc1),PetscErrorCode,(MatCoarsen,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function MatCoarsenDestroy(arg0::Type{Float64},arg1::Union(Ptr{MatCoarsen},StridedArray{MatCoarsen},Ptr{Void}))
    ccall((:MatCoarsenDestroy,petsc1),PetscErrorCode,(Ptr{MatCoarsen},),arg1)
end 
=#
function MatCoarsenRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:MatCoarsenRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function MatCoarsenView(arg1::MatCoarsen,arg2::PetscViewer{Float64})
    ccall((:MatCoarsenView,petsc1),PetscErrorCode,(MatCoarsen,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatCoarsenSetFromOptions(arg0::Type{Float64},arg1::MatCoarsen)
    ccall((:MatCoarsenSetFromOptions,petsc1),PetscErrorCode,(MatCoarsen,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function MatCoarsenGetType(arg0::Type{Float64},arg1::MatCoarsen,arg2::Union(Ptr{MatCoarsenType},StridedArray{MatCoarsenType},Ptr{Void}))
    ccall((:MatCoarsenGetType,petsc1),PetscErrorCode,(MatCoarsen,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
function MatMeshToCellGraph(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatMeshToCellGraph,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{Mat{Float64}}),arg1,arg2,arg3)
    return err
end

function MatHasOperation(arg1::Mat{Float64},arg2::MatOperation,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:MatHasOperation,petsc1),PetscErrorCode,(Mat{Float64},MatOperation,Ptr{PetscBool}),arg1,arg2,arg3)
    return err
end

function MatShellSetOperation(arg1::Mat{Float64},arg2::MatOperation,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:MatShellSetOperation,petsc1),PetscErrorCode,(Mat{Float64},MatOperation,Ptr{Void}),arg1,arg2,arg3)
    return err
end

function MatShellGetOperation(arg1::Mat{Float64},arg2::MatOperation,arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    err = ccall((:MatShellGetOperation,petsc1),PetscErrorCode,(Mat{Float64},MatOperation,Ptr{Ptr{Void}}),arg1,arg2,arg3)
    return err
end

function MatShellSetContext(arg1::Mat{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:MatShellSetContext,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Void}),arg1,arg2)
    return err
end

function MatMPIBAIJSetHashTableFactor(arg1::Mat{Float64},PetscReal::Integer)
    err = ccall((:MatMPIBAIJSetHashTableFactor,petsc1),PetscErrorCode,(Mat{Float64},Cint),arg1,PetscReal)
    return err
end

function MatISGetLocalMat(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatISGetLocalMat,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2)
    return err
end

function MatISSetLocalMat(arg1::Mat{Float64},arg2::Mat{Float64})
    err = ccall((:MatISSetLocalMat,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64}),arg1,arg2)
    return err
end

function MatISGetMPIXAIJ(arg1::Mat{Float64},arg2::MatReuse,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatISGetMPIXAIJ,petsc1),PetscErrorCode,(Mat{Float64},MatReuse,Ptr{Mat{Float64}}),arg1,arg2,arg3)
    return err
end

#= skipping function with undefined symbols: 
 function MatNullSpaceCreate(arg1::MPI_Comm,arg2::PetscBool,arg3::Integer,arg4::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg5::Union(Ptr{MatNullSpace},StridedArray{MatNullSpace},Ptr{Void}))
    ccall((:MatNullSpaceCreate,petsc1),PetscErrorCode,(comm_type,PetscBool,PetscInt,Ptr{Vec{Float64}},Ptr{MatNullSpace}),arg1.val,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function MatNullSpaceSetFunction(arg0::Type{Float64},arg1::MatNullSpace,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:MatNullSpaceSetFunction,petsc1),PetscErrorCode,(MatNullSpace,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function MatNullSpaceDestroy(arg0::Type{Float64},arg1::Union(Ptr{MatNullSpace},StridedArray{MatNullSpace},Ptr{Void}))
    ccall((:MatNullSpaceDestroy,petsc1),PetscErrorCode,(Ptr{MatNullSpace},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function MatNullSpaceRemove(arg1::MatNullSpace,arg2::Vec{Float64})
    ccall((:MatNullSpaceRemove,petsc1),PetscErrorCode,(MatNullSpace,Vec{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatGetNullSpace(arg1::Mat{Float64},arg2::Union(Ptr{MatNullSpace},StridedArray{MatNullSpace},Ptr{Void}))
    ccall((:MatGetNullSpace,petsc1),PetscErrorCode,(Mat{Float64},Ptr{MatNullSpace}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatGetTransposeNullSpace(arg1::Mat{Float64},arg2::Union(Ptr{MatNullSpace},StridedArray{MatNullSpace},Ptr{Void}))
    ccall((:MatGetTransposeNullSpace,petsc1),PetscErrorCode,(Mat{Float64},Ptr{MatNullSpace}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatSetTransposeNullSpace(arg1::Mat{Float64},arg2::MatNullSpace)
    ccall((:MatSetTransposeNullSpace,petsc1),PetscErrorCode,(Mat{Float64},MatNullSpace),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatSetNullSpace(arg1::Mat{Float64},arg2::MatNullSpace)
    ccall((:MatSetNullSpace,petsc1),PetscErrorCode,(Mat{Float64},MatNullSpace),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatSetNearNullSpace(arg1::Mat{Float64},arg2::MatNullSpace)
    ccall((:MatSetNearNullSpace,petsc1),PetscErrorCode,(Mat{Float64},MatNullSpace),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatGetNearNullSpace(arg1::Mat{Float64},arg2::Union(Ptr{MatNullSpace},StridedArray{MatNullSpace},Ptr{Void}))
    ccall((:MatGetNearNullSpace,petsc1),PetscErrorCode,(Mat{Float64},Ptr{MatNullSpace}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatNullSpaceTest(arg1::MatNullSpace,arg2::Mat{Float64},arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:MatNullSpaceTest,petsc1),PetscErrorCode,(MatNullSpace,Mat{Float64},Ptr{PetscBool}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function MatNullSpaceView(arg1::MatNullSpace,arg2::PetscViewer{Float64})
    ccall((:MatNullSpaceView,petsc1),PetscErrorCode,(MatNullSpace,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatNullSpaceGetVecs(arg1::MatNullSpace,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Ptr{Vec{Float64}}},StridedArray{Ptr{Vec{Float64}}},Ptr{Void}))
    ccall((:MatNullSpaceGetVecs,petsc1),PetscErrorCode,(MatNullSpace,Ptr{PetscBool},Ptr{PetscInt},Ptr{Ptr{Vec{Float64}}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function MatNullSpaceCreateRigidBody(arg1::Vec{Float64},arg2::Union(Ptr{MatNullSpace},StridedArray{MatNullSpace},Ptr{Void}))
    ccall((:MatNullSpaceCreateRigidBody,petsc1),PetscErrorCode,(Vec{Float64},Ptr{MatNullSpace}),arg1,arg2)
end 
=#
function MatReorderingSeqSBAIJ(arg1::Mat{Float64},arg2::IS{Float64})
    err = ccall((:MatReorderingSeqSBAIJ,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64}),arg1,arg2)
    return err
end

function MatMPISBAIJSetHashTableFactor(arg1::Mat{Float64},PetscReal::Integer)
    err = ccall((:MatMPISBAIJSetHashTableFactor,petsc1),PetscErrorCode,(Mat{Float64},Cint),arg1,PetscReal)
    return err
end

function MatSeqSBAIJSetColumnIndices(arg1::Mat{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatSeqSBAIJSetColumnIndices,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function MatSeqBAIJInvertBlockDiagonal(arg1::Mat{Float64})
    err = ccall((:MatSeqBAIJInvertBlockDiagonal,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
    return err
end

function MatCreateMAIJ(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateMAIJ,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{Mat{Float64}}),arg1,arg2,arg3)
    return err
end

function MatMAIJRedimension(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatMAIJRedimension,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{Mat{Float64}}),arg1,arg2,arg3)
    return err
end

function MatMAIJGetAIJ(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatMAIJGetAIJ,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2)
    return err
end

function MatComputeExplicitOperator(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatComputeExplicitOperator,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2)
    return err
end

function MatDiagonalScaleLocal(arg1::Mat{Float64},arg2::Vec{Float64})
    err = ccall((:MatDiagonalScaleLocal,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64}),arg1,arg2)
    return err
end

function MatCreateMFFD(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateMFFD,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatMFFDSetBase(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:MatMFFDSetBase,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function MatMFFDSetFunction(arg1::Mat{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:MatMFFDSetFunction,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
    return err
end

function MatMFFDSetFunctioni(arg1::Mat{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:MatMFFDSetFunctioni,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Void}),arg1,arg2)
    return err
end

function MatMFFDSetFunctioniBase(arg1::Mat{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:MatMFFDSetFunctioniBase,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Void}),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function MatMFFDAddNullSpace(arg1::Mat{Float64},arg2::MatNullSpace)
    ccall((:MatMFFDAddNullSpace,petsc1),PetscErrorCode,(Mat{Float64},MatNullSpace),arg1,arg2)
end 
=#
function MatMFFDSetHHistory(arg1::Mat{Float64},arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg3::Integer)
    err = ccall((:MatMFFDSetHHistory,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Float64},PetscInt),arg1,arg2,arg3)
    return err
end

function MatMFFDResetHHistory(arg1::Mat{Float64})
    err = ccall((:MatMFFDResetHHistory,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
    return err
end

function MatMFFDSetFunctionError(arg1::Mat{Float64},PetscReal::Integer)
    err = ccall((:MatMFFDSetFunctionError,petsc1),PetscErrorCode,(Mat{Float64},Cint),arg1,PetscReal)
    return err
end

function MatMFFDSetPeriod(arg1::Mat{Float64},arg2::Integer)
    err = ccall((:MatMFFDSetPeriod,petsc1),PetscErrorCode,(Mat{Float64},PetscInt),arg1,arg2)
    return err
end

function MatMFFDGetH(arg1::Mat{Float64},arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:MatMFFDGetH,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Float64}),arg1,arg2)
    return err
end

function MatMFFDSetOptionsPrefix(arg1::Mat{Float64},arg2::Union(ByteString,Symbol))
    err = ccall((:MatMFFDSetOptionsPrefix,petsc1),PetscErrorCode,(Mat{Float64},Cstring),arg1,arg2)
    return err
end

function MatMFFDCheckPositivity(arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:MatMFFDCheckPositivity,petsc1),PetscErrorCode,(Ptr{Void},Vec{Float64},Vec{Float64},Ptr{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function MatMFFDSetCheckh(arg1::Mat{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:MatMFFDSetCheckh,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
    return err
end

function MatMFFDSetType(arg1::Mat{Float64},arg2::MatMFFDType)
    err = ccall((:MatMFFDSetType,petsc1),PetscErrorCode,(Mat{Float64},Cstring),arg1,arg2)
    return err
end

function MatMFFDRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:MatMFFDRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

function MatMFFDDSSetUmin(arg1::Mat{Float64},PetscReal::Integer)
    err = ccall((:MatMFFDDSSetUmin,petsc1),PetscErrorCode,(Mat{Float64},Cint),arg1,PetscReal)
    return err
end

function MatMFFDWPSetComputeNormU(arg1::Mat{Float64},arg2::PetscBool)
    err = ccall((:MatMFFDWPSetComputeNormU,petsc1),PetscErrorCode,(Mat{Float64},PetscBool),arg1,arg2)
    return err
end

function PetscViewerMathematicaPutMatrix(arg1::PetscViewer{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:PetscViewerMathematicaPutMatrix,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscInt,PetscInt,Ptr{Cint}),arg1,arg2,arg3,arg4)
    return err
end

function PetscViewerMathematicaPutCSRMatrix(arg1::PetscViewer{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:PetscViewerMathematicaPutCSRMatrix,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatCreateNest(arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg6::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateNest,petsc1),PetscErrorCode,(comm_type,PetscInt,Ptr{IS{Float64}},PetscInt,Ptr{IS{Float64}},Ptr{Mat{Float64}},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
    return err
end

function MatNestGetSize(arg1::Mat{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatNestGetSize,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function MatNestGetISs(arg1::Mat{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:MatNestGetISs,petsc1),PetscErrorCode,(Mat{Float64},Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3)
    return err
end

function MatNestGetLocalISs(arg1::Mat{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:MatNestGetLocalISs,petsc1),PetscErrorCode,(Mat{Float64},Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3)
    return err
end

function MatNestGetSubMats(arg1::Mat{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Ptr{Ptr{Mat{Float64}}}},StridedArray{Ptr{Ptr{Mat{Float64}}}},Ptr{Void}))
    err = ccall((:MatNestGetSubMats,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{Ptr{Mat{Float64}}}}),arg1,arg2,arg3,arg4)
    return err
end

function MatNestGetSubMat(arg1::Mat{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatNestGetSubMat,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,PetscInt,Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function MatNestSetVecType(arg1::Mat{Float64},arg2::VecType)
    err = ccall((:MatNestSetVecType,petsc1),PetscErrorCode,(Mat{Float64},Cstring),arg1,arg2)
    return err
end

function MatNestSetSubMats(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg6::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatNestSetSubMats,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{IS{Float64}},PetscInt,Ptr{IS{Float64}},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatNestSetSubMat(arg1::Mat{Float64},arg2::Integer,arg3::Integer,arg4::Mat{Float64})
    err = ccall((:MatNestSetSubMat,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,PetscInt,Mat{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function MatChop(arg1::Mat{Float64},PetscReal::Integer)
    err = ccall((:MatChop,petsc1),PetscErrorCode,(Mat{Float64},Cint),arg1,PetscReal)
    return err
end

function MatComputeBandwidth(arg1::Mat{Float64},PetscReal::Integer,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:MatComputeBandwidth,petsc1),PetscErrorCode,(Mat{Float64},Cint,Ptr{PetscInt}),arg1,PetscReal,arg2)
    return err
end

function MatSubdomainsCreateCoalesce(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    err = ccall((:MatSubdomainsCreateCoalesce,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},Ptr{Ptr{IS{Float64}}}),arg1,arg2,arg3,arg4)
    return err
end

function PCExoticSetType(arg1::PC{Float64},arg2::PCExoticType)
    err = ccall((:PCExoticSetType,petsc1),PetscErrorCode,(PC{Float64},PCExoticType),arg1,arg2)
    return err
end

function PCInitializePackage(arg0::Type{Float64})
    err = ccall((:PCInitializePackage,petsc1),PetscErrorCode,())
    return err
end

function PCCreate(arg1::MPI_Comm,arg2::Union(Ptr{PC{Float64}},StridedArray{PC{Float64}},Ptr{Void}))
    err = ccall((:PCCreate,petsc1),PetscErrorCode,(comm_type,Ptr{PC{Float64}}),arg1.val,arg2)
    return err
end

function PCSetType(arg1::PC{Float64},arg2::PCType)
    err = ccall((:PCSetType,petsc1),PetscErrorCode,(PC{Float64},Cstring),arg1,arg2)
    return err
end

function PCGetType(arg1::PC{Float64},arg2::Union(Ptr{PCType},StridedArray{PCType},Ptr{Void}))
    (arg2_,tmp) = symbol_get_before(arg2)
    err = ccall((:PCGetType,petsc1),PetscErrorCode,(PC{Float64},Ptr{Ptr{Uint8}}),arg1,arg2_)
    symbol_get_after(arg2_,arg2)
    return err
end

function PCSetUp(arg1::PC{Float64})
    err = ccall((:PCSetUp,petsc1),PetscErrorCode,(PC{Float64},),arg1)
    return err
end

function PCGetSetUpFailedReason(arg1::PC{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PCGetSetUpFailedReason,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function PCSetUpOnBlocks(arg1::PC{Float64})
    err = ccall((:PCSetUpOnBlocks,petsc1),PetscErrorCode,(PC{Float64},),arg1)
    return err
end

function PCApply(arg1::PC{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:PCApply,petsc1),PetscErrorCode,(PC{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function PCApplySymmetricLeft(arg1::PC{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:PCApplySymmetricLeft,petsc1),PetscErrorCode,(PC{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function PCApplySymmetricRight(arg1::PC{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:PCApplySymmetricRight,petsc1),PetscErrorCode,(PC{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function PCApplyBAorAB(arg1::PC{Float64},arg2::PCSide,arg3::Vec{Float64},arg4::Vec{Float64},arg5::Vec{Float64})
    err = ccall((:PCApplyBAorAB,petsc1),PetscErrorCode,(PC{Float64},PCSide,Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function PCApplyTranspose(arg1::PC{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:PCApplyTranspose,petsc1),PetscErrorCode,(PC{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function PCApplyTransposeExists(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PCApplyTransposeExists,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function PCApplyBAorABTranspose(arg1::PC{Float64},arg2::PCSide,arg3::Vec{Float64},arg4::Vec{Float64},arg5::Vec{Float64})
    err = ccall((:PCApplyBAorABTranspose,petsc1),PetscErrorCode,(PC{Float64},PCSide,Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function PCSetReusePreconditioner(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCSetReusePreconditioner,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCGetReusePreconditioner(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PCGetReusePreconditioner,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function PCSetErrorIfFailure(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCSetErrorIfFailure,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCApplyRichardson(arg1::PC{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64},PetscReal::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::PetscBool,arg9::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg10::Union(Ptr{PCRichardsonConvergedReason},StridedArray{PCRichardsonConvergedReason},Ptr{Void}))
    err = ccall((:PCApplyRichardson,petsc1),PetscErrorCode,(PC{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Cint,Cint,Cint,PetscInt,PetscBool,Ptr{PetscInt},Ptr{PCRichardsonConvergedReason}),arg1,arg2,arg3,arg4,PetscReal,arg5,arg6,arg7,arg8,arg9,arg10)
    return err
end

function PCApplyRichardsonExists(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PCApplyRichardsonExists,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function PCSetInitialGuessNonzero(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCSetInitialGuessNonzero,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCGetInitialGuessNonzero(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PCGetInitialGuessNonzero,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function PCSetUseAmat(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCSetUseAmat,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCGetUseAmat(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PCGetUseAmat,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function PCRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PCRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

function PCReset(arg1::PC{Float64})
    err = ccall((:PCReset,petsc1),PetscErrorCode,(PC{Float64},),arg1)
    return err
end

function PCDestroy(arg1::Union(Ptr{PC{Float64}},StridedArray{PC{Float64}},Ptr{Void}))
    err = ccall((:PCDestroy,petsc1),PetscErrorCode,(Ptr{PC{Float64}},),arg1)
    return err
end

function PCSetFromOptions(arg1::PC{Float64})
    err = ccall((:PCSetFromOptions,petsc1),PetscErrorCode,(PC{Float64},),arg1)
    return err
end

function PCFactorGetMatrix(arg1::PC{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:PCFactorGetMatrix,petsc1),PetscErrorCode,(PC{Float64},Ptr{Mat{Float64}}),arg1,arg2)
    return err
end

function PCSetModifySubMatrices(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PCSetModifySubMatrices,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
    return err
end

function PCModifySubMatrices(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg5::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg6::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PCModifySubMatrices,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Ptr{IS{Float64}},Ptr{IS{Float64}},Ptr{Mat{Float64}},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function PCSetOperators(arg1::PC{Float64},arg2::Mat{Float64},arg3::Mat{Float64})
    err = ccall((:PCSetOperators,petsc1),PetscErrorCode,(PC{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3)
    return err
end

function PCGetOperators(arg1::PC{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:PCGetOperators,petsc1),PetscErrorCode,(PC{Float64},Ptr{Mat{Float64}},Ptr{Mat{Float64}}),arg1,arg2,arg3)
    return err
end

function PCGetOperatorsSet(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PCGetOperatorsSet,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3)
    return err
end

function PCView(arg1::PC{Float64},arg2::PetscViewer{Float64})
    err = ccall((:PCView,petsc1),PetscErrorCode,(PC{Float64},PetscViewer{Float64}),arg1,arg2)
    return err
end

function PCLoad(arg1::PC{Float64},arg2::PetscViewer{Float64})
    err = ccall((:PCLoad,petsc1),PetscErrorCode,(PC{Float64},PetscViewer{Float64}),arg1,arg2)
    return err
end

function PCAppendOptionsPrefix(arg1::PC{Float64},arg2::Union(ByteString,Symbol))
    err = ccall((:PCAppendOptionsPrefix,petsc1),PetscErrorCode,(PC{Float64},Cstring),arg1,arg2)
    return err
end

function PCGetOptionsPrefix(arg1::PC{Float64},arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:PCGetOptionsPrefix,petsc1),PetscErrorCode,(PC{Float64},Ptr{Ptr{Uint8}}),arg1,arg2)
    return err
end

function PCComputeExplicitOperator(arg1::PC{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:PCComputeExplicitOperator,petsc1),PetscErrorCode,(PC{Float64},Ptr{Mat{Float64}}),arg1,arg2)
    return err
end

function PCGetDiagonalScale(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PCGetDiagonalScale,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function PCDiagonalScaleLeft(arg1::PC{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:PCDiagonalScaleLeft,petsc1),PetscErrorCode,(PC{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function PCDiagonalScaleRight(arg1::PC{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:PCDiagonalScaleRight,petsc1),PetscErrorCode,(PC{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function PCSetDiagonalScale(arg1::PC{Float64},arg2::Vec{Float64})
    err = ccall((:PCSetDiagonalScale,petsc1),PetscErrorCode,(PC{Float64},Vec{Float64}),arg1,arg2)
    return err
end

function PCJacobiSetType(arg1::PC{Float64},arg2::PCJacobiType)
    err = ccall((:PCJacobiSetType,petsc1),PetscErrorCode,(PC{Float64},PCJacobiType),arg1,arg2)
    return err
end

function PCJacobiGetType(arg1::PC{Float64},arg2::Union(Ptr{PCJacobiType},StridedArray{PCJacobiType},Ptr{Void}))
    err = ccall((:PCJacobiGetType,petsc1),PetscErrorCode,(PC{Float64},Ptr{PCJacobiType}),arg1,arg2)
    return err
end

function PCJacobiSetUseAbs(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCJacobiSetUseAbs,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCJacobiGetUseAbs(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PCJacobiGetUseAbs,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function PCSORSetSymmetric(arg1::PC{Float64},arg2::MatSORType)
    err = ccall((:PCSORSetSymmetric,petsc1),PetscErrorCode,(PC{Float64},MatSORType),arg1,arg2)
    return err
end

function PCSORGetSymmetric(arg1::PC{Float64},arg2::Union(Ptr{MatSORType},StridedArray{MatSORType},Ptr{Void}))
    err = ccall((:PCSORGetSymmetric,petsc1),PetscErrorCode,(PC{Float64},Ptr{MatSORType}),arg1,arg2)
    return err
end

function PCSORSetOmega(arg1::PC{Float64},PetscReal::Integer)
    err = ccall((:PCSORSetOmega,petsc1),PetscErrorCode,(PC{Float64},Cint),arg1,PetscReal)
    return err
end

function PCSORGetOmega(arg1::PC{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:PCSORGetOmega,petsc1),PetscErrorCode,(PC{Float64},Ptr{Cint}),arg1,arg2)
    return err
end

function PCSORSetIterations(arg1::PC{Float64},arg2::Integer,arg3::Integer)
    err = ccall((:PCSORSetIterations,petsc1),PetscErrorCode,(PC{Float64},PetscInt,PetscInt),arg1,arg2,arg3)
    return err
end

function PCSORGetIterations(arg1::PC{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PCSORGetIterations,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function PCEisenstatSetOmega(arg1::PC{Float64},PetscReal::Integer)
    err = ccall((:PCEisenstatSetOmega,petsc1),PetscErrorCode,(PC{Float64},Cint),arg1,PetscReal)
    return err
end

function PCEisenstatGetOmega(arg1::PC{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:PCEisenstatGetOmega,petsc1),PetscErrorCode,(PC{Float64},Ptr{Cint}),arg1,arg2)
    return err
end

function PCEisenstatSetNoDiagonalScaling(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCEisenstatSetNoDiagonalScaling,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCEisenstatGetNoDiagonalScaling(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PCEisenstatGetNoDiagonalScaling,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function PCBJacobiSetTotalBlocks(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PCBJacobiSetTotalBlocks,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function PCBJacobiGetTotalBlocks(arg1::PC{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:PCBJacobiGetTotalBlocks,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3)
    return err
end

function PCBJacobiSetLocalBlocks(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PCBJacobiSetLocalBlocks,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function PCBJacobiGetLocalBlocks(arg1::PC{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:PCBJacobiGetLocalBlocks,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3)
    return err
end

function PCShellSetApply(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PCShellSetApply,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
    return err
end

function PCShellSetApplyBA(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PCShellSetApplyBA,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
    return err
end

function PCShellSetApplyTranspose(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PCShellSetApplyTranspose,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
    return err
end

function PCShellSetSetUp(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PCShellSetSetUp,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
    return err
end

function PCShellSetApplyRichardson(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PCShellSetApplyRichardson,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
    return err
end

function PCShellSetView(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PCShellSetView,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
    return err
end

function PCShellSetDestroy(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PCShellSetDestroy,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
    return err
end

function PCShellSetContext(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PCShellSetContext,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
    return err
end

function PCShellGetContext(arg1::PC{Float64},arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    err = ccall((:PCShellGetContext,petsc1),PetscErrorCode,(PC{Float64},Ptr{Ptr{Void}}),arg1,arg2)
    return err
end

function PCShellSetName(arg1::PC{Float64},arg2::Union(ByteString,Symbol))
    err = ccall((:PCShellSetName,petsc1),PetscErrorCode,(PC{Float64},Cstring),arg1,arg2)
    return err
end

function PCShellGetName(arg1::PC{Float64},arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:PCShellGetName,petsc1),PetscErrorCode,(PC{Float64},Ptr{Ptr{Uint8}}),arg1,arg2)
    return err
end

function PCFactorSetZeroPivot(arg1::PC{Float64},PetscReal::Integer)
    err = ccall((:PCFactorSetZeroPivot,petsc1),PetscErrorCode,(PC{Float64},Cint),arg1,PetscReal)
    return err
end

function PCFactorSetShiftType(arg1::PC{Float64},arg2::MatFactorShiftType)
    err = ccall((:PCFactorSetShiftType,petsc1),PetscErrorCode,(PC{Float64},MatFactorShiftType),arg1,arg2)
    return err
end

function PCFactorSetShiftAmount(arg1::PC{Float64},PetscReal::Integer)
    err = ccall((:PCFactorSetShiftAmount,petsc1),PetscErrorCode,(PC{Float64},Cint),arg1,PetscReal)
    return err
end

function PCFactorSetMatSolverPackage(arg1::PC{Float64},arg2::Union(ByteString,Symbol))
    err = ccall((:PCFactorSetMatSolverPackage,petsc1),PetscErrorCode,(PC{Float64},Cstring),arg1,arg2)
    return err
end

function PCFactorGetMatSolverPackage(arg1::PC{Float64},arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:PCFactorGetMatSolverPackage,petsc1),PetscErrorCode,(PC{Float64},Ptr{Ptr{Uint8}}),arg1,arg2)
    return err
end

function PCFactorSetUpMatSolverPackage(arg1::PC{Float64})
    err = ccall((:PCFactorSetUpMatSolverPackage,petsc1),PetscErrorCode,(PC{Float64},),arg1)
    return err
end

function PCFactorSetFill(arg1::PC{Float64},PetscReal::Integer)
    err = ccall((:PCFactorSetFill,petsc1),PetscErrorCode,(PC{Float64},Cint),arg1,PetscReal)
    return err
end

function PCFactorSetColumnPivot(arg1::PC{Float64},PetscReal::Integer)
    err = ccall((:PCFactorSetColumnPivot,petsc1),PetscErrorCode,(PC{Float64},Cint),arg1,PetscReal)
    return err
end

function PCFactorReorderForNonzeroDiagonal(arg1::PC{Float64},PetscReal::Integer)
    err = ccall((:PCFactorReorderForNonzeroDiagonal,petsc1),PetscErrorCode,(PC{Float64},Cint),arg1,PetscReal)
    return err
end

function PCFactorSetMatOrderingType(arg1::PC{Float64},arg2::MatOrderingType)
    err = ccall((:PCFactorSetMatOrderingType,petsc1),PetscErrorCode,(PC{Float64},Cstring),arg1,arg2)
    return err
end

function PCFactorSetReuseOrdering(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCFactorSetReuseOrdering,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCFactorSetReuseFill(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCFactorSetReuseFill,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCFactorSetUseInPlace(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCFactorSetUseInPlace,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCFactorGetUseInPlace(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PCFactorGetUseInPlace,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function PCFactorSetAllowDiagonalFill(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCFactorSetAllowDiagonalFill,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCFactorGetAllowDiagonalFill(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PCFactorGetAllowDiagonalFill,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function PCFactorSetPivotInBlocks(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCFactorSetPivotInBlocks,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCFactorSetLevels(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCFactorSetLevels,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCFactorGetLevels(arg1::PC{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PCFactorGetLevels,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function PCFactorSetDropTolerance(arg1::PC{Float64},PetscReal::Integer,arg2::Integer,arg3::Integer)
    err = ccall((:PCFactorSetDropTolerance,petsc1),PetscErrorCode,(PC{Float64},Cint,Cint,PetscInt),arg1,PetscReal,arg2,arg3)
    return err
end

function PCASMSetLocalSubdomains(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:PCASMSetLocalSubdomains,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function PCASMSetTotalSubdomains(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:PCASMSetTotalSubdomains,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function PCASMSetOverlap(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCASMSetOverlap,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCASMSetDMSubdomains(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCASMSetDMSubdomains,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCASMGetDMSubdomains(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PCASMGetDMSubdomains,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function PCASMSetSortIndices(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCASMSetSortIndices,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCASMSetType(arg1::PC{Float64},arg2::PCASMType)
    err = ccall((:PCASMSetType,petsc1),PetscErrorCode,(PC{Float64},PCASMType),arg1,arg2)
    return err
end

function PCASMGetType(arg1::PC{Float64},arg2::Union(Ptr{PCASMType},StridedArray{PCASMType},Ptr{Void}))
    err = ccall((:PCASMGetType,petsc1),PetscErrorCode,(PC{Float64},Ptr{PCASMType}),arg1,arg2)
    return err
end

function PCASMSetLocalType(arg1::PC{Float64},arg2::PCCompositeType)
    err = ccall((:PCASMSetLocalType,petsc1),PetscErrorCode,(PC{Float64},PCCompositeType),arg1,arg2)
    return err
end

function PCASMGetLocalType(arg1::PC{Float64},arg2::Union(Ptr{PCCompositeType},StridedArray{PCCompositeType},Ptr{Void}))
    err = ccall((:PCASMGetLocalType,petsc1),PetscErrorCode,(PC{Float64},Ptr{PCCompositeType}),arg1,arg2)
    return err
end

function PCASMCreateSubdomains(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    err = ccall((:PCASMCreateSubdomains,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{Ptr{IS{Float64}}}),arg1,arg2,arg3)
    return err
end

function PCASMDestroySubdomains(arg1::Integer,arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:PCASMDestroySubdomains,petsc1),PetscErrorCode,(PetscInt,Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3)
    return err
end

function PCASMCreateSubdomains2D(arg1::Integer,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg8::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}),arg9::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    err = ccall((:PCASMCreateSubdomains2D,petsc1),PetscErrorCode,(PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Ptr{IS{Float64}}},Ptr{Ptr{IS{Float64}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
    return err
end

function PCASMGetLocalSubdomains(arg1::PC{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}),arg4::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    err = ccall((:PCASMGetLocalSubdomains,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscInt},Ptr{Ptr{IS{Float64}}},Ptr{Ptr{IS{Float64}}}),arg1,arg2,arg3,arg4)
    return err
end

function PCASMGetLocalSubmatrices(arg1::PC{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{Mat{Float64}}},StridedArray{Ptr{Mat{Float64}}},Ptr{Void}))
    err = ccall((:PCASMGetLocalSubmatrices,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscInt},Ptr{Ptr{Mat{Float64}}}),arg1,arg2,arg3)
    return err
end

function PCGASMSetTotalSubdomains(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCGASMSetTotalSubdomains,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCGASMSetSubdomains(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:PCGASMSetSubdomains,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function PCGASMSetOverlap(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCGASMSetOverlap,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCGASMSetUseDMSubdomains(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCGASMSetUseDMSubdomains,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCGASMGetUseDMSubdomains(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PCGASMGetUseDMSubdomains,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function PCGASMSetSortIndices(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCGASMSetSortIndices,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCGASMSetType(arg1::PC{Float64},arg2::PCGASMType)
    err = ccall((:PCGASMSetType,petsc1),PetscErrorCode,(PC{Float64},PCGASMType),arg1,arg2)
    return err
end

function PCGASMCreateSubdomains(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    err = ccall((:PCGASMCreateSubdomains,petsc1),PetscErrorCode,(Mat{Float64},PetscInt,Ptr{PetscInt},Ptr{Ptr{IS{Float64}}}),arg1,arg2,arg3,arg4)
    return err
end

function PCGASMDestroySubdomains(arg1::Integer,arg2::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}),arg3::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    err = ccall((:PCGASMDestroySubdomains,petsc1),PetscErrorCode,(PetscInt,Ptr{Ptr{IS{Float64}}},Ptr{Ptr{IS{Float64}}}),arg1,arg2,arg3)
    return err
end

function PCGASMCreateSubdomains2D(arg1::PC{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg9::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}),arg10::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    err = ccall((:PCGASMCreateSubdomains2D,petsc1),PetscErrorCode,(PC{Float64},PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Ptr{IS{Float64}}},Ptr{Ptr{IS{Float64}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
    return err
end

function PCGASMGetSubdomains(arg1::PC{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}),arg4::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    err = ccall((:PCGASMGetSubdomains,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscInt},Ptr{Ptr{IS{Float64}}},Ptr{Ptr{IS{Float64}}}),arg1,arg2,arg3,arg4)
    return err
end

function PCGASMGetSubmatrices(arg1::PC{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{Mat{Float64}}},StridedArray{Ptr{Mat{Float64}}},Ptr{Void}))
    err = ccall((:PCGASMGetSubmatrices,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscInt},Ptr{Ptr{Mat{Float64}}}),arg1,arg2,arg3)
    return err
end

function PCCompositeSetType(arg1::PC{Float64},arg2::PCCompositeType)
    err = ccall((:PCCompositeSetType,petsc1),PetscErrorCode,(PC{Float64},PCCompositeType),arg1,arg2)
    return err
end

function PCCompositeGetType(arg1::PC{Float64},arg2::Union(Ptr{PCCompositeType},StridedArray{PCCompositeType},Ptr{Void}))
    err = ccall((:PCCompositeGetType,petsc1),PetscErrorCode,(PC{Float64},Ptr{PCCompositeType}),arg1,arg2)
    return err
end

function PCCompositeAddPC(arg1::PC{Float64},arg2::PCType)
    err = ccall((:PCCompositeAddPC,petsc1),PetscErrorCode,(PC{Float64},Cstring),arg1,arg2)
    return err
end

function PCCompositeGetPC(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{PC{Float64}},StridedArray{PC{Float64}},Ptr{Void}))
    err = ccall((:PCCompositeGetPC,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Ptr{PC{Float64}}),arg1,arg2,arg3)
    return err
end

function PCCompositeSpecialSetAlpha(arg1::PC{Float64},arg2::Float64)
    err = ccall((:PCCompositeSpecialSetAlpha,petsc1),PetscErrorCode,(PC{Float64},Float64),arg1,arg2)
    return err
end

function PCRedundantSetNumber(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCRedundantSetNumber,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCRedundantSetScatter(arg1::PC{Float64},arg2::VecScatter{Float64},arg3::VecScatter{Float64})
    err = ccall((:PCRedundantSetScatter,petsc1),PetscErrorCode,(PC{Float64},VecScatter{Float64},VecScatter{Float64}),arg1,arg2,arg3)
    return err
end

function PCRedundantGetOperators(arg1::PC{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:PCRedundantGetOperators,petsc1),PetscErrorCode,(PC{Float64},Ptr{Mat{Float64}},Ptr{Mat{Float64}}),arg1,arg2,arg3)
    return err
end

function PCSPAISetEpsilon(arg1::PC{Float64},arg2::Cdouble)
    err = ccall((:PCSPAISetEpsilon,petsc1),PetscErrorCode,(PC{Float64},Cdouble),arg1,arg2)
    return err
end

function PCSPAISetNBSteps(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCSPAISetNBSteps,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCSPAISetMax(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCSPAISetMax,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCSPAISetMaxNew(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCSPAISetMaxNew,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCSPAISetBlockSize(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCSPAISetBlockSize,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCSPAISetCacheSize(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCSPAISetCacheSize,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCSPAISetVerbose(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCSPAISetVerbose,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCSPAISetSp(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCSPAISetSp,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCHYPRESetType(arg1::PC{Float64},arg2::Union(ByteString,Symbol))
    err = ccall((:PCHYPRESetType,petsc1),PetscErrorCode,(PC{Float64},Cstring),arg1,arg2)
    return err
end

function PCHYPREGetType(arg1::PC{Float64},arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:PCHYPREGetType,petsc1),PetscErrorCode,(PC{Float64},Ptr{Ptr{Uint8}}),arg1,arg2)
    return err
end

function PCHYPRESetDiscreteGradient(arg1::PC{Float64},arg2::Mat{Float64})
    err = ccall((:PCHYPRESetDiscreteGradient,petsc1),PetscErrorCode,(PC{Float64},Mat{Float64}),arg1,arg2)
    return err
end

function PCHYPRESetDiscreteCurl(arg1::PC{Float64},arg2::Mat{Float64})
    err = ccall((:PCHYPRESetDiscreteCurl,petsc1),PetscErrorCode,(PC{Float64},Mat{Float64}),arg1,arg2)
    return err
end

function PCHYPRESetEdgeConstantVectors(arg1::PC{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64})
    err = ccall((:PCHYPRESetEdgeConstantVectors,petsc1),PetscErrorCode,(PC{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function PCHYPRESetAlphaPoissonMatrix(arg1::PC{Float64},arg2::Mat{Float64})
    err = ccall((:PCHYPRESetAlphaPoissonMatrix,petsc1),PetscErrorCode,(PC{Float64},Mat{Float64}),arg1,arg2)
    return err
end

function PCHYPRESetBetaPoissonMatrix(arg1::PC{Float64},arg2::Mat{Float64})
    err = ccall((:PCHYPRESetBetaPoissonMatrix,petsc1),PetscErrorCode,(PC{Float64},Mat{Float64}),arg1,arg2)
    return err
end

function PCBJacobiGetLocalBlocks(arg1::PC{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:PCBJacobiGetLocalBlocks,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3)
    return err
end

function PCBJacobiGetTotalBlocks(arg1::PC{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    err = ccall((:PCBJacobiGetTotalBlocks,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3)
    return err
end

function PCFieldSplitSetFields(arg1::PC{Float64},arg2::Union(ByteString,Symbol),arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PCFieldSplitSetFields,petsc1),PetscErrorCode,(PC{Float64},Cstring,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function PCFieldSplitSetType(arg1::PC{Float64},arg2::PCCompositeType)
    err = ccall((:PCFieldSplitSetType,petsc1),PetscErrorCode,(PC{Float64},PCCompositeType),arg1,arg2)
    return err
end

function PCFieldSplitGetType(arg1::PC{Float64},arg2::Union(Ptr{PCCompositeType},StridedArray{PCCompositeType},Ptr{Void}))
    err = ccall((:PCFieldSplitGetType,petsc1),PetscErrorCode,(PC{Float64},Ptr{PCCompositeType}),arg1,arg2)
    return err
end

function PCFieldSplitSetBlockSize(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCFieldSplitSetBlockSize,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCFieldSplitSetIS(arg1::PC{Float64},arg2::Union(ByteString,Symbol),arg3::IS{Float64})
    err = ccall((:PCFieldSplitSetIS,petsc1),PetscErrorCode,(PC{Float64},Cstring,IS{Float64}),arg1,arg2,arg3)
    return err
end

function PCFieldSplitGetIS(arg1::PC{Float64},arg2::Union(ByteString,Symbol),arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:PCFieldSplitGetIS,petsc1),PetscErrorCode,(PC{Float64},Cstring,Ptr{IS{Float64}}),arg1,arg2,arg3)
    return err
end

function PCFieldSplitSetDMSplits(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCFieldSplitSetDMSplits,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCFieldSplitGetDMSplits(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PCFieldSplitGetDMSplits,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function PCFieldSplitSetDiagUseAmat(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCFieldSplitSetDiagUseAmat,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCFieldSplitGetDiagUseAmat(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PCFieldSplitGetDiagUseAmat,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function PCFieldSplitSetOffDiagUseAmat(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCFieldSplitSetOffDiagUseAmat,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCFieldSplitGetOffDiagUseAmat(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PCFieldSplitGetOffDiagUseAmat,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function PETSC_DEPRECATED(arg0::Type{Float64})
    err = ccall((:PETSC_DEPRECATED,petsc1),Cint,())
    return err
end

function PCFieldSplitSchurPrecondition(arg1::PC{Float64},arg2::PCFieldSplitSchurPreType,arg3::Mat{Float64})
    err = ccall((:PCFieldSplitSchurPrecondition,petsc1),PetscErrorCode,(PC{Float64},PCFieldSplitSchurPreType,Mat{Float64}),arg1,arg2,arg3)
    return err
end

function PCFieldSplitSetSchurPre(arg1::PC{Float64},arg2::PCFieldSplitSchurPreType,arg3::Mat{Float64})
    err = ccall((:PCFieldSplitSetSchurPre,petsc1),PetscErrorCode,(PC{Float64},PCFieldSplitSchurPreType,Mat{Float64}),arg1,arg2,arg3)
    return err
end

function PCFieldSplitGetSchurPre(arg1::PC{Float64},arg2::Union(Ptr{PCFieldSplitSchurPreType},StridedArray{PCFieldSplitSchurPreType},Ptr{Void}),arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:PCFieldSplitGetSchurPre,petsc1),PetscErrorCode,(PC{Float64},Ptr{PCFieldSplitSchurPreType},Ptr{Mat{Float64}}),arg1,arg2,arg3)
    return err
end

function PCFieldSplitSetSchurFactType(arg1::PC{Float64},arg2::PCFieldSplitSchurFactType)
    err = ccall((:PCFieldSplitSetSchurFactType,petsc1),PetscErrorCode,(PC{Float64},PCFieldSplitSchurFactType),arg1,arg2)
    return err
end

function PCFieldSplitGetSchurBlocks(arg1::PC{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg5::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:PCFieldSplitGetSchurBlocks,petsc1),PetscErrorCode,(PC{Float64},Ptr{Mat{Float64}},Ptr{Mat{Float64}},Ptr{Mat{Float64}},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function PCFieldSplitSchurGetS(arg1::PC{Float64},S::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:PCFieldSplitSchurGetS,petsc1),PetscErrorCode,(PC{Float64},Ptr{Mat{Float64}}),arg1,S)
    return err
end

function PCFieldSplitSchurRestoreS(arg1::PC{Float64},S::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:PCFieldSplitSchurRestoreS,petsc1),PetscErrorCode,(PC{Float64},Ptr{Mat{Float64}}),arg1,S)
    return err
end

function PCGalerkinSetRestriction(arg1::PC{Float64},arg2::Mat{Float64})
    err = ccall((:PCGalerkinSetRestriction,petsc1),PetscErrorCode,(PC{Float64},Mat{Float64}),arg1,arg2)
    return err
end

function PCGalerkinSetInterpolation(arg1::PC{Float64},arg2::Mat{Float64})
    err = ccall((:PCGalerkinSetInterpolation,petsc1),PetscErrorCode,(PC{Float64},Mat{Float64}),arg1,arg2)
    return err
end

function PCSetCoordinates(arg1::PC{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:PCSetCoordinates,petsc1),PetscErrorCode,(PC{Float64},PetscInt,PetscInt,Ptr{Cint}),arg1,arg2,arg3,arg4)
    return err
end

function PCPythonSetType(arg1::PC{Float64},arg2::Union(ByteString,Symbol))
    err = ccall((:PCPythonSetType,petsc1),PetscErrorCode,(PC{Float64},Cstring),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function PCSetDM(arg1::PC{Float64},arg2::DM)
    ccall((:PCSetDM,petsc1),PetscErrorCode,(PC{Float64},DM),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PCGetDM(arg1::PC{Float64},arg2::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:PCGetDM,petsc1),PetscErrorCode,(PC{Float64},Ptr{DM}),arg1,arg2)
end 
=#
function PCSetApplicationContext(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PCSetApplicationContext,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
    return err
end

function PCGetApplicationContext(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PCGetApplicationContext,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
    return err
end

function PCBiCGStabCUSPSetTolerance(arg1::PC{Float64},PetscReal::Integer)
    err = ccall((:PCBiCGStabCUSPSetTolerance,petsc1),PetscErrorCode,(PC{Float64},Cint),arg1,PetscReal)
    return err
end

function PCBiCGStabCUSPSetIterations(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCBiCGStabCUSPSetIterations,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCBiCGStabCUSPSetUseVerboseMonitor(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCBiCGStabCUSPSetUseVerboseMonitor,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCAINVCUSPSetDropTolerance(arg1::PC{Float64},PetscReal::Integer)
    err = ccall((:PCAINVCUSPSetDropTolerance,petsc1),PetscErrorCode,(PC{Float64},Cint),arg1,PetscReal)
    return err
end

function PCAINVCUSPUseScaling(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCAINVCUSPUseScaling,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCAINVCUSPSetNonzeros(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCAINVCUSPSetNonzeros,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCAINVCUSPSetLinParameter(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCAINVCUSPSetLinParameter,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCPARMSSetGlobal(arg1::PC{Float64},arg2::PCPARMSGlobalType)
    err = ccall((:PCPARMSSetGlobal,petsc1),PetscErrorCode,(PC{Float64},PCPARMSGlobalType),arg1,arg2)
    return err
end

function PCPARMSSetLocal(arg1::PC{Float64},arg2::PCPARMSLocalType)
    err = ccall((:PCPARMSSetLocal,petsc1),PetscErrorCode,(PC{Float64},PCPARMSLocalType),arg1,arg2)
    return err
end

function PCPARMSSetSolveTolerances(arg1::PC{Float64},PetscReal::Integer,arg2::Integer)
    err = ccall((:PCPARMSSetSolveTolerances,petsc1),PetscErrorCode,(PC{Float64},Cint,PetscInt),arg1,PetscReal,arg2)
    return err
end

function PCPARMSSetSolveRestart(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCPARMSSetSolveRestart,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCPARMSSetNonsymPerm(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCPARMSSetNonsymPerm,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCPARMSSetFill(arg1::PC{Float64},arg2::Integer,arg3::Integer,arg4::Integer)
    err = ccall((:PCPARMSSetFill,petsc1),PetscErrorCode,(PC{Float64},PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
    return err
end

function PCGAMGSetType(arg1::PC{Float64},arg2::PCGAMGType)
    err = ccall((:PCGAMGSetType,petsc1),PetscErrorCode,(PC{Float64},Cstring),arg1,arg2)
    return err
end

function PCGAMGGetType(arg1::PC{Float64},arg2::Union(Ptr{PCGAMGType},StridedArray{PCGAMGType},Ptr{Void}))
    (arg2_,tmp) = symbol_get_before(arg2)
    err = ccall((:PCGAMGGetType,petsc1),PetscErrorCode,(PC{Float64},Ptr{Ptr{Uint8}}),arg1,arg2_)
    symbol_get_after(arg2_,arg2)
    return err
end

function PCGAMGSetProcEqLim(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCGAMGSetProcEqLim,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCGAMGSetRepartitioning(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCGAMGSetRepartitioning,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCGAMGSetUseASMAggs(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCGAMGSetUseASMAggs,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCGAMGSetSolverType(arg1::PC{Float64},arg2::Union(ByteString,Symbol),arg3::Integer)
    err = ccall((:PCGAMGSetSolverType,petsc1),PetscErrorCode,(PC{Float64},Cstring,PetscInt),arg1,arg2,arg3)
    return err
end

function PCGAMGSetThreshold(arg1::PC{Float64},PetscReal::Integer)
    err = ccall((:PCGAMGSetThreshold,petsc1),PetscErrorCode,(PC{Float64},Cint),arg1,PetscReal)
    return err
end

function PCGAMGSetCoarseEqLim(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCGAMGSetCoarseEqLim,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCGAMGSetNlevels(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCGAMGSetNlevels,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCGAMGSetNSmooths(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCGAMGSetNSmooths,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCGAMGSetSymGraph(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCGAMGSetSymGraph,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCGAMGSetSquareGraph(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCGAMGSetSquareGraph,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCGAMGSetReuseInterpolation(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCGAMGSetReuseInterpolation,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCGAMGFinalizePackage(arg0::Type{Float64})
    err = ccall((:PCGAMGFinalizePackage,petsc1),PetscErrorCode,())
    return err
end

function PCGAMGInitializePackage(arg0::Type{Float64})
    err = ccall((:PCGAMGInitializePackage,petsc1),PetscErrorCode,())
    return err
end

function PCGAMGRegister(arg0::Type{Float64},arg1::PCGAMGType,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PCGAMGRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

function PCGAMGClassicalSetType(arg1::PC{Float64},arg2::PCGAMGClassicalType)
    err = ccall((:PCGAMGClassicalSetType,petsc1),PetscErrorCode,(PC{Float64},Cstring),arg1,arg2)
    return err
end

function PCGAMGClassicalGetType(arg1::PC{Float64},arg2::Union(Ptr{PCGAMGClassicalType},StridedArray{PCGAMGClassicalType},Ptr{Void}))
    (arg2_,tmp) = symbol_get_before(arg2)
    err = ccall((:PCGAMGClassicalGetType,petsc1),PetscErrorCode,(PC{Float64},Ptr{Ptr{Uint8}}),arg1,arg2_)
    symbol_get_after(arg2_,arg2)
    return err
end

function PCBDDCSetChangeOfBasisMat(arg1::PC{Float64},arg2::Mat{Float64})
    err = ccall((:PCBDDCSetChangeOfBasisMat,petsc1),PetscErrorCode,(PC{Float64},Mat{Float64}),arg1,arg2)
    return err
end

function PCBDDCSetPrimalVerticesLocalIS(arg1::PC{Float64},arg2::IS{Float64})
    err = ccall((:PCBDDCSetPrimalVerticesLocalIS,petsc1),PetscErrorCode,(PC{Float64},IS{Float64}),arg1,arg2)
    return err
end

function PCBDDCSetCoarseningRatio(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCBDDCSetCoarseningRatio,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCBDDCSetLevels(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCBDDCSetLevels,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function PCBDDCSetNullSpace(arg1::PC{Float64},arg2::MatNullSpace)
    ccall((:PCBDDCSetNullSpace,petsc1),PetscErrorCode,(PC{Float64},MatNullSpace),arg1,arg2)
end 
=#
function PCBDDCSetDirichletBoundaries(arg1::PC{Float64},arg2::IS{Float64})
    err = ccall((:PCBDDCSetDirichletBoundaries,petsc1),PetscErrorCode,(PC{Float64},IS{Float64}),arg1,arg2)
    return err
end

function PCBDDCSetDirichletBoundariesLocal(arg1::PC{Float64},arg2::IS{Float64})
    err = ccall((:PCBDDCSetDirichletBoundariesLocal,petsc1),PetscErrorCode,(PC{Float64},IS{Float64}),arg1,arg2)
    return err
end

function PCBDDCGetDirichletBoundaries(arg1::PC{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:PCBDDCGetDirichletBoundaries,petsc1),PetscErrorCode,(PC{Float64},Ptr{IS{Float64}}),arg1,arg2)
    return err
end

function PCBDDCGetDirichletBoundariesLocal(arg1::PC{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:PCBDDCGetDirichletBoundariesLocal,petsc1),PetscErrorCode,(PC{Float64},Ptr{IS{Float64}}),arg1,arg2)
    return err
end

function PCBDDCSetNeumannBoundaries(arg1::PC{Float64},arg2::IS{Float64})
    err = ccall((:PCBDDCSetNeumannBoundaries,petsc1),PetscErrorCode,(PC{Float64},IS{Float64}),arg1,arg2)
    return err
end

function PCBDDCSetNeumannBoundariesLocal(arg1::PC{Float64},arg2::IS{Float64})
    err = ccall((:PCBDDCSetNeumannBoundariesLocal,petsc1),PetscErrorCode,(PC{Float64},IS{Float64}),arg1,arg2)
    return err
end

function PCBDDCGetNeumannBoundaries(arg1::PC{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:PCBDDCGetNeumannBoundaries,petsc1),PetscErrorCode,(PC{Float64},Ptr{IS{Float64}}),arg1,arg2)
    return err
end

function PCBDDCGetNeumannBoundariesLocal(arg1::PC{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:PCBDDCGetNeumannBoundariesLocal,petsc1),PetscErrorCode,(PC{Float64},Ptr{IS{Float64}}),arg1,arg2)
    return err
end

function PCBDDCSetDofsSplitting(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:PCBDDCSetDofsSplitting,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Ptr{IS{Float64}}),arg1,arg2,arg3)
    return err
end

function PCBDDCSetDofsSplittingLocal(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    err = ccall((:PCBDDCSetDofsSplittingLocal,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Ptr{IS{Float64}}),arg1,arg2,arg3)
    return err
end

function PCBDDCSetLocalAdjacencyGraph(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),PetscCopyMode::Integer)
    err = ccall((:PCBDDCSetLocalAdjacencyGraph,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Ptr{PetscInt},Ptr{PetscInt},Cint),arg1,arg2,arg3,arg4,PetscCopyMode)
    return err
end

function PCBDDCCreateFETIDPOperators(arg1::PC{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg3::Union(Ptr{PC{Float64}},StridedArray{PC{Float64}},Ptr{Void}))
    err = ccall((:PCBDDCCreateFETIDPOperators,petsc1),PetscErrorCode,(PC{Float64},Ptr{Mat{Float64}},Ptr{PC{Float64}}),arg1,arg2,arg3)
    return err
end

function PCBDDCMatFETIDPGetRHS(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:PCBDDCMatFETIDPGetRHS,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function PCBDDCMatFETIDPGetSolution(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:PCBDDCMatFETIDPGetSolution,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function PCISSetUseStiffnessScaling(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCISSetUseStiffnessScaling,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCISSetSubdomainScalingFactor(arg1::PC{Float64},arg2::Float64)
    err = ccall((:PCISSetSubdomainScalingFactor,petsc1),PetscErrorCode,(PC{Float64},Float64),arg1,arg2)
    return err
end

function PCISSetSubdomainDiagonalScaling(arg1::PC{Float64},arg2::Vec{Float64})
    err = ccall((:PCISSetSubdomainDiagonalScaling,petsc1),PetscErrorCode,(PC{Float64},Vec{Float64}),arg1,arg2)
    return err
end

function PCMGSetType(arg1::PC{Float64},arg2::PCMGType)
    err = ccall((:PCMGSetType,petsc1),PetscErrorCode,(PC{Float64},PCMGType),arg1,arg2)
    return err
end

function PCMGGetType(arg1::PC{Float64},arg2::Union(Ptr{PCMGType},StridedArray{PCMGType},Ptr{Void}))
    err = ccall((:PCMGGetType,petsc1),PetscErrorCode,(PC{Float64},Ptr{PCMGType}),arg1,arg2)
    return err
end

function PCMGSetLevels(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{MPI_Comm},StridedArray{MPI_Comm},Ptr{Void}))
    err = ccall((:PCMGSetLevels,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Ptr{comm_type}),arg1,arg2,arg3)
    return err
end

function PCMGGetLevels(arg1::PC{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:PCMGGetLevels,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function PCMGSetNumberSmoothUp(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCMGSetNumberSmoothUp,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCMGSetNumberSmoothDown(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCMGSetNumberSmoothDown,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCMGSetCycleType(arg1::PC{Float64},arg2::PCMGCycleType)
    err = ccall((:PCMGSetCycleType,petsc1),PetscErrorCode,(PC{Float64},PCMGCycleType),arg1,arg2)
    return err
end

function PCMGSetCycleTypeOnLevel(arg1::PC{Float64},arg2::Integer,arg3::PCMGCycleType)
    err = ccall((:PCMGSetCycleTypeOnLevel,petsc1),PetscErrorCode,(PC{Float64},PetscInt,PCMGCycleType),arg1,arg2,arg3)
    return err
end

function PCMGSetCyclesOnLevel(arg1::PC{Float64},arg2::Integer,arg3::Integer)
    err = ccall((:PCMGSetCyclesOnLevel,petsc1),PetscErrorCode,(PC{Float64},PetscInt,PetscInt),arg1,arg2,arg3)
    return err
end

function PCMGMultiplicativeSetCycles(arg1::PC{Float64},arg2::Integer)
    err = ccall((:PCMGMultiplicativeSetCycles,petsc1),PetscErrorCode,(PC{Float64},PetscInt),arg1,arg2)
    return err
end

function PCMGSetGalerkin(arg1::PC{Float64},arg2::PetscBool)
    err = ccall((:PCMGSetGalerkin,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
    return err
end

function PCMGGetGalerkin(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:PCMGGetGalerkin,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function PCMGSetRhs(arg1::PC{Float64},arg2::Integer,arg3::Vec{Float64})
    err = ccall((:PCMGSetRhs,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Vec{Float64}),arg1,arg2,arg3)
    return err
end

function PCMGSetX(arg1::PC{Float64},arg2::Integer,arg3::Vec{Float64})
    err = ccall((:PCMGSetX,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Vec{Float64}),arg1,arg2,arg3)
    return err
end

function PCMGSetR(arg1::PC{Float64},arg2::Integer,arg3::Vec{Float64})
    err = ccall((:PCMGSetR,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Vec{Float64}),arg1,arg2,arg3)
    return err
end

function PCMGSetRestriction(arg1::PC{Float64},arg2::Integer,arg3::Mat{Float64})
    err = ccall((:PCMGSetRestriction,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Mat{Float64}),arg1,arg2,arg3)
    return err
end

function PCMGGetRestriction(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:PCMGGetRestriction,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Ptr{Mat{Float64}}),arg1,arg2,arg3)
    return err
end

function PCMGSetInterpolation(arg1::PC{Float64},arg2::Integer,arg3::Mat{Float64})
    err = ccall((:PCMGSetInterpolation,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Mat{Float64}),arg1,arg2,arg3)
    return err
end

function PCMGGetInterpolation(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:PCMGGetInterpolation,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Ptr{Mat{Float64}}),arg1,arg2,arg3)
    return err
end

function PCMGSetRScale(arg1::PC{Float64},arg2::Integer,arg3::Vec{Float64})
    err = ccall((:PCMGSetRScale,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Vec{Float64}),arg1,arg2,arg3)
    return err
end

function PCMGGetRScale(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:PCMGGetRScale,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Ptr{Vec{Float64}}),arg1,arg2,arg3)
    return err
end

function PCMGSetResidual(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Mat{Float64})
    err = ccall((:PCMGSetResidual,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Ptr{Void},Mat{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function PCMGResidualDefault(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64})
    err = ccall((:PCMGResidualDefault,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function KSPInitializePackage(arg0::Type{Float64})
    err = ccall((:KSPInitializePackage,petsc1),PetscErrorCode,())
    return err
end

function KSPCreate(arg1::MPI_Comm,arg2::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    err = ccall((:KSPCreate,petsc1),PetscErrorCode,(comm_type,Ptr{KSP{Float64}}),arg1.val,arg2)
    return err
end

function KSPSetType(arg1::KSP{Float64},arg2::KSPType)
    err = ccall((:KSPSetType,petsc1),PetscErrorCode,(KSP{Float64},Cstring),arg1,arg2)
    return err
end

function KSPGetType(arg1::KSP{Float64},arg2::Union(Ptr{KSPType},StridedArray{KSPType},Ptr{Void}))
    (arg2_,tmp) = symbol_get_before(arg2)
    err = ccall((:KSPGetType,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Ptr{Uint8}}),arg1,arg2_)
    symbol_get_after(arg2_,arg2)
    return err
end

function KSPSetUp(arg1::KSP{Float64})
    err = ccall((:KSPSetUp,petsc1),PetscErrorCode,(KSP{Float64},),arg1)
    return err
end

function KSPSetUpOnBlocks(arg1::KSP{Float64})
    err = ccall((:KSPSetUpOnBlocks,petsc1),PetscErrorCode,(KSP{Float64},),arg1)
    return err
end

function KSPSolve(arg1::KSP{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:KSPSolve,petsc1),PetscErrorCode,(KSP{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function KSPSolveTranspose(arg1::KSP{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:KSPSolveTranspose,petsc1),PetscErrorCode,(KSP{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function KSPReset(arg1::KSP{Float64})
    err = ccall((:KSPReset,petsc1),PetscErrorCode,(KSP{Float64},),arg1)
    return err
end

function KSPDestroy(arg1::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    err = ccall((:KSPDestroy,petsc1),PetscErrorCode,(Ptr{KSP{Float64}},),arg1)
    return err
end

function KSPSetReusePreconditioner(arg1::KSP{Float64},arg2::PetscBool)
    err = ccall((:KSPSetReusePreconditioner,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
    return err
end

function KSPSetSkipPCSetFromOptions(arg1::KSP{Float64},arg2::PetscBool)
    err = ccall((:KSPSetSkipPCSetFromOptions,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
    return err
end

function KSPRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

function KSPSetPCSide(arg1::KSP{Float64},arg2::PCSide)
    err = ccall((:KSPSetPCSide,petsc1),PetscErrorCode,(KSP{Float64},PCSide),arg1,arg2)
    return err
end

function KSPGetPCSide(arg1::KSP{Float64},arg2::Union(Ptr{PCSide},StridedArray{PCSide},Ptr{Void}))
    err = ccall((:KSPGetPCSide,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PCSide}),arg1,arg2)
    return err
end

function KSPSetTolerances(arg1::KSP{Float64},PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Integer)
    err = ccall((:KSPSetTolerances,petsc1),PetscErrorCode,(KSP{Float64},Cint,Cint,Cint,PetscInt),arg1,PetscReal,arg2,arg3,arg4)
    return err
end

function KSPGetTolerances(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:KSPGetTolerances,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function KSPSetInitialGuessNonzero(arg1::KSP{Float64},arg2::PetscBool)
    err = ccall((:KSPSetInitialGuessNonzero,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
    return err
end

function KSPGetInitialGuessNonzero(arg1::KSP{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:KSPGetInitialGuessNonzero,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function KSPSetInitialGuessKnoll(arg1::KSP{Float64},arg2::PetscBool)
    err = ccall((:KSPSetInitialGuessKnoll,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
    return err
end

function KSPGetInitialGuessKnoll(arg1::KSP{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:KSPGetInitialGuessKnoll,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function KSPSetErrorIfNotConverged(arg1::KSP{Float64},arg2::PetscBool)
    err = ccall((:KSPSetErrorIfNotConverged,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
    return err
end

function KSPGetErrorIfNotConverged(arg1::KSP{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:KSPGetErrorIfNotConverged,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function KSPSetComputeEigenvalues(arg1::KSP{Float64},arg2::PetscBool)
    err = ccall((:KSPSetComputeEigenvalues,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
    return err
end

function KSPGetComputeEigenvalues(arg1::KSP{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:KSPGetComputeEigenvalues,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function KSPSetComputeSingularValues(arg1::KSP{Float64},arg2::PetscBool)
    err = ccall((:KSPSetComputeSingularValues,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
    return err
end

function KSPGetComputeSingularValues(arg1::KSP{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:KSPGetComputeSingularValues,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function KSPGetRhs(arg1::KSP{Float64},arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:KSPGetRhs,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Vec{Float64}}),arg1,arg2)
    return err
end

function KSPGetSolution(arg1::KSP{Float64},arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:KSPGetSolution,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Vec{Float64}}),arg1,arg2)
    return err
end

function KSPGetResidualNorm(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:KSPGetResidualNorm,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint}),arg1,arg2)
    return err
end

function KSPGetIterationNumber(arg1::KSP{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:KSPGetIterationNumber,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function KSPGetTotalIterations(arg1::KSP{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:KSPGetTotalIterations,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function KSPCreateVecs(arg1::KSP{Float64},arg2::Integer,arg3::Union(Ptr{Ptr{Vec{Float64}}},StridedArray{Ptr{Vec{Float64}}},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{Ptr{Vec{Float64}}},StridedArray{Ptr{Vec{Float64}}},Ptr{Void}))
    err = ccall((:KSPCreateVecs,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,Ptr{Ptr{Vec{Float64}}},PetscInt,Ptr{Ptr{Vec{Float64}}}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function KSPSetPostSolve(arg1::KSP{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPSetPostSolve,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
    return err
end

function KSPSetPC(arg1::KSP{Float64},arg2::PC{Float64})
    err = ccall((:KSPSetPC,petsc1),PetscErrorCode,(KSP{Float64},PC{Float64}),arg1,arg2)
    return err
end

function KSPGetPC(arg1::KSP{Float64},arg2::Union(Ptr{PC{Float64}},StridedArray{PC{Float64}},Ptr{Void}))
    err = ccall((:KSPGetPC,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PC{Float64}}),arg1,arg2)
    return err
end

function KSPMonitor(arg1::KSP{Float64},arg2::Integer,PetscReal::Integer)
    err = ccall((:KSPMonitor,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,Cint),arg1,arg2,PetscReal)
    return err
end

function KSPMonitorSet(arg1::KSP{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPMonitorSet,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
    return err
end

function KSPMonitorCancel(arg1::KSP{Float64})
    err = ccall((:KSPMonitorCancel,petsc1),PetscErrorCode,(KSP{Float64},),arg1)
    return err
end

function KSPGetMonitorContext(arg1::KSP{Float64},arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    err = ccall((:KSPGetMonitorContext,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Ptr{Void}}),arg1,arg2)
    return err
end

function KSPGetResidualHistory(arg1::KSP{Float64},arg2::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:KSPGetResidualHistory,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Ptr{Cint}},Ptr{PetscInt}),arg1,arg2,arg3)
    return err
end

function KSPSetResidualHistory(arg1::KSP{Float64},PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg2::Integer,arg3::PetscBool)
    err = ccall((:KSPSetResidualHistory,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint},PetscInt,PetscBool),arg1,PetscReal,arg2,arg3)
    return err
end

function KSPBuildSolutionDefault(arg1::KSP{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:KSPBuildSolutionDefault,petsc1),PetscErrorCode,(KSP{Float64},Vec{Float64},Ptr{Vec{Float64}}),arg1,arg2,arg3)
    return err
end

function KSPBuildResidualDefault(arg1::KSP{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:KSPBuildResidualDefault,petsc1),PetscErrorCode,(KSP{Float64},Vec{Float64},Vec{Float64},Ptr{Vec{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function KSPDestroyDefault(arg1::KSP{Float64})
    err = ccall((:KSPDestroyDefault,petsc1),PetscErrorCode,(KSP{Float64},),arg1)
    return err
end

function KSPSetWorkVecs(arg1::KSP{Float64},arg2::Integer)
    err = ccall((:KSPSetWorkVecs,petsc1),PetscErrorCode,(KSP{Float64},PetscInt),arg1,arg2)
    return err
end

function PCKSPGetKSP(arg1::PC{Float64},arg2::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    err = ccall((:PCKSPGetKSP,petsc1),PetscErrorCode,(PC{Float64},Ptr{KSP{Float64}}),arg1,arg2)
    return err
end

function PCBJacobiGetSubKSP(arg1::PC{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Ptr{KSP{Float64}}},StridedArray{Ptr{KSP{Float64}}},Ptr{Void}))
    err = ccall((:PCBJacobiGetSubKSP,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{KSP{Float64}}}),arg1,arg2,arg3,arg4)
    return err
end

function PCASMGetSubKSP(arg1::PC{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Ptr{KSP{Float64}}},StridedArray{Ptr{KSP{Float64}}},Ptr{Void}))
    err = ccall((:PCASMGetSubKSP,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{KSP{Float64}}}),arg1,arg2,arg3,arg4)
    return err
end

function PCGASMGetSubKSP(arg1::PC{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Ptr{KSP{Float64}}},StridedArray{Ptr{KSP{Float64}}},Ptr{Void}))
    err = ccall((:PCGASMGetSubKSP,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{KSP{Float64}}}),arg1,arg2,arg3,arg4)
    return err
end

function PCFieldSplitGetSubKSP(arg1::PC{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{KSP{Float64}}},StridedArray{Ptr{KSP{Float64}}},Ptr{Void}))
    err = ccall((:PCFieldSplitGetSubKSP,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscInt},Ptr{Ptr{KSP{Float64}}}),arg1,arg2,arg3)
    return err
end

function PCMGGetSmoother(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    err = ccall((:PCMGGetSmoother,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Ptr{KSP{Float64}}),arg1,arg2,arg3)
    return err
end

function PCMGGetSmootherDown(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    err = ccall((:PCMGGetSmootherDown,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Ptr{KSP{Float64}}),arg1,arg2,arg3)
    return err
end

function PCMGGetSmootherUp(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    err = ccall((:PCMGGetSmootherUp,petsc1),PetscErrorCode,(PC{Float64},PetscInt,Ptr{KSP{Float64}}),arg1,arg2,arg3)
    return err
end

function PCMGGetCoarseSolve(arg1::PC{Float64},arg2::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    err = ccall((:PCMGGetCoarseSolve,petsc1),PetscErrorCode,(PC{Float64},Ptr{KSP{Float64}}),arg1,arg2)
    return err
end

function PCGalerkinGetKSP(arg1::PC{Float64},arg2::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    err = ccall((:PCGalerkinGetKSP,petsc1),PetscErrorCode,(PC{Float64},Ptr{KSP{Float64}}),arg1,arg2)
    return err
end

function KSPBuildSolution(arg1::KSP{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:KSPBuildSolution,petsc1),PetscErrorCode,(KSP{Float64},Vec{Float64},Ptr{Vec{Float64}}),arg1,arg2,arg3)
    return err
end

function KSPBuildResidual(arg1::KSP{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:KSPBuildResidual,petsc1),PetscErrorCode,(KSP{Float64},Vec{Float64},Vec{Float64},Ptr{Vec{Float64}}),arg1,arg2,arg3,arg4)
    return err
end

function KSPRichardsonSetScale(arg1::KSP{Float64},PetscReal::Integer)
    err = ccall((:KSPRichardsonSetScale,petsc1),PetscErrorCode,(KSP{Float64},Cint),arg1,PetscReal)
    return err
end

function KSPRichardsonSetSelfScale(arg1::KSP{Float64},arg2::PetscBool)
    err = ccall((:KSPRichardsonSetSelfScale,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
    return err
end

function KSPChebyshevSetEigenvalues(arg1::KSP{Float64},PetscReal::Integer,arg2::Integer)
    err = ccall((:KSPChebyshevSetEigenvalues,petsc1),PetscErrorCode,(KSP{Float64},Cint,Cint),arg1,PetscReal,arg2)
    return err
end

function KSPChebyshevEstEigSet(arg1::KSP{Float64},PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Integer)
    err = ccall((:KSPChebyshevEstEigSet,petsc1),PetscErrorCode,(KSP{Float64},Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4)
    return err
end

#= skipping function with undefined symbols: 
 function KSPChebyshevEstEigSetRandom(arg1::KSP{Float64},arg2::PetscRandom)
    ccall((:KSPChebyshevEstEigSetRandom,petsc1),PetscErrorCode,(KSP{Float64},PetscRandom),arg1,arg2)
end 
=#
function KSPChebyshevEstEigGetKSP(arg1::KSP{Float64},arg2::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    err = ccall((:KSPChebyshevEstEigGetKSP,petsc1),PetscErrorCode,(KSP{Float64},Ptr{KSP{Float64}}),arg1,arg2)
    return err
end

function KSPComputeExtremeSingularValues(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:KSPComputeExtremeSingularValues,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3)
    return err
end

function KSPComputeEigenvalues(arg1::KSP{Float64},arg2::Integer,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:KSPComputeEigenvalues,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,Ptr{Cint},Ptr{Cint},Ptr{PetscInt}),arg1,arg2,PetscReal,arg3,arg4)
    return err
end

function KSPComputeEigenvaluesExplicitly(arg1::KSP{Float64},arg2::Integer,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:KSPComputeEigenvaluesExplicitly,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,Ptr{Cint},Ptr{Cint}),arg1,arg2,PetscReal,arg3)
    return err
end

function KSPFCGSetMmax(arg1::KSP{Float64},arg2::Integer)
    err = ccall((:KSPFCGSetMmax,petsc1),PetscErrorCode,(KSP{Float64},PetscInt),arg1,arg2)
    return err
end

function KSPFCGGetMmax(arg1::KSP{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:KSPFCGGetMmax,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function KSPFCGSetNprealloc(arg1::KSP{Float64},arg2::Integer)
    err = ccall((:KSPFCGSetNprealloc,petsc1),PetscErrorCode,(KSP{Float64},PetscInt),arg1,arg2)
    return err
end

function KSPFCGGetNprealloc(arg1::KSP{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:KSPFCGGetNprealloc,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function KSPFCGSetTruncationType(arg1::KSP{Float64},arg2::KSPFCGTruncationType)
    err = ccall((:KSPFCGSetTruncationType,petsc1),PetscErrorCode,(KSP{Float64},KSPFCGTruncationType),arg1,arg2)
    return err
end

function KSPFCGGetTruncationType(arg1::KSP{Float64},arg2::Union(Ptr{KSPFCGTruncationType},StridedArray{KSPFCGTruncationType},Ptr{Void}))
    err = ccall((:KSPFCGGetTruncationType,petsc1),PetscErrorCode,(KSP{Float64},Ptr{KSPFCGTruncationType}),arg1,arg2)
    return err
end

function KSPGMRESSetRestart(arg1::KSP{Float64},arg2::Integer)
    err = ccall((:KSPGMRESSetRestart,petsc1),PetscErrorCode,(KSP{Float64},PetscInt),arg1,arg2)
    return err
end

function KSPGMRESGetRestart(arg1::KSP{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:KSPGMRESGetRestart,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function KSPGMRESSetHapTol(arg1::KSP{Float64},PetscReal::Integer)
    err = ccall((:KSPGMRESSetHapTol,petsc1),PetscErrorCode,(KSP{Float64},Cint),arg1,PetscReal)
    return err
end

function KSPGMRESSetPreAllocateVectors(arg1::KSP{Float64})
    err = ccall((:KSPGMRESSetPreAllocateVectors,petsc1),PetscErrorCode,(KSP{Float64},),arg1)
    return err
end

function KSPGMRESSetOrthogonalization(arg1::KSP{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPGMRESSetOrthogonalization,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Void}),arg1,arg2)
    return err
end

function KSPGMRESGetOrthogonalization(arg1::KSP{Float64},arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    err = ccall((:KSPGMRESGetOrthogonalization,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Ptr{Void}}),arg1,arg2)
    return err
end

function KSPGMRESModifiedGramSchmidtOrthogonalization(arg1::KSP{Float64},arg2::Integer)
    err = ccall((:KSPGMRESModifiedGramSchmidtOrthogonalization,petsc1),PetscErrorCode,(KSP{Float64},PetscInt),arg1,arg2)
    return err
end

function KSPGMRESClassicalGramSchmidtOrthogonalization(arg1::KSP{Float64},arg2::Integer)
    err = ccall((:KSPGMRESClassicalGramSchmidtOrthogonalization,petsc1),PetscErrorCode,(KSP{Float64},PetscInt),arg1,arg2)
    return err
end

function KSPLGMRESSetAugDim(arg1::KSP{Float64},arg2::Integer)
    err = ccall((:KSPLGMRESSetAugDim,petsc1),PetscErrorCode,(KSP{Float64},PetscInt),arg1,arg2)
    return err
end

function KSPLGMRESSetConstant(arg1::KSP{Float64})
    err = ccall((:KSPLGMRESSetConstant,petsc1),PetscErrorCode,(KSP{Float64},),arg1)
    return err
end

function KSPGCRSetRestart(arg1::KSP{Float64},arg2::Integer)
    err = ccall((:KSPGCRSetRestart,petsc1),PetscErrorCode,(KSP{Float64},PetscInt),arg1,arg2)
    return err
end

function KSPGCRGetRestart(arg1::KSP{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:KSPGCRGetRestart,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function KSPGCRSetModifyPC(arg1::KSP{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPGCRSetModifyPC,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
    return err
end

function KSPGMRESSetCGSRefinementType(arg1::KSP{Float64},arg2::KSPGMRESCGSRefinementType)
    err = ccall((:KSPGMRESSetCGSRefinementType,petsc1),PetscErrorCode,(KSP{Float64},KSPGMRESCGSRefinementType),arg1,arg2)
    return err
end

function KSPGMRESGetCGSRefinementType(arg1::KSP{Float64},arg2::Union(Ptr{KSPGMRESCGSRefinementType},StridedArray{KSPGMRESCGSRefinementType},Ptr{Void}))
    err = ccall((:KSPGMRESGetCGSRefinementType,petsc1),PetscErrorCode,(KSP{Float64},Ptr{KSPGMRESCGSRefinementType}),arg1,arg2)
    return err
end

function KSPFGMRESModifyPCNoChange(arg1::KSP{Float64},arg2::Integer,arg3::Integer,PetscReal::Integer,arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPFGMRESModifyPCNoChange,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,PetscInt,Cint,Ptr{Void}),arg1,arg2,arg3,PetscReal,arg4)
    return err
end

function KSPFGMRESModifyPCKSP(arg1::KSP{Float64},arg2::Integer,arg3::Integer,PetscReal::Integer,arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPFGMRESModifyPCKSP,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,PetscInt,Cint,Ptr{Void}),arg1,arg2,arg3,PetscReal,arg4)
    return err
end

function KSPFGMRESSetModifyPC(arg1::KSP{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPFGMRESSetModifyPC,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
    return err
end

function KSPQCGSetTrustRegionRadius(arg1::KSP{Float64},PetscReal::Integer)
    err = ccall((:KSPQCGSetTrustRegionRadius,petsc1),PetscErrorCode,(KSP{Float64},Cint),arg1,PetscReal)
    return err
end

function KSPQCGGetQuadratic(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:KSPQCGGetQuadratic,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint}),arg1,arg2)
    return err
end

function KSPQCGGetTrialStepNorm(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:KSPQCGGetTrialStepNorm,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint}),arg1,arg2)
    return err
end

function KSPBCGSLSetXRes(arg1::KSP{Float64},PetscReal::Integer)
    err = ccall((:KSPBCGSLSetXRes,petsc1),PetscErrorCode,(KSP{Float64},Cint),arg1,PetscReal)
    return err
end

function KSPBCGSLSetPol(arg1::KSP{Float64},arg2::PetscBool)
    err = ccall((:KSPBCGSLSetPol,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
    return err
end

function KSPBCGSLSetEll(arg1::KSP{Float64},arg2::Integer)
    err = ccall((:KSPBCGSLSetEll,petsc1),PetscErrorCode,(KSP{Float64},PetscInt),arg1,arg2)
    return err
end

function KSPBCGSLSetUsePseudoinverse(arg1::KSP{Float64},arg2::PetscBool)
    err = ccall((:KSPBCGSLSetUsePseudoinverse,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
    return err
end

function KSPSetFromOptions(arg1::KSP{Float64})
    err = ccall((:KSPSetFromOptions,petsc1),PetscErrorCode,(KSP{Float64},),arg1)
    return err
end

function KSPAddOptionsChecker(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPAddOptionsChecker,petsc1),PetscErrorCode,(Ptr{Void},),arg1)
    return err
end

function KSPMonitorSingularValue(arg1::KSP{Float64},arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPMonitorSingularValue,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
    return err
end

function KSPMonitorDefault(arg1::KSP{Float64},arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPMonitorDefault,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
    return err
end

function KSPLSQRMonitorDefault(arg1::KSP{Float64},arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPLSQRMonitorDefault,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
    return err
end

function KSPMonitorRange(arg1::KSP{Float64},arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPMonitorRange,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
    return err
end

function KSPMonitorDynamicTolerance(ksp::KSP{Float64},its::Integer,fnorm::Integer,dummy::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPMonitorDynamicTolerance,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,Cint,Ptr{Void}),ksp,its,fnorm,dummy)
    return err
end

function KSPMonitorDynamicToleranceDestroy(arg0::Type{Float64},dummy::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    err = ccall((:KSPMonitorDynamicToleranceDestroy,petsc1),PetscErrorCode,(Ptr{Ptr{Void}},),dummy)
    return err
end

function KSPMonitorTrueResidualNorm(arg1::KSP{Float64},arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPMonitorTrueResidualNorm,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
    return err
end

function KSPMonitorTrueResidualMaxNorm(arg1::KSP{Float64},arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPMonitorTrueResidualMaxNorm,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
    return err
end

function KSPMonitorDefaultShort(arg1::KSP{Float64},arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPMonitorDefaultShort,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
    return err
end

function KSPMonitorSolution(arg1::KSP{Float64},arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPMonitorSolution,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
    return err
end

function KSPMonitorSAWs(arg1::KSP{Float64},arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPMonitorSAWs,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
    return err
end

function KSPMonitorSAWsCreate(arg1::KSP{Float64},arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    err = ccall((:KSPMonitorSAWsCreate,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Ptr{Void}}),arg1,arg2)
    return err
end

function KSPMonitorSAWsDestroy(arg0::Type{Float64},arg1::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    err = ccall((:KSPMonitorSAWsDestroy,petsc1),PetscErrorCode,(Ptr{Ptr{Void}},),arg1)
    return err
end

function KSPGMRESMonitorKrylov(arg1::KSP{Float64},arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPGMRESMonitorKrylov,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
    return err
end

function KSPUnwindPreconditioner(arg1::KSP{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    err = ccall((:KSPUnwindPreconditioner,petsc1),PetscErrorCode,(KSP{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
    return err
end

function KSPInitialResidual(arg1::KSP{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64},arg5::Vec{Float64},arg6::Vec{Float64})
    err = ccall((:KSPInitialResidual,petsc1),PetscErrorCode,(KSP{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function KSPSetOperators(arg1::KSP{Float64},arg2::Mat{Float64},arg3::Mat{Float64})
    err = ccall((:KSPSetOperators,petsc1),PetscErrorCode,(KSP{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3)
    return err
end

function KSPGetOperators(arg1::KSP{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:KSPGetOperators,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Mat{Float64}},Ptr{Mat{Float64}}),arg1,arg2,arg3)
    return err
end

function KSPGetOperatorsSet(arg1::KSP{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:KSPGetOperatorsSet,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3)
    return err
end

function KSPSetOptionsPrefix(arg1::KSP{Float64},arg2::Union(ByteString,Symbol))
    err = ccall((:KSPSetOptionsPrefix,petsc1),PetscErrorCode,(KSP{Float64},Cstring),arg1,arg2)
    return err
end

function KSPAppendOptionsPrefix(arg1::KSP{Float64},arg2::Union(ByteString,Symbol))
    err = ccall((:KSPAppendOptionsPrefix,petsc1),PetscErrorCode,(KSP{Float64},Cstring),arg1,arg2)
    return err
end

function KSPGetOptionsPrefix(arg1::KSP{Float64},arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    err = ccall((:KSPGetOptionsPrefix,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Ptr{Uint8}}),arg1,arg2)
    return err
end

function KSPSetTabLevel(arg1::KSP{Float64},arg2::Integer)
    err = ccall((:KSPSetTabLevel,petsc1),PetscErrorCode,(KSP{Float64},PetscInt),arg1,arg2)
    return err
end

function KSPGetTabLevel(arg1::KSP{Float64},arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    err = ccall((:KSPGetTabLevel,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PetscInt}),arg1,arg2)
    return err
end

function KSPSetDiagonalScale(arg1::KSP{Float64},arg2::PetscBool)
    err = ccall((:KSPSetDiagonalScale,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
    return err
end

function KSPGetDiagonalScale(arg1::KSP{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:KSPGetDiagonalScale,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function KSPSetDiagonalScaleFix(arg1::KSP{Float64},arg2::PetscBool)
    err = ccall((:KSPSetDiagonalScaleFix,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
    return err
end

function KSPGetDiagonalScaleFix(arg1::KSP{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    err = ccall((:KSPGetDiagonalScaleFix,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PetscBool}),arg1,arg2)
    return err
end

function KSPView(arg1::KSP{Float64},arg2::PetscViewer{Float64})
    err = ccall((:KSPView,petsc1),PetscErrorCode,(KSP{Float64},PetscViewer{Float64}),arg1,arg2)
    return err
end

function KSPLoad(arg1::KSP{Float64},arg2::PetscViewer{Float64})
    err = ccall((:KSPLoad,petsc1),PetscErrorCode,(KSP{Float64},PetscViewer{Float64}),arg1,arg2)
    return err
end

function KSPReasonViewFromOptions(arg1::KSP{Float64})
    err = ccall((:KSPReasonViewFromOptions,petsc1),PetscErrorCode,(KSP{Float64},),arg1)
    return err
end

function KSPLSQRSetStandardErrorVec(arg1::KSP{Float64},arg2::Vec{Float64})
    err = ccall((:KSPLSQRSetStandardErrorVec,petsc1),PetscErrorCode,(KSP{Float64},Vec{Float64}),arg1,arg2)
    return err
end

function KSPLSQRGetStandardErrorVec(arg1::KSP{Float64},arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    err = ccall((:KSPLSQRGetStandardErrorVec,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Vec{Float64}}),arg1,arg2)
    return err
end

function PCRedundantGetKSP(arg1::PC{Float64},arg2::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    err = ccall((:PCRedundantGetKSP,petsc1),PetscErrorCode,(PC{Float64},Ptr{KSP{Float64}}),arg1,arg2)
    return err
end

function PCRedistributeGetKSP(arg1::PC{Float64},arg2::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    err = ccall((:PCRedistributeGetKSP,petsc1),PetscErrorCode,(PC{Float64},Ptr{KSP{Float64}}),arg1,arg2)
    return err
end

function KSPSetNormType(arg1::KSP{Float64},arg2::KSPNormType)
    err = ccall((:KSPSetNormType,petsc1),PetscErrorCode,(KSP{Float64},KSPNormType),arg1,arg2)
    return err
end

function KSPGetNormType(arg1::KSP{Float64},arg2::Union(Ptr{KSPNormType},StridedArray{KSPNormType},Ptr{Void}))
    err = ccall((:KSPGetNormType,petsc1),PetscErrorCode,(KSP{Float64},Ptr{KSPNormType}),arg1,arg2)
    return err
end

function KSPSetSupportedNorm(ksp::KSP{Float64},arg1::KSPNormType,arg2::PCSide,arg3::Integer)
    err = ccall((:KSPSetSupportedNorm,petsc1),PetscErrorCode,(KSP{Float64},KSPNormType,PCSide,PetscInt),ksp,arg1,arg2,arg3)
    return err
end

function KSPSetCheckNormIteration(arg1::KSP{Float64},arg2::Integer)
    err = ccall((:KSPSetCheckNormIteration,petsc1),PetscErrorCode,(KSP{Float64},PetscInt),arg1,arg2)
    return err
end

function KSPSetLagNorm(arg1::KSP{Float64},arg2::PetscBool)
    err = ccall((:KSPSetLagNorm,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
    return err
end

function KSPSetConvergenceTest(arg1::KSP{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPSetConvergenceTest,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
    return err
end

function KSPGetConvergenceContext(arg1::KSP{Float64},arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    err = ccall((:KSPGetConvergenceContext,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Ptr{Void}}),arg1,arg2)
    return err
end

function KSPConvergedDefault(arg1::KSP{Float64},arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{KSPConvergedReason},StridedArray{KSPConvergedReason},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPConvergedDefault,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,Cint,Ptr{KSPConvergedReason},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
    return err
end

function KSPConvergedLSQR(arg1::KSP{Float64},arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{KSPConvergedReason},StridedArray{KSPConvergedReason},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPConvergedLSQR,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,Cint,Ptr{KSPConvergedReason},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
    return err
end

function KSPConvergedDefaultDestroy(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPConvergedDefaultDestroy,petsc1),PetscErrorCode,(Ptr{Void},),arg1)
    return err
end

function KSPConvergedDefaultCreate(arg0::Type{Float64},arg1::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    err = ccall((:KSPConvergedDefaultCreate,petsc1),PetscErrorCode,(Ptr{Ptr{Void}},),arg1)
    return err
end

function KSPConvergedDefaultSetUIRNorm(arg1::KSP{Float64})
    err = ccall((:KSPConvergedDefaultSetUIRNorm,petsc1),PetscErrorCode,(KSP{Float64},),arg1)
    return err
end

function KSPConvergedDefaultSetUMIRNorm(arg1::KSP{Float64})
    err = ccall((:KSPConvergedDefaultSetUMIRNorm,petsc1),PetscErrorCode,(KSP{Float64},),arg1)
    return err
end

function KSPConvergedSkip(arg1::KSP{Float64},arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{KSPConvergedReason},StridedArray{KSPConvergedReason},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPConvergedSkip,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,Cint,Ptr{KSPConvergedReason},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
    return err
end

function KSPGetConvergedReason(arg1::KSP{Float64},arg2::Union(Ptr{KSPConvergedReason},StridedArray{KSPConvergedReason},Ptr{Void}))
    err = ccall((:KSPGetConvergedReason,petsc1),PetscErrorCode,(KSP{Float64},Ptr{KSPConvergedReason}),arg1,arg2)
    return err
end

function KSPCGSetType(arg1::KSP{Float64},arg2::KSPCGType)
    err = ccall((:KSPCGSetType,petsc1),PetscErrorCode,(KSP{Float64},KSPCGType),arg1,arg2)
    return err
end

function KSPCGUseSingleReduction(arg1::KSP{Float64},arg2::PetscBool)
    err = ccall((:KSPCGUseSingleReduction,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
    return err
end

function KSPNASHSetRadius(arg1::KSP{Float64},PetscReal::Integer)
    err = ccall((:KSPNASHSetRadius,petsc1),PetscErrorCode,(KSP{Float64},Cint),arg1,PetscReal)
    return err
end

function KSPNASHGetNormD(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:KSPNASHGetNormD,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint}),arg1,arg2)
    return err
end

function KSPNASHGetObjFcn(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:KSPNASHGetObjFcn,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint}),arg1,arg2)
    return err
end

function KSPSTCGSetRadius(arg1::KSP{Float64},PetscReal::Integer)
    err = ccall((:KSPSTCGSetRadius,petsc1),PetscErrorCode,(KSP{Float64},Cint),arg1,PetscReal)
    return err
end

function KSPSTCGGetNormD(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:KSPSTCGGetNormD,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint}),arg1,arg2)
    return err
end

function KSPSTCGGetObjFcn(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:KSPSTCGGetObjFcn,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint}),arg1,arg2)
    return err
end

function KSPGLTRSetRadius(arg1::KSP{Float64},PetscReal::Integer)
    err = ccall((:KSPGLTRSetRadius,petsc1),PetscErrorCode,(KSP{Float64},Cint),arg1,PetscReal)
    return err
end

function KSPGLTRGetNormD(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:KSPGLTRGetNormD,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint}),arg1,arg2)
    return err
end

function KSPGLTRGetObjFcn(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:KSPGLTRGetObjFcn,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint}),arg1,arg2)
    return err
end

function KSPGLTRGetMinEig(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:KSPGLTRGetMinEig,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint}),arg1,arg2)
    return err
end

function KSPGLTRGetLambda(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:KSPGLTRGetLambda,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint}),arg1,arg2)
    return err
end

function KSPPythonSetType(arg1::KSP{Float64},arg2::Union(ByteString,Symbol))
    err = ccall((:KSPPythonSetType,petsc1),PetscErrorCode,(KSP{Float64},Cstring),arg1,arg2)
    return err
end

function PCPreSolve(arg1::PC{Float64},arg2::KSP{Float64})
    err = ccall((:PCPreSolve,petsc1),PetscErrorCode,(PC{Float64},KSP{Float64}),arg1,arg2)
    return err
end

function PCPostSolve(arg1::PC{Float64},arg2::KSP{Float64})
    err = ccall((:PCPostSolve,petsc1),PetscErrorCode,(PC{Float64},KSP{Float64}),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function KSPMonitorLGResidualNormCreate(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Union(Ptr{Ptr{PetscObject}},StridedArray{Ptr{PetscObject}},Ptr{Void}))
    ccall((:KSPMonitorLGResidualNormCreate,petsc1),PetscErrorCode,(Cstring,Cstring,Cint,Cint,Cint,Cint,Ptr{Ptr{PetscObject}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function KSPMonitorLGResidualNorm(arg1::KSP{Float64},arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{PetscObject},StridedArray{PetscObject},Ptr{Void}))
    ccall((:KSPMonitorLGResidualNorm,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,Cint,Ptr{PetscObject}),arg1,arg2,PetscReal,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function KSPMonitorLGResidualNormDestroy(arg0::Type{Float64},arg1::Union(Ptr{Ptr{PetscObject}},StridedArray{Ptr{PetscObject}},Ptr{Void}))
    ccall((:KSPMonitorLGResidualNormDestroy,petsc1),PetscErrorCode,(Ptr{Ptr{PetscObject}},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function KSPMonitorLGTrueResidualNormCreate(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Union(Ptr{Ptr{PetscObject}},StridedArray{Ptr{PetscObject}},Ptr{Void}))
    ccall((:KSPMonitorLGTrueResidualNormCreate,petsc1),PetscErrorCode,(Cstring,Cstring,Cint,Cint,Cint,Cint,Ptr{Ptr{PetscObject}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function KSPMonitorLGTrueResidualNorm(arg1::KSP{Float64},arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{PetscObject},StridedArray{PetscObject},Ptr{Void}))
    ccall((:KSPMonitorLGTrueResidualNorm,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,Cint,Ptr{PetscObject}),arg1,arg2,PetscReal,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function KSPMonitorLGTrueResidualNormDestroy(arg0::Type{Float64},arg1::Union(Ptr{Ptr{PetscObject}},StridedArray{Ptr{PetscObject}},Ptr{Void}))
    ccall((:KSPMonitorLGTrueResidualNormDestroy,petsc1),PetscErrorCode,(Ptr{Ptr{PetscObject}},),arg1)
end 
=#
function KSPMonitorLGRange(arg1::KSP{Float64},arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPMonitorLGRange,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
    return err
end

function PCShellSetPreSolve(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PCShellSetPreSolve,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
    return err
end

function PCShellSetPostSolve(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PCShellSetPostSolve,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function KSPFischerGuessCreate(arg1::KSP{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{KSPFischerGuess},StridedArray{KSPFischerGuess},Ptr{Void}))
    ccall((:KSPFischerGuessCreate,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,PetscInt,Ptr{KSPFischerGuess}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function KSPFischerGuessDestroy(arg0::Type{Float64},arg1::Union(Ptr{KSPFischerGuess},StridedArray{KSPFischerGuess},Ptr{Void}))
    ccall((:KSPFischerGuessDestroy,petsc1),PetscErrorCode,(Ptr{KSPFischerGuess},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function KSPFischerGuessReset(arg0::Type{Float64},arg1::KSPFischerGuess)
    ccall((:KSPFischerGuessReset,petsc1),PetscErrorCode,(KSPFischerGuess,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function KSPFischerGuessUpdate(arg1::KSPFischerGuess,arg2::Vec{Float64})
    ccall((:KSPFischerGuessUpdate,petsc1),PetscErrorCode,(KSPFischerGuess,Vec{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function KSPFischerGuessFormGuess(arg1::KSPFischerGuess,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:KSPFischerGuessFormGuess,petsc1),PetscErrorCode,(KSPFischerGuess,Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function KSPFischerGuessSetFromOptions(arg0::Type{Float64},arg1::KSPFischerGuess)
    ccall((:KSPFischerGuessSetFromOptions,petsc1),PetscErrorCode,(KSPFischerGuess,),arg1)
end 
=#
function KSPSetUseFischerGuess(arg1::KSP{Float64},arg2::Integer,arg3::Integer)
    err = ccall((:KSPSetUseFischerGuess,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,PetscInt),arg1,arg2,arg3)
    return err
end

#= skipping function with undefined symbols: 
 function KSPSetFischerGuess(arg1::KSP{Float64},arg2::KSPFischerGuess)
    ccall((:KSPSetFischerGuess,petsc1),PetscErrorCode,(KSP{Float64},KSPFischerGuess),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function KSPGetFischerGuess(arg1::KSP{Float64},arg2::Union(Ptr{KSPFischerGuess},StridedArray{KSPFischerGuess},Ptr{Void}))
    ccall((:KSPGetFischerGuess,petsc1),PetscErrorCode,(KSP{Float64},Ptr{KSPFischerGuess}),arg1,arg2)
end 
=#
function MatCreateSchurComplement(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64},arg4::Mat{Float64},arg5::Mat{Float64},arg6::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateSchurComplement,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64},Mat{Float64},Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatSchurComplementGetKSP(arg1::Mat{Float64},arg2::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    err = ccall((:MatSchurComplementGetKSP,petsc1),PetscErrorCode,(Mat{Float64},Ptr{KSP{Float64}}),arg1,arg2)
    return err
end

function MatSchurComplementSetKSP(arg1::Mat{Float64},arg2::KSP{Float64})
    err = ccall((:MatSchurComplementSetKSP,petsc1),PetscErrorCode,(Mat{Float64},KSP{Float64}),arg1,arg2)
    return err
end

function MatSchurComplementSetSubMatrices(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64},arg4::Mat{Float64},arg5::Mat{Float64},arg6::Mat{Float64})
    err = ccall((:MatSchurComplementSetSubMatrices,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64},Mat{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatSchurComplementUpdateSubMatrices(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64},arg4::Mat{Float64},arg5::Mat{Float64},arg6::Mat{Float64})
    err = ccall((:MatSchurComplementUpdateSubMatrices,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64},Mat{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatSchurComplementGetSubMatrices(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg5::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg6::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatSchurComplementGetSubMatrices,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}},Ptr{Mat{Float64}},Ptr{Mat{Float64}},Ptr{Mat{Float64}},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

function MatSchurComplementSetAinvType(arg1::Mat{Float64},arg2::MatSchurComplementAinvType)
    err = ccall((:MatSchurComplementSetAinvType,petsc1),PetscErrorCode,(Mat{Float64},MatSchurComplementAinvType),arg1,arg2)
    return err
end

function MatSchurComplementGetAinvType(arg1::Mat{Float64},arg2::Union(Ptr{MatSchurComplementAinvType},StridedArray{MatSchurComplementAinvType},Ptr{Void}))
    err = ccall((:MatSchurComplementGetAinvType,petsc1),PetscErrorCode,(Mat{Float64},Ptr{MatSchurComplementAinvType}),arg1,arg2)
    return err
end

function MatSchurComplementGetPmat(arg1::Mat{Float64},arg2::MatReuse,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatSchurComplementGetPmat,petsc1),PetscErrorCode,(Mat{Float64},MatReuse,Ptr{Mat{Float64}}),arg1,arg2,arg3)
    return err
end

function MatSchurComplementComputeExplicitOperator(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatSchurComplementComputeExplicitOperator,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2)
    return err
end

function MatGetSchurComplement(arg1::Mat{Float64},arg2::IS{Float64},arg3::IS{Float64},arg4::IS{Float64},arg5::IS{Float64},arg6::MatReuse,arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg8::MatSchurComplementAinvType,arg9::MatReuse,arg10::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatGetSchurComplement,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},IS{Float64},IS{Float64},IS{Float64},MatReuse,Ptr{Mat{Float64}},MatSchurComplementAinvType,MatReuse,Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
    return err
end

function MatCreateSchurComplementPmat(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64},arg4::Mat{Float64},arg5::MatSchurComplementAinvType,arg6::MatReuse,arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    err = ccall((:MatCreateSchurComplementPmat,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64},Mat{Float64},MatSchurComplementAinvType,MatReuse,Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
    return err
end

#= skipping function with undefined symbols: 
 function KSPSetDM(arg1::KSP{Float64},arg2::DM)
    ccall((:KSPSetDM,petsc1),PetscErrorCode,(KSP{Float64},DM),arg1,arg2)
end 
=#
function KSPSetDMActive(arg1::KSP{Float64},arg2::PetscBool)
    err = ccall((:KSPSetDMActive,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function KSPGetDM(arg1::KSP{Float64},arg2::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:KSPGetDM,petsc1),PetscErrorCode,(KSP{Float64},Ptr{DM}),arg1,arg2)
end 
=#
function KSPSetApplicationContext(arg1::KSP{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPSetApplicationContext,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Void}),arg1,arg2)
    return err
end

function KSPGetApplicationContext(arg1::KSP{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPGetApplicationContext,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Void}),arg1,arg2)
    return err
end

function KSPSetComputeRHS(arg1::KSP{Float64},func::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPSetComputeRHS,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Void},Ptr{Void}),arg1,func,arg2)
    return err
end

function KSPSetComputeOperators(arg1::KSP{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPSetComputeOperators,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
    return err
end

function KSPSetComputeInitialGuess(arg1::KSP{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPSetComputeInitialGuess,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
    return err
end

#= skipping function with undefined symbols: 
 function DMKSPSetComputeOperators(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMKSPSetComputeOperators,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMKSPGetComputeOperators(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMKSPGetComputeOperators,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMKSPSetComputeRHS(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMKSPSetComputeRHS,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMKSPGetComputeRHS(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMKSPGetComputeRHS,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMKSPSetComputeInitialGuess(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMKSPSetComputeInitialGuess,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMKSPGetComputeInitialGuess(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMKSPGetComputeInitialGuess,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMGlobalToLocalSolve(arg1::DM,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:DMGlobalToLocalSolve,petsc1),PetscErrorCode,(DM,Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexProjectField(arg1::DM,arg2::Vec{Float64},arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg4::InsertMode,arg5::Vec{Float64})
    ccall((:DMPlexProjectField,petsc1),PetscErrorCode,(DM,Vec{Float64},Ptr{Ptr{Void}},InsertMode,Vec{Float64}),arg1,arg2,arg3,arg4,arg5)
end 
=#
