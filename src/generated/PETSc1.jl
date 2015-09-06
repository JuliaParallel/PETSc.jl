# Julia wrapper for header: ../../deps/petsc-3.6.0/include/petsc.h
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

#= skipping function with undefined symbols: 
 function PetscBagCreate(arg0::Type{Float64},arg1::MPI_Comm,size_t::Integer,arg2::Union(Ptr{PetscBag},StridedArray{PetscBag},Ptr{Void}))
    ccall((:PetscBagCreate,petsc1),PetscErrorCode,(comm_type,Cint,Ptr{PetscBag}),arg1.val,size_t,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscBagDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscBag},StridedArray{PetscBag},Ptr{Void}))
    ccall((:PetscBagDestroy,petsc1),PetscErrorCode,(Ptr{PetscBag},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscBagGetData(arg0::Type{Float64},arg1::PetscBag,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:PetscBagGetData,petsc1),PetscErrorCode,(PetscBag,Ptr{Ptr{Void}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscBagRegisterReal(arg0::Type{Float64},arg1::PetscBag,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),PetscReal::Integer,arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol))
    ccall((:PetscBagRegisterReal,petsc1),PetscErrorCode,(PetscBag,Ptr{Void},Cint,Cstring,Cstring),arg1,arg2,PetscReal,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscBagRegisterRealArray(arg0::Type{Float64},arg1::PetscBag,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Integer,arg4::Union(ByteString,Symbol),arg5::Union(ByteString,Symbol))
    ccall((:PetscBagRegisterRealArray,petsc1),PetscErrorCode,(PetscBag,Ptr{Void},PetscInt,Cstring,Cstring),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscBagRegisterString(arg0::Type{Float64},arg1::PetscBag,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Integer,arg4::Union(ByteString,Symbol),arg5::Union(ByteString,Symbol),arg6::Union(ByteString,Symbol))
    ccall((:PetscBagRegisterString,petsc1),PetscErrorCode,(PetscBag,Ptr{Void},PetscInt,Cstring,Cstring,Cstring),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function PetscBagRegisterScalar(arg0::Type{Float64},arg1::PetscBag,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Float64,arg4::Union(ByteString,Symbol),arg5::Union(ByteString,Symbol))
    ccall((:PetscBagRegisterScalar,petsc1),PetscErrorCode,(PetscBag,Ptr{Void},Float64,Cstring,Cstring),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscBagRegisterInt(arg0::Type{Float64},arg1::PetscBag,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Integer,arg4::Union(ByteString,Symbol),arg5::Union(ByteString,Symbol))
    ccall((:PetscBagRegisterInt,petsc1),PetscErrorCode,(PetscBag,Ptr{Void},PetscInt,Cstring,Cstring),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscBagRegister64bitInt(arg0::Type{Float64},arg1::PetscBag,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Petsc64bitInt,arg4::Union(ByteString,Symbol),arg5::Union(ByteString,Symbol))
    ccall((:PetscBagRegister64bitInt,petsc1),PetscErrorCode,(PetscBag,Ptr{Void},Petsc64bitInt,Cstring,Cstring),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscBagRegisterIntArray(arg0::Type{Float64},arg1::PetscBag,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Integer,arg4::Union(ByteString,Symbol),arg5::Union(ByteString,Symbol))
    ccall((:PetscBagRegisterIntArray,petsc1),PetscErrorCode,(PetscBag,Ptr{Void},PetscInt,Cstring,Cstring),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscBagRegisterEnum(arg0::Type{Float64},arg1::PetscBag,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}),arg4::PetscEnum,arg5::Union(ByteString,Symbol),arg6::Union(ByteString,Symbol))
    ccall((:PetscBagRegisterEnum,petsc1),PetscErrorCode,(PetscBag,Ptr{Void},Ptr{Ptr{Uint8}},PetscEnum,Cstring,Cstring),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function PetscBagRegisterBool(arg0::Type{Float64},arg1::PetscBag,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::PetscBool,arg4::Union(ByteString,Symbol),arg5::Union(ByteString,Symbol))
    ccall((:PetscBagRegisterBool,petsc1),PetscErrorCode,(PetscBag,Ptr{Void},PetscBool,Cstring,Cstring),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscBagRegisterBoolArray(arg0::Type{Float64},arg1::PetscBag,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Integer,arg4::Union(ByteString,Symbol),arg5::Union(ByteString,Symbol))
    ccall((:PetscBagRegisterBoolArray,petsc1),PetscErrorCode,(PetscBag,Ptr{Void},PetscInt,Cstring,Cstring),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscBagGetNames(arg0::Type{Float64},arg1::PetscBag,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:PetscBagGetNames,petsc1),PetscErrorCode,(PetscBag,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscBagSetFromOptions(arg0::Type{Float64},arg1::PetscBag)
    ccall((:PetscBagSetFromOptions,petsc1),PetscErrorCode,(PetscBag,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscBagGetName(arg0::Type{Float64},arg1::PetscBag,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:PetscBagGetName,petsc1),PetscErrorCode,(PetscBag,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscBagSetName(arg0::Type{Float64},arg1::PetscBag,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol))
    ccall((:PetscBagSetName,petsc1),PetscErrorCode,(PetscBag,Cstring,Cstring),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscBagSetOptionsPrefix(arg0::Type{Float64},arg1::PetscBag,arg2::Union(ByteString,Symbol))
    ccall((:PetscBagSetOptionsPrefix,petsc1),PetscErrorCode,(PetscBag,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscBagView(arg1::PetscBag,arg2::PetscViewer{Float64})
    ccall((:PetscBagView,petsc1),PetscErrorCode,(PetscBag,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscBagLoad(arg1::PetscViewer{Float64},arg2::PetscBag)
    ccall((:PetscBagLoad,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscBag),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscBagSetViewer(arg0::Type{Float64},arg1::PetscBag,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscBagSetViewer,petsc1),PetscErrorCode,(PetscBag,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscBagSetLoader(arg0::Type{Float64},arg1::PetscBag,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscBagSetLoader,petsc1),PetscErrorCode,(PetscBag,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscBagSetDestroy(arg0::Type{Float64},arg1::PetscBag,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscBagSetDestroy,petsc1),PetscErrorCode,(PetscBag,Ptr{Void}),arg1,arg2)
end 
=#
function PetscGetCPUTime(arg0::Type{Float64},arg1::Union(Ptr{PetscLogDouble},StridedArray{PetscLogDouble},Ptr{Void}))
    err = ccall((:PetscGetCPUTime,petsc1),PetscErrorCode,(Ptr{PetscLogDouble},),arg1)
    return err
end

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
#= skipping function with undefined symbols: 
 function PetscTableCreate(arg0::Type{Float64},arg1::Integer,arg2::Integer,arg3::Union(Ptr{PetscTable},StridedArray{PetscTable},Ptr{Void}))
    ccall((:PetscTableCreate,petsc1),PetscErrorCode,(PetscInt,PetscInt,Ptr{PetscTable}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscTableCreateCopy(arg0::Type{Float64},arg1::PetscTable,arg2::Union(Ptr{PetscTable},StridedArray{PetscTable},Ptr{Void}))
    ccall((:PetscTableCreateCopy,petsc1),PetscErrorCode,(PetscTable,Ptr{PetscTable}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscTableDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscTable},StridedArray{PetscTable},Ptr{Void}))
    ccall((:PetscTableDestroy,petsc1),PetscErrorCode,(Ptr{PetscTable},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscTableGetCount(arg0::Type{Float64},arg1::PetscTable,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscTableGetCount,petsc1),PetscErrorCode,(PetscTable,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscTableIsEmpty(arg0::Type{Float64},arg1::PetscTable,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscTableIsEmpty,petsc1),PetscErrorCode,(PetscTable,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscTableAddExpand(arg0::Type{Float64},arg1::PetscTable,arg2::Integer,arg3::Integer,arg4::InsertMode)
    ccall((:PetscTableAddExpand,petsc1),PetscErrorCode,(PetscTable,PetscInt,PetscInt,InsertMode),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscTableAddCountExpand(arg0::Type{Float64},arg1::PetscTable,arg2::Integer)
    ccall((:PetscTableAddCountExpand,petsc1),PetscErrorCode,(PetscTable,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscTableGetHeadPosition(arg0::Type{Float64},arg1::PetscTable,arg2::Union(Ptr{PetscTablePosition},StridedArray{PetscTablePosition},Ptr{Void}))
    ccall((:PetscTableGetHeadPosition,petsc1),PetscErrorCode,(PetscTable,Ptr{PetscTablePosition}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscTableGetNext(arg0::Type{Float64},arg1::PetscTable,arg2::Union(Ptr{PetscTablePosition},StridedArray{PetscTablePosition},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscTableGetNext,petsc1),PetscErrorCode,(PetscTable,Ptr{PetscTablePosition},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscTableRemoveAll(arg0::Type{Float64},arg1::PetscTable)
    ccall((:PetscTableRemoveAll,petsc1),PetscErrorCode,(PetscTable,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscMatlabEngineCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{PetscMatlabEngine},StridedArray{PetscMatlabEngine},Ptr{Void}))
    ccall((:PetscMatlabEngineCreate,petsc1),PetscErrorCode,(comm_type,Cstring,Ptr{PetscMatlabEngine}),arg1.val,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscMatlabEngineDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscMatlabEngine},StridedArray{PetscMatlabEngine},Ptr{Void}))
    ccall((:PetscMatlabEngineDestroy,petsc1),PetscErrorCode,(Ptr{PetscMatlabEngine},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscMatlabEngineGetOutput(arg0::Type{Float64},arg1::PetscMatlabEngine,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:PetscMatlabEngineGetOutput,petsc1),PetscErrorCode,(PetscMatlabEngine,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscMatlabEnginePrintOutput(arg0::Type{Float64},arg1::PetscMatlabEngine,arg2::Union(Ptr{FILE},StridedArray{FILE},Ptr{Void}))
    ccall((:PetscMatlabEnginePrintOutput,petsc1),PetscErrorCode,(PetscMatlabEngine,Ptr{FILE}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscMatlabEnginePut(arg0::Type{Float64},arg1::PetscMatlabEngine,arg2::PetscObject)
    ccall((:PetscMatlabEnginePut,petsc1),PetscErrorCode,(PetscMatlabEngine,PetscObject),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscMatlabEngineGet(arg0::Type{Float64},arg1::PetscMatlabEngine,arg2::PetscObject)
    ccall((:PetscMatlabEngineGet,petsc1),PetscErrorCode,(PetscMatlabEngine,PetscObject),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscMatlabEnginePutArray(arg0::Type{Float64},arg1::PetscMatlabEngine,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::Union(ByteString,Symbol))
    ccall((:PetscMatlabEnginePutArray,petsc1),PetscErrorCode,(PetscMatlabEngine,Cint,Cint,Ptr{Float64},Cstring),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscMatlabEngineGetArray(arg0::Type{Float64},arg1::PetscMatlabEngine,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::Union(ByteString,Symbol))
    ccall((:PetscMatlabEngineGetArray,petsc1),PetscErrorCode,(PetscMatlabEngine,Cint,Cint,Ptr{Float64},Cstring),arg1,arg2,arg3,arg4,arg5)
end 
=#
function PETSC_MATLAB_ENGINE_(arg0::Type{Float64},arg1::MPI_Comm)
    err = ccall((:PETSC_MATLAB_ENGINE_,petsc1),PetscMatlabEngine,(comm_type,),arg1.val)
    return err
end

function PetscDrawInitializePackage(arg0::Type{Float64})
    err = ccall((:PetscDrawInitializePackage,petsc1),PetscErrorCode,())
    return err
end

function PetscDrawRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscDrawRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function PetscDrawGetType(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(Ptr{PetscDrawType},StridedArray{PetscDrawType},Ptr{Void}))
    ccall((:PetscDrawGetType,petsc1),PetscErrorCode,(PetscDraw,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSetType(arg0::Type{Float64},arg1::PetscDraw,arg2::PetscDrawType)
    ccall((:PetscDrawSetType,petsc1),PetscErrorCode,(PetscDraw,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Union(Ptr{PetscDraw},StridedArray{PetscDraw},Ptr{Void}))
    ccall((:PetscDrawCreate,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint,Cint,Cint,Cint,Ptr{PetscDraw}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSetFromOptions(arg0::Type{Float64},arg1::PetscDraw)
    ccall((:PetscDrawSetFromOptions,petsc1),PetscErrorCode,(PetscDraw,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSetSave(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(ByteString,Symbol),arg3::PetscBool)
    ccall((:PetscDrawSetSave,petsc1),PetscErrorCode,(PetscDraw,Cstring,PetscBool),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSetSaveFinalImage(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(ByteString,Symbol))
    ccall((:PetscDrawSetSaveFinalImage,petsc1),PetscErrorCode,(PetscDraw,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawView(arg1::PetscDraw,arg2::PetscViewer{Float64})
    ccall((:PetscDrawView,petsc1),PetscErrorCode,(PetscDraw,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawOpenGLUT(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Union(Ptr{PetscDraw},StridedArray{PetscDraw},Ptr{Void}))
    ccall((:PetscDrawOpenGLUT,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint,Cint,Cint,Cint,Ptr{PetscDraw}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawOpenNull(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscDraw},StridedArray{PetscDraw},Ptr{Void}))
    ccall((:PetscDrawOpenNull,petsc1),PetscErrorCode,(comm_type,Ptr{PetscDraw}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscDraw},StridedArray{PetscDraw},Ptr{Void}))
    ccall((:PetscDrawDestroy,petsc1),PetscErrorCode,(Ptr{PetscDraw},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawIsNull(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscDrawIsNull,petsc1),PetscErrorCode,(PetscDraw,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawGetPopup(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(Ptr{PetscDraw},StridedArray{PetscDraw},Ptr{Void}))
    ccall((:PetscDrawGetPopup,petsc1),PetscErrorCode,(PetscDraw,Ptr{PetscDraw}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawCheckResizedWindow(arg0::Type{Float64},arg1::PetscDraw)
    ccall((:PetscDrawCheckResizedWindow,petsc1),PetscErrorCode,(PetscDraw,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawResizeWindow(arg0::Type{Float64},arg1::PetscDraw,arg2::Integer,arg3::Integer)
    ccall((:PetscDrawResizeWindow,petsc1),PetscErrorCode,(PetscDraw,Cint,Cint),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawScalePopup(arg0::Type{Float64},arg1::PetscDraw,PetscReal::Integer,arg2::Integer)
    ccall((:PetscDrawScalePopup,petsc1),PetscErrorCode,(PetscDraw,Cint,Cint),arg1,PetscReal,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawPixelToCoordinate(arg0::Type{Float64},arg1::PetscDraw,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscDrawPixelToCoordinate,petsc1),PetscErrorCode,(PetscDraw,PetscInt,PetscInt,Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawCoordinateToPixel(arg0::Type{Float64},arg1::PetscDraw,PetscReal::Integer,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscDrawCoordinateToPixel,petsc1),PetscErrorCode,(PetscDraw,Cint,Cint,Ptr{PetscInt},Ptr{PetscInt}),arg1,PetscReal,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawIndicatorFunction(arg0::Type{Float64},arg1::PetscDraw,PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg7::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscDrawIndicatorFunction,petsc1),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Cint,Cint,Ptr{Void},Ptr{Void}),arg1,PetscReal,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawLine(arg0::Type{Float64},arg1::PetscDraw,PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer)
    ccall((:PetscDrawLine,petsc1),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawArrow(arg0::Type{Float64},arg1::PetscDraw,PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer)
    ccall((:PetscDrawArrow,petsc1),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawLineSetWidth(arg0::Type{Float64},arg1::PetscDraw,PetscReal::Integer)
    ccall((:PetscDrawLineSetWidth,petsc1),PetscErrorCode,(PetscDraw,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawLineGetWidth(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscDrawLineGetWidth,petsc1),PetscErrorCode,(PetscDraw,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawMarker(arg0::Type{Float64},arg1::PetscDraw,PetscReal::Integer,arg2::Integer,arg3::Integer)
    ccall((:PetscDrawMarker,petsc1),PetscErrorCode,(PetscDraw,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSetMarkerType(arg0::Type{Float64},arg1::PetscDraw,arg2::PetscDrawMarkerType)
    ccall((:PetscDrawSetMarkerType,petsc1),PetscErrorCode,(PetscDraw,PetscDrawMarkerType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawGetMarkerType(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(Ptr{PetscDrawMarkerType},StridedArray{PetscDrawMarkerType},Ptr{Void}))
    ccall((:PetscDrawGetMarkerType,petsc1),PetscErrorCode,(PetscDraw,Ptr{PetscDrawMarkerType}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawPoint(arg0::Type{Float64},arg1::PetscDraw,PetscReal::Integer,arg2::Integer,arg3::Integer)
    ccall((:PetscDrawPoint,petsc1),PetscErrorCode,(PetscDraw,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawPointPixel(arg0::Type{Float64},arg1::PetscDraw,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:PetscDrawPointPixel,petsc1),PetscErrorCode,(PetscDraw,PetscInt,PetscInt,Cint),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawPointSetSize(arg0::Type{Float64},arg1::PetscDraw,PetscReal::Integer)
    ccall((:PetscDrawPointSetSize,petsc1),PetscErrorCode,(PetscDraw,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawRectangle(arg0::Type{Float64},arg1::PetscDraw,PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Integer)
    ccall((:PetscDrawRectangle,petsc1),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Cint,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawTriangle(arg0::Type{Float64},arg1::PetscDraw,PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Integer,arg9::Integer)
    ccall((:PetscDrawTriangle,petsc1),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Cint,Cint,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawEllipse(arg0::Type{Float64},arg1::PetscDraw,PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer)
    ccall((:PetscDrawEllipse,petsc1),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawTensorContourPatch(arg0::Type{Float64},arg1::PetscDraw,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),PetscReal::Integer,arg6::Integer,arg7::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscDrawTensorContourPatch,petsc1),PetscErrorCode,(PetscDraw,Cint,Cint,Ptr{Cint},Ptr{Cint},Cint,Cint,Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,PetscReal,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawTensorContour(arg0::Type{Float64},arg1::PetscDraw,arg2::Integer,arg3::Integer,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscDrawTensorContour,petsc1),PetscErrorCode,(PetscDraw,Cint,Cint,Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,PetscReal,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawString(arg0::Type{Float64},arg1::PetscDraw,PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Union(ByteString,Symbol))
    ccall((:PetscDrawString,petsc1),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Cstring),arg1,PetscReal,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawStringCentered(arg0::Type{Float64},arg1::PetscDraw,PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Union(ByteString,Symbol))
    ccall((:PetscDrawStringCentered,petsc1),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Cstring),arg1,PetscReal,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawStringBoxed(arg0::Type{Float64},arg1::PetscDraw,PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(ByteString,Symbol),arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg7::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscDrawStringBoxed,petsc1),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Cint,Cstring,Ptr{Cint},Ptr{Cint}),arg1,PetscReal,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawStringBoxedSize(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscDrawStringBoxedSize,petsc1),PetscErrorCode,(PetscDraw,Cstring,Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawStringVertical(arg0::Type{Float64},arg1::PetscDraw,PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Union(ByteString,Symbol))
    ccall((:PetscDrawStringVertical,petsc1),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Cstring),arg1,PetscReal,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawStringSetSize(arg0::Type{Float64},arg1::PetscDraw,PetscReal::Integer,arg2::Integer)
    ccall((:PetscDrawStringSetSize,petsc1),PetscErrorCode,(PetscDraw,Cint,Cint),arg1,PetscReal,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawStringGetSize(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscDrawStringGetSize,petsc1),PetscErrorCode,(PetscDraw,Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSetViewPort(arg0::Type{Float64},arg1::PetscDraw,PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:PetscDrawSetViewPort,petsc1),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawGetViewPort(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscDrawGetViewPort,petsc1),PetscErrorCode,(PetscDraw,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSplitViewPort(arg0::Type{Float64},arg1::PetscDraw)
    ccall((:PetscDrawSplitViewPort,petsc1),PetscErrorCode,(PetscDraw,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSetCoordinates(arg0::Type{Float64},arg1::PetscDraw,PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:PetscDrawSetCoordinates,petsc1),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawGetCoordinates(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscDrawGetCoordinates,petsc1),PetscErrorCode,(PetscDraw,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSetTitle(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(ByteString,Symbol))
    ccall((:PetscDrawSetTitle,petsc1),PetscErrorCode,(PetscDraw,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawAppendTitle(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(ByteString,Symbol))
    ccall((:PetscDrawAppendTitle,petsc1),PetscErrorCode,(PetscDraw,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawGetTitle(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:PetscDrawGetTitle,petsc1),PetscErrorCode,(PetscDraw,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSetPause(arg0::Type{Float64},arg1::PetscDraw,PetscReal::Integer)
    ccall((:PetscDrawSetPause,petsc1),PetscErrorCode,(PetscDraw,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawGetPause(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscDrawGetPause,petsc1),PetscErrorCode,(PetscDraw,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawPause(arg0::Type{Float64},arg1::PetscDraw)
    ccall((:PetscDrawPause,petsc1),PetscErrorCode,(PetscDraw,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSetDoubleBuffer(arg0::Type{Float64},arg1::PetscDraw)
    ccall((:PetscDrawSetDoubleBuffer,petsc1),PetscErrorCode,(PetscDraw,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawFlush(arg0::Type{Float64},arg1::PetscDraw)
    ccall((:PetscDrawFlush,petsc1),PetscErrorCode,(PetscDraw,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSynchronizedFlush(arg0::Type{Float64},arg1::PetscDraw)
    ccall((:PetscDrawSynchronizedFlush,petsc1),PetscErrorCode,(PetscDraw,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawClear(arg0::Type{Float64},arg1::PetscDraw)
    ccall((:PetscDrawClear,petsc1),PetscErrorCode,(PetscDraw,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSave(arg0::Type{Float64},arg1::PetscDraw)
    ccall((:PetscDrawSave,petsc1),PetscErrorCode,(PetscDraw,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSynchronizedClear(arg0::Type{Float64},arg1::PetscDraw)
    ccall((:PetscDrawSynchronizedClear,petsc1),PetscErrorCode,(PetscDraw,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawBOP(arg0::Type{Float64},arg1::PetscDraw)
    ccall((:PetscDrawBOP,petsc1),PetscErrorCode,(PetscDraw,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawEOP(arg0::Type{Float64},arg1::PetscDraw)
    ccall((:PetscDrawEOP,petsc1),PetscErrorCode,(PetscDraw,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSetDisplay(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(ByteString,Symbol))
    ccall((:PetscDrawSetDisplay,petsc1),PetscErrorCode,(PetscDraw,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawGetSingleton(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(Ptr{PetscDraw},StridedArray{PetscDraw},Ptr{Void}))
    ccall((:PetscDrawGetSingleton,petsc1),PetscErrorCode,(PetscDraw,Ptr{PetscDraw}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawRestoreSingleton(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(Ptr{PetscDraw},StridedArray{PetscDraw},Ptr{Void}))
    ccall((:PetscDrawRestoreSingleton,petsc1),PetscErrorCode,(PetscDraw,Ptr{PetscDraw}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawGetCurrentPoint(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscDrawGetCurrentPoint,petsc1),PetscErrorCode,(PetscDraw,Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSetCurrentPoint(arg0::Type{Float64},arg1::PetscDraw,PetscReal::Integer,arg2::Integer)
    ccall((:PetscDrawSetCurrentPoint,petsc1),PetscErrorCode,(PetscDraw,Cint,Cint),arg1,PetscReal,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawPushCurrentPoint(arg0::Type{Float64},arg1::PetscDraw,PetscReal::Integer,arg2::Integer)
    ccall((:PetscDrawPushCurrentPoint,petsc1),PetscErrorCode,(PetscDraw,Cint,Cint),arg1,PetscReal,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawPopCurrentPoint(arg0::Type{Float64},arg1::PetscDraw)
    ccall((:PetscDrawPopCurrentPoint,petsc1),PetscErrorCode,(PetscDraw,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawGetBoundingBox(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscDrawGetBoundingBox,petsc1),PetscErrorCode,(PetscDraw,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawGetMouseButton(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(Ptr{PetscDrawButton},StridedArray{PetscDrawButton},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscDrawGetMouseButton,petsc1),PetscErrorCode,(PetscDraw,Ptr{PetscDrawButton},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSynchronizedGetMouseButton(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(Ptr{PetscDrawButton},StridedArray{PetscDrawButton},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscDrawSynchronizedGetMouseButton,petsc1),PetscErrorCode,(PetscDraw,Ptr{PetscDrawButton},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawZoom(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscDrawZoom,petsc1),PetscErrorCode,(PetscDraw,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawViewPortsCreate(arg0::Type{Float64},arg1::PetscDraw,arg2::Integer,arg3::Union(Ptr{Ptr{PetscDrawViewPorts}},StridedArray{Ptr{PetscDrawViewPorts}},Ptr{Void}))
    ccall((:PetscDrawViewPortsCreate,petsc1),PetscErrorCode,(PetscDraw,PetscInt,Ptr{Ptr{PetscDrawViewPorts}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawViewPortsCreateRect(arg0::Type{Float64},arg1::PetscDraw,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Ptr{PetscDrawViewPorts}},StridedArray{Ptr{PetscDrawViewPorts}},Ptr{Void}))
    ccall((:PetscDrawViewPortsCreateRect,petsc1),PetscErrorCode,(PetscDraw,PetscInt,PetscInt,Ptr{Ptr{PetscDrawViewPorts}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawViewPortsDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscDrawViewPorts},StridedArray{PetscDrawViewPorts},Ptr{Void}))
    ccall((:PetscDrawViewPortsDestroy,petsc1),PetscErrorCode,(Ptr{PetscDrawViewPorts},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawViewPortsSet(arg0::Type{Float64},arg1::Union(Ptr{PetscDrawViewPorts},StridedArray{PetscDrawViewPorts},Ptr{Void}),arg2::Integer)
    ccall((:PetscDrawViewPortsSet,petsc1),PetscErrorCode,(Ptr{PetscDrawViewPorts},PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawAxisCreate(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(Ptr{PetscDrawAxis},StridedArray{PetscDrawAxis},Ptr{Void}))
    ccall((:PetscDrawAxisCreate,petsc1),PetscErrorCode,(PetscDraw,Ptr{PetscDrawAxis}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawAxisDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscDrawAxis},StridedArray{PetscDrawAxis},Ptr{Void}))
    ccall((:PetscDrawAxisDestroy,petsc1),PetscErrorCode,(Ptr{PetscDrawAxis},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawAxisDraw(arg0::Type{Float64},arg1::PetscDrawAxis)
    ccall((:PetscDrawAxisDraw,petsc1),PetscErrorCode,(PetscDrawAxis,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawAxisSetLimits(arg0::Type{Float64},arg1::PetscDrawAxis,PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:PetscDrawAxisSetLimits,petsc1),PetscErrorCode,(PetscDrawAxis,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawAxisGetLimits(arg0::Type{Float64},arg1::PetscDrawAxis,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscDrawAxisGetLimits,petsc1),PetscErrorCode,(PetscDrawAxis,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawAxisSetHoldLimits(arg0::Type{Float64},arg1::PetscDrawAxis,arg2::PetscBool)
    ccall((:PetscDrawAxisSetHoldLimits,petsc1),PetscErrorCode,(PetscDrawAxis,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawAxisSetColors(arg0::Type{Float64},arg1::PetscDrawAxis,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:PetscDrawAxisSetColors,petsc1),PetscErrorCode,(PetscDrawAxis,Cint,Cint,Cint),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawAxisSetLabels(arg0::Type{Float64},arg1::PetscDrawAxis,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol))
    ccall((:PetscDrawAxisSetLabels,petsc1),PetscErrorCode,(PetscDrawAxis,Cstring,Cstring,Cstring),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawLGCreate(arg0::Type{Float64},arg1::PetscDraw,arg2::Integer,arg3::Union(Ptr{PetscDrawLG},StridedArray{PetscDrawLG},Ptr{Void}))
    ccall((:PetscDrawLGCreate,petsc1),PetscErrorCode,(PetscDraw,PetscInt,Ptr{PetscDrawLG}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawLGDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscDrawLG},StridedArray{PetscDrawLG},Ptr{Void}))
    ccall((:PetscDrawLGDestroy,petsc1),PetscErrorCode,(Ptr{PetscDrawLG},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawLGAddPoint(arg0::Type{Float64},arg1::PetscDrawLG,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscDrawLGAddPoint,petsc1),PetscErrorCode,(PetscDrawLG,Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawLGAddCommonPoint(arg0::Type{Float64},arg1::PetscDrawLG,PetscReal::Integer,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscDrawLGAddCommonPoint,petsc1),PetscErrorCode,(PetscDrawLG,Cint,Ptr{Cint}),arg1,PetscReal,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawLGAddPoints(arg0::Type{Float64},arg1::PetscDrawLG,arg2::Integer,arg3::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg4::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}))
    ccall((:PetscDrawLGAddPoints,petsc1),PetscErrorCode,(PetscDrawLG,PetscInt,Ptr{Ptr{Cint}},Ptr{Ptr{Cint}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawLGDraw(arg0::Type{Float64},arg1::PetscDrawLG)
    ccall((:PetscDrawLGDraw,petsc1),PetscErrorCode,(PetscDrawLG,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawLGView(arg1::PetscDrawLG,arg2::PetscViewer{Float64})
    ccall((:PetscDrawLGView,petsc1),PetscErrorCode,(PetscDrawLG,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawLGReset(arg0::Type{Float64},arg1::PetscDrawLG)
    ccall((:PetscDrawLGReset,petsc1),PetscErrorCode,(PetscDrawLG,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawLGSetDimension(arg0::Type{Float64},arg1::PetscDrawLG,arg2::Integer)
    ccall((:PetscDrawLGSetDimension,petsc1),PetscErrorCode,(PetscDrawLG,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawLGGetDimension(arg0::Type{Float64},arg1::PetscDrawLG,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscDrawLGGetDimension,petsc1),PetscErrorCode,(PetscDrawLG,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawLGSetLegend(arg0::Type{Float64},arg1::PetscDrawLG,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:PetscDrawLGSetLegend,petsc1),PetscErrorCode,(PetscDrawLG,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawLGGetAxis(arg0::Type{Float64},arg1::PetscDrawLG,arg2::Union(Ptr{PetscDrawAxis},StridedArray{PetscDrawAxis},Ptr{Void}))
    ccall((:PetscDrawLGGetAxis,petsc1),PetscErrorCode,(PetscDrawLG,Ptr{PetscDrawAxis}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawLGGetDraw(arg0::Type{Float64},arg1::PetscDrawLG,arg2::Union(Ptr{PetscDraw},StridedArray{PetscDraw},Ptr{Void}))
    ccall((:PetscDrawLGGetDraw,petsc1),PetscErrorCode,(PetscDrawLG,Ptr{PetscDraw}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawLGSetUseMarkers(arg0::Type{Float64},arg1::PetscDrawLG,arg2::PetscBool)
    ccall((:PetscDrawLGSetUseMarkers,petsc1),PetscErrorCode,(PetscDrawLG,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawLGSetLimits(arg0::Type{Float64},arg1::PetscDrawLG,PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:PetscDrawLGSetLimits,petsc1),PetscErrorCode,(PetscDrawLG,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawLGSetColors(arg0::Type{Float64},arg1::PetscDrawLG,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscDrawLGSetColors,petsc1),PetscErrorCode,(PetscDrawLG,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawLGSetFromOptions(arg0::Type{Float64},arg1::PetscDrawLG)
    ccall((:PetscDrawLGSetFromOptions,petsc1),PetscErrorCode,(PetscDrawLG,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSPCreate(arg0::Type{Float64},arg1::PetscDraw,arg2::Integer,arg3::Union(Ptr{PetscDrawSP},StridedArray{PetscDrawSP},Ptr{Void}))
    ccall((:PetscDrawSPCreate,petsc1),PetscErrorCode,(PetscDraw,Cint,Ptr{PetscDrawSP}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSPDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscDrawSP},StridedArray{PetscDrawSP},Ptr{Void}))
    ccall((:PetscDrawSPDestroy,petsc1),PetscErrorCode,(Ptr{PetscDrawSP},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSPAddPoint(arg0::Type{Float64},arg1::PetscDrawSP,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscDrawSPAddPoint,petsc1),PetscErrorCode,(PetscDrawSP,Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSPAddPoints(arg0::Type{Float64},arg1::PetscDrawSP,arg2::Integer,arg3::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg4::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}))
    ccall((:PetscDrawSPAddPoints,petsc1),PetscErrorCode,(PetscDrawSP,Cint,Ptr{Ptr{Cint}},Ptr{Ptr{Cint}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSPDraw(arg0::Type{Float64},arg1::PetscDrawSP,arg2::PetscBool)
    ccall((:PetscDrawSPDraw,petsc1),PetscErrorCode,(PetscDrawSP,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSPReset(arg0::Type{Float64},arg1::PetscDrawSP)
    ccall((:PetscDrawSPReset,petsc1),PetscErrorCode,(PetscDrawSP,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSPSetDimension(arg0::Type{Float64},arg1::PetscDrawSP,arg2::Integer)
    ccall((:PetscDrawSPSetDimension,petsc1),PetscErrorCode,(PetscDrawSP,Cint),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSPGetAxis(arg0::Type{Float64},arg1::PetscDrawSP,arg2::Union(Ptr{PetscDrawAxis},StridedArray{PetscDrawAxis},Ptr{Void}))
    ccall((:PetscDrawSPGetAxis,petsc1),PetscErrorCode,(PetscDrawSP,Ptr{PetscDrawAxis}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSPGetDraw(arg0::Type{Float64},arg1::PetscDrawSP,arg2::Union(Ptr{PetscDraw},StridedArray{PetscDraw},Ptr{Void}))
    ccall((:PetscDrawSPGetDraw,petsc1),PetscErrorCode,(PetscDrawSP,Ptr{PetscDraw}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawSPSetLimits(arg0::Type{Float64},arg1::PetscDrawSP,PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:PetscDrawSPSetLimits,petsc1),PetscErrorCode,(PetscDrawSP,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawLGSPDraw(arg0::Type{Float64},arg1::PetscDrawLG,arg2::PetscDrawSP)
    ccall((:PetscDrawLGSPDraw,petsc1),PetscErrorCode,(PetscDrawLG,PetscDrawSP),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawHGCreate(arg0::Type{Float64},arg1::PetscDraw,arg2::Integer,arg3::Union(Ptr{PetscDrawHG},StridedArray{PetscDrawHG},Ptr{Void}))
    ccall((:PetscDrawHGCreate,petsc1),PetscErrorCode,(PetscDraw,Cint,Ptr{PetscDrawHG}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawHGDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscDrawHG},StridedArray{PetscDrawHG},Ptr{Void}))
    ccall((:PetscDrawHGDestroy,petsc1),PetscErrorCode,(Ptr{PetscDrawHG},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawHGAddValue(arg0::Type{Float64},arg1::PetscDrawHG,PetscReal::Integer)
    ccall((:PetscDrawHGAddValue,petsc1),PetscErrorCode,(PetscDrawHG,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawHGDraw(arg0::Type{Float64},arg1::PetscDrawHG)
    ccall((:PetscDrawHGDraw,petsc1),PetscErrorCode,(PetscDrawHG,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawHGView(arg1::PetscDrawHG,arg2::PetscViewer{Float64})
    ccall((:PetscDrawHGView,petsc1),PetscErrorCode,(PetscDrawHG,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawHGReset(arg0::Type{Float64},arg1::PetscDrawHG)
    ccall((:PetscDrawHGReset,petsc1),PetscErrorCode,(PetscDrawHG,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawHGGetAxis(arg0::Type{Float64},arg1::PetscDrawHG,arg2::Union(Ptr{PetscDrawAxis},StridedArray{PetscDrawAxis},Ptr{Void}))
    ccall((:PetscDrawHGGetAxis,petsc1),PetscErrorCode,(PetscDrawHG,Ptr{PetscDrawAxis}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawHGGetDraw(arg0::Type{Float64},arg1::PetscDrawHG,arg2::Union(Ptr{PetscDraw},StridedArray{PetscDraw},Ptr{Void}))
    ccall((:PetscDrawHGGetDraw,petsc1),PetscErrorCode,(PetscDrawHG,Ptr{PetscDraw}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawHGSetLimits(arg0::Type{Float64},arg1::PetscDrawHG,PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:PetscDrawHGSetLimits,petsc1),PetscErrorCode,(PetscDrawHG,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawHGSetNumberBins(arg0::Type{Float64},arg1::PetscDrawHG,arg2::Integer)
    ccall((:PetscDrawHGSetNumberBins,petsc1),PetscErrorCode,(PetscDrawHG,Cint),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawHGSetColor(arg0::Type{Float64},arg1::PetscDrawHG,arg2::Integer)
    ccall((:PetscDrawHGSetColor,petsc1),PetscErrorCode,(PetscDrawHG,Cint),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawHGCalcStats(arg0::Type{Float64},arg1::PetscDrawHG,arg2::PetscBool)
    ccall((:PetscDrawHGCalcStats,petsc1),PetscErrorCode,(PetscDrawHG,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawHGIntegerBins(arg0::Type{Float64},arg1::PetscDrawHG,arg2::PetscBool)
    ccall((:PetscDrawHGIntegerBins,petsc1),PetscErrorCode,(PetscDrawHG,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawBarCreate(arg0::Type{Float64},arg1::PetscDraw,arg2::Union(Ptr{PetscDrawBar},StridedArray{PetscDrawBar},Ptr{Void}))
    ccall((:PetscDrawBarCreate,petsc1),PetscErrorCode,(PetscDraw,Ptr{PetscDrawBar}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawBarSetData(arg0::Type{Float64},arg1::PetscDrawBar,arg2::Integer,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:PetscDrawBarSetData,petsc1),PetscErrorCode,(PetscDrawBar,PetscInt,Ptr{Cint},Ptr{Ptr{Uint8}}),arg1,arg2,PetscReal,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawBarDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscDrawBar},StridedArray{PetscDrawBar},Ptr{Void}))
    ccall((:PetscDrawBarDestroy,petsc1),PetscErrorCode,(Ptr{PetscDrawBar},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawBarDraw(arg0::Type{Float64},arg1::PetscDrawBar)
    ccall((:PetscDrawBarDraw,petsc1),PetscErrorCode,(PetscDrawBar,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawBarSetColor(arg0::Type{Float64},arg1::PetscDrawBar,arg2::Integer)
    ccall((:PetscDrawBarSetColor,petsc1),PetscErrorCode,(PetscDrawBar,Cint),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawBarSetLimits(arg0::Type{Float64},arg1::PetscDrawBar,PetscReal::Integer,arg2::Integer)
    ccall((:PetscDrawBarSetLimits,petsc1),PetscErrorCode,(PetscDrawBar,Cint,Cint),arg1,PetscReal,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawBarSort(arg0::Type{Float64},arg1::PetscDrawBar,arg2::PetscBool,PetscReal::Integer)
    ccall((:PetscDrawBarSort,petsc1),PetscErrorCode,(PetscDrawBar,PetscBool,Cint),arg1,arg2,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawBarSetFromOptions(arg0::Type{Float64},arg1::PetscDrawBar)
    ccall((:PetscDrawBarSetFromOptions,petsc1),PetscErrorCode,(PetscDrawBar,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawBarGetAxis(arg0::Type{Float64},arg1::PetscDrawBar,arg2::Union(Ptr{PetscDrawAxis},StridedArray{PetscDrawAxis},Ptr{Void}))
    ccall((:PetscDrawBarGetAxis,petsc1),PetscErrorCode,(PetscDrawBar,Ptr{PetscDrawAxis}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDrawBarGetDraw(arg0::Type{Float64},arg1::PetscDrawBar,arg2::Union(Ptr{PetscDraw},StridedArray{PetscDraw},Ptr{Void}))
    ccall((:PetscDrawBarGetDraw,petsc1),PetscErrorCode,(PetscDrawBar,Ptr{PetscDraw}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscViewerDrawGetDraw(arg1::PetscViewer{Float64},arg2::Integer,arg3::Union(Ptr{PetscDraw},StridedArray{PetscDraw},Ptr{Void}))
    ccall((:PetscViewerDrawGetDraw,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscInt,Ptr{PetscDraw}),arg1,arg2,arg3)
end 
=#
function PetscViewerDrawBaseAdd(arg1::PetscViewer{Float64},arg2::Integer)
    err = ccall((:PetscViewerDrawBaseAdd,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscInt),arg1,arg2)
    return err
end

function PetscViewerDrawBaseSet(arg1::PetscViewer{Float64},arg2::Integer)
    err = ccall((:PetscViewerDrawBaseSet,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscInt),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function PetscViewerDrawGetDrawLG(arg1::PetscViewer{Float64},arg2::Integer,arg3::Union(Ptr{PetscDrawLG},StridedArray{PetscDrawLG},Ptr{Void}))
    ccall((:PetscViewerDrawGetDrawLG,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscInt,Ptr{PetscDrawLG}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscViewerDrawGetDrawAxis(arg1::PetscViewer{Float64},arg2::Integer,arg3::Union(Ptr{PetscDrawAxis},StridedArray{PetscDrawAxis},Ptr{Void}))
    ccall((:PetscViewerDrawGetDrawAxis,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscInt,Ptr{PetscDrawAxis}),arg1,arg2,arg3)
end 
=#
function PetscDrawUtilitySetCmapHue(arg0::Type{Float64},arg1::Union(Ptr{Cuchar},StridedArray{Cuchar},Ptr{Void}),arg2::Union(Ptr{Cuchar},StridedArray{Cuchar},Ptr{Void}),arg3::Union(Ptr{Cuchar},StridedArray{Cuchar},Ptr{Void}),arg4::Integer)
    err = ccall((:PetscDrawUtilitySetCmapHue,petsc1),PetscErrorCode,(Ptr{Cuchar},Ptr{Cuchar},Ptr{Cuchar},Cint),arg1,arg2,arg3,arg4)
    return err
end

function PetscDrawUtilitySetGamma(arg0::Type{Float64})
    err = ccall((:PetscDrawUtilitySetGamma,petsc1),PetscErrorCode,())
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

function DMInitializePackage(arg0::Type{Float64})
    err = ccall((:DMInitializePackage,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function DMCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMCreate,petsc1),PetscErrorCode,(comm_type,Ptr{DM}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMClone(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMClone,petsc1),PetscErrorCode,(DM,Ptr{DM}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetType(arg0::Type{Float64},arg1::DM,arg2::DMType)
    ccall((:DMSetType,petsc1),PetscErrorCode,(DM,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetType(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{DMType},StridedArray{DMType},Ptr{Void}))
    ccall((:DMGetType,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
function DMRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:DMRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

function DMRegisterDestroy(arg0::Type{Float64})
    err = ccall((:DMRegisterDestroy,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function DMView(arg1::DM,arg2::PetscViewer{Float64})
    ccall((:DMView,petsc1),PetscErrorCode,(DM,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMLoad(arg1::DM,arg2::PetscViewer{Float64})
    ccall((:DMLoad,petsc1),PetscErrorCode,(DM,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDestroy(arg0::Type{Float64},arg1::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMDestroy,petsc1),PetscErrorCode,(Ptr{DM},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function DMCreateGlobalVector(arg1::DM,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMCreateGlobalVector,petsc1),PetscErrorCode,(DM,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMCreateLocalVector(arg1::DM,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMCreateLocalVector,petsc1),PetscErrorCode,(DM,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetLocalVector(arg1::DM,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMGetLocalVector,petsc1),PetscErrorCode,(DM,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMRestoreLocalVector(arg1::DM,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMRestoreLocalVector,petsc1),PetscErrorCode,(DM,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetGlobalVector(arg1::DM,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMGetGlobalVector,petsc1),PetscErrorCode,(DM,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMRestoreGlobalVector(arg1::DM,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMRestoreGlobalVector,petsc1),PetscErrorCode,(DM,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMClearGlobalVectors(arg0::Type{Float64},arg1::DM)
    ccall((:DMClearGlobalVectors,petsc1),PetscErrorCode,(DM,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetNamedGlobalVector(arg1::DM,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMGetNamedGlobalVector,petsc1),PetscErrorCode,(DM,Cstring,Ptr{Vec{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMRestoreNamedGlobalVector(arg1::DM,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMRestoreNamedGlobalVector,petsc1),PetscErrorCode,(DM,Cstring,Ptr{Vec{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetNamedLocalVector(arg1::DM,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMGetNamedLocalVector,petsc1),PetscErrorCode,(DM,Cstring,Ptr{Vec{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMRestoreNamedLocalVector(arg1::DM,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMRestoreNamedLocalVector,petsc1),PetscErrorCode,(DM,Cstring,Ptr{Vec{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetLocalToGlobalMapping(arg1::DM,arg2::Union(Ptr{ISLocalToGlobalMapping{Float64}},StridedArray{ISLocalToGlobalMapping{Float64}},Ptr{Void}))
    ccall((:DMGetLocalToGlobalMapping,petsc1),PetscErrorCode,(DM,Ptr{ISLocalToGlobalMapping{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMCreateFieldIS(arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{Ptr{Uint8}}},StridedArray{Ptr{Ptr{Uint8}}},Ptr{Void}),arg4::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    ccall((:DMCreateFieldIS,petsc1),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{Ptr{Ptr{Uint8}}},Ptr{Ptr{IS{Float64}}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetBlockSize(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMGetBlockSize,petsc1),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMCreateColoring(arg1::DM,arg2::ISColoringType,arg3::Union(Ptr{ISColoring{Float64}},StridedArray{ISColoring{Float64}},Ptr{Void}))
    ccall((:DMCreateColoring,petsc1),PetscErrorCode,(DM,ISColoringType,Ptr{ISColoring{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMCreateMatrix(arg1::DM,arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:DMCreateMatrix,petsc1),PetscErrorCode,(DM,Ptr{Mat{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetMatrixPreallocateOnly(arg0::Type{Float64},arg1::DM,arg2::PetscBool)
    ccall((:DMSetMatrixPreallocateOnly,petsc1),PetscErrorCode,(DM,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMCreateInterpolation(arg1::DM,arg2::DM,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg4::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMCreateInterpolation,petsc1),PetscErrorCode,(DM,DM,Ptr{Mat{Float64}},Ptr{Vec{Float64}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMCreateInjection(arg1::DM,arg2::DM,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:DMCreateInjection,petsc1),PetscErrorCode,(DM,DM,Ptr{Mat{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetWorkArray(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::PetscDataType,arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMGetWorkArray,petsc1),PetscErrorCode,(DM,PetscInt,PetscDataType,Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMRestoreWorkArray(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::PetscDataType,arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMRestoreWorkArray,petsc1),PetscErrorCode,(DM,PetscInt,PetscDataType,Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMRefine(arg0::Type{Float64},arg1::DM,arg2::MPI_Comm,arg3::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMRefine,petsc1),PetscErrorCode,(DM,comm_type,Ptr{DM}),arg1,arg2.val,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMCoarsen(arg0::Type{Float64},arg1::DM,arg2::MPI_Comm,arg3::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMCoarsen,petsc1),PetscErrorCode,(DM,comm_type,Ptr{DM}),arg1,arg2.val,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMRefineHierarchy(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMRefineHierarchy,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{DM}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMCoarsenHierarchy(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMCoarsenHierarchy,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{DM}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMCoarsenHookAdd(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMCoarsenHookAdd,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMRefineHookAdd(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMRefineHookAdd,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMRestrict(arg1::DM,arg2::Mat{Float64},arg3::Vec{Float64},arg4::Mat{Float64},arg5::DM)
    ccall((:DMRestrict,petsc1),PetscErrorCode,(DM,Mat{Float64},Vec{Float64},Mat{Float64},DM),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMInterpolate(arg1::DM,arg2::Mat{Float64},arg3::DM)
    ccall((:DMInterpolate,petsc1),PetscErrorCode,(DM,Mat{Float64},DM),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetFromOptions(arg0::Type{Float64},arg1::DM)
    ccall((:DMSetFromOptions,petsc1),PetscErrorCode,(DM,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function DMCreateInterpolationScale(arg1::DM,arg2::DM,arg3::Mat{Float64},arg4::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMCreateInterpolationScale,petsc1),PetscErrorCode,(DM,DM,Mat{Float64},Ptr{Vec{Float64}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMCreateAggregates(arg1::DM,arg2::DM,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:DMCreateAggregates,petsc1),PetscErrorCode,(DM,DM,Ptr{Mat{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMGlobalToLocalHookAdd(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMGlobalToLocalHookAdd,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMLocalToGlobalHookAdd(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMLocalToGlobalHookAdd,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMGlobalToLocalBegin(arg1::DM,arg2::Vec{Float64},arg3::InsertMode,arg4::Vec{Float64})
    ccall((:DMGlobalToLocalBegin,petsc1),PetscErrorCode,(DM,Vec{Float64},InsertMode,Vec{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMGlobalToLocalEnd(arg1::DM,arg2::Vec{Float64},arg3::InsertMode,arg4::Vec{Float64})
    ccall((:DMGlobalToLocalEnd,petsc1),PetscErrorCode,(DM,Vec{Float64},InsertMode,Vec{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMLocalToGlobalBegin(arg1::DM,arg2::Vec{Float64},arg3::InsertMode,arg4::Vec{Float64})
    ccall((:DMLocalToGlobalBegin,petsc1),PetscErrorCode,(DM,Vec{Float64},InsertMode,Vec{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMLocalToGlobalEnd(arg1::DM,arg2::Vec{Float64},arg3::InsertMode,arg4::Vec{Float64})
    ccall((:DMLocalToGlobalEnd,petsc1),PetscErrorCode,(DM,Vec{Float64},InsertMode,Vec{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMLocalToLocalBegin(arg1::DM,arg2::Vec{Float64},arg3::InsertMode,arg4::Vec{Float64})
    ccall((:DMLocalToLocalBegin,petsc1),PetscErrorCode,(DM,Vec{Float64},InsertMode,Vec{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMLocalToLocalEnd(arg1::DM,arg2::Vec{Float64},arg3::InsertMode,arg4::Vec{Float64})
    ccall((:DMLocalToLocalEnd,petsc1),PetscErrorCode,(DM,Vec{Float64},InsertMode,Vec{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMConvert(arg0::Type{Float64},arg1::DM,arg2::DMType,arg3::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMConvert,petsc1),PetscErrorCode,(DM,Cstring,Ptr{DM}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetDimension(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMGetDimension,petsc1),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetDimension(arg0::Type{Float64},arg1::DM,arg2::Integer)
    ccall((:DMSetDimension,petsc1),PetscErrorCode,(DM,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetDimPoints(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMGetDimPoints,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetCoordinateDM(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMGetCoordinateDM,petsc1),PetscErrorCode,(DM,Ptr{DM}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetCoordinateDM(arg0::Type{Float64},arg1::DM,arg2::DM)
    ccall((:DMSetCoordinateDM,petsc1),PetscErrorCode,(DM,DM),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetCoordinateDim(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMGetCoordinateDim,petsc1),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetCoordinateDim(arg0::Type{Float64},arg1::DM,arg2::Integer)
    ccall((:DMSetCoordinateDim,petsc1),PetscErrorCode,(DM,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetCoordinateSection(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}))
    ccall((:DMGetCoordinateSection,petsc1),PetscErrorCode,(DM,Ptr{PetscSection}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetCoordinateSection(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::PetscSection)
    ccall((:DMSetCoordinateSection,petsc1),PetscErrorCode,(DM,PetscInt,PetscSection),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetCoordinates(arg1::DM,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMGetCoordinates,petsc1),PetscErrorCode,(DM,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetCoordinates(arg1::DM,arg2::Vec{Float64})
    ccall((:DMSetCoordinates,petsc1),PetscErrorCode,(DM,Vec{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetCoordinatesLocal(arg1::DM,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMGetCoordinatesLocal,petsc1),PetscErrorCode,(DM,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetCoordinatesLocal(arg1::DM,arg2::Vec{Float64})
    ccall((:DMSetCoordinatesLocal,petsc1),PetscErrorCode,(DM,Vec{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMLocatePoints(arg1::DM,arg2::Vec{Float64},arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:DMLocatePoints,petsc1),PetscErrorCode,(DM,Vec{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetPeriodicity(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg3::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg4::Union(Ptr{Ptr{DMBoundaryType}},StridedArray{Ptr{DMBoundaryType}},Ptr{Void}))
    ccall((:DMGetPeriodicity,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Cint}},Ptr{Ptr{Cint}},Ptr{Ptr{DMBoundaryType}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetPeriodicity(arg0::Type{Float64},arg1::DM,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{DMBoundaryType},StridedArray{DMBoundaryType},Ptr{Void}))
    ccall((:DMSetPeriodicity,petsc1),PetscErrorCode,(DM,Ptr{Cint},Ptr{Cint},Ptr{DMBoundaryType}),arg1,PetscReal,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMSubDomainHookAdd(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMSubDomainHookAdd,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMSubDomainRestrict(arg1::DM,arg2::VecScatter{Float64},arg3::VecScatter{Float64},arg4::DM)
    ccall((:DMSubDomainRestrict,petsc1),PetscErrorCode,(DM,VecScatter{Float64},VecScatter{Float64},DM),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetOptionsPrefix(arg0::Type{Float64},arg1::DM,arg2::Union(ByteString,Symbol))
    ccall((:DMSetOptionsPrefix,petsc1),PetscErrorCode,(DM,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetVecType(arg0::Type{Float64},arg1::DM,arg2::VecType)
    ccall((:DMSetVecType,petsc1),PetscErrorCode,(DM,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetVecType(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{VecType},StridedArray{VecType},Ptr{Void}))
    ccall((:DMGetVecType,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetMatType(arg0::Type{Float64},arg1::DM,arg2::MatType)
    ccall((:DMSetMatType,petsc1),PetscErrorCode,(DM,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetMatType(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{MatType},StridedArray{MatType},Ptr{Void}))
    ccall((:DMGetMatType,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetApplicationContext(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMSetApplicationContext,petsc1),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetApplicationContextDestroy(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMSetApplicationContextDestroy,petsc1),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetApplicationContext(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMGetApplicationContext,petsc1),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetVariableBounds(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMSetVariableBounds,petsc1),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMHasVariableBounds(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:DMHasVariableBounds,petsc1),PetscErrorCode,(DM,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMHasColoring(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:DMHasColoring,petsc1),PetscErrorCode,(DM,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMComputeVariableBounds(arg1::DM,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:DMComputeVariableBounds,petsc1),PetscErrorCode,(DM,Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMCreateSubDM(arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg5::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMCreateSubDM,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{IS{Float64}},Ptr{DM}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMCreateFieldDecomposition(arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{Ptr{Uint8}}},StridedArray{Ptr{Ptr{Uint8}}},Ptr{Void}),arg4::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}),arg5::Union(Ptr{Ptr{DM}},StridedArray{Ptr{DM}},Ptr{Void}))
    ccall((:DMCreateFieldDecomposition,petsc1),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{Ptr{Ptr{Uint8}}},Ptr{Ptr{IS{Float64}}},Ptr{Ptr{DM}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMCreateDomainDecomposition(arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{Ptr{Uint8}}},StridedArray{Ptr{Ptr{Uint8}}},Ptr{Void}),arg4::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}),arg5::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}),arg6::Union(Ptr{Ptr{DM}},StridedArray{Ptr{DM}},Ptr{Void}))
    ccall((:DMCreateDomainDecomposition,petsc1),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{Ptr{Ptr{Uint8}}},Ptr{Ptr{IS{Float64}}},Ptr{Ptr{IS{Float64}}},Ptr{Ptr{DM}}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function DMCreateDomainDecompositionScatters(arg1::DM,arg2::Integer,arg3::Union(Ptr{DM},StridedArray{DM},Ptr{Void}),arg4::Union(Ptr{Ptr{VecScatter{Float64}}},StridedArray{Ptr{VecScatter{Float64}}},Ptr{Void}),arg5::Union(Ptr{Ptr{VecScatter{Float64}}},StridedArray{Ptr{VecScatter{Float64}}},Ptr{Void}),arg6::Union(Ptr{Ptr{VecScatter{Float64}}},StridedArray{Ptr{VecScatter{Float64}}},Ptr{Void}))
    ccall((:DMCreateDomainDecompositionScatters,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{DM},Ptr{Ptr{VecScatter{Float64}}},Ptr{Ptr{VecScatter{Float64}}},Ptr{Ptr{VecScatter{Float64}}}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetRefineLevel(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMGetRefineLevel,petsc1),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetCoarsenLevel(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMGetCoarsenLevel,petsc1),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end 
=#
function DMFinalizePackage(arg0::Type{Float64})
    err = ccall((:DMFinalizePackage,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function VecGetDM(arg1::Vec{Float64},arg2::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:VecGetDM,petsc1),PetscErrorCode,(Vec{Float64},Ptr{DM}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function VecSetDM(arg1::Vec{Float64},arg2::DM)
    ccall((:VecSetDM,petsc1),PetscErrorCode,(Vec{Float64},DM),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatGetDM(arg1::Mat{Float64},arg2::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:MatGetDM,petsc1),PetscErrorCode,(Mat{Float64},Ptr{DM}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatSetDM(arg1::Mat{Float64},arg2::DM)
    ccall((:MatSetDM,petsc1),PetscErrorCode,(Mat{Float64},DM),arg1,arg2)
end 
=#
function DMPrintCellVector(arg0::Type{Float64},arg1::Integer,arg2::Union(ByteString,Symbol),arg3::Integer,arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:DMPrintCellVector,petsc1),PetscErrorCode,(PetscInt,Cstring,PetscInt,Ptr{Float64}),arg1,arg2,arg3,arg4)
    return err
end

function DMPrintCellMatrix(arg0::Type{Float64},arg1::Integer,arg2::Union(ByteString,Symbol),arg3::Integer,arg4::Integer,arg5::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    err = ccall((:DMPrintCellMatrix,petsc1),PetscErrorCode,(PetscInt,Cstring,PetscInt,PetscInt,Ptr{Float64}),arg1,arg2,arg3,arg4,arg5)
    return err
end

#= skipping function with undefined symbols: 
 function DMPrintLocalVec(arg1::DM,arg2::Union(ByteString,Symbol),PetscReal::Integer,arg3::Vec{Float64})
    ccall((:DMPrintLocalVec,petsc1),PetscErrorCode,(DM,Cstring,Cint,Vec{Float64}),arg1,arg2,PetscReal,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetDefaultSection(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}))
    ccall((:DMGetDefaultSection,petsc1),PetscErrorCode,(DM,Ptr{PetscSection}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetDefaultSection(arg0::Type{Float64},arg1::DM,arg2::PetscSection)
    ccall((:DMSetDefaultSection,petsc1),PetscErrorCode,(DM,PetscSection),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetDefaultConstraints(arg1::DM,arg2::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}),arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:DMGetDefaultConstraints,petsc1),PetscErrorCode,(DM,Ptr{PetscSection},Ptr{Mat{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetDefaultConstraints(arg1::DM,arg2::PetscSection,arg3::Mat{Float64})
    ccall((:DMSetDefaultConstraints,petsc1),PetscErrorCode,(DM,PetscSection,Mat{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetDefaultGlobalSection(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}))
    ccall((:DMGetDefaultGlobalSection,petsc1),PetscErrorCode,(DM,Ptr{PetscSection}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetDefaultGlobalSection(arg0::Type{Float64},arg1::DM,arg2::PetscSection)
    ccall((:DMSetDefaultGlobalSection,petsc1),PetscErrorCode,(DM,PetscSection),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetDefaultSF(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscSF},StridedArray{PetscSF},Ptr{Void}))
    ccall((:DMGetDefaultSF,petsc1),PetscErrorCode,(DM,Ptr{PetscSF}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetDefaultSF(arg0::Type{Float64},arg1::DM,arg2::PetscSF)
    ccall((:DMSetDefaultSF,petsc1),PetscErrorCode,(DM,PetscSF),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMCreateDefaultSF(arg0::Type{Float64},arg1::DM,arg2::PetscSection,arg3::PetscSection)
    ccall((:DMCreateDefaultSF,petsc1),PetscErrorCode,(DM,PetscSection,PetscSection),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetPointSF(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscSF},StridedArray{PetscSF},Ptr{Void}))
    ccall((:DMGetPointSF,petsc1),PetscErrorCode,(DM,Ptr{PetscSF}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetPointSF(arg0::Type{Float64},arg1::DM,arg2::PetscSF)
    ccall((:DMSetPointSF,petsc1),PetscErrorCode,(DM,PetscSF),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetOutputDM(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMGetOutputDM,petsc1),PetscErrorCode,(DM,Ptr{DM}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetOutputSequenceNumber(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:DMGetOutputSequenceNumber,petsc1),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{Cint}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetOutputSequenceNumber(arg0::Type{Float64},arg1::DM,arg2::Integer,PetscReal::Integer)
    ccall((:DMSetOutputSequenceNumber,petsc1),PetscErrorCode,(DM,PetscInt,Cint),arg1,arg2,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function DMOutputSequenceLoad(arg1::DM,arg2::PetscViewer{Float64},arg3::Union(ByteString,Symbol),arg4::Integer,arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:DMOutputSequenceLoad,petsc1),PetscErrorCode,(DM,PetscViewer{Float64},Cstring,PetscInt,Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetDS(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscDS},StridedArray{PetscDS},Ptr{Void}))
    ccall((:DMGetDS,petsc1),PetscErrorCode,(DM,Ptr{PetscDS}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetDS(arg0::Type{Float64},arg1::DM,arg2::PetscDS)
    ccall((:DMSetDS,petsc1),PetscErrorCode,(DM,PetscDS),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetNumFields(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMGetNumFields,petsc1),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetNumFields(arg0::Type{Float64},arg1::DM,arg2::Integer)
    ccall((:DMSetNumFields,petsc1),PetscErrorCode,(DM,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMGetField(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscObject},StridedArray{PetscObject},Ptr{Void}))
    ccall((:DMGetField,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscObject}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMSetField(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::PetscObject)
    ccall((:DMSetField,petsc1),PetscErrorCode,(DM,PetscInt,PetscObject),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMInterpolationCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{DMInterpolationInfo},StridedArray{DMInterpolationInfo},Ptr{Void}))
    ccall((:DMInterpolationCreate,petsc1),PetscErrorCode,(comm_type,Ptr{DMInterpolationInfo}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMInterpolationSetDim(arg0::Type{Float64},arg1::DMInterpolationInfo,arg2::Integer)
    ccall((:DMInterpolationSetDim,petsc1),PetscErrorCode,(DMInterpolationInfo,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMInterpolationGetDim(arg0::Type{Float64},arg1::DMInterpolationInfo,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMInterpolationGetDim,petsc1),PetscErrorCode,(DMInterpolationInfo,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMInterpolationSetDof(arg0::Type{Float64},arg1::DMInterpolationInfo,arg2::Integer)
    ccall((:DMInterpolationSetDof,petsc1),PetscErrorCode,(DMInterpolationInfo,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMInterpolationGetDof(arg0::Type{Float64},arg1::DMInterpolationInfo,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMInterpolationGetDof,petsc1),PetscErrorCode,(DMInterpolationInfo,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMInterpolationAddPoints(arg0::Type{Float64},arg1::DMInterpolationInfo,arg2::Integer,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:DMInterpolationAddPoints,petsc1),PetscErrorCode,(DMInterpolationInfo,PetscInt,Ptr{Cint}),arg1,arg2,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function DMInterpolationSetUp(arg0::Type{Float64},arg1::DMInterpolationInfo,arg2::DM,arg3::PetscBool)
    ccall((:DMInterpolationSetUp,petsc1),PetscErrorCode,(DMInterpolationInfo,DM,PetscBool),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMInterpolationGetCoordinates(arg1::DMInterpolationInfo,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMInterpolationGetCoordinates,petsc1),PetscErrorCode,(DMInterpolationInfo,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMInterpolationGetVector(arg1::DMInterpolationInfo,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMInterpolationGetVector,petsc1),PetscErrorCode,(DMInterpolationInfo,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMInterpolationRestoreVector(arg1::DMInterpolationInfo,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMInterpolationRestoreVector,petsc1),PetscErrorCode,(DMInterpolationInfo,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMInterpolationEvaluate(arg1::DMInterpolationInfo,arg2::DM,arg3::Vec{Float64},arg4::Vec{Float64})
    ccall((:DMInterpolationEvaluate,petsc1),PetscErrorCode,(DMInterpolationInfo,DM,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMInterpolationDestroy(arg0::Type{Float64},arg1::Union(Ptr{DMInterpolationInfo},StridedArray{DMInterpolationInfo},Ptr{Void}))
    ccall((:DMInterpolationDestroy,petsc1),PetscErrorCode,(Ptr{DMInterpolationInfo},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PFCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Union(Ptr{PF},StridedArray{PF},Ptr{Void}))
    ccall((:PFCreate,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,Ptr{PF}),arg1.val,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PFSetType(arg0::Type{Float64},arg1::PF,arg2::PFType,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PFSetType,petsc1),PetscErrorCode,(PF,Cstring,Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PFSet(arg0::Type{Float64},arg1::PF,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg6::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PFSet,petsc1),PetscErrorCode,(PF,Ptr{Void},Ptr{Void},Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function PFApply(arg0::Type{Float64},arg1::PF,arg2::Integer,arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:PFApply,petsc1),PetscErrorCode,(PF,PetscInt,Ptr{Float64},Ptr{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PFApplyVec(arg1::PF,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:PFApplyVec,petsc1),PetscErrorCode,(PF,Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end 
=#
function PFInitializePackage(arg0::Type{Float64})
    err = ccall((:PFInitializePackage,petsc1),PetscErrorCode,())
    return err
end

function PFRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PFRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function PFDestroy(arg0::Type{Float64},arg1::Union(Ptr{PF},StridedArray{PF},Ptr{Void}))
    ccall((:PFDestroy,petsc1),PetscErrorCode,(Ptr{PF},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PFSetFromOptions(arg0::Type{Float64},arg1::PF)
    ccall((:PFSetFromOptions,petsc1),PetscErrorCode,(PF,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PFGetType(arg0::Type{Float64},arg1::PF,arg2::Union(Ptr{PFType},StridedArray{PFType},Ptr{Void}))
    ccall((:PFGetType,petsc1),PetscErrorCode,(PF,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PFView(arg1::PF,arg2::PetscViewer{Float64})
    ccall((:PFView,petsc1),PetscErrorCode,(PF,PetscViewer{Float64}),arg1,arg2)
end 
=#
function AOInitializePackage(arg0::Type{Float64})
    err = ccall((:AOInitializePackage,petsc1),PetscErrorCode,())
    return err
end

function AOCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:AOCreate,petsc1),PetscErrorCode,(comm_type,Ptr{Cint}),arg1.val,arg2)
    return err
end

function AOSetIS(arg0::Type{Float64})
    err = ccall((:AOSetIS,petsc1),PetscErrorCode,())
    return err
end

function AOSetFromOptions(arg0::Type{Float64})
    err = ccall((:AOSetFromOptions,petsc1),PetscErrorCode,())
    return err
end

function AOCreateBasic(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:AOCreateBasic,petsc1),PetscErrorCode,(comm_type,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Cint}),arg1.val,arg2,arg3,arg4,arg5)
    return err
end

function AOCreateBasicIS(arg1::IS{Float64},arg2::IS{Float64},arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:AOCreateBasicIS,petsc1),PetscErrorCode,(IS{Float64},IS{Float64},Ptr{Cint}),arg1,arg2,arg3)
    return err
end

function AOCreateMemoryScalable(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:AOCreateMemoryScalable,petsc1),PetscErrorCode,(comm_type,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Cint}),arg1.val,arg2,arg3,arg4,arg5)
    return err
end

function AOCreateMemoryScalableIS(arg1::IS{Float64},arg2::IS{Float64},arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:AOCreateMemoryScalableIS,petsc1),PetscErrorCode,(IS{Float64},IS{Float64},Ptr{Cint}),arg1,arg2,arg3)
    return err
end

function AOCreateMapping(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:AOCreateMapping,petsc1),PetscErrorCode,(comm_type,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Cint}),arg1.val,arg2,arg3,arg4,arg5)
    return err
end

function AOCreateMappingIS(arg1::IS{Float64},arg2::IS{Float64},arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:AOCreateMappingIS,petsc1),PetscErrorCode,(IS{Float64},IS{Float64},Ptr{Cint}),arg1,arg2,arg3)
    return err
end

function AOView(arg0::Type{Float64})
    err = ccall((:AOView,petsc1),PetscErrorCode,())
    return err
end

function AOSetType(arg0::Type{Float64})
    err = ccall((:AOSetType,petsc1),PetscErrorCode,())
    return err
end

function AOGetType(arg0::Type{Float64})
    err = ccall((:AOGetType,petsc1),PetscErrorCode,())
    return err
end

function AORegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:AORegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

function AOPetscToApplication(arg0::Type{Float64})
    err = ccall((:AOPetscToApplication,petsc1),PetscErrorCode,())
    return err
end

function AOApplicationToPetsc(arg0::Type{Float64})
    err = ccall((:AOApplicationToPetsc,petsc1),PetscErrorCode,())
    return err
end

function AOPetscToApplicationIS(arg0::Type{Float64})
    err = ccall((:AOPetscToApplicationIS,petsc1),PetscErrorCode,())
    return err
end

function AOApplicationToPetscIS(arg0::Type{Float64})
    err = ccall((:AOApplicationToPetscIS,petsc1),PetscErrorCode,())
    return err
end

function AOPetscToApplicationPermuteInt(arg0::Type{Float64})
    err = ccall((:AOPetscToApplicationPermuteInt,petsc1),PetscErrorCode,())
    return err
end

function AOApplicationToPetscPermuteInt(arg0::Type{Float64})
    err = ccall((:AOApplicationToPetscPermuteInt,petsc1),PetscErrorCode,())
    return err
end

function AOPetscToApplicationPermuteReal(arg0::Type{Float64})
    err = ccall((:AOPetscToApplicationPermuteReal,petsc1),PetscErrorCode,())
    return err
end

function AOApplicationToPetscPermuteReal(arg0::Type{Float64})
    err = ccall((:AOApplicationToPetscPermuteReal,petsc1),PetscErrorCode,())
    return err
end

function AOMappingHasApplicationIndex(arg0::Type{Float64})
    err = ccall((:AOMappingHasApplicationIndex,petsc1),PetscErrorCode,())
    return err
end

function AOMappingHasPetscIndex(arg0::Type{Float64})
    err = ccall((:AOMappingHasPetscIndex,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function PetscQuadratureCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscQuadrature},StridedArray{PetscQuadrature},Ptr{Void}))
    ccall((:PetscQuadratureCreate,petsc1),PetscErrorCode,(comm_type,Ptr{PetscQuadrature}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscQuadratureDuplicate(arg0::Type{Float64},arg1::PetscQuadrature,arg2::Union(Ptr{PetscQuadrature},StridedArray{PetscQuadrature},Ptr{Void}))
    ccall((:PetscQuadratureDuplicate,petsc1),PetscErrorCode,(PetscQuadrature,Ptr{PetscQuadrature}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscQuadratureGetOrder(arg0::Type{Float64},arg1::PetscQuadrature,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscQuadratureGetOrder,petsc1),PetscErrorCode,(PetscQuadrature,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscQuadratureSetOrder(arg0::Type{Float64},arg1::PetscQuadrature,arg2::Integer)
    ccall((:PetscQuadratureSetOrder,petsc1),PetscErrorCode,(PetscQuadrature,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscQuadratureGetData(arg0::Type{Float64},arg1::PetscQuadrature,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg5::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}))
    ccall((:PetscQuadratureGetData,petsc1),PetscErrorCode,(PetscQuadrature,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscQuadratureSetData(arg0::Type{Float64},arg1::PetscQuadrature,arg2::Integer,arg3::Integer,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscQuadratureSetData,petsc1),PetscErrorCode,(PetscQuadrature,PetscInt,PetscInt,Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,PetscReal,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscQuadratureView(arg1::PetscQuadrature,arg2::PetscViewer{Float64})
    ccall((:PetscQuadratureView,petsc1),PetscErrorCode,(PetscQuadrature,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscQuadratureDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscQuadrature},StridedArray{PetscQuadrature},Ptr{Void}))
    ccall((:PetscQuadratureDestroy,petsc1),PetscErrorCode,(Ptr{PetscQuadrature},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscQuadratureExpandComposite(arg0::Type{Float64},arg1::PetscQuadrature,arg2::Integer,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{PetscQuadrature},StridedArray{PetscQuadrature},Ptr{Void}))
    ccall((:PetscQuadratureExpandComposite,petsc1),PetscErrorCode,(PetscQuadrature,PetscInt,Ptr{Cint},Ptr{Cint},Ptr{PetscQuadrature}),arg1,arg2,PetscReal,arg3,arg4)
end 
=#
function PetscDTLegendreEval(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg7::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:PetscDTLegendreEval,petsc1),PetscErrorCode,(PetscInt,Ptr{Cint},PetscInt,Ptr{PetscInt},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
    return err
end

function PetscDTGaussQuadrature(arg0::Type{Float64},arg1::Integer,PetscReal::Integer,arg2::Integer,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:PetscDTGaussQuadrature,petsc1),PetscErrorCode,(PetscInt,Cint,Cint,Ptr{Cint},Ptr{Cint}),arg1,PetscReal,arg2,arg3,arg4)
    return err
end

function PetscDTReconstructPoly(arg0::Type{Float64},arg1::Integer,arg2::Integer,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:PetscDTReconstructPoly,petsc1),PetscErrorCode,(PetscInt,PetscInt,Ptr{Cint},PetscInt,Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6)
    return err
end

#= skipping function with undefined symbols: 
 function PetscDTGaussTensorQuadrature(arg0::Type{Float64},arg1::Integer,arg2::Integer,PetscReal::Integer,arg3::Integer,arg4::Union(Ptr{PetscQuadrature},StridedArray{PetscQuadrature},Ptr{Void}))
    ccall((:PetscDTGaussTensorQuadrature,petsc1),PetscErrorCode,(PetscInt,PetscInt,Cint,Cint,Ptr{PetscQuadrature}),arg1,arg2,PetscReal,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDTGaussJacobiQuadrature(arg0::Type{Float64},arg1::Integer,arg2::Integer,PetscReal::Integer,arg3::Integer,arg4::Union(Ptr{PetscQuadrature},StridedArray{PetscQuadrature},Ptr{Void}))
    ccall((:PetscDTGaussJacobiQuadrature,petsc1),PetscErrorCode,(PetscInt,PetscInt,Cint,Cint,Ptr{PetscQuadrature}),arg1,arg2,PetscReal,arg3,arg4)
end 
=#
function PetscFEInitializePackage(arg0::Type{Float64})
    err = ccall((:PetscFEInitializePackage,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function PetscSpaceCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscSpace},StridedArray{PetscSpace},Ptr{Void}))
    ccall((:PetscSpaceCreate,petsc1),PetscErrorCode,(comm_type,Ptr{PetscSpace}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSpaceDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscSpace},StridedArray{PetscSpace},Ptr{Void}))
    ccall((:PetscSpaceDestroy,petsc1),PetscErrorCode,(Ptr{PetscSpace},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSpaceSetType(arg0::Type{Float64},arg1::PetscSpace,arg2::PetscSpaceType)
    ccall((:PetscSpaceSetType,petsc1),PetscErrorCode,(PetscSpace,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSpaceGetType(arg0::Type{Float64},arg1::PetscSpace,arg2::Union(Ptr{PetscSpaceType},StridedArray{PetscSpaceType},Ptr{Void}))
    ccall((:PetscSpaceGetType,petsc1),PetscErrorCode,(PetscSpace,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSpaceSetUp(arg0::Type{Float64},arg1::PetscSpace)
    ccall((:PetscSpaceSetUp,petsc1),PetscErrorCode,(PetscSpace,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSpaceSetFromOptions(arg0::Type{Float64},arg1::PetscSpace)
    ccall((:PetscSpaceSetFromOptions,petsc1),PetscErrorCode,(PetscSpace,),arg1)
end 
=#
function PetscSpaceRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscSpaceRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

function PetscSpaceRegisterDestroy(arg0::Type{Float64})
    err = ccall((:PetscSpaceRegisterDestroy,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function PetscSpaceGetDimension(arg0::Type{Float64},arg1::PetscSpace,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscSpaceGetDimension,petsc1),PetscErrorCode,(PetscSpace,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSpaceSetOrder(arg0::Type{Float64},arg1::PetscSpace,arg2::Integer)
    ccall((:PetscSpaceSetOrder,petsc1),PetscErrorCode,(PetscSpace,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSpaceGetOrder(arg0::Type{Float64},arg1::PetscSpace,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscSpaceGetOrder,petsc1),PetscErrorCode,(PetscSpace,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSpaceEvaluate(arg0::Type{Float64},arg1::PetscSpace,arg2::Integer,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscSpaceEvaluate,petsc1),PetscErrorCode,(PetscSpace,PetscInt,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,PetscReal,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSpacePolynomialSetNumVariables(arg0::Type{Float64},arg1::PetscSpace,arg2::Integer)
    ccall((:PetscSpacePolynomialSetNumVariables,petsc1),PetscErrorCode,(PetscSpace,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSpacePolynomialGetNumVariables(arg0::Type{Float64},arg1::PetscSpace,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscSpacePolynomialGetNumVariables,petsc1),PetscErrorCode,(PetscSpace,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSpacePolynomialSetSymmetric(arg0::Type{Float64},arg1::PetscSpace,arg2::PetscBool)
    ccall((:PetscSpacePolynomialSetSymmetric,petsc1),PetscErrorCode,(PetscSpace,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSpacePolynomialGetSymmetric(arg0::Type{Float64},arg1::PetscSpace,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscSpacePolynomialGetSymmetric,petsc1),PetscErrorCode,(PetscSpace,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSpacePolynomialSetTensor(arg0::Type{Float64},arg1::PetscSpace,arg2::PetscBool)
    ccall((:PetscSpacePolynomialSetTensor,petsc1),PetscErrorCode,(PetscSpace,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSpacePolynomialGetTensor(arg0::Type{Float64},arg1::PetscSpace,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscSpacePolynomialGetTensor,petsc1),PetscErrorCode,(PetscSpace,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSpaceDGSetQuadrature(arg0::Type{Float64},arg1::PetscSpace,arg2::PetscQuadrature)
    ccall((:PetscSpaceDGSetQuadrature,petsc1),PetscErrorCode,(PetscSpace,PetscQuadrature),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSpaceDGGetQuadrature(arg0::Type{Float64},arg1::PetscSpace,arg2::Union(Ptr{PetscQuadrature},StridedArray{PetscQuadrature},Ptr{Void}))
    ccall((:PetscSpaceDGGetQuadrature,petsc1),PetscErrorCode,(PetscSpace,Ptr{PetscQuadrature}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDualSpaceCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscDualSpace},StridedArray{PetscDualSpace},Ptr{Void}))
    ccall((:PetscDualSpaceCreate,petsc1),PetscErrorCode,(comm_type,Ptr{PetscDualSpace}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDualSpaceDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscDualSpace},StridedArray{PetscDualSpace},Ptr{Void}))
    ccall((:PetscDualSpaceDestroy,petsc1),PetscErrorCode,(Ptr{PetscDualSpace},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDualSpaceDuplicate(arg0::Type{Float64},arg1::PetscDualSpace,arg2::Union(Ptr{PetscDualSpace},StridedArray{PetscDualSpace},Ptr{Void}))
    ccall((:PetscDualSpaceDuplicate,petsc1),PetscErrorCode,(PetscDualSpace,Ptr{PetscDualSpace}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDualSpaceSetType(arg0::Type{Float64},arg1::PetscDualSpace,arg2::PetscDualSpaceType)
    ccall((:PetscDualSpaceSetType,petsc1),PetscErrorCode,(PetscDualSpace,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDualSpaceGetType(arg0::Type{Float64},arg1::PetscDualSpace,arg2::Union(Ptr{PetscDualSpaceType},StridedArray{PetscDualSpaceType},Ptr{Void}))
    ccall((:PetscDualSpaceGetType,petsc1),PetscErrorCode,(PetscDualSpace,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDualSpaceSetUp(arg0::Type{Float64},arg1::PetscDualSpace)
    ccall((:PetscDualSpaceSetUp,petsc1),PetscErrorCode,(PetscDualSpace,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDualSpaceSetFromOptions(arg0::Type{Float64},arg1::PetscDualSpace)
    ccall((:PetscDualSpaceSetFromOptions,petsc1),PetscErrorCode,(PetscDualSpace,),arg1)
end 
=#
function PetscDualSpaceRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscDualSpaceRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

function PetscDualSpaceRegisterDestroy(arg0::Type{Float64})
    err = ccall((:PetscDualSpaceRegisterDestroy,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function PetscDualSpaceGetDimension(arg0::Type{Float64},arg1::PetscDualSpace,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscDualSpaceGetDimension,petsc1),PetscErrorCode,(PetscDualSpace,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDualSpaceSetOrder(arg0::Type{Float64},arg1::PetscDualSpace,arg2::Integer)
    ccall((:PetscDualSpaceSetOrder,petsc1),PetscErrorCode,(PetscDualSpace,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDualSpaceGetOrder(arg0::Type{Float64},arg1::PetscDualSpace,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscDualSpaceGetOrder,petsc1),PetscErrorCode,(PetscDualSpace,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDualSpaceSetDM(arg0::Type{Float64},arg1::PetscDualSpace,arg2::DM)
    ccall((:PetscDualSpaceSetDM,petsc1),PetscErrorCode,(PetscDualSpace,DM),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDualSpaceGetDM(arg0::Type{Float64},arg1::PetscDualSpace,arg2::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:PetscDualSpaceGetDM,petsc1),PetscErrorCode,(PetscDualSpace,Ptr{DM}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDualSpaceGetFunctional(arg0::Type{Float64},arg1::PetscDualSpace,arg2::Integer,arg3::Union(Ptr{PetscQuadrature},StridedArray{PetscQuadrature},Ptr{Void}))
    ccall((:PetscDualSpaceGetFunctional,petsc1),PetscErrorCode,(PetscDualSpace,PetscInt,Ptr{PetscQuadrature}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDualSpaceCreateReferenceCell(arg0::Type{Float64},arg1::PetscDualSpace,arg2::Integer,arg3::PetscBool,arg4::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:PetscDualSpaceCreateReferenceCell,petsc1),PetscErrorCode,(PetscDualSpace,PetscInt,PetscBool,Ptr{DM}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDualSpaceApply(arg0::Type{Float64},arg1::PetscDualSpace,arg2::Integer,arg3::Union(Ptr{PetscFECellGeom},StridedArray{PetscFECellGeom},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg6::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg7::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:PetscDualSpaceApply,petsc1),PetscErrorCode,(PetscDualSpace,PetscInt,Ptr{PetscFECellGeom},PetscInt,Ptr{Void},Ptr{Void},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDualSpaceLagrangeGetContinuity(arg0::Type{Float64},arg1::PetscDualSpace,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscDualSpaceLagrangeGetContinuity,petsc1),PetscErrorCode,(PetscDualSpace,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDualSpaceLagrangeSetContinuity(arg0::Type{Float64},arg1::PetscDualSpace,arg2::PetscBool)
    ccall((:PetscDualSpaceLagrangeSetContinuity,petsc1),PetscErrorCode,(PetscDualSpace,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDualSpaceGetHeightSubspace(arg0::Type{Float64},arg1::PetscDualSpace,arg2::Integer,arg3::Union(Ptr{PetscDualSpace},StridedArray{PetscDualSpace},Ptr{Void}))
    ccall((:PetscDualSpaceGetHeightSubspace,petsc1),PetscErrorCode,(PetscDualSpace,PetscInt,Ptr{PetscDualSpace}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDualSpaceSimpleSetDimension(arg0::Type{Float64},arg1::PetscDualSpace,arg2::Integer)
    ccall((:PetscDualSpaceSimpleSetDimension,petsc1),PetscErrorCode,(PetscDualSpace,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDualSpaceSimpleSetFunctional(arg0::Type{Float64},arg1::PetscDualSpace,arg2::Integer,arg3::PetscQuadrature)
    ccall((:PetscDualSpaceSimpleSetFunctional,petsc1),PetscErrorCode,(PetscDualSpace,PetscInt,PetscQuadrature),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFECreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscFE},StridedArray{PetscFE},Ptr{Void}))
    ccall((:PetscFECreate,petsc1),PetscErrorCode,(comm_type,Ptr{PetscFE}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFEDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscFE},StridedArray{PetscFE},Ptr{Void}))
    ccall((:PetscFEDestroy,petsc1),PetscErrorCode,(Ptr{PetscFE},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFESetType(arg0::Type{Float64},arg1::PetscFE,arg2::PetscFEType)
    ccall((:PetscFESetType,petsc1),PetscErrorCode,(PetscFE,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFEGetType(arg0::Type{Float64},arg1::PetscFE,arg2::Union(Ptr{PetscFEType},StridedArray{PetscFEType},Ptr{Void}))
    ccall((:PetscFEGetType,petsc1),PetscErrorCode,(PetscFE,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFESetUp(arg0::Type{Float64},arg1::PetscFE)
    ccall((:PetscFESetUp,petsc1),PetscErrorCode,(PetscFE,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFESetFromOptions(arg0::Type{Float64},arg1::PetscFE)
    ccall((:PetscFESetFromOptions,petsc1),PetscErrorCode,(PetscFE,),arg1)
end 
=#
function PetscFERegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscFERegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

function PetscFERegisterDestroy(arg0::Type{Float64})
    err = ccall((:PetscFERegisterDestroy,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function PetscFECreateDefault(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer,arg4::PetscBool,arg5::Union(ByteString,Symbol),arg6::Integer,arg7::Union(Ptr{PetscFE},StridedArray{PetscFE},Ptr{Void}))
    ccall((:PetscFECreateDefault,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,PetscBool,Cstring,PetscInt,Ptr{PetscFE}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFEGetDimension(arg0::Type{Float64},arg1::PetscFE,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscFEGetDimension,petsc1),PetscErrorCode,(PetscFE,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFEGetSpatialDimension(arg0::Type{Float64},arg1::PetscFE,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscFEGetSpatialDimension,petsc1),PetscErrorCode,(PetscFE,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFESetNumComponents(arg0::Type{Float64},arg1::PetscFE,arg2::Integer)
    ccall((:PetscFESetNumComponents,petsc1),PetscErrorCode,(PetscFE,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFEGetNumComponents(arg0::Type{Float64},arg1::PetscFE,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscFEGetNumComponents,petsc1),PetscErrorCode,(PetscFE,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFEGetTileSizes(arg0::Type{Float64},arg1::PetscFE,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscFEGetTileSizes,petsc1),PetscErrorCode,(PetscFE,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFESetTileSizes(arg0::Type{Float64},arg1::PetscFE,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer)
    ccall((:PetscFESetTileSizes,petsc1),PetscErrorCode,(PetscFE,PetscInt,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFESetBasisSpace(arg0::Type{Float64},arg1::PetscFE,arg2::PetscSpace)
    ccall((:PetscFESetBasisSpace,petsc1),PetscErrorCode,(PetscFE,PetscSpace),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFEGetBasisSpace(arg0::Type{Float64},arg1::PetscFE,arg2::Union(Ptr{PetscSpace},StridedArray{PetscSpace},Ptr{Void}))
    ccall((:PetscFEGetBasisSpace,petsc1),PetscErrorCode,(PetscFE,Ptr{PetscSpace}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFESetDualSpace(arg0::Type{Float64},arg1::PetscFE,arg2::PetscDualSpace)
    ccall((:PetscFESetDualSpace,petsc1),PetscErrorCode,(PetscFE,PetscDualSpace),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFEGetDualSpace(arg0::Type{Float64},arg1::PetscFE,arg2::Union(Ptr{PetscDualSpace},StridedArray{PetscDualSpace},Ptr{Void}))
    ccall((:PetscFEGetDualSpace,petsc1),PetscErrorCode,(PetscFE,Ptr{PetscDualSpace}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFESetQuadrature(arg0::Type{Float64},arg1::PetscFE,arg2::PetscQuadrature)
    ccall((:PetscFESetQuadrature,petsc1),PetscErrorCode,(PetscFE,PetscQuadrature),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFEGetQuadrature(arg0::Type{Float64},arg1::PetscFE,arg2::Union(Ptr{PetscQuadrature},StridedArray{PetscQuadrature},Ptr{Void}))
    ccall((:PetscFEGetQuadrature,petsc1),PetscErrorCode,(PetscFE,Ptr{PetscQuadrature}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFEGetNumDof(arg0::Type{Float64},arg1::PetscFE,arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:PetscFEGetNumDof,petsc1),PetscErrorCode,(PetscFE,Ptr{Ptr{PetscInt}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFEGetDefaultTabulation(arg0::Type{Float64},arg1::PetscFE,arg2::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg3::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg4::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}))
    ccall((:PetscFEGetDefaultTabulation,petsc1),PetscErrorCode,(PetscFE,Ptr{Ptr{Cint}},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFEGetFaceTabulation(arg0::Type{Float64},arg1::PetscFE,arg2::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}))
    ccall((:PetscFEGetFaceTabulation,petsc1),PetscErrorCode,(PetscFE,Ptr{Ptr{Cint}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFEGetTabulation(arg0::Type{Float64},arg1::PetscFE,arg2::Integer,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg4::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg5::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}))
    ccall((:PetscFEGetTabulation,petsc1),PetscErrorCode,(PetscFE,PetscInt,Ptr{Cint},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}}),arg1,arg2,PetscReal,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFERestoreTabulation(arg0::Type{Float64},arg1::PetscFE,arg2::Integer,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg4::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg5::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}))
    ccall((:PetscFERestoreTabulation,petsc1),PetscErrorCode,(PetscFE,PetscInt,Ptr{Cint},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}}),arg1,arg2,PetscReal,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFERefine(arg0::Type{Float64},arg1::PetscFE,arg2::Union(Ptr{PetscFE},StridedArray{PetscFE},Ptr{Void}))
    ccall((:PetscFERefine,petsc1),PetscErrorCode,(PetscFE,Ptr{PetscFE}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFEIntegrate(arg0::Type{Float64},arg1::PetscFE,arg2::PetscDS,arg3::Integer,arg4::Integer,arg5::Union(Ptr{PetscFECellGeom},StridedArray{PetscFECellGeom},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::PetscDS,arg8::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscFEIntegrate,petsc1),PetscErrorCode,(PetscFE,PetscDS,PetscInt,PetscInt,Ptr{PetscFECellGeom},Ptr{Float64},PetscDS,Ptr{Float64},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFEIntegrateResidual(arg0::Type{Float64},arg1::PetscFE,arg2::PetscDS,arg3::Integer,arg4::Integer,arg5::Union(Ptr{PetscFECellGeom},StridedArray{PetscFECellGeom},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg8::PetscDS,arg9::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg10::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:PetscFEIntegrateResidual,petsc1),PetscErrorCode,(PetscFE,PetscDS,PetscInt,PetscInt,Ptr{PetscFECellGeom},Ptr{Float64},Ptr{Float64},PetscDS,Ptr{Float64},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFEIntegrateBdResidual(arg0::Type{Float64},arg1::PetscFE,arg2::PetscDS,arg3::Integer,arg4::Integer,arg5::Union(Ptr{PetscFECellGeom},StridedArray{PetscFECellGeom},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg8::PetscDS,arg9::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg10::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:PetscFEIntegrateBdResidual,petsc1),PetscErrorCode,(PetscFE,PetscDS,PetscInt,PetscInt,Ptr{PetscFECellGeom},Ptr{Float64},Ptr{Float64},PetscDS,Ptr{Float64},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFEIntegrateJacobian(arg0::Type{Float64},arg1::PetscFE,arg2::PetscDS,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{PetscFECellGeom},StridedArray{PetscFECellGeom},Ptr{Void}),arg7::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg8::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg9::PetscDS,arg10::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg11::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:PetscFEIntegrateJacobian,petsc1),PetscErrorCode,(PetscFE,PetscDS,PetscInt,PetscInt,PetscInt,Ptr{PetscFECellGeom},Ptr{Float64},Ptr{Float64},PetscDS,Ptr{Float64},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFEIntegrateBdJacobian(arg0::Type{Float64},arg1::PetscFE,arg2::PetscDS,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{PetscFECellGeom},StridedArray{PetscFECellGeom},Ptr{Void}),arg7::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg8::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg9::PetscDS,arg10::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg11::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:PetscFEIntegrateBdJacobian,petsc1),PetscErrorCode,(PetscFE,PetscDS,PetscInt,PetscInt,PetscInt,Ptr{PetscFECellGeom},Ptr{Float64},Ptr{Float64},PetscDS,Ptr{Float64},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFECompositeGetMapping(arg0::Type{Float64},arg1::PetscFE,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg4::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg5::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}))
    ccall((:PetscFECompositeGetMapping,petsc1),PetscErrorCode,(PetscFE,Ptr{PetscInt},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFEOpenCLSetRealType(arg0::Type{Float64},arg1::PetscFE,arg2::PetscDataType)
    ccall((:PetscFEOpenCLSetRealType,petsc1),PetscErrorCode,(PetscFE,PetscDataType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFEOpenCLGetRealType(arg0::Type{Float64},arg1::PetscFE,arg2::Union(Ptr{PetscDataType},StridedArray{PetscDataType},Ptr{Void}))
    ccall((:PetscFEOpenCLGetRealType,petsc1),PetscErrorCode,(PetscFE,Ptr{PetscDataType}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetInterpolationType(arg0::Type{Float64},arg1::DM,arg2::DMDAInterpolationType)
    ccall((:DMDASetInterpolationType,petsc1),PetscErrorCode,(DM,DMDAInterpolationType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetInterpolationType(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{DMDAInterpolationType},StridedArray{DMDAInterpolationType},Ptr{Void}))
    ccall((:DMDAGetInterpolationType,petsc1),PetscErrorCode,(DM,Ptr{DMDAInterpolationType}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetElementType(arg0::Type{Float64},arg1::DM,arg2::DMDAElementType)
    ccall((:DMDASetElementType,petsc1),PetscErrorCode,(DM,DMDAElementType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetElementType(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{DMDAElementType},StridedArray{DMDAElementType},Ptr{Void}))
    ccall((:DMDAGetElementType,petsc1),PetscErrorCode,(DM,Ptr{DMDAElementType}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetElements(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMDAGetElements,petsc1),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMDARestoreElements(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMDARestoreElements,petsc1),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMDACreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMDACreate,petsc1),PetscErrorCode,(comm_type,Ptr{DM}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetSizes(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:DMDASetSizes,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMDACreate1d(arg0::Type{Float64},arg1::MPI_Comm,arg2::DMBoundaryType,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMDACreate1d,petsc1),PetscErrorCode,(comm_type,DMBoundaryType,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{DM}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function DMDACreate2d(arg0::Type{Float64},arg1::MPI_Comm,arg2::DMBoundaryType,arg3::DMBoundaryType,arg4::DMDAStencilType,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Integer,arg9::Integer,arg10::Integer,arg11::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg12::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg13::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMDACreate2d,petsc1),PetscErrorCode,(comm_type,DMBoundaryType,DMBoundaryType,DMDAStencilType,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{DM}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13)
end 
=#
#= skipping function with undefined symbols: 
 function DMDACreate3d(arg0::Type{Float64},arg1::MPI_Comm,arg2::DMBoundaryType,arg3::DMBoundaryType,arg4::DMBoundaryType,arg5::DMDAStencilType,arg6::Integer,arg7::Integer,arg8::Integer,arg9::Integer,arg10::Integer,arg11::Integer,arg12::Integer,arg13::Integer,arg14::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg15::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg16::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg17::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMDACreate3d,petsc1),PetscErrorCode,(comm_type,DMBoundaryType,DMBoundaryType,DMBoundaryType,DMDAStencilType,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{DM}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGlobalToNaturalBegin(arg1::DM,arg2::Vec{Float64},arg3::InsertMode,arg4::Vec{Float64})
    ccall((:DMDAGlobalToNaturalBegin,petsc1),PetscErrorCode,(DM,Vec{Float64},InsertMode,Vec{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGlobalToNaturalEnd(arg1::DM,arg2::Vec{Float64},arg3::InsertMode,arg4::Vec{Float64})
    ccall((:DMDAGlobalToNaturalEnd,petsc1),PetscErrorCode,(DM,Vec{Float64},InsertMode,Vec{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMDANaturalToGlobalBegin(arg1::DM,arg2::Vec{Float64},arg3::InsertMode,arg4::Vec{Float64})
    ccall((:DMDANaturalToGlobalBegin,petsc1),PetscErrorCode,(DM,Vec{Float64},InsertMode,Vec{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMDANaturalToGlobalEnd(arg1::DM,arg2::Vec{Float64},arg3::InsertMode,arg4::Vec{Float64})
    ccall((:DMDANaturalToGlobalEnd,petsc1),PetscErrorCode,(DM,Vec{Float64},InsertMode,Vec{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetCorners(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMDAGetCorners,petsc1),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetGhostCorners(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMDAGetGhostCorners,petsc1),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetInfo(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg8::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg9::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg10::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg11::Union(Ptr{DMBoundaryType},StridedArray{DMBoundaryType},Ptr{Void}),arg12::Union(Ptr{DMBoundaryType},StridedArray{DMBoundaryType},Ptr{Void}),arg13::Union(Ptr{DMBoundaryType},StridedArray{DMBoundaryType},Ptr{Void}),arg14::Union(Ptr{DMDAStencilType},StridedArray{DMDAStencilType},Ptr{Void}))
    ccall((:DMDAGetInfo,petsc1),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{DMBoundaryType},Ptr{DMBoundaryType},Ptr{DMBoundaryType},Ptr{DMDAStencilType}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetProcessorSubset(arg0::Type{Float64},arg1::DM,arg2::DMDADirection,arg3::Integer,arg4::Union(Ptr{MPI_Comm},StridedArray{MPI_Comm},Ptr{Void}))
    ccall((:DMDAGetProcessorSubset,petsc1),PetscErrorCode,(DM,DMDADirection,PetscInt,Ptr{comm_type}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetProcessorSubsets(arg0::Type{Float64},arg1::DM,arg2::DMDADirection,arg3::Union(Ptr{MPI_Comm},StridedArray{MPI_Comm},Ptr{Void}))
    ccall((:DMDAGetProcessorSubsets,petsc1),PetscErrorCode,(DM,DMDADirection,Ptr{comm_type}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetRay(arg1::DM,arg2::DMDADirection,arg3::Integer,arg4::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg5::Union(Ptr{VecScatter{Float64}},StridedArray{VecScatter{Float64}},Ptr{Void}))
    ccall((:DMDAGetRay,petsc1),PetscErrorCode,(DM,DMDADirection,PetscInt,Ptr{Vec{Float64}},Ptr{VecScatter{Float64}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGlobalToNaturalAllCreate(arg1::DM,arg2::Union(Ptr{VecScatter{Float64}},StridedArray{VecScatter{Float64}},Ptr{Void}))
    ccall((:DMDAGlobalToNaturalAllCreate,petsc1),PetscErrorCode,(DM,Ptr{VecScatter{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDANaturalAllToGlobalCreate(arg1::DM,arg2::Union(Ptr{VecScatter{Float64}},StridedArray{VecScatter{Float64}},Ptr{Void}))
    ccall((:DMDANaturalAllToGlobalCreate,petsc1),PetscErrorCode,(DM,Ptr{VecScatter{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetScatter(arg1::DM,arg2::Union(Ptr{VecScatter{Float64}},StridedArray{VecScatter{Float64}},Ptr{Void}),arg3::Union(Ptr{VecScatter{Float64}},StridedArray{VecScatter{Float64}},Ptr{Void}))
    ccall((:DMDAGetScatter,petsc1),PetscErrorCode,(DM,Ptr{VecScatter{Float64}},Ptr{VecScatter{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetNeighbors(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Ptr{PetscMPIInt}},StridedArray{Ptr{PetscMPIInt}},Ptr{Void}))
    ccall((:DMDAGetNeighbors,petsc1),PetscErrorCode,(DM,Ptr{Ptr{PetscMPIInt}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetAOType(arg0::Type{Float64},arg1::DM,arg2::AOType)
    ccall((:DMDASetAOType,petsc1),PetscErrorCode,(DM,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetAO(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:DMDAGetAO,petsc1),PetscErrorCode,(DM,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetUniformCoordinates(arg0::Type{Float64},arg1::DM,PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer)
    ccall((:DMDASetUniformCoordinates,petsc1),PetscErrorCode,(DM,Cint,Cint,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetCoordinateArray(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMDAGetCoordinateArray,petsc1),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDARestoreCoordinateArray(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMDARestoreCoordinateArray,petsc1),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetBoundingBox(arg0::Type{Float64},arg1::DM,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:DMDAGetBoundingBox,petsc1),PetscErrorCode,(DM,Ptr{Cint},Ptr{Cint}),arg1,PetscReal,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetLocalBoundingBox(arg0::Type{Float64},arg1::DM,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:DMDAGetLocalBoundingBox,petsc1),PetscErrorCode,(DM,Ptr{Cint},Ptr{Cint}),arg1,PetscReal,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetLogicalCoordinate(arg0::Type{Float64},arg1::DM,arg2::Float64,arg3::Float64,arg4::Float64,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg8::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg9::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg10::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:DMDAGetLogicalCoordinate,petsc1),PetscErrorCode,(DM,Float64,Float64,Float64,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64},Ptr{Float64},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAMapCoordsToPeriodicDomain(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:DMDAMapCoordsToPeriodicDomain,petsc1),PetscErrorCode,(DM,Ptr{Float64},Ptr{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetReducedDMDA(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMDAGetReducedDMDA,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{DM}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetFieldName(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(ByteString,Symbol))
    ccall((:DMDASetFieldName,petsc1),PetscErrorCode,(DM,PetscInt,Cstring),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetFieldName(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:DMDAGetFieldName,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{Ptr{Uint8}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetFieldNames(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:DMDASetFieldNames,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetFieldNames(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Ptr{Ptr{Uint8}}},StridedArray{Ptr{Ptr{Uint8}}},Ptr{Void}))
    ccall((:DMDAGetFieldNames,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Ptr{Uint8}}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetCoordinateName(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(ByteString,Symbol))
    ccall((:DMDASetCoordinateName,petsc1),PetscErrorCode,(DM,PetscInt,Cstring),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetCoordinateName(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:DMDAGetCoordinateName,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{Ptr{Uint8}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetBoundaryType(arg0::Type{Float64},arg1::DM,arg2::DMBoundaryType,arg3::DMBoundaryType,arg4::DMBoundaryType)
    ccall((:DMDASetBoundaryType,petsc1),PetscErrorCode,(DM,DMBoundaryType,DMBoundaryType,DMBoundaryType),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetDof(arg0::Type{Float64},arg1::DM,arg2::Integer)
    ccall((:DMDASetDof,petsc1),PetscErrorCode,(DM,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetOverlap(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:DMDASetOverlap,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetOverlap(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMDAGetOverlap,petsc1),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetNumLocalSubDomains(arg0::Type{Float64},arg1::DM,arg2::Integer)
    ccall((:DMDASetNumLocalSubDomains,petsc1),PetscErrorCode,(DM,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetNumLocalSubDomains(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMDAGetNumLocalSubDomains,petsc1),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetOffset(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMDAGetOffset,petsc1),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetOffset(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer)
    ccall((:DMDASetOffset,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetNonOverlappingRegion(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMDAGetNonOverlappingRegion,petsc1),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetNonOverlappingRegion(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer)
    ccall((:DMDASetNonOverlappingRegion,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetStencilWidth(arg0::Type{Float64},arg1::DM,arg2::Integer)
    ccall((:DMDASetStencilWidth,petsc1),PetscErrorCode,(DM,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetOwnershipRanges(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMDASetOwnershipRanges,petsc1),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetOwnershipRanges(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg3::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg4::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMDAGetOwnershipRanges,petsc1),PetscErrorCode,(DM,Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetNumProcs(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:DMDASetNumProcs,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetStencilType(arg0::Type{Float64},arg1::DM,arg2::DMDAStencilType)
    ccall((:DMDASetStencilType,petsc1),PetscErrorCode,(DM,DMDAStencilType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAVecGetArray(arg1::DM,arg2::Vec{Float64},arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMDAVecGetArray,petsc1),PetscErrorCode,(DM,Vec{Float64},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAVecRestoreArray(arg1::DM,arg2::Vec{Float64},arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMDAVecRestoreArray,petsc1),PetscErrorCode,(DM,Vec{Float64},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAVecGetArrayDOF(arg1::DM,arg2::Vec{Float64},arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMDAVecGetArrayDOF,petsc1),PetscErrorCode,(DM,Vec{Float64},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAVecRestoreArrayDOF(arg1::DM,arg2::Vec{Float64},arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMDAVecRestoreArrayDOF,petsc1),PetscErrorCode,(DM,Vec{Float64},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAVecGetArrayRead(arg1::DM,arg2::Vec{Float64},arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMDAVecGetArrayRead,petsc1),PetscErrorCode,(DM,Vec{Float64},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAVecRestoreArrayRead(arg1::DM,arg2::Vec{Float64},arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMDAVecRestoreArrayRead,petsc1),PetscErrorCode,(DM,Vec{Float64},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAVecGetArrayDOFRead(arg1::DM,arg2::Vec{Float64},arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMDAVecGetArrayDOFRead,petsc1),PetscErrorCode,(DM,Vec{Float64},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAVecRestoreArrayDOFRead(arg1::DM,arg2::Vec{Float64},arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMDAVecRestoreArrayDOFRead,petsc1),PetscErrorCode,(DM,Vec{Float64},Ptr{Void}),arg1,arg2,arg3)
end 
=#
function DMDASplitComm2d(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{MPI_Comm},StridedArray{MPI_Comm},Ptr{Void}))
    err = ccall((:DMDASplitComm2d,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,Ptr{comm_type}),arg1.val,arg2,arg3,arg4,arg5)
    return err
end

#= skipping function with undefined symbols: 
 function DMDACreatePatchIS(arg1::DM,arg2::Union(Ptr{MatStencil},StridedArray{MatStencil},Ptr{Void}),arg3::Union(Ptr{MatStencil},StridedArray{MatStencil},Ptr{Void}),arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:DMDACreatePatchIS,petsc1),PetscErrorCode,(DM,Ptr{MatStencil},Ptr{MatStencil},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetLocalInfo(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{DMDALocalInfo},StridedArray{DMDALocalInfo},Ptr{Void}))
    ccall((:DMDAGetLocalInfo,petsc1),PetscErrorCode,(DM,Ptr{DMDALocalInfo}),arg1,arg2)
end 
=#
function MatRegisterDAAD(arg0::Type{Float64})
    err = ccall((:MatRegisterDAAD,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function MatCreateDAAD(arg1::DM,arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateDAAD,petsc1),PetscErrorCode,(DM,Ptr{Mat{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatCreateSeqUSFFT(arg1::Vec{Float64},arg2::DM,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateSeqUSFFT,petsc1),PetscErrorCode,(Vec{Float64},DM,Ptr{Mat{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetGetMatrix(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMDASetGetMatrix,petsc1),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetBlockFills(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMDASetBlockFills,petsc1),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetRefinementFactor(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:DMDASetRefinementFactor,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetRefinementFactor(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMDAGetRefinementFactor,petsc1),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetArray(arg0::Type{Float64},arg1::DM,arg2::PetscBool,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMDAGetArray,petsc1),PetscErrorCode,(DM,PetscBool,Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDARestoreArray(arg0::Type{Float64},arg1::DM,arg2::PetscBool,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMDARestoreArray,petsc1),PetscErrorCode,(DM,PetscBool,Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDACreatePF(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PF},StridedArray{PF},Ptr{Void}))
    ccall((:DMDACreatePF,petsc1),PetscErrorCode,(DM,Ptr{PF}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetNumCells(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMDAGetNumCells,petsc1),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetCellPoint(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMDAGetCellPoint,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetNumVertices(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMDAGetNumVertices,petsc1),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetNumFaces(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMDAGetNumFaces,petsc1),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetHeightStratum(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMDAGetHeightStratum,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetDepthStratum(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMDAGetDepthStratum,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMDACreateSection(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}))
    ccall((:DMDACreateSection,petsc1),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscSection}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAComputeCellGeometryFEM(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::PetscQuadrature,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:DMDAComputeCellGeometryFEM,petsc1),PetscErrorCode,(DM,PetscInt,PetscQuadrature,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,PetscReal,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetTransitiveClosure(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::PetscBool,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMDAGetTransitiveClosure,petsc1),PetscErrorCode,(DM,PetscInt,PetscBool,Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMDARestoreTransitiveClosure(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::PetscBool,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMDARestoreTransitiveClosure,petsc1),PetscErrorCode,(DM,PetscInt,PetscBool,Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAVecGetClosure(arg1::DM,arg2::PetscSection,arg3::Vec{Float64},arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:DMDAVecGetClosure,petsc1),PetscErrorCode,(DM,PetscSection,Vec{Float64},PetscInt,Ptr{PetscInt},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAVecRestoreClosure(arg1::DM,arg2::PetscSection,arg3::Vec{Float64},arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:DMDAVecRestoreClosure,petsc1),PetscErrorCode,(DM,PetscSection,Vec{Float64},PetscInt,Ptr{PetscInt},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAVecSetClosure(arg1::DM,arg2::PetscSection,arg3::Vec{Float64},arg4::Integer,arg5::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg6::InsertMode)
    ccall((:DMDAVecSetClosure,petsc1),PetscErrorCode,(DM,PetscSection,Vec{Float64},PetscInt,Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetClosure(arg0::Type{Float64},arg1::DM,arg2::PetscSection,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMDAGetClosure,petsc1),PetscErrorCode,(DM,PetscSection,PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMDARestoreClosure(arg0::Type{Float64},arg1::DM,arg2::PetscSection,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMDARestoreClosure,petsc1),PetscErrorCode,(DM,PetscSection,PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetClosureScalar(arg0::Type{Float64},arg1::DM,arg2::PetscSection,arg3::Integer,arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:DMDAGetClosureScalar,petsc1),PetscErrorCode,(DM,PetscSection,PetscInt,Ptr{Float64},Ptr{PetscInt},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function DMDARestoreClosureScalar(arg0::Type{Float64},arg1::DM,arg2::PetscSection,arg3::Integer,arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:DMDARestoreClosureScalar,petsc1),PetscErrorCode,(DM,PetscSection,PetscInt,Ptr{Float64},Ptr{PetscInt},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetClosureScalar(arg0::Type{Float64},arg1::DM,arg2::PetscSection,arg3::Integer,arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg6::InsertMode)
    ccall((:DMDASetClosureScalar,petsc1),PetscErrorCode,(DM,PetscSection,PetscInt,Ptr{Float64},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAConvertToCell(arg0::Type{Float64},arg1::DM,arg2::MatStencil,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMDAConvertToCell,petsc1),PetscErrorCode,(DM,MatStencil,Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetVertexCoordinates(arg0::Type{Float64},arg1::DM,PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer)
    ccall((:DMDASetVertexCoordinates,petsc1),PetscErrorCode,(DM,Cint,Cint,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASetPreallocationCenterDimension(arg0::Type{Float64},arg1::DM,arg2::Integer)
    ccall((:DMDASetPreallocationCenterDimension,petsc1),PetscErrorCode,(DM,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAGetPreallocationCenterDimension(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMDAGetPreallocationCenterDimension,petsc1),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAProjectFunction(arg1::DM,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg4::InsertMode,arg5::Vec{Float64})
    ccall((:DMDAProjectFunction,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}},InsertMode,Vec{Float64}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAProjectFunctionLocal(arg1::DM,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg4::InsertMode,arg5::Vec{Float64})
    ccall((:DMDAProjectFunctionLocal,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}},InsertMode,Vec{Float64}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAComputeL2Diff(arg1::DM,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg4::Vec{Float64},arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:DMDAComputeL2Diff,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}},Vec{Float64},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMDAComputeL2GradientDiff(arg1::DM,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg4::Vec{Float64},PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:DMDAComputeL2GradientDiff,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}},Vec{Float64},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,PetscReal,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMCompositeCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMCompositeCreate,petsc1),PetscErrorCode,(comm_type,Ptr{DM}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMCompositeAddDM(arg0::Type{Float64},arg1::DM,arg2::DM)
    ccall((:DMCompositeAddDM,petsc1),PetscErrorCode,(DM,DM),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMCompositeSetCoupling(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMCompositeSetCoupling,petsc1),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMCompositeAddVecScatter(arg1::DM,arg2::VecScatter{Float64})
    ccall((:DMCompositeAddVecScatter,petsc1),PetscErrorCode,(DM,VecScatter{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMCompositeScatterArray(arg1::DM,arg2::Vec{Float64},arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMCompositeScatterArray,petsc1),PetscErrorCode,(DM,Vec{Float64},Ptr{Vec{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMCompositeGatherArray(arg1::DM,arg2::Vec{Float64},arg3::InsertMode,arg4::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMCompositeGatherArray,petsc1),PetscErrorCode,(DM,Vec{Float64},InsertMode,Ptr{Vec{Float64}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMCompositeGetNumberDM(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMCompositeGetNumberDM,petsc1),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMCompositeGetAccessArray(arg1::DM,arg2::Vec{Float64},arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMCompositeGetAccessArray,petsc1),PetscErrorCode,(DM,Vec{Float64},PetscInt,Ptr{PetscInt},Ptr{Vec{Float64}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMCompositeRestoreAccessArray(arg1::DM,arg2::Vec{Float64},arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMCompositeRestoreAccessArray,petsc1),PetscErrorCode,(DM,Vec{Float64},PetscInt,Ptr{PetscInt},Ptr{Vec{Float64}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMCompositeGetEntriesArray(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMCompositeGetEntriesArray,petsc1),PetscErrorCode,(DM,Ptr{DM}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMCompositeGetGlobalISs(arg1::DM,arg2::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    ccall((:DMCompositeGetGlobalISs,petsc1),PetscErrorCode,(DM,Ptr{Ptr{IS{Float64}}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMCompositeGetLocalISs(arg1::DM,arg2::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    ccall((:DMCompositeGetLocalISs,petsc1),PetscErrorCode,(DM,Ptr{Ptr{IS{Float64}}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMCompositeGetISLocalToGlobalMappings(arg1::DM,arg2::Union(Ptr{Ptr{ISLocalToGlobalMapping{Float64}}},StridedArray{Ptr{ISLocalToGlobalMapping{Float64}}},Ptr{Void}))
    ccall((:DMCompositeGetISLocalToGlobalMappings,petsc1),PetscErrorCode,(DM,Ptr{Ptr{ISLocalToGlobalMapping{Float64}}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPatchCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPatchCreate,petsc1),PetscErrorCode,(comm_type,Ptr{DM}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPatchZoom(arg1::DM,arg2::Vec{Float64},arg3::MatStencil,arg4::MatStencil,arg5::MPI_Comm,arg6::Union(Ptr{DM},StridedArray{DM},Ptr{Void}),arg7::Union(Ptr{PetscSF},StridedArray{PetscSF},Ptr{Void}),arg8::Union(Ptr{PetscSF},StridedArray{PetscSF},Ptr{Void}))
    ccall((:DMPatchZoom,petsc1),PetscErrorCode,(DM,Vec{Float64},MatStencil,MatStencil,comm_type,Ptr{DM},Ptr{PetscSF},Ptr{PetscSF}),arg1,arg2,arg3,arg4,arg5.val,arg6,arg7,arg8)
end 
=#
#= skipping function with undefined symbols: 
 function DMPatchSolve(arg0::Type{Float64},arg1::DM)
    ccall((:DMPatchSolve,petsc1),PetscErrorCode,(DM,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function DMPatchGetPatchSize(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{MatStencil},StridedArray{MatStencil},Ptr{Void}))
    ccall((:DMPatchGetPatchSize,petsc1),PetscErrorCode,(DM,Ptr{MatStencil}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPatchSetPatchSize(arg0::Type{Float64},arg1::DM,arg2::MatStencil)
    ccall((:DMPatchSetPatchSize,petsc1),PetscErrorCode,(DM,MatStencil),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPatchGetCommSize(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{MatStencil},StridedArray{MatStencil},Ptr{Void}))
    ccall((:DMPatchGetCommSize,petsc1),PetscErrorCode,(DM,Ptr{MatStencil}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPatchSetCommSize(arg0::Type{Float64},arg1::DM,arg2::MatStencil)
    ccall((:DMPatchSetCommSize,petsc1),PetscErrorCode,(DM,MatStencil),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPatchGetCoarse(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPatchGetCoarse,petsc1),PetscErrorCode,(DM,Ptr{DM}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPatchCreateGrid(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::MatStencil,arg4::MatStencil,arg5::MatStencil,arg6::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPatchCreateGrid,petsc1),PetscErrorCode,(comm_type,PetscInt,MatStencil,MatStencil,MatStencil,Ptr{DM}),arg1.val,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function PetscLimiterCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscLimiter},StridedArray{PetscLimiter},Ptr{Void}))
    ccall((:PetscLimiterCreate,petsc1),PetscErrorCode,(comm_type,Ptr{PetscLimiter}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscLimiterDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscLimiter},StridedArray{PetscLimiter},Ptr{Void}))
    ccall((:PetscLimiterDestroy,petsc1),PetscErrorCode,(Ptr{PetscLimiter},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscLimiterSetType(arg0::Type{Float64},arg1::PetscLimiter,arg2::PetscLimiterType)
    ccall((:PetscLimiterSetType,petsc1),PetscErrorCode,(PetscLimiter,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscLimiterGetType(arg0::Type{Float64},arg1::PetscLimiter,arg2::Union(Ptr{PetscLimiterType},StridedArray{PetscLimiterType},Ptr{Void}))
    ccall((:PetscLimiterGetType,petsc1),PetscErrorCode,(PetscLimiter,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscLimiterSetUp(arg0::Type{Float64},arg1::PetscLimiter)
    ccall((:PetscLimiterSetUp,petsc1),PetscErrorCode,(PetscLimiter,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscLimiterSetFromOptions(arg0::Type{Float64},arg1::PetscLimiter)
    ccall((:PetscLimiterSetFromOptions,petsc1),PetscErrorCode,(PetscLimiter,),arg1)
end 
=#
function PetscLimiterRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscLimiterRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

function PetscLimiterRegisterDestroy(arg0::Type{Float64})
    err = ccall((:PetscLimiterRegisterDestroy,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function PetscLimiterLimit(arg0::Type{Float64},arg1::PetscLimiter,PetscReal::Integer,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscLimiterLimit,petsc1),PetscErrorCode,(PetscLimiter,Cint,Ptr{Cint}),arg1,PetscReal,arg2)
end 
=#
function PetscFVInitializePackage(arg0::Type{Float64})
    err = ccall((:PetscFVInitializePackage,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function PetscFVCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscFV},StridedArray{PetscFV},Ptr{Void}))
    ccall((:PetscFVCreate,petsc1),PetscErrorCode,(comm_type,Ptr{PetscFV}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFVDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscFV},StridedArray{PetscFV},Ptr{Void}))
    ccall((:PetscFVDestroy,petsc1),PetscErrorCode,(Ptr{PetscFV},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFVSetType(arg0::Type{Float64},arg1::PetscFV,arg2::PetscFVType)
    ccall((:PetscFVSetType,petsc1),PetscErrorCode,(PetscFV,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFVGetType(arg0::Type{Float64},arg1::PetscFV,arg2::Union(Ptr{PetscFVType},StridedArray{PetscFVType},Ptr{Void}))
    ccall((:PetscFVGetType,petsc1),PetscErrorCode,(PetscFV,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFVSetUp(arg0::Type{Float64},arg1::PetscFV)
    ccall((:PetscFVSetUp,petsc1),PetscErrorCode,(PetscFV,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFVSetFromOptions(arg0::Type{Float64},arg1::PetscFV)
    ccall((:PetscFVSetFromOptions,petsc1),PetscErrorCode,(PetscFV,),arg1)
end 
=#
function PetscFVRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscFVRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

function PetscFVRegisterDestroy(arg0::Type{Float64})
    err = ccall((:PetscFVRegisterDestroy,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function PetscFVSetLimiter(arg0::Type{Float64},arg1::PetscFV,arg2::PetscLimiter)
    ccall((:PetscFVSetLimiter,petsc1),PetscErrorCode,(PetscFV,PetscLimiter),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFVGetLimiter(arg0::Type{Float64},arg1::PetscFV,arg2::Union(Ptr{PetscLimiter},StridedArray{PetscLimiter},Ptr{Void}))
    ccall((:PetscFVGetLimiter,petsc1),PetscErrorCode,(PetscFV,Ptr{PetscLimiter}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFVSetNumComponents(arg0::Type{Float64},arg1::PetscFV,arg2::Integer)
    ccall((:PetscFVSetNumComponents,petsc1),PetscErrorCode,(PetscFV,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFVGetNumComponents(arg0::Type{Float64},arg1::PetscFV,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscFVGetNumComponents,petsc1),PetscErrorCode,(PetscFV,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFVSetSpatialDimension(arg0::Type{Float64},arg1::PetscFV,arg2::Integer)
    ccall((:PetscFVSetSpatialDimension,petsc1),PetscErrorCode,(PetscFV,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFVGetSpatialDimension(arg0::Type{Float64},arg1::PetscFV,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscFVGetSpatialDimension,petsc1),PetscErrorCode,(PetscFV,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFVSetComputeGradients(arg0::Type{Float64},arg1::PetscFV,arg2::PetscBool)
    ccall((:PetscFVSetComputeGradients,petsc1),PetscErrorCode,(PetscFV,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFVGetComputeGradients(arg0::Type{Float64},arg1::PetscFV,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscFVGetComputeGradients,petsc1),PetscErrorCode,(PetscFV,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFVSetQuadrature(arg0::Type{Float64},arg1::PetscFV,arg2::PetscQuadrature)
    ccall((:PetscFVSetQuadrature,petsc1),PetscErrorCode,(PetscFV,PetscQuadrature),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFVGetQuadrature(arg0::Type{Float64},arg1::PetscFV,arg2::Union(Ptr{PetscQuadrature},StridedArray{PetscQuadrature},Ptr{Void}))
    ccall((:PetscFVGetQuadrature,petsc1),PetscErrorCode,(PetscFV,Ptr{PetscQuadrature}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFVSetDualSpace(arg0::Type{Float64},arg1::PetscFV,arg2::PetscDualSpace)
    ccall((:PetscFVSetDualSpace,petsc1),PetscErrorCode,(PetscFV,PetscDualSpace),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFVGetDualSpace(arg0::Type{Float64},arg1::PetscFV,arg2::Union(Ptr{PetscDualSpace},StridedArray{PetscDualSpace},Ptr{Void}))
    ccall((:PetscFVGetDualSpace,petsc1),PetscErrorCode,(PetscFV,Ptr{PetscDualSpace}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFVRefine(arg0::Type{Float64},arg1::PetscFV,arg2::Union(Ptr{PetscFV},StridedArray{PetscFV},Ptr{Void}))
    ccall((:PetscFVRefine,petsc1),PetscErrorCode,(PetscFV,Ptr{PetscFV}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFVGetDefaultTabulation(arg0::Type{Float64},arg1::PetscFV,arg2::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg3::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg4::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}))
    ccall((:PetscFVGetDefaultTabulation,petsc1),PetscErrorCode,(PetscFV,Ptr{Ptr{Cint}},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFVGetTabulation(arg0::Type{Float64},arg1::PetscFV,arg2::Integer,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg4::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg5::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}))
    ccall((:PetscFVGetTabulation,petsc1),PetscErrorCode,(PetscFV,PetscInt,Ptr{Cint},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}}),arg1,arg2,PetscReal,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFVRestoreTabulation(arg0::Type{Float64},arg1::PetscFV,arg2::Integer,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg4::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg5::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}))
    ccall((:PetscFVRestoreTabulation,petsc1),PetscErrorCode,(PetscFV,PetscInt,Ptr{Cint},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}}),arg1,arg2,PetscReal,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFVComputeGradient(arg0::Type{Float64},arg1::PetscFV,arg2::Integer,arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:PetscFVComputeGradient,petsc1),PetscErrorCode,(PetscFV,PetscInt,Ptr{Float64},Ptr{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFVIntegrateRHSFunction(arg0::Type{Float64},arg1::PetscFV,arg2::PetscDS,arg3::Integer,arg4::Integer,arg5::Union(Ptr{PetscFVFaceGeom},StridedArray{PetscFVFaceGeom},Ptr{Void}),arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg7::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg8::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg9::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg10::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:PetscFVIntegrateRHSFunction,petsc1),PetscErrorCode,(PetscFV,PetscDS,PetscInt,PetscInt,Ptr{PetscFVFaceGeom},Ptr{Cint},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFVLeastSquaresSetMaxFaces(arg0::Type{Float64},arg1::PetscFV,arg2::Integer)
    ccall((:PetscFVLeastSquaresSetMaxFaces,petsc1),PetscErrorCode,(PetscFV,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscPartitionerCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscPartitioner},StridedArray{PetscPartitioner},Ptr{Void}))
    ccall((:PetscPartitionerCreate,petsc1),PetscErrorCode,(comm_type,Ptr{PetscPartitioner}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscPartitionerDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscPartitioner},StridedArray{PetscPartitioner},Ptr{Void}))
    ccall((:PetscPartitionerDestroy,petsc1),PetscErrorCode,(Ptr{PetscPartitioner},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscPartitionerSetType(arg0::Type{Float64},arg1::PetscPartitioner,arg2::PetscPartitionerType)
    ccall((:PetscPartitionerSetType,petsc1),PetscErrorCode,(PetscPartitioner,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscPartitionerGetType(arg0::Type{Float64},arg1::PetscPartitioner,arg2::Union(Ptr{PetscPartitionerType},StridedArray{PetscPartitionerType},Ptr{Void}))
    ccall((:PetscPartitionerGetType,petsc1),PetscErrorCode,(PetscPartitioner,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscPartitionerSetUp(arg0::Type{Float64},arg1::PetscPartitioner)
    ccall((:PetscPartitionerSetUp,petsc1),PetscErrorCode,(PetscPartitioner,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscPartitionerSetFromOptions(arg0::Type{Float64},arg1::PetscPartitioner)
    ccall((:PetscPartitionerSetFromOptions,petsc1),PetscErrorCode,(PetscPartitioner,),arg1)
end 
=#
function PetscPartitionerRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscPartitionerRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

function PetscPartitionerRegisterDestroy(arg0::Type{Float64})
    err = ccall((:PetscPartitionerRegisterDestroy,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function PetscPartitionerPartition(arg1::PetscPartitioner,arg2::DM,arg3::PetscSection,arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:PetscPartitionerPartition,petsc1),PetscErrorCode,(PetscPartitioner,DM,PetscSection,Ptr{IS{Float64}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscPartitionerShellSetPartition(arg0::Type{Float64},arg1::PetscPartitioner,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscPartitionerShellSetPartition,petsc1),PetscErrorCode,(PetscPartitioner,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexCreate,petsc1),PetscErrorCode,(comm_type,Ptr{DM}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateCohesiveSubmesh(arg0::Type{Float64},arg1::DM,arg2::PetscBool,arg3::Union(ByteString,Symbol),arg4::Integer,arg5::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexCreateCohesiveSubmesh,petsc1),PetscErrorCode,(DM,PetscBool,Cstring,PetscInt,Ptr{DM}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateFromCellList(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::PetscBool,arg7::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg8::Integer,arg9::Union(Ptr{Cdouble},StridedArray{Cdouble},Ptr{Void}),arg10::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexCreateFromCellList,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,PetscInt,PetscBool,Ptr{Cint},PetscInt,Ptr{Cdouble},Ptr{DM}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateFromDAG(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:DMPlexCreateFromDAG,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateReferenceCell(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::PetscBool,arg4::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexCreateReferenceCell,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscBool,Ptr{DM}),arg1.val,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetChart(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexGetChart,petsc1),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetChart(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer)
    ccall((:DMPlexSetChart,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetConeSize(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexGetConeSize,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetConeSize(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer)
    ccall((:DMPlexSetConeSize,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexAddConeSize(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer)
    ccall((:DMPlexAddConeSize,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetCone(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMPlexGetCone,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{Ptr{PetscInt}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetCone(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexSetCone,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexInsertCone(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:DMPlexInsertCone,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexInsertConeOrientation(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:DMPlexInsertConeOrientation,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetConeOrientation(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMPlexGetConeOrientation,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{Ptr{PetscInt}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetConeOrientation(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexSetConeOrientation,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetSupportSize(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexGetSupportSize,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetSupportSize(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer)
    ccall((:DMPlexSetSupportSize,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetSupport(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMPlexGetSupport,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{Ptr{PetscInt}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetSupport(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexSetSupport,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexInsertSupport(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:DMPlexInsertSupport,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetConeSection(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}))
    ccall((:DMPlexGetConeSection,petsc1),PetscErrorCode,(DM,Ptr{PetscSection}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetSupportSection(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}))
    ccall((:DMPlexGetSupportSection,petsc1),PetscErrorCode,(DM,Ptr{PetscSection}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetCones(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMPlexGetCones,petsc1),PetscErrorCode,(DM,Ptr{Ptr{PetscInt}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetConeOrientations(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMPlexGetConeOrientations,petsc1),PetscErrorCode,(DM,Ptr{Ptr{PetscInt}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetMaxSizes(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexGetMaxSizes,petsc1),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSymmetrize(arg0::Type{Float64},arg1::DM)
    ccall((:DMPlexSymmetrize,petsc1),PetscErrorCode,(DM,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexStratify(arg0::Type{Float64},arg1::DM)
    ccall((:DMPlexStratify,petsc1),PetscErrorCode,(DM,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexEqual(arg0::Type{Float64},arg1::DM,arg2::DM,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:DMPlexEqual,petsc1),PetscErrorCode,(DM,DM,Ptr{PetscBool}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexReverseCell(arg0::Type{Float64},arg1::DM,arg2::Integer)
    ccall((:DMPlexReverseCell,petsc1),PetscErrorCode,(DM,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexOrient(arg0::Type{Float64},arg1::DM)
    ccall((:DMPlexOrient,petsc1),PetscErrorCode,(DM,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexInterpolate(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexInterpolate,petsc1),PetscErrorCode,(DM,Ptr{DM}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexUninterpolate(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexUninterpolate,petsc1),PetscErrorCode,(DM,Ptr{DM}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexLoad(arg1::PetscViewer{Float64},arg2::DM)
    ccall((:DMPlexLoad,petsc1),PetscErrorCode,(PetscViewer{Float64},DM),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexPreallocateOperator(arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Mat{Float64},arg8::PetscBool)
    ccall((:DMPlexPreallocateOperator,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Mat{Float64},PetscBool),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetPointLocal(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexGetPointLocal,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexPointLocalRead(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMPlexPointLocalRead,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{Float64},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexPointLocalRef(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMPlexPointLocalRef,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{Float64},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetPointLocalField(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexGetPointLocalField,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexPointLocalFieldRef(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMPlexPointLocalFieldRef,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,Ptr{Float64},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexPointLocalFieldRead(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMPlexPointLocalFieldRead,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,Ptr{Float64},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetPointGlobal(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexGetPointGlobal,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexPointGlobalRead(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMPlexPointGlobalRead,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{Float64},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexPointGlobalRef(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMPlexPointGlobalRef,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{Float64},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetPointGlobalField(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexGetPointGlobalField,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexPointGlobalFieldRef(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMPlexPointGlobalFieldRef,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,Ptr{Float64},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexPointGlobalFieldRead(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMPlexPointGlobalFieldRead,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,Ptr{Float64},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMLabelCreate(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{DMLabel},StridedArray{DMLabel},Ptr{Void}))
    ccall((:DMLabelCreate,petsc1),PetscErrorCode,(Cstring,Ptr{DMLabel}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMLabelView(arg1::DMLabel,arg2::PetscViewer{Float64})
    ccall((:DMLabelView,petsc1),PetscErrorCode,(DMLabel,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMLabelDestroy(arg0::Type{Float64},arg1::Union(Ptr{DMLabel},StridedArray{DMLabel},Ptr{Void}))
    ccall((:DMLabelDestroy,petsc1),PetscErrorCode,(Ptr{DMLabel},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function DMLabelDuplicate(arg0::Type{Float64},arg1::DMLabel,arg2::Union(Ptr{DMLabel},StridedArray{DMLabel},Ptr{Void}))
    ccall((:DMLabelDuplicate,petsc1),PetscErrorCode,(DMLabel,Ptr{DMLabel}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMLabelGetName(arg0::Type{Float64},arg1::DMLabel,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:DMLabelGetName,petsc1),PetscErrorCode,(DMLabel,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMLabelGetValue(arg0::Type{Float64},arg1::DMLabel,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMLabelGetValue,petsc1),PetscErrorCode,(DMLabel,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMLabelSetValue(arg0::Type{Float64},arg1::DMLabel,arg2::Integer,arg3::Integer)
    ccall((:DMLabelSetValue,petsc1),PetscErrorCode,(DMLabel,PetscInt,PetscInt),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMLabelClearValue(arg0::Type{Float64},arg1::DMLabel,arg2::Integer,arg3::Integer)
    ccall((:DMLabelClearValue,petsc1),PetscErrorCode,(DMLabel,PetscInt,PetscInt),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMLabelInsertIS(arg1::DMLabel,arg2::IS{Float64},arg3::Integer)
    ccall((:DMLabelInsertIS,petsc1),PetscErrorCode,(DMLabel,IS{Float64},PetscInt),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMLabelGetNumValues(arg0::Type{Float64},arg1::DMLabel,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMLabelGetNumValues,petsc1),PetscErrorCode,(DMLabel,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMLabelGetStratumBounds(arg0::Type{Float64},arg1::DMLabel,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMLabelGetStratumBounds,petsc1),PetscErrorCode,(DMLabel,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMLabelGetValueIS(arg1::DMLabel,arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:DMLabelGetValueIS,petsc1),PetscErrorCode,(DMLabel,Ptr{IS{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMLabelStratumHasPoint(arg0::Type{Float64},arg1::DMLabel,arg2::Integer,arg3::Integer,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:DMLabelStratumHasPoint,petsc1),PetscErrorCode,(DMLabel,PetscInt,PetscInt,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMLabelGetStratumSize(arg0::Type{Float64},arg1::DMLabel,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMLabelGetStratumSize,petsc1),PetscErrorCode,(DMLabel,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMLabelGetStratumIS(arg1::DMLabel,arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:DMLabelGetStratumIS,petsc1),PetscErrorCode,(DMLabel,PetscInt,Ptr{IS{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMLabelClearStratum(arg0::Type{Float64},arg1::DMLabel,arg2::Integer)
    ccall((:DMLabelClearStratum,petsc1),PetscErrorCode,(DMLabel,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMLabelCreateIndex(arg0::Type{Float64},arg1::DMLabel,arg2::Integer,arg3::Integer)
    ccall((:DMLabelCreateIndex,petsc1),PetscErrorCode,(DMLabel,PetscInt,PetscInt),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMLabelDestroyIndex(arg0::Type{Float64},arg1::DMLabel)
    ccall((:DMLabelDestroyIndex,petsc1),PetscErrorCode,(DMLabel,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function DMLabelHasValue(arg0::Type{Float64},arg1::DMLabel,arg2::Integer,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:DMLabelHasValue,petsc1),PetscErrorCode,(DMLabel,PetscInt,Ptr{PetscBool}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMLabelHasPoint(arg0::Type{Float64},arg1::DMLabel,arg2::Integer,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:DMLabelHasPoint,petsc1),PetscErrorCode,(DMLabel,PetscInt,Ptr{PetscBool}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMLabelFilter(arg0::Type{Float64},arg1::DMLabel,arg2::Integer,arg3::Integer)
    ccall((:DMLabelFilter,petsc1),PetscErrorCode,(DMLabel,PetscInt,PetscInt),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMLabelPermute(arg1::DMLabel,arg2::IS{Float64},arg3::Union(Ptr{DMLabel},StridedArray{DMLabel},Ptr{Void}))
    ccall((:DMLabelPermute,petsc1),PetscErrorCode,(DMLabel,IS{Float64},Ptr{DMLabel}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMLabelDistribute(arg0::Type{Float64},arg1::DMLabel,arg2::PetscSF,arg3::Union(Ptr{DMLabel},StridedArray{DMLabel},Ptr{Void}))
    ccall((:DMLabelDistribute,petsc1),PetscErrorCode,(DMLabel,PetscSF,Ptr{DMLabel}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateLabel(arg0::Type{Float64},arg1::DM,arg2::Union(ByteString,Symbol))
    ccall((:DMPlexCreateLabel,petsc1),PetscErrorCode,(DM,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetLabelValue(arg0::Type{Float64},arg1::DM,arg2::Union(ByteString,Symbol),arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexGetLabelValue,petsc1),PetscErrorCode,(DM,Cstring,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetLabelValue(arg0::Type{Float64},arg1::DM,arg2::Union(ByteString,Symbol),arg3::Integer,arg4::Integer)
    ccall((:DMPlexSetLabelValue,petsc1),PetscErrorCode,(DM,Cstring,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexClearLabelValue(arg0::Type{Float64},arg1::DM,arg2::Union(ByteString,Symbol),arg3::Integer,arg4::Integer)
    ccall((:DMPlexClearLabelValue,petsc1),PetscErrorCode,(DM,Cstring,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetLabelSize(arg0::Type{Float64},arg1::DM,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexGetLabelSize,petsc1),PetscErrorCode,(DM,Cstring,Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetLabelIdIS(arg1::DM,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:DMPlexGetLabelIdIS,petsc1),PetscErrorCode,(DM,Cstring,Ptr{IS{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetStratumSize(arg0::Type{Float64},arg1::DM,arg2::Union(ByteString,Symbol),arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexGetStratumSize,petsc1),PetscErrorCode,(DM,Cstring,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetStratumIS(arg1::DM,arg2::Union(ByteString,Symbol),arg3::Integer,arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:DMPlexGetStratumIS,petsc1),PetscErrorCode,(DM,Cstring,PetscInt,Ptr{IS{Float64}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexClearLabelStratum(arg0::Type{Float64},arg1::DM,arg2::Union(ByteString,Symbol),arg3::Integer)
    ccall((:DMPlexClearLabelStratum,petsc1),PetscErrorCode,(DM,Cstring,PetscInt),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetLabelOutput(arg0::Type{Float64},arg1::DM,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:DMPlexGetLabelOutput,petsc1),PetscErrorCode,(DM,Cstring,Ptr{PetscBool}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetLabelOutput(arg0::Type{Float64},arg1::DM,arg2::Union(ByteString,Symbol),arg3::PetscBool)
    ccall((:DMPlexSetLabelOutput,petsc1),PetscErrorCode,(DM,Cstring,PetscBool),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionCreateGlobalSectionLabel(arg0::Type{Float64},arg1::PetscSection,arg2::PetscSF,arg3::PetscBool,arg4::DMLabel,arg5::Integer,arg6::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}))
    ccall((:PetscSectionCreateGlobalSectionLabel,petsc1),PetscErrorCode,(PetscSection,PetscSF,PetscBool,DMLabel,PetscInt,Ptr{PetscSection}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetNumLabels(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexGetNumLabels,petsc1),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetLabelName(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:DMPlexGetLabelName,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{Ptr{Uint8}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexHasLabel(arg0::Type{Float64},arg1::DM,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:DMPlexHasLabel,petsc1),PetscErrorCode,(DM,Cstring,Ptr{PetscBool}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetLabel(arg0::Type{Float64},arg1::DM,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{DMLabel},StridedArray{DMLabel},Ptr{Void}))
    ccall((:DMPlexGetLabel,petsc1),PetscErrorCode,(DM,Cstring,Ptr{DMLabel}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetLabelByNum(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{DMLabel},StridedArray{DMLabel},Ptr{Void}))
    ccall((:DMPlexGetLabelByNum,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{DMLabel}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexAddLabel(arg0::Type{Float64},arg1::DM,arg2::DMLabel)
    ccall((:DMPlexAddLabel,petsc1),PetscErrorCode,(DM,DMLabel),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexRemoveLabel(arg0::Type{Float64},arg1::DM,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{DMLabel},StridedArray{DMLabel},Ptr{Void}))
    ccall((:DMPlexRemoveLabel,petsc1),PetscErrorCode,(DM,Cstring,Ptr{DMLabel}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetCellNumbering(arg1::DM,arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:DMPlexGetCellNumbering,petsc1),PetscErrorCode,(DM,Ptr{IS{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetVertexNumbering(arg1::DM,arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:DMPlexGetVertexNumbering,petsc1),PetscErrorCode,(DM,Ptr{IS{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreatePointNumbering(arg1::DM,arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:DMPlexCreatePointNumbering,petsc1),PetscErrorCode,(DM,Ptr{IS{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetDepth(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexGetDepth,petsc1),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetDepthLabel(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{DMLabel},StridedArray{DMLabel},Ptr{Void}))
    ccall((:DMPlexGetDepthLabel,petsc1),PetscErrorCode,(DM,Ptr{DMLabel}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetDepthStratum(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexGetDepthStratum,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetHeightStratum(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexGetHeightStratum,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetMeet(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMPlexGetMeet,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetFullMeet(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMPlexGetFullMeet,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexRestoreMeet(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMPlexRestoreMeet,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetJoin(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMPlexGetJoin,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetFullJoin(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMPlexGetFullJoin,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexRestoreJoin(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMPlexRestoreJoin,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetTransitiveClosure(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::PetscBool,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMPlexGetTransitiveClosure,petsc1),PetscErrorCode,(DM,PetscInt,PetscBool,Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexRestoreTransitiveClosure(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::PetscBool,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMPlexRestoreTransitiveClosure,petsc1),PetscErrorCode,(DM,PetscInt,PetscBool,Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGenerate(arg0::Type{Float64},arg1::DM,arg2::Union(ByteString,Symbol),arg3::PetscBool,arg4::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexGenerate,petsc1),PetscErrorCode,(DM,Cstring,PetscBool,Ptr{DM}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCopyCoordinates(arg0::Type{Float64},arg1::DM,arg2::DM)
    ccall((:DMPlexCopyCoordinates,petsc1),PetscErrorCode,(DM,DM),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCopyLabels(arg0::Type{Float64},arg1::DM,arg2::DM)
    ccall((:DMPlexCopyLabels,petsc1),PetscErrorCode,(DM,DM),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateDoublet(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::PetscBool,arg4::PetscBool,arg5::PetscBool,PetscReal::Integer,arg6::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexCreateDoublet,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscBool,PetscBool,PetscBool,Cint,Ptr{DM}),arg1.val,arg2,arg3,arg4,arg5,PetscReal,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateSquareBoundary(arg0::Type{Float64},arg1::DM,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexCreateSquareBoundary,petsc1),PetscErrorCode,(DM,Ptr{Cint},Ptr{Cint},Ptr{PetscInt}),arg1,PetscReal,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateCubeBoundary(arg0::Type{Float64},arg1::DM,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexCreateCubeBoundary,petsc1),PetscErrorCode,(DM,Ptr{Cint},Ptr{Cint},Ptr{PetscInt}),arg1,PetscReal,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateSquareMesh(arg0::Type{Float64},arg1::DM,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::DMBoundaryType,arg5::DMBoundaryType)
    ccall((:DMPlexCreateSquareMesh,petsc1),PetscErrorCode,(DM,Ptr{Cint},Ptr{Cint},Ptr{PetscInt},DMBoundaryType,DMBoundaryType),arg1,PetscReal,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateBoxMesh(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::PetscBool,arg4::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexCreateBoxMesh,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscBool,Ptr{DM}),arg1.val,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateHexBoxMesh(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::DMBoundaryType,arg5::DMBoundaryType,arg6::DMBoundaryType,arg7::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexCreateHexBoxMesh,petsc1),PetscErrorCode,(comm_type,PetscInt,Ptr{PetscInt},DMBoundaryType,DMBoundaryType,DMBoundaryType,Ptr{DM}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateConeSection(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}))
    ccall((:DMPlexCreateConeSection,petsc1),PetscErrorCode,(DM,Ptr{PetscSection}),arg1,arg2)
end 
=#
function DMPlexInvertCell(arg0::Type{Float64},arg1::Integer,arg2::Integer,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:DMPlexInvertCell,petsc1),PetscErrorCode,(PetscInt,PetscInt,Ptr{Cint}),arg1,arg2,arg3)
    return err
end

#= skipping function with undefined symbols: 
 function DMPlexLocalizeCoordinate(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:DMPlexLocalizeCoordinate,petsc1),PetscErrorCode,(DM,Ptr{Float64},Ptr{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexLocalizeCoordinates(arg0::Type{Float64},arg1::DM)
    ccall((:DMPlexLocalizeCoordinates,petsc1),PetscErrorCode,(DM,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCheckSymmetry(arg0::Type{Float64},arg1::DM)
    ccall((:DMPlexCheckSymmetry,petsc1),PetscErrorCode,(DM,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCheckSkeleton(arg0::Type{Float64},arg1::DM,arg2::PetscBool,arg3::Integer)
    ccall((:DMPlexCheckSkeleton,petsc1),PetscErrorCode,(DM,PetscBool,PetscInt),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCheckFaces(arg0::Type{Float64},arg1::DM,arg2::PetscBool,arg3::Integer)
    ccall((:DMPlexCheckFaces,petsc1),PetscErrorCode,(DM,PetscBool,PetscInt),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexTriangleSetOptions(arg0::Type{Float64},arg1::DM,arg2::Union(ByteString,Symbol))
    ccall((:DMPlexTriangleSetOptions,petsc1),PetscErrorCode,(DM,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexTetgenSetOptions(arg0::Type{Float64},arg1::DM,arg2::Union(ByteString,Symbol))
    ccall((:DMPlexTetgenSetOptions,petsc1),PetscErrorCode,(DM,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateNeighborCSR(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg5::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMPlexCreateNeighborCSR,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetPartitioner(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscPartitioner},StridedArray{PetscPartitioner},Ptr{Void}))
    ccall((:DMPlexGetPartitioner,petsc1),PetscErrorCode,(DM,Ptr{PetscPartitioner}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetPartitioner(arg0::Type{Float64},arg1::DM,arg2::PetscPartitioner)
    ccall((:DMPlexSetPartitioner,petsc1),PetscErrorCode,(DM,PetscPartitioner),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreatePartition(arg1::DM,arg2::Union(ByteString,Symbol),arg3::Integer,arg4::PetscBool,arg5::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}),arg6::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg7::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}),arg8::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:DMPlexCreatePartition,petsc1),PetscErrorCode,(DM,Cstring,PetscInt,PetscBool,Ptr{PetscSection},Ptr{IS{Float64}},Ptr{PetscSection},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreatePartitionerGraph(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg5::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMPlexCreatePartitionerGraph,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreatePartitionClosure(arg1::DM,arg2::PetscSection,arg3::IS{Float64},arg4::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}),arg5::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:DMPlexCreatePartitionClosure,petsc1),PetscErrorCode,(DM,PetscSection,IS{Float64},Ptr{PetscSection},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexPartitionLabelInvert(arg0::Type{Float64},arg1::DM,arg2::DMLabel,arg3::PetscSF,arg4::DMLabel)
    ccall((:DMPlexPartitionLabelInvert,petsc1),PetscErrorCode,(DM,DMLabel,PetscSF,DMLabel),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexPartitionLabelClosure(arg0::Type{Float64},arg1::DM,arg2::DMLabel)
    ccall((:DMPlexPartitionLabelClosure,petsc1),PetscErrorCode,(DM,DMLabel),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexPartitionLabelAdjacency(arg0::Type{Float64},arg1::DM,arg2::DMLabel)
    ccall((:DMPlexPartitionLabelAdjacency,petsc1),PetscErrorCode,(DM,DMLabel),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexPartitionLabelCreateSF(arg0::Type{Float64},arg1::DM,arg2::DMLabel,arg3::Union(Ptr{PetscSF},StridedArray{PetscSF},Ptr{Void}))
    ccall((:DMPlexPartitionLabelCreateSF,petsc1),PetscErrorCode,(DM,DMLabel,Ptr{PetscSF}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexDistribute(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscSF},StridedArray{PetscSF},Ptr{Void}),arg4::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexDistribute,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscSF},Ptr{DM}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexDistributeOverlap(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscSF},StridedArray{PetscSF},Ptr{Void}),arg4::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexDistributeOverlap,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscSF},Ptr{DM}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexDistributeField(arg1::DM,arg2::PetscSF,arg3::PetscSection,arg4::Vec{Float64},arg5::PetscSection,arg6::Vec{Float64})
    ccall((:DMPlexDistributeField,petsc1),PetscErrorCode,(DM,PetscSF,PetscSection,Vec{Float64},PetscSection,Vec{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexDistributeFieldIS(arg1::DM,arg2::PetscSF,arg3::PetscSection,arg4::IS{Float64},arg5::PetscSection,arg6::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:DMPlexDistributeFieldIS,petsc1),PetscErrorCode,(DM,PetscSF,PetscSection,IS{Float64},PetscSection,Ptr{IS{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexDistributeData(arg0::Type{Float64},arg1::DM,arg2::PetscSF,arg3::PetscSection,arg4::MPI_Datatype,arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg6::PetscSection,arg7::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:DMPlexDistributeData,petsc1),PetscErrorCode,(DM,PetscSF,PetscSection,MPI_Datatype,Ptr{Void},PetscSection,Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexMigrate(arg0::Type{Float64},arg1::DM,arg2::PetscSF,arg3::DM)
    ccall((:DMPlexMigrate,petsc1),PetscErrorCode,(DM,PetscSF,DM),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetAdjacencyUseCone(arg0::Type{Float64},arg1::DM,arg2::PetscBool)
    ccall((:DMPlexSetAdjacencyUseCone,petsc1),PetscErrorCode,(DM,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetAdjacencyUseCone(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:DMPlexGetAdjacencyUseCone,petsc1),PetscErrorCode,(DM,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetAdjacencyUseClosure(arg0::Type{Float64},arg1::DM,arg2::PetscBool)
    ccall((:DMPlexSetAdjacencyUseClosure,petsc1),PetscErrorCode,(DM,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetAdjacencyUseClosure(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:DMPlexGetAdjacencyUseClosure,petsc1),PetscErrorCode,(DM,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetAdjacencyUseAnchors(arg0::Type{Float64},arg1::DM,arg2::PetscBool)
    ccall((:DMPlexSetAdjacencyUseAnchors,petsc1),PetscErrorCode,(DM,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetAdjacencyUseAnchors(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:DMPlexGetAdjacencyUseAnchors,petsc1),PetscErrorCode,(DM,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetAdjacency(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMPlexGetAdjacency,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetOrdering(arg1::DM,arg2::MatOrderingType,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:DMPlexGetOrdering,petsc1),PetscErrorCode,(DM,Cstring,Ptr{IS{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexPermute(arg1::DM,arg2::IS{Float64},arg3::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexPermute,petsc1),PetscErrorCode,(DM,IS{Float64},Ptr{DM}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateProcessSF(arg1::DM,arg2::PetscSF,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Union(Ptr{PetscSF},StridedArray{PetscSF},Ptr{Void}))
    ccall((:DMPlexCreateProcessSF,petsc1),PetscErrorCode,(DM,PetscSF,Ptr{IS{Float64}},Ptr{PetscSF}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateTwoSidedProcessSF(arg1::DM,arg2::PetscSF,arg3::PetscSection,arg4::IS{Float64},arg5::PetscSection,arg6::IS{Float64},arg7::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg8::Union(Ptr{PetscSF},StridedArray{PetscSF},Ptr{Void}))
    ccall((:DMPlexCreateTwoSidedProcessSF,petsc1),PetscErrorCode,(DM,PetscSF,PetscSection,IS{Float64},PetscSection,IS{Float64},Ptr{IS{Float64}},Ptr{PetscSF}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexDistributeOwnership(arg1::DM,arg2::PetscSection,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::PetscSection,arg5::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:DMPlexDistributeOwnership,petsc1),PetscErrorCode,(DM,PetscSection,Ptr{IS{Float64}},PetscSection,Ptr{IS{Float64}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateOverlap(arg1::DM,arg2::Integer,arg3::PetscSection,arg4::IS{Float64},arg5::PetscSection,arg6::IS{Float64},arg7::Union(Ptr{DMLabel},StridedArray{DMLabel},Ptr{Void}))
    ccall((:DMPlexCreateOverlap,petsc1),PetscErrorCode,(DM,PetscInt,PetscSection,IS{Float64},PetscSection,IS{Float64},Ptr{DMLabel}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateOverlapMigrationSF(arg0::Type{Float64},arg1::DM,arg2::PetscSF,arg3::Union(Ptr{PetscSF},StridedArray{PetscSF},Ptr{Void}))
    ccall((:DMPlexCreateOverlapMigrationSF,petsc1),PetscErrorCode,(DM,PetscSF,Ptr{PetscSF}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexStratifyMigrationSF(arg0::Type{Float64},arg1::DM,arg2::PetscSF,arg3::Union(Ptr{PetscSF},StridedArray{PetscSF},Ptr{Void}))
    ccall((:DMPlexStratifyMigrationSF,petsc1),PetscErrorCode,(DM,PetscSF,Ptr{PetscSF}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateSubmesh(arg0::Type{Float64},arg1::DM,arg2::DMLabel,arg3::Integer,arg4::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexCreateSubmesh,petsc1),PetscErrorCode,(DM,DMLabel,PetscInt,Ptr{DM}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateHybridMesh(arg0::Type{Float64},arg1::DM,arg2::DMLabel,arg3::Union(Ptr{DMLabel},StridedArray{DMLabel},Ptr{Void}),arg4::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexCreateHybridMesh,petsc1),PetscErrorCode,(DM,DMLabel,Ptr{DMLabel},Ptr{DM}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetSubpointMap(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{DMLabel},StridedArray{DMLabel},Ptr{Void}))
    ccall((:DMPlexGetSubpointMap,petsc1),PetscErrorCode,(DM,Ptr{DMLabel}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetSubpointMap(arg0::Type{Float64},arg1::DM,arg2::DMLabel)
    ccall((:DMPlexSetSubpointMap,petsc1),PetscErrorCode,(DM,DMLabel),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateSubpointIS(arg1::DM,arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:DMPlexCreateSubpointIS,petsc1),PetscErrorCode,(DM,Ptr{IS{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexMarkBoundaryFaces(arg0::Type{Float64},arg1::DM,arg2::DMLabel)
    ccall((:DMPlexMarkBoundaryFaces,petsc1),PetscErrorCode,(DM,DMLabel),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexLabelComplete(arg0::Type{Float64},arg1::DM,arg2::DMLabel)
    ccall((:DMPlexLabelComplete,petsc1),PetscErrorCode,(DM,DMLabel),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexLabelCohesiveComplete(arg0::Type{Float64},arg1::DM,arg2::DMLabel,arg3::DMLabel,arg4::PetscBool,arg5::DM)
    ccall((:DMPlexLabelCohesiveComplete,petsc1),PetscErrorCode,(DM,DMLabel,DMLabel,PetscBool,DM),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexLabelAddCells(arg0::Type{Float64},arg1::DM,arg2::DMLabel)
    ccall((:DMPlexLabelAddCells,petsc1),PetscErrorCode,(DM,DMLabel),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetRefinementLimit(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:DMPlexGetRefinementLimit,petsc1),PetscErrorCode,(DM,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetRefinementLimit(arg0::Type{Float64},arg1::DM,PetscReal::Integer)
    ccall((:DMPlexSetRefinementLimit,petsc1),PetscErrorCode,(DM,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetRefinementUniform(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:DMPlexGetRefinementUniform,petsc1),PetscErrorCode,(DM,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetRefinementUniform(arg0::Type{Float64},arg1::DM,arg2::PetscBool)
    ccall((:DMPlexSetRefinementUniform,petsc1),PetscErrorCode,(DM,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetCoarseDM(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexGetCoarseDM,petsc1),PetscErrorCode,(DM,Ptr{DM}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetCoarseDM(arg0::Type{Float64},arg1::DM,arg2::DM)
    ccall((:DMPlexSetCoarseDM,petsc1),PetscErrorCode,(DM,DM),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateCoarsePointIS(arg1::DM,arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:DMPlexCreateCoarsePointIS,petsc1),PetscErrorCode,(DM,Ptr{IS{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetNumFaceVertices(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexGetNumFaceVertices,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetOrientedFace(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Integer,arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg8::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg9::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:DMPlexGetOrientedFace,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetMinRadius(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:DMPlexGetMinRadius,petsc1),PetscErrorCode,(DM,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetMinRadius(arg0::Type{Float64},arg1::DM,PetscReal::Integer)
    ccall((:DMPlexSetMinRadius,petsc1),PetscErrorCode,(DM,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexComputeCellGeometryFVM(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:DMPlexComputeCellGeometryFVM,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,PetscReal,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexComputeGeometryFVM(arg1::DM,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMPlexComputeGeometryFVM,petsc1),PetscErrorCode,(DM,Ptr{Vec{Float64}},Ptr{Vec{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexComputeGradientFVM(arg1::DM,arg2::PetscFV,arg3::Vec{Float64},arg4::Vec{Float64},arg5::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexComputeGradientFVM,petsc1),PetscErrorCode,(DM,PetscFV,Vec{Float64},Vec{Float64},Ptr{DM}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexInsertBoundaryValues(arg1::DM,arg2::Vec{Float64},PetscReal::Integer,arg3::Vec{Float64},arg4::Vec{Float64},arg5::Vec{Float64})
    ccall((:DMPlexInsertBoundaryValues,petsc1),PetscErrorCode,(DM,Vec{Float64},Cint,Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,PetscReal,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateSection(arg1::DM,arg2::Integer,arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Integer,arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg8::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg9::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg10::IS{Float64},arg11::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}))
    ccall((:DMPlexCreateSection,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{IS{Float64}},Ptr{IS{Float64}},IS{Float64},Ptr{PetscSection}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexComputeCellGeometryAffineFEM(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:DMPlexComputeCellGeometryAffineFEM,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexComputeCellGeometryFEM(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::PetscFE,arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg7::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:DMPlexComputeCellGeometryFEM,petsc1),PetscErrorCode,(DM,PetscInt,PetscFE,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexComputeGeometryFEM(arg1::DM,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMPlexComputeGeometryFEM,petsc1),PetscErrorCode,(DM,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexVecGetClosure(arg1::DM,arg2::PetscSection,arg3::Vec{Float64},arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:DMPlexVecGetClosure,petsc1),PetscErrorCode,(DM,PetscSection,Vec{Float64},PetscInt,Ptr{PetscInt},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexVecRestoreClosure(arg1::DM,arg2::PetscSection,arg3::Vec{Float64},arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:DMPlexVecRestoreClosure,petsc1),PetscErrorCode,(DM,PetscSection,Vec{Float64},PetscInt,Ptr{PetscInt},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexVecSetClosure(arg1::DM,arg2::PetscSection,arg3::Vec{Float64},arg4::Integer,arg5::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg6::InsertMode)
    ccall((:DMPlexVecSetClosure,petsc1),PetscErrorCode,(DM,PetscSection,Vec{Float64},PetscInt,Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexMatSetClosure(arg1::DM,arg2::PetscSection,arg3::PetscSection,arg4::Mat{Float64},arg5::Integer,arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::InsertMode)
    ccall((:DMPlexMatSetClosure,petsc1),PetscErrorCode,(DM,PetscSection,PetscSection,Mat{Float64},PetscInt,Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexMatSetClosureRefined(arg1::DM,arg2::PetscSection,arg3::PetscSection,arg4::DM,arg5::PetscSection,arg6::PetscSection,arg7::Mat{Float64},arg8::Integer,arg9::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg10::InsertMode)
    ccall((:DMPlexMatSetClosureRefined,petsc1),PetscErrorCode,(DM,PetscSection,PetscSection,DM,PetscSection,PetscSection,Mat{Float64},PetscInt,Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexMatGetClosureIndicesRefined(arg0::Type{Float64},arg1::DM,arg2::PetscSection,arg3::PetscSection,arg4::DM,arg5::PetscSection,arg6::PetscSection,arg7::Integer,arg8::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg9::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexMatGetClosureIndicesRefined,petsc1),PetscErrorCode,(DM,PetscSection,PetscSection,DM,PetscSection,PetscSection,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateClosureIndex(arg0::Type{Float64},arg1::DM,arg2::PetscSection)
    ccall((:DMPlexCreateClosureIndex,petsc1),PetscErrorCode,(DM,PetscSection),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateFromFile(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::PetscBool,arg4::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexCreateFromFile,petsc1),PetscErrorCode,(comm_type,Cstring,PetscBool,Ptr{DM}),arg1.val,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateExodus(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::PetscBool,arg4::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexCreateExodus,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscBool,Ptr{DM}),arg1.val,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateExodusFromFile(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::PetscBool,arg4::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexCreateExodusFromFile,petsc1),PetscErrorCode,(comm_type,Cstring,PetscBool,Ptr{DM}),arg1.val,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateCGNS(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::PetscBool,arg4::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexCreateCGNS,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscBool,Ptr{DM}),arg1.val,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateCGNSFromFile(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::PetscBool,arg4::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexCreateCGNSFromFile,petsc1),PetscErrorCode,(comm_type,Cstring,PetscBool,Ptr{DM}),arg1.val,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateGmsh(arg1::MPI_Comm,arg2::PetscViewer{Float64},arg3::PetscBool,arg4::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexCreateGmsh,petsc1),PetscErrorCode,(comm_type,PetscViewer{Float64},PetscBool,Ptr{DM}),arg1.val,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateGmshFromFile(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::PetscBool,arg4::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexCreateGmshFromFile,petsc1),PetscErrorCode,(comm_type,Cstring,PetscBool,Ptr{DM}),arg1.val,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateFluent(arg1::MPI_Comm,arg2::PetscViewer{Float64},arg3::PetscBool,arg4::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexCreateFluent,petsc1),PetscErrorCode,(comm_type,PetscViewer{Float64},PetscBool,Ptr{DM}),arg1.val,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateFluentFromFile(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::PetscBool,arg4::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexCreateFluentFromFile,petsc1),PetscErrorCode,(comm_type,Cstring,PetscBool,Ptr{DM}),arg1.val,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexConstructGhostCells(arg0::Type{Float64},arg1::DM,arg2::Union(ByteString,Symbol),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexConstructGhostCells,petsc1),PetscErrorCode,(DM,Cstring,Ptr{PetscInt},Ptr{DM}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexConstructCohesiveCells(arg0::Type{Float64},arg1::DM,arg2::DMLabel,arg3::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexConstructCohesiveCells,petsc1),PetscErrorCode,(DM,DMLabel,Ptr{DM}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetHybridBounds(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexGetHybridBounds,petsc1),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetHybridBounds(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer)
    ccall((:DMPlexSetHybridBounds,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetVTKCellHeight(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexGetVTKCellHeight,petsc1),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetVTKCellHeight(arg0::Type{Float64},arg1::DM,arg2::Integer)
    ccall((:DMPlexSetVTKCellHeight,petsc1),PetscErrorCode,(DM,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexVTKWriteAll(arg1::PetscObject,arg2::PetscViewer{Float64})
    ccall((:DMPlexVTKWriteAll,petsc1),PetscErrorCode,(PetscObject,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetScale(arg0::Type{Float64},arg1::DM,arg2::PetscUnit,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:DMPlexGetScale,petsc1),PetscErrorCode,(DM,PetscUnit,Ptr{Cint}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetScale(arg0::Type{Float64},arg1::DM,arg2::PetscUnit,PetscReal::Integer)
    ccall((:DMPlexSetScale,petsc1),PetscErrorCode,(DM,PetscUnit,Cint),arg1,arg2,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexAddBoundary(arg0::Type{Float64},arg1::DM,arg2::PetscBool,arg3::Union(ByteString,Symbol),arg4::Union(ByteString,Symbol),arg5::Integer,arg6::Integer,arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg8::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg9::Integer,arg10::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg11::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMPlexAddBoundary,petsc1),PetscErrorCode,(DM,PetscBool,Cstring,Cstring,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Void},PetscInt,Ptr{PetscInt},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetNumBoundary(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexGetNumBoundary,petsc1),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetBoundary(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg4::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}),arg5::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg8::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg9::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg10::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg11::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg12::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:DMPlexGetBoundary,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscBool},Ptr{Ptr{Uint8}},Ptr{Ptr{Uint8}},Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{Void}},Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexIsBoundaryPoint(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:DMPlexIsBoundaryPoint,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscBool}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCopyBoundary(arg0::Type{Float64},arg1::DM,arg2::DM)
    ccall((:DMPlexCopyBoundary,petsc1),PetscErrorCode,(DM,DM),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexInsertBoundaryValuesFEM(arg1::DM,arg2::Vec{Float64})
    ccall((:DMPlexInsertBoundaryValuesFEM,petsc1),PetscErrorCode,(DM,Vec{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetMaxProjectionHeight(arg0::Type{Float64},arg1::DM,arg2::Integer)
    ccall((:DMPlexSetMaxProjectionHeight,petsc1),PetscErrorCode,(DM,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetMaxProjectionHeight(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexGetMaxProjectionHeight,petsc1),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexProjectFunction(arg1::DM,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg4::InsertMode,arg5::Vec{Float64})
    ccall((:DMPlexProjectFunction,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}},InsertMode,Vec{Float64}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexProjectFunctionLocal(arg1::DM,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg4::InsertMode,arg5::Vec{Float64})
    ccall((:DMPlexProjectFunctionLocal,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}},InsertMode,Vec{Float64}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexProjectFieldLocal(arg1::DM,arg2::Vec{Float64},arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg4::InsertMode,arg5::Vec{Float64})
    ccall((:DMPlexProjectFieldLocal,petsc1),PetscErrorCode,(DM,Vec{Float64},Ptr{Ptr{Void}},InsertMode,Vec{Float64}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexComputeL2Diff(arg1::DM,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg4::Vec{Float64},arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:DMPlexComputeL2Diff,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}},Vec{Float64},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexComputeL2GradientDiff(arg1::DM,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg4::Vec{Float64},PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:DMPlexComputeL2GradientDiff,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}},Vec{Float64},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,PetscReal,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexComputeL2FieldDiff(arg1::DM,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg4::Vec{Float64},PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:DMPlexComputeL2FieldDiff,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}},Vec{Float64},Ptr{Cint}),arg1,arg2,arg3,arg4,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexComputeIntegralFEM(arg1::DM,arg2::Vec{Float64},arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMPlexComputeIntegralFEM,petsc1),PetscErrorCode,(DM,Vec{Float64},Ptr{Cint},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexComputeInterpolatorFEM(arg1::DM,arg2::DM,arg3::Mat{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMPlexComputeInterpolatorFEM,petsc1),PetscErrorCode,(DM,DM,Mat{Float64},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexComputeInjectorFEM(arg1::DM,arg2::DM,arg3::Union(Ptr{VecScatter{Float64}},StridedArray{VecScatter{Float64}},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMPlexComputeInjectorFEM,petsc1),PetscErrorCode,(DM,DM,Ptr{VecScatter{Float64}},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateRigidBody(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{MatNullSpace},StridedArray{MatNullSpace},Ptr{Void}))
    ccall((:DMPlexCreateRigidBody,petsc1),PetscErrorCode,(DM,Ptr{MatNullSpace}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSNESComputeResidualFEM(arg1::DM,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMPlexSNESComputeResidualFEM,petsc1),PetscErrorCode,(DM,Vec{Float64},Vec{Float64},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSNESComputeJacobianFEM(arg1::DM,arg2::Vec{Float64},arg3::Mat{Float64},arg4::Mat{Float64},arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMPlexSNESComputeJacobianFEM,petsc1),PetscErrorCode,(DM,Vec{Float64},Mat{Float64},Mat{Float64},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexTSComputeRHSFunctionFVM(arg1::DM,PetscReal::Integer,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMPlexTSComputeRHSFunctionFVM,petsc1),PetscErrorCode,(DM,Cint,Vec{Float64},Vec{Float64},Ptr{Void}),arg1,PetscReal,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexTSComputeIFunctionFEM(arg1::DM,PetscReal::Integer,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64},arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMPlexTSComputeIFunctionFEM,petsc1),PetscErrorCode,(DM,Cint,Vec{Float64},Vec{Float64},Vec{Float64},Ptr{Void}),arg1,PetscReal,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexComputeRHSFunctionFVM(arg1::DM,PetscReal::Integer,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMPlexComputeRHSFunctionFVM,petsc1),PetscErrorCode,(DM,Cint,Vec{Float64},Vec{Float64},Ptr{Void}),arg1,PetscReal,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetAnchors(arg1::DM,arg2::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}),arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:DMPlexGetAnchors,petsc1),PetscErrorCode,(DM,Ptr{PetscSection},Ptr{IS{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetAnchors(arg1::DM,arg2::PetscSection,arg3::IS{Float64})
    ccall((:DMPlexSetAnchors,petsc1),PetscErrorCode,(DM,PetscSection,IS{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetReferenceTree(arg0::Type{Float64},arg1::DM,arg2::DM)
    ccall((:DMPlexSetReferenceTree,petsc1),PetscErrorCode,(DM,DM),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetReferenceTree(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexGetReferenceTree,petsc1),PetscErrorCode,(DM,Ptr{DM}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexReferenceTreeGetChildSymmetry(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg8::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexReferenceTreeGetChildSymmetry,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexCreateDefaultReferenceTree(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::PetscBool,arg4::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexCreateDefaultReferenceTree,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscBool,Ptr{DM}),arg1.val,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSetTree(arg0::Type{Float64},arg1::DM,arg2::PetscSection,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexSetTree,petsc1),PetscErrorCode,(DM,PetscSection,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetTree(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}),arg3::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg4::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg5::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}),arg6::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMPlexGetTree,petsc1),PetscErrorCode,(DM,Ptr{PetscSection},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{PetscSection},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetTreeParent(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMPlexGetTreeParent,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetTreeChildren(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:DMPlexGetTreeChildren,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexTreeRefineCell(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexTreeRefineCell,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{DM}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMRedundantCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::PetscMPIInt,arg3::Integer,arg4::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMRedundantCreate,petsc1),PetscErrorCode,(comm_type,PetscMPIInt,PetscInt,Ptr{DM}),arg1.val,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMRedundantSetSize(arg0::Type{Float64},arg1::DM,arg2::PetscMPIInt,arg3::Integer)
    ccall((:DMRedundantSetSize,petsc1),PetscErrorCode,(DM,PetscMPIInt,PetscInt),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMRedundantGetSize(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMRedundantGetSize,petsc1),PetscErrorCode,(DM,Ptr{PetscMPIInt},Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMShellCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMShellCreate,petsc1),PetscErrorCode,(comm_type,Ptr{DM}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMShellSetMatrix(arg1::DM,arg2::Mat{Float64})
    ccall((:DMShellSetMatrix,petsc1),PetscErrorCode,(DM,Mat{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMShellSetGlobalVector(arg1::DM,arg2::Vec{Float64})
    ccall((:DMShellSetGlobalVector,petsc1),PetscErrorCode,(DM,Vec{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMShellSetLocalVector(arg1::DM,arg2::Vec{Float64})
    ccall((:DMShellSetLocalVector,petsc1),PetscErrorCode,(DM,Vec{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMShellSetCreateGlobalVector(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMShellSetCreateGlobalVector,petsc1),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMShellSetCreateLocalVector(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMShellSetCreateLocalVector,petsc1),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMShellSetGlobalToLocal(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMShellSetGlobalToLocal,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMShellSetGlobalToLocalVecScatter(arg1::DM,arg2::VecScatter{Float64})
    ccall((:DMShellSetGlobalToLocalVecScatter,petsc1),PetscErrorCode,(DM,VecScatter{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMShellSetLocalToGlobal(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMShellSetLocalToGlobal,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMShellSetLocalToGlobalVecScatter(arg1::DM,arg2::VecScatter{Float64})
    ccall((:DMShellSetLocalToGlobalVecScatter,petsc1),PetscErrorCode,(DM,VecScatter{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMShellSetLocalToLocal(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMShellSetLocalToLocal,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMShellSetLocalToLocalVecScatter(arg1::DM,arg2::VecScatter{Float64})
    ccall((:DMShellSetLocalToLocalVecScatter,petsc1),PetscErrorCode,(DM,VecScatter{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMShellSetCreateMatrix(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMShellSetCreateMatrix,petsc1),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMShellSetCoarsen(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMShellSetCoarsen,petsc1),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMShellSetRefine(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMShellSetRefine,petsc1),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMShellSetCreateInterpolation(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMShellSetCreateInterpolation,petsc1),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMShellSetCreateInjection(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMShellSetCreateInjection,petsc1),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMShellSetCreateFieldDecomposition(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMShellSetCreateFieldDecomposition,petsc1),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMGlobalToLocalBeginDefaultShell(arg1::DM,arg2::Vec{Float64},arg3::InsertMode,arg4::Vec{Float64})
    ccall((:DMGlobalToLocalBeginDefaultShell,petsc1),PetscErrorCode,(DM,Vec{Float64},InsertMode,Vec{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMGlobalToLocalEndDefaultShell(arg1::DM,arg2::Vec{Float64},arg3::InsertMode,arg4::Vec{Float64})
    ccall((:DMGlobalToLocalEndDefaultShell,petsc1),PetscErrorCode,(DM,Vec{Float64},InsertMode,Vec{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMLocalToGlobalBeginDefaultShell(arg1::DM,arg2::Vec{Float64},arg3::InsertMode,arg4::Vec{Float64})
    ccall((:DMLocalToGlobalBeginDefaultShell,petsc1),PetscErrorCode,(DM,Vec{Float64},InsertMode,Vec{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMLocalToGlobalEndDefaultShell(arg1::DM,arg2::Vec{Float64},arg3::InsertMode,arg4::Vec{Float64})
    ccall((:DMLocalToGlobalEndDefaultShell,petsc1),PetscErrorCode,(DM,Vec{Float64},InsertMode,Vec{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMLocalToLocalBeginDefaultShell(arg1::DM,arg2::Vec{Float64},arg3::InsertMode,arg4::Vec{Float64})
    ccall((:DMLocalToLocalBeginDefaultShell,petsc1),PetscErrorCode,(DM,Vec{Float64},InsertMode,Vec{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMLocalToLocalEndDefaultShell(arg1::DM,arg2::Vec{Float64},arg3::InsertMode,arg4::Vec{Float64})
    ccall((:DMLocalToLocalEndDefaultShell,petsc1),PetscErrorCode,(DM,Vec{Float64},InsertMode,Vec{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMSlicedCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg8::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMSlicedCreate,petsc1),PetscErrorCode,(comm_type,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{DM}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end 
=#
#= skipping function with undefined symbols: 
 function DMSlicedSetPreallocation(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMSlicedSetPreallocation,petsc1),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMSlicedSetBlockFills(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMSlicedSetBlockFills,petsc1),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMSlicedSetGhosts(arg0::Type{Float64},arg1::DM,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:DMSlicedSetGhosts,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end 
=#
function PetscDSInitializePackage(arg0::Type{Float64})
    err = ccall((:PetscDSInitializePackage,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function PetscDSCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscDS},StridedArray{PetscDS},Ptr{Void}))
    ccall((:PetscDSCreate,petsc1),PetscErrorCode,(comm_type,Ptr{PetscDS}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscDS},StridedArray{PetscDS},Ptr{Void}))
    ccall((:PetscDSDestroy,petsc1),PetscErrorCode,(Ptr{PetscDS},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSSetType(arg0::Type{Float64},arg1::PetscDS,arg2::PetscDSType)
    ccall((:PetscDSSetType,petsc1),PetscErrorCode,(PetscDS,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetType(arg0::Type{Float64},arg1::PetscDS,arg2::Union(Ptr{PetscDSType},StridedArray{PetscDSType},Ptr{Void}))
    ccall((:PetscDSGetType,petsc1),PetscErrorCode,(PetscDS,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSSetUp(arg0::Type{Float64},arg1::PetscDS)
    ccall((:PetscDSSetUp,petsc1),PetscErrorCode,(PetscDS,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSSetFromOptions(arg0::Type{Float64},arg1::PetscDS)
    ccall((:PetscDSSetFromOptions,petsc1),PetscErrorCode,(PetscDS,),arg1)
end 
=#
function PetscDSRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:PetscDSRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

function PetscDSRegisterDestroy(arg0::Type{Float64})
    err = ccall((:PetscDSRegisterDestroy,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function PetscDSGetSpatialDimension(arg0::Type{Float64},arg1::PetscDS,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscDSGetSpatialDimension,petsc1),PetscErrorCode,(PetscDS,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetNumFields(arg0::Type{Float64},arg1::PetscDS,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscDSGetNumFields,petsc1),PetscErrorCode,(PetscDS,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetTotalDimension(arg0::Type{Float64},arg1::PetscDS,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscDSGetTotalDimension,petsc1),PetscErrorCode,(PetscDS,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetTotalBdDimension(arg0::Type{Float64},arg1::PetscDS,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscDSGetTotalBdDimension,petsc1),PetscErrorCode,(PetscDS,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetTotalComponents(arg0::Type{Float64},arg1::PetscDS,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscDSGetTotalComponents,petsc1),PetscErrorCode,(PetscDS,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetFieldOffset(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscDSGetFieldOffset,petsc1),PetscErrorCode,(PetscDS,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetBdFieldOffset(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscDSGetBdFieldOffset,petsc1),PetscErrorCode,(PetscDS,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetComponentOffset(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:PetscDSGetComponentOffset,petsc1),PetscErrorCode,(PetscDS,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetComponentOffsets(arg0::Type{Float64},arg1::PetscDS,arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:PetscDSGetComponentOffsets,petsc1),PetscErrorCode,(PetscDS,Ptr{Ptr{PetscInt}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetComponentBdOffsets(arg0::Type{Float64},arg1::PetscDS,arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:PetscDSGetComponentBdOffsets,petsc1),PetscErrorCode,(PetscDS,Ptr{Ptr{PetscInt}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetComponentDerivativeOffsets(arg0::Type{Float64},arg1::PetscDS,arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:PetscDSGetComponentDerivativeOffsets,petsc1),PetscErrorCode,(PetscDS,Ptr{Ptr{PetscInt}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetComponentBdDerivativeOffsets(arg0::Type{Float64},arg1::PetscDS,arg2::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}))
    ccall((:PetscDSGetComponentBdDerivativeOffsets,petsc1),PetscErrorCode,(PetscDS,Ptr{Ptr{PetscInt}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetDiscretization(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::Union(Ptr{PetscObject},StridedArray{PetscObject},Ptr{Void}))
    ccall((:PetscDSGetDiscretization,petsc1),PetscErrorCode,(PetscDS,PetscInt,Ptr{PetscObject}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSSetDiscretization(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::PetscObject)
    ccall((:PetscDSSetDiscretization,petsc1),PetscErrorCode,(PetscDS,PetscInt,PetscObject),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSAddDiscretization(arg0::Type{Float64},arg1::PetscDS,arg2::PetscObject)
    ccall((:PetscDSAddDiscretization,petsc1),PetscErrorCode,(PetscDS,PetscObject),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetBdDiscretization(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::Union(Ptr{PetscObject},StridedArray{PetscObject},Ptr{Void}))
    ccall((:PetscDSGetBdDiscretization,petsc1),PetscErrorCode,(PetscDS,PetscInt,Ptr{PetscObject}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSSetBdDiscretization(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::PetscObject)
    ccall((:PetscDSSetBdDiscretization,petsc1),PetscErrorCode,(PetscDS,PetscInt,PetscObject),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSAddBdDiscretization(arg0::Type{Float64},arg1::PetscDS,arg2::PetscObject)
    ccall((:PetscDSAddBdDiscretization,petsc1),PetscErrorCode,(PetscDS,PetscObject),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetImplicit(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscDSGetImplicit,petsc1),PetscErrorCode,(PetscDS,PetscInt,Ptr{PetscBool}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSSetImplicit(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::PetscBool)
    ccall((:PetscDSSetImplicit,petsc1),PetscErrorCode,(PetscDS,PetscInt,PetscBool),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetAdjacency(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscDSGetAdjacency,petsc1),PetscErrorCode,(PetscDS,PetscInt,Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSSetAdjacency(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::PetscBool,arg4::PetscBool)
    ccall((:PetscDSSetAdjacency,petsc1),PetscErrorCode,(PetscDS,PetscInt,PetscBool,PetscBool),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetObjective(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:PetscDSGetObjective,petsc1),PetscErrorCode,(PetscDS,PetscInt,Ptr{Ptr{Void}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSSetObjective(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscDSSetObjective,petsc1),PetscErrorCode,(PetscDS,PetscInt,Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetResidual(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg4::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:PetscDSGetResidual,petsc1),PetscErrorCode,(PetscDS,PetscInt,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSSetResidual(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscDSSetResidual,petsc1),PetscErrorCode,(PetscDS,PetscInt,Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetJacobian(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg5::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg6::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg7::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:PetscDSGetJacobian,petsc1),PetscErrorCode,(PetscDS,PetscInt,PetscInt,Ptr{Ptr{Void}},Ptr{Ptr{Void}},Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSSetJacobian(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg6::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg7::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscDSSetJacobian,petsc1),PetscErrorCode,(PetscDS,PetscInt,PetscInt,Ptr{Void},Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetRiemannSolver(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:PetscDSGetRiemannSolver,petsc1),PetscErrorCode,(PetscDS,PetscInt,Ptr{Ptr{Void}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSSetRiemannSolver(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscDSSetRiemannSolver,petsc1),PetscErrorCode,(PetscDS,PetscInt,Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetContext(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:PetscDSGetContext,petsc1),PetscErrorCode,(PetscDS,PetscInt,Ptr{Ptr{Void}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSSetContext(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscDSSetContext,petsc1),PetscErrorCode,(PetscDS,PetscInt,Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetBdResidual(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg4::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:PetscDSGetBdResidual,petsc1),PetscErrorCode,(PetscDS,PetscInt,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSSetBdResidual(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscDSSetBdResidual,petsc1),PetscErrorCode,(PetscDS,PetscInt,Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetBdJacobian(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg5::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg6::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg7::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:PetscDSGetBdJacobian,petsc1),PetscErrorCode,(PetscDS,PetscInt,PetscInt,Ptr{Ptr{Void}},Ptr{Ptr{Void}},Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSSetBdJacobian(arg0::Type{Float64},arg1::PetscDS,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg6::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg7::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscDSSetBdJacobian,petsc1),PetscErrorCode,(PetscDS,PetscInt,PetscInt,Ptr{Void},Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetTabulation(arg0::Type{Float64},arg1::PetscDS,arg2::Union(Ptr{Ptr{Ptr{Cint}}},StridedArray{Ptr{Ptr{Cint}}},Ptr{Void}),arg3::Union(Ptr{Ptr{Ptr{Cint}}},StridedArray{Ptr{Ptr{Cint}}},Ptr{Void}))
    ccall((:PetscDSGetTabulation,petsc1),PetscErrorCode,(PetscDS,Ptr{Ptr{Ptr{Cint}}},Ptr{Ptr{Ptr{Cint}}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetBdTabulation(arg0::Type{Float64},arg1::PetscDS,arg2::Union(Ptr{Ptr{Ptr{Cint}}},StridedArray{Ptr{Ptr{Cint}}},Ptr{Void}),arg3::Union(Ptr{Ptr{Ptr{Cint}}},StridedArray{Ptr{Ptr{Cint}}},Ptr{Void}))
    ccall((:PetscDSGetBdTabulation,petsc1),PetscErrorCode,(PetscDS,Ptr{Ptr{Ptr{Cint}}},Ptr{Ptr{Ptr{Cint}}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetEvaluationArrays(arg0::Type{Float64},arg1::PetscDS,arg2::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}),arg3::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}),arg4::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:PetscDSGetEvaluationArrays,petsc1),PetscErrorCode,(PetscDS,Ptr{Ptr{Float64}},Ptr{Ptr{Float64}},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetWeakFormArrays(arg0::Type{Float64},arg1::PetscDS,arg2::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}),arg3::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}),arg4::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}),arg5::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}),arg6::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}),arg7::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:PetscDSGetWeakFormArrays,petsc1),PetscErrorCode,(PetscDS,Ptr{Ptr{Float64}},Ptr{Ptr{Float64}},Ptr{Ptr{Float64}},Ptr{Ptr{Float64}},Ptr{Ptr{Float64}},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDSGetRefCoordArrays(arg0::Type{Float64},arg1::PetscDS,arg2::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg3::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:PetscDSGetRefCoordArrays,petsc1),PetscErrorCode,(PetscDS,Ptr{Ptr{Cint}},Ptr{Ptr{Float64}}),arg1,arg2,arg3)
end 
=#
function CharacteristicInitializePackage(arg0::Type{Float64})
    err = ccall((:CharacteristicInitializePackage,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function CharacteristicCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{Characteristic},StridedArray{Characteristic},Ptr{Void}))
    ccall((:CharacteristicCreate,petsc1),PetscErrorCode,(comm_type,Ptr{Characteristic}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function CharacteristicSetType(arg0::Type{Float64},arg1::Characteristic,arg2::CharacteristicType)
    ccall((:CharacteristicSetType,petsc1),PetscErrorCode,(Characteristic,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function CharacteristicSetUp(arg0::Type{Float64},arg1::Characteristic)
    ccall((:CharacteristicSetUp,petsc1),PetscErrorCode,(Characteristic,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function CharacteristicSetVelocityInterpolation(arg1::Characteristic,arg2::DM,arg3::Vec{Float64},arg4::Vec{Float64},arg5::Integer,arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg8::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:CharacteristicSetVelocityInterpolation,petsc1),PetscErrorCode,(Characteristic,DM,Vec{Float64},Vec{Float64},PetscInt,Ptr{PetscInt},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end 
=#
#= skipping function with undefined symbols: 
 function CharacteristicSetVelocityInterpolationLocal(arg1::Characteristic,arg2::DM,arg3::Vec{Float64},arg4::Vec{Float64},arg5::Integer,arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg7::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg8::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:CharacteristicSetVelocityInterpolationLocal,petsc1),PetscErrorCode,(Characteristic,DM,Vec{Float64},Vec{Float64},PetscInt,Ptr{PetscInt},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end 
=#
#= skipping function with undefined symbols: 
 function CharacteristicSetFieldInterpolation(arg1::Characteristic,arg2::DM,arg3::Vec{Float64},arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg7::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:CharacteristicSetFieldInterpolation,petsc1),PetscErrorCode,(Characteristic,DM,Vec{Float64},PetscInt,Ptr{PetscInt},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function CharacteristicSetFieldInterpolationLocal(arg1::Characteristic,arg2::DM,arg3::Vec{Float64},arg4::Integer,arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg7::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:CharacteristicSetFieldInterpolationLocal,petsc1),PetscErrorCode,(Characteristic,DM,Vec{Float64},PetscInt,Ptr{PetscInt},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function CharacteristicSolve(arg1::Characteristic,PetscReal::Integer,arg2::Vec{Float64})
    ccall((:CharacteristicSolve,petsc1),PetscErrorCode,(Characteristic,Cint,Vec{Float64}),arg1,PetscReal,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function CharacteristicDestroy(arg0::Type{Float64},arg1::Union(Ptr{Characteristic},StridedArray{Characteristic},Ptr{Void}))
    ccall((:CharacteristicDestroy,petsc1),PetscErrorCode,(Ptr{Characteristic},),arg1)
end 
=#
function CharacteristicRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:CharacteristicRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
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
function SNESInitializePackage(arg0::Type{Float64})
    err = ccall((:SNESInitializePackage,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function SNESCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{SNES},StridedArray{SNES},Ptr{Void}))
    ccall((:SNESCreate,petsc1),PetscErrorCode,(comm_type,Ptr{SNES}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESReset(arg0::Type{Float64},arg1::SNES)
    ccall((:SNESReset,petsc1),PetscErrorCode,(SNES,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function SNESDestroy(arg0::Type{Float64},arg1::Union(Ptr{SNES},StridedArray{SNES},Ptr{Void}))
    ccall((:SNESDestroy,petsc1),PetscErrorCode,(Ptr{SNES},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetType(arg0::Type{Float64},arg1::SNES,arg2::SNESType)
    ccall((:SNESSetType,petsc1),PetscErrorCode,(SNES,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMonitor(arg0::Type{Float64},arg1::SNES,arg2::Integer,PetscReal::Integer)
    ccall((:SNESMonitor,petsc1),PetscErrorCode,(SNES,PetscInt,Cint),arg1,arg2,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMonitorSet(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESMonitorSet,petsc1),PetscErrorCode,(SNES,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMonitorCancel(arg0::Type{Float64},arg1::SNES)
    ccall((:SNESMonitorCancel,petsc1),PetscErrorCode,(SNES,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMonitorSAWs(arg0::Type{Float64},arg1::SNES,arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESMonitorSAWs,petsc1),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMonitorSAWsCreate(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:SNESMonitorSAWsCreate,petsc1),PetscErrorCode,(SNES,Ptr{Ptr{Void}}),arg1,arg2)
end 
=#
function SNESMonitorSAWsDestroy(arg0::Type{Float64},arg1::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    err = ccall((:SNESMonitorSAWsDestroy,petsc1),PetscErrorCode,(Ptr{Ptr{Void}},),arg1)
    return err
end

#= skipping function with undefined symbols: 
 function SNESSetConvergenceHistory(arg0::Type{Float64},arg1::SNES,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Integer,arg4::PetscBool)
    ccall((:SNESSetConvergenceHistory,petsc1),PetscErrorCode,(SNES,Ptr{Cint},Ptr{PetscInt},PetscInt,PetscBool),arg1,PetscReal,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetConvergenceHistory(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg3::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:SNESGetConvergenceHistory,petsc1),PetscErrorCode,(SNES,Ptr{Ptr{Cint}},Ptr{Ptr{PetscInt}},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetUp(arg0::Type{Float64},arg1::SNES)
    ccall((:SNESSetUp,petsc1),PetscErrorCode,(SNES,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSolve(arg1::SNES,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:SNESSolve,petsc1),PetscErrorCode,(SNES,Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetErrorIfNotConverged(arg0::Type{Float64},arg1::SNES,arg2::PetscBool)
    ccall((:SNESSetErrorIfNotConverged,petsc1),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetErrorIfNotConverged(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:SNESGetErrorIfNotConverged,petsc1),PetscErrorCode,(SNES,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetWorkVecs(arg0::Type{Float64},arg1::SNES,arg2::Integer)
    ccall((:SNESSetWorkVecs,petsc1),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end 
=#
function SNESAddOptionsChecker(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:SNESAddOptionsChecker,petsc1),PetscErrorCode,(Ptr{Void},),arg1)
    return err
end

#= skipping function with undefined symbols: 
 function SNESSetUpdate(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESSetUpdate,petsc1),PetscErrorCode,(SNES,Ptr{Void}),arg1,arg2)
end 
=#
function SNESRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:SNESRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function SNESGetKSP(arg1::SNES,arg2::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    ccall((:SNESGetKSP,petsc1),PetscErrorCode,(SNES,Ptr{KSP{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetKSP(arg1::SNES,arg2::KSP{Float64})
    ccall((:SNESSetKSP,petsc1),PetscErrorCode,(SNES,KSP{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetSolution(arg1::SNES,arg2::Vec{Float64})
    ccall((:SNESSetSolution,petsc1),PetscErrorCode,(SNES,Vec{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetSolution(arg1::SNES,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:SNESGetSolution,petsc1),PetscErrorCode,(SNES,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetSolutionUpdate(arg1::SNES,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:SNESGetSolutionUpdate,petsc1),PetscErrorCode,(SNES,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetRhs(arg1::SNES,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:SNESGetRhs,petsc1),PetscErrorCode,(SNES,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESView(arg1::SNES,arg2::PetscViewer{Float64})
    ccall((:SNESView,petsc1),PetscErrorCode,(SNES,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLoad(arg1::SNES,arg2::PetscViewer{Float64})
    ccall((:SNESLoad,petsc1),PetscErrorCode,(SNES,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESReasonViewFromOptions(arg0::Type{Float64},arg1::SNES)
    ccall((:SNESReasonViewFromOptions,petsc1),PetscErrorCode,(SNES,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetOptionsPrefix(arg0::Type{Float64},arg1::SNES,arg2::Union(ByteString,Symbol))
    ccall((:SNESSetOptionsPrefix,petsc1),PetscErrorCode,(SNES,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESAppendOptionsPrefix(arg0::Type{Float64},arg1::SNES,arg2::Union(ByteString,Symbol))
    ccall((:SNESAppendOptionsPrefix,petsc1),PetscErrorCode,(SNES,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetOptionsPrefix(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:SNESGetOptionsPrefix,petsc1),PetscErrorCode,(SNES,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetFromOptions(arg0::Type{Float64},arg1::SNES)
    ccall((:SNESSetFromOptions,petsc1),PetscErrorCode,(SNES,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function MatCreateSNESMF(arg1::SNES,arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateSNESMF,petsc1),PetscErrorCode,(SNES,Ptr{Mat{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatMFFDComputeJacobian(arg1::SNES,arg2::Vec{Float64},arg3::Mat{Float64},arg4::Mat{Float64},arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:MatMFFDComputeJacobian,petsc1),PetscErrorCode,(SNES,Vec{Float64},Mat{Float64},Mat{Float64},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function MatDAADSetSNES(arg1::Mat{Float64},arg2::SNES)
    ccall((:MatDAADSetSNES,petsc1),PetscErrorCode,(Mat{Float64},SNES),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetType(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{SNESType},StridedArray{SNESType},Ptr{Void}))
    ccall((:SNESGetType,petsc1),PetscErrorCode,(SNES,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMonitorDefault(arg0::Type{Float64},arg1::SNES,arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESMonitorDefault,petsc1),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMonitorRange(arg0::Type{Float64},arg1::SNES,arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESMonitorRange,petsc1),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMonitorRatio(arg0::Type{Float64},arg1::SNES,arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESMonitorRatio,petsc1),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMonitorSetRatio(arg1::SNES,arg2::PetscViewer{Float64})
    ccall((:SNESMonitorSetRatio,petsc1),PetscErrorCode,(SNES,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMonitorSolution(arg0::Type{Float64},arg1::SNES,arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESMonitorSolution,petsc1),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMonitorResidual(arg0::Type{Float64},arg1::SNES,arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESMonitorResidual,petsc1),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMonitorSolutionUpdate(arg0::Type{Float64},arg1::SNES,arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESMonitorSolutionUpdate,petsc1),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMonitorDefaultShort(arg0::Type{Float64},arg1::SNES,arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESMonitorDefaultShort,petsc1),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMonitorDefaultField(arg0::Type{Float64},arg1::SNES,arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESMonitorDefaultField,petsc1),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMonitorJacUpdateSpectrum(arg0::Type{Float64},arg1::SNES,arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESMonitorJacUpdateSpectrum,petsc1),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMonitorFields(arg0::Type{Float64},arg1::SNES,arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESMonitorFields,petsc1),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end 
=#
function KSPMonitorSNES(arg1::KSP{Float64},arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:KSPMonitorSNES,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
    return err
end

#= skipping function with undefined symbols: 
 function KSPMonitorSNESLGResidualNormCreate(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Union(Ptr{Ptr{PetscObject}},StridedArray{Ptr{PetscObject}},Ptr{Void}))
    ccall((:KSPMonitorSNESLGResidualNormCreate,petsc1),PetscErrorCode,(Cstring,Cstring,Cint,Cint,Cint,Cint,Ptr{Ptr{PetscObject}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function KSPMonitorSNESLGResidualNorm(arg1::KSP{Float64},arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{PetscObject},StridedArray{PetscObject},Ptr{Void}))
    ccall((:KSPMonitorSNESLGResidualNorm,petsc1),PetscErrorCode,(KSP{Float64},PetscInt,Cint,Ptr{PetscObject}),arg1,arg2,PetscReal,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function KSPMonitorSNESLGResidualNormDestroy(arg0::Type{Float64},arg1::Union(Ptr{Ptr{PetscObject}},StridedArray{Ptr{PetscObject}},Ptr{Void}))
    ccall((:KSPMonitorSNESLGResidualNormDestroy,petsc1),PetscErrorCode,(Ptr{Ptr{PetscObject}},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetTolerances(arg0::Type{Float64},arg1::SNES,PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer)
    ccall((:SNESSetTolerances,petsc1),PetscErrorCode,(SNES,Cint,Cint,Cint,PetscInt,PetscInt),arg1,PetscReal,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetTolerances(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:SNESGetTolerances,petsc1),PetscErrorCode,(SNES,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetTrustRegionTolerance(arg0::Type{Float64},arg1::SNES,PetscReal::Integer)
    ccall((:SNESSetTrustRegionTolerance,petsc1),PetscErrorCode,(SNES,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetIterationNumber(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:SNESGetIterationNumber,petsc1),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetIterationNumber(arg0::Type{Float64},arg1::SNES,arg2::Integer)
    ccall((:SNESSetIterationNumber,petsc1),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetNonlinearStepFailures(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:SNESGetNonlinearStepFailures,petsc1),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetMaxNonlinearStepFailures(arg0::Type{Float64},arg1::SNES,arg2::Integer)
    ccall((:SNESSetMaxNonlinearStepFailures,petsc1),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetMaxNonlinearStepFailures(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:SNESGetMaxNonlinearStepFailures,petsc1),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetNumberFunctionEvals(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:SNESGetNumberFunctionEvals,petsc1),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetLagPreconditioner(arg0::Type{Float64},arg1::SNES,arg2::Integer)
    ccall((:SNESSetLagPreconditioner,petsc1),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetLagPreconditioner(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:SNESGetLagPreconditioner,petsc1),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetLagJacobian(arg0::Type{Float64},arg1::SNES,arg2::Integer)
    ccall((:SNESSetLagJacobian,petsc1),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetLagJacobian(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:SNESGetLagJacobian,petsc1),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetLagPreconditionerPersists(arg0::Type{Float64},arg1::SNES,arg2::PetscBool)
    ccall((:SNESSetLagPreconditionerPersists,petsc1),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetLagJacobianPersists(arg0::Type{Float64},arg1::SNES,arg2::PetscBool)
    ccall((:SNESSetLagJacobianPersists,petsc1),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetGridSequence(arg0::Type{Float64},arg1::SNES,arg2::Integer)
    ccall((:SNESSetGridSequence,petsc1),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetGridSequence(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:SNESGetGridSequence,petsc1),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetLinearSolveIterations(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:SNESGetLinearSolveIterations,petsc1),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetLinearSolveFailures(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:SNESGetLinearSolveFailures,petsc1),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetMaxLinearSolveFailures(arg0::Type{Float64},arg1::SNES,arg2::Integer)
    ccall((:SNESSetMaxLinearSolveFailures,petsc1),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetMaxLinearSolveFailures(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:SNESGetMaxLinearSolveFailures,petsc1),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetCountersReset(arg0::Type{Float64},arg1::SNES,arg2::PetscBool)
    ccall((:SNESSetCountersReset,petsc1),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESKSPSetUseEW(arg0::Type{Float64},arg1::SNES,arg2::PetscBool)
    ccall((:SNESKSPSetUseEW,petsc1),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESKSPGetUseEW(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:SNESKSPGetUseEW,petsc1),PetscErrorCode,(SNES,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESKSPSetParametersEW(arg0::Type{Float64},arg1::SNES,arg2::Integer,PetscReal::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer)
    ccall((:SNESKSPSetParametersEW,petsc1),PetscErrorCode,(SNES,PetscInt,Cint,Cint,Cint,Cint,Cint,Cint),arg1,arg2,PetscReal,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function SNESKSPGetParametersEW(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg7::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg8::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:SNESKSPGetParametersEW,petsc1),PetscErrorCode,(SNES,Ptr{PetscInt},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMonitorLGCreate(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(ByteString,Symbol),arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Union(Ptr{Ptr{PetscObject}},StridedArray{Ptr{PetscObject}},Ptr{Void}))
    ccall((:SNESMonitorLGCreate,petsc1),PetscErrorCode,(Cstring,Cstring,Cint,Cint,Cint,Cint,Ptr{Ptr{PetscObject}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMonitorLGResidualNorm(arg0::Type{Float64},arg1::SNES,arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{PetscObject},StridedArray{PetscObject},Ptr{Void}))
    ccall((:SNESMonitorLGResidualNorm,petsc1),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{PetscObject}),arg1,arg2,PetscReal,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMonitorLGDestroy(arg0::Type{Float64},arg1::Union(Ptr{Ptr{PetscObject}},StridedArray{Ptr{PetscObject}},Ptr{Void}))
    ccall((:SNESMonitorLGDestroy,petsc1),PetscErrorCode,(Ptr{Ptr{PetscObject}},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMonitorLGRange(arg0::Type{Float64},arg1::SNES,arg2::Integer,PetscReal::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESMonitorLGRange,petsc1),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetApplicationContext(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESSetApplicationContext,petsc1),PetscErrorCode,(SNES,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetApplicationContext(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESGetApplicationContext,petsc1),PetscErrorCode,(SNES,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetComputeApplicationContext(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESSetComputeApplicationContext,petsc1),PetscErrorCode,(SNES,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESPythonSetType(arg0::Type{Float64},arg1::SNES,arg2::Union(ByteString,Symbol))
    ccall((:SNESPythonSetType,petsc1),PetscErrorCode,(SNES,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetFunctionDomainError(arg0::Type{Float64},arg1::SNES)
    ccall((:SNESSetFunctionDomainError,petsc1),PetscErrorCode,(SNES,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetFunctionDomainError(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:SNESGetFunctionDomainError,petsc1),PetscErrorCode,(SNES,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetConvergenceTest(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESSetConvergenceTest,petsc1),PetscErrorCode,(SNES,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function SNESConvergedDefault(arg0::Type{Float64},arg1::SNES,arg2::Integer,PetscReal::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{SNESConvergedReason},StridedArray{SNESConvergedReason},Ptr{Void}),arg6::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESConvergedDefault,petsc1),PetscErrorCode,(SNES,PetscInt,Cint,Cint,Cint,Ptr{SNESConvergedReason},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function SNESConvergedSkip(arg0::Type{Float64},arg1::SNES,arg2::Integer,PetscReal::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{SNESConvergedReason},StridedArray{SNESConvergedReason},Ptr{Void}),arg6::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESConvergedSkip,petsc1),PetscErrorCode,(SNES,PetscInt,Cint,Cint,Cint,Ptr{SNESConvergedReason},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetConvergedReason(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{SNESConvergedReason},StridedArray{SNESConvergedReason},Ptr{Void}))
    ccall((:SNESGetConvergedReason,petsc1),PetscErrorCode,(SNES,Ptr{SNESConvergedReason}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetFunction(arg1::SNES,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg4::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:SNESGetFunction,petsc1),PetscErrorCode,(SNES,Ptr{Vec{Float64}},Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function SNESComputeFunction(arg1::SNES,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:SNESComputeFunction,petsc1),PetscErrorCode,(SNES,Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetJacobian(arg1::SNES,arg2::Mat{Float64},arg3::Mat{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESSetJacobian,petsc1),PetscErrorCode,(SNES,Mat{Float64},Mat{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetJacobian(arg1::SNES,arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg4::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg5::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:SNESGetJacobian,petsc1),PetscErrorCode,(SNES,Ptr{Mat{Float64}},Ptr{Mat{Float64}},Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function SNESObjectiveComputeFunctionDefaultFD(arg1::SNES,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESObjectiveComputeFunctionDefaultFD,petsc1),PetscErrorCode,(SNES,Vec{Float64},Vec{Float64},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function SNESComputeJacobianDefault(arg1::SNES,arg2::Vec{Float64},arg3::Mat{Float64},arg4::Mat{Float64},arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESComputeJacobianDefault,petsc1),PetscErrorCode,(SNES,Vec{Float64},Mat{Float64},Mat{Float64},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function SNESComputeJacobianDefaultColor(arg1::SNES,arg2::Vec{Float64},arg3::Mat{Float64},arg4::Mat{Float64},arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESComputeJacobianDefaultColor,petsc1),PetscErrorCode,(SNES,Vec{Float64},Mat{Float64},Mat{Float64},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetComputeInitialGuess(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESSetComputeInitialGuess,petsc1),PetscErrorCode,(SNES,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetPicard(arg1::SNES,arg2::Vec{Float64},arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Mat{Float64},arg5::Mat{Float64},arg6::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg7::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESSetPicard,petsc1),PetscErrorCode,(SNES,Vec{Float64},Ptr{Void},Mat{Float64},Mat{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetPicard(arg1::SNES,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg5::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg6::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg7::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:SNESGetPicard,petsc1),PetscErrorCode,(SNES,Ptr{Vec{Float64}},Ptr{Ptr{Void}},Ptr{Mat{Float64}},Ptr{Mat{Float64}},Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetInitialFunction(arg1::SNES,arg2::Vec{Float64})
    ccall((:SNESSetInitialFunction,petsc1),PetscErrorCode,(SNES,Vec{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetObjective(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESSetObjective,petsc1),PetscErrorCode,(SNES,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetObjective(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:SNESGetObjective,petsc1),PetscErrorCode,(SNES,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESComputeObjective(arg1::SNES,arg2::Vec{Float64},arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:SNESComputeObjective,petsc1),PetscErrorCode,(SNES,Vec{Float64},Ptr{Cint}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetNormSchedule(arg0::Type{Float64},arg1::SNES,arg2::SNESNormSchedule)
    ccall((:SNESSetNormSchedule,petsc1),PetscErrorCode,(SNES,SNESNormSchedule),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetNormSchedule(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{SNESNormSchedule},StridedArray{SNESNormSchedule},Ptr{Void}))
    ccall((:SNESGetNormSchedule,petsc1),PetscErrorCode,(SNES,Ptr{SNESNormSchedule}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetFunctionType(arg0::Type{Float64},arg1::SNES,arg2::SNESFunctionType)
    ccall((:SNESSetFunctionType,petsc1),PetscErrorCode,(SNES,SNESFunctionType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetFunctionType(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{SNESFunctionType},StridedArray{SNESFunctionType},Ptr{Void}))
    ccall((:SNESGetFunctionType,petsc1),PetscErrorCode,(SNES,Ptr{SNESFunctionType}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetNGS(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESSetNGS,petsc1),PetscErrorCode,(SNES,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetNGS(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:SNESGetNGS,petsc1),PetscErrorCode,(SNES,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetUseNGS(arg0::Type{Float64},arg1::SNES,arg2::PetscBool)
    ccall((:SNESSetUseNGS,petsc1),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetUseNGS(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:SNESGetUseNGS,petsc1),PetscErrorCode,(SNES,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESComputeNGS(arg1::SNES,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:SNESComputeNGS,petsc1),PetscErrorCode,(SNES,Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESNGSSetSweeps(arg0::Type{Float64},arg1::SNES,arg2::Integer)
    ccall((:SNESNGSSetSweeps,petsc1),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESNGSGetSweeps(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:SNESNGSGetSweeps,petsc1),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESNGSSetTolerances(arg0::Type{Float64},arg1::SNES,PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:SNESNGSSetTolerances,petsc1),PetscErrorCode,(SNES,Cint,Cint,Cint,PetscInt),arg1,PetscReal,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function SNESNGSGetTolerances(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:SNESNGSGetTolerances,petsc1),PetscErrorCode,(SNES,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function SNESUpdateCheckJacobian(arg0::Type{Float64},arg1::SNES,arg2::Integer)
    ccall((:SNESUpdateCheckJacobian,petsc1),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESShellGetContext(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:SNESShellGetContext,petsc1),PetscErrorCode,(SNES,Ptr{Ptr{Void}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESShellSetContext(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESShellSetContext,petsc1),PetscErrorCode,(SNES,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESShellSetSolve(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESShellSetSolve,petsc1),PetscErrorCode,(SNES,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{SNESLineSearch},StridedArray{SNESLineSearch},Ptr{Void}))
    ccall((:SNESLineSearchCreate,petsc1),PetscErrorCode,(comm_type,Ptr{SNESLineSearch}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchReset(arg0::Type{Float64},arg1::SNESLineSearch)
    ccall((:SNESLineSearchReset,petsc1),PetscErrorCode,(SNESLineSearch,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchView(arg1::SNESLineSearch,arg2::PetscViewer{Float64})
    ccall((:SNESLineSearchView,petsc1),PetscErrorCode,(SNESLineSearch,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchDestroy(arg0::Type{Float64},arg1::Union(Ptr{SNESLineSearch},StridedArray{SNESLineSearch},Ptr{Void}))
    ccall((:SNESLineSearchDestroy,petsc1),PetscErrorCode,(Ptr{SNESLineSearch},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchSetType(arg0::Type{Float64},arg1::SNESLineSearch,arg2::SNESLineSearchType)
    ccall((:SNESLineSearchSetType,petsc1),PetscErrorCode,(SNESLineSearch,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchSetFromOptions(arg0::Type{Float64},arg1::SNESLineSearch)
    ccall((:SNESLineSearchSetFromOptions,petsc1),PetscErrorCode,(SNESLineSearch,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchSetFunction(arg0::Type{Float64},arg1::SNESLineSearch,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESLineSearchSetFunction,petsc1),PetscErrorCode,(SNESLineSearch,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchSetUp(arg0::Type{Float64},arg1::SNESLineSearch)
    ccall((:SNESLineSearchSetUp,petsc1),PetscErrorCode,(SNESLineSearch,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchApply(arg1::SNESLineSearch,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Vec{Float64})
    ccall((:SNESLineSearchApply,petsc1),PetscErrorCode,(SNESLineSearch,Vec{Float64},Vec{Float64},Ptr{Cint},Vec{Float64}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchPreCheck(arg1::SNESLineSearch,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:SNESLineSearchPreCheck,petsc1),PetscErrorCode,(SNESLineSearch,Vec{Float64},Vec{Float64},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchPostCheck(arg1::SNESLineSearch,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64},arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg6::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:SNESLineSearchPostCheck,petsc1),PetscErrorCode,(SNESLineSearch,Vec{Float64},Vec{Float64},Vec{Float64},Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchSetWorkVecs(arg0::Type{Float64},arg1::SNESLineSearch,arg2::Integer)
    ccall((:SNESLineSearchSetWorkVecs,petsc1),PetscErrorCode,(SNESLineSearch,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchSetPreCheck(arg0::Type{Float64},arg1::SNESLineSearch,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),ctx::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESLineSearchSetPreCheck,petsc1),PetscErrorCode,(SNESLineSearch,Ptr{Void},Ptr{Void}),arg1,arg2,ctx)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchSetPostCheck(arg0::Type{Float64},arg1::SNESLineSearch,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),ctx::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESLineSearchSetPostCheck,petsc1),PetscErrorCode,(SNESLineSearch,Ptr{Void},Ptr{Void}),arg1,arg2,ctx)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchGetPreCheck(arg0::Type{Float64},arg1::SNESLineSearch,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),ctx::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:SNESLineSearchGetPreCheck,petsc1),PetscErrorCode,(SNESLineSearch,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,ctx)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchGetPostCheck(arg0::Type{Float64},arg1::SNESLineSearch,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),ctx::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:SNESLineSearchGetPostCheck,petsc1),PetscErrorCode,(SNESLineSearch,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,ctx)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchSetVIFunctions(arg0::Type{Float64},arg1::SNESLineSearch,arg2::SNESLineSearchVIProjectFunc,arg3::SNESLineSearchVINormFunc)
    ccall((:SNESLineSearchSetVIFunctions,petsc1),PetscErrorCode,(SNESLineSearch,SNESLineSearchVIProjectFunc,SNESLineSearchVINormFunc),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchGetVIFunctions(arg0::Type{Float64},arg1::SNESLineSearch,arg2::Union(Ptr{SNESLineSearchVIProjectFunc},StridedArray{SNESLineSearchVIProjectFunc},Ptr{Void}),arg3::Union(Ptr{SNESLineSearchVINormFunc},StridedArray{SNESLineSearchVINormFunc},Ptr{Void}))
    ccall((:SNESLineSearchGetVIFunctions,petsc1),PetscErrorCode,(SNESLineSearch,Ptr{SNESLineSearchVIProjectFunc},Ptr{SNESLineSearchVINormFunc}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchSetSNES(arg0::Type{Float64},arg1::SNESLineSearch,arg2::SNES)
    ccall((:SNESLineSearchSetSNES,petsc1),PetscErrorCode,(SNESLineSearch,SNES),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchGetSNES(arg0::Type{Float64},arg1::SNESLineSearch,arg2::Union(Ptr{SNES},StridedArray{SNES},Ptr{Void}))
    ccall((:SNESLineSearchGetSNES,petsc1),PetscErrorCode,(SNESLineSearch,Ptr{SNES}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchGetTolerances(arg0::Type{Float64},arg1::SNESLineSearch,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg7::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:SNESLineSearchGetTolerances,petsc1),PetscErrorCode,(SNESLineSearch,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchSetTolerances(arg0::Type{Float64},arg1::SNESLineSearch,PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer)
    ccall((:SNESLineSearchSetTolerances,petsc1),PetscErrorCode,(SNESLineSearch,Cint,Cint,Cint,Cint,Cint,PetscInt),arg1,PetscReal,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchPreCheckPicard(arg1::SNESLineSearch,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESLineSearchPreCheckPicard,petsc1),PetscErrorCode,(SNESLineSearch,Vec{Float64},Vec{Float64},Ptr{PetscBool},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchGetLambda(arg0::Type{Float64},arg1::SNESLineSearch,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:SNESLineSearchGetLambda,petsc1),PetscErrorCode,(SNESLineSearch,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchSetLambda(arg0::Type{Float64},arg1::SNESLineSearch,PetscReal::Integer)
    ccall((:SNESLineSearchSetLambda,petsc1),PetscErrorCode,(SNESLineSearch,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchGetDamping(arg0::Type{Float64},arg1::SNESLineSearch,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:SNESLineSearchGetDamping,petsc1),PetscErrorCode,(SNESLineSearch,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchSetDamping(arg0::Type{Float64},arg1::SNESLineSearch,PetscReal::Integer)
    ccall((:SNESLineSearchSetDamping,petsc1),PetscErrorCode,(SNESLineSearch,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchGetOrder(arg0::Type{Float64},arg1::SNESLineSearch,order::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:SNESLineSearchGetOrder,petsc1),PetscErrorCode,(SNESLineSearch,Ptr{PetscInt}),arg1,order)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchSetOrder(arg0::Type{Float64},arg1::SNESLineSearch,order::Integer)
    ccall((:SNESLineSearchSetOrder,petsc1),PetscErrorCode,(SNESLineSearch,PetscInt),arg1,order)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchGetReason(arg0::Type{Float64},arg1::SNESLineSearch,arg2::Union(Ptr{SNESLineSearchReason},StridedArray{SNESLineSearchReason},Ptr{Void}))
    ccall((:SNESLineSearchGetReason,petsc1),PetscErrorCode,(SNESLineSearch,Ptr{SNESLineSearchReason}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchSetReason(arg0::Type{Float64},arg1::SNESLineSearch,arg2::SNESLineSearchReason)
    ccall((:SNESLineSearchSetReason,petsc1),PetscErrorCode,(SNESLineSearch,SNESLineSearchReason),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchGetVecs(arg1::SNESLineSearch,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg4::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg5::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg6::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:SNESLineSearchGetVecs,petsc1),PetscErrorCode,(SNESLineSearch,Ptr{Vec{Float64}},Ptr{Vec{Float64}},Ptr{Vec{Float64}},Ptr{Vec{Float64}},Ptr{Vec{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchSetVecs(arg1::SNESLineSearch,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64},arg5::Vec{Float64},arg6::Vec{Float64})
    ccall((:SNESLineSearchSetVecs,petsc1),PetscErrorCode,(SNESLineSearch,Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchGetNorms(arg0::Type{Float64},arg1::SNESLineSearch,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:SNESLineSearchGetNorms,petsc1),PetscErrorCode,(SNESLineSearch,Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchSetNorms(arg0::Type{Float64},arg1::SNESLineSearch,PetscReal::Integer,arg2::Integer,arg3::Integer)
    ccall((:SNESLineSearchSetNorms,petsc1),PetscErrorCode,(SNESLineSearch,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchComputeNorms(arg0::Type{Float64},arg1::SNESLineSearch)
    ccall((:SNESLineSearchComputeNorms,petsc1),PetscErrorCode,(SNESLineSearch,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchSetComputeNorms(arg0::Type{Float64},arg1::SNESLineSearch,arg2::PetscBool)
    ccall((:SNESLineSearchSetComputeNorms,petsc1),PetscErrorCode,(SNESLineSearch,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchSetMonitor(arg0::Type{Float64},arg1::SNESLineSearch,arg2::PetscBool)
    ccall((:SNESLineSearchSetMonitor,petsc1),PetscErrorCode,(SNESLineSearch,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchGetMonitor(arg1::SNESLineSearch,arg2::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    ccall((:SNESLineSearchGetMonitor,petsc1),PetscErrorCode,(SNESLineSearch,Ptr{PetscViewer{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchAppendOptionsPrefix(arg0::Type{Float64},arg1::SNESLineSearch,prefix::Union(ByteString,Symbol))
    ccall((:SNESLineSearchAppendOptionsPrefix,petsc1),PetscErrorCode,(SNESLineSearch,Cstring),arg1,prefix)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchGetOptionsPrefix(arg0::Type{Float64},arg1::SNESLineSearch,prefix::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:SNESLineSearchGetOptionsPrefix,petsc1),PetscErrorCode,(SNESLineSearch,Ptr{Ptr{Uint8}}),arg1,prefix)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchShellSetUserFunc(arg0::Type{Float64},arg1::SNESLineSearch,arg2::SNESLineSearchUserFunc,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESLineSearchShellSetUserFunc,petsc1),PetscErrorCode,(SNESLineSearch,SNESLineSearchUserFunc,Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchShellGetUserFunc(arg0::Type{Float64},arg1::SNESLineSearch,arg2::Union(Ptr{SNESLineSearchUserFunc},StridedArray{SNESLineSearchUserFunc},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:SNESLineSearchShellGetUserFunc,petsc1),PetscErrorCode,(SNESLineSearch,Ptr{SNESLineSearchUserFunc},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchBTSetAlpha(arg0::Type{Float64},arg1::SNESLineSearch,PetscReal::Integer)
    ccall((:SNESLineSearchBTSetAlpha,petsc1),PetscErrorCode,(SNESLineSearch,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function SNESLineSearchBTGetAlpha(arg0::Type{Float64},arg1::SNESLineSearch,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:SNESLineSearchBTGetAlpha,petsc1),PetscErrorCode,(SNESLineSearch,Ptr{Cint}),arg1,arg2)
end 
=#
function SNESLineSearchRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:SNESLineSearchRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function SNESVISetVariableBounds(arg1::SNES,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:SNESVISetVariableBounds,petsc1),PetscErrorCode,(SNES,Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESVISetComputeVariableBounds(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESVISetComputeVariableBounds,petsc1),PetscErrorCode,(SNES,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESVIGetInactiveSet(arg1::SNES,arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:SNESVIGetInactiveSet,petsc1),PetscErrorCode,(SNES,Ptr{IS{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESVIGetActiveSetIS(arg1::SNES,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:SNESVIGetActiveSetIS,petsc1),PetscErrorCode,(SNES,Vec{Float64},Vec{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function SNESVIComputeInactiveSetFnorm(arg1::SNES,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:SNESVIComputeInactiveSetFnorm,petsc1),PetscErrorCode,(SNES,Vec{Float64},Vec{Float64},Ptr{Cint}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function SNESVISetRedundancyCheck(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESVISetRedundancyCheck,petsc1),PetscErrorCode,(SNES,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESTestLocalMin(arg0::Type{Float64},arg1::SNES)
    ccall((:SNESTestLocalMin,petsc1),PetscErrorCode,(SNES,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function SNESComputeJacobian(arg1::SNES,arg2::Vec{Float64},arg3::Mat{Float64},arg4::Mat{Float64})
    ccall((:SNESComputeJacobian,petsc1),PetscErrorCode,(SNES,Vec{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetDM(arg0::Type{Float64},arg1::SNES,arg2::DM)
    ccall((:SNESSetDM,petsc1),PetscErrorCode,(SNES,DM),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetDM(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:SNESGetDM,petsc1),PetscErrorCode,(SNES,Ptr{DM}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetNPC(arg0::Type{Float64},arg1::SNES,arg2::SNES)
    ccall((:SNESSetNPC,petsc1),PetscErrorCode,(SNES,SNES),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetNPC(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{SNES},StridedArray{SNES},Ptr{Void}))
    ccall((:SNESGetNPC,petsc1),PetscErrorCode,(SNES,Ptr{SNES}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESHasNPC(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:SNESHasNPC,petsc1),PetscErrorCode,(SNES,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESApplyNPC(arg1::SNES,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64})
    ccall((:SNESApplyNPC,petsc1),PetscErrorCode,(SNES,Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetNPCFunction(arg1::SNES,arg2::Vec{Float64},arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:SNESGetNPCFunction,petsc1),PetscErrorCode,(SNES,Vec{Float64},Ptr{Cint}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESComputeFunctionDefaultNPC(arg1::SNES,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:SNESComputeFunctionDefaultNPC,petsc1),PetscErrorCode,(SNES,Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetNPCSide(arg0::Type{Float64},arg1::SNES,arg2::PCSide)
    ccall((:SNESSetNPCSide,petsc1),PetscErrorCode,(SNES,PCSide),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetNPCSide(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{PCSide},StridedArray{PCSide},Ptr{Void}))
    ccall((:SNESGetNPCSide,petsc1),PetscErrorCode,(SNES,Ptr{PCSide}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESSetLineSearch(arg0::Type{Float64},arg1::SNES,arg2::SNESLineSearch)
    ccall((:SNESSetLineSearch,petsc1),PetscErrorCode,(SNES,SNESLineSearch),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESGetLineSearch(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{SNESLineSearch},StridedArray{SNESLineSearch},Ptr{Void}))
    ccall((:SNESGetLineSearch,petsc1),PetscErrorCode,(SNES,Ptr{SNESLineSearch}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESRestrictHookAdd(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESRestrictHookAdd,petsc1),PetscErrorCode,(SNES,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESRestrictHooksRun(arg0::Type{Float64},arg1::SNES,arg2::SNES)
    ccall((:SNESRestrictHooksRun,petsc1),PetscErrorCode,(SNES,SNES),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMSNESSetFunction(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMSNESSetFunction,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMSNESGetFunction(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:DMSNESGetFunction,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMSNESSetNGS(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMSNESSetNGS,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMSNESGetNGS(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:DMSNESGetNGS,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMSNESSetJacobian(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMSNESSetJacobian,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMSNESGetJacobian(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:DMSNESGetJacobian,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMSNESSetPicard(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMSNESSetPicard,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMSNESGetPicard(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg4::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:DMSNESGetPicard,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMSNESSetObjective(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMSNESSetObjective,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMSNESGetObjective(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:DMSNESGetObjective,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASNESSetFunctionLocal(arg0::Type{Float64},arg1::DM,arg2::InsertMode,arg3::DMDASNESFunction,arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMDASNESSetFunctionLocal,petsc1),PetscErrorCode,(DM,InsertMode,DMDASNESFunction,Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASNESSetJacobianLocal(arg0::Type{Float64},arg1::DM,arg2::DMDASNESJacobian,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMDASNESSetJacobianLocal,petsc1),PetscErrorCode,(DM,DMDASNESJacobian,Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASNESSetObjectiveLocal(arg0::Type{Float64},arg1::DM,arg2::DMDASNESObjective,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMDASNESSetObjectiveLocal,petsc1),PetscErrorCode,(DM,DMDASNESObjective,Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDASNESSetPicardLocal(arg0::Type{Float64},arg1::DM,arg2::InsertMode,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMDASNESSetPicardLocal,petsc1),PetscErrorCode,(DM,InsertMode,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSNESGetGeometryFEM(arg1::DM,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:DMPlexSNESGetGeometryFEM,petsc1),PetscErrorCode,(DM,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSNESGetGeometryFVM(arg1::DM,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:DMPlexSNESGetGeometryFVM,petsc1),PetscErrorCode,(DM,Ptr{Vec{Float64}},Ptr{Vec{Float64}},Ptr{Cint}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexSNESGetGradientDM(arg0::Type{Float64},arg1::DM,arg2::PetscFV,arg3::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:DMPlexSNESGetGradientDM,petsc1),PetscErrorCode,(DM,PetscFV,Ptr{DM}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetCellFields(arg1::DM,arg2::Integer,arg3::Integer,arg4::Vec{Float64},arg5::Vec{Float64},arg6::Vec{Float64},arg7::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}),arg8::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}),arg9::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:DMPlexGetCellFields,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,Vec{Float64},Vec{Float64},Vec{Float64},Ptr{Ptr{Float64}},Ptr{Ptr{Float64}},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexRestoreCellFields(arg1::DM,arg2::Integer,arg3::Integer,arg4::Vec{Float64},arg5::Vec{Float64},arg6::Vec{Float64},arg7::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}),arg8::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}),arg9::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:DMPlexRestoreCellFields,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,Vec{Float64},Vec{Float64},Vec{Float64},Ptr{Ptr{Float64}},Ptr{Ptr{Float64}},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetFaceFields(arg1::DM,arg2::Integer,arg3::Integer,arg4::Vec{Float64},arg5::Vec{Float64},arg6::Vec{Float64},arg7::Vec{Float64},arg8::Vec{Float64},arg9::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}),arg10::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:DMPlexGetFaceFields,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Ptr{Ptr{Float64}},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexRestoreFaceFields(arg1::DM,arg2::Integer,arg3::Integer,arg4::Vec{Float64},arg5::Vec{Float64},arg6::Vec{Float64},arg7::Vec{Float64},arg8::Vec{Float64},arg9::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}),arg10::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:DMPlexRestoreFaceFields,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Ptr{Ptr{Float64}},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexGetFaceGeometry(arg1::DM,arg2::Integer,arg3::Integer,arg4::Vec{Float64},arg5::Vec{Float64},arg6::Union(Ptr{Ptr{PetscFVFaceGeom}},StridedArray{Ptr{PetscFVFaceGeom}},Ptr{Void}),arg7::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}))
    ccall((:DMPlexGetFaceGeometry,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,Vec{Float64},Vec{Float64},Ptr{Ptr{PetscFVFaceGeom}},Ptr{Ptr{Cint}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexRestoreFaceGeometry(arg1::DM,arg2::Integer,arg3::Integer,arg4::Vec{Float64},arg5::Vec{Float64},arg6::Union(Ptr{Ptr{PetscFVFaceGeom}},StridedArray{Ptr{PetscFVFaceGeom}},Ptr{Void}),arg7::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}))
    ccall((:DMPlexRestoreFaceGeometry,petsc1),PetscErrorCode,(DM,PetscInt,PetscInt,Vec{Float64},Vec{Float64},Ptr{Ptr{PetscFVFaceGeom}},Ptr{Ptr{Cint}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function DMSNESSetFunctionLocal(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMSNESSetFunctionLocal,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMSNESSetJacobianLocal(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMSNESSetJacobianLocal,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMultiblockSetFields(arg0::Type{Float64},arg1::SNES,arg2::Union(ByteString,Symbol),arg3::Integer,arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:SNESMultiblockSetFields,petsc1),PetscErrorCode,(SNES,Cstring,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMultiblockSetIS(arg1::SNES,arg2::Union(ByteString,Symbol),arg3::IS{Float64})
    ccall((:SNESMultiblockSetIS,petsc1),PetscErrorCode,(SNES,Cstring,IS{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMultiblockSetBlockSize(arg0::Type{Float64},arg1::SNES,arg2::Integer)
    ccall((:SNESMultiblockSetBlockSize,petsc1),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESMultiblockSetType(arg0::Type{Float64},arg1::SNES,arg2::PCCompositeType)
    ccall((:SNESMultiblockSetType,petsc1),PetscErrorCode,(SNES,PCCompositeType),arg1,arg2)
end 
=#
function SNESMSRegister(arg0::Type{Float64},arg1::SNESMSType,arg2::Integer,arg3::Integer,PetscReal::Integer,arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:SNESMSRegister,petsc1),PetscErrorCode,(Cstring,PetscInt,PetscInt,Cint,Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,PetscReal,arg4,arg5,arg6)
    return err
end

#= skipping function with undefined symbols: 
 function SNESMSSetType(arg0::Type{Float64},arg1::SNES,arg2::SNESMSType)
    ccall((:SNESMSSetType,petsc1),PetscErrorCode,(SNES,Cstring),arg1,arg2)
end 
=#
function SNESMSFinalizePackage(arg0::Type{Float64})
    err = ccall((:SNESMSFinalizePackage,petsc1),PetscErrorCode,())
    return err
end

function SNESMSInitializePackage(arg0::Type{Float64})
    err = ccall((:SNESMSInitializePackage,petsc1),PetscErrorCode,())
    return err
end

function SNESMSRegisterDestroy(arg0::Type{Float64})
    err = ccall((:SNESMSRegisterDestroy,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function SNESNGMRESSetRestartType(arg0::Type{Float64},arg1::SNES,arg2::SNESNGMRESRestartType)
    ccall((:SNESNGMRESSetRestartType,petsc1),PetscErrorCode,(SNES,SNESNGMRESRestartType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESNGMRESSetSelectType(arg0::Type{Float64},arg1::SNES,arg2::SNESNGMRESSelectType)
    ccall((:SNESNGMRESSetSelectType,petsc1),PetscErrorCode,(SNES,SNESNGMRESSelectType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESNCGSetType(arg0::Type{Float64},arg1::SNES,arg2::SNESNCGType)
    ccall((:SNESNCGSetType,petsc1),PetscErrorCode,(SNES,SNESNCGType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESQNSetType(arg0::Type{Float64},arg1::SNES,arg2::SNESQNType)
    ccall((:SNESQNSetType,petsc1),PetscErrorCode,(SNES,SNESQNType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESQNSetScaleType(arg0::Type{Float64},arg1::SNES,arg2::SNESQNScaleType)
    ccall((:SNESQNSetScaleType,petsc1),PetscErrorCode,(SNES,SNESQNScaleType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESQNSetRestartType(arg0::Type{Float64},arg1::SNES,arg2::SNESQNRestartType)
    ccall((:SNESQNSetRestartType,petsc1),PetscErrorCode,(SNES,SNESQNRestartType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESNASMGetType(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{PCASMType},StridedArray{PCASMType},Ptr{Void}))
    ccall((:SNESNASMGetType,petsc1),PetscErrorCode,(SNES,Ptr{PCASMType}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESNASMSetType(arg0::Type{Float64},arg1::SNES,arg2::PCASMType)
    ccall((:SNESNASMSetType,petsc1),PetscErrorCode,(SNES,PCASMType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESNASMGetSubdomains(arg1::SNES,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{SNES}},StridedArray{Ptr{SNES}},Ptr{Void}),arg4::Union(Ptr{Ptr{VecScatter{Float64}}},StridedArray{Ptr{VecScatter{Float64}}},Ptr{Void}),arg5::Union(Ptr{Ptr{VecScatter{Float64}}},StridedArray{Ptr{VecScatter{Float64}}},Ptr{Void}),arg6::Union(Ptr{Ptr{VecScatter{Float64}}},StridedArray{Ptr{VecScatter{Float64}}},Ptr{Void}))
    ccall((:SNESNASMGetSubdomains,petsc1),PetscErrorCode,(SNES,Ptr{PetscInt},Ptr{Ptr{SNES}},Ptr{Ptr{VecScatter{Float64}}},Ptr{Ptr{VecScatter{Float64}}},Ptr{Ptr{VecScatter{Float64}}}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function SNESNASMSetSubdomains(arg1::SNES,arg2::Integer,arg3::Union(Ptr{SNES},StridedArray{SNES},Ptr{Void}),arg4::Union(Ptr{VecScatter{Float64}},StridedArray{VecScatter{Float64}},Ptr{Void}),arg5::Union(Ptr{VecScatter{Float64}},StridedArray{VecScatter{Float64}},Ptr{Void}),arg6::Union(Ptr{VecScatter{Float64}},StridedArray{VecScatter{Float64}},Ptr{Void}))
    ccall((:SNESNASMSetSubdomains,petsc1),PetscErrorCode,(SNES,PetscInt,Ptr{SNES},Ptr{VecScatter{Float64}},Ptr{VecScatter{Float64}},Ptr{VecScatter{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function SNESNASMSetDamping(arg0::Type{Float64},arg1::SNES,PetscReal::Integer)
    ccall((:SNESNASMSetDamping,petsc1),PetscErrorCode,(SNES,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function SNESNASMGetDamping(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:SNESNASMGetDamping,petsc1),PetscErrorCode,(SNES,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESNASMGetSubdomainVecs(arg1::SNES,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{Vec{Float64}}},StridedArray{Ptr{Vec{Float64}}},Ptr{Void}),arg4::Union(Ptr{Ptr{Vec{Float64}}},StridedArray{Ptr{Vec{Float64}}},Ptr{Void}),arg5::Union(Ptr{Ptr{Vec{Float64}}},StridedArray{Ptr{Vec{Float64}}},Ptr{Void}),arg6::Union(Ptr{Ptr{Vec{Float64}}},StridedArray{Ptr{Vec{Float64}}},Ptr{Void}))
    ccall((:SNESNASMGetSubdomainVecs,petsc1),PetscErrorCode,(SNES,Ptr{PetscInt},Ptr{Ptr{Vec{Float64}}},Ptr{Ptr{Vec{Float64}}},Ptr{Ptr{Vec{Float64}}},Ptr{Ptr{Vec{Float64}}}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function SNESNASMSetComputeFinalJacobian(arg0::Type{Float64},arg1::SNES,arg2::PetscBool)
    ccall((:SNESNASMSetComputeFinalJacobian,petsc1),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESCompositeSetType(arg0::Type{Float64},arg1::SNES,arg2::SNESCompositeType)
    ccall((:SNESCompositeSetType,petsc1),PetscErrorCode,(SNES,SNESCompositeType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESCompositeAddSNES(arg0::Type{Float64},arg1::SNES,arg2::SNESType)
    ccall((:SNESCompositeAddSNES,petsc1),PetscErrorCode,(SNES,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESCompositeGetSNES(arg0::Type{Float64},arg1::SNES,arg2::Integer,arg3::Union(Ptr{SNES},StridedArray{SNES},Ptr{Void}))
    ccall((:SNESCompositeGetSNES,petsc1),PetscErrorCode,(SNES,PetscInt,Ptr{SNES}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESCompositeGetNumber(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:SNESCompositeGetNumber,petsc1),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESCompositeSetDamping(arg0::Type{Float64},arg1::SNES,arg2::Integer,PetscReal::Integer)
    ccall((:SNESCompositeSetDamping,petsc1),PetscErrorCode,(SNES,PetscInt,Cint),arg1,arg2,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASSetType(arg0::Type{Float64},arg1::SNES,arg2::SNESFASType)
    ccall((:SNESFASSetType,petsc1),PetscErrorCode,(SNES,SNESFASType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASGetType(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{SNESFASType},StridedArray{SNESFASType},Ptr{Void}))
    ccall((:SNESFASGetType,petsc1),PetscErrorCode,(SNES,Ptr{SNESFASType}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASSetLevels(arg0::Type{Float64},arg1::SNES,arg2::Integer,arg3::Union(Ptr{MPI_Comm},StridedArray{MPI_Comm},Ptr{Void}))
    ccall((:SNESFASSetLevels,petsc1),PetscErrorCode,(SNES,PetscInt,Ptr{comm_type}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASGetLevels(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:SNESFASGetLevels,petsc1),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASGetCycleSNES(arg0::Type{Float64},arg1::SNES,arg2::Integer,arg3::Union(Ptr{SNES},StridedArray{SNES},Ptr{Void}))
    ccall((:SNESFASGetCycleSNES,petsc1),PetscErrorCode,(SNES,PetscInt,Ptr{SNES}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASSetNumberSmoothUp(arg0::Type{Float64},arg1::SNES,arg2::Integer)
    ccall((:SNESFASSetNumberSmoothUp,petsc1),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASSetNumberSmoothDown(arg0::Type{Float64},arg1::SNES,arg2::Integer)
    ccall((:SNESFASSetNumberSmoothDown,petsc1),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASSetCycles(arg0::Type{Float64},arg1::SNES,arg2::Integer)
    ccall((:SNESFASSetCycles,petsc1),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASSetMonitor(arg0::Type{Float64},arg1::SNES,arg2::PetscBool)
    ccall((:SNESFASSetMonitor,petsc1),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASSetLog(arg0::Type{Float64},arg1::SNES,arg2::PetscBool)
    ccall((:SNESFASSetLog,petsc1),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASSetGalerkin(arg0::Type{Float64},arg1::SNES,arg2::PetscBool)
    ccall((:SNESFASSetGalerkin,petsc1),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASGetGalerkin(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:SNESFASGetGalerkin,petsc1),PetscErrorCode,(SNES,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASCycleGetSmoother(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{SNES},StridedArray{SNES},Ptr{Void}))
    ccall((:SNESFASCycleGetSmoother,petsc1),PetscErrorCode,(SNES,Ptr{SNES}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASCycleGetSmootherUp(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{SNES},StridedArray{SNES},Ptr{Void}))
    ccall((:SNESFASCycleGetSmootherUp,petsc1),PetscErrorCode,(SNES,Ptr{SNES}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASCycleGetSmootherDown(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{SNES},StridedArray{SNES},Ptr{Void}))
    ccall((:SNESFASCycleGetSmootherDown,petsc1),PetscErrorCode,(SNES,Ptr{SNES}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASCycleGetCorrection(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{SNES},StridedArray{SNES},Ptr{Void}))
    ccall((:SNESFASCycleGetCorrection,petsc1),PetscErrorCode,(SNES,Ptr{SNES}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASCycleGetInterpolation(arg1::SNES,arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:SNESFASCycleGetInterpolation,petsc1),PetscErrorCode,(SNES,Ptr{Mat{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASCycleGetRestriction(arg1::SNES,arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:SNESFASCycleGetRestriction,petsc1),PetscErrorCode,(SNES,Ptr{Mat{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASCycleGetInjection(arg1::SNES,arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:SNESFASCycleGetInjection,petsc1),PetscErrorCode,(SNES,Ptr{Mat{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASCycleGetRScale(arg1::SNES,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:SNESFASCycleGetRScale,petsc1),PetscErrorCode,(SNES,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASCycleSetCycles(arg0::Type{Float64},arg1::SNES,arg2::Integer)
    ccall((:SNESFASCycleSetCycles,petsc1),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASCycleIsFine(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:SNESFASCycleIsFine,petsc1),PetscErrorCode,(SNES,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASSetInterpolation(arg1::SNES,arg2::Integer,arg3::Mat{Float64})
    ccall((:SNESFASSetInterpolation,petsc1),PetscErrorCode,(SNES,PetscInt,Mat{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASGetInterpolation(arg1::SNES,arg2::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:SNESFASGetInterpolation,petsc1),PetscErrorCode,(SNES,PetscInt,Ptr{Mat{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASSetRestriction(arg1::SNES,arg2::Integer,arg3::Mat{Float64})
    ccall((:SNESFASSetRestriction,petsc1),PetscErrorCode,(SNES,PetscInt,Mat{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASGetRestriction(arg1::SNES,arg2::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:SNESFASGetRestriction,petsc1),PetscErrorCode,(SNES,PetscInt,Ptr{Mat{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASSetInjection(arg1::SNES,arg2::Integer,arg3::Mat{Float64})
    ccall((:SNESFASSetInjection,petsc1),PetscErrorCode,(SNES,PetscInt,Mat{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASGetInjection(arg1::SNES,arg2::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:SNESFASGetInjection,petsc1),PetscErrorCode,(SNES,PetscInt,Ptr{Mat{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASSetRScale(arg1::SNES,arg2::Integer,arg3::Vec{Float64})
    ccall((:SNESFASSetRScale,petsc1),PetscErrorCode,(SNES,PetscInt,Vec{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASGetRScale(arg1::SNES,arg2::Integer,arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:SNESFASGetRScale,petsc1),PetscErrorCode,(SNES,PetscInt,Ptr{Vec{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASSetContinuation(arg0::Type{Float64},arg1::SNES,arg2::PetscBool)
    ccall((:SNESFASSetContinuation,petsc1),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASGetSmoother(arg0::Type{Float64},arg1::SNES,arg2::Integer,arg3::Union(Ptr{SNES},StridedArray{SNES},Ptr{Void}))
    ccall((:SNESFASGetSmoother,petsc1),PetscErrorCode,(SNES,PetscInt,Ptr{SNES}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASGetSmootherUp(arg0::Type{Float64},arg1::SNES,arg2::Integer,arg3::Union(Ptr{SNES},StridedArray{SNES},Ptr{Void}))
    ccall((:SNESFASGetSmootherUp,petsc1),PetscErrorCode,(SNES,PetscInt,Ptr{SNES}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASGetSmootherDown(arg0::Type{Float64},arg1::SNES,arg2::Integer,arg3::Union(Ptr{SNES},StridedArray{SNES},Ptr{Void}))
    ccall((:SNESFASGetSmootherDown,petsc1),PetscErrorCode,(SNES,PetscInt,Ptr{SNES}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASGetCoarseSolve(arg0::Type{Float64},arg1::SNES,arg2::Union(Ptr{SNES},StridedArray{SNES},Ptr{Void}))
    ccall((:SNESFASGetCoarseSolve,petsc1),PetscErrorCode,(SNES,Ptr{SNES}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASFullSetDownSweep(arg0::Type{Float64},arg1::SNES,arg2::PetscBool)
    ccall((:SNESFASFullSetDownSweep,petsc1),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASCreateCoarseVec(arg1::SNES,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:SNESFASCreateCoarseVec,petsc1),PetscErrorCode,(SNES,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESFASRestrict(arg1::SNES,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:SNESFASRestrict,petsc1),PetscErrorCode,(SNES,Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMSNESCheckFromOptions(arg1::SNES,arg2::Vec{Float64},arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg4::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:DMSNESCheckFromOptions,petsc1),PetscErrorCode,(SNES,Vec{Float64},Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4)
end 
=#
function TSInitializePackage(arg0::Type{Float64})
    err = ccall((:TSInitializePackage,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function TSCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{TS},StridedArray{TS},Ptr{Void}))
    ccall((:TSCreate,petsc1),PetscErrorCode,(comm_type,Ptr{TS}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSClone(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{TS},StridedArray{TS},Ptr{Void}))
    ccall((:TSClone,petsc1),PetscErrorCode,(TS,Ptr{TS}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSDestroy(arg0::Type{Float64},arg1::Union(Ptr{TS},StridedArray{TS},Ptr{Void}))
    ccall((:TSDestroy,petsc1),PetscErrorCode,(Ptr{TS},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetProblemType(arg0::Type{Float64},arg1::TS,arg2::TSProblemType)
    ccall((:TSSetProblemType,petsc1),PetscErrorCode,(TS,TSProblemType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetProblemType(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{TSProblemType},StridedArray{TSProblemType},Ptr{Void}))
    ccall((:TSGetProblemType,petsc1),PetscErrorCode,(TS,Ptr{TSProblemType}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitor(arg1::TS,arg2::Integer,PetscReal::Integer,arg3::Vec{Float64})
    ccall((:TSMonitor,petsc1),PetscErrorCode,(TS,PetscInt,Cint,Vec{Float64}),arg1,arg2,PetscReal,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorSet(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSMonitorSet,petsc1),PetscErrorCode,(TS,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorCancel(arg0::Type{Float64},arg1::TS)
    ccall((:TSMonitorCancel,petsc1),PetscErrorCode,(TS,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetOptionsPrefix(arg0::Type{Float64},arg1::TS,arg2::Union(ByteString,Symbol))
    ccall((:TSSetOptionsPrefix,petsc1),PetscErrorCode,(TS,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSAppendOptionsPrefix(arg0::Type{Float64},arg1::TS,arg2::Union(ByteString,Symbol))
    ccall((:TSAppendOptionsPrefix,petsc1),PetscErrorCode,(TS,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetOptionsPrefix(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:TSGetOptionsPrefix,petsc1),PetscErrorCode,(TS,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetFromOptions(arg0::Type{Float64},arg1::TS)
    ccall((:TSSetFromOptions,petsc1),PetscErrorCode,(TS,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetUp(arg0::Type{Float64},arg1::TS)
    ccall((:TSSetUp,petsc1),PetscErrorCode,(TS,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSReset(arg0::Type{Float64},arg1::TS)
    ccall((:TSReset,petsc1),PetscErrorCode,(TS,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetSolution(arg1::TS,arg2::Vec{Float64})
    ccall((:TSSetSolution,petsc1),PetscErrorCode,(TS,Vec{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetSolution(arg1::TS,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:TSGetSolution,petsc1),PetscErrorCode,(TS,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetSaveTrajectory(arg0::Type{Float64},arg1::TS)
    ccall((:TSSetSaveTrajectory,petsc1),PetscErrorCode,(TS,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSTrajectoryCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{TSTrajectory},StridedArray{TSTrajectory},Ptr{Void}))
    ccall((:TSTrajectoryCreate,petsc1),PetscErrorCode,(comm_type,Ptr{TSTrajectory}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSTrajectoryDestroy(arg0::Type{Float64},arg1::Union(Ptr{TSTrajectory},StridedArray{TSTrajectory},Ptr{Void}))
    ccall((:TSTrajectoryDestroy,petsc1),PetscErrorCode,(Ptr{TSTrajectory},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSTrajectorySetType(arg0::Type{Float64},arg1::TSTrajectory,arg2::TSTrajectoryType)
    ccall((:TSTrajectorySetType,petsc1),PetscErrorCode,(TSTrajectory,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSTrajectorySet(arg1::TSTrajectory,arg2::TS,arg3::Integer,PetscReal::Integer,arg4::Vec{Float64})
    ccall((:TSTrajectorySet,petsc1),PetscErrorCode,(TSTrajectory,TS,PetscInt,Cint,Vec{Float64}),arg1,arg2,arg3,PetscReal,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSTrajectoryGet(arg0::Type{Float64},arg1::TSTrajectory,arg2::TS,arg3::Integer,arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TSTrajectoryGet,petsc1),PetscErrorCode,(TSTrajectory,TS,PetscInt,Ptr{Cint}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSTrajectorySetFromOptions(arg0::Type{Float64},arg1::TSTrajectory)
    ccall((:TSTrajectorySetFromOptions,petsc1),PetscErrorCode,(TSTrajectory,),arg1)
end 
=#
function TSTrajectoryRegisterAll(arg0::Type{Float64})
    err = ccall((:TSTrajectoryRegisterAll,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function TSSetCostGradients(arg1::TS,arg2::Integer,arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg4::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:TSSetCostGradients,petsc1),PetscErrorCode,(TS,PetscInt,Ptr{Vec{Float64}},Ptr{Vec{Float64}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetCostGradients(arg1::TS,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{Vec{Float64}}},StridedArray{Ptr{Vec{Float64}}},Ptr{Void}),arg4::Union(Ptr{Ptr{Vec{Float64}}},StridedArray{Ptr{Vec{Float64}}},Ptr{Void}))
    ccall((:TSGetCostGradients,petsc1),PetscErrorCode,(TS,Ptr{PetscInt},Ptr{Ptr{Vec{Float64}}},Ptr{Ptr{Vec{Float64}}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetCostIntegrand(arg0::Type{Float64},arg1::TS,arg2::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg6::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSSetCostIntegrand,petsc1),PetscErrorCode,(TS,PetscInt,Ptr{Void},Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetCostIntegral(arg1::TS,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:TSGetCostIntegral,petsc1),PetscErrorCode,(TS,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdjointSetRHSJacobian(arg1::TS,arg2::Mat{Float64},arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSAdjointSetRHSJacobian,petsc1),PetscErrorCode,(TS,Mat{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdjointSolve(arg0::Type{Float64},arg1::TS)
    ccall((:TSAdjointSolve,petsc1),PetscErrorCode,(TS,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdjointSetSteps(arg0::Type{Float64},arg1::TS,arg2::Integer)
    ccall((:TSAdjointSetSteps,petsc1),PetscErrorCode,(TS,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdjointComputeRHSJacobian(arg1::TS,PetscReal::Integer,arg2::Vec{Float64},arg3::Mat{Float64})
    ccall((:TSAdjointComputeRHSJacobian,petsc1),PetscErrorCode,(TS,Cint,Vec{Float64},Mat{Float64}),arg1,PetscReal,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdjointStep(arg0::Type{Float64},arg1::TS)
    ccall((:TSAdjointStep,petsc1),PetscErrorCode,(TS,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdjointSetUp(arg0::Type{Float64},arg1::TS)
    ccall((:TSAdjointSetUp,petsc1),PetscErrorCode,(TS,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdjointComputeDRDPFunction(arg1::TS,PetscReal::Integer,arg2::Vec{Float64},arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:TSAdjointComputeDRDPFunction,petsc1),PetscErrorCode,(TS,Cint,Vec{Float64},Ptr{Vec{Float64}}),arg1,PetscReal,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdjointComputeDRDYFunction(arg1::TS,PetscReal::Integer,arg2::Vec{Float64},arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:TSAdjointComputeDRDYFunction,petsc1),PetscErrorCode,(TS,Cint,Vec{Float64},Ptr{Vec{Float64}}),arg1,PetscReal,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdjointComputeCostIntegrand(arg1::TS,PetscReal::Integer,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:TSAdjointComputeCostIntegrand,petsc1),PetscErrorCode,(TS,Cint,Vec{Float64},Vec{Float64}),arg1,PetscReal,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetDuration(arg0::Type{Float64},arg1::TS,arg2::Integer,PetscReal::Integer)
    ccall((:TSSetDuration,petsc1),PetscErrorCode,(TS,PetscInt,Cint),arg1,arg2,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetDuration(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TSGetDuration,petsc1),PetscErrorCode,(TS,Ptr{PetscInt},Ptr{Cint}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetExactFinalTime(arg0::Type{Float64},arg1::TS,arg2::TSExactFinalTimeOption)
    ccall((:TSSetExactFinalTime,petsc1),PetscErrorCode,(TS,TSExactFinalTimeOption),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorDefault(arg1::TS,arg2::Integer,PetscReal::Integer,arg3::Vec{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSMonitorDefault,petsc1),PetscErrorCode,(TS,PetscInt,Cint,Vec{Float64},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorDrawCtxCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Integer,arg9::Union(Ptr{TSMonitorDrawCtx},StridedArray{TSMonitorDrawCtx},Ptr{Void}))
    ccall((:TSMonitorDrawCtxCreate,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint,Cint,Cint,Cint,PetscInt,Ptr{TSMonitorDrawCtx}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorDrawCtxDestroy(arg0::Type{Float64},arg1::Union(Ptr{TSMonitorDrawCtx},StridedArray{TSMonitorDrawCtx},Ptr{Void}))
    ccall((:TSMonitorDrawCtxDestroy,petsc1),PetscErrorCode,(Ptr{TSMonitorDrawCtx},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorDrawSolution(arg1::TS,arg2::Integer,PetscReal::Integer,arg3::Vec{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSMonitorDrawSolution,petsc1),PetscErrorCode,(TS,PetscInt,Cint,Vec{Float64},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorDrawSolutionPhase(arg1::TS,arg2::Integer,PetscReal::Integer,arg3::Vec{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSMonitorDrawSolutionPhase,petsc1),PetscErrorCode,(TS,PetscInt,Cint,Vec{Float64},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorDrawError(arg1::TS,arg2::Integer,PetscReal::Integer,arg3::Vec{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSMonitorDrawError,petsc1),PetscErrorCode,(TS,PetscInt,Cint,Vec{Float64},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorSolutionBinary(arg1::TS,arg2::Integer,PetscReal::Integer,arg3::Vec{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSMonitorSolutionBinary,petsc1),PetscErrorCode,(TS,PetscInt,Cint,Vec{Float64},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorSolutionVTK(arg1::TS,arg2::Integer,PetscReal::Integer,arg3::Vec{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSMonitorSolutionVTK,petsc1),PetscErrorCode,(TS,PetscInt,Cint,Vec{Float64},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end 
=#
function TSMonitorSolutionVTKDestroy(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:TSMonitorSolutionVTKDestroy,petsc1),PetscErrorCode,(Ptr{Void},),arg1)
    return err
end

#= skipping function with undefined symbols: 
 function TSStep(arg0::Type{Float64},arg1::TS)
    ccall((:TSStep,petsc1),PetscErrorCode,(TS,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSEvaluateStep(arg1::TS,arg2::Integer,arg3::Vec{Float64},arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:TSEvaluateStep,petsc1),PetscErrorCode,(TS,PetscInt,Vec{Float64},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSSolve(arg1::TS,arg2::Vec{Float64})
    ccall((:TSSolve,petsc1),PetscErrorCode,(TS,Vec{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetEquationType(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{TSEquationType},StridedArray{TSEquationType},Ptr{Void}))
    ccall((:TSGetEquationType,petsc1),PetscErrorCode,(TS,Ptr{TSEquationType}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetEquationType(arg0::Type{Float64},arg1::TS,arg2::TSEquationType)
    ccall((:TSSetEquationType,petsc1),PetscErrorCode,(TS,TSEquationType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetConvergedReason(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{TSConvergedReason},StridedArray{TSConvergedReason},Ptr{Void}))
    ccall((:TSGetConvergedReason,petsc1),PetscErrorCode,(TS,Ptr{TSConvergedReason}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetConvergedReason(arg0::Type{Float64},arg1::TS,arg2::TSConvergedReason)
    ccall((:TSSetConvergedReason,petsc1),PetscErrorCode,(TS,TSConvergedReason),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetSolveTime(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TSGetSolveTime,petsc1),PetscErrorCode,(TS,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetSNESIterations(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:TSGetSNESIterations,petsc1),PetscErrorCode,(TS,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetKSPIterations(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:TSGetKSPIterations,petsc1),PetscErrorCode,(TS,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetStepRejections(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:TSGetStepRejections,petsc1),PetscErrorCode,(TS,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetMaxStepRejections(arg0::Type{Float64},arg1::TS,arg2::Integer)
    ccall((:TSSetMaxStepRejections,petsc1),PetscErrorCode,(TS,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetSNESFailures(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:TSGetSNESFailures,petsc1),PetscErrorCode,(TS,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetMaxSNESFailures(arg0::Type{Float64},arg1::TS,arg2::Integer)
    ccall((:TSSetMaxSNESFailures,petsc1),PetscErrorCode,(TS,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetErrorIfStepFails(arg0::Type{Float64},arg1::TS,arg2::PetscBool)
    ccall((:TSSetErrorIfStepFails,petsc1),PetscErrorCode,(TS,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSRollBack(arg0::Type{Float64},arg1::TS)
    ccall((:TSRollBack,petsc1),PetscErrorCode,(TS,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetTotalSteps(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:TSGetTotalSteps,petsc1),PetscErrorCode,(TS,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetStages(arg1::TS,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{Vec{Float64}}},StridedArray{Ptr{Vec{Float64}}},Ptr{Void}))
    ccall((:TSGetStages,petsc1),PetscErrorCode,(TS,Ptr{PetscInt},Ptr{Ptr{Vec{Float64}}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetInitialTimeStep(arg0::Type{Float64},arg1::TS,PetscReal::Integer,arg2::Integer)
    ccall((:TSSetInitialTimeStep,petsc1),PetscErrorCode,(TS,Cint,Cint),arg1,PetscReal,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetTimeStep(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TSGetTimeStep,petsc1),PetscErrorCode,(TS,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetTime(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TSGetTime,petsc1),PetscErrorCode,(TS,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetTime(arg0::Type{Float64},arg1::TS,PetscReal::Integer)
    ccall((:TSSetTime,petsc1),PetscErrorCode,(TS,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetTimeStepNumber(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:TSGetTimeStepNumber,petsc1),PetscErrorCode,(TS,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetTimeStep(arg0::Type{Float64},arg1::TS,PetscReal::Integer)
    ccall((:TSSetTimeStep,petsc1),PetscErrorCode,(TS,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetPrevTime(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TSGetPrevTime,petsc1),PetscErrorCode,(TS,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetRHSFunction(arg1::TS,arg2::Vec{Float64},arg3::TSRHSFunction,arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSSetRHSFunction,petsc1),PetscErrorCode,(TS,Vec{Float64},TSRHSFunction,Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetRHSFunction(arg1::TS,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg3::Union(Ptr{TSRHSFunction},StridedArray{TSRHSFunction},Ptr{Void}),arg4::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:TSGetRHSFunction,petsc1),PetscErrorCode,(TS,Ptr{Vec{Float64}},Ptr{TSRHSFunction},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetRHSJacobian(arg1::TS,arg2::Mat{Float64},arg3::Mat{Float64},arg4::TSRHSJacobian,arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSSetRHSJacobian,petsc1),PetscErrorCode,(TS,Mat{Float64},Mat{Float64},TSRHSJacobian,Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetRHSJacobian(arg1::TS,arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg4::Union(Ptr{TSRHSJacobian},StridedArray{TSRHSJacobian},Ptr{Void}),arg5::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:TSGetRHSJacobian,petsc1),PetscErrorCode,(TS,Ptr{Mat{Float64}},Ptr{Mat{Float64}},Ptr{TSRHSJacobian},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function TSRHSJacobianSetReuse(arg0::Type{Float64},arg1::TS,arg2::PetscBool)
    ccall((:TSRHSJacobianSetReuse,petsc1),PetscErrorCode,(TS,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetSolutionFunction(arg0::Type{Float64},arg1::TS,arg2::TSSolutionFunction,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSSetSolutionFunction,petsc1),PetscErrorCode,(TS,TSSolutionFunction,Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetForcingFunction(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSSetForcingFunction,petsc1),PetscErrorCode,(TS,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetIFunction(arg1::TS,arg2::Vec{Float64},arg3::TSIFunction,arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSSetIFunction,petsc1),PetscErrorCode,(TS,Vec{Float64},TSIFunction,Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetIFunction(arg1::TS,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg3::Union(Ptr{TSIFunction},StridedArray{TSIFunction},Ptr{Void}),arg4::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:TSGetIFunction,petsc1),PetscErrorCode,(TS,Ptr{Vec{Float64}},Ptr{TSIFunction},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetIJacobian(arg1::TS,arg2::Mat{Float64},arg3::Mat{Float64},arg4::TSIJacobian,arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSSetIJacobian,petsc1),PetscErrorCode,(TS,Mat{Float64},Mat{Float64},TSIJacobian,Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetIJacobian(arg1::TS,arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg4::Union(Ptr{TSIJacobian},StridedArray{TSIJacobian},Ptr{Void}),arg5::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:TSGetIJacobian,petsc1),PetscErrorCode,(TS,Ptr{Mat{Float64}},Ptr{Mat{Float64}},Ptr{TSIJacobian},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function TSComputeRHSFunctionLinear(arg1::TS,PetscReal::Integer,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSComputeRHSFunctionLinear,petsc1),PetscErrorCode,(TS,Cint,Vec{Float64},Vec{Float64},Ptr{Void}),arg1,PetscReal,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSComputeRHSJacobianConstant(arg1::TS,PetscReal::Integer,arg2::Vec{Float64},arg3::Mat{Float64},arg4::Mat{Float64},arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSComputeRHSJacobianConstant,petsc1),PetscErrorCode,(TS,Cint,Vec{Float64},Mat{Float64},Mat{Float64},Ptr{Void}),arg1,PetscReal,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function TSComputeIFunctionLinear(arg1::TS,PetscReal::Integer,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64},arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSComputeIFunctionLinear,petsc1),PetscErrorCode,(TS,Cint,Vec{Float64},Vec{Float64},Vec{Float64},Ptr{Void}),arg1,PetscReal,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function TSComputeIJacobianConstant(arg1::TS,PetscReal::Integer,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Integer,arg5::Mat{Float64},arg6::Mat{Float64},arg7::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSComputeIJacobianConstant,petsc1),PetscErrorCode,(TS,Cint,Vec{Float64},Vec{Float64},Cint,Mat{Float64},Mat{Float64},Ptr{Void}),arg1,PetscReal,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function TSComputeSolutionFunction(arg1::TS,PetscReal::Integer,arg2::Vec{Float64})
    ccall((:TSComputeSolutionFunction,petsc1),PetscErrorCode,(TS,Cint,Vec{Float64}),arg1,PetscReal,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSComputeForcingFunction(arg1::TS,PetscReal::Integer,arg2::Vec{Float64})
    ccall((:TSComputeForcingFunction,petsc1),PetscErrorCode,(TS,Cint,Vec{Float64}),arg1,PetscReal,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSComputeIJacobianDefaultColor(arg1::TS,PetscReal::Integer,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Integer,arg5::Mat{Float64},arg6::Mat{Float64},arg7::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSComputeIJacobianDefaultColor,petsc1),PetscErrorCode,(TS,Cint,Vec{Float64},Vec{Float64},Cint,Mat{Float64},Mat{Float64},Ptr{Void}),arg1,PetscReal,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetPreStep(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSSetPreStep,petsc1),PetscErrorCode,(TS,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetPreStage(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSSetPreStage,petsc1),PetscErrorCode,(TS,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetPostStage(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSSetPostStage,petsc1),PetscErrorCode,(TS,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetPostStep(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSSetPostStep,petsc1),PetscErrorCode,(TS,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSPreStep(arg0::Type{Float64},arg1::TS)
    ccall((:TSPreStep,petsc1),PetscErrorCode,(TS,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSPreStage(arg0::Type{Float64},arg1::TS,PetscReal::Integer)
    ccall((:TSPreStage,petsc1),PetscErrorCode,(TS,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function TSPostStage(arg1::TS,PetscReal::Integer,arg2::Integer,arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:TSPostStage,petsc1),PetscErrorCode,(TS,Cint,PetscInt,Ptr{Vec{Float64}}),arg1,PetscReal,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TSPostStep(arg0::Type{Float64},arg1::TS)
    ccall((:TSPostStep,petsc1),PetscErrorCode,(TS,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetRetainStages(arg0::Type{Float64},arg1::TS,arg2::PetscBool)
    ccall((:TSSetRetainStages,petsc1),PetscErrorCode,(TS,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSInterpolate(arg1::TS,PetscReal::Integer,arg2::Vec{Float64})
    ccall((:TSInterpolate,petsc1),PetscErrorCode,(TS,Cint,Vec{Float64}),arg1,PetscReal,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetTolerances(arg1::TS,PetscReal::Integer,arg2::Vec{Float64},arg3::Integer,arg4::Vec{Float64})
    ccall((:TSSetTolerances,petsc1),PetscErrorCode,(TS,Cint,Vec{Float64},Cint,Vec{Float64}),arg1,PetscReal,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetTolerances(arg1::TS,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:TSGetTolerances,petsc1),PetscErrorCode,(TS,Ptr{Cint},Ptr{Vec{Float64}},Ptr{Cint},Ptr{Vec{Float64}}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function TSErrorWeightedNormInfinity(arg1::TS,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TSErrorWeightedNormInfinity,petsc1),PetscErrorCode,(TS,Vec{Float64},Vec{Float64},Ptr{Cint}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSErrorWeightedNorm2(arg1::TS,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TSErrorWeightedNorm2,petsc1),PetscErrorCode,(TS,Vec{Float64},Vec{Float64},Ptr{Cint}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSErrorWeightedNorm(arg1::TS,arg2::Vec{Float64},arg3::Vec{Float64},arg4::NormType,arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TSErrorWeightedNorm,petsc1),PetscErrorCode,(TS,Vec{Float64},Vec{Float64},NormType,Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetCFLTimeLocal(arg0::Type{Float64},arg1::TS,PetscReal::Integer)
    ccall((:TSSetCFLTimeLocal,petsc1),PetscErrorCode,(TS,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetCFLTime(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TSGetCFLTime,petsc1),PetscErrorCode,(TS,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSPseudoSetTimeStep(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSPseudoSetTimeStep,petsc1),PetscErrorCode,(TS,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TSPseudoTimeStepDefault(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSPseudoTimeStepDefault,petsc1),PetscErrorCode,(TS,Ptr{Cint},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TSPseudoComputeTimeStep(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TSPseudoComputeTimeStep,petsc1),PetscErrorCode,(TS,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSPseudoSetMaxTimeStep(arg0::Type{Float64},arg1::TS,PetscReal::Integer)
    ccall((:TSPseudoSetMaxTimeStep,petsc1),PetscErrorCode,(TS,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function TSPseudoSetVerifyTimeStep(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSPseudoSetVerifyTimeStep,petsc1),PetscErrorCode,(TS,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TSPseudoVerifyTimeStepDefault(arg1::TS,arg2::Vec{Float64},arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:TSPseudoVerifyTimeStepDefault,petsc1),PetscErrorCode,(TS,Vec{Float64},Ptr{Void},Ptr{Cint},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function TSPseudoVerifyTimeStep(arg1::TS,arg2::Vec{Float64},arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:TSPseudoVerifyTimeStep,petsc1),PetscErrorCode,(TS,Vec{Float64},Ptr{Cint},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSPseudoSetTimeStepIncrement(arg0::Type{Float64},arg1::TS,PetscReal::Integer)
    ccall((:TSPseudoSetTimeStepIncrement,petsc1),PetscErrorCode,(TS,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function TSPseudoIncrementDtFromInitialDt(arg0::Type{Float64},arg1::TS)
    ccall((:TSPseudoIncrementDtFromInitialDt,petsc1),PetscErrorCode,(TS,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSPythonSetType(arg0::Type{Float64},arg1::TS,arg2::Union(ByteString,Symbol))
    ccall((:TSPythonSetType,petsc1),PetscErrorCode,(TS,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSComputeRHSFunction(arg1::TS,PetscReal::Integer,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:TSComputeRHSFunction,petsc1),PetscErrorCode,(TS,Cint,Vec{Float64},Vec{Float64}),arg1,PetscReal,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TSComputeRHSJacobian(arg1::TS,PetscReal::Integer,arg2::Vec{Float64},arg3::Mat{Float64},arg4::Mat{Float64})
    ccall((:TSComputeRHSJacobian,petsc1),PetscErrorCode,(TS,Cint,Vec{Float64},Mat{Float64},Mat{Float64}),arg1,PetscReal,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSComputeIFunction(arg1::TS,PetscReal::Integer,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64},arg5::PetscBool)
    ccall((:TSComputeIFunction,petsc1),PetscErrorCode,(TS,Cint,Vec{Float64},Vec{Float64},Vec{Float64},PetscBool),arg1,PetscReal,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function TSComputeIJacobian(arg1::TS,PetscReal::Integer,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Integer,arg5::Mat{Float64},arg6::Mat{Float64},arg7::PetscBool)
    ccall((:TSComputeIJacobian,petsc1),PetscErrorCode,(TS,Cint,Vec{Float64},Vec{Float64},Cint,Mat{Float64},Mat{Float64},PetscBool),arg1,PetscReal,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function TSComputeLinearStability(arg0::Type{Float64},arg1::TS,PetscReal::Integer,arg2::Integer,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TSComputeLinearStability,petsc1),PetscErrorCode,(TS,Cint,Cint,Ptr{Cint},Ptr{Cint}),arg1,PetscReal,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSVISetVariableBounds(arg1::TS,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:TSVISetVariableBounds,petsc1),PetscErrorCode,(TS,Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMTSSetRHSFunction(arg0::Type{Float64},arg1::DM,arg2::TSRHSFunction,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMTSSetRHSFunction,petsc1),PetscErrorCode,(DM,TSRHSFunction,Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMTSGetRHSFunction(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{TSRHSFunction},StridedArray{TSRHSFunction},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:DMTSGetRHSFunction,petsc1),PetscErrorCode,(DM,Ptr{TSRHSFunction},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMTSSetRHSJacobian(arg0::Type{Float64},arg1::DM,arg2::TSRHSJacobian,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMTSSetRHSJacobian,petsc1),PetscErrorCode,(DM,TSRHSJacobian,Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMTSGetRHSJacobian(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{TSRHSJacobian},StridedArray{TSRHSJacobian},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:DMTSGetRHSJacobian,petsc1),PetscErrorCode,(DM,Ptr{TSRHSJacobian},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMTSSetIFunction(arg0::Type{Float64},arg1::DM,arg2::TSIFunction,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMTSSetIFunction,petsc1),PetscErrorCode,(DM,TSIFunction,Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMTSGetIFunction(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{TSIFunction},StridedArray{TSIFunction},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:DMTSGetIFunction,petsc1),PetscErrorCode,(DM,Ptr{TSIFunction},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMTSSetIJacobian(arg0::Type{Float64},arg1::DM,arg2::TSIJacobian,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMTSSetIJacobian,petsc1),PetscErrorCode,(DM,TSIJacobian,Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMTSGetIJacobian(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{TSIJacobian},StridedArray{TSIJacobian},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:DMTSGetIJacobian,petsc1),PetscErrorCode,(DM,Ptr{TSIJacobian},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMTSSetSolutionFunction(arg0::Type{Float64},arg1::DM,arg2::TSSolutionFunction,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMTSSetSolutionFunction,petsc1),PetscErrorCode,(DM,TSSolutionFunction,Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMTSGetSolutionFunction(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{TSSolutionFunction},StridedArray{TSSolutionFunction},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:DMTSGetSolutionFunction,petsc1),PetscErrorCode,(DM,Ptr{TSSolutionFunction},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMTSSetForcingFunction(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMTSSetForcingFunction,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMTSGetForcingFunction(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:DMTSGetForcingFunction,petsc1),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMTSGetMinRadius(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:DMTSGetMinRadius,petsc1),PetscErrorCode,(DM,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function DMTSSetMinRadius(arg0::Type{Float64},arg1::DM,PetscReal::Integer)
    ccall((:DMTSSetMinRadius,petsc1),PetscErrorCode,(DM,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function DMTSCheckFromOptions(arg1::TS,arg2::Vec{Float64},arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}),arg4::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:DMTSCheckFromOptions,petsc1),PetscErrorCode,(TS,Vec{Float64},Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMTSSetIFunctionLocal(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMTSSetIFunctionLocal,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMTSSetIJacobianLocal(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMTSSetIJacobianLocal,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMTSSetRHSFunctionLocal(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMTSSetRHSFunctionLocal,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMTSSetIFunctionSerialize(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMTSSetIFunctionSerialize,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMTSSetIJacobianSerialize(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMTSSetIJacobianSerialize,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDATSSetRHSFunctionLocal(arg0::Type{Float64},arg1::DM,arg2::InsertMode,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMDATSSetRHSFunctionLocal,petsc1),PetscErrorCode,(DM,InsertMode,Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMDATSSetRHSJacobianLocal(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMDATSSetRHSJacobianLocal,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMDATSSetIFunctionLocal(arg0::Type{Float64},arg1::DM,arg2::InsertMode,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMDATSSetIFunctionLocal,petsc1),PetscErrorCode,(DM,InsertMode,Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function DMDATSSetIJacobianLocal(arg0::Type{Float64},arg1::DM,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:DMDATSSetIJacobianLocal,petsc1),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function DMPlexTSGetGeometryFVM(arg1::DM,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:DMPlexTSGetGeometryFVM,petsc1),PetscErrorCode,(DM,Ptr{Vec{Float64}},Ptr{Vec{Float64}},Ptr{Cint}),arg1,arg2,arg3,arg4)
end 
=#
function TSMonitorDMDARayDestroy(arg0::Type{Float64},arg1::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    err = ccall((:TSMonitorDMDARayDestroy,petsc1),PetscErrorCode,(Ptr{Ptr{Void}},),arg1)
    return err
end

#= skipping function with undefined symbols: 
 function TSMonitorDMDARay(arg1::TS,arg2::Integer,PetscReal::Integer,arg3::Vec{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSMonitorDMDARay,petsc1),PetscErrorCode,(TS,PetscInt,Cint,Vec{Float64},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorLGDMDARay(arg1::TS,arg2::Integer,PetscReal::Integer,arg3::Vec{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSMonitorLGDMDARay,petsc1),PetscErrorCode,(TS,PetscInt,Cint,Vec{Float64},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetType(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{TSType},StridedArray{TSType},Ptr{Void}))
    ccall((:TSGetType,petsc1),PetscErrorCode,(TS,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetType(arg0::Type{Float64},arg1::TS,arg2::TSType)
    ccall((:TSSetType,petsc1),PetscErrorCode,(TS,Cstring),arg1,arg2)
end 
=#
function TSRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:TSRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function TSGetSNES(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{SNES},StridedArray{SNES},Ptr{Void}))
    ccall((:TSGetSNES,petsc1),PetscErrorCode,(TS,Ptr{SNES}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetSNES(arg0::Type{Float64},arg1::TS,arg2::SNES)
    ccall((:TSSetSNES,petsc1),PetscErrorCode,(TS,SNES),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetKSP(arg1::TS,arg2::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    ccall((:TSGetKSP,petsc1),PetscErrorCode,(TS,Ptr{KSP{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSView(arg1::TS,arg2::PetscViewer{Float64})
    ccall((:TSView,petsc1),PetscErrorCode,(TS,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSLoad(arg1::TS,arg2::PetscViewer{Float64})
    ccall((:TSLoad,petsc1),PetscErrorCode,(TS,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetApplicationContext(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSGetApplicationContext,petsc1),PetscErrorCode,(TS,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorLGCtxCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Integer,arg9::Union(Ptr{TSMonitorLGCtx},StridedArray{TSMonitorLGCtx},Ptr{Void}))
    ccall((:TSMonitorLGCtxCreate,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint,Cint,Cint,Cint,PetscInt,Ptr{TSMonitorLGCtx}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorLGCtxDestroy(arg0::Type{Float64},arg1::Union(Ptr{TSMonitorLGCtx},StridedArray{TSMonitorLGCtx},Ptr{Void}))
    ccall((:TSMonitorLGCtxDestroy,petsc1),PetscErrorCode,(Ptr{TSMonitorLGCtx},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorLGTimeStep(arg1::TS,arg2::Integer,PetscReal::Integer,arg3::Vec{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSMonitorLGTimeStep,petsc1),PetscErrorCode,(TS,PetscInt,Cint,Vec{Float64},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorLGSolution(arg1::TS,arg2::Integer,PetscReal::Integer,arg3::Vec{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSMonitorLGSolution,petsc1),PetscErrorCode,(TS,PetscInt,Cint,Vec{Float64},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorLGSetVariableNames(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:TSMonitorLGSetVariableNames,petsc1),PetscErrorCode,(TS,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorLGGetVariableNames(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Ptr{Ptr{Uint8}}},StridedArray{Ptr{Ptr{Uint8}}},Ptr{Void}))
    ccall((:TSMonitorLGGetVariableNames,petsc1),PetscErrorCode,(TS,Ptr{Ptr{Ptr{Uint8}}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorLGCtxSetVariableNames(arg0::Type{Float64},arg1::TSMonitorLGCtx,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:TSMonitorLGCtxSetVariableNames,petsc1),PetscErrorCode,(TSMonitorLGCtx,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorLGSetDisplayVariables(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:TSMonitorLGSetDisplayVariables,petsc1),PetscErrorCode,(TS,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorLGCtxSetDisplayVariables(arg0::Type{Float64},arg1::TSMonitorLGCtx,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:TSMonitorLGCtxSetDisplayVariables,petsc1),PetscErrorCode,(TSMonitorLGCtx,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorLGSetTransform(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSMonitorLGSetTransform,petsc1),PetscErrorCode,(TS,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorLGCtxSetTransform(arg0::Type{Float64},arg1::TSMonitorLGCtx,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSMonitorLGCtxSetTransform,petsc1),PetscErrorCode,(TSMonitorLGCtx,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorLGError(arg1::TS,arg2::Integer,PetscReal::Integer,arg3::Vec{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSMonitorLGError,petsc1),PetscErrorCode,(TS,PetscInt,Cint,Vec{Float64},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorLGSNESIterations(arg1::TS,arg2::Integer,PetscReal::Integer,arg3::Vec{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSMonitorLGSNESIterations,petsc1),PetscErrorCode,(TS,PetscInt,Cint,Vec{Float64},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorLGKSPIterations(arg1::TS,arg2::Integer,PetscReal::Integer,arg3::Vec{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSMonitorLGKSPIterations,petsc1),PetscErrorCode,(TS,PetscInt,Cint,Vec{Float64},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorEnvelopeCtxCreate(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{TSMonitorEnvelopeCtx},StridedArray{TSMonitorEnvelopeCtx},Ptr{Void}))
    ccall((:TSMonitorEnvelopeCtxCreate,petsc1),PetscErrorCode,(TS,Ptr{TSMonitorEnvelopeCtx}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorEnvelope(arg1::TS,arg2::Integer,PetscReal::Integer,arg3::Vec{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSMonitorEnvelope,petsc1),PetscErrorCode,(TS,PetscInt,Cint,Vec{Float64},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorEnvelopeGetBounds(arg1::TS,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:TSMonitorEnvelopeGetBounds,petsc1),PetscErrorCode,(TS,Ptr{Vec{Float64}},Ptr{Vec{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorEnvelopeCtxDestroy(arg0::Type{Float64},arg1::Union(Ptr{TSMonitorEnvelopeCtx},StridedArray{TSMonitorEnvelopeCtx},Ptr{Void}))
    ccall((:TSMonitorEnvelopeCtxDestroy,petsc1),PetscErrorCode,(Ptr{TSMonitorEnvelopeCtx},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorSPEigCtxCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(ByteString,Symbol),arg3::Union(ByteString,Symbol),arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Integer,arg9::Union(Ptr{TSMonitorSPEigCtx},StridedArray{TSMonitorSPEigCtx},Ptr{Void}))
    ccall((:TSMonitorSPEigCtxCreate,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint,Cint,Cint,Cint,PetscInt,Ptr{TSMonitorSPEigCtx}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorSPEigCtxDestroy(arg0::Type{Float64},arg1::Union(Ptr{TSMonitorSPEigCtx},StridedArray{TSMonitorSPEigCtx},Ptr{Void}))
    ccall((:TSMonitorSPEigCtxDestroy,petsc1),PetscErrorCode,(Ptr{TSMonitorSPEigCtx},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSMonitorSPEig(arg1::TS,arg2::Integer,PetscReal::Integer,arg3::Vec{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSMonitorSPEig,petsc1),PetscErrorCode,(TS,PetscInt,Cint,Vec{Float64},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetEventMonitor(arg0::Type{Float64},arg1::TS,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg6::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg7::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSSetEventMonitor,petsc1),PetscErrorCode,(TS,PetscInt,Ptr{PetscInt},Ptr{PetscBool},Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetEventTolerances(arg0::Type{Float64},arg1::TS,PetscReal::Integer,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TSSetEventTolerances,petsc1),PetscErrorCode,(TS,Cint,Ptr{Cint}),arg1,PetscReal,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSSSPSetType(arg0::Type{Float64},arg1::TS,arg2::TSSSPType)
    ccall((:TSSSPSetType,petsc1),PetscErrorCode,(TS,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSSSPGetType(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{TSSSPType},StridedArray{TSSSPType},Ptr{Void}))
    ccall((:TSSSPGetType,petsc1),PetscErrorCode,(TS,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSSSPSetNumStages(arg0::Type{Float64},arg1::TS,arg2::Integer)
    ccall((:TSSSPSetNumStages,petsc1),PetscErrorCode,(TS,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSSSPGetNumStages(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:TSSSPGetNumStages,petsc1),PetscErrorCode,(TS,Ptr{PetscInt}),arg1,arg2)
end 
=#
function TSSSPFinalizePackage(arg0::Type{Float64})
    err = ccall((:TSSSPFinalizePackage,petsc1),PetscErrorCode,())
    return err
end

function TSSSPInitializePackage(arg0::Type{Float64})
    err = ccall((:TSSSPInitializePackage,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function TSGetAdapt(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{TSAdapt},StridedArray{TSAdapt},Ptr{Void}))
    ccall((:TSGetAdapt,petsc1),PetscErrorCode,(TS,Ptr{TSAdapt}),arg1,arg2)
end 
=#
function TSAdaptRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:TSAdaptRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

function TSAdaptInitializePackage(arg0::Type{Float64})
    err = ccall((:TSAdaptInitializePackage,petsc1),PetscErrorCode,())
    return err
end

function TSAdaptFinalizePackage(arg0::Type{Float64})
    err = ccall((:TSAdaptFinalizePackage,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function TSAdaptCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{TSAdapt},StridedArray{TSAdapt},Ptr{Void}))
    ccall((:TSAdaptCreate,petsc1),PetscErrorCode,(comm_type,Ptr{TSAdapt}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdaptSetType(arg0::Type{Float64},arg1::TSAdapt,arg2::TSAdaptType)
    ccall((:TSAdaptSetType,petsc1),PetscErrorCode,(TSAdapt,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdaptSetOptionsPrefix(arg0::Type{Float64},arg1::TSAdapt,arg2::Union(ByteString,Symbol))
    ccall((:TSAdaptSetOptionsPrefix,petsc1),PetscErrorCode,(TSAdapt,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdaptCandidatesClear(arg0::Type{Float64},arg1::TSAdapt)
    ccall((:TSAdaptCandidatesClear,petsc1),PetscErrorCode,(TSAdapt,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdaptCandidateAdd(arg0::Type{Float64},arg1::TSAdapt,arg2::Union(ByteString,Symbol),arg3::Integer,arg4::Integer,PetscReal::Integer,arg5::Integer,arg6::PetscBool)
    ccall((:TSAdaptCandidateAdd,petsc1),PetscErrorCode,(TSAdapt,Cstring,PetscInt,PetscInt,Cint,Cint,PetscBool),arg1,arg2,arg3,arg4,PetscReal,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdaptCandidatesGet(arg0::Type{Float64},arg1::TSAdapt,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg4::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg5::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg6::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}))
    ccall((:TSAdaptCandidatesGet,petsc1),PetscErrorCode,(TSAdapt,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdaptChoose(arg0::Type{Float64},arg1::TSAdapt,arg2::TS,PetscReal::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:TSAdaptChoose,petsc1),PetscErrorCode,(TSAdapt,TS,Cint,Ptr{PetscInt},Ptr{Cint},Ptr{PetscBool}),arg1,arg2,PetscReal,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdaptCheckStage(arg0::Type{Float64},arg1::TSAdapt,arg2::TS,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:TSAdaptCheckStage,petsc1),PetscErrorCode,(TSAdapt,TS,Ptr{PetscBool}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdaptView(arg1::TSAdapt,arg2::PetscViewer{Float64})
    ccall((:TSAdaptView,petsc1),PetscErrorCode,(TSAdapt,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdaptLoad(arg1::TSAdapt,arg2::PetscViewer{Float64})
    ccall((:TSAdaptLoad,petsc1),PetscErrorCode,(TSAdapt,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdaptSetFromOptions(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::TSAdapt)
    ccall((:TSAdaptSetFromOptions,petsc1),PetscErrorCode,(Ptr{PetscOptions},TSAdapt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdaptReset(arg0::Type{Float64},arg1::TSAdapt)
    ccall((:TSAdaptReset,petsc1),PetscErrorCode,(TSAdapt,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdaptDestroy(arg0::Type{Float64},arg1::Union(Ptr{TSAdapt},StridedArray{TSAdapt},Ptr{Void}))
    ccall((:TSAdaptDestroy,petsc1),PetscErrorCode,(Ptr{TSAdapt},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdaptSetMonitor(arg0::Type{Float64},arg1::TSAdapt,arg2::PetscBool)
    ccall((:TSAdaptSetMonitor,petsc1),PetscErrorCode,(TSAdapt,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdaptSetStepLimits(arg0::Type{Float64},arg1::TSAdapt,PetscReal::Integer,arg2::Integer)
    ccall((:TSAdaptSetStepLimits,petsc1),PetscErrorCode,(TSAdapt,Cint,Cint),arg1,PetscReal,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSAdaptSetCheckStage(arg0::Type{Float64},arg1::TSAdapt,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSAdaptSetCheckStage,petsc1),PetscErrorCode,(TSAdapt,Ptr{Void}),arg1,arg2)
end 
=#
function TSGLAdaptRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:TSGLAdaptRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

function TSGLAdaptInitializePackage(arg0::Type{Float64})
    err = ccall((:TSGLAdaptInitializePackage,petsc1),PetscErrorCode,())
    return err
end

function TSGLAdaptFinalizePackage(arg0::Type{Float64})
    err = ccall((:TSGLAdaptFinalizePackage,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function TSGLAdaptCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{TSGLAdapt},StridedArray{TSGLAdapt},Ptr{Void}))
    ccall((:TSGLAdaptCreate,petsc1),PetscErrorCode,(comm_type,Ptr{TSGLAdapt}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSGLAdaptSetType(arg0::Type{Float64},arg1::TSGLAdapt,arg2::TSGLAdaptType)
    ccall((:TSGLAdaptSetType,petsc1),PetscErrorCode,(TSGLAdapt,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSGLAdaptSetOptionsPrefix(arg0::Type{Float64},arg1::TSGLAdapt,arg2::Union(ByteString,Symbol))
    ccall((:TSGLAdaptSetOptionsPrefix,petsc1),PetscErrorCode,(TSGLAdapt,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSGLAdaptChoose(arg0::Type{Float64},arg1::TSGLAdapt,arg2::Integer,arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Integer,arg6::Integer,arg7::Integer,arg8::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg9::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg10::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:TSGLAdaptChoose,petsc1),PetscErrorCode,(TSGLAdapt,PetscInt,Ptr{PetscInt},Ptr{Cint},Ptr{Cint},PetscInt,Cint,Cint,Ptr{PetscInt},Ptr{Cint},Ptr{PetscBool}),arg1,arg2,arg3,PetscReal,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end 
=#
#= skipping function with undefined symbols: 
 function TSGLAdaptView(arg1::TSGLAdapt,arg2::PetscViewer{Float64})
    ccall((:TSGLAdaptView,petsc1),PetscErrorCode,(TSGLAdapt,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSGLAdaptSetFromOptions(arg0::Type{Float64},arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::TSGLAdapt)
    ccall((:TSGLAdaptSetFromOptions,petsc1),PetscErrorCode,(Ptr{PetscOptions},TSGLAdapt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSGLAdaptDestroy(arg0::Type{Float64},arg1::Union(Ptr{TSGLAdapt},StridedArray{TSGLAdapt},Ptr{Void}))
    ccall((:TSGLAdaptDestroy,petsc1),PetscErrorCode,(Ptr{TSGLAdapt},),arg1)
end 
=#
function TSGLAcceptRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::TSGLAcceptFunction)
    err = ccall((:TSGLAcceptRegister,petsc1),PetscErrorCode,(Cstring,TSGLAcceptFunction),arg1,arg2)
    return err
end

function TSGLRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:TSGLRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

function TSGLInitializePackage(arg0::Type{Float64})
    err = ccall((:TSGLInitializePackage,petsc1),PetscErrorCode,())
    return err
end

function TSGLFinalizePackage(arg0::Type{Float64})
    err = ccall((:TSGLFinalizePackage,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function TSGLSetType(arg0::Type{Float64},arg1::TS,arg2::TSGLType)
    ccall((:TSGLSetType,petsc1),PetscErrorCode,(TS,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSGLGetAdapt(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{TSGLAdapt},StridedArray{TSGLAdapt},Ptr{Void}))
    ccall((:TSGLGetAdapt,petsc1),PetscErrorCode,(TS,Ptr{TSGLAdapt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSGLSetAcceptType(arg0::Type{Float64},arg1::TS,arg2::TSGLAcceptType)
    ccall((:TSGLSetAcceptType,petsc1),PetscErrorCode,(TS,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSEIMEXSetMaxRows(arg0::Type{Float64},ts::TS,arg1::Integer)
    ccall((:TSEIMEXSetMaxRows,petsc1),PetscErrorCode,(TS,PetscInt),ts,arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSEIMEXSetRowCol(arg0::Type{Float64},ts::TS,arg1::Integer,arg2::Integer)
    ccall((:TSEIMEXSetRowCol,petsc1),PetscErrorCode,(TS,PetscInt,PetscInt),ts,arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSEIMEXSetOrdAdapt(arg0::Type{Float64},arg1::TS,arg2::PetscBool)
    ccall((:TSEIMEXSetOrdAdapt,petsc1),PetscErrorCode,(TS,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSRKGetType(arg0::Type{Float64},ts::TS,arg1::Union(Ptr{TSRKType},StridedArray{TSRKType},Ptr{Void}))
    ccall((:TSRKGetType,petsc1),PetscErrorCode,(TS,Ptr{Ptr{Uint8}}),ts,arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSRKSetType(arg0::Type{Float64},ts::TS,arg1::TSRKType)
    ccall((:TSRKSetType,petsc1),PetscErrorCode,(TS,Cstring),ts,arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSRKSetFullyImplicit(arg0::Type{Float64},arg1::TS,arg2::PetscBool)
    ccall((:TSRKSetFullyImplicit,petsc1),PetscErrorCode,(TS,PetscBool),arg1,arg2)
end 
=#
function TSRKRegister(arg0::Type{Float64},arg1::TSRKType,arg2::Integer,arg3::Integer,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg7::Integer,arg8::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:TSRKRegister,petsc1),PetscErrorCode,(Cstring,PetscInt,PetscInt,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},PetscInt,Ptr{Cint}),arg1,arg2,arg3,PetscReal,arg4,arg5,arg6,arg7,arg8)
    return err
end

function TSRKFinalizePackage(arg0::Type{Float64})
    err = ccall((:TSRKFinalizePackage,petsc1),PetscErrorCode,())
    return err
end

function TSRKInitializePackage(arg0::Type{Float64})
    err = ccall((:TSRKInitializePackage,petsc1),PetscErrorCode,())
    return err
end

function TSRKRegisterDestroy(arg0::Type{Float64})
    err = ccall((:TSRKRegisterDestroy,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function TSARKIMEXGetType(arg0::Type{Float64},ts::TS,arg1::Union(Ptr{TSARKIMEXType},StridedArray{TSARKIMEXType},Ptr{Void}))
    ccall((:TSARKIMEXGetType,petsc1),PetscErrorCode,(TS,Ptr{Ptr{Uint8}}),ts,arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSARKIMEXSetType(arg0::Type{Float64},ts::TS,arg1::TSARKIMEXType)
    ccall((:TSARKIMEXSetType,petsc1),PetscErrorCode,(TS,Cstring),ts,arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSARKIMEXSetFullyImplicit(arg0::Type{Float64},arg1::TS,arg2::PetscBool)
    ccall((:TSARKIMEXSetFullyImplicit,petsc1),PetscErrorCode,(TS,PetscBool),arg1,arg2)
end 
=#
function TSARKIMEXRegister(arg0::Type{Float64},arg1::TSARKIMEXType,arg2::Integer,arg3::Integer,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg7::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg8::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg9::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg10::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg11::Integer,arg12::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg13::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:TSARKIMEXRegister,petsc1),PetscErrorCode,(Cstring,PetscInt,PetscInt,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},PetscInt,Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,PetscReal,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13)
    return err
end

function TSARKIMEXFinalizePackage(arg0::Type{Float64})
    err = ccall((:TSARKIMEXFinalizePackage,petsc1),PetscErrorCode,())
    return err
end

function TSARKIMEXInitializePackage(arg0::Type{Float64})
    err = ccall((:TSARKIMEXInitializePackage,petsc1),PetscErrorCode,())
    return err
end

function TSARKIMEXRegisterDestroy(arg0::Type{Float64})
    err = ccall((:TSARKIMEXRegisterDestroy,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function TSRosWGetType(arg0::Type{Float64},ts::TS,arg1::Union(Ptr{TSRosWType},StridedArray{TSRosWType},Ptr{Void}))
    ccall((:TSRosWGetType,petsc1),PetscErrorCode,(TS,Ptr{Ptr{Uint8}}),ts,arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSRosWSetType(arg0::Type{Float64},ts::TS,arg1::TSRosWType)
    ccall((:TSRosWSetType,petsc1),PetscErrorCode,(TS,Cstring),ts,arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TSRosWSetRecomputeJacobian(arg0::Type{Float64},arg1::TS,arg2::PetscBool)
    ccall((:TSRosWSetRecomputeJacobian,petsc1),PetscErrorCode,(TS,PetscBool),arg1,arg2)
end 
=#
function TSRosWRegister(arg0::Type{Float64},arg1::TSRosWType,arg2::Integer,arg3::Integer,PetscReal::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg7::Integer,arg8::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    err = ccall((:TSRosWRegister,petsc1),PetscErrorCode,(Cstring,PetscInt,PetscInt,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},PetscInt,Ptr{Cint}),arg1,arg2,arg3,PetscReal,arg4,arg5,arg6,arg7,arg8)
    return err
end

function TSRosWRegisterRos4(arg0::Type{Float64},arg1::TSRosWType,PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer)
    err = ccall((:TSRosWRegisterRos4,petsc1),PetscErrorCode,(Cstring,Cint,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4,arg5)
    return err
end

function TSRosWFinalizePackage(arg0::Type{Float64})
    err = ccall((:TSRosWFinalizePackage,petsc1),PetscErrorCode,())
    return err
end

function TSRosWInitializePackage(arg0::Type{Float64})
    err = ccall((:TSRosWInitializePackage,petsc1),PetscErrorCode,())
    return err
end

function TSRosWRegisterDestroy(arg0::Type{Float64})
    err = ccall((:TSRosWRegisterDestroy,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function TSThetaSetTheta(arg0::Type{Float64},arg1::TS,PetscReal::Integer)
    ccall((:TSThetaSetTheta,petsc1),PetscErrorCode,(TS,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function TSThetaGetTheta(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TSThetaGetTheta,petsc1),PetscErrorCode,(TS,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSThetaGetEndpoint(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:TSThetaGetEndpoint,petsc1),PetscErrorCode,(TS,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSThetaSetEndpoint(arg0::Type{Float64},arg1::TS,arg2::PetscBool)
    ccall((:TSThetaSetEndpoint,petsc1),PetscErrorCode,(TS,PetscBool),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSAlphaSetAdapt(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSAlphaSetAdapt,petsc1),PetscErrorCode,(TS,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TSAlphaAdaptDefault(arg1::TS,PetscReal::Integer,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg6::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TSAlphaAdaptDefault,petsc1),PetscErrorCode,(TS,Cint,Vec{Float64},Vec{Float64},Ptr{Cint},Ptr{PetscBool},Ptr{Void}),arg1,PetscReal,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function TSAlphaSetRadius(arg0::Type{Float64},arg1::TS,PetscReal::Integer)
    ccall((:TSAlphaSetRadius,petsc1),PetscErrorCode,(TS,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function TSAlphaSetParams(arg0::Type{Float64},arg1::TS,PetscReal::Integer,arg2::Integer,arg3::Integer)
    ccall((:TSAlphaSetParams,petsc1),PetscErrorCode,(TS,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TSAlphaGetParams(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TSAlphaGetParams,petsc1),PetscErrorCode,(TS,Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TSSetDM(arg0::Type{Float64},arg1::TS,arg2::DM)
    ccall((:TSSetDM,petsc1),PetscErrorCode,(TS,DM),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TSGetDM(arg0::Type{Float64},arg1::TS,arg2::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:TSGetDM,petsc1),PetscErrorCode,(TS,Ptr{DM}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function SNESTSFormFunction(arg1::SNES,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESTSFormFunction,petsc1),PetscErrorCode,(SNES,Vec{Float64},Vec{Float64},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function SNESTSFormJacobian(arg1::SNES,arg2::Vec{Float64},arg3::Mat{Float64},arg4::Mat{Float64},arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:SNESTSFormJacobian,petsc1),PetscErrorCode,(SNES,Vec{Float64},Mat{Float64},Mat{Float64},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end 
=#
function VecFischer(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64},arg5::Vec{Float64})
    err = ccall((:VecFischer,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5)
    return err
end

function VecSFischer(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64},PetscReal::Integer,arg5::Vec{Float64})
    err = ccall((:VecSFischer,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Cint,Vec{Float64}),arg1,arg2,arg3,arg4,PetscReal,arg5)
    return err
end

function MatDFischer(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64},arg5::Vec{Float64},arg6::Vec{Float64},arg7::Vec{Float64},arg8::Vec{Float64},arg9::Vec{Float64})
    err = ccall((:MatDFischer,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
    return err
end

function MatDSFischer(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64},arg5::Vec{Float64},PetscReal::Integer,arg6::Vec{Float64},arg7::Vec{Float64},arg8::Vec{Float64},arg9::Vec{Float64},arg10::Vec{Float64})
    err = ccall((:MatDSFischer,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Cint,Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5,PetscReal,arg6,arg7,arg8,arg9,arg10)
    return err
end

function TaoInitializePackage(arg0::Type{Float64})
    err = ccall((:TaoInitializePackage,petsc1),PetscErrorCode,())
    return err
end

function TaoFinalizePackage(arg0::Type{Float64})
    err = ccall((:TaoFinalizePackage,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function TaoCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{Tao},StridedArray{Tao},Ptr{Void}))
    ccall((:TaoCreate,petsc1),PetscErrorCode,(comm_type,Ptr{Tao}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetFromOptions(arg0::Type{Float64},arg1::Tao)
    ccall((:TaoSetFromOptions,petsc1),PetscErrorCode,(Tao,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetUp(arg0::Type{Float64},arg1::Tao)
    ccall((:TaoSetUp,petsc1),PetscErrorCode,(Tao,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetType(arg0::Type{Float64},arg1::Tao,arg2::Union(ByteString,Symbol))
    ccall((:TaoSetType,petsc1),PetscErrorCode,(Tao,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGetType(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:TaoGetType,petsc1),PetscErrorCode,(Tao,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetApplicationContext(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoSetApplicationContext,petsc1),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGetApplicationContext(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoGetApplicationContext,petsc1),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoDestroy(arg0::Type{Float64},arg1::Union(Ptr{Tao},StridedArray{Tao},Ptr{Void}))
    ccall((:TaoDestroy,petsc1),PetscErrorCode,(Ptr{Tao},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetOptionsPrefix(arg0::Type{Float64},arg1::Tao,arg2::Union(ByteString,Symbol))
    ccall((:TaoSetOptionsPrefix,petsc1),PetscErrorCode,(Tao,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoView(arg1::Tao,arg2::PetscViewer{Float64})
    ccall((:TaoView,petsc1),PetscErrorCode,(Tao,PetscViewer{Float64}),arg1,arg2)
end 
=#
function TaoRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:TaoRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

function TaoRegisterDestroy(arg0::Type{Float64})
    err = ccall((:TaoRegisterDestroy,petsc1),PetscErrorCode,())
    return err
end

#= skipping function with undefined symbols: 
 function TaoGetConvergedReason(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{TaoConvergedReason},StridedArray{TaoConvergedReason},Ptr{Void}))
    ccall((:TaoGetConvergedReason,petsc1),PetscErrorCode,(Tao,Ptr{TaoConvergedReason}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGetSolutionStatus(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg7::Union(Ptr{TaoConvergedReason},StridedArray{TaoConvergedReason},Ptr{Void}))
    ccall((:TaoGetSolutionStatus,petsc1),PetscErrorCode,(Tao,Ptr{PetscInt},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{TaoConvergedReason}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetConvergedReason(arg0::Type{Float64},arg1::Tao,arg2::TaoConvergedReason)
    ccall((:TaoSetConvergedReason,petsc1),PetscErrorCode,(Tao,TaoConvergedReason),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetInitialVector(arg1::Tao,arg2::Vec{Float64})
    ccall((:TaoSetInitialVector,petsc1),PetscErrorCode,(Tao,Vec{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGetSolutionVector(arg1::Tao,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:TaoGetSolutionVector,petsc1),PetscErrorCode,(Tao,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGetGradientVector(arg1::Tao,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:TaoGetGradientVector,petsc1),PetscErrorCode,(Tao,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetObjectiveRoutine(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoSetObjectiveRoutine,petsc1),PetscErrorCode,(Tao,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetGradientRoutine(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoSetGradientRoutine,petsc1),PetscErrorCode,(Tao,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetObjectiveAndGradientRoutine(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoSetObjectiveAndGradientRoutine,petsc1),PetscErrorCode,(Tao,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetHessianRoutine(arg1::Tao,arg2::Mat{Float64},arg3::Mat{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoSetHessianRoutine,petsc1),PetscErrorCode,(Tao,Mat{Float64},Mat{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetSeparableObjectiveRoutine(arg1::Tao,arg2::Vec{Float64},arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoSetSeparableObjectiveRoutine,petsc1),PetscErrorCode,(Tao,Vec{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetConstraintsRoutine(arg1::Tao,arg2::Vec{Float64},arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoSetConstraintsRoutine,petsc1),PetscErrorCode,(Tao,Vec{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetInequalityConstraintsRoutine(arg1::Tao,arg2::Vec{Float64},arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoSetInequalityConstraintsRoutine,petsc1),PetscErrorCode,(Tao,Vec{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetEqualityConstraintsRoutine(arg1::Tao,arg2::Vec{Float64},arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoSetEqualityConstraintsRoutine,petsc1),PetscErrorCode,(Tao,Vec{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetJacobianRoutine(arg1::Tao,arg2::Mat{Float64},arg3::Mat{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoSetJacobianRoutine,petsc1),PetscErrorCode,(Tao,Mat{Float64},Mat{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetJacobianStateRoutine(arg1::Tao,arg2::Mat{Float64},arg3::Mat{Float64},arg4::Mat{Float64},arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg6::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoSetJacobianStateRoutine,petsc1),PetscErrorCode,(Tao,Mat{Float64},Mat{Float64},Mat{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetJacobianDesignRoutine(arg1::Tao,arg2::Mat{Float64},arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoSetJacobianDesignRoutine,petsc1),PetscErrorCode,(Tao,Mat{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetJacobianInequalityRoutine(arg1::Tao,arg2::Mat{Float64},arg3::Mat{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoSetJacobianInequalityRoutine,petsc1),PetscErrorCode,(Tao,Mat{Float64},Mat{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetJacobianEqualityRoutine(arg1::Tao,arg2::Mat{Float64},arg3::Mat{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoSetJacobianEqualityRoutine,petsc1),PetscErrorCode,(Tao,Mat{Float64},Mat{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetStateDesignIS(arg1::Tao,arg2::IS{Float64},arg3::IS{Float64})
    ccall((:TaoSetStateDesignIS,petsc1),PetscErrorCode,(Tao,IS{Float64},IS{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoComputeObjective(arg1::Tao,arg2::Vec{Float64},arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TaoComputeObjective,petsc1),PetscErrorCode,(Tao,Vec{Float64},Ptr{Cint}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoComputeSeparableObjective(arg1::Tao,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:TaoComputeSeparableObjective,petsc1),PetscErrorCode,(Tao,Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoComputeGradient(arg1::Tao,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:TaoComputeGradient,petsc1),PetscErrorCode,(Tao,Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoComputeObjectiveAndGradient(arg1::Tao,arg2::Vec{Float64},arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Vec{Float64})
    ccall((:TaoComputeObjectiveAndGradient,petsc1),PetscErrorCode,(Tao,Vec{Float64},Ptr{Cint},Vec{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TaoComputeConstraints(arg1::Tao,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:TaoComputeConstraints,petsc1),PetscErrorCode,(Tao,Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoComputeInequalityConstraints(arg1::Tao,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:TaoComputeInequalityConstraints,petsc1),PetscErrorCode,(Tao,Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoComputeEqualityConstraints(arg1::Tao,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:TaoComputeEqualityConstraints,petsc1),PetscErrorCode,(Tao,Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoDefaultComputeGradient(arg1::Tao,arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoDefaultComputeGradient,petsc1),PetscErrorCode,(Tao,Vec{Float64},Vec{Float64},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TaoIsObjectiveDefined(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:TaoIsObjectiveDefined,petsc1),PetscErrorCode,(Tao,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoIsGradientDefined(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:TaoIsGradientDefined,petsc1),PetscErrorCode,(Tao,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoIsObjectiveAndGradientDefined(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:TaoIsObjectiveAndGradientDefined,petsc1),PetscErrorCode,(Tao,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoComputeHessian(arg1::Tao,arg2::Vec{Float64},arg3::Mat{Float64},arg4::Mat{Float64})
    ccall((:TaoComputeHessian,petsc1),PetscErrorCode,(Tao,Vec{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TaoComputeJacobian(arg1::Tao,arg2::Vec{Float64},arg3::Mat{Float64},arg4::Mat{Float64})
    ccall((:TaoComputeJacobian,petsc1),PetscErrorCode,(Tao,Vec{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TaoComputeJacobianState(arg1::Tao,arg2::Vec{Float64},arg3::Mat{Float64},arg4::Mat{Float64},arg5::Mat{Float64})
    ccall((:TaoComputeJacobianState,petsc1),PetscErrorCode,(Tao,Vec{Float64},Mat{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function TaoComputeJacobianEquality(arg1::Tao,arg2::Vec{Float64},arg3::Mat{Float64},arg4::Mat{Float64})
    ccall((:TaoComputeJacobianEquality,petsc1),PetscErrorCode,(Tao,Vec{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TaoComputeJacobianInequality(arg1::Tao,arg2::Vec{Float64},arg3::Mat{Float64},arg4::Mat{Float64})
    ccall((:TaoComputeJacobianInequality,petsc1),PetscErrorCode,(Tao,Vec{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TaoComputeJacobianDesign(arg1::Tao,arg2::Vec{Float64},arg3::Mat{Float64})
    ccall((:TaoComputeJacobianDesign,petsc1),PetscErrorCode,(Tao,Vec{Float64},Mat{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoDefaultComputeHessian(arg1::Tao,arg2::Vec{Float64},arg3::Mat{Float64},arg4::Mat{Float64},arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoDefaultComputeHessian,petsc1),PetscErrorCode,(Tao,Vec{Float64},Mat{Float64},Mat{Float64},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function TaoDefaultComputeHessianColor(arg1::Tao,arg2::Vec{Float64},arg3::Mat{Float64},arg4::Mat{Float64},arg5::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoDefaultComputeHessianColor,petsc1),PetscErrorCode,(Tao,Vec{Float64},Mat{Float64},Mat{Float64},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function TaoComputeDualVariables(arg1::Tao,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:TaoComputeDualVariables,petsc1),PetscErrorCode,(Tao,Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoComputeDualVariables(arg1::Tao,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:TaoComputeDualVariables,petsc1),PetscErrorCode,(Tao,Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetVariableBounds(arg1::Tao,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:TaoSetVariableBounds,petsc1),PetscErrorCode,(Tao,Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGetVariableBounds(arg1::Tao,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:TaoGetVariableBounds,petsc1),PetscErrorCode,(Tao,Ptr{Vec{Float64}},Ptr{Vec{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGetDualVariables(arg1::Tao,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:TaoGetDualVariables,petsc1),PetscErrorCode,(Tao,Ptr{Vec{Float64}},Ptr{Vec{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetInequalityBounds(arg1::Tao,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:TaoSetInequalityBounds,petsc1),PetscErrorCode,(Tao,Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGetInequalityBounds(arg1::Tao,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:TaoGetInequalityBounds,petsc1),PetscErrorCode,(Tao,Ptr{Vec{Float64}},Ptr{Vec{Float64}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetVariableBoundsRoutine(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoSetVariableBoundsRoutine,petsc1),PetscErrorCode,(Tao,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoComputeVariableBounds(arg0::Type{Float64},arg1::Tao)
    ccall((:TaoComputeVariableBounds,petsc1),PetscErrorCode,(Tao,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGetTolerances(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TaoGetTolerances,petsc1),PetscErrorCode,(Tao,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetTolerances(arg0::Type{Float64},arg1::Tao,PetscReal::Integer,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer)
    ccall((:TaoSetTolerances,petsc1),PetscErrorCode,(Tao,Cint,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGetConstraintTolerances(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TaoGetConstraintTolerances,petsc1),PetscErrorCode,(Tao,Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetConstraintTolerances(arg0::Type{Float64},arg1::Tao,PetscReal::Integer,arg2::Integer)
    ccall((:TaoSetConstraintTolerances,petsc1),PetscErrorCode,(Tao,Cint,Cint),arg1,PetscReal,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetFunctionLowerBound(arg0::Type{Float64},arg1::Tao,PetscReal::Integer)
    ccall((:TaoSetFunctionLowerBound,petsc1),PetscErrorCode,(Tao,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetInitialTrustRegionRadius(arg0::Type{Float64},arg1::Tao,PetscReal::Integer)
    ccall((:TaoSetInitialTrustRegionRadius,petsc1),PetscErrorCode,(Tao,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetMaximumIterations(arg0::Type{Float64},arg1::Tao,arg2::Integer)
    ccall((:TaoSetMaximumIterations,petsc1),PetscErrorCode,(Tao,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetMaximumFunctionEvaluations(arg0::Type{Float64},arg1::Tao,arg2::Integer)
    ccall((:TaoSetMaximumFunctionEvaluations,petsc1),PetscErrorCode,(Tao,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGetFunctionLowerBound(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TaoGetFunctionLowerBound,petsc1),PetscErrorCode,(Tao,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGetInitialTrustRegionRadius(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TaoGetInitialTrustRegionRadius,petsc1),PetscErrorCode,(Tao,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGetCurrentTrustRegionRadius(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TaoGetCurrentTrustRegionRadius,petsc1),PetscErrorCode,(Tao,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGetMaximumIterations(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:TaoGetMaximumIterations,petsc1),PetscErrorCode,(Tao,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGetCurrentFunctionEvaluations(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:TaoGetCurrentFunctionEvaluations,petsc1),PetscErrorCode,(Tao,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGetMaximumFunctionEvaluations(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:TaoGetMaximumFunctionEvaluations,petsc1),PetscErrorCode,(Tao,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGetIterationNumber(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:TaoGetIterationNumber,petsc1),PetscErrorCode,(Tao,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetIterationNumber(arg0::Type{Float64},arg1::Tao,arg2::Integer)
    ccall((:TaoSetIterationNumber,petsc1),PetscErrorCode,(Tao,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGetTotalIterationNumber(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:TaoGetTotalIterationNumber,petsc1),PetscErrorCode,(Tao,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetTotalIterationNumber(arg0::Type{Float64},arg1::Tao,arg2::Integer)
    ccall((:TaoSetTotalIterationNumber,petsc1),PetscErrorCode,(Tao,PetscInt),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetOptionsPrefix(arg0::Type{Float64},arg1::Tao,p::Union(ByteString,Symbol))
    ccall((:TaoSetOptionsPrefix,petsc1),PetscErrorCode,(Tao,Cstring),arg1,p)
end 
=#
#= skipping function with undefined symbols: 
 function TaoAppendOptionsPrefix(arg0::Type{Float64},arg1::Tao,p::Union(ByteString,Symbol))
    ccall((:TaoAppendOptionsPrefix,petsc1),PetscErrorCode,(Tao,Cstring),arg1,p)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGetOptionsPrefix(arg0::Type{Float64},arg1::Tao,p::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:TaoGetOptionsPrefix,petsc1),PetscErrorCode,(Tao,Ptr{Ptr{Uint8}}),arg1,p)
end 
=#
#= skipping function with undefined symbols: 
 function TaoResetStatistics(arg0::Type{Float64},arg1::Tao)
    ccall((:TaoResetStatistics,petsc1),PetscErrorCode,(Tao,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGetKSP(arg1::Tao,arg2::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    ccall((:TaoGetKSP,petsc1),PetscErrorCode,(Tao,Ptr{KSP{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGetLinearSolveIterations(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:TaoGetLinearSolveIterations,petsc1),PetscErrorCode,(Tao,Ptr{PetscInt}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{TaoLineSearch},StridedArray{TaoLineSearch},Ptr{Void}))
    ccall((:TaoLineSearchCreate,petsc1),PetscErrorCode,(comm_type,Ptr{TaoLineSearch}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchSetFromOptions(arg0::Type{Float64},arg1::TaoLineSearch)
    ccall((:TaoLineSearchSetFromOptions,petsc1),PetscErrorCode,(TaoLineSearch,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchSetUp(arg0::Type{Float64},arg1::TaoLineSearch)
    ccall((:TaoLineSearchSetUp,petsc1),PetscErrorCode,(TaoLineSearch,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchDestroy(arg0::Type{Float64},arg1::Union(Ptr{TaoLineSearch},StridedArray{TaoLineSearch},Ptr{Void}))
    ccall((:TaoLineSearchDestroy,petsc1),PetscErrorCode,(Ptr{TaoLineSearch},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchView(arg1::TaoLineSearch,arg2::PetscViewer{Float64})
    ccall((:TaoLineSearchView,petsc1),PetscErrorCode,(TaoLineSearch,PetscViewer{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchReset(arg0::Type{Float64},arg1::TaoLineSearch)
    ccall((:TaoLineSearchReset,petsc1),PetscErrorCode,(TaoLineSearch,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchAppendOptionsPrefix(arg0::Type{Float64},arg1::TaoLineSearch,prefix::Union(ByteString,Symbol))
    ccall((:TaoLineSearchAppendOptionsPrefix,petsc1),PetscErrorCode,(TaoLineSearch,Cstring),arg1,prefix)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchGetOptionsPrefix(arg0::Type{Float64},arg1::TaoLineSearch,prefix::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:TaoLineSearchGetOptionsPrefix,petsc1),PetscErrorCode,(TaoLineSearch,Ptr{Ptr{Uint8}}),arg1,prefix)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchApply(arg1::TaoLineSearch,arg2::Vec{Float64},arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Vec{Float64},arg5::Vec{Float64},arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg7::Union(Ptr{TaoLineSearchConvergedReason},StridedArray{TaoLineSearchConvergedReason},Ptr{Void}))
    ccall((:TaoLineSearchApply,petsc1),PetscErrorCode,(TaoLineSearch,Vec{Float64},Ptr{Cint},Vec{Float64},Vec{Float64},Ptr{Cint},Ptr{TaoLineSearchConvergedReason}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchGetStepLength(arg0::Type{Float64},arg1::TaoLineSearch,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TaoLineSearchGetStepLength,petsc1),PetscErrorCode,(TaoLineSearch,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchGetStartingVector(arg1::TaoLineSearch,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:TaoLineSearchGetStartingVector,petsc1),PetscErrorCode,(TaoLineSearch,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchGetStepDirection(arg1::TaoLineSearch,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:TaoLineSearchGetStepDirection,petsc1),PetscErrorCode,(TaoLineSearch,Ptr{Vec{Float64}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchSetInitialStepLength(arg0::Type{Float64},arg1::TaoLineSearch,PetscReal::Integer)
    ccall((:TaoLineSearchSetInitialStepLength,petsc1),PetscErrorCode,(TaoLineSearch,Cint),arg1,PetscReal)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchGetSolution(arg1::TaoLineSearch,arg2::Vec{Float64},arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Vec{Float64},arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg6::Union(Ptr{TaoLineSearchConvergedReason},StridedArray{TaoLineSearchConvergedReason},Ptr{Void}))
    ccall((:TaoLineSearchGetSolution,petsc1),PetscErrorCode,(TaoLineSearch,Vec{Float64},Ptr{Cint},Vec{Float64},Ptr{Cint},Ptr{TaoLineSearchConvergedReason}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchGetFullStepObjective(arg0::Type{Float64},arg1::TaoLineSearch,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TaoLineSearchGetFullStepObjective,petsc1),PetscErrorCode,(TaoLineSearch,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchGetNumberFunctionEvaluations(arg0::Type{Float64},arg1::TaoLineSearch,arg2::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg3::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg4::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:TaoLineSearchGetNumberFunctionEvaluations,petsc1),PetscErrorCode,(TaoLineSearch,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchGetType(arg0::Type{Float64},arg1::TaoLineSearch,arg2::Union(Ptr{Ptr{Uint8}},StridedArray{Ptr{Uint8}},Ptr{Void}))
    ccall((:TaoLineSearchGetType,petsc1),PetscErrorCode,(TaoLineSearch,Ptr{Ptr{Uint8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchSetType(arg0::Type{Float64},arg1::TaoLineSearch,arg2::Union(ByteString,Symbol))
    ccall((:TaoLineSearchSetType,petsc1),PetscErrorCode,(TaoLineSearch,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchIsUsingTaoRoutines(arg0::Type{Float64},arg1::TaoLineSearch,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:TaoLineSearchIsUsingTaoRoutines,petsc1),PetscErrorCode,(TaoLineSearch,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchSetObjectiveAndGTSRoutine(arg0::Type{Float64},arg1::TaoLineSearch,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoLineSearchSetObjectiveAndGTSRoutine,petsc1),PetscErrorCode,(TaoLineSearch,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchSetObjectiveRoutine(arg0::Type{Float64},arg1::TaoLineSearch,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoLineSearchSetObjectiveRoutine,petsc1),PetscErrorCode,(TaoLineSearch,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchSetGradientRoutine(arg0::Type{Float64},arg1::TaoLineSearch,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoLineSearchSetGradientRoutine,petsc1),PetscErrorCode,(TaoLineSearch,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchSetObjectiveAndGradientRoutine(arg0::Type{Float64},arg1::TaoLineSearch,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoLineSearchSetObjectiveAndGradientRoutine,petsc1),PetscErrorCode,(TaoLineSearch,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchComputeObjective(arg1::TaoLineSearch,arg2::Vec{Float64},arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TaoLineSearchComputeObjective,petsc1),PetscErrorCode,(TaoLineSearch,Vec{Float64},Ptr{Cint}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchComputeGradient(arg1::TaoLineSearch,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:TaoLineSearchComputeGradient,petsc1),PetscErrorCode,(TaoLineSearch,Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchComputeObjectiveAndGradient(arg1::TaoLineSearch,arg2::Vec{Float64},arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Vec{Float64})
    ccall((:TaoLineSearchComputeObjectiveAndGradient,petsc1),PetscErrorCode,(TaoLineSearch,Vec{Float64},Ptr{Cint},Vec{Float64}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchComputeObjectiveAndGTS(arg1::TaoLineSearch,arg2::Vec{Float64},arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:TaoLineSearchComputeObjectiveAndGTS,petsc1),PetscErrorCode,(TaoLineSearch,Vec{Float64},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLineSearchSetVariableBounds(arg1::TaoLineSearch,arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:TaoLineSearchSetVariableBounds,petsc1),PetscErrorCode,(TaoLineSearch,Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end 
=#
function TaoLineSearchInitializePackage(arg0::Type{Float64})
    err = ccall((:TaoLineSearchInitializePackage,petsc1),PetscErrorCode,())
    return err
end

function TaoLineSearchFinalizePackage(arg0::Type{Float64})
    err = ccall((:TaoLineSearchFinalizePackage,petsc1),PetscErrorCode,())
    return err
end

function TaoLineSearchRegister(arg0::Type{Float64},arg1::Union(ByteString,Symbol),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    err = ccall((:TaoLineSearchRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
    return err
end

#= skipping function with undefined symbols: 
 function TaoLineSearchUseTaoRoutines(arg0::Type{Float64},arg1::TaoLineSearch,arg2::Tao)
    ccall((:TaoLineSearchUseTaoRoutines,petsc1),PetscErrorCode,(TaoLineSearch,Tao),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGetLineSearch(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{TaoLineSearch},StridedArray{TaoLineSearch},Ptr{Void}))
    ccall((:TaoGetLineSearch,petsc1),PetscErrorCode,(Tao,Ptr{TaoLineSearch}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetConvergenceHistory(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}),arg6::Integer,arg7::PetscBool)
    ccall((:TaoSetConvergenceHistory,petsc1),PetscErrorCode,(Tao,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{PetscInt},PetscInt,PetscBool),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGetConvergenceHistory(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg3::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg4::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg5::Union(Ptr{Ptr{PetscInt}},StridedArray{Ptr{PetscInt}},Ptr{Void}),arg6::Union(Ptr{PetscInt},StridedArray{PetscInt},Ptr{Void}))
    ccall((:TaoGetConvergenceHistory,petsc1),PetscErrorCode,(Tao,Ptr{Ptr{Cint}},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}},Ptr{Ptr{PetscInt}},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetMonitor(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoSetMonitor,petsc1),PetscErrorCode,(Tao,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function TaoCancelMonitors(arg0::Type{Float64},arg1::Tao)
    ccall((:TaoCancelMonitors,petsc1),PetscErrorCode,(Tao,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TaoDefaultMonitor(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoDefaultMonitor,petsc1),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoDefaultSMonitor(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoDefaultSMonitor,petsc1),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoDefaultCMonitor(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoDefaultCMonitor,petsc1),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSolutionMonitor(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoSolutionMonitor,petsc1),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSeparableObjectiveMonitor(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoSeparableObjectiveMonitor,petsc1),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoGradientMonitor(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoGradientMonitor,petsc1),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoStepDirectionMonitor(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoStepDirectionMonitor,petsc1),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoDrawSolutionMonitor(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoDrawSolutionMonitor,petsc1),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoDrawStepMonitor(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoDrawStepMonitor,petsc1),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoDrawGradientMonitor(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoDrawGradientMonitor,petsc1),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoAddLineSearchCounts(arg0::Type{Float64},arg1::Tao)
    ccall((:TaoAddLineSearchCounts,petsc1),PetscErrorCode,(Tao,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function TaoDefaultConvergenceTest(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoDefaultConvergenceTest,petsc1),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSetConvergenceTest(arg0::Type{Float64},arg1::Tao,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:TaoSetConvergenceTest,petsc1),PetscErrorCode,(Tao,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoSQPCONSetStateDesignIS(arg1::Tao,arg2::IS{Float64},arg3::IS{Float64})
    ccall((:TaoSQPCONSetStateDesignIS,petsc1),PetscErrorCode,(Tao,IS{Float64},IS{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoLCLSetStateDesignIS(arg1::Tao,arg2::IS{Float64},arg3::IS{Float64})
    ccall((:TaoLCLSetStateDesignIS,petsc1),PetscErrorCode,(Tao,IS{Float64},IS{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function TaoMonitor(arg0::Type{Float64},arg1::Tao,arg2::Integer,PetscReal::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{TaoConvergedReason},StridedArray{TaoConvergedReason},Ptr{Void}))
    ccall((:TaoMonitor,petsc1),PetscErrorCode,(Tao,PetscInt,Cint,Cint,Cint,Cint,Ptr{TaoConvergedReason}),arg1,arg2,PetscReal,arg3,arg4,arg5,arg6)
end 
=#
