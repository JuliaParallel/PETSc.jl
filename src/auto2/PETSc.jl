# Julia wrapper for header: /home/jared/build/petsc-3.6.0/include/petscksp.h
# Automatically generated using Clang.jl wrap_c, version 0.0.0


function PetscIsInfOrNanReal(arg0::Type{Float64})
    ccall((:PetscIsInfOrNanReal,petsc1),PetscErrorCode,())
end

function PetscIsNormalReal(arg0::Type{Float64})
    ccall((:PetscIsNormalReal,petsc1),PetscBool,())
end

function PetscSetHelpVersionFunctions(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscSetHelpVersionFunctions,petsc1),PetscErrorCode,(Ptr{Void},Ptr{Void}),arg1,arg2)
end

function PetscCommDuplicate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{MPI_Comm},StridedArray{MPI_Comm},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscCommDuplicate,petsc1),PetscErrorCode,(comm_type,Ptr{comm_type},Ptr{Cint}),arg1.val,arg2,arg3)
end

function PetscCommDestroy(arg0::Type{Float64},arg1::Union(Ptr{MPI_Comm},StridedArray{MPI_Comm},Ptr{Void}))
    ccall((:PetscCommDestroy,petsc1),PetscErrorCode,(Ptr{comm_type},),arg1)
end

function PetscMallocSet(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscMallocSet,petsc1),PetscErrorCode,(Ptr{Void},Ptr{Void}),arg1,arg2)
end

function PetscMallocClear(arg0::Type{Float64})
    ccall((:PetscMallocClear,petsc1),PetscErrorCode,())
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
    ccall((:PetscMallocGetCurrentUsage,petsc1),PetscErrorCode,(Ptr{PetscLogDouble},),arg1)
end

function PetscMallocGetMaximumUsage(arg0::Type{Float64},arg1::Union(Ptr{PetscLogDouble},StridedArray{PetscLogDouble},Ptr{Void}))
    ccall((:PetscMallocGetMaximumUsage,petsc1),PetscErrorCode,(Ptr{PetscLogDouble},),arg1)
end

function PetscMallocDebug(arg0::Type{Float64},arg1::PetscBool)
    ccall((:PetscMallocDebug,petsc1),PetscErrorCode,(PetscBool,),arg1)
end

function PetscMallocGetDebug(arg0::Type{Float64},arg1::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscMallocGetDebug,petsc1),PetscErrorCode,(Ptr{PetscBool},),arg1)
end

function PetscMallocValidate(arg1::Integer,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscMallocValidate,petsc1),PetscErrorCode,(Cint,Cstring,Cstring),arg1,arg2,arg3)
end

function PetscMallocSetDumpLog(arg0::Type{Float64})
    ccall((:PetscMallocSetDumpLog,petsc1),PetscErrorCode,())
end

function PetscMallocSetDumpLogThreshold(arg0::Type{Float64},arg1::PetscLogDouble)
    ccall((:PetscMallocSetDumpLogThreshold,petsc1),PetscErrorCode,(PetscLogDouble,),arg1)
end

function PetscMallocGetDumpLog(arg0::Type{Float64},arg1::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscMallocGetDumpLog,petsc1),PetscErrorCode,(Ptr{PetscBool},),arg1)
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
    ccall((:PetscDataTypeGetSize,petsc1),PetscErrorCode,(PetscDataType,Ptr{Cint}),arg1,arg2)
end

function PetscDataTypeFromString(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(Ptr{PetscDataType},StridedArray{PetscDataType},Ptr{Void}),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscDataTypeFromString,petsc1),PetscErrorCode,(Cstring,Ptr{PetscDataType},Ptr{PetscBool}),arg1,arg2,arg3)
end

function PetscBitMemcpy(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Integer,arg5::Integer,arg6::PetscDataType)
    ccall((:PetscBitMemcpy,petsc1),PetscErrorCode,(Ptr{Void},Int32,Ptr{Void},Int32,Int32,PetscDataType),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscMemmove(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),size_t::Integer)
    ccall((:PetscMemmove,petsc1),PetscErrorCode,(Ptr{Void},Ptr{Void},Cint),arg1,arg2,size_t)
end

function PetscMemcmp(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),size_t::Integer,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscMemcmp,petsc1),PetscErrorCode,(Ptr{Void},Ptr{Void},Cint,Ptr{PetscBool}),arg1,arg2,size_t,arg3)
end

function PetscStrlen(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscStrlen,petsc1),PetscErrorCode,(Cstring,Ptr{Cint}),arg1,arg2)
end

function PetscStrToArray(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Uint8,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Ptr{ASCIIString}},StridedArray{Ptr{ASCIIString}},Ptr{Void}))
    ccall((:PetscStrToArray,petsc1),PetscErrorCode,(Cstring,Uint8,Ptr{Cint},Ptr{Ptr{Ptr{UInt8}}}),arg1,arg2,arg3,arg4)
end

function PetscStrToArrayDestroy(arg1::Integer,arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PetscStrToArrayDestroy,petsc1),PetscErrorCode,(Cint,Ptr{Ptr{UInt8}}),arg1,arg2)
end

function PetscStrcmp(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscStrcmp,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{PetscBool}),arg1,arg2,arg3)
end

function PetscStrgrt(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscStrgrt,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{PetscBool}),arg1,arg2,arg3)
end

function PetscStrcasecmp(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscStrcasecmp,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{PetscBool}),arg1,arg2,arg3)
end

function PetscStrncmp(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscStrncmp,petsc1),PetscErrorCode,(Cstring,Cstring,Cint,Ptr{PetscBool}),arg1,arg2,size_t,arg3)
end

function PetscStrcpy(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscStrcpy,petsc1),PetscErrorCode,(Cstring,Cstring),arg1,arg2)
end

function PetscStrcat(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscStrcat,petsc1),PetscErrorCode,(Cstring,Cstring),arg1,arg2)
end

function PetscStrncat(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer)
    ccall((:PetscStrncat,petsc1),PetscErrorCode,(Cstring,Cstring,Cint),arg1,arg2,size_t)
end

function PetscStrncpy(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer)
    ccall((:PetscStrncpy,petsc1),PetscErrorCode,(Cstring,Cstring,Cint),arg1,arg2,size_t)
end

function PetscStrchr(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Uint8,arg3::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PetscStrchr,petsc1),PetscErrorCode,(Cstring,Uint8,Ptr{Ptr{UInt8}}),arg1,arg2,arg3)
end

function PetscStrtolower(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscStrtolower,petsc1),PetscErrorCode,(Cstring,),arg1)
end

function PetscStrtoupper(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscStrtoupper,petsc1),PetscErrorCode,(Cstring,),arg1)
end

function PetscStrrchr(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Uint8,arg3::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PetscStrrchr,petsc1),PetscErrorCode,(Cstring,Uint8,Ptr{Ptr{UInt8}}),arg1,arg2,arg3)
end

function PetscStrstr(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PetscStrstr,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Ptr{UInt8}}),arg1,arg2,arg3)
end

function PetscStrrstr(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PetscStrrstr,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Ptr{UInt8}}),arg1,arg2,arg3)
end

function PetscStrendswith(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscStrendswith,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{PetscBool}),arg1,arg2,arg3)
end

function PetscStrbeginswith(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscStrbeginswith,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{PetscBool}),arg1,arg2,arg3)
end

function PetscStrendswithwhich(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscStrendswithwhich,petsc1),PetscErrorCode,(Cstring,Ptr{Ptr{UInt8}},Ptr{Int32}),arg1,arg2,arg3)
end

function PetscStrallocpy(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PetscStrallocpy,petsc1),PetscErrorCode,(Cstring,Ptr{Ptr{UInt8}}),arg1,arg2)
end

function PetscStrArrayallocpy(arg1::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}),arg2::Union(Ptr{Ptr{ASCIIString}},StridedArray{Ptr{ASCIIString}},Ptr{Void}))
    ccall((:PetscStrArrayallocpy,petsc1),PetscErrorCode,(Ptr{Ptr{UInt8}},Ptr{Ptr{Ptr{UInt8}}}),arg1,arg2)
end

function PetscStrArrayDestroy(arg1::Union(Ptr{Ptr{ASCIIString}},StridedArray{Ptr{ASCIIString}},Ptr{Void}))
    ccall((:PetscStrArrayDestroy,petsc1),PetscErrorCode,(Ptr{Ptr{Ptr{UInt8}}},),arg1)
end

function PetscStrNArrayallocpy(arg1::Integer,arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}),arg3::Union(Ptr{Ptr{ASCIIString}},StridedArray{Ptr{ASCIIString}},Ptr{Void}))
    ccall((:PetscStrNArrayallocpy,petsc1),PetscErrorCode,(Int32,Ptr{Ptr{UInt8}},Ptr{Ptr{Ptr{UInt8}}}),arg1,arg2,arg3)
end

function PetscStrNArrayDestroy(arg1::Integer,arg2::Union(Ptr{Ptr{ASCIIString}},StridedArray{Ptr{ASCIIString}},Ptr{Void}))
    ccall((:PetscStrNArrayDestroy,petsc1),PetscErrorCode,(Int32,Ptr{Ptr{Ptr{UInt8}}}),arg1,arg2)
end

function PetscStrreplace(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer)
    ccall((:PetscStrreplace,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint),arg1.val,arg2,arg3,size_t)
end

function PetscStrcmpNoError(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscStrcmpNoError,petsc1),Void,(Cstring,Cstring,Ptr{PetscBool}),arg1,arg2,arg3)
end

#= skipping function with undefined symbols: 
 function PetscTokenCreate(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Uint8,arg3::Union(Ptr{PetscToken},StridedArray{PetscToken},Ptr{Void}))
    ccall((:PetscTokenCreate,petsc1),PetscErrorCode,(Cstring,Uint8,Ptr{PetscToken}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscTokenFind(arg1::PetscToken,arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PetscTokenFind,petsc1),PetscErrorCode,(PetscToken,Ptr{Ptr{UInt8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscTokenDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscToken},StridedArray{PetscToken},Ptr{Void}))
    ccall((:PetscTokenDestroy,petsc1),PetscErrorCode,(Ptr{PetscToken},),arg1)
end 
=#
function PetscEListFind(arg1::Integer,arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscEListFind,petsc1),PetscErrorCode,(Int32,Ptr{Ptr{UInt8}},Cstring,Ptr{Int32},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function PetscEnumFind(arg1::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{PetscEnum},StridedArray{PetscEnum},Ptr{Void}),arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscEnumFind,petsc1),PetscErrorCode,(Ptr{Ptr{UInt8}},Cstring,Ptr{PetscEnum},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function PetscMaxSum(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscMaxSum,petsc1),PetscErrorCode,(comm_type,Ptr{Int32},Ptr{Int32},Ptr{Int32}),arg1.val,arg2,arg3,arg4)
end

#= skipping function with undefined symbols: 
 function MPIULong_Send(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::Integer,arg3::MPI_Datatype,arg4::PetscMPIInt,arg5::PetscMPIInt,arg6::MPI_Comm)
    ccall((:MPIULong_Send,petsc1),PetscErrorCode,(Ptr{Void},Int32,MPI_Datatype,PetscMPIInt,PetscMPIInt,comm_type),arg1,arg2,arg3,arg4,arg5,arg6.val)
end 
=#
#= skipping function with undefined symbols: 
 function MPIULong_Recv(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::Integer,arg3::MPI_Datatype,arg4::PetscMPIInt,arg5::PetscMPIInt,arg6::MPI_Comm)
    ccall((:MPIULong_Recv,petsc1),PetscErrorCode,(Ptr{Void},Int32,MPI_Datatype,PetscMPIInt,PetscMPIInt,comm_type),arg1,arg2,arg3,arg4,arg5,arg6.val)
end 
=#
function PetscErrorPrintfInitialize(arg0::Type{Float64})
    ccall((:PetscErrorPrintfInitialize,petsc1),PetscErrorCode,())
end

function PetscErrorMessage(arg1::Integer,arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}),arg3::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PetscErrorMessage,petsc1),PetscErrorCode,(Cint,Ptr{Ptr{UInt8}},Ptr{Ptr{UInt8}}),arg1,arg2,arg3)
end

function PetscTraceBackErrorHandler(arg1::MPI_Comm,arg2::Integer,arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::PetscErrorCode,arg6::PetscErrorType,arg7::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg8::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscTraceBackErrorHandler,petsc1),PetscErrorCode,(comm_type,Cint,Cstring,Cstring,PetscErrorCode,PetscErrorType,Cstring,Ptr{Void}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscIgnoreErrorHandler(arg1::MPI_Comm,arg2::Integer,arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::PetscErrorCode,arg6::PetscErrorType,arg7::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg8::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscIgnoreErrorHandler,petsc1),PetscErrorCode,(comm_type,Cint,Cstring,Cstring,PetscErrorCode,PetscErrorType,Cstring,Ptr{Void}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscEmacsClientErrorHandler(arg1::MPI_Comm,arg2::Integer,arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::PetscErrorCode,arg6::PetscErrorType,arg7::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg8::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscEmacsClientErrorHandler,petsc1),PetscErrorCode,(comm_type,Cint,Cstring,Cstring,PetscErrorCode,PetscErrorType,Cstring,Ptr{Void}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscMPIAbortErrorHandler(arg1::MPI_Comm,arg2::Integer,arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::PetscErrorCode,arg6::PetscErrorType,arg7::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg8::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscMPIAbortErrorHandler,petsc1),PetscErrorCode,(comm_type,Cint,Cstring,Cstring,PetscErrorCode,PetscErrorType,Cstring,Ptr{Void}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscAbortErrorHandler(arg1::MPI_Comm,arg2::Integer,arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::PetscErrorCode,arg6::PetscErrorType,arg7::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg8::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscAbortErrorHandler,petsc1),PetscErrorCode,(comm_type,Cint,Cstring,Cstring,PetscErrorCode,PetscErrorType,Cstring,Ptr{Void}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscAttachDebuggerErrorHandler(arg1::MPI_Comm,arg2::Integer,arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::PetscErrorCode,arg6::PetscErrorType,arg7::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg8::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscAttachDebuggerErrorHandler,petsc1),PetscErrorCode,(comm_type,Cint,Cstring,Cstring,PetscErrorCode,PetscErrorType,Cstring,Ptr{Void}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscReturnErrorHandler(arg1::MPI_Comm,arg2::Integer,arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::PetscErrorCode,arg6::PetscErrorType,arg7::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg8::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscReturnErrorHandler,petsc1),PetscErrorCode,(comm_type,Cint,Cstring,Cstring,PetscErrorCode,PetscErrorType,Cstring,Ptr{Void}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscPushErrorHandler(arg0::Type{Float64},handler::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscPushErrorHandler,petsc1),PetscErrorCode,(Ptr{Void},Ptr{Void}),handler,arg1)
end

function PetscPopErrorHandler(arg0::Type{Float64})
    ccall((:PetscPopErrorHandler,petsc1),PetscErrorCode,())
end

function PetscSignalHandlerDefault(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscSignalHandlerDefault,petsc1),PetscErrorCode,(Cint,Ptr{Void}),arg1,arg2)
end

function PetscPushSignalHandler(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscPushSignalHandler,petsc1),PetscErrorCode,(Ptr{Void},Ptr{Void}),arg1,arg2)
end

function PetscPopSignalHandler(arg0::Type{Float64})
    ccall((:PetscPopSignalHandler,petsc1),PetscErrorCode,())
end

function PetscCheckPointerSetIntensity(arg0::Type{Float64},arg1::Integer)
    ccall((:PetscCheckPointerSetIntensity,petsc1),PetscErrorCode,(Int32,),arg1)
end

function PetscSetFPTrap(arg0::Type{Float64},arg1::PetscFPTrap)
    ccall((:PetscSetFPTrap,petsc1),PetscErrorCode,(PetscFPTrap,),arg1)
end

function PetscFPTrapPush(arg0::Type{Float64},arg1::PetscFPTrap)
    ccall((:PetscFPTrapPush,petsc1),PetscErrorCode,(PetscFPTrap,),arg1)
end

function PetscFPTrapPop(arg0::Type{Float64})
    ccall((:PetscFPTrapPop,petsc1),PetscErrorCode,())
end

function PetscStackCopy(arg0::Type{Float64},arg1::Union(Ptr{PetscStack},StridedArray{PetscStack},Ptr{Void}),arg2::Union(Ptr{PetscStack},StridedArray{PetscStack},Ptr{Void}))
    ccall((:PetscStackCopy,petsc1),PetscErrorCode,(Ptr{PetscStack},Ptr{PetscStack}),arg1,arg2)
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
    ccall((:PetscStackDestroy,petsc1),PetscErrorCode,())
end

function PetscClassIdRegister(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(Ptr{PetscClassId},StridedArray{PetscClassId},Ptr{Void}))
    ccall((:PetscClassIdRegister,petsc1),PetscErrorCode,(Cstring,Ptr{PetscClassId}),arg1,arg2)
end

function PetscMemoryGetCurrentUsage(arg0::Type{Float64},arg1::Union(Ptr{PetscLogDouble},StridedArray{PetscLogDouble},Ptr{Void}))
    ccall((:PetscMemoryGetCurrentUsage,petsc1),PetscErrorCode,(Ptr{PetscLogDouble},),arg1)
end

function PetscMemoryGetMaximumUsage(arg0::Type{Float64},arg1::Union(Ptr{PetscLogDouble},StridedArray{PetscLogDouble},Ptr{Void}))
    ccall((:PetscMemoryGetMaximumUsage,petsc1),PetscErrorCode,(Ptr{PetscLogDouble},),arg1)
end

function PetscMemorySetGetMaximumUsage(arg0::Type{Float64})
    ccall((:PetscMemorySetGetMaximumUsage,petsc1),PetscErrorCode,())
end

function PetscMemoryTrace(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscMemoryTrace,petsc1),PetscErrorCode,(Cstring,),arg1)
end

function PetscInfoAllow(arg1::PetscBool,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscInfoAllow,petsc1),PetscErrorCode,(PetscBool,Cstring),arg1,arg2)
end

function PetscSleep(arg0::Type{Float64})
    ccall((:PetscSleep,petsc1),PetscErrorCode,())
end

function PetscInitialize(arg1::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg2::Union(Ptr{Ptr{ASCIIString}},StridedArray{Ptr{ASCIIString}},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscInitialize,petsc1),PetscErrorCode,(Ptr{Cint},Ptr{Ptr{Ptr{UInt8}}},Cstring,Cstring),arg1,arg2,arg3,arg4)
end

function PetscInitializeNoPointers(arg1::Integer,arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscInitializeNoPointers,petsc1),PetscErrorCode,(Cint,Ptr{Ptr{UInt8}},Cstring,Cstring),arg1,arg2,arg3,arg4)
end

function PetscInitializeNoArguments(arg0::Type{Float64})
    ccall((:PetscInitializeNoArguments,petsc1),PetscErrorCode,())
end

function PetscInitialized(arg0::Type{Float64},arg1::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscInitialized,petsc1),PetscErrorCode,(Ptr{PetscBool},),arg1)
end

function PetscFinalized(arg0::Type{Float64},arg1::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscFinalized,petsc1),PetscErrorCode,(Ptr{PetscBool},),arg1)
end

function PetscFinalize(arg0::Type{Float64})
    ccall((:PetscFinalize,petsc1),PetscErrorCode,())
end

function PetscInitializeFortran(arg0::Type{Float64})
    ccall((:PetscInitializeFortran,petsc1),PetscErrorCode,())
end

function PetscGetArgs(arg1::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg2::Union(Ptr{Ptr{ASCIIString}},StridedArray{Ptr{ASCIIString}},Ptr{Void}))
    ccall((:PetscGetArgs,petsc1),PetscErrorCode,(Ptr{Cint},Ptr{Ptr{Ptr{UInt8}}}),arg1,arg2)
end

function PetscGetArguments(arg1::Union(Ptr{Ptr{ASCIIString}},StridedArray{Ptr{ASCIIString}},Ptr{Void}))
    ccall((:PetscGetArguments,petsc1),PetscErrorCode,(Ptr{Ptr{Ptr{UInt8}}},),arg1)
end

function PetscFreeArguments(arg1::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PetscFreeArguments,petsc1),PetscErrorCode,(Ptr{Ptr{UInt8}},),arg1)
end

function PetscEnd(arg0::Type{Float64})
    ccall((:PetscEnd,petsc1),PetscErrorCode,())
end

function PetscSysInitializePackage(arg0::Type{Float64})
    ccall((:PetscSysInitializePackage,petsc1),PetscErrorCode,())
end

function PetscPythonInitialize(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscPythonInitialize,petsc1),PetscErrorCode,(Cstring,Cstring),arg1,arg2)
end

function PetscPythonFinalize(arg0::Type{Float64})
    ccall((:PetscPythonFinalize,petsc1),PetscErrorCode,())
end

function PetscPythonPrintError(arg0::Type{Float64})
    ccall((:PetscPythonPrintError,petsc1),PetscErrorCode,())
end

#= skipping function with undefined symbols: 
 function PetscPythonMonitorSet(arg1::PetscObject,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
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
 function PetscObjectGetClassName(arg1::PetscObject,arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PetscObjectGetClassName,petsc1),PetscErrorCode,(PetscObject,Ptr{Ptr{UInt8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectSetType(arg1::PetscObject,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscObjectSetType,petsc1),PetscErrorCode,(PetscObject,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectSetPrecision(arg0::Type{Float64},arg1::PetscObject,arg2::PetscPrecision)
    ccall((:PetscObjectSetPrecision,petsc1),PetscErrorCode,(PetscObject,PetscPrecision),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectGetType(arg1::PetscObject,arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PetscObjectGetType,petsc1),PetscErrorCode,(PetscObject,Ptr{Ptr{UInt8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectSetName(arg1::PetscObject,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscObjectSetName,petsc1),PetscErrorCode,(PetscObject,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectGetName(arg1::PetscObject,arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PetscObjectGetName,petsc1),PetscErrorCode,(PetscObject,Ptr{Ptr{UInt8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectSetTabLevel(arg0::Type{Float64},arg1::PetscObject,arg2::Integer)
    ccall((:PetscObjectSetTabLevel,petsc1),PetscErrorCode,(PetscObject,Int32),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectGetTabLevel(arg0::Type{Float64},arg1::PetscObject,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscObjectGetTabLevel,petsc1),PetscErrorCode,(PetscObject,Ptr{Int32}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectIncrementTabLevel(arg0::Type{Float64},arg1::PetscObject,arg2::PetscObject,arg3::Integer)
    ccall((:PetscObjectIncrementTabLevel,petsc1),PetscErrorCode,(PetscObject,PetscObject,Int32),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectReference(arg0::Type{Float64},arg1::PetscObject)
    ccall((:PetscObjectReference,petsc1),PetscErrorCode,(PetscObject,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectGetReference(arg0::Type{Float64},arg1::PetscObject,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscObjectGetReference,petsc1),PetscErrorCode,(PetscObject,Ptr{Int32}),arg1,arg2)
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
 function PetscObjectCompose(arg1::PetscObject,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::PetscObject)
    ccall((:PetscObjectCompose,petsc1),PetscErrorCode,(PetscObject,Cstring,PetscObject),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectRemoveReference(arg1::PetscObject,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscObjectRemoveReference,petsc1),PetscErrorCode,(PetscObject,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectQuery(arg1::PetscObject,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{PetscObject},StridedArray{PetscObject},Ptr{Void}))
    ccall((:PetscObjectQuery,petsc1),PetscErrorCode,(PetscObject,Cstring,Ptr{PetscObject}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectComposeFunction_Private(arg1::PetscObject,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
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
    ccall((:PetscCommGetNewTag,petsc1),PetscErrorCode,(comm_type,Ptr{PetscMPIInt}),arg1.val,arg2)
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
 function PetscObjectsListGetGlobalNumbering(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{PetscObject},StridedArray{PetscObject},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscObjectsListGetGlobalNumbering,petsc1),PetscErrorCode,(comm_type,Int32,Ptr{PetscObject},Ptr{Int32},Ptr{Int32}),arg1.val,arg2,arg3,arg4,arg5)
end 
=#
function PetscOptionsHasName(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsHasName,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{PetscBool}),arg1,arg2,arg3)
end

function PetscOptionsGetInt(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsGetInt,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Int32},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function PetscOptionsGetBool(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsGetBool,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function PetscOptionsGetReal(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsGetReal,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Cint},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function PetscOptionsGetScalar(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsGetScalar,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Float64},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function PetscOptionsGetIntArray(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsGetIntArray,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Int32},Ptr{Int32},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function PetscOptionsGetRealArray(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),Float64::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsGetRealArray,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Cint},Ptr{Int32},Ptr{PetscBool}),arg1,arg2,PetscReal,arg3,arg4)
end

function PetscOptionsGetScalarArray(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsGetScalarArray,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Float64},Ptr{Int32},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function PetscOptionsGetBoolArray(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsGetBoolArray,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{PetscBool},Ptr{Int32},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function PetscOptionsGetString(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsGetString,petsc1),PetscErrorCode,(Cstring,Cstring,Cstring,Cint,Ptr{PetscBool}),arg1,arg2,arg3,size_t,arg4)
end

function PetscOptionsGetStringArray(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsGetStringArray,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Ptr{UInt8}},Ptr{Int32},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function PetscOptionsGetEList(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsGetEList,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Ptr{UInt8}},Int32,Ptr{Int32},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscOptionsGetEnum(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}),arg4::Union(Ptr{PetscEnum},StridedArray{PetscEnum},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsGetEnum,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Ptr{UInt8}},Ptr{PetscEnum},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function PetscOptionsGetEnumArray(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}),arg4::Union(Ptr{PetscEnum},StridedArray{PetscEnum},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsGetEnumArray,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Ptr{UInt8}},Ptr{PetscEnum},Ptr{Int32},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscOptionsValidKey(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsValidKey,petsc1),PetscErrorCode,(Cstring,Ptr{PetscBool}),arg1,arg2)
end

function PetscOptionsSetAlias(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscOptionsSetAlias,petsc1),PetscErrorCode,(Cstring,Cstring),arg1,arg2)
end

function PetscOptionsSetValue(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscOptionsSetValue,petsc1),PetscErrorCode,(Cstring,Cstring),arg1,arg2)
end

function PetscOptionsClearValue(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscOptionsClearValue,petsc1),PetscErrorCode,(Cstring,),arg1)
end

function PetscOptionsAllUsed(arg0::Type{Float64},arg1::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscOptionsAllUsed,petsc1),PetscErrorCode,(Ptr{Int32},),arg1)
end

function PetscOptionsUsed(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsUsed,petsc1),PetscErrorCode,(Cstring,Ptr{PetscBool}),arg1,arg2)
end

function PetscOptionsLeft(arg0::Type{Float64})
    ccall((:PetscOptionsLeft,petsc1),PetscErrorCode,())
end

function PetscOptionsView(arg1::PetscViewer{Float64})
    ccall((:PetscOptionsView,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
end

function PetscOptionsCreate(arg0::Type{Float64})
    ccall((:PetscOptionsCreate,petsc1),PetscErrorCode,())
end

function PetscOptionsInsert(arg1::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg2::Union(Ptr{Ptr{ASCIIString}},StridedArray{Ptr{ASCIIString}},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscOptionsInsert,petsc1),PetscErrorCode,(Ptr{Cint},Ptr{Ptr{Ptr{UInt8}}},Cstring),arg1,arg2,arg3)
end

function PetscOptionsInsertFile(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::PetscBool)
    ccall((:PetscOptionsInsertFile,petsc1),PetscErrorCode,(comm_type,Cstring,PetscBool),arg1.val,arg2,arg3)
end

function PetscOptionsInsertString(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscOptionsInsertString,petsc1),PetscErrorCode,(Cstring,),arg1)
end

function PetscOptionsDestroy(arg0::Type{Float64})
    ccall((:PetscOptionsDestroy,petsc1),PetscErrorCode,())
end

function PetscOptionsClear(arg0::Type{Float64})
    ccall((:PetscOptionsClear,petsc1),PetscErrorCode,())
end

function PetscOptionsPrefixPush(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscOptionsPrefixPush,petsc1),PetscErrorCode,(Cstring,),arg1)
end

function PetscOptionsPrefixPop(arg0::Type{Float64})
    ccall((:PetscOptionsPrefixPop,petsc1),PetscErrorCode,())
end

function PetscOptionsReject(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscOptionsReject,petsc1),PetscErrorCode,(Cstring,Cstring),arg1,arg2)
end

function PetscOptionsGetAll(arg1::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PetscOptionsGetAll,petsc1),PetscErrorCode,(Ptr{Ptr{UInt8}},),arg1)
end

function PetscOptionsGetenv(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsGetenv,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint,Ptr{PetscBool}),arg1.val,arg2,arg3,size_t,arg4)
end

function PetscOptionsStringToInt(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscOptionsStringToInt,petsc1),PetscErrorCode,(Cstring,Ptr{Int32}),arg1,arg2)
end

function PetscOptionsStringToReal(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscOptionsStringToReal,petsc1),PetscErrorCode,(Cstring,Ptr{Cint}),arg1,arg2)
end

function PetscOptionsStringToBool(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsStringToBool,petsc1),PetscErrorCode,(Cstring,Ptr{PetscBool}),arg1,arg2)
end

function PetscOptionsMonitorSet(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscOptionsMonitorSet,petsc1),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function PetscOptionsMonitorCancel(arg0::Type{Float64})
    ccall((:PetscOptionsMonitorCancel,petsc1),PetscErrorCode,())
end

function PetscOptionsMonitorDefault(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscOptionsMonitorDefault,petsc1),PetscErrorCode,(Cstring,Cstring,Ptr{Void}),arg1,arg2,arg3)
end

#= skipping function with undefined symbols: 
 function PetscOptionsBegin_Private(arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::MPI_Comm,arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
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
 function PetscOptionsHead(arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscOptionsHead,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsEnum_Private(arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}),arg6::PetscEnum,arg7::Union(Ptr{PetscEnum},StridedArray{PetscEnum},Ptr{Void}),arg8::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsEnum_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{Ptr{UInt8}},PetscEnum,Ptr{PetscEnum},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsInt_Private(arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::Integer,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsInt_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Int32,Ptr{Int32},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsReal_Private(arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),Float64::Integer,arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg6::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsReal_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Cint,Ptr{Cint},Ptr{PetscBool}),arg1,arg2,arg3,arg4,PetscReal,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsScalar_Private(arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::Float64,arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsScalar_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Float64,Ptr{Float64},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsName_Private(arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsName_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsString_Private(arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg6::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer,arg7::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsString_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Cstring,Cstring,Cint,Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,size_t,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsBool_Private(arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::PetscBool,arg6::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg7::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsBool_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,PetscBool,Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsBoolGroupBegin_Private(arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsBoolGroupBegin_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsBoolGroup_Private(arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsBoolGroup_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsBoolGroupEnd_Private(arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsBoolGroupEnd_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsFList_Private(arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::PetscFunctionList,arg6::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg7::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer,arg8::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsFList_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,PetscFunctionList,Cstring,Cstring,Cint,Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,size_t,arg8)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsEList_Private(arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}),arg6::Integer,arg7::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg8::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg9::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsEList_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{Ptr{UInt8}},Int32,Cstring,Ptr{Int32},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsRealArray_Private(arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),Float64::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsRealArray_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{Cint},Ptr{Int32},Ptr{PetscBool}),arg1,arg2,arg3,arg4,PetscReal,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsScalarArray_Private(arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsScalarArray_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{Float64},Ptr{Int32},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsIntArray_Private(arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsIntArray_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{Int32},Ptr{Int32},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsStringArray_Private(arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}),arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsStringArray_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{Ptr{UInt8}},Ptr{Int32},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsBoolArray_Private(arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsBoolArray_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{PetscBool},Ptr{Int32},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscOptionsEnumArray_Private(arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}),arg6::Union(Ptr{PetscEnum},StridedArray{PetscEnum},Ptr{Void}),arg7::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg8::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsEnumArray_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{Ptr{UInt8}},Ptr{PetscEnum},Ptr{Int32},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end 
=#
function PetscOptionsSetFromOptions(arg0::Type{Float64})
    ccall((:PetscOptionsSetFromOptions,petsc1),PetscErrorCode,())
end

function PetscOptionsSAWsDestroy(arg0::Type{Float64})
    ccall((:PetscOptionsSAWsDestroy,petsc1),PetscErrorCode,())
end

function PetscMemoryShowUsage(arg1::PetscViewer{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscMemoryShowUsage,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring),arg1,arg2)
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
 function PetscObjectQueryFunction_Private(arg1::PetscObject,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:PetscObjectQueryFunction_Private,petsc1),PetscErrorCode,(PetscObject,Cstring,Ptr{Ptr{Void}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectSetOptionsPrefix(arg1::PetscObject,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscObjectSetOptionsPrefix,petsc1),PetscErrorCode,(PetscObject,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectAppendOptionsPrefix(arg1::PetscObject,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscObjectAppendOptionsPrefix,petsc1),PetscErrorCode,(PetscObject,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectPrependOptionsPrefix(arg1::PetscObject,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscObjectPrependOptionsPrefix,petsc1),PetscErrorCode,(PetscObject,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectGetOptionsPrefix(arg1::PetscObject,arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PetscObjectGetOptionsPrefix,petsc1),PetscErrorCode,(PetscObject,Ptr{Ptr{UInt8}}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectChangeTypeName(arg1::PetscObject,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscObjectChangeTypeName,petsc1),PetscErrorCode,(PetscObject,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectRegisterDestroy(arg0::Type{Float64},arg1::PetscObject)
    ccall((:PetscObjectRegisterDestroy,petsc1),PetscErrorCode,(PetscObject,),arg1)
end 
=#
function PetscObjectRegisterDestroyAll(arg0::Type{Float64})
    ccall((:PetscObjectRegisterDestroyAll,petsc1),PetscErrorCode,())
end

#= skipping function with undefined symbols: 
 function PetscObjectViewFromOptions(arg1::PetscObject,arg2::PetscObject,arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscObjectViewFromOptions,petsc1),PetscErrorCode,(PetscObject,PetscObject,Cstring),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectName(arg0::Type{Float64},arg1::PetscObject)
    ccall((:PetscObjectName,petsc1),PetscErrorCode,(PetscObject,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectTypeCompare(arg1::PetscObject,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscObjectTypeCompare,petsc1),PetscErrorCode,(PetscObject,Cstring,Ptr{PetscBool}),arg1,arg2,arg3)
end 
=#
function PetscRegisterFinalize(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscRegisterFinalize,petsc1),PetscErrorCode,(Ptr{Void},),arg1)
end

function PetscRegisterFinalizeAll(arg0::Type{Float64})
    ccall((:PetscRegisterFinalizeAll,petsc1),PetscErrorCode,())
end

function PetscDLOpen(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::PetscDLMode,arg3::Union(Ptr{PetscDLHandle},StridedArray{PetscDLHandle},Ptr{Void}))
    ccall((:PetscDLOpen,petsc1),PetscErrorCode,(Cstring,PetscDLMode,Ptr{PetscDLHandle}),arg1,arg2,arg3)
end

function PetscDLClose(arg0::Type{Float64},arg1::Union(Ptr{PetscDLHandle},StridedArray{PetscDLHandle},Ptr{Void}))
    ccall((:PetscDLClose,petsc1),PetscErrorCode,(Ptr{PetscDLHandle},),arg1)
end

function PetscDLSym(arg1::PetscDLHandle,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:PetscDLSym,petsc1),PetscErrorCode,(PetscDLHandle,Cstring,Ptr{Ptr{Void}}),arg1,arg2,arg3)
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
 function PetscObjectListFind(arg1::PetscObjectList,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{PetscObject},StridedArray{PetscObject},Ptr{Void}))
    ccall((:PetscObjectListFind,petsc1),PetscErrorCode,(PetscObjectList,Cstring,Ptr{PetscObject}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectListReverseFind(arg1::PetscObjectList,arg2::PetscObject,arg3::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}),arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscObjectListReverseFind,petsc1),PetscErrorCode,(PetscObjectList,PetscObject,Ptr{Ptr{UInt8}},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectListAdd(arg1::Union(Ptr{PetscObjectList},StridedArray{PetscObjectList},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::PetscObject)
    ccall((:PetscObjectListAdd,petsc1),PetscErrorCode,(Ptr{PetscObjectList},Cstring,PetscObject),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectListRemoveReference(arg1::Union(Ptr{PetscObjectList},StridedArray{PetscObjectList},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscObjectListRemoveReference,petsc1),PetscErrorCode,(Ptr{PetscObjectList},Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscObjectListDuplicate(arg0::Type{Float64},arg1::PetscObjectList,arg2::Union(Ptr{PetscObjectList},StridedArray{PetscObjectList},Ptr{Void}))
    ccall((:PetscObjectListDuplicate,petsc1),PetscErrorCode,(PetscObjectList,Ptr{PetscObjectList}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFunctionListAdd_Private(arg1::Union(Ptr{PetscFunctionList},StridedArray{PetscFunctionList},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscFunctionListAdd_Private,petsc1),PetscErrorCode,(Ptr{PetscFunctionList},Cstring,Ptr{Void}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFunctionListDestroy(arg0::Type{Float64},arg1::Union(Ptr{PetscFunctionList},StridedArray{PetscFunctionList},Ptr{Void}))
    ccall((:PetscFunctionListDestroy,petsc1),PetscErrorCode,(Ptr{PetscFunctionList},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFunctionListFind_Private(arg1::PetscFunctionList,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:PetscFunctionListFind_Private,petsc1),PetscErrorCode,(PetscFunctionList,Cstring,Ptr{Ptr{Void}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFunctionListPrintTypes(arg1::MPI_Comm,arg2::Union(Ptr{FILE},StridedArray{FILE},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg6::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg7::PetscFunctionList,arg8::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
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
 function PetscFunctionListGet(arg1::PetscFunctionList,arg2::Union(Ptr{Ptr{ASCIIString}},StridedArray{Ptr{ASCIIString}},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscFunctionListGet,petsc1),PetscErrorCode,(PetscFunctionList,Ptr{Ptr{Ptr{UInt8}}},Ptr{Cint}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDLLibraryAppend(arg1::MPI_Comm,arg2::Union(Ptr{PetscDLLibrary},StridedArray{PetscDLLibrary},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscDLLibraryAppend,petsc1),PetscErrorCode,(comm_type,Ptr{PetscDLLibrary},Cstring),arg1.val,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDLLibraryPrepend(arg1::MPI_Comm,arg2::Union(Ptr{PetscDLLibrary},StridedArray{PetscDLLibrary},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscDLLibraryPrepend,petsc1),PetscErrorCode,(comm_type,Ptr{PetscDLLibrary},Cstring),arg1.val,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDLLibrarySym(arg1::MPI_Comm,arg2::Union(Ptr{PetscDLLibrary},StridedArray{PetscDLLibrary},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:PetscDLLibrarySym,petsc1),PetscErrorCode,(comm_type,Ptr{PetscDLLibrary},Cstring,Cstring,Ptr{Ptr{Void}}),arg1.val,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDLLibraryPrintPath(arg0::Type{Float64},arg1::PetscDLLibrary)
    ccall((:PetscDLLibraryPrintPath,petsc1),PetscErrorCode,(PetscDLLibrary,),arg1)
end 
=#
function PetscDLLibraryRetrieve(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscDLLibraryRetrieve,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint,Ptr{PetscBool}),arg1.val,arg2,arg3,size_t,arg4)
end

#= skipping function with undefined symbols: 
 function PetscDLLibraryOpen(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{PetscDLLibrary},StridedArray{PetscDLLibrary},Ptr{Void}))
    ccall((:PetscDLLibraryOpen,petsc1),PetscErrorCode,(comm_type,Cstring,Ptr{PetscDLLibrary}),arg1.val,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscDLLibraryClose(arg0::Type{Float64},arg1::PetscDLLibrary)
    ccall((:PetscDLLibraryClose,petsc1),PetscErrorCode,(PetscDLLibrary,),arg1)
end 
=#
function PetscSplitOwnership(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSplitOwnership,petsc1),PetscErrorCode,(comm_type,Ptr{Int32},Ptr{Int32}),arg1.val,arg2,arg3)
end

function PetscSplitOwnershipBlock(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSplitOwnershipBlock,petsc1),PetscErrorCode,(comm_type,Int32,Ptr{Int32},Ptr{Int32}),arg1.val,arg2,arg3,arg4)
end

function PetscSequentialPhaseBegin(arg0::Type{Float64},arg1::MPI_Comm,arg2::PetscMPIInt)
    ccall((:PetscSequentialPhaseBegin,petsc1),PetscErrorCode,(comm_type,PetscMPIInt),arg1.val,arg2)
end

function PetscSequentialPhaseEnd(arg0::Type{Float64},arg1::MPI_Comm,arg2::PetscMPIInt)
    ccall((:PetscSequentialPhaseEnd,petsc1),PetscErrorCode,(comm_type,PetscMPIInt),arg1.val,arg2)
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
    ccall((:PetscInfoDeactivateClass,petsc1),PetscErrorCode,(PetscClassId,),arg1)
end

function PetscInfoActivateClass(arg0::Type{Float64},arg1::PetscClassId)
    ccall((:PetscInfoActivateClass,petsc1),PetscErrorCode,(PetscClassId,),arg1)
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
function PetscFixFilename(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscFixFilename,petsc1),PetscErrorCode,(Cstring,Cstring),arg1,arg2)
end

#= skipping function with undefined symbols: 
 function PetscFOpen(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(Ptr{Ptr{FILE}},StridedArray{Ptr{FILE}},Ptr{Void}))
    ccall((:PetscFOpen,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Ptr{Ptr{FILE}}),arg1.val,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscFClose(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{FILE},StridedArray{FILE},Ptr{Void}))
    ccall((:PetscFClose,petsc1),PetscErrorCode,(comm_type,Ptr{FILE}),arg1.val,arg2)
end 
=#
function PetscVSNPrintf(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),va_list::Integer)
    ccall((:PetscVSNPrintf,petsc1),PetscErrorCode,(Cstring,Cint,Cstring,Ptr{Cint},Cint),arg1,size_t,arg2,arg3,va_list)
end

#= skipping function with undefined symbols: 
 function PetscVFPrintfDefault(arg1::Union(Ptr{FILE},StridedArray{FILE},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),va_list::Integer)
    ccall((:PetscVFPrintfDefault,petsc1),PetscErrorCode,(Ptr{FILE},Cstring,Cint),arg1,arg2,va_list)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSynchronizedFlush(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{FILE},StridedArray{FILE},Ptr{Void}))
    ccall((:PetscSynchronizedFlush,petsc1),PetscErrorCode,(comm_type,Ptr{FILE}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSynchronizedFGets(arg1::MPI_Comm,arg2::Union(Ptr{FILE},StridedArray{FILE},Ptr{Void}),size_t::Integer,arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscSynchronizedFGets,petsc1),PetscErrorCode,(comm_type,Ptr{FILE},Cint,Cstring),arg1.val,arg2,size_t,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscStartMatlab(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(Ptr{Ptr{FILE}},StridedArray{Ptr{FILE}},Ptr{Void}))
    ccall((:PetscStartMatlab,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Ptr{Ptr{FILE}}),arg1.val,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscStartJava(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(Ptr{Ptr{FILE}},StridedArray{Ptr{FILE}},Ptr{Void}))
    ccall((:PetscStartJava,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Ptr{Ptr{FILE}}),arg1.val,arg2,arg3,arg4)
end 
=#
function PetscGetPetscDir(arg1::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PetscGetPetscDir,petsc1),PetscErrorCode,(Ptr{Ptr{UInt8}},),arg1)
end

function PetscPopUpSelect(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}),arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscPopUpSelect,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint,Ptr{Ptr{UInt8}},Ptr{Cint}),arg1.val,arg2,arg3,arg4,arg5,arg6)
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
function PetscIntView(arg1::Integer,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::PetscViewer{Float64})
    ccall((:PetscIntView,petsc1),PetscErrorCode,(Int32,Ptr{Int32},PetscViewer{Float64}),arg1,arg2,arg3)
end

function PetscRealView(arg1::Integer,Float64::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg2::PetscViewer{Float64})
    ccall((:PetscRealView,petsc1),PetscErrorCode,(Int32,Ptr{Cint},PetscViewer{Float64}),arg1,PetscReal,arg2)
end

function PetscScalarView(arg1::Integer,arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg3::PetscViewer{Float64})
    ccall((:PetscScalarView,petsc1),PetscErrorCode,(Int32,Ptr{Float64},PetscViewer{Float64}),arg1,arg2,arg3)
end

function PetscGetHostName(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer)
    ccall((:PetscGetHostName,petsc1),PetscErrorCode,(Cstring,Cint),arg1,size_t)
end

function PetscGetUserName(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer)
    ccall((:PetscGetUserName,petsc1),PetscErrorCode,(Cstring,Cint),arg1,size_t)
end

function PetscGetProgramName(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer)
    ccall((:PetscGetProgramName,petsc1),PetscErrorCode,(Cstring,Cint),arg1,size_t)
end

function PetscSetProgramName(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscSetProgramName,petsc1),PetscErrorCode,(Cstring,),arg1)
end

function PetscGetDate(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer)
    ccall((:PetscGetDate,petsc1),PetscErrorCode,(Cstring,Cint),arg1,size_t)
end

function PetscGetVersion(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer)
    ccall((:PetscGetVersion,petsc1),PetscErrorCode,(Cstring,Cint),arg1,size_t)
end

function PetscSortInt(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSortInt,petsc1),PetscErrorCode,(Int32,Ptr{Int32}),arg1,arg2)
end

function PetscSortRemoveDupsInt(arg0::Type{Float64},arg1::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSortRemoveDupsInt,petsc1),PetscErrorCode,(Ptr{Int32},Ptr{Int32}),arg1,arg2)
end

function PetscFindInt(arg0::Type{Float64},arg1::Integer,arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscFindInt,petsc1),PetscErrorCode,(Int32,Int32,Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3,arg4)
end

function PetscSortIntWithPermutation(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSortIntWithPermutation,petsc1),PetscErrorCode,(Int32,Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3)
end

function PetscSortStrWithPermutation(arg1::Integer,arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSortStrWithPermutation,petsc1),PetscErrorCode,(Int32,Ptr{Ptr{UInt8}},Ptr{Int32}),arg1,arg2,arg3)
end

function PetscSortIntWithArray(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSortIntWithArray,petsc1),PetscErrorCode,(Int32,Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3)
end

function PetscSortIntWithArrayPair(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSortIntWithArrayPair,petsc1),PetscErrorCode,(Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3,arg4)
end

function PetscSortMPIInt(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}))
    ccall((:PetscSortMPIInt,petsc1),PetscErrorCode,(Int32,Ptr{PetscMPIInt}),arg1,arg2)
end

function PetscSortRemoveDupsMPIInt(arg0::Type{Float64},arg1::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg2::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}))
    ccall((:PetscSortRemoveDupsMPIInt,petsc1),PetscErrorCode,(Ptr{Int32},Ptr{PetscMPIInt}),arg1,arg2)
end

function PetscSortMPIIntWithArray(arg0::Type{Float64},arg1::PetscMPIInt,arg2::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}),arg3::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}))
    ccall((:PetscSortMPIIntWithArray,petsc1),PetscErrorCode,(PetscMPIInt,Ptr{PetscMPIInt},Ptr{PetscMPIInt}),arg1,arg2,arg3)
end

function PetscSortIntWithScalarArray(arg1::Integer,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:PetscSortIntWithScalarArray,petsc1),PetscErrorCode,(Int32,Ptr{Int32},Ptr{Float64}),arg1,arg2,arg3)
end

function PetscSortIntWithDataArray(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),size_t::Integer,arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscSortIntWithDataArray,petsc1),PetscErrorCode,(Int32,Ptr{Int32},Ptr{Void},Cint,Ptr{Void}),arg1,arg2,arg3,size_t,arg4)
end

function PetscSortReal(arg0::Type{Float64},arg1::Integer,Float64::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscSortReal,petsc1),PetscErrorCode,(Int32,Ptr{Cint}),arg1,PetscReal)
end

function PetscSortRealWithPermutation(arg0::Type{Float64},arg1::Integer,Float64::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSortRealWithPermutation,petsc1),PetscErrorCode,(Int32,Ptr{Cint},Ptr{Int32}),arg1,PetscReal,arg2)
end

function PetscSortRemoveDupsReal(arg0::Type{Float64},arg1::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),Float64::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscSortRemoveDupsReal,petsc1),PetscErrorCode,(Ptr{Int32},Ptr{Cint}),arg1,PetscReal)
end

function PetscSortSplit(arg1::Integer,arg2::Integer,arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSortSplit,petsc1),PetscErrorCode,(Int32,Int32,Ptr{Float64},Ptr{Int32}),arg1,arg2,arg3,arg4)
end

function PetscSortSplitReal(arg0::Type{Float64},arg1::Integer,arg2::Integer,Float64::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSortSplitReal,petsc1),PetscErrorCode,(Int32,Int32,Ptr{Cint},Ptr{Int32}),arg1,arg2,PetscReal,arg3)
end

function PetscProcessTree(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg6::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg7::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg8::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:PetscProcessTree,petsc1),PetscErrorCode,(Int32,Ptr{PetscBool},Ptr{Int32},Ptr{Int32},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscMergeIntArrayPair(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg8::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg9::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:PetscMergeIntArrayPair,petsc1),PetscErrorCode,(Int32,Ptr{Int32},Ptr{Int32},Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function PetscMergeIntArray(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:PetscMergeIntArray,petsc1),PetscErrorCode,(Int32,Ptr{Int32},Int32,Ptr{Int32},Ptr{Int32},Ptr{Ptr{Int32}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscSetDisplay(arg0::Type{Float64})
    ccall((:PetscSetDisplay,petsc1),PetscErrorCode,())
end

function PetscGetDisplay(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer)
    ccall((:PetscGetDisplay,petsc1),PetscErrorCode,(Cstring,Cint),arg1,size_t)
end

function PetscRandomInitializePackage(arg0::Type{Float64})
    ccall((:PetscRandomInitializePackage,petsc1),PetscErrorCode,())
end

function PetscRandomRegister(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscRandomRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
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
    ccall((:PetscRandomGetType,petsc1),PetscErrorCode,(PetscRandom,Ptr{Cstring}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscRandomCreate(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscRandom},StridedArray{PetscRandom},Ptr{Void}))
    ccall((:PetscRandomCreate,petsc1),PetscErrorCode,(comm_type,Ptr{PetscRandom}),arg1.val,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscRandomGetValue(arg1::PetscRandom,arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:PetscRandomGetValue,petsc1),PetscErrorCode,(PetscRandom,Ptr{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscRandomGetValueReal(arg0::Type{Float64},arg1::PetscRandom,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscRandomGetValueReal,petsc1),PetscErrorCode,(PetscRandom,Ptr{Cint}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscRandomGetInterval(arg1::PetscRandom,arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:PetscRandomGetInterval,petsc1),PetscErrorCode,(PetscRandom,Ptr{Float64},Ptr{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscRandomSetInterval(arg1::PetscRandom,arg2::Float64,arg3::Float64)
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
function PetscGetFullPath(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer)
    ccall((:PetscGetFullPath,petsc1),PetscErrorCode,(Cstring,Cstring,Cint),arg1,arg2,size_t)
end

function PetscGetRelativePath(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer)
    ccall((:PetscGetRelativePath,petsc1),PetscErrorCode,(Cstring,Cstring,Cint),arg1,arg2,size_t)
end

function PetscGetWorkingDirectory(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer)
    ccall((:PetscGetWorkingDirectory,petsc1),PetscErrorCode,(Cstring,Cint),arg1,size_t)
end

function PetscGetRealPath(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscGetRealPath,petsc1),PetscErrorCode,(Cstring,Cstring),arg1,arg2)
end

function PetscGetHomeDirectory(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer)
    ccall((:PetscGetHomeDirectory,petsc1),PetscErrorCode,(Cstring,Cint),arg1,size_t)
end

function PetscTestFile(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Uint8,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscTestFile,petsc1),PetscErrorCode,(Cstring,Uint8,Ptr{PetscBool}),arg1,arg2,arg3)
end

function PetscTestDirectory(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Uint8,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscTestDirectory,petsc1),PetscErrorCode,(Cstring,Uint8,Ptr{PetscBool}),arg1,arg2,arg3)
end

function PetscBinaryRead(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Integer,arg4::PetscDataType)
    ccall((:PetscBinaryRead,petsc1),PetscErrorCode,(Cint,Ptr{Void},Int32,PetscDataType),arg1,arg2,arg3,arg4)
end

function PetscBinarySynchronizedRead(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Integer,arg5::PetscDataType)
    ccall((:PetscBinarySynchronizedRead,petsc1),PetscErrorCode,(comm_type,Cint,Ptr{Void},Int32,PetscDataType),arg1.val,arg2,arg3,arg4,arg5)
end

function PetscBinarySynchronizedWrite(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Integer,arg5::PetscDataType,arg6::PetscBool)
    ccall((:PetscBinarySynchronizedWrite,petsc1),PetscErrorCode,(comm_type,Cint,Ptr{Void},Int32,PetscDataType,PetscBool),arg1.val,arg2,arg3,arg4,arg5,arg6)
end

function PetscBinaryWrite(arg0::Type{Float64},arg1::Integer,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Integer,arg4::PetscDataType,arg5::PetscBool)
    ccall((:PetscBinaryWrite,petsc1),PetscErrorCode,(Cint,Ptr{Void},Int32,PetscDataType,PetscBool),arg1,arg2,arg3,arg4,arg5)
end

function PetscBinaryOpen(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::PetscFileMode,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscBinaryOpen,petsc1),PetscErrorCode,(Cstring,PetscFileMode,Ptr{Cint}),arg1,arg2,arg3)
end

function PetscBinaryClose(arg0::Type{Float64},arg1::Integer)
    ccall((:PetscBinaryClose,petsc1),PetscErrorCode,(Cint,),arg1)
end

function PetscSharedTmp(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscSharedTmp,petsc1),PetscErrorCode,(comm_type,Ptr{PetscBool}),arg1.val,arg2)
end

function PetscSharedWorkingDirectory(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscSharedWorkingDirectory,petsc1),PetscErrorCode,(comm_type,Ptr{PetscBool}),arg1.val,arg2)
end

function PetscGetTmp(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer)
    ccall((:PetscGetTmp,petsc1),PetscErrorCode,(comm_type,Cstring,Cint),arg1.val,arg2,size_t)
end

function PetscFileRetrieve(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscFileRetrieve,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint,Ptr{PetscBool}),arg1.val,arg2,arg3,size_t,arg4)
end

function PetscLs(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscLs,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint,Ptr{PetscBool}),arg1.val,arg2,arg3,size_t,arg4)
end

function PetscOpenSocket(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Integer,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscOpenSocket,petsc1),PetscErrorCode,(Cstring,Cint,Ptr{Cint}),arg1,arg2,arg3)
end

function PetscBinarySeek(arg0::Type{Float64},arg1::Integer,off_t::Integer,arg2::PetscBinarySeekType,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscBinarySeek,petsc1),PetscErrorCode,(Cint,Cint,PetscBinarySeekType,Ptr{Cint}),arg1,off_t,arg2,arg3)
end

function PetscBinarySynchronizedSeek(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,off_t::Integer,arg3::PetscBinarySeekType,arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscBinarySynchronizedSeek,petsc1),PetscErrorCode,(comm_type,Cint,Cint,PetscBinarySeekType,Ptr{Cint}),arg1.val,arg2,off_t,arg3,arg4)
end

function PetscByteSwap(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::PetscDataType,arg3::Integer)
    ccall((:PetscByteSwap,petsc1),PetscErrorCode,(Ptr{Void},PetscDataType,Int32),arg1,arg2,arg3)
end

function PetscSetDebugTerminal(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscSetDebugTerminal,petsc1),PetscErrorCode,(Cstring,),arg1)
end

function PetscSetDebugger(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::PetscBool)
    ccall((:PetscSetDebugger,petsc1),PetscErrorCode,(Cstring,PetscBool),arg1,arg2)
end

function PetscSetDefaultDebugger(arg0::Type{Float64})
    ccall((:PetscSetDefaultDebugger,petsc1),PetscErrorCode,())
end

function PetscSetDebuggerFromString(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscSetDebuggerFromString,petsc1),PetscErrorCode,(Cstring,),arg1)
end

function PetscAttachDebugger(arg0::Type{Float64})
    ccall((:PetscAttachDebugger,petsc1),PetscErrorCode,())
end

function PetscStopForDebugger(arg0::Type{Float64})
    ccall((:PetscStopForDebugger,petsc1),PetscErrorCode,())
end

function PetscGatherNumberOfMessages(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}),arg3::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}),arg4::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}))
    ccall((:PetscGatherNumberOfMessages,petsc1),PetscErrorCode,(comm_type,Ptr{PetscMPIInt},Ptr{PetscMPIInt},Ptr{PetscMPIInt}),arg1.val,arg2,arg3,arg4)
end

function PetscGatherMessageLengths(arg0::Type{Float64},arg1::MPI_Comm,arg2::PetscMPIInt,arg3::PetscMPIInt,arg4::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}),arg5::Union(Ptr{Ptr{PetscMPIInt}},StridedArray{Ptr{PetscMPIInt}},Ptr{Void}),arg6::Union(Ptr{Ptr{PetscMPIInt}},StridedArray{Ptr{PetscMPIInt}},Ptr{Void}))
    ccall((:PetscGatherMessageLengths,petsc1),PetscErrorCode,(comm_type,PetscMPIInt,PetscMPIInt,Ptr{PetscMPIInt},Ptr{Ptr{PetscMPIInt}},Ptr{Ptr{PetscMPIInt}}),arg1.val,arg2,arg3,arg4,arg5,arg6)
end

function PetscGatherMessageLengths2(arg0::Type{Float64},arg1::MPI_Comm,arg2::PetscMPIInt,arg3::PetscMPIInt,arg4::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}),arg5::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}),arg6::Union(Ptr{Ptr{PetscMPIInt}},StridedArray{Ptr{PetscMPIInt}},Ptr{Void}),arg7::Union(Ptr{Ptr{PetscMPIInt}},StridedArray{Ptr{PetscMPIInt}},Ptr{Void}),arg8::Union(Ptr{Ptr{PetscMPIInt}},StridedArray{Ptr{PetscMPIInt}},Ptr{Void}))
    ccall((:PetscGatherMessageLengths2,petsc1),PetscErrorCode,(comm_type,PetscMPIInt,PetscMPIInt,Ptr{PetscMPIInt},Ptr{PetscMPIInt},Ptr{Ptr{PetscMPIInt}},Ptr{Ptr{PetscMPIInt}},Ptr{Ptr{PetscMPIInt}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

#= skipping function with undefined symbols: 
 function PetscPostIrecvInt(arg0::Type{Float64},arg1::MPI_Comm,arg2::PetscMPIInt,arg3::PetscMPIInt,arg4::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}),arg5::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}),arg6::Union(Ptr{Ptr{Ptr{Int32}}},StridedArray{Ptr{Ptr{Int32}}},Ptr{Void}),arg7::Union(Ptr{Ptr{MPI_Request}},StridedArray{Ptr{MPI_Request}},Ptr{Void}))
    ccall((:PetscPostIrecvInt,petsc1),PetscErrorCode,(comm_type,PetscMPIInt,PetscMPIInt,Ptr{PetscMPIInt},Ptr{PetscMPIInt},Ptr{Ptr{Ptr{Int32}}},Ptr{Ptr{MPI_Request}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscPostIrecvScalar(arg1::MPI_Comm,arg2::PetscMPIInt,arg3::PetscMPIInt,arg4::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}),arg5::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}),arg6::Union(Ptr{Ptr{Ptr{Float64}}},StridedArray{Ptr{Ptr{Float64}}},Ptr{Void}),arg7::Union(Ptr{Ptr{MPI_Request}},StridedArray{Ptr{MPI_Request}},Ptr{Void}))
    ccall((:PetscPostIrecvScalar,petsc1),PetscErrorCode,(comm_type,PetscMPIInt,PetscMPIInt,Ptr{PetscMPIInt},Ptr{PetscMPIInt},Ptr{Ptr{Ptr{Float64}}},Ptr{Ptr{MPI_Request}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscCommBuildTwoSided(arg0::Type{Float64},arg1::MPI_Comm,arg2::PetscMPIInt,arg3::MPI_Datatype,arg4::Integer,arg5::Union(Ptr{PetscMPIInt},StridedArray{PetscMPIInt},Ptr{Void}),arg6::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg7::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg8::Union(Ptr{Ptr{PetscMPIInt}},StridedArray{Ptr{PetscMPIInt}},Ptr{Void}),arg9::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscCommBuildTwoSided,petsc1),PetscErrorCode,(comm_type,PetscMPIInt,MPI_Datatype,Int32,Ptr{PetscMPIInt},Ptr{Void},Ptr{Int32},Ptr{Ptr{PetscMPIInt}},Ptr{Void}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end 
=#
function PetscCommBuildTwoSidedSetType(arg0::Type{Float64},arg1::MPI_Comm,arg2::PetscBuildTwoSidedType)
    ccall((:PetscCommBuildTwoSidedSetType,petsc1),PetscErrorCode,(comm_type,PetscBuildTwoSidedType),arg1.val,arg2)
end

function PetscCommBuildTwoSidedGetType(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscBuildTwoSidedType},StridedArray{PetscBuildTwoSidedType},Ptr{Void}))
    ccall((:PetscCommBuildTwoSidedGetType,petsc1),PetscErrorCode,(comm_type,Ptr{PetscBuildTwoSidedType}),arg1.val,arg2)
end

function PetscSSEIsEnabled(arg0::Type{Float64},arg1::MPI_Comm,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscSSEIsEnabled,petsc1),PetscErrorCode,(comm_type,Ptr{PetscBool},Ptr{PetscBool}),arg1.val,arg2,arg3)
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
    ccall((:PetscSubcommSetNumber,petsc1),PetscErrorCode,(PetscSubcomm,Int32),arg1,arg2)
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
    ccall((:PetscSegBufferCreate,petsc1),PetscErrorCode,())
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
function PetscGoogleDriveAuthorize(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer)
    ccall((:PetscGoogleDriveAuthorize,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint),arg1.val,arg2,arg3,size_t)
end

function PetscGoogleDriveRefresh(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer)
    ccall((:PetscGoogleDriveRefresh,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint),arg1.val,arg2,arg3,size_t)
end

function PetscGoogleDriveUpload(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscGoogleDriveUpload,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring),arg1.val,arg2,arg3)
end

function PetscBoxAuthorize(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer)
    ccall((:PetscBoxAuthorize,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint),arg1.val,arg2,arg3,size_t)
end

function PetscBoxRefresh(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer)
    ccall((:PetscBoxRefresh,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cstring,Cint),arg1.val,arg2,arg3,arg4,size_t)
end

function PetscTextBelt(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscTextBelt,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Ptr{PetscBool}),arg1.val,arg2,arg3,arg4)
end

function PetscPullJSONValue(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscPullJSONValue,petsc1),PetscErrorCode,(Cstring,Cstring,Cstring,Cint,Ptr{PetscBool}),arg1,arg2,arg3,size_t,arg4)
end

function PetscPushJSONValue(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer)
    ccall((:PetscPushJSONValue,petsc1),PetscErrorCode,(Cstring,Cstring,Cstring,Cint),arg1,arg2,arg3,size_t)
end

function ISInitializePackage(arg0::Type{Float64})
    ccall((:ISInitializePackage,petsc1),PetscErrorCode,())
end

function ISSetType(arg1::IS{Float64},arg2::ISType)
    ccall((:ISSetType,petsc1),PetscErrorCode,(IS{Float64},Cstring),arg1,arg2)
end

function ISGetType(arg1::IS{Float64},arg2::Union(Ptr{ISType},StridedArray{ISType},Ptr{Void}))
    ccall((:ISGetType,petsc1),PetscErrorCode,(IS{Float64},Ptr{Ptr{UInt8}}),arg1,arg2)
end

function ISRegister(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:ISRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
end

function ISCreate(arg1::MPI_Comm,arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISCreate,petsc1),PetscErrorCode,(comm_type,Ptr{IS{Float64}}),arg1.val,arg2)
end

function ISCreateGeneral(arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),PetscCopyMode::Integer,arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISCreateGeneral,petsc1),PetscErrorCode,(comm_type,Int32,Ptr{Int32},Cint,Ptr{IS{Float64}}),arg1.val,arg2,arg3,PetscCopyMode,arg4)
end

function ISGeneralSetIndices(arg1::IS{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),PetscCopyMode::Integer)
    ccall((:ISGeneralSetIndices,petsc1),PetscErrorCode,(IS{Float64},Int32,Ptr{Int32},Cint),arg1,arg2,arg3,PetscCopyMode)
end

function ISCreateBlock(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),PetscCopyMode::Integer,arg5::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISCreateBlock,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Ptr{Int32},Cint,Ptr{IS{Float64}}),arg1.val,arg2,arg3,arg4,PetscCopyMode,arg5)
end

function ISBlockSetIndices(arg1::IS{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),PetscCopyMode::Integer)
    ccall((:ISBlockSetIndices,petsc1),PetscErrorCode,(IS{Float64},Int32,Int32,Ptr{Int32},Cint),arg1,arg2,arg3,arg4,PetscCopyMode)
end

function ISCreateStride(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISCreateStride,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Ptr{IS{Float64}}),arg1.val,arg2,arg3,arg4,arg5)
end

function ISStrideSetStride(arg1::IS{Float64},arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:ISStrideSetStride,petsc1),PetscErrorCode,(IS{Float64},Int32,Int32,Int32),arg1,arg2,arg3,arg4)
end

function ISDestroy(arg1::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISDestroy,petsc1),PetscErrorCode,(Ptr{IS{Float64}},),arg1)
end

function ISSetPermutation(arg1::IS{Float64})
    ccall((:ISSetPermutation,petsc1),PetscErrorCode,(IS{Float64},),arg1)
end

function ISPermutation(arg1::IS{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:ISPermutation,petsc1),PetscErrorCode,(IS{Float64},Ptr{PetscBool}),arg1,arg2)
end

function ISSetIdentity(arg1::IS{Float64})
    ccall((:ISSetIdentity,petsc1),PetscErrorCode,(IS{Float64},),arg1)
end

function ISIdentity(arg1::IS{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:ISIdentity,petsc1),PetscErrorCode,(IS{Float64},Ptr{PetscBool}),arg1,arg2)
end

function ISContiguousLocal(arg1::IS{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:ISContiguousLocal,petsc1),PetscErrorCode,(IS{Float64},Int32,Int32,Ptr{Int32},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function ISGetIndices(arg1::IS{Float64},arg2::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:ISGetIndices,petsc1),PetscErrorCode,(IS{Float64},Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISRestoreIndices(arg1::IS{Float64},arg2::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:ISRestoreIndices,petsc1),PetscErrorCode,(IS{Float64},Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISGetTotalIndices(arg1::IS{Float64},arg2::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:ISGetTotalIndices,petsc1),PetscErrorCode,(IS{Float64},Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISRestoreTotalIndices(arg1::IS{Float64},arg2::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:ISRestoreTotalIndices,petsc1),PetscErrorCode,(IS{Float64},Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISGetNonlocalIndices(arg1::IS{Float64},arg2::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:ISGetNonlocalIndices,petsc1),PetscErrorCode,(IS{Float64},Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISRestoreNonlocalIndices(arg1::IS{Float64},arg2::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:ISRestoreNonlocalIndices,petsc1),PetscErrorCode,(IS{Float64},Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISGetNonlocalIS(arg1::IS{Float64},is::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISGetNonlocalIS,petsc1),PetscErrorCode,(IS{Float64},Ptr{IS{Float64}}),arg1,is)
end

function ISRestoreNonlocalIS(arg1::IS{Float64},is::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISRestoreNonlocalIS,petsc1),PetscErrorCode,(IS{Float64},Ptr{IS{Float64}}),arg1,is)
end

function ISGetSize(arg1::IS{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISGetSize,petsc1),PetscErrorCode,(IS{Float64},Ptr{Int32}),arg1,arg2)
end

function ISGetLocalSize(arg1::IS{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISGetLocalSize,petsc1),PetscErrorCode,(IS{Float64},Ptr{Int32}),arg1,arg2)
end

function ISInvertPermutation(arg1::IS{Float64},arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISInvertPermutation,petsc1),PetscErrorCode,(IS{Float64},Int32,Ptr{IS{Float64}}),arg1,arg2,arg3)
end

function ISView(arg1::IS{Float64},arg2::PetscViewer{Float64})
    ccall((:ISView,petsc1),PetscErrorCode,(IS{Float64},PetscViewer{Float64}),arg1,arg2)
end

function ISEqual(arg1::IS{Float64},arg2::IS{Float64},arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:ISEqual,petsc1),PetscErrorCode,(IS{Float64},IS{Float64},Ptr{PetscBool}),arg1,arg2,arg3)
end

function ISSort(arg1::IS{Float64})
    ccall((:ISSort,petsc1),PetscErrorCode,(IS{Float64},),arg1)
end

function ISSortRemoveDups(arg1::IS{Float64})
    ccall((:ISSortRemoveDups,petsc1),PetscErrorCode,(IS{Float64},),arg1)
end

function ISSorted(arg1::IS{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:ISSorted,petsc1),PetscErrorCode,(IS{Float64},Ptr{PetscBool}),arg1,arg2)
end

function ISDifference(arg1::IS{Float64},arg2::IS{Float64},arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISDifference,petsc1),PetscErrorCode,(IS{Float64},IS{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3)
end

function ISSum(arg1::IS{Float64},arg2::IS{Float64},arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISSum,petsc1),PetscErrorCode,(IS{Float64},IS{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3)
end

function ISExpand(arg1::IS{Float64},arg2::IS{Float64},arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISExpand,petsc1),PetscErrorCode,(IS{Float64},IS{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3)
end

function ISGetMinMax(arg1::IS{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISGetMinMax,petsc1),PetscErrorCode,(IS{Float64},Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3)
end

function ISBlockGetIndices(arg1::IS{Float64},arg2::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:ISBlockGetIndices,petsc1),PetscErrorCode,(IS{Float64},Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISBlockRestoreIndices(arg1::IS{Float64},arg2::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:ISBlockRestoreIndices,petsc1),PetscErrorCode,(IS{Float64},Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISBlockGetLocalSize(arg1::IS{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISBlockGetLocalSize,petsc1),PetscErrorCode,(IS{Float64},Ptr{Int32}),arg1,arg2)
end

function ISBlockGetSize(arg1::IS{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISBlockGetSize,petsc1),PetscErrorCode,(IS{Float64},Ptr{Int32}),arg1,arg2)
end

function ISGetBlockSize(arg1::IS{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISGetBlockSize,petsc1),PetscErrorCode,(IS{Float64},Ptr{Int32}),arg1,arg2)
end

function ISSetBlockSize(arg1::IS{Float64},arg2::Integer)
    ccall((:ISSetBlockSize,petsc1),PetscErrorCode,(IS{Float64},Int32),arg1,arg2)
end

function ISStrideGetInfo(arg1::IS{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISStrideGetInfo,petsc1),PetscErrorCode,(IS{Float64},Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3)
end

function ISToGeneral(arg1::IS{Float64})
    ccall((:ISToGeneral,petsc1),PetscErrorCode,(IS{Float64},),arg1)
end

function ISDuplicate(arg1::IS{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISDuplicate,petsc1),PetscErrorCode,(IS{Float64},Ptr{IS{Float64}}),arg1,arg2)
end

function ISCopy(arg1::IS{Float64},arg2::IS{Float64})
    ccall((:ISCopy,petsc1),PetscErrorCode,(IS{Float64},IS{Float64}),arg1,arg2)
end

function ISAllGather(arg1::IS{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISAllGather,petsc1),PetscErrorCode,(IS{Float64},Ptr{IS{Float64}}),arg1,arg2)
end

function ISComplement(arg1::IS{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISComplement,petsc1),PetscErrorCode,(IS{Float64},Int32,Int32,Ptr{IS{Float64}}),arg1,arg2,arg3,arg4)
end

function ISConcatenate(arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISConcatenate,petsc1),PetscErrorCode,(comm_type,Int32,Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1.val,arg2,arg3,arg4)
end

function ISListToPair(arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg5::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISListToPair,petsc1),PetscErrorCode,(comm_type,Int32,Ptr{IS{Float64}},Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1.val,arg2,arg3,arg4,arg5)
end

function ISPairToList(arg1::IS{Float64},arg2::IS{Float64},arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    ccall((:ISPairToList,petsc1),PetscErrorCode,(IS{Float64},IS{Float64},Ptr{Int32},Ptr{Ptr{IS{Float64}}}),arg1,arg2,arg3,arg4)
end

function ISEmbed(arg1::IS{Float64},arg2::IS{Float64},arg3::PetscBool,arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISEmbed,petsc1),PetscErrorCode,(IS{Float64},IS{Float64},PetscBool,Ptr{IS{Float64}}),arg1,arg2,arg3,arg4)
end

function ISSortPermutation(arg1::IS{Float64},arg2::PetscBool,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISSortPermutation,petsc1),PetscErrorCode,(IS{Float64},PetscBool,Ptr{IS{Float64}}),arg1,arg2,arg3)
end

function ISOnComm(arg1::IS{Float64},arg2::MPI_Comm,PetscCopyMode::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISOnComm,petsc1),PetscErrorCode,(IS{Float64},comm_type,Cint,Ptr{IS{Float64}}),arg1,arg2.val,PetscCopyMode,arg3)
end

function ISLocalToGlobalMappingCreate(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),PetscCopyMode::Integer,arg5::Union(Ptr{ISLocalToGlobalMapping{Float64}},StridedArray{ISLocalToGlobalMapping{Float64}},Ptr{Void}))
    ccall((:ISLocalToGlobalMappingCreate,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Ptr{Int32},Cint,Ptr{ISLocalToGlobalMapping{Float64}}),arg1.val,arg2,arg3,arg4,PetscCopyMode,arg5)
end

function ISLocalToGlobalMappingCreateIS(arg1::IS{Float64},arg2::Union(Ptr{ISLocalToGlobalMapping{Float64}},StridedArray{ISLocalToGlobalMapping{Float64}},Ptr{Void}))
    ccall((:ISLocalToGlobalMappingCreateIS,petsc1),PetscErrorCode,(IS{Float64},Ptr{ISLocalToGlobalMapping{Float64}}),arg1,arg2)
end

#= skipping function with undefined symbols: 
 function ISLocalToGlobalMappingCreateSF(arg1::PetscSF,arg2::Integer,arg3::Union(Ptr{ISLocalToGlobalMapping{Float64}},StridedArray{ISLocalToGlobalMapping{Float64}},Ptr{Void}))
    ccall((:ISLocalToGlobalMappingCreateSF,petsc1),PetscErrorCode,(PetscSF,Int32,Ptr{ISLocalToGlobalMapping{Float64}}),arg1,arg2,arg3)
end 
=#
function ISLocalToGlobalMappingView(arg1::ISLocalToGlobalMapping{Float64},arg2::PetscViewer{Float64})
    ccall((:ISLocalToGlobalMappingView,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},PetscViewer{Float64}),arg1,arg2)
end

function ISLocalToGlobalMappingDestroy(arg1::Union(Ptr{ISLocalToGlobalMapping{Float64}},StridedArray{ISLocalToGlobalMapping{Float64}},Ptr{Void}))
    ccall((:ISLocalToGlobalMappingDestroy,petsc1),PetscErrorCode,(Ptr{ISLocalToGlobalMapping{Float64}},),arg1)
end

function ISLocalToGlobalMappingApply(arg1::ISLocalToGlobalMapping{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISLocalToGlobalMappingApply,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},Int32,Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3,arg4)
end

function ISLocalToGlobalMappingApplyBlock(arg1::ISLocalToGlobalMapping{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISLocalToGlobalMappingApplyBlock,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},Int32,Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3,arg4)
end

function ISLocalToGlobalMappingApplyIS(arg1::ISLocalToGlobalMapping{Float64},arg2::IS{Float64},arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISLocalToGlobalMappingApplyIS,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},IS{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3)
end

function ISGlobalToLocalMappingApply(arg1::ISLocalToGlobalMapping{Float64},arg2::ISGlobalToLocalMappingType,arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISGlobalToLocalMappingApply,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},ISGlobalToLocalMappingType,Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function ISGlobalToLocalMappingApplyBlock(arg1::ISLocalToGlobalMapping{Float64},arg2::ISGlobalToLocalMappingType,arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISGlobalToLocalMappingApplyBlock,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},ISGlobalToLocalMappingType,Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function ISGlobalToLocalMappingApplyIS(arg1::ISLocalToGlobalMapping{Float64},arg2::ISGlobalToLocalMappingType,arg3::IS{Float64},arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISGlobalToLocalMappingApplyIS,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},ISGlobalToLocalMappingType,IS{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4)
end

function ISLocalToGlobalMappingGetSize(arg1::ISLocalToGlobalMapping{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISLocalToGlobalMappingGetSize,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},Ptr{Int32}),arg1,arg2)
end

function ISLocalToGlobalMappingGetInfo(arg1::ISLocalToGlobalMapping{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg4::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg5::Union(Ptr{Ptr{Ptr{Int32}}},StridedArray{Ptr{Ptr{Int32}}},Ptr{Void}))
    ccall((:ISLocalToGlobalMappingGetInfo,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},Ptr{Int32},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}},Ptr{Ptr{Ptr{Int32}}}),arg1,arg2,arg3,arg4,arg5)
end

function ISLocalToGlobalMappingRestoreInfo(arg1::ISLocalToGlobalMapping{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg4::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg5::Union(Ptr{Ptr{Ptr{Int32}}},StridedArray{Ptr{Ptr{Int32}}},Ptr{Void}))
    ccall((:ISLocalToGlobalMappingRestoreInfo,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},Ptr{Int32},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}},Ptr{Ptr{Ptr{Int32}}}),arg1,arg2,arg3,arg4,arg5)
end

function ISLocalToGlobalMappingGetBlockInfo(arg1::ISLocalToGlobalMapping{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg4::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg5::Union(Ptr{Ptr{Ptr{Int32}}},StridedArray{Ptr{Ptr{Int32}}},Ptr{Void}))
    ccall((:ISLocalToGlobalMappingGetBlockInfo,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},Ptr{Int32},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}},Ptr{Ptr{Ptr{Int32}}}),arg1,arg2,arg3,arg4,arg5)
end

function ISLocalToGlobalMappingRestoreBlockInfo(arg1::ISLocalToGlobalMapping{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg4::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg5::Union(Ptr{Ptr{Ptr{Int32}}},StridedArray{Ptr{Ptr{Int32}}},Ptr{Void}))
    ccall((:ISLocalToGlobalMappingRestoreBlockInfo,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},Ptr{Int32},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}},Ptr{Ptr{Ptr{Int32}}}),arg1,arg2,arg3,arg4,arg5)
end

function ISLocalToGlobalMappingGetIndices(arg1::ISLocalToGlobalMapping{Float64},arg2::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:ISLocalToGlobalMappingGetIndices,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISLocalToGlobalMappingRestoreIndices(arg1::ISLocalToGlobalMapping{Float64},arg2::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:ISLocalToGlobalMappingRestoreIndices,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISLocalToGlobalMappingGetBlockIndices(arg1::ISLocalToGlobalMapping{Float64},arg2::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:ISLocalToGlobalMappingGetBlockIndices,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISLocalToGlobalMappingRestoreBlockIndices(arg1::ISLocalToGlobalMapping{Float64},arg2::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:ISLocalToGlobalMappingRestoreBlockIndices,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISLocalToGlobalMappingConcatenate(arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{ISLocalToGlobalMapping{Float64}},StridedArray{ISLocalToGlobalMapping{Float64}},Ptr{Void}),arg4::Union(Ptr{ISLocalToGlobalMapping{Float64}},StridedArray{ISLocalToGlobalMapping{Float64}},Ptr{Void}))
    ccall((:ISLocalToGlobalMappingConcatenate,petsc1),PetscErrorCode,(comm_type,Int32,Ptr{ISLocalToGlobalMapping{Float64}},Ptr{ISLocalToGlobalMapping{Float64}}),arg1.val,arg2,arg3,arg4)
end

function ISG2LMapApply(arg1::ISLocalToGlobalMapping{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISG2LMapApply,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},Int32,Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3,arg4)
end

function ISLocalToGlobalMappingGetBlockSize(arg1::ISLocalToGlobalMapping{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISLocalToGlobalMappingGetBlockSize,petsc1),PetscErrorCode,(ISLocalToGlobalMapping{Float64},Ptr{Int32}),arg1,arg2)
end

function ISAllGatherColors(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}))
    ccall((:ISAllGatherColors,petsc1),PetscErrorCode,(comm_type,Int32,Ptr{Cint},Ptr{Int32},Ptr{Ptr{Cint}}),arg1.val,arg2,arg3,arg4,arg5)
end

function ISColoringCreate(arg1::MPI_Comm,arg2::Integer,arg3::Integer,ISColoringValue::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),PetscCopyMode::Integer,arg4::Union(Ptr{ISColoring{Float64}},StridedArray{ISColoring{Float64}},Ptr{Void}))
    ccall((:ISColoringCreate,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Ptr{Cint},Cint,Ptr{ISColoring{Float64}}),arg1.val,arg2,arg3,ISColoringValue,PetscCopyMode,arg4)
end

function ISColoringDestroy(arg1::Union(Ptr{ISColoring{Float64}},StridedArray{ISColoring{Float64}},Ptr{Void}))
    ccall((:ISColoringDestroy,petsc1),PetscErrorCode,(Ptr{ISColoring{Float64}},),arg1)
end

function ISColoringView(arg1::ISColoring{Float64},arg2::PetscViewer{Float64})
    ccall((:ISColoringView,petsc1),PetscErrorCode,(ISColoring{Float64},PetscViewer{Float64}),arg1,arg2)
end

#= skipping function with undefined symbols: 
 function ISColoringViewFromOptions(arg1::ISColoring{Float64},arg2::PetscObject,arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:ISColoringViewFromOptions,petsc1),PetscErrorCode,(ISColoring{Float64},PetscObject,Cstring),arg1,arg2,arg3)
end 
=#
function ISColoringGetIS(arg1::ISColoring{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    ccall((:ISColoringGetIS,petsc1),PetscErrorCode,(ISColoring{Float64},Ptr{Int32},Ptr{Ptr{IS{Float64}}}),arg1,arg2,arg3)
end

function ISColoringRestoreIS(arg1::ISColoring{Float64},arg2::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    ccall((:ISColoringRestoreIS,petsc1),PetscErrorCode,(ISColoring{Float64},Ptr{Ptr{IS{Float64}}}),arg1,arg2)
end

function ISColoringReference(arg1::ISColoring{Float64})
    ccall((:ISColoringReference,petsc1),PetscErrorCode,(ISColoring{Float64},),arg1)
end

function ISColoringSetType(arg1::ISColoring{Float64},arg2::ISColoringType)
    ccall((:ISColoringSetType,petsc1),PetscErrorCode,(ISColoring{Float64},ISColoringType),arg1,arg2)
end

function ISPartitioningToNumbering(arg1::IS{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISPartitioningToNumbering,petsc1),PetscErrorCode,(IS{Float64},Ptr{IS{Float64}}),arg1,arg2)
end

function ISPartitioningCount(arg1::IS{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISPartitioningCount,petsc1),PetscErrorCode,(IS{Float64},Int32,Ptr{Int32}),arg1,arg2,arg3)
end

function ISCompressIndicesGeneral(arg1::Integer,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg6::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISCompressIndicesGeneral,petsc1),PetscErrorCode,(Int32,Int32,Int32,Int32,Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function ISCompressIndicesSorted(arg1::Integer,arg2::Integer,arg3::Integer,arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg5::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISCompressIndicesSorted,petsc1),PetscErrorCode,(Int32,Int32,Int32,Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4,arg5)
end

function ISExpandIndicesGeneral(arg1::Integer,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg6::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISExpandIndicesGeneral,petsc1),PetscErrorCode,(Int32,Int32,Int32,Int32,Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscLayoutSetUp(arg1::PetscLayout{Float64})
    ccall((:PetscLayoutSetUp,petsc1),PetscErrorCode,(PetscLayout{Float64},),arg1)
end

function PetscLayoutDestroy(arg1::Union(Ptr{PetscLayout{Float64}},StridedArray{PetscLayout{Float64}},Ptr{Void}))
    ccall((:PetscLayoutDestroy,petsc1),PetscErrorCode,(Ptr{PetscLayout{Float64}{Float64}},),arg1)
end

function PetscLayoutDuplicate(arg1::PetscLayout{Float64},arg2::Union(Ptr{PetscLayout{Float64}},StridedArray{PetscLayout{Float64}},Ptr{Void}))
    ccall((:PetscLayoutDuplicate,petsc1),PetscErrorCode,(PetscLayout{Float64},Ptr{PetscLayout{Float64}{Float64}}),arg1,arg2)
end

function PetscLayoutReference(arg1::PetscLayout{Float64},arg2::Union(Ptr{PetscLayout{Float64}},StridedArray{PetscLayout{Float64}},Ptr{Void}))
    ccall((:PetscLayoutReference,petsc1),PetscErrorCode,(PetscLayout{Float64},Ptr{PetscLayout{Float64}{Float64}}),arg1,arg2)
end

function PetscLayoutSetLocalSize(arg1::PetscLayout{Float64},arg2::Integer)
    ccall((:PetscLayoutSetLocalSize,petsc1),PetscErrorCode,(PetscLayout{Float64},Int32),arg1,arg2)
end

function PetscLayoutGetLocalSize(arg1::PetscLayout{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscLayoutGetLocalSize,petsc1),PetscErrorCode,(PetscLayout{Float64},Ptr{Int32}),arg1,arg2)
end

function PetscLayoutSetSize(arg1::PetscLayout{Float64},arg2::Integer)
    ccall((:PetscLayoutSetSize,petsc1),PetscErrorCode,(PetscLayout{Float64},Int32),arg1,arg2)
end

function PetscLayoutGetSize(arg1::PetscLayout{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscLayoutGetSize,petsc1),PetscErrorCode,(PetscLayout{Float64},Ptr{Int32}),arg1,arg2)
end

function PetscLayoutSetBlockSize(arg1::PetscLayout{Float64},arg2::Integer)
    ccall((:PetscLayoutSetBlockSize,petsc1),PetscErrorCode,(PetscLayout{Float64},Int32),arg1,arg2)
end

function PetscLayoutGetBlockSize(arg1::PetscLayout{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscLayoutGetBlockSize,petsc1),PetscErrorCode,(PetscLayout{Float64},Ptr{Int32}),arg1,arg2)
end

function PetscLayoutGetRange(arg1::PetscLayout{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscLayoutGetRange,petsc1),PetscErrorCode,(PetscLayout{Float64},Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3)
end

function PetscLayoutGetRanges(arg1::PetscLayout{Float64},arg2::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:PetscLayoutGetRanges,petsc1),PetscErrorCode,(PetscLayout{Float64},Ptr{Ptr{Int32}}),arg1,arg2)
end

function PetscLayoutSetISLocalToGlobalMapping(arg1::PetscLayout{Float64},arg2::ISLocalToGlobalMapping{Float64})
    ccall((:PetscLayoutSetISLocalToGlobalMapping,petsc1),PetscErrorCode,(PetscLayout{Float64},ISLocalToGlobalMapping{Float64}),arg1,arg2)
end

#= skipping function with undefined symbols: 
 function PetscSFSetGraphLayout(arg1::PetscSF,arg2::PetscLayout{Float64},arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),PetscCopyMode::Integer,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSFSetGraphLayout,petsc1),PetscErrorCode,(PetscSF,PetscLayout{Float64},Int32,Ptr{Int32},Cint,Ptr{Int32}),arg1,arg2,arg3,arg4,PetscCopyMode,arg5)
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
 function PetscSectionGetNumFields(arg0::Type{Float64},arg1::PetscSection,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetNumFields,petsc1),PetscErrorCode,(PetscSection,Ptr{Int32}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetNumFields(arg0::Type{Float64},arg1::PetscSection,arg2::Integer)
    ccall((:PetscSectionSetNumFields,petsc1),PetscErrorCode,(PetscSection,Int32),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetFieldName(arg1::PetscSection,arg2::Integer,arg3::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PetscSectionGetFieldName,petsc1),PetscErrorCode,(PetscSection,Int32,Ptr{Ptr{UInt8}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetFieldName(arg1::PetscSection,arg2::Integer,arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscSectionSetFieldName,petsc1),PetscErrorCode,(PetscSection,Int32,Cstring),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetFieldComponents(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetFieldComponents,petsc1),PetscErrorCode,(PetscSection,Int32,Ptr{Int32}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetFieldComponents(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer)
    ccall((:PetscSectionSetFieldComponents,petsc1),PetscErrorCode,(PetscSection,Int32,Int32),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetChart(arg0::Type{Float64},arg1::PetscSection,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetChart,petsc1),PetscErrorCode,(PetscSection,Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetChart(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer)
    ccall((:PetscSectionSetChart,petsc1),PetscErrorCode,(PetscSection,Int32,Int32),arg1,arg2,arg3)
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
 function PetscSectionGetDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetDof,petsc1),PetscErrorCode,(PetscSection,Int32,Ptr{Int32}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer)
    ccall((:PetscSectionSetDof,petsc1),PetscErrorCode,(PetscSection,Int32,Int32),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionAddDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer)
    ccall((:PetscSectionAddDof,petsc1),PetscErrorCode,(PetscSection,Int32,Int32),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetFieldDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetFieldDof,petsc1),PetscErrorCode,(PetscSection,Int32,Int32,Ptr{Int32}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetFieldDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:PetscSectionSetFieldDof,petsc1),PetscErrorCode,(PetscSection,Int32,Int32,Int32),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionAddFieldDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:PetscSectionAddFieldDof,petsc1),PetscErrorCode,(PetscSection,Int32,Int32,Int32),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionHasConstraints(arg0::Type{Float64},arg1::PetscSection,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscSectionHasConstraints,petsc1),PetscErrorCode,(PetscSection,Ptr{PetscBool}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetConstraintDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetConstraintDof,petsc1),PetscErrorCode,(PetscSection,Int32,Ptr{Int32}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetConstraintDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer)
    ccall((:PetscSectionSetConstraintDof,petsc1),PetscErrorCode,(PetscSection,Int32,Int32),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionAddConstraintDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer)
    ccall((:PetscSectionAddConstraintDof,petsc1),PetscErrorCode,(PetscSection,Int32,Int32),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetFieldConstraintDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetFieldConstraintDof,petsc1),PetscErrorCode,(PetscSection,Int32,Int32,Ptr{Int32}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetFieldConstraintDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:PetscSectionSetFieldConstraintDof,petsc1),PetscErrorCode,(PetscSection,Int32,Int32,Int32),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionAddFieldConstraintDof(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:PetscSectionAddFieldConstraintDof,petsc1),PetscErrorCode,(PetscSection,Int32,Int32,Int32),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetConstraintIndices(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:PetscSectionGetConstraintIndices,petsc1),PetscErrorCode,(PetscSection,Int32,Ptr{Ptr{Int32}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetConstraintIndices(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionSetConstraintIndices,petsc1),PetscErrorCode,(PetscSection,Int32,Ptr{Int32}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetFieldConstraintIndices(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:PetscSectionGetFieldConstraintIndices,petsc1),PetscErrorCode,(PetscSection,Int32,Int32,Ptr{Ptr{Int32}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetFieldConstraintIndices(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionSetFieldConstraintIndices,petsc1),PetscErrorCode,(PetscSection,Int32,Int32,Ptr{Int32}),arg1,arg2,arg3,arg4)
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
 function PetscSectionGetMaxDof(arg0::Type{Float64},arg1::PetscSection,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetMaxDof,petsc1),PetscErrorCode,(PetscSection,Ptr{Int32}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetStorageSize(arg0::Type{Float64},arg1::PetscSection,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetStorageSize,petsc1),PetscErrorCode,(PetscSection,Ptr{Int32}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetConstrainedStorageSize(arg0::Type{Float64},arg1::PetscSection,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetConstrainedStorageSize,petsc1),PetscErrorCode,(PetscSection,Ptr{Int32}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetOffset(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetOffset,petsc1),PetscErrorCode,(PetscSection,Int32,Ptr{Int32}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetOffset(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer)
    ccall((:PetscSectionSetOffset,petsc1),PetscErrorCode,(PetscSection,Int32,Int32),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetFieldOffset(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetFieldOffset,petsc1),PetscErrorCode,(PetscSection,Int32,Int32,Ptr{Int32}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionSetFieldOffset(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:PetscSectionSetFieldOffset,petsc1),PetscErrorCode,(PetscSection,Int32,Int32,Int32),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionGetOffsetRange(arg0::Type{Float64},arg1::PetscSection,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetOffsetRange,petsc1),PetscErrorCode,(PetscSection,Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3)
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
 function PetscSectionCreateGlobalSectionCensored(arg0::Type{Float64},arg1::PetscSection,arg2::PetscSF,arg3::PetscBool,arg4::Integer,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}))
    ccall((:PetscSectionCreateGlobalSectionCensored,petsc1),PetscErrorCode,(PetscSection,PetscSF,PetscBool,Int32,Ptr{Int32},Ptr{PetscSection}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionCreateSubsection(arg0::Type{Float64},arg1::PetscSection,arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{PetscSection},StridedArray{PetscSection},Ptr{Void}))
    ccall((:PetscSectionCreateSubsection,petsc1),PetscErrorCode,(PetscSection,Int32,Ptr{Int32},Ptr{PetscSection}),arg1,arg2,arg3,arg4)
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
    ccall((:PetscSectionGetField,petsc1),PetscErrorCode,(PetscSection,Int32,Ptr{PetscSection}),arg1,arg2,arg3)
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
 function PetscSFCreateRemoteOffsets(arg0::Type{Float64},arg1::PetscSF,arg2::PetscSection,arg3::PetscSection,arg4::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:PetscSFCreateRemoteOffsets,petsc1),PetscErrorCode,(PetscSF,PetscSection,PetscSection,Ptr{Ptr{Int32}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSFDistributeSection(arg0::Type{Float64},arg1::PetscSF,arg2::PetscSection,arg3::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg4::PetscSection)
    ccall((:PetscSFDistributeSection,petsc1),PetscErrorCode,(PetscSF,PetscSection,Ptr{Ptr{Int32}},PetscSection),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSFCreateSectionSF(arg0::Type{Float64},arg1::PetscSF,arg2::PetscSection,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::PetscSection,arg5::Union(Ptr{PetscSF},StridedArray{PetscSF},Ptr{Void}))
    ccall((:PetscSFCreateSectionSF,petsc1),PetscErrorCode,(PetscSF,PetscSection,Ptr{Int32},PetscSection,Ptr{PetscSF}),arg1,arg2,arg3,arg4,arg5)
end 
=#
function PetscViewerInitializePackage(arg0::Type{Float64})
    ccall((:PetscViewerInitializePackage,petsc1),PetscErrorCode,())
end

function PetscViewerRegister(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscViewerRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
end

function PetscViewerCreate(arg1::MPI_Comm,arg2::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    ccall((:PetscViewerCreate,petsc1),PetscErrorCode,(comm_type,Ptr{PetscViewer{Float64}}),arg1.val,arg2)
end

function PetscViewerSetFromOptions(arg1::PetscViewer{Float64})
    ccall((:PetscViewerSetFromOptions,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
end

#= skipping function with undefined symbols: 
 function PetscViewerASCIIOpenWithFILE(arg1::MPI_Comm,arg2::Union(Ptr{FILE},StridedArray{FILE},Ptr{Void}),arg3::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    ccall((:PetscViewerASCIIOpenWithFILE,petsc1),PetscErrorCode,(comm_type,Ptr{FILE},Ptr{PetscViewer{Float64}}),arg1.val,arg2,arg3)
end 
=#
function PetscViewerASCIIOpen(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    ccall((:PetscViewerASCIIOpen,petsc1),PetscErrorCode,(comm_type,Cstring,Ptr{PetscViewer{Float64}}),arg1.val,arg2,arg3)
end

#= skipping function with undefined symbols: 
 function PetscViewerASCIISetFILE(arg1::PetscViewer{Float64},arg2::Union(Ptr{FILE},StridedArray{FILE},Ptr{Void}))
    ccall((:PetscViewerASCIISetFILE,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{FILE}),arg1,arg2)
end 
=#
function PetscViewerBinaryOpen(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::PetscFileMode,arg4::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    ccall((:PetscViewerBinaryOpen,petsc1),PetscErrorCode,(comm_type,Cstring,PetscFileMode,Ptr{PetscViewer{Float64}}),arg1.val,arg2,arg3,arg4)
end

function PetscViewerBinaryGetFlowControl(arg1::PetscViewer{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscViewerBinaryGetFlowControl,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Int32}),arg1,arg2)
end

function PetscViewerBinarySetFlowControl(arg1::PetscViewer{Float64},arg2::Integer)
    ccall((:PetscViewerBinarySetFlowControl,petsc1),PetscErrorCode,(PetscViewer{Float64},Int32),arg1,arg2)
end

function PetscViewerBinarySetUseMPIIO(arg1::PetscViewer{Float64},arg2::PetscBool)
    ccall((:PetscViewerBinarySetUseMPIIO,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscBool),arg1,arg2)
end

function PetscViewerBinaryGetUseMPIIO(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscViewerBinaryGetUseMPIIO,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{PetscBool}),arg1,arg2)
end

function PetscViewerSocketOpen(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Integer,arg4::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    ccall((:PetscViewerSocketOpen,petsc1),PetscErrorCode,(comm_type,Cstring,Cint,Ptr{PetscViewer{Float64}}),arg1.val,arg2,arg3,arg4)
end

function PetscViewerStringOpen(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),size_t::Integer,arg3::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    ccall((:PetscViewerStringOpen,petsc1),PetscErrorCode,(comm_type,Cstring,Cint,Ptr{PetscViewer{Float64}}),arg1.val,arg2,size_t,arg3)
end

function PetscViewerDrawOpen(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    ccall((:PetscViewerDrawOpen,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Cint,Cint,Cint,Cint,Ptr{PetscViewer{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscViewerDrawSetDrawType(arg1::PetscViewer{Float64},arg2::PetscDrawType)
    ccall((:PetscViewerDrawSetDrawType,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring),arg1,arg2)
end

function PetscViewerMathematicaOpen(arg1::MPI_Comm,arg2::Integer,arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    ccall((:PetscViewerMathematicaOpen,petsc1),PetscErrorCode,(comm_type,Cint,Cstring,Cstring,Ptr{PetscViewer{Float64}}),arg1.val,arg2,arg3,arg4,arg5)
end

function PetscViewerSiloOpen(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    ccall((:PetscViewerSiloOpen,petsc1),PetscErrorCode,(comm_type,Cstring,Ptr{PetscViewer{Float64}}),arg1.val,arg2,arg3)
end

function PetscViewerMatlabOpen(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::PetscFileMode,arg4::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    ccall((:PetscViewerMatlabOpen,petsc1),PetscErrorCode,(comm_type,Cstring,PetscFileMode,Ptr{PetscViewer{Float64}}),arg1.val,arg2,arg3,arg4)
end

function PetscViewerGetType(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscViewerType},StridedArray{PetscViewerType},Ptr{Void}))
    ccall((:PetscViewerGetType,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Ptr{UInt8}}),arg1,arg2)
end

function PetscViewerSetType(arg1::PetscViewer{Float64},arg2::PetscViewerType)
    ccall((:PetscViewerSetType,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring),arg1,arg2)
end

function PetscViewerDestroy(arg1::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    ccall((:PetscViewerDestroy,petsc1),PetscErrorCode,(Ptr{PetscViewer{Float64}},),arg1)
end

function PetscViewerGetSingleton(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    ccall((:PetscViewerGetSingleton,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{PetscViewer{Float64}}),arg1,arg2)
end

function PetscViewerRestoreSingleton(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    ccall((:PetscViewerRestoreSingleton,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{PetscViewer{Float64}}),arg1,arg2)
end

function PetscViewerGetSubcomm(arg1::PetscViewer{Float64},arg2::MPI_Comm,arg3::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    ccall((:PetscViewerGetSubcomm,petsc1),PetscErrorCode,(PetscViewer{Float64},comm_type,Ptr{PetscViewer{Float64}}),arg1,arg2.val,arg3)
end

function PetscViewerRestoreSubcomm(arg1::PetscViewer{Float64},arg2::MPI_Comm,arg3::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    ccall((:PetscViewerRestoreSubcomm,petsc1),PetscErrorCode,(PetscViewer{Float64},comm_type,Ptr{PetscViewer{Float64}}),arg1,arg2.val,arg3)
end

function PetscViewerSetUp(arg1::PetscViewer{Float64})
    ccall((:PetscViewerSetUp,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
end

function PetscViewerView(arg1::PetscViewer{Float64},arg2::PetscViewer{Float64})
    ccall((:PetscViewerView,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscViewer{Float64}),arg1,arg2)
end

function PetscViewerAppendOptionsPrefix(arg1::PetscViewer{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscViewerAppendOptionsPrefix,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring),arg1,arg2)
end

function PetscViewerGetOptionsPrefix(arg1::PetscViewer{Float64},arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PetscViewerGetOptionsPrefix,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Ptr{UInt8}}),arg1,arg2)
end

function PetscViewerSetFormat(arg1::PetscViewer{Float64},arg2::PetscViewerFormat)
    ccall((:PetscViewerSetFormat,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscViewerFormat),arg1,arg2)
end

function PetscViewerPushFormat(arg1::PetscViewer{Float64},arg2::PetscViewerFormat)
    ccall((:PetscViewerPushFormat,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscViewerFormat),arg1,arg2)
end

function PetscViewerPopFormat(arg1::PetscViewer{Float64})
    ccall((:PetscViewerPopFormat,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
end

function PetscViewerGetFormat(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscViewerFormat},StridedArray{PetscViewerFormat},Ptr{Void}))
    ccall((:PetscViewerGetFormat,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{PetscViewerFormat}),arg1,arg2)
end

function PetscViewerFlush(arg1::PetscViewer{Float64})
    ccall((:PetscViewerFlush,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
end

function PetscOptionsGetViewer(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}),arg5::Union(Ptr{PetscViewerFormat},StridedArray{PetscViewerFormat},Ptr{Void}),arg6::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsGetViewer,petsc1),PetscErrorCode,(comm_type,Cstring,Cstring,Ptr{PetscViewer{Float64}},Ptr{PetscViewerFormat},Ptr{PetscBool}),arg1.val,arg2,arg3,arg4,arg5,arg6)
end

#= skipping function with undefined symbols: 
 function PetscOptionsViewer_Private(arg1::Union(Ptr{PetscOptions},StridedArray{PetscOptions},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg5::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}),arg6::Union(Ptr{PetscViewerFormat},StridedArray{PetscViewerFormat},Ptr{Void}),arg7::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsViewer_Private,petsc1),PetscErrorCode,(Ptr{PetscOptions},Cstring,Cstring,Cstring,Ptr{PetscViewer{Float64}},Ptr{PetscViewerFormat},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function PetscViewerASCIIGetPointer(arg1::PetscViewer{Float64},arg2::Union(Ptr{Ptr{FILE}},StridedArray{Ptr{FILE}},Ptr{Void}))
    ccall((:PetscViewerASCIIGetPointer,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Ptr{FILE}}),arg1,arg2)
end 
=#
function PetscViewerFileGetMode(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscFileMode},StridedArray{PetscFileMode},Ptr{Void}))
    ccall((:PetscViewerFileGetMode,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{PetscFileMode}),arg1,arg2)
end

function PetscViewerFileSetMode(arg1::PetscViewer{Float64},arg2::PetscFileMode)
    ccall((:PetscViewerFileSetMode,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscFileMode),arg1,arg2)
end

function PetscViewerRead(arg1::PetscViewer{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::PetscDataType)
    ccall((:PetscViewerRead,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Void},Int32,Ptr{Int32},PetscDataType),arg1,arg2,arg3,arg4,arg5)
end

function PetscViewerASCIISynchronizedAllow(arg1::PetscViewer{Float64},arg2::PetscBool)
    ccall((:PetscViewerASCIISynchronizedAllow,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscBool),arg1,arg2)
end

function PetscViewerASCIIPushTab(arg1::PetscViewer{Float64})
    ccall((:PetscViewerASCIIPushTab,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
end

function PetscViewerASCIIPopTab(arg1::PetscViewer{Float64})
    ccall((:PetscViewerASCIIPopTab,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
end

function PetscViewerASCIIUseTabs(arg1::PetscViewer{Float64},arg2::PetscBool)
    ccall((:PetscViewerASCIIUseTabs,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscBool),arg1,arg2)
end

function PetscViewerASCIISetTab(arg1::PetscViewer{Float64},arg2::Integer)
    ccall((:PetscViewerASCIISetTab,petsc1),PetscErrorCode,(PetscViewer{Float64},Int32),arg1,arg2)
end

function PetscViewerASCIIGetTab(arg1::PetscViewer{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscViewerASCIIGetTab,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Int32}),arg1,arg2)
end

function PetscViewerASCIIAddTab(arg1::PetscViewer{Float64},arg2::Integer)
    ccall((:PetscViewerASCIIAddTab,petsc1),PetscErrorCode,(PetscViewer{Float64},Int32),arg1,arg2)
end

function PetscViewerASCIISubtractTab(arg1::PetscViewer{Float64},arg2::Integer)
    ccall((:PetscViewerASCIISubtractTab,petsc1),PetscErrorCode,(PetscViewer{Float64},Int32),arg1,arg2)
end

function PetscViewerASCIIRead(arg1::PetscViewer{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::PetscDataType)
    ccall((:PetscViewerASCIIRead,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Void},Int32,Ptr{Int32},PetscDataType),arg1,arg2,arg3,arg4,arg5)
end

function PetscViewerBinaryGetDescriptor(arg1::PetscViewer{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscViewerBinaryGetDescriptor,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Cint}),arg1,arg2)
end

#= skipping function with undefined symbols: 
 function PetscViewerBinaryGetInfoPointer(arg1::PetscViewer{Float64},arg2::Union(Ptr{Ptr{FILE}},StridedArray{Ptr{FILE}},Ptr{Void}))
    ccall((:PetscViewerBinaryGetInfoPointer,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Ptr{FILE}}),arg1,arg2)
end 
=#
function PetscViewerBinaryRead(arg1::PetscViewer{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::PetscDataType)
    ccall((:PetscViewerBinaryRead,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Void},Int32,Ptr{Int32},PetscDataType),arg1,arg2,arg3,arg4,arg5)
end

function PetscViewerBinaryWrite(arg1::PetscViewer{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Integer,arg4::PetscDataType,arg5::PetscBool)
    ccall((:PetscViewerBinaryWrite,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Void},Int32,PetscDataType,PetscBool),arg1,arg2,arg3,arg4,arg5)
end

function PetscViewerStringSetString(arg1::PetscViewer{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Integer)
    ccall((:PetscViewerStringSetString,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring,Int32),arg1,arg2,arg3)
end

function PetscViewerDrawClear(arg1::PetscViewer{Float64})
    ccall((:PetscViewerDrawClear,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
end

function PetscViewerDrawSetHold(arg1::PetscViewer{Float64},arg2::PetscBool)
    ccall((:PetscViewerDrawSetHold,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscBool),arg1,arg2)
end

function PetscViewerDrawGetHold(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscViewerDrawGetHold,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{PetscBool}),arg1,arg2)
end

function PetscViewerDrawSetPause(arg1::PetscViewer{Float64},Float64::Integer)
    ccall((:PetscViewerDrawSetPause,petsc1),PetscErrorCode,(PetscViewer{Float64},Cint),arg1,PetscReal)
end

function PetscViewerDrawGetPause(arg1::PetscViewer{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscViewerDrawGetPause,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Cint}),arg1,arg2)
end

function PetscViewerDrawSetInfo(arg1::PetscViewer{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer)
    ccall((:PetscViewerDrawSetInfo,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring,Cstring,Cint,Cint,Cint,Cint),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscViewerDrawResize(arg1::PetscViewer{Float64},arg2::Integer,arg3::Integer)
    ccall((:PetscViewerDrawResize,petsc1),PetscErrorCode,(PetscViewer{Float64},Cint,Cint),arg1,arg2,arg3)
end

function PetscViewerDrawSetBounds(arg1::PetscViewer{Float64},arg2::Integer,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscViewerDrawSetBounds,petsc1),PetscErrorCode,(PetscViewer{Float64},Int32,Ptr{Cint}),arg1,arg2,arg3)
end

function PetscViewerDrawGetBounds(arg1::PetscViewer{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}))
    ccall((:PetscViewerDrawGetBounds,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Int32},Ptr{Ptr{Cint}}),arg1,arg2,arg3)
end

function PetscViewerSocketSetConnection(arg1::PetscViewer{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Integer)
    ccall((:PetscViewerSocketSetConnection,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring,Cint),arg1,arg2,arg3)
end

function PetscViewerBinarySkipInfo(arg1::PetscViewer{Float64})
    ccall((:PetscViewerBinarySkipInfo,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
end

function PetscViewerBinarySetSkipInfo(arg1::PetscViewer{Float64},arg2::PetscBool)
    ccall((:PetscViewerBinarySetSkipInfo,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscBool),arg1,arg2)
end

function PetscViewerBinaryGetSkipInfo(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscViewerBinaryGetSkipInfo,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{PetscBool}),arg1,arg2)
end

function PetscViewerBinarySetSkipOptions(arg1::PetscViewer{Float64},arg2::PetscBool)
    ccall((:PetscViewerBinarySetSkipOptions,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscBool),arg1,arg2)
end

function PetscViewerBinaryGetSkipOptions(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscViewerBinaryGetSkipOptions,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{PetscBool}),arg1,arg2)
end

function PetscViewerBinarySetSkipHeader(arg1::PetscViewer{Float64},arg2::PetscBool)
    ccall((:PetscViewerBinarySetSkipHeader,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscBool),arg1,arg2)
end

function PetscViewerBinaryGetSkipHeader(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscViewerBinaryGetSkipHeader,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{PetscBool}),arg1,arg2)
end

function PetscViewerBinaryReadStringArray(arg1::PetscViewer{Float64},arg2::Union(Ptr{Ptr{ASCIIString}},StridedArray{Ptr{ASCIIString}},Ptr{Void}))
    ccall((:PetscViewerBinaryReadStringArray,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Ptr{Ptr{UInt8}}}),arg1,arg2)
end

function PetscViewerBinaryWriteStringArray(arg1::PetscViewer{Float64},arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PetscViewerBinaryWriteStringArray,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Ptr{UInt8}}),arg1,arg2)
end

function PetscViewerFileSetName(arg1::PetscViewer{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscViewerFileSetName,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring),arg1,arg2)
end

function PetscViewerFileGetName(arg1::PetscViewer{Float64},arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PetscViewerFileGetName,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Ptr{UInt8}}),arg1,arg2)
end

#= skipping function with undefined symbols: 
 function PetscViewerVUGetPointer(arg1::PetscViewer{Float64},arg2::Union(Ptr{Ptr{FILE}},StridedArray{Ptr{FILE}},Ptr{Void}))
    ccall((:PetscViewerVUGetPointer,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Ptr{FILE}}),arg1,arg2)
end 
=#
function PetscViewerVUSetVecSeen(arg1::PetscViewer{Float64},arg2::PetscBool)
    ccall((:PetscViewerVUSetVecSeen,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscBool),arg1,arg2)
end

function PetscViewerVUGetVecSeen(arg1::PetscViewer{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscViewerVUGetVecSeen,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{PetscBool}),arg1,arg2)
end

function PetscViewerVUFlushDeferred(arg1::PetscViewer{Float64})
    ccall((:PetscViewerVUFlushDeferred,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
end

function PetscViewerMathematicaInitializePackage(arg0::Type{Float64})
    ccall((:PetscViewerMathematicaInitializePackage,petsc1),PetscErrorCode,())
end

function PetscViewerMathematicaFinalizePackage(arg0::Type{Float64})
    ccall((:PetscViewerMathematicaFinalizePackage,petsc1),PetscErrorCode,())
end

function PetscViewerMathematicaGetName(arg1::PetscViewer{Float64},arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PetscViewerMathematicaGetName,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Ptr{UInt8}}),arg1,arg2)
end

function PetscViewerMathematicaSetName(arg1::PetscViewer{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscViewerMathematicaSetName,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring),arg1,arg2)
end

function PetscViewerMathematicaClearName(arg1::PetscViewer{Float64})
    ccall((:PetscViewerMathematicaClearName,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
end

function PetscViewerMathematicaSkipPackets(arg1::PetscViewer{Float64},arg2::Integer)
    ccall((:PetscViewerMathematicaSkipPackets,petsc1),PetscErrorCode,(PetscViewer{Float64},Cint),arg1,arg2)
end

function PetscViewerSiloGetName(arg1::PetscViewer{Float64},arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PetscViewerSiloGetName,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Ptr{UInt8}}),arg1,arg2)
end

function PetscViewerSiloSetName(arg1::PetscViewer{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscViewerSiloSetName,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring),arg1,arg2)
end

function PetscViewerSiloClearName(arg1::PetscViewer{Float64})
    ccall((:PetscViewerSiloClearName,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
end

function PetscViewerSiloGetMeshName(arg1::PetscViewer{Float64},arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PetscViewerSiloGetMeshName,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Ptr{UInt8}}),arg1,arg2)
end

function PetscViewerSiloSetMeshName(arg1::PetscViewer{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscViewerSiloSetMeshName,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring),arg1,arg2)
end

function PetscViewerSiloClearMeshName(arg1::PetscViewer{Float64})
    ccall((:PetscViewerSiloClearMeshName,petsc1),PetscErrorCode,(PetscViewer{Float64},),arg1)
end

function PetscViewerNetcdfOpen(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::PetscFileMode,arg4::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    ccall((:PetscViewerNetcdfOpen,petsc1),PetscErrorCode,(comm_type,Cstring,PetscFileMode,Ptr{PetscViewer{Float64}}),arg1.val,arg2,arg3,arg4)
end

function PetscViewerNetcdfGetID(arg1::PetscViewer{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscViewerNetcdfGetID,petsc1),PetscErrorCode,(PetscViewer{Float64},Ptr{Cint}),arg1,arg2)
end

#= skipping function with undefined symbols: 
 function PetscViewerVTKAddField(arg1::PetscViewer{Float64},arg2::PetscObject,PetscViewerVTKWriteFunction::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::PetscViewerVTKFieldType,arg4::PetscObject)
    ccall((:PetscViewerVTKAddField,petsc1),PetscErrorCode,(PetscViewer{Float64},PetscObject,Ptr{Void},PetscViewerVTKFieldType,PetscObject),arg1,arg2,PetscViewerVTKWriteFunction,arg3,arg4)
end 
=#
function PetscViewerVTKOpen(arg1::MPI_Comm,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::PetscFileMode,arg4::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    ccall((:PetscViewerVTKOpen,petsc1),PetscErrorCode,(comm_type,Cstring,PetscFileMode,Ptr{PetscViewer{Float64}}),arg1.val,arg2,arg3,arg4)
end

function PETSC_VIEWER_STDOUT_(arg0::Type{Float64},arg1::MPI_Comm)
    ccall((:PETSC_VIEWER_STDOUT_,petsc1),PetscViewer,(comm_type,),arg1.val)
end

function PetscViewerASCIIGetStdout(arg1::MPI_Comm,arg2::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    ccall((:PetscViewerASCIIGetStdout,petsc1),PetscErrorCode,(comm_type,Ptr{PetscViewer{Float64}}),arg1.val,arg2)
end

function PETSC_VIEWER_STDERR_(arg0::Type{Float64},arg1::MPI_Comm)
    ccall((:PETSC_VIEWER_STDERR_,petsc1),PetscViewer,(comm_type,),arg1.val)
end

function PetscViewerASCIIGetStderr(arg1::MPI_Comm,arg2::Union(Ptr{PetscViewer{Float64}},StridedArray{PetscViewer{Float64}},Ptr{Void}))
    ccall((:PetscViewerASCIIGetStderr,petsc1),PetscErrorCode,(comm_type,Ptr{PetscViewer{Float64}}),arg1.val,arg2)
end

function PETSC_VIEWER_DRAW_(arg0::Type{Float64},arg1::MPI_Comm)
    ccall((:PETSC_VIEWER_DRAW_,petsc1),PetscViewer,(comm_type,),arg1.val)
end

function PETSC_VIEWER_SOCKET_(arg0::Type{Float64},arg1::MPI_Comm)
    ccall((:PETSC_VIEWER_SOCKET_,petsc1),PetscViewer,(comm_type,),arg1.val)
end

function PETSC_VIEWER_BINARY_(arg0::Type{Float64},arg1::MPI_Comm)
    ccall((:PETSC_VIEWER_BINARY_,petsc1),PetscViewer,(comm_type,),arg1.val)
end

function PETSC_VIEWER_MATLAB_(arg0::Type{Float64},arg1::MPI_Comm)
    ccall((:PETSC_VIEWER_MATLAB_,petsc1),PetscViewer,(comm_type,),arg1.val)
end

function PETSC_VIEWER_HDF5_(arg0::Type{Float64},arg1::MPI_Comm)
    ccall((:PETSC_VIEWER_HDF5_,petsc1),PetscViewer,(comm_type,),arg1.val)
end

function PetscViewerMatlabGetArray(arg1::PetscViewer{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PetscViewerMatlabGetArray,petsc1),PetscErrorCode,(PetscViewer{Float64},Cint,Cint,Ptr{Float64},Cstring),arg1,arg2,arg3,arg4,arg5)
end

function PetscViewerMatlabPutVariable(arg1::PetscViewer{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PetscViewerMatlabPutVariable,petsc1),PetscErrorCode,(PetscViewer{Float64},Cstring,Ptr{Void}),arg1,arg2,arg3)
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
    ccall((:PetscViewersGetViewer,petsc1),PetscErrorCode,(PetscViewers,Int32,Ptr{PetscViewer{Float64}}),arg1,arg2,arg3)
end 
=#
function VecInitializePackage(arg0::Type{Float64})
    ccall((:VecInitializePackage,petsc1),PetscErrorCode,())
end

function VecFinalizePackage(arg0::Type{Float64})
    ccall((:VecFinalizePackage,petsc1),PetscErrorCode,())
end

function VecCreate(arg1::MPI_Comm,arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:VecCreate,petsc1),PetscErrorCode,(comm_type,Ptr{Vec{Float64}}),arg1.val,arg2)
end

function VecCreateSeq(arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:VecCreateSeq,petsc1),PetscErrorCode,(comm_type,Int32,Ptr{Vec{Float64}}),arg1.val,arg2,arg3)
end

function VecCreateMPI(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:VecCreateMPI,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Ptr{Vec{Float64}}),arg1.val,arg2,arg3,arg4)
end

function VecCreateSeqWithArray(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:VecCreateSeqWithArray,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Ptr{Float64},Ptr{Vec{Float64}}),arg1.val,arg2,arg3,arg4,arg5)
end

function VecCreateMPIWithArray(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg6::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:VecCreateMPIWithArray,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Ptr{Float64},Ptr{Vec{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6)
end

function VecCreateShared(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:VecCreateShared,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Ptr{Vec{Float64}}),arg1.val,arg2,arg3,arg4)
end

function VecSetFromOptions(arg1::Vec{Float64})
    ccall((:VecSetFromOptions,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
end

function VecDestroy(arg1::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:VecDestroy,petsc1),PetscErrorCode,(Ptr{Vec{Float64}},),arg1)
end

function VecZeroEntries(arg1::Vec{Float64})
    ccall((:VecZeroEntries,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
end

function VecSetOptionsPrefix(arg1::Vec{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:VecSetOptionsPrefix,petsc1),PetscErrorCode,(Vec{Float64},Cstring),arg1,arg2)
end

function VecAppendOptionsPrefix(arg1::Vec{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:VecAppendOptionsPrefix,petsc1),PetscErrorCode,(Vec{Float64},Cstring),arg1,arg2)
end

function VecGetOptionsPrefix(arg1::Vec{Float64},arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:VecGetOptionsPrefix,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Ptr{UInt8}}),arg1,arg2)
end

function VecSetSizes(arg1::Vec{Float64},arg2::Integer,arg3::Integer)
    ccall((:VecSetSizes,petsc1),PetscErrorCode,(Vec{Float64},Int32,Int32),arg1,arg2,arg3)
end

function VecDotNorm2(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:VecDotNorm2,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{Float64},Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function VecDot(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:VecDot,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{Float64}),arg1,arg2,arg3)
end

function VecDotRealPart(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:VecDotRealPart,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{Cint}),arg1,arg2,arg3)
end

function VecTDot(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:VecTDot,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{Float64}),arg1,arg2,arg3)
end

function VecMDot(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:VecMDot,petsc1),PetscErrorCode,(Vec{Float64},Int32,Ptr{Vec{Float64}},Ptr{Float64}),arg1,arg2,arg3,arg4)
end

function VecMTDot(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:VecMTDot,petsc1),PetscErrorCode,(Vec{Float64},Int32,Ptr{Vec{Float64}},Ptr{Float64}),arg1,arg2,arg3,arg4)
end

function VecGetSubVector(arg1::Vec{Float64},arg2::IS{Float64},arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:VecGetSubVector,petsc1),PetscErrorCode,(Vec{Float64},IS{Float64},Ptr{Vec{Float64}}),arg1,arg2,arg3)
end

function VecRestoreSubVector(arg1::Vec{Float64},arg2::IS{Float64},arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:VecRestoreSubVector,petsc1),PetscErrorCode,(Vec{Float64},IS{Float64},Ptr{Vec{Float64}}),arg1,arg2,arg3)
end

function VecNorm(arg1::Vec{Float64},arg2::NormType,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:VecNorm,petsc1),PetscErrorCode,(Vec{Float64},NormType,Ptr{Cint}),arg1,arg2,arg3)
end

function VecNormAvailable(arg1::Vec{Float64},arg2::NormType,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:VecNormAvailable,petsc1),PetscErrorCode,(Vec{Float64},NormType,Ptr{PetscBool},Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function VecNormalize(arg1::Vec{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:VecNormalize,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Cint}),arg1,arg2)
end

function VecSum(arg1::Vec{Float64},arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:VecSum,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Float64}),arg1,arg2)
end

function VecMax(arg1::Vec{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:VecMax,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Int32},Ptr{Cint}),arg1,arg2,arg3)
end

function VecMin(arg1::Vec{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:VecMin,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Int32},Ptr{Cint}),arg1,arg2,arg3)
end

function VecScale(arg1::Vec{Float64},arg2::Float64)
    ccall((:VecScale,petsc1),PetscErrorCode,(Vec{Float64},Float64),arg1,arg2)
end

function VecCopy(arg1::Vec{Float64},arg2::Vec{Float64})
    ccall((:VecCopy,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64}),arg1,arg2)
end

#= skipping function with undefined symbols: 
 function VecSetRandom(arg1::Vec{Float64},arg2::PetscRandom)
    ccall((:VecSetRandom,petsc1),PetscErrorCode,(Vec{Float64},PetscRandom),arg1,arg2)
end 
=#
function VecSet(arg1::Vec{Float64},arg2::Float64)
    ccall((:VecSet,petsc1),PetscErrorCode,(Vec{Float64},Float64),arg1,arg2)
end

function VecSetInf(arg1::Vec{Float64})
    ccall((:VecSetInf,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
end

function VecSwap(arg1::Vec{Float64},arg2::Vec{Float64})
    ccall((:VecSwap,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64}),arg1,arg2)
end

function VecAXPY(arg1::Vec{Float64},arg2::Float64,arg3::Vec{Float64})
    ccall((:VecAXPY,petsc1),PetscErrorCode,(Vec{Float64},Float64,Vec{Float64}),arg1,arg2,arg3)
end

function VecAXPBY(arg1::Vec{Float64},arg2::Float64,arg3::Float64,arg4::Vec{Float64})
    ccall((:VecAXPBY,petsc1),PetscErrorCode,(Vec{Float64},Float64,Float64,Vec{Float64}),arg1,arg2,arg3,arg4)
end

function VecMAXPY(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg4::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:VecMAXPY,petsc1),PetscErrorCode,(Vec{Float64},Int32,Ptr{Float64},Ptr{Vec{Float64}}),arg1,arg2,arg3,arg4)
end

function VecAYPX(arg1::Vec{Float64},arg2::Float64,arg3::Vec{Float64})
    ccall((:VecAYPX,petsc1),PetscErrorCode,(Vec{Float64},Float64,Vec{Float64}),arg1,arg2,arg3)
end

function VecWAXPY(arg1::Vec{Float64},arg2::Float64,arg3::Vec{Float64},arg4::Vec{Float64})
    ccall((:VecWAXPY,petsc1),PetscErrorCode,(Vec{Float64},Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
end

function VecAXPBYPCZ(arg1::Vec{Float64},arg2::Float64,arg3::Float64,arg4::Float64,arg5::Vec{Float64},arg6::Vec{Float64})
    ccall((:VecAXPBYPCZ,petsc1),PetscErrorCode,(Vec{Float64},Float64,Float64,Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecPointwiseMax(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:VecPointwiseMax,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function VecPointwiseMaxAbs(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:VecPointwiseMaxAbs,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function VecPointwiseMin(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:VecPointwiseMin,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function VecPointwiseMult(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:VecPointwiseMult,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function VecPointwiseDivide(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:VecPointwiseDivide,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function VecMaxPointwiseDivide(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:VecMaxPointwiseDivide,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{Cint}),arg1,arg2,arg3)
end

function VecShift(arg1::Vec{Float64},arg2::Float64)
    ccall((:VecShift,petsc1),PetscErrorCode,(Vec{Float64},Float64),arg1,arg2)
end

function VecReciprocal(arg1::Vec{Float64})
    ccall((:VecReciprocal,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
end

function VecPermute(arg1::Vec{Float64},arg2::IS{Float64},arg3::PetscBool)
    ccall((:VecPermute,petsc1),PetscErrorCode,(Vec{Float64},IS{Float64},PetscBool),arg1,arg2,arg3)
end

function VecSqrtAbs(arg1::Vec{Float64})
    ccall((:VecSqrtAbs,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
end

function VecLog(arg1::Vec{Float64})
    ccall((:VecLog,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
end

function VecExp(arg1::Vec{Float64})
    ccall((:VecExp,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
end

function VecAbs(arg1::Vec{Float64})
    ccall((:VecAbs,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
end

function VecDuplicate(arg1::Vec{Float64},arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:VecDuplicate,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Vec{Float64}}),arg1,arg2)
end

function VecDuplicateVecs(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Ptr{Vec{Float64}}},StridedArray{Ptr{Vec{Float64}}},Ptr{Void}))
    ccall((:VecDuplicateVecs,petsc1),PetscErrorCode,(Vec{Float64},Int32,Ptr{Ptr{Vec{Float64}}}),arg1,arg2,arg3)
end

function VecDestroyVecs(arg1::Integer,arg2::Union(Ptr{Ptr{Vec{Float64}}},StridedArray{Ptr{Vec{Float64}}},Ptr{Void}))
    ccall((:VecDestroyVecs,petsc1),PetscErrorCode,(Int32,Ptr{Ptr{Vec{Float64}}}),arg1,arg2)
end

function VecStrideNormAll(arg1::Vec{Float64},arg2::NormType,Float64::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:VecStrideNormAll,petsc1),PetscErrorCode,(Vec{Float64},NormType,Ptr{Cint}),arg1,arg2,PetscReal)
end

function VecStrideMaxAll(arg1::Vec{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),Float64::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:VecStrideMaxAll,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Int32},Ptr{Cint}),arg1,arg2,PetscReal)
end

function VecStrideMinAll(arg1::Vec{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),Float64::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:VecStrideMinAll,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Int32},Ptr{Cint}),arg1,arg2,PetscReal)
end

function VecStrideScaleAll(arg1::Vec{Float64},arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:VecStrideScaleAll,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Float64}),arg1,arg2)
end

function VecUniqueEntries(arg1::Vec{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:VecUniqueEntries,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Int32},Ptr{Ptr{Float64}}),arg1,arg2,arg3)
end

function VecStrideNorm(arg1::Vec{Float64},arg2::Integer,arg3::NormType,arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:VecStrideNorm,petsc1),PetscErrorCode,(Vec{Float64},Int32,NormType,Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function VecStrideMax(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:VecStrideMax,petsc1),PetscErrorCode,(Vec{Float64},Int32,Ptr{Int32},Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function VecStrideMin(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:VecStrideMin,petsc1),PetscErrorCode,(Vec{Float64},Int32,Ptr{Int32},Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function VecStrideScale(arg1::Vec{Float64},arg2::Integer,arg3::Float64)
    ccall((:VecStrideScale,petsc1),PetscErrorCode,(Vec{Float64},Int32,Float64),arg1,arg2,arg3)
end

function VecStrideSet(arg1::Vec{Float64},arg2::Integer,arg3::Float64)
    ccall((:VecStrideSet,petsc1),PetscErrorCode,(Vec{Float64},Int32,Float64),arg1,arg2,arg3)
end

function VecStrideGather(arg1::Vec{Float64},arg2::Integer,arg3::Vec{Float64},arg4::InsertMode)
    ccall((:VecStrideGather,petsc1),PetscErrorCode,(Vec{Float64},Int32,Vec{Float64},InsertMode),arg1,arg2,arg3,arg4)
end

function VecStrideScatter(arg1::Vec{Float64},arg2::Integer,arg3::Vec{Float64},arg4::InsertMode)
    ccall((:VecStrideScatter,petsc1),PetscErrorCode,(Vec{Float64},Int32,Vec{Float64},InsertMode),arg1,arg2,arg3,arg4)
end

function VecStrideGatherAll(arg1::Vec{Float64},arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg3::InsertMode)
    ccall((:VecStrideGatherAll,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Vec{Float64}},InsertMode),arg1,arg2,arg3)
end

function VecStrideScatterAll(arg1::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg2::Vec{Float64},arg3::InsertMode)
    ccall((:VecStrideScatterAll,petsc1),PetscErrorCode,(Ptr{Vec{Float64}},Vec{Float64},InsertMode),arg1,arg2,arg3)
end

function VecStrideSubSetScatter(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Vec{Float64},arg6::InsertMode)
    ccall((:VecStrideSubSetScatter,petsc1),PetscErrorCode,(Vec{Float64},Int32,Ptr{Int32},Ptr{Int32},Vec{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecStrideSubSetGather(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Vec{Float64},arg6::InsertMode)
    ccall((:VecStrideSubSetGather,petsc1),PetscErrorCode,(Vec{Float64},Int32,Ptr{Int32},Ptr{Int32},Vec{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecSetValues(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::InsertMode)
    ccall((:VecSetValues,petsc1),PetscErrorCode,(Vec{Float64},Int32,Ptr{Int32},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5)
end

function VecGetValues(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:VecGetValues,petsc1),PetscErrorCode,(Vec{Float64},Int32,Ptr{Int32},Ptr{Float64}),arg1,arg2,arg3,arg4)
end

function VecAssemblyBegin(arg1::Vec{Float64})
    ccall((:VecAssemblyBegin,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
end

function VecAssemblyEnd(arg1::Vec{Float64})
    ccall((:VecAssemblyEnd,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
end

function VecStashSetInitialSize(arg1::Vec{Float64},arg2::Integer,arg3::Integer)
    ccall((:VecStashSetInitialSize,petsc1),PetscErrorCode,(Vec{Float64},Int32,Int32),arg1,arg2,arg3)
end

function VecStashView(arg1::Vec{Float64},arg2::PetscViewer{Float64})
    ccall((:VecStashView,petsc1),PetscErrorCode,(Vec{Float64},PetscViewer{Float64}),arg1,arg2)
end

#= skipping function with undefined symbols: 
 function VecStashViewFromOptions(arg1::Vec{Float64},arg2::PetscObject,arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:VecStashViewFromOptions,petsc1),PetscErrorCode,(Vec{Float64},PetscObject,Cstring),arg1,arg2,arg3)
end 
=#
function VecStashGetInfo(arg1::Vec{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:VecStashGetInfo,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3,arg4,arg5)
end

function VecGetBlockSize(arg1::Vec{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:VecGetBlockSize,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Int32}),arg1,arg2)
end

function VecSetValuesBlocked(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::InsertMode)
    ccall((:VecSetValuesBlocked,petsc1),PetscErrorCode,(Vec{Float64},Int32,Ptr{Int32},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5)
end

function VecSetType(arg1::Vec{Float64},arg2::VecType)
    ccall((:VecSetType,petsc1),PetscErrorCode,(Vec{Float64},Cstring),arg1,arg2)
end

function VecGetType(arg1::Vec{Float64},arg2::Union(Ptr{VecType},StridedArray{VecType},Ptr{Void}))
    ccall((:VecGetType,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Cstring}),arg1,arg2)
end

function VecRegister(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:VecRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
end

function VecScatterCreate(arg1::Vec{Float64},arg2::IS{Float64},arg3::Vec{Float64},arg4::IS{Float64},arg5::Union(Ptr{VecScatter{Float64}},StridedArray{VecScatter{Float64}},Ptr{Void}))
    ccall((:VecScatterCreate,petsc1),PetscErrorCode,(Vec{Float64},IS{Float64},Vec{Float64},IS{Float64},Ptr{VecScatter{Float64}}),arg1,arg2,arg3,arg4,arg5)
end

function VecScatterCreateEmpty(arg1::MPI_Comm,arg2::Union(Ptr{VecScatter{Float64}},StridedArray{VecScatter{Float64}},Ptr{Void}))
    ccall((:VecScatterCreateEmpty,petsc1),PetscErrorCode,(comm_type,Ptr{VecScatter{Float64}}),arg1.val,arg2)
end

function VecScatterCreateLocal(arg1::VecScatter{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Integer,arg7::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg8::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg9::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg10::Integer)
    ccall((:VecScatterCreateLocal,petsc1),PetscErrorCode,(VecScatter{Float64},Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32},Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32},Int32),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function VecScatterBegin(arg1::VecScatter{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::InsertMode,arg5::ScatterMode)
    ccall((:VecScatterBegin,petsc1),PetscErrorCode,(VecScatter{Float64},Vec{Float64},Vec{Float64},InsertMode,ScatterMode),arg1,arg2,arg3,arg4,arg5)
end

function VecScatterEnd(arg1::VecScatter{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::InsertMode,arg5::ScatterMode)
    ccall((:VecScatterEnd,petsc1),PetscErrorCode,(VecScatter{Float64},Vec{Float64},Vec{Float64},InsertMode,ScatterMode),arg1,arg2,arg3,arg4,arg5)
end

function VecScatterDestroy(arg1::Union(Ptr{VecScatter{Float64}},StridedArray{VecScatter{Float64}},Ptr{Void}))
    ccall((:VecScatterDestroy,petsc1),PetscErrorCode,(Ptr{VecScatter{Float64}},),arg1)
end

function VecScatterCopy(arg1::VecScatter{Float64},arg2::Union(Ptr{VecScatter{Float64}},StridedArray{VecScatter{Float64}},Ptr{Void}))
    ccall((:VecScatterCopy,petsc1),PetscErrorCode,(VecScatter{Float64},Ptr{VecScatter{Float64}}),arg1,arg2)
end

function VecScatterView(arg1::VecScatter{Float64},arg2::PetscViewer{Float64})
    ccall((:VecScatterView,petsc1),PetscErrorCode,(VecScatter{Float64},PetscViewer{Float64}),arg1,arg2)
end

function VecScatterGetMerged(arg1::VecScatter{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:VecScatterGetMerged,petsc1),PetscErrorCode,(VecScatter{Float64},Ptr{PetscBool}),arg1,arg2)
end

function VecGetArray4d(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Integer,arg9::Integer,arg10::Union(Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}},StridedArray{Ptr{Ptr{Ptr{Ptr{Float64}}}}},Ptr{Void}))
    ccall((:VecGetArray4d,petsc1),PetscErrorCode,(Vec{Float64},Int32,Int32,Int32,Int32,Int32,Int32,Int32,Int32,Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function VecRestoreArray4d(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Integer,arg9::Integer,arg10::Union(Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}},StridedArray{Ptr{Ptr{Ptr{Ptr{Float64}}}}},Ptr{Void}))
    ccall((:VecRestoreArray4d,petsc1),PetscErrorCode,(Vec{Float64},Int32,Int32,Int32,Int32,Int32,Int32,Int32,Int32,Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function VecGetArray3d(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Union(Ptr{Ptr{Ptr{Ptr{Float64}}}},StridedArray{Ptr{Ptr{Ptr{Float64}}}},Ptr{Void}))
    ccall((:VecGetArray3d,petsc1),PetscErrorCode,(Vec{Float64},Int32,Int32,Int32,Int32,Int32,Int32,Ptr{Ptr{Ptr{Ptr{Float64}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function VecRestoreArray3d(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Union(Ptr{Ptr{Ptr{Ptr{Float64}}}},StridedArray{Ptr{Ptr{Ptr{Float64}}}},Ptr{Void}))
    ccall((:VecRestoreArray3d,petsc1),PetscErrorCode,(Vec{Float64},Int32,Int32,Int32,Int32,Int32,Int32,Ptr{Ptr{Ptr{Ptr{Float64}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function VecGetArray2d(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Ptr{Ptr{Float64}}},StridedArray{Ptr{Ptr{Float64}}},Ptr{Void}))
    ccall((:VecGetArray2d,petsc1),PetscErrorCode,(Vec{Float64},Int32,Int32,Int32,Int32,Ptr{Ptr{Ptr{Float64}}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecRestoreArray2d(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Ptr{Ptr{Float64}}},StridedArray{Ptr{Ptr{Float64}}},Ptr{Void}))
    ccall((:VecRestoreArray2d,petsc1),PetscErrorCode,(Vec{Float64},Int32,Int32,Int32,Int32,Ptr{Ptr{Ptr{Float64}}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecGetArray1d(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:VecGetArray1d,petsc1),PetscErrorCode,(Vec{Float64},Int32,Int32,Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4)
end

function VecRestoreArray1d(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:VecRestoreArray1d,petsc1),PetscErrorCode,(Vec{Float64},Int32,Int32,Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4)
end

function VecGetArray4dRead(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Integer,arg9::Integer,arg10::Union(Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}},StridedArray{Ptr{Ptr{Ptr{Ptr{Float64}}}}},Ptr{Void}))
    ccall((:VecGetArray4dRead,petsc1),PetscErrorCode,(Vec{Float64},Int32,Int32,Int32,Int32,Int32,Int32,Int32,Int32,Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function VecRestoreArray4dRead(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Integer,arg9::Integer,arg10::Union(Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}},StridedArray{Ptr{Ptr{Ptr{Ptr{Float64}}}}},Ptr{Void}))
    ccall((:VecRestoreArray4dRead,petsc1),PetscErrorCode,(Vec{Float64},Int32,Int32,Int32,Int32,Int32,Int32,Int32,Int32,Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function VecGetArray3dRead(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Union(Ptr{Ptr{Ptr{Ptr{Float64}}}},StridedArray{Ptr{Ptr{Ptr{Float64}}}},Ptr{Void}))
    ccall((:VecGetArray3dRead,petsc1),PetscErrorCode,(Vec{Float64},Int32,Int32,Int32,Int32,Int32,Int32,Ptr{Ptr{Ptr{Ptr{Float64}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function VecRestoreArray3dRead(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Union(Ptr{Ptr{Ptr{Ptr{Float64}}}},StridedArray{Ptr{Ptr{Ptr{Float64}}}},Ptr{Void}))
    ccall((:VecRestoreArray3dRead,petsc1),PetscErrorCode,(Vec{Float64},Int32,Int32,Int32,Int32,Int32,Int32,Ptr{Ptr{Ptr{Ptr{Float64}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function VecGetArray2dRead(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Ptr{Ptr{Float64}}},StridedArray{Ptr{Ptr{Float64}}},Ptr{Void}))
    ccall((:VecGetArray2dRead,petsc1),PetscErrorCode,(Vec{Float64},Int32,Int32,Int32,Int32,Ptr{Ptr{Ptr{Float64}}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecRestoreArray2dRead(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Ptr{Ptr{Float64}}},StridedArray{Ptr{Ptr{Float64}}},Ptr{Void}))
    ccall((:VecRestoreArray2dRead,petsc1),PetscErrorCode,(Vec{Float64},Int32,Int32,Int32,Int32,Ptr{Ptr{Ptr{Float64}}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecGetArray1dRead(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:VecGetArray1dRead,petsc1),PetscErrorCode,(Vec{Float64},Int32,Int32,Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4)
end

function VecRestoreArray1dRead(arg1::Vec{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:VecRestoreArray1dRead,petsc1),PetscErrorCode,(Vec{Float64},Int32,Int32,Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4)
end

function VecPlaceArray(arg1::Vec{Float64},arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:VecPlaceArray,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Float64}),arg1,arg2)
end

function VecResetArray(arg1::Vec{Float64})
    ccall((:VecResetArray,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
end

function VecReplaceArray(arg1::Vec{Float64},arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:VecReplaceArray,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Float64}),arg1,arg2)
end

function VecGetArrays(arg1::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg2::Integer,arg3::Union(Ptr{Ptr{Ptr{Float64}}},StridedArray{Ptr{Ptr{Float64}}},Ptr{Void}))
    ccall((:VecGetArrays,petsc1),PetscErrorCode,(Ptr{Vec{Float64}},Int32,Ptr{Ptr{Ptr{Float64}}}),arg1,arg2,arg3)
end

function VecRestoreArrays(arg1::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg2::Integer,arg3::Union(Ptr{Ptr{Ptr{Float64}}},StridedArray{Ptr{Ptr{Float64}}},Ptr{Void}))
    ccall((:VecRestoreArrays,petsc1),PetscErrorCode,(Ptr{Vec{Float64}},Int32,Ptr{Ptr{Ptr{Float64}}}),arg1,arg2,arg3)
end

function VecView(arg1::Vec{Float64},arg2::PetscViewer{Float64})
    ccall((:VecView,petsc1),PetscErrorCode,(Vec{Float64},PetscViewer{Float64}),arg1,arg2)
end

function VecEqual(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:VecEqual,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{PetscBool}),arg1,arg2,arg3)
end

function VecLoad(arg1::Vec{Float64},arg2::PetscViewer{Float64})
    ccall((:VecLoad,petsc1),PetscErrorCode,(Vec{Float64},PetscViewer{Float64}),arg1,arg2)
end

function VecGetSize(arg1::Vec{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:VecGetSize,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Int32}),arg1,arg2)
end

function VecGetLocalSize(arg1::Vec{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:VecGetLocalSize,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Int32}),arg1,arg2)
end

function VecGetOwnershipRange(arg1::Vec{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:VecGetOwnershipRange,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3)
end

function VecGetOwnershipRanges(arg1::Vec{Float64},arg2::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:VecGetOwnershipRanges,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Ptr{Int32}}),arg1,arg2)
end

function VecSetLocalToGlobalMapping(arg1::Vec{Float64},arg2::ISLocalToGlobalMapping{Float64})
    ccall((:VecSetLocalToGlobalMapping,petsc1),PetscErrorCode,(Vec{Float64},ISLocalToGlobalMapping{Float64}),arg1,arg2)
end

function VecSetValuesLocal(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::InsertMode)
    ccall((:VecSetValuesLocal,petsc1),PetscErrorCode,(Vec{Float64},Int32,Ptr{Int32},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5)
end

function VecGetLocalToGlobalMapping(arg1::Vec{Float64},arg2::Union(Ptr{ISLocalToGlobalMapping{Float64}},StridedArray{ISLocalToGlobalMapping{Float64}},Ptr{Void}))
    ccall((:VecGetLocalToGlobalMapping,petsc1),PetscErrorCode,(Vec{Float64},Ptr{ISLocalToGlobalMapping{Float64}}),arg1,arg2)
end

function VecDotBegin(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:VecDotBegin,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{Float64}),arg1,arg2,arg3)
end

function VecDotEnd(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:VecDotEnd,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{Float64}),arg1,arg2,arg3)
end

function VecTDotBegin(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:VecTDotBegin,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{Float64}),arg1,arg2,arg3)
end

function VecTDotEnd(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:VecTDotEnd,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{Float64}),arg1,arg2,arg3)
end

function VecNormBegin(arg1::Vec{Float64},arg2::NormType,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:VecNormBegin,petsc1),PetscErrorCode,(Vec{Float64},NormType,Ptr{Cint}),arg1,arg2,arg3)
end

function VecNormEnd(arg1::Vec{Float64},arg2::NormType,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:VecNormEnd,petsc1),PetscErrorCode,(Vec{Float64},NormType,Ptr{Cint}),arg1,arg2,arg3)
end

function VecMDotBegin(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:VecMDotBegin,petsc1),PetscErrorCode,(Vec{Float64},Int32,Ptr{Vec{Float64}},Ptr{Float64}),arg1,arg2,arg3,arg4)
end

function VecMDotEnd(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:VecMDotEnd,petsc1),PetscErrorCode,(Vec{Float64},Int32,Ptr{Vec{Float64}},Ptr{Float64}),arg1,arg2,arg3,arg4)
end

function VecMTDotBegin(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:VecMTDotBegin,petsc1),PetscErrorCode,(Vec{Float64},Int32,Ptr{Vec{Float64}},Ptr{Float64}),arg1,arg2,arg3,arg4)
end

function VecMTDotEnd(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:VecMTDotEnd,petsc1),PetscErrorCode,(Vec{Float64},Int32,Ptr{Vec{Float64}},Ptr{Float64}),arg1,arg2,arg3,arg4)
end

function PetscCommSplitReductionBegin(arg0::Type{Float64},arg1::MPI_Comm)
    ccall((:PetscCommSplitReductionBegin,petsc1),PetscErrorCode,(comm_type,),arg1.val)
end

function VecSetOption(arg1::Vec{Float64},arg2::VecOption,arg3::PetscBool)
    ccall((:VecSetOption,petsc1),PetscErrorCode,(Vec{Float64},VecOption,PetscBool),arg1,arg2,arg3)
end

function VecGetArray(arg1::Vec{Float64},arg2::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:VecGetArray,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Ptr{Float64}}),arg1,arg2)
end

function VecGetArrayRead(arg1::Vec{Float64},arg2::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:VecGetArrayRead,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Ptr{Float64}}),arg1,arg2)
end

function VecRestoreArray(arg1::Vec{Float64},arg2::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:VecRestoreArray,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Ptr{Float64}}),arg1,arg2)
end

function VecRestoreArrayRead(arg1::Vec{Float64},arg2::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:VecRestoreArrayRead,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Ptr{Float64}}),arg1,arg2)
end

function VecGetLocalVector(arg1::Vec{Float64},arg2::Vec{Float64})
    ccall((:VecGetLocalVector,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64}),arg1,arg2)
end

function VecRestoreLocalVector(arg1::Vec{Float64},arg2::Vec{Float64})
    ccall((:VecRestoreLocalVector,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64}),arg1,arg2)
end

function VecGetLocalVectorRead(arg1::Vec{Float64},arg2::Vec{Float64})
    ccall((:VecGetLocalVectorRead,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64}),arg1,arg2)
end

function VecRestoreLocalVectorRead(arg1::Vec{Float64},arg2::Vec{Float64})
    ccall((:VecRestoreLocalVectorRead,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64}),arg1,arg2)
end

function VecContourScale(arg1::Vec{Float64},Float64::Integer,arg2::Integer)
    ccall((:VecContourScale,petsc1),PetscErrorCode,(Vec{Float64},Cint,Cint),arg1,PetscReal,arg2)
end

function VecSetOperation(arg1::Vec{Float64},arg2::VecOperation,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:VecSetOperation,petsc1),PetscErrorCode,(Vec{Float64},VecOperation,Ptr{Void}),arg1,arg2,arg3)
end

function VecMPISetGhost(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:VecMPISetGhost,petsc1),PetscErrorCode,(Vec{Float64},Int32,Ptr{Int32}),arg1,arg2,arg3)
end

function VecCreateGhost(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:VecCreateGhost,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Ptr{Int32},Ptr{Vec{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6)
end

function VecCreateGhostWithArray(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:VecCreateGhostWithArray,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Ptr{Int32},Ptr{Float64},Ptr{Vec{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
end

function VecCreateGhostBlock(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:VecCreateGhostBlock,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Int32,Ptr{Int32},Ptr{Vec{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
end

function VecCreateGhostBlockWithArray(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg8::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:VecCreateGhostBlockWithArray,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Int32,Ptr{Int32},Ptr{Float64},Ptr{Vec{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function VecGhostGetLocalForm(arg1::Vec{Float64},arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:VecGhostGetLocalForm,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Vec{Float64}}),arg1,arg2)
end

function VecGhostRestoreLocalForm(arg1::Vec{Float64},arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:VecGhostRestoreLocalForm,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Vec{Float64}}),arg1,arg2)
end

function VecGhostIsLocalForm(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:VecGhostIsLocalForm,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{PetscBool}),arg1,arg2,arg3)
end

function VecGhostUpdateBegin(arg1::Vec{Float64},arg2::InsertMode,arg3::ScatterMode)
    ccall((:VecGhostUpdateBegin,petsc1),PetscErrorCode,(Vec{Float64},InsertMode,ScatterMode),arg1,arg2,arg3)
end

function VecGhostUpdateEnd(arg1::Vec{Float64},arg2::InsertMode,arg3::ScatterMode)
    ccall((:VecGhostUpdateEnd,petsc1),PetscErrorCode,(Vec{Float64},InsertMode,ScatterMode),arg1,arg2,arg3)
end

function VecConjugate(arg1::Vec{Float64})
    ccall((:VecConjugate,petsc1),PetscErrorCode,(Vec{Float64},),arg1)
end

function VecScatterCreateToAll(arg1::Vec{Float64},arg2::Union(Ptr{VecScatter{Float64}},StridedArray{VecScatter{Float64}},Ptr{Void}),arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:VecScatterCreateToAll,petsc1),PetscErrorCode,(Vec{Float64},Ptr{VecScatter{Float64}},Ptr{Vec{Float64}}),arg1,arg2,arg3)
end

function VecScatterCreateToZero(arg1::Vec{Float64},arg2::Union(Ptr{VecScatter{Float64}},StridedArray{VecScatter{Float64}},Ptr{Void}),arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:VecScatterCreateToZero,petsc1),PetscErrorCode,(Vec{Float64},Ptr{VecScatter{Float64}},Ptr{Vec{Float64}}),arg1,arg2,arg3)
end

function ISComplementVec(arg1::IS{Float64},arg2::Vec{Float64},arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:ISComplementVec,petsc1),PetscErrorCode,(IS{Float64},Vec{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3)
end

function VecPow(arg1::Vec{Float64},arg2::Float64)
    ccall((:VecPow,petsc1),PetscErrorCode,(Vec{Float64},Float64),arg1,arg2)
end

function VecMedian(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64})
    ccall((:VecMedian,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
end

function VecWhichBetween(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:VecWhichBetween,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4)
end

function VecWhichBetweenOrEqual(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:VecWhichBetweenOrEqual,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4)
end

function VecWhichGreaterThan(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:VecWhichGreaterThan,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3)
end

function VecWhichLessThan(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:VecWhichLessThan,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3)
end

function VecWhichEqual(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:VecWhichEqual,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{IS{Float64}}),arg1,arg2,arg3)
end

function VecISAXPY(arg1::Vec{Float64},arg2::IS{Float64},arg3::Float64,arg4::Vec{Float64})
    ccall((:VecISAXPY,petsc1),PetscErrorCode,(Vec{Float64},IS{Float64},Float64,Vec{Float64}),arg1,arg2,arg3,arg4)
end

function VecISSet(arg1::Vec{Float64},arg2::IS{Float64},arg3::Float64)
    ccall((:VecISSet,petsc1),PetscErrorCode,(Vec{Float64},IS{Float64},Float64),arg1,arg2,arg3)
end

function VecBoundGradientProjection(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64},arg5::Vec{Float64})
    ccall((:VecBoundGradientProjection,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5)
end

function VecStepBoundInfo(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64},arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg7::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:VecStepBoundInfo,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function VecStepMax(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:VecStepMax,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Ptr{Cint}),arg1,arg2,arg3)
end

function VecStepMaxBounded(arg1::Vec{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64},arg5::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:VecStepMaxBounded,petsc1),PetscErrorCode,(Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
end

function PetscViewerMathematicaGetVector(arg1::PetscViewer{Float64},arg2::Vec{Float64})
    ccall((:PetscViewerMathematicaGetVector,petsc1),PetscErrorCode,(PetscViewer{Float64},Vec{Float64}),arg1,arg2)
end

function PetscViewerMathematicaPutVector(arg1::PetscViewer{Float64},arg2::Vec{Float64})
    ccall((:PetscViewerMathematicaPutVector,petsc1),PetscErrorCode,(PetscViewer{Float64},Vec{Float64}),arg1,arg2)
end

#= skipping function with undefined symbols: 
 function VecsDestroy(arg0::Type{Float64},arg1::Vecs)
    ccall((:VecsDestroy,petsc1),PetscErrorCode,(Vecs,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function VecsCreateSeq(arg0::Type{Float64},arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Vecs},StridedArray{Vecs},Ptr{Void}))
    ccall((:VecsCreateSeq,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Ptr{Vecs}),arg1.val,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function VecsCreateSeqWithArray(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::Union(Ptr{Vecs},StridedArray{Vecs},Ptr{Void}))
    ccall((:VecsCreateSeqWithArray,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Ptr{Float64},Ptr{Vecs}),arg1.val,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function VecsDuplicate(arg0::Type{Float64},arg1::Vecs,arg2::Union(Ptr{Vecs},StridedArray{Vecs},Ptr{Void}))
    ccall((:VecsDuplicate,petsc1),PetscErrorCode,(Vecs,Ptr{Vecs}),arg1,arg2)
end 
=#
function VecNestGetSubVecs(arg1::Vec{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Ptr{Vec{Float64}}},StridedArray{Ptr{Vec{Float64}}},Ptr{Void}))
    ccall((:VecNestGetSubVecs,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Int32},Ptr{Ptr{Vec{Float64}}}),arg1,arg2,arg3)
end

function VecNestGetSubVec(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:VecNestGetSubVec,petsc1),PetscErrorCode,(Vec{Float64},Int32,Ptr{Vec{Float64}}),arg1,arg2,arg3)
end

function VecNestSetSubVecs(arg1::Vec{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:VecNestSetSubVecs,petsc1),PetscErrorCode,(Vec{Float64},Int32,Ptr{Int32},Ptr{Vec{Float64}}),arg1,arg2,arg3,arg4)
end

function VecNestSetSubVec(arg1::Vec{Float64},arg2::Integer,arg3::Vec{Float64})
    ccall((:VecNestSetSubVec,petsc1),PetscErrorCode,(Vec{Float64},Int32,Vec{Float64}),arg1,arg2,arg3)
end

function VecCreateNest(arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg5::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:VecCreateNest,petsc1),PetscErrorCode,(comm_type,Int32,Ptr{IS{Float64}},Ptr{Vec{Float64}},Ptr{Vec{Float64}}),arg1.val,arg2,arg3,arg4,arg5)
end

function VecNestGetSize(arg1::Vec{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:VecNestGetSize,petsc1),PetscErrorCode,(Vec{Float64},Ptr{Int32}),arg1,arg2)
end

function PetscOptionsGetVec(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Vec{Float64},arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PetscOptionsGetVec,petsc1),PetscErrorCode,(Cstring,Cstring,Vec{Float64},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function VecChop(arg1::Vec{Float64},Float64::Integer)
    ccall((:VecChop,petsc1),PetscErrorCode,(Vec{Float64},Cint),arg1,PetscReal)
end

function VecGetLayout(arg1::Vec{Float64},arg2::Union(Ptr{PetscLayout{Float64}},StridedArray{PetscLayout{Float64}},Ptr{Void}))
    ccall((:VecGetLayout,petsc1),PetscErrorCode,(Vec{Float64},Ptr{PetscLayout{Float64}{Float64}}),arg1,arg2)
end

function VecSetLayout(arg1::Vec{Float64},arg2::PetscLayout{Float64})
    ccall((:VecSetLayout,petsc1),PetscErrorCode,(Vec{Float64},PetscLayout{Float64}),arg1,arg2)
end

#= skipping function with undefined symbols: 
 function PetscSectionVecView(arg1::PetscSection,arg2::Vec{Float64},arg3::PetscViewer{Float64})
    ccall((:PetscSectionVecView,petsc1),PetscErrorCode,(PetscSection,Vec{Float64},PetscViewer{Float64}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function VecGetValuesSection(arg1::Vec{Float64},arg2::PetscSection,arg3::Integer,arg4::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:VecGetValuesSection,petsc1),PetscErrorCode,(Vec{Float64},PetscSection,Int32,Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function VecSetValuesSection(arg1::Vec{Float64},arg2::PetscSection,arg3::Integer,arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::InsertMode)
    ccall((:VecSetValuesSection,petsc1),PetscErrorCode,(Vec{Float64},PetscSection,Int32,Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5)
end 
=#
#= skipping function with undefined symbols: 
 function PetscSectionVecNorm(arg1::PetscSection,arg2::PetscSection,arg3::Vec{Float64},arg4::NormType,Float64::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscSectionVecNorm,petsc1),PetscErrorCode,(PetscSection,PetscSection,Vec{Float64},NormType,Ptr{Cint}),arg1,arg2,arg3,arg4,PetscReal)
end 
=#
function MatGetFactor(arg1::Mat{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::MatFactorType,arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatGetFactor,petsc1),PetscErrorCode,(Mat{Float64},Cstring,MatFactorType,Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4)
end

function MatGetFactorAvailable(arg1::Mat{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::MatFactorType,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:MatGetFactorAvailable,petsc1),PetscErrorCode,(Mat{Float64},Cstring,MatFactorType,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function MatFactorGetSolverPackage(arg1::Mat{Float64},arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:MatFactorGetSolverPackage,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Ptr{UInt8}}),arg1,arg2)
end

function MatGetFactorType(arg1::Mat{Float64},arg2::Union(Ptr{MatFactorType},StridedArray{MatFactorType},Ptr{Void}))
    ccall((:MatGetFactorType,petsc1),PetscErrorCode,(Mat{Float64},Ptr{MatFactorType}),arg1,arg2)
end

function MatSolverPackageRegister(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::MatType,arg3::MatFactorType,arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:MatSolverPackageRegister,petsc1),PetscErrorCode,(Cstring,Cstring,MatFactorType,Ptr{Void}),arg1,arg2,arg3,arg4)
end

function MatSolverPackageGet(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::MatType,arg3::MatFactorType,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg5::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg6::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:MatSolverPackageGet,petsc1),PetscErrorCode,(Cstring,Cstring,MatFactorType,Ptr{PetscBool},Ptr{PetscBool},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatInitializePackage(arg0::Type{Float64})
    ccall((:MatInitializePackage,petsc1),PetscErrorCode,())
end

function MatCreate(arg1::MPI_Comm,arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreate,petsc1),PetscErrorCode,(comm_type,Ptr{Mat{Float64}}),arg1.val,arg2)
end

function MatSetSizes(arg1::Mat{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer)
    ccall((:MatSetSizes,petsc1),PetscErrorCode,(Mat{Float64},Int32,Int32,Int32,Int32),arg1,arg2,arg3,arg4,arg5)
end

function MatSetType(arg1::Mat{Float64},arg2::MatType)
    ccall((:MatSetType,petsc1),PetscErrorCode,(Mat{Float64},Cstring),arg1,arg2)
end

function MatSetFromOptions(arg1::Mat{Float64})
    ccall((:MatSetFromOptions,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
end

function MatRegisterBaseName(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:MatRegisterBaseName,petsc1),PetscErrorCode,(Cstring,Cstring,Cstring),arg1,arg2,arg3)
end

function MatSetOptionsPrefix(arg1::Mat{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:MatSetOptionsPrefix,petsc1),PetscErrorCode,(Mat{Float64},Cstring),arg1,arg2)
end

function MatAppendOptionsPrefix(arg1::Mat{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:MatAppendOptionsPrefix,petsc1),PetscErrorCode,(Mat{Float64},Cstring),arg1,arg2)
end

function MatGetOptionsPrefix(arg1::Mat{Float64},arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:MatGetOptionsPrefix,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Ptr{UInt8}}),arg1,arg2)
end

function MatSetErrorIfFPE(arg1::Mat{Float64},arg2::PetscBool)
    ccall((:MatSetErrorIfFPE,petsc1),PetscErrorCode,(Mat{Float64},PetscBool),arg1,arg2)
end

function MatCreateSeqDense(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg5::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateSeqDense,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Ptr{Float64},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5)
end

function MatCreateDense(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateDense,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Int32,Ptr{Float64},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateSeqAIJ(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateSeqAIJ,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Ptr{Int32},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6)
end

function MatCreateAIJ(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg8::Integer,arg9::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg10::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateAIJ,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Int32,Int32,Ptr{Int32},Int32,Ptr{Int32},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function MatCreateMPIAIJWithArrays(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg8::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg9::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateMPIAIJWithArrays,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function MatCreateMPIAIJWithSplitArrays(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg8::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg9::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg10::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg11::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg12::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateMPIAIJWithSplitArrays,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12)
end

function MatCreateSeqBAIJ(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateSeqBAIJ,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Int32,Ptr{Int32},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateBAIJ(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg9::Integer,arg10::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg11::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateBAIJ,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Int32,Int32,Int32,Ptr{Int32},Int32,Ptr{Int32},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
end

function MatCreateMPIBAIJWithArrays(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg8::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg9::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg10::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateMPIBAIJWithArrays,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Int32,Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function MatCreateMPIAdj(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateMPIAdj,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateSeqSBAIJ(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateSeqSBAIJ,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Int32,Ptr{Int32},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateSBAIJ(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg9::Integer,arg10::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg11::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateSBAIJ,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Int32,Int32,Int32,Ptr{Int32},Int32,Ptr{Int32},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
end

function MatCreateMPISBAIJWithArrays(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg8::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg9::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg10::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateMPISBAIJWithArrays,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Int32,Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function MatSeqSBAIJSetPreallocationCSR(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:MatSeqSBAIJSetPreallocationCSR,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5)
end

function MatMPISBAIJSetPreallocationCSR(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:MatMPISBAIJSetPreallocationCSR,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5)
end

function MatXAIJSetPreallocation(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatXAIJSetPreallocation,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatCreateShell(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateShell,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Int32,Ptr{Void},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateNormal(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateNormal,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2)
end

function MatCreateLRC(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64},arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateLRC,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4)
end

function MatCreateIS(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::ISLocalToGlobalMapping{Float64},arg8::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateIS,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Int32,Int32,ISLocalToGlobalMapping{Float64},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatCreateSeqAIJCRL(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateSeqAIJCRL,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Ptr{Int32},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6)
end

function MatCreateMPIAIJCRL(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Integer,arg7::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg8::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateMPIAIJCRL,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Ptr{Int32},Int32,Ptr{Int32},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatCreateSeqBSTRM(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateSeqBSTRM,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Int32,Ptr{Int32},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateMPIBSTRM(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg9::Integer,arg10::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg11::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateMPIBSTRM,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Int32,Int32,Int32,Ptr{Int32},Int32,Ptr{Int32},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
end

function MatCreateSeqSBSTRM(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateSeqSBSTRM,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Int32,Ptr{Int32},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateMPISBSTRM(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg9::Integer,arg10::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg11::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateMPISBSTRM,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Int32,Int32,Int32,Ptr{Int32},Int32,Ptr{Int32},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
end

function MatCreateScatter(arg1::MPI_Comm,arg2::VecScatter{Float64},arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateScatter,petsc1),PetscErrorCode,(comm_type,VecScatter{Float64},Ptr{Mat{Float64}}),arg1.val,arg2,arg3)
end

function MatScatterSetVecScatter(arg1::Mat{Float64},arg2::VecScatter{Float64})
    ccall((:MatScatterSetVecScatter,petsc1),PetscErrorCode,(Mat{Float64},VecScatter{Float64}),arg1,arg2)
end

function MatScatterGetVecScatter(arg1::Mat{Float64},arg2::Union(Ptr{VecScatter{Float64}},StridedArray{VecScatter{Float64}},Ptr{Void}))
    ccall((:MatScatterGetVecScatter,petsc1),PetscErrorCode,(Mat{Float64},Ptr{VecScatter{Float64}}),arg1,arg2)
end

function MatCreateBlockMat(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateBlockMat,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Int32,Ptr{Int32},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCompositeAddMat(arg1::Mat{Float64},arg2::Mat{Float64})
    ccall((:MatCompositeAddMat,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64}),arg1,arg2)
end

function MatCompositeMerge(arg1::Mat{Float64})
    ccall((:MatCompositeMerge,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
end

function MatCreateComposite(arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateComposite,petsc1),PetscErrorCode,(comm_type,Int32,Ptr{Mat{Float64}},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4)
end

function MatCompositeSetType(arg1::Mat{Float64},arg2::MatCompositeType)
    ccall((:MatCompositeSetType,petsc1),PetscErrorCode,(Mat{Float64},MatCompositeType),arg1,arg2)
end

function MatCreateFFT(arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::MatType,arg5::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateFFT,petsc1),PetscErrorCode,(comm_type,Int32,Ptr{Int32},Cstring,Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5)
end

function MatCreateSeqCUFFT(arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateSeqCUFFT,petsc1),PetscErrorCode,(comm_type,Int32,Ptr{Int32},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4)
end

function MatCreateTranspose(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateTranspose,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2)
end

function MatCreateHermitianTranspose(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateHermitianTranspose,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2)
end

function MatCreateSubMatrix(arg1::Mat{Float64},arg2::IS{Float64},arg3::IS{Float64},arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateSubMatrix,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},IS{Float64},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4)
end

function MatSubMatrixUpdate(arg1::Mat{Float64},arg2::Mat{Float64},arg3::IS{Float64},arg4::IS{Float64})
    ccall((:MatSubMatrixUpdate,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},IS{Float64},IS{Float64}),arg1,arg2,arg3,arg4)
end

function MatCreateLocalRef(arg1::Mat{Float64},arg2::IS{Float64},arg3::IS{Float64},arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateLocalRef,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},IS{Float64},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4)
end

function MatPythonSetType(arg1::Mat{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:MatPythonSetType,petsc1),PetscErrorCode,(Mat{Float64},Cstring),arg1,arg2)
end

function MatSetUp(arg1::Mat{Float64})
    ccall((:MatSetUp,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
end

function MatDestroy(arg1::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatDestroy,petsc1),PetscErrorCode,(Ptr{Mat{Float64}},),arg1)
end

function MatGetNonzeroState(arg1::Mat{Float64},arg2::Union(Ptr{PetscObjectState},StridedArray{PetscObjectState},Ptr{Void}))
    ccall((:MatGetNonzeroState,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscObjectState}),arg1,arg2)
end

function MatConjugate(arg1::Mat{Float64})
    ccall((:MatConjugate,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
end

function MatRealPart(arg1::Mat{Float64})
    ccall((:MatRealPart,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
end

function MatImaginaryPart(arg1::Mat{Float64})
    ccall((:MatImaginaryPart,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
end

function MatGetDiagonalBlock(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatGetDiagonalBlock,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2)
end

function MatGetTrace(arg1::Mat{Float64},arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:MatGetTrace,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Float64}),arg1,arg2)
end

function MatInvertBlockDiagonal(arg1::Mat{Float64},arg2::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:MatInvertBlockDiagonal,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Ptr{Float64}}),arg1,arg2)
end

function MatSetValues(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::InsertMode)
    ccall((:MatSetValues,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Int32,Ptr{Int32},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatSetValuesBlocked(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::InsertMode)
    ccall((:MatSetValuesBlocked,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Int32,Ptr{Int32},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatSetValuesRow(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:MatSetValuesRow,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Float64}),arg1,arg2,arg3)
end

function MatSetValuesRowLocal(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:MatSetValuesRowLocal,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Float64}),arg1,arg2,arg3)
end

function MatSetValuesBatch(arg1::Mat{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:MatSetValuesBatch,petsc1),PetscErrorCode,(Mat{Float64},Int32,Int32,Ptr{Int32},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5)
end

#= skipping function with undefined symbols: 
 function MatSetRandom(arg1::Mat{Float64},arg2::PetscRandom)
    ccall((:MatSetRandom,petsc1),PetscErrorCode,(Mat{Float64},PetscRandom),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatSetValuesStencil(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{MatStencil},StridedArray{MatStencil},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{MatStencil},StridedArray{MatStencil},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::InsertMode)
    ccall((:MatSetValuesStencil,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{MatStencil},Int32,Ptr{MatStencil},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function MatSetValuesBlockedStencil(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{MatStencil},StridedArray{MatStencil},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{MatStencil},StridedArray{MatStencil},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::InsertMode)
    ccall((:MatSetValuesBlockedStencil,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{MatStencil},Int32,Ptr{MatStencil},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
function MatSetStencil(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Integer)
    ccall((:MatSetStencil,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Ptr{Int32},Int32),arg1,arg2,arg3,arg4,arg5)
end

function MatSetColoring(arg1::Mat{Float64},arg2::ISColoring{Float64})
    ccall((:MatSetColoring,petsc1),PetscErrorCode,(Mat{Float64},ISColoring{Float64}),arg1,arg2)
end

function MatSetValuesAdifor(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:MatSetValuesAdifor,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Void}),arg1,arg2,arg3)
end

function MatAssemblyBegin(arg1::Mat{Float64},arg2::MatAssemblyType)
    ccall((:MatAssemblyBegin,petsc1),PetscErrorCode,(Mat{Float64},MatAssemblyType),arg1,arg2)
end

function MatAssemblyEnd(arg1::Mat{Float64},arg2::MatAssemblyType)
    ccall((:MatAssemblyEnd,petsc1),PetscErrorCode,(Mat{Float64},MatAssemblyType),arg1,arg2)
end

function MatAssembled(arg1::Mat{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:MatAssembled,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscBool}),arg1,arg2)
end

function MatSetOption(arg1::Mat{Float64},arg2::MatOption,arg3::PetscBool)
    ccall((:MatSetOption,petsc1),PetscErrorCode,(Mat{Float64},MatOption,PetscBool),arg1,arg2,arg3)
end

function MatGetOption(arg1::Mat{Float64},arg2::MatOption,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:MatGetOption,petsc1),PetscErrorCode,(Mat{Float64},MatOption,Ptr{PetscBool}),arg1,arg2,arg3)
end

function MatGetType(arg1::Mat{Float64},arg2::Union(Ptr{MatType},StridedArray{MatType},Ptr{Void}))
    ccall((:MatGetType,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Cstring}),arg1,arg2)
end

function MatGetValues(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:MatGetValues,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Int32,Ptr{Int32},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatGetRow(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg5::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:MatGetRow,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Ptr{Ptr{Int32}},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4,arg5)
end

function MatRestoreRow(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg5::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:MatRestoreRow,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Ptr{Ptr{Int32}},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4,arg5)
end

function MatGetRowUpperTriangular(arg1::Mat{Float64})
    ccall((:MatGetRowUpperTriangular,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
end

function MatRestoreRowUpperTriangular(arg1::Mat{Float64})
    ccall((:MatRestoreRowUpperTriangular,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
end

function MatGetColumn(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg5::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:MatGetColumn,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Ptr{Ptr{Int32}},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4,arg5)
end

function MatRestoreColumn(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg5::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:MatRestoreColumn,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Ptr{Ptr{Int32}},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4,arg5)
end

function MatGetColumnVector(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Integer)
    ccall((:MatGetColumnVector,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Int32),arg1,arg2,arg3)
end

function MatSeqAIJGetArray(arg1::Mat{Float64},arg2::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:MatSeqAIJGetArray,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Ptr{Float64}}),arg1,arg2)
end

function MatSeqAIJRestoreArray(arg1::Mat{Float64},arg2::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:MatSeqAIJRestoreArray,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Ptr{Float64}}),arg1,arg2)
end

function MatSeqAIJGetMaxRowNonzeros(arg1::Mat{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatSeqAIJGetMaxRowNonzeros,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Int32}),arg1,arg2)
end

function MatSeqAIJSetValuesLocalFast(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::InsertMode)
    ccall((:MatSeqAIJSetValuesLocalFast,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Int32,Ptr{Int32},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatDenseGetArray(arg1::Mat{Float64},arg2::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:MatDenseGetArray,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Ptr{Float64}}),arg1,arg2)
end

function MatDenseRestoreArray(arg1::Mat{Float64},arg2::Union(Ptr{Ptr{Float64}},StridedArray{Ptr{Float64}},Ptr{Void}))
    ccall((:MatDenseRestoreArray,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Ptr{Float64}}),arg1,arg2)
end

function MatGetBlockSize(arg1::Mat{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatGetBlockSize,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Int32}),arg1,arg2)
end

function MatSetBlockSize(arg1::Mat{Float64},arg2::Integer)
    ccall((:MatSetBlockSize,petsc1),PetscErrorCode,(Mat{Float64},Int32),arg1,arg2)
end

function MatGetBlockSizes(arg1::Mat{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatGetBlockSizes,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3)
end

function MatSetBlockSizes(arg1::Mat{Float64},arg2::Integer,arg3::Integer)
    ccall((:MatSetBlockSizes,petsc1),PetscErrorCode,(Mat{Float64},Int32,Int32),arg1,arg2,arg3)
end

function MatSetBlockSizesFromMats(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64})
    ccall((:MatSetBlockSizesFromMats,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3)
end

function MatSetNThreads(arg1::Mat{Float64},arg2::Integer)
    ccall((:MatSetNThreads,petsc1),PetscErrorCode,(Mat{Float64},Int32),arg1,arg2)
end

function MatGetNThreads(arg1::Mat{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatGetNThreads,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Int32}),arg1,arg2)
end

function MatMult(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:MatMult,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function MatMultDiagonalBlock(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:MatMultDiagonalBlock,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function MatMultAdd(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64})
    ccall((:MatMultAdd,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
end

function MatMultTranspose(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:MatMultTranspose,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function MatMultHermitianTranspose(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:MatMultHermitianTranspose,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function MatIsTranspose(arg1::Mat{Float64},arg2::Mat{Float64},Float64::Integer,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:MatIsTranspose,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Cint,Ptr{PetscBool}),arg1,arg2,PetscReal,arg3)
end

function MatIsHermitianTranspose(arg1::Mat{Float64},arg2::Mat{Float64},Float64::Integer,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:MatIsHermitianTranspose,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Cint,Ptr{PetscBool}),arg1,arg2,PetscReal,arg3)
end

function MatMultTransposeAdd(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64})
    ccall((:MatMultTransposeAdd,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
end

function MatMultHermitianTransposeAdd(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64})
    ccall((:MatMultHermitianTransposeAdd,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
end

function MatMultConstrained(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:MatMultConstrained,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function MatMultTransposeConstrained(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:MatMultTransposeConstrained,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function MatMatSolve(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64})
    ccall((:MatMatSolve,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3)
end

function MatResidual(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64})
    ccall((:MatResidual,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
end

function MatConvert(arg1::Mat{Float64},arg2::MatType,arg3::MatReuse,arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatConvert,petsc1),PetscErrorCode,(Mat{Float64},Cstring,MatReuse,Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4)
end

function MatDuplicate(arg1::Mat{Float64},arg2::MatDuplicateOption,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatDuplicate,petsc1),PetscErrorCode,(Mat{Float64},MatDuplicateOption,Ptr{Mat{Float64}}),arg1,arg2,arg3)
end

function MatCopy(arg1::Mat{Float64},arg2::Mat{Float64},arg3::MatStructure)
    ccall((:MatCopy,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},MatStructure),arg1,arg2,arg3)
end

function MatView(arg1::Mat{Float64},arg2::PetscViewer{Float64})
    ccall((:MatView,petsc1),PetscErrorCode,(Mat{Float64},PetscViewer{Float64}),arg1,arg2)
end

function MatIsSymmetric(arg1::Mat{Float64},Float64::Integer,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:MatIsSymmetric,petsc1),PetscErrorCode,(Mat{Float64},Cint,Ptr{PetscBool}),arg1,PetscReal,arg2)
end

function MatIsStructurallySymmetric(arg1::Mat{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:MatIsStructurallySymmetric,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscBool}),arg1,arg2)
end

function MatIsHermitian(arg1::Mat{Float64},Float64::Integer,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:MatIsHermitian,petsc1),PetscErrorCode,(Mat{Float64},Cint,Ptr{PetscBool}),arg1,PetscReal,arg2)
end

function MatIsSymmetricKnown(arg1::Mat{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:MatIsSymmetricKnown,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3)
end

function MatIsHermitianKnown(arg1::Mat{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:MatIsHermitianKnown,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3)
end

function MatMissingDiagonal(arg1::Mat{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatMissingDiagonal,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscBool},Ptr{Int32}),arg1,arg2,arg3)
end

function MatLoad(arg1::Mat{Float64},arg2::PetscViewer{Float64})
    ccall((:MatLoad,petsc1),PetscErrorCode,(Mat{Float64},PetscViewer{Float64}),arg1,arg2)
end

function MatGetRowIJ(arg1::Mat{Float64},arg2::Integer,arg3::PetscBool,arg4::PetscBool,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg7::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg8::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:MatGetRowIJ,petsc1),PetscErrorCode,(Mat{Float64},Int32,PetscBool,PetscBool,Ptr{Int32},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatRestoreRowIJ(arg1::Mat{Float64},arg2::Integer,arg3::PetscBool,arg4::PetscBool,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg7::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg8::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:MatRestoreRowIJ,petsc1),PetscErrorCode,(Mat{Float64},Int32,PetscBool,PetscBool,Ptr{Int32},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatGetColumnIJ(arg1::Mat{Float64},arg2::Integer,arg3::PetscBool,arg4::PetscBool,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg7::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg8::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:MatGetColumnIJ,petsc1),PetscErrorCode,(Mat{Float64},Int32,PetscBool,PetscBool,Ptr{Int32},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatRestoreColumnIJ(arg1::Mat{Float64},arg2::Integer,arg3::PetscBool,arg4::PetscBool,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg7::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg8::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:MatRestoreColumnIJ,petsc1),PetscErrorCode,(Mat{Float64},Int32,PetscBool,PetscBool,Ptr{Int32},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatGetInfo(arg1::Mat{Float64},arg2::MatInfoType,arg3::Union(Ptr{MatInfo},StridedArray{MatInfo},Ptr{Void}))
    ccall((:MatGetInfo,petsc1),PetscErrorCode,(Mat{Float64},MatInfoType,Ptr{MatInfo}),arg1,arg2,arg3)
end

function MatGetDiagonal(arg1::Mat{Float64},arg2::Vec{Float64})
    ccall((:MatGetDiagonal,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64}),arg1,arg2)
end

function MatGetRowMax(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatGetRowMax,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Ptr{Int32}),arg1,arg2,arg3)
end

function MatGetRowMin(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatGetRowMin,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Ptr{Int32}),arg1,arg2,arg3)
end

function MatGetRowMaxAbs(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatGetRowMaxAbs,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Ptr{Int32}),arg1,arg2,arg3)
end

function MatGetRowMinAbs(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatGetRowMinAbs,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Ptr{Int32}),arg1,arg2,arg3)
end

function MatGetRowSum(arg1::Mat{Float64},arg2::Vec{Float64})
    ccall((:MatGetRowSum,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64}),arg1,arg2)
end

function MatTranspose(arg1::Mat{Float64},arg2::MatReuse,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatTranspose,petsc1),PetscErrorCode,(Mat{Float64},MatReuse,Ptr{Mat{Float64}}),arg1,arg2,arg3)
end

function MatHermitianTranspose(arg1::Mat{Float64},arg2::MatReuse,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatHermitianTranspose,petsc1),PetscErrorCode,(Mat{Float64},MatReuse,Ptr{Mat{Float64}}),arg1,arg2,arg3)
end

function MatPermute(arg1::Mat{Float64},arg2::IS{Float64},arg3::IS{Float64},arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatPermute,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},IS{Float64},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4)
end

function MatDiagonalScale(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:MatDiagonalScale,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function MatDiagonalSet(arg1::Mat{Float64},arg2::Vec{Float64},arg3::InsertMode)
    ccall((:MatDiagonalSet,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},InsertMode),arg1,arg2,arg3)
end

function MatEqual(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:MatEqual,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Ptr{PetscBool}),arg1,arg2,arg3)
end

function MatMultEqual(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Integer,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:MatMultEqual,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Int32,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function MatMultAddEqual(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Integer,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:MatMultAddEqual,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Int32,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function MatMultTransposeEqual(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Integer,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:MatMultTransposeEqual,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Int32,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function MatMultTransposeAddEqual(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Integer,arg4::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:MatMultTransposeAddEqual,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Int32,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function MatNorm(arg1::Mat{Float64},arg2::NormType,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:MatNorm,petsc1),PetscErrorCode,(Mat{Float64},NormType,Ptr{Cint}),arg1,arg2,arg3)
end

function MatGetColumnNorms(arg1::Mat{Float64},arg2::NormType,arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:MatGetColumnNorms,petsc1),PetscErrorCode,(Mat{Float64},NormType,Ptr{Cint}),arg1,arg2,arg3)
end

function MatZeroEntries(arg1::Mat{Float64})
    ccall((:MatZeroEntries,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
end

function MatZeroRows(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Float64,arg5::Vec{Float64},arg6::Vec{Float64})
    ccall((:MatZeroRows,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatZeroRowsIS(arg1::Mat{Float64},arg2::IS{Float64},arg3::Float64,arg4::Vec{Float64},arg5::Vec{Float64})
    ccall((:MatZeroRowsIS,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5)
end

#= skipping function with undefined symbols: 
 function MatZeroRowsStencil(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{MatStencil},StridedArray{MatStencil},Ptr{Void}),arg4::Float64,arg5::Vec{Float64},arg6::Vec{Float64})
    ccall((:MatZeroRowsStencil,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{MatStencil},Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
#= skipping function with undefined symbols: 
 function MatZeroRowsColumnsStencil(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{MatStencil},StridedArray{MatStencil},Ptr{Void}),arg4::Float64,arg5::Vec{Float64},arg6::Vec{Float64})
    ccall((:MatZeroRowsColumnsStencil,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{MatStencil},Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
end 
=#
function MatZeroRowsColumns(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Float64,arg5::Vec{Float64},arg6::Vec{Float64})
    ccall((:MatZeroRowsColumns,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatZeroRowsColumnsIS(arg1::Mat{Float64},arg2::IS{Float64},arg3::Float64,arg4::Vec{Float64},arg5::Vec{Float64})
    ccall((:MatZeroRowsColumnsIS,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5)
end

function MatGetSize(arg1::Mat{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatGetSize,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3)
end

function MatGetLocalSize(arg1::Mat{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatGetLocalSize,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3)
end

function MatGetOwnershipRange(arg1::Mat{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatGetOwnershipRange,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3)
end

function MatGetOwnershipRanges(arg1::Mat{Float64},arg2::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:MatGetOwnershipRanges,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Ptr{Int32}}),arg1,arg2)
end

function MatGetOwnershipRangeColumn(arg1::Mat{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatGetOwnershipRangeColumn,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3)
end

function MatGetOwnershipRangesColumn(arg1::Mat{Float64},arg2::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:MatGetOwnershipRangesColumn,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Ptr{Int32}}),arg1,arg2)
end

function MatGetOwnershipIS(arg1::Mat{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:MatGetOwnershipIS,petsc1),PetscErrorCode,(Mat{Float64},Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3)
end

function MatGetSubMatrices(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg5::MatReuse,arg6::Union(Ptr{Ptr{Mat{Float64}}},StridedArray{Ptr{Mat{Float64}}},Ptr{Void}))
    ccall((:MatGetSubMatrices,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{IS{Float64}},Ptr{IS{Float64}},MatReuse,Ptr{Ptr{Mat{Float64}}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatGetSubMatricesMPI(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg5::MatReuse,arg6::Union(Ptr{Ptr{Mat{Float64}}},StridedArray{Ptr{Mat{Float64}}},Ptr{Void}))
    ccall((:MatGetSubMatricesMPI,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{IS{Float64}},Ptr{IS{Float64}},MatReuse,Ptr{Ptr{Mat{Float64}}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatDestroyMatrices(arg1::Integer,arg2::Union(Ptr{Ptr{Mat{Float64}}},StridedArray{Ptr{Mat{Float64}}},Ptr{Void}))
    ccall((:MatDestroyMatrices,petsc1),PetscErrorCode,(Int32,Ptr{Ptr{Mat{Float64}}}),arg1,arg2)
end

function MatGetSubMatrix(arg1::Mat{Float64},arg2::IS{Float64},arg3::IS{Float64},arg4::MatReuse,arg5::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatGetSubMatrix,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},IS{Float64},MatReuse,Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4,arg5)
end

function MatGetLocalSubMatrix(arg1::Mat{Float64},arg2::IS{Float64},arg3::IS{Float64},arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatGetLocalSubMatrix,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},IS{Float64},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4)
end

function MatRestoreLocalSubMatrix(arg1::Mat{Float64},arg2::IS{Float64},arg3::IS{Float64},arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatRestoreLocalSubMatrix,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},IS{Float64},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4)
end

function MatGetSeqNonzeroStructure(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatGetSeqNonzeroStructure,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2)
end

function MatDestroySeqNonzeroStructure(arg1::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatDestroySeqNonzeroStructure,petsc1),PetscErrorCode,(Ptr{Mat{Float64}},),arg1)
end

function MatCreateMPIAIJSumSeqAIJ(arg1::MPI_Comm,arg2::Mat{Float64},arg3::Integer,arg4::Integer,arg5::MatReuse,arg6::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateMPIAIJSumSeqAIJ,petsc1),PetscErrorCode,(comm_type,Mat{Float64},Int32,Int32,MatReuse,Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6)
end

function MatCreateMPIAIJSumSeqAIJSymbolic(arg1::MPI_Comm,arg2::Mat{Float64},arg3::Integer,arg4::Integer,arg5::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateMPIAIJSumSeqAIJSymbolic,petsc1),PetscErrorCode,(comm_type,Mat{Float64},Int32,Int32,Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5)
end

function MatCreateMPIAIJSumSeqAIJNumeric(arg1::Mat{Float64},arg2::Mat{Float64})
    ccall((:MatCreateMPIAIJSumSeqAIJNumeric,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64}),arg1,arg2)
end

function MatMPIAIJGetLocalMat(arg1::Mat{Float64},arg2::MatReuse,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatMPIAIJGetLocalMat,petsc1),PetscErrorCode,(Mat{Float64},MatReuse,Ptr{Mat{Float64}}),arg1,arg2,arg3)
end

function MatMPIAIJGetLocalMatCondensed(arg1::Mat{Float64},arg2::MatReuse,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg5::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatMPIAIJGetLocalMatCondensed,petsc1),PetscErrorCode,(Mat{Float64},MatReuse,Ptr{IS{Float64}},Ptr{IS{Float64}},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4,arg5)
end

function MatGetBrowsOfAcols(arg1::Mat{Float64},arg2::Mat{Float64},arg3::MatReuse,arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg5::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg6::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatGetBrowsOfAcols,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},MatReuse,Ptr{IS{Float64}},Ptr{IS{Float64}},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatGetGhosts(arg1::Mat{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:MatGetGhosts,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Int32},Ptr{Ptr{Int32}}),arg1,arg2,arg3)
end

function MatIncreaseOverlap(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Integer)
    ccall((:MatIncreaseOverlap,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{IS{Float64}},Int32),arg1,arg2,arg3,arg4)
end

function MatMatMult(arg1::Mat{Float64},arg2::Mat{Float64},arg3::MatReuse,Float64::Integer,arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatMatMult,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},MatReuse,Cint,Ptr{Mat{Float64}}),arg1,arg2,arg3,PetscReal,arg4)
end

function MatMatMultSymbolic(arg1::Mat{Float64},arg2::Mat{Float64},Float64::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatMatMultSymbolic,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Cint,Ptr{Mat{Float64}}),arg1,arg2,PetscReal,arg3)
end

function MatMatMultNumeric(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64})
    ccall((:MatMatMultNumeric,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3)
end

function MatMatMatMult(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64},arg4::MatReuse,Float64::Integer,arg5::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatMatMatMult,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64},MatReuse,Cint,Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4,PetscReal,arg5)
end

function MatMatMatMultSymbolic(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64},Float64::Integer,arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatMatMatMultSymbolic,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64},Cint,Ptr{Mat{Float64}}),arg1,arg2,arg3,PetscReal,arg4)
end

function MatMatMatMultNumeric(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64},arg4::Mat{Float64})
    ccall((:MatMatMatMultNumeric,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3,arg4)
end

function MatPtAP(arg1::Mat{Float64},arg2::Mat{Float64},arg3::MatReuse,Float64::Integer,arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatPtAP,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},MatReuse,Cint,Ptr{Mat{Float64}}),arg1,arg2,arg3,PetscReal,arg4)
end

function MatPtAPSymbolic(arg1::Mat{Float64},arg2::Mat{Float64},Float64::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatPtAPSymbolic,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Cint,Ptr{Mat{Float64}}),arg1,arg2,PetscReal,arg3)
end

function MatPtAPNumeric(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64})
    ccall((:MatPtAPNumeric,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3)
end

function MatRARt(arg1::Mat{Float64},arg2::Mat{Float64},arg3::MatReuse,Float64::Integer,arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatRARt,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},MatReuse,Cint,Ptr{Mat{Float64}}),arg1,arg2,arg3,PetscReal,arg4)
end

function MatRARtSymbolic(arg1::Mat{Float64},arg2::Mat{Float64},Float64::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatRARtSymbolic,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Cint,Ptr{Mat{Float64}}),arg1,arg2,PetscReal,arg3)
end

function MatRARtNumeric(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64})
    ccall((:MatRARtNumeric,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3)
end

function MatTransposeMatMult(arg1::Mat{Float64},arg2::Mat{Float64},arg3::MatReuse,Float64::Integer,arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatTransposeMatMult,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},MatReuse,Cint,Ptr{Mat{Float64}}),arg1,arg2,arg3,PetscReal,arg4)
end

function MatTransposetMatMultSymbolic(arg1::Mat{Float64},arg2::Mat{Float64},Float64::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatTransposetMatMultSymbolic,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Cint,Ptr{Mat{Float64}}),arg1,arg2,PetscReal,arg3)
end

function MatTransposetMatMultNumeric(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64})
    ccall((:MatTransposetMatMultNumeric,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3)
end

function MatMatTransposeMult(arg1::Mat{Float64},arg2::Mat{Float64},arg3::MatReuse,Float64::Integer,arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatMatTransposeMult,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},MatReuse,Cint,Ptr{Mat{Float64}}),arg1,arg2,arg3,PetscReal,arg4)
end

function MatMatTransposeMultSymbolic(arg1::Mat{Float64},arg2::Mat{Float64},Float64::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatMatTransposeMultSymbolic,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Cint,Ptr{Mat{Float64}}),arg1,arg2,PetscReal,arg3)
end

function MatMatTransposeMultNumeric(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64})
    ccall((:MatMatTransposeMultNumeric,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3)
end

function MatAXPY(arg1::Mat{Float64},arg2::Float64,arg3::Mat{Float64},arg4::MatStructure)
    ccall((:MatAXPY,petsc1),PetscErrorCode,(Mat{Float64},Float64,Mat{Float64},MatStructure),arg1,arg2,arg3,arg4)
end

function MatAYPX(arg1::Mat{Float64},arg2::Float64,arg3::Mat{Float64},arg4::MatStructure)
    ccall((:MatAYPX,petsc1),PetscErrorCode,(Mat{Float64},Float64,Mat{Float64},MatStructure),arg1,arg2,arg3,arg4)
end

function MatScale(arg1::Mat{Float64},arg2::Float64)
    ccall((:MatScale,petsc1),PetscErrorCode,(Mat{Float64},Float64),arg1,arg2)
end

function MatShift(arg1::Mat{Float64},arg2::Float64)
    ccall((:MatShift,petsc1),PetscErrorCode,(Mat{Float64},Float64),arg1,arg2)
end

function MatSetLocalToGlobalMapping(arg1::Mat{Float64},arg2::ISLocalToGlobalMapping{Float64},arg3::ISLocalToGlobalMapping{Float64})
    ccall((:MatSetLocalToGlobalMapping,petsc1),PetscErrorCode,(Mat{Float64},ISLocalToGlobalMapping{Float64},ISLocalToGlobalMapping{Float64}),arg1,arg2,arg3)
end

function MatGetLocalToGlobalMapping(arg1::Mat{Float64},arg2::Union(Ptr{ISLocalToGlobalMapping{Float64}},StridedArray{ISLocalToGlobalMapping{Float64}},Ptr{Void}),arg3::Union(Ptr{ISLocalToGlobalMapping{Float64}},StridedArray{ISLocalToGlobalMapping{Float64}},Ptr{Void}))
    ccall((:MatGetLocalToGlobalMapping,petsc1),PetscErrorCode,(Mat{Float64},Ptr{ISLocalToGlobalMapping{Float64}},Ptr{ISLocalToGlobalMapping{Float64}}),arg1,arg2,arg3)
end

function MatGetLayouts(arg1::Mat{Float64},arg2::Union(Ptr{PetscLayout{Float64}},StridedArray{PetscLayout{Float64}},Ptr{Void}),arg3::Union(Ptr{PetscLayout{Float64}},StridedArray{PetscLayout{Float64}},Ptr{Void}))
    ccall((:MatGetLayouts,petsc1),PetscErrorCode,(Mat{Float64},Ptr{PetscLayout{Float64}{Float64}},Ptr{PetscLayout{Float64}{Float64}}),arg1,arg2,arg3)
end

function MatZeroRowsLocal(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Float64,arg5::Vec{Float64},arg6::Vec{Float64})
    ccall((:MatZeroRowsLocal,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatZeroRowsLocalIS(arg1::Mat{Float64},arg2::IS{Float64},arg3::Float64,arg4::Vec{Float64},arg5::Vec{Float64})
    ccall((:MatZeroRowsLocalIS,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5)
end

function MatZeroRowsColumnsLocal(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Float64,arg5::Vec{Float64},arg6::Vec{Float64})
    ccall((:MatZeroRowsColumnsLocal,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatZeroRowsColumnsLocalIS(arg1::Mat{Float64},arg2::IS{Float64},arg3::Float64,arg4::Vec{Float64},arg5::Vec{Float64})
    ccall((:MatZeroRowsColumnsLocalIS,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},Float64,Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5)
end

function MatSetValuesLocal(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::InsertMode)
    ccall((:MatSetValuesLocal,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Int32,Ptr{Int32},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatSetValuesBlockedLocal(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::InsertMode)
    ccall((:MatSetValuesBlockedLocal,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Int32,Ptr{Int32},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatStashSetInitialSize(arg1::Mat{Float64},arg2::Integer,arg3::Integer)
    ccall((:MatStashSetInitialSize,petsc1),PetscErrorCode,(Mat{Float64},Int32,Int32),arg1,arg2,arg3)
end

function MatStashGetInfo(arg1::Mat{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatStashGetInfo,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3,arg4,arg5)
end

function MatInterpolate(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:MatInterpolate,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function MatInterpolateAdd(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64})
    ccall((:MatInterpolateAdd,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
end

function MatRestrict(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:MatRestrict,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function MatCreateVecs(arg1::Mat{Float64},arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:MatCreateVecs,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Vec{Float64}},Ptr{Vec{Float64}}),arg1,arg2,arg3)
end

function MatGetMultiProcBlock(arg1::Mat{Float64},arg2::MPI_Comm,arg3::MatReuse,arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatGetMultiProcBlock,petsc1),PetscErrorCode,(Mat{Float64},comm_type,MatReuse,Ptr{Mat{Float64}}),arg1,arg2.val,arg3,arg4)
end

function MatFindZeroDiagonals(arg1::Mat{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:MatFindZeroDiagonals,petsc1),PetscErrorCode,(Mat{Float64},Ptr{IS{Float64}}),arg1,arg2)
end

function MatFindOffBlockDiagonalEntries(arg1::Mat{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:MatFindOffBlockDiagonalEntries,petsc1),PetscErrorCode,(Mat{Float64},Ptr{IS{Float64}}),arg1,arg2)
end

function MatCreateMPIMatConcatenateSeqMat(arg1::MPI_Comm,arg2::Mat{Float64},arg3::Integer,arg4::MatReuse,arg5::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateMPIMatConcatenateSeqMat,petsc1),PetscErrorCode,(comm_type,Mat{Float64},Int32,MatReuse,Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5)
end

function MatInodeAdjustForInodes(arg1::Mat{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:MatInodeAdjustForInodes,petsc1),PetscErrorCode,(Mat{Float64},Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3)
end

function MatInodeGetInodeSizes(arg1::Mat{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatInodeGetInodeSizes,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Int32},Ptr{Ptr{Int32}},Ptr{Int32}),arg1,arg2,arg3,arg4)
end

function MatSeqAIJSetColumnIndices(arg1::Mat{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatSeqAIJSetColumnIndices,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Int32}),arg1,arg2)
end

function MatSeqBAIJSetColumnIndices(arg1::Mat{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatSeqBAIJSetColumnIndices,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Int32}),arg1,arg2)
end

function MatCreateSeqAIJWithArrays(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateSeqAIJWithArrays,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateSeqBAIJWithArrays(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg8::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateSeqBAIJWithArrays,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatCreateSeqSBAIJWithArrays(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg8::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateSeqSBAIJWithArrays,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatCreateSeqAIJFromTriple(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg8::Integer,arg9::PetscBool)
    ccall((:MatCreateSeqAIJFromTriple,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Mat{Float64}},Int32,PetscBool),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function MatSeqBAIJSetPreallocation(arg1::Mat{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatSeqBAIJSetPreallocation,petsc1),PetscErrorCode,(Mat{Float64},Int32,Int32,Ptr{Int32}),arg1,arg2,arg3,arg4)
end

function MatSeqSBAIJSetPreallocation(arg1::Mat{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatSeqSBAIJSetPreallocation,petsc1),PetscErrorCode,(Mat{Float64},Int32,Int32,Ptr{Int32}),arg1,arg2,arg3,arg4)
end

function MatSeqAIJSetPreallocation(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatSeqAIJSetPreallocation,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32}),arg1,arg2,arg3)
end

function MatMPIBAIJSetPreallocation(arg1::Mat{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Integer,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatMPIBAIJSetPreallocation,petsc1),PetscErrorCode,(Mat{Float64},Int32,Int32,Ptr{Int32},Int32,Ptr{Int32}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatMPISBAIJSetPreallocation(arg1::Mat{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Integer,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatMPISBAIJSetPreallocation,petsc1),PetscErrorCode,(Mat{Float64},Int32,Int32,Ptr{Int32},Int32,Ptr{Int32}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatMPIAIJSetPreallocation(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatMPIAIJSetPreallocation,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Int32,Ptr{Int32}),arg1,arg2,arg3,arg4,arg5)
end

function MatSeqAIJSetPreallocationCSR(arg1::Mat{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:MatSeqAIJSetPreallocationCSR,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Int32},Ptr{Int32},Ptr{Float64}),arg1,arg2,arg3,arg4)
end

function MatSeqBAIJSetPreallocationCSR(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:MatSeqBAIJSetPreallocationCSR,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5)
end

function MatMPIAIJSetPreallocationCSR(arg1::Mat{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:MatMPIAIJSetPreallocationCSR,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Int32},Ptr{Int32},Ptr{Float64}),arg1,arg2,arg3,arg4)
end

function MatMPIBAIJSetPreallocationCSR(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:MatMPIBAIJSetPreallocationCSR,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5)
end

function MatMPIAdjSetPreallocation(arg1::Mat{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatMPIAdjSetPreallocation,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Int32},Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3,arg4)
end

function MatMPIDenseSetPreallocation(arg1::Mat{Float64},arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:MatMPIDenseSetPreallocation,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Float64}),arg1,arg2)
end

function MatSeqDenseSetPreallocation(arg1::Mat{Float64},arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:MatSeqDenseSetPreallocation,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Float64}),arg1,arg2)
end

function MatMPIAIJGetSeqAIJ(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg4::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:MatMPIAIJGetSeqAIJ,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}},Ptr{Mat{Float64}},Ptr{Ptr{Int32}}),arg1,arg2,arg3,arg4)
end

function MatMPIBAIJGetSeqBAIJ(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg4::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:MatMPIBAIJGetSeqBAIJ,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}},Ptr{Mat{Float64}},Ptr{Ptr{Int32}}),arg1,arg2,arg3,arg4)
end

function MatMPIAdjCreateNonemptySubcommMat(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatMPIAdjCreateNonemptySubcommMat,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2)
end

function MatISSetPreallocation(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatISSetPreallocation,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Int32,Ptr{Int32}),arg1,arg2,arg3,arg4,arg5)
end

function MatSeqDenseSetLDA(arg1::Mat{Float64},arg2::Integer)
    ccall((:MatSeqDenseSetLDA,petsc1),PetscErrorCode,(Mat{Float64},Int32),arg1,arg2)
end

function MatDenseGetLocalMatrix(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatDenseGetLocalMatrix,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2)
end

function MatStoreValues(arg1::Mat{Float64})
    ccall((:MatStoreValues,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
end

function MatRetrieveValues(arg1::Mat{Float64})
    ccall((:MatRetrieveValues,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
end

function MatDAADSetCtx(arg1::Mat{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:MatDAADSetCtx,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Void}),arg1,arg2)
end

function MatFindNonzeroRows(arg1::Mat{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:MatFindNonzeroRows,petsc1),PetscErrorCode,(Mat{Float64},Ptr{IS{Float64}}),arg1,arg2)
end

function MatGetOrdering(arg1::Mat{Float64},arg2::MatOrderingType,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:MatGetOrdering,petsc1),PetscErrorCode,(Mat{Float64},Cstring,Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4)
end

#= skipping function with undefined symbols: 
 function MatGetOrderingList(arg0::Type{Float64},arg1::Union(Ptr{PetscFunctionList},StridedArray{PetscFunctionList},Ptr{Void}))
    ccall((:MatGetOrderingList,petsc1),PetscErrorCode,(Ptr{PetscFunctionList},),arg1)
end 
=#
function MatOrderingRegister(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:MatOrderingRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
end

function MatReorderForNonzeroDiagonal(arg1::Mat{Float64},Float64::Integer,arg2::IS{Float64},arg3::IS{Float64})
    ccall((:MatReorderForNonzeroDiagonal,petsc1),PetscErrorCode,(Mat{Float64},Cint,IS{Float64},IS{Float64}),arg1,PetscReal,arg2,arg3)
end

function MatCreateLaplacian(arg1::Mat{Float64},Float64::Integer,arg2::PetscBool,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateLaplacian,petsc1),PetscErrorCode,(Mat{Float64},Cint,PetscBool,Ptr{Mat{Float64}}),arg1,PetscReal,arg2,arg3)
end

function MatFactorInfoInitialize(arg0::Type{Float64},arg1::Union(Ptr{MatFactorInfo},StridedArray{MatFactorInfo},Ptr{Void}))
    ccall((:MatFactorInfoInitialize,petsc1),PetscErrorCode,(Ptr{MatFactorInfo},),arg1)
end

function MatCholeskyFactor(arg1::Mat{Float64},arg2::IS{Float64},arg3::Union(Ptr{MatFactorInfo},StridedArray{MatFactorInfo},Ptr{Void}))
    ccall((:MatCholeskyFactor,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},Ptr{MatFactorInfo}),arg1,arg2,arg3)
end

function MatCholeskyFactorSymbolic(arg1::Mat{Float64},arg2::Mat{Float64},arg3::IS{Float64},arg4::Union(Ptr{MatFactorInfo},StridedArray{MatFactorInfo},Ptr{Void}))
    ccall((:MatCholeskyFactorSymbolic,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},IS{Float64},Ptr{MatFactorInfo}),arg1,arg2,arg3,arg4)
end

function MatCholeskyFactorNumeric(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Union(Ptr{MatFactorInfo},StridedArray{MatFactorInfo},Ptr{Void}))
    ccall((:MatCholeskyFactorNumeric,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Ptr{MatFactorInfo}),arg1,arg2,arg3)
end

function MatLUFactor(arg1::Mat{Float64},arg2::IS{Float64},arg3::IS{Float64},arg4::Union(Ptr{MatFactorInfo},StridedArray{MatFactorInfo},Ptr{Void}))
    ccall((:MatLUFactor,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},IS{Float64},Ptr{MatFactorInfo}),arg1,arg2,arg3,arg4)
end

function MatILUFactor(arg1::Mat{Float64},arg2::IS{Float64},arg3::IS{Float64},arg4::Union(Ptr{MatFactorInfo},StridedArray{MatFactorInfo},Ptr{Void}))
    ccall((:MatILUFactor,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},IS{Float64},Ptr{MatFactorInfo}),arg1,arg2,arg3,arg4)
end

function MatLUFactorSymbolic(arg1::Mat{Float64},arg2::Mat{Float64},arg3::IS{Float64},arg4::IS{Float64},arg5::Union(Ptr{MatFactorInfo},StridedArray{MatFactorInfo},Ptr{Void}))
    ccall((:MatLUFactorSymbolic,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},IS{Float64},IS{Float64},Ptr{MatFactorInfo}),arg1,arg2,arg3,arg4,arg5)
end

function MatILUFactorSymbolic(arg1::Mat{Float64},arg2::Mat{Float64},arg3::IS{Float64},arg4::IS{Float64},arg5::Union(Ptr{MatFactorInfo},StridedArray{MatFactorInfo},Ptr{Void}))
    ccall((:MatILUFactorSymbolic,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},IS{Float64},IS{Float64},Ptr{MatFactorInfo}),arg1,arg2,arg3,arg4,arg5)
end

function MatICCFactorSymbolic(arg1::Mat{Float64},arg2::Mat{Float64},arg3::IS{Float64},arg4::Union(Ptr{MatFactorInfo},StridedArray{MatFactorInfo},Ptr{Void}))
    ccall((:MatICCFactorSymbolic,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},IS{Float64},Ptr{MatFactorInfo}),arg1,arg2,arg3,arg4)
end

function MatICCFactor(arg1::Mat{Float64},arg2::IS{Float64},arg3::Union(Ptr{MatFactorInfo},StridedArray{MatFactorInfo},Ptr{Void}))
    ccall((:MatICCFactor,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},Ptr{MatFactorInfo}),arg1,arg2,arg3)
end

function MatLUFactorNumeric(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Union(Ptr{MatFactorInfo},StridedArray{MatFactorInfo},Ptr{Void}))
    ccall((:MatLUFactorNumeric,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Ptr{MatFactorInfo}),arg1,arg2,arg3)
end

function MatGetInertia(arg1::Mat{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatGetInertia,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Int32},Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3,arg4)
end

function MatSolve(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:MatSolve,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function MatForwardSolve(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:MatForwardSolve,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function MatBackwardSolve(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:MatBackwardSolve,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function MatSolveAdd(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64})
    ccall((:MatSolveAdd,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
end

function MatSolveTranspose(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:MatSolveTranspose,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function MatSolveTransposeAdd(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64})
    ccall((:MatSolveTransposeAdd,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
end

#= skipping function with undefined symbols: 
 function MatSolves(arg1::Mat{Float64},arg2::Vecs,arg3::Vecs)
    ccall((:MatSolves,petsc1),PetscErrorCode,(Mat{Float64},Vecs,Vecs),arg1,arg2,arg3)
end 
=#
function MatSetUnfactored(arg1::Mat{Float64})
    ccall((:MatSetUnfactored,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
end

function MatSOR(arg1::Mat{Float64},arg2::Vec{Float64},Float64::Integer,arg3::MatSORType,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Vec{Float64})
    ccall((:MatSOR,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Cint,MatSORType,Cint,Int32,Int32,Vec{Float64}),arg1,arg2,PetscReal,arg3,arg4,arg5,arg6,arg7)
end

#= skipping function with undefined symbols: 
 function MatColoringCreate(arg1::Mat{Float64},arg2::Union(Ptr{MatColoring},StridedArray{MatColoring},Ptr{Void}))
    ccall((:MatColoringCreate,petsc1),PetscErrorCode,(Mat{Float64},Ptr{MatColoring}),arg1,arg2)
end 
=#
function MatColoringGetDegrees(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatColoringGetDegrees,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32}),arg1,arg2,arg3)
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
    ccall((:MatColoringSetDistance,petsc1),PetscErrorCode,(MatColoring,Int32),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatColoringGetDistance(arg0::Type{Float64},arg1::MatColoring,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatColoringGetDistance,petsc1),PetscErrorCode,(MatColoring,Ptr{Int32}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatColoringSetMaxColors(arg0::Type{Float64},arg1::MatColoring,arg2::Integer)
    ccall((:MatColoringSetMaxColors,petsc1),PetscErrorCode,(MatColoring,Int32),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatColoringGetMaxColors(arg0::Type{Float64},arg1::MatColoring,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatColoringGetMaxColors,petsc1),PetscErrorCode,(MatColoring,Ptr{Int32}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatColoringApply(arg1::MatColoring,arg2::Union(Ptr{ISColoring{Float64}},StridedArray{ISColoring{Float64}},Ptr{Void}))
    ccall((:MatColoringApply,petsc1),PetscErrorCode,(MatColoring,Ptr{ISColoring{Float64}}),arg1,arg2)
end 
=#
function MatColoringRegister(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:MatColoringRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
end

function MatColoringPatch(arg1::Mat{Float64},arg2::Integer,arg3::Integer,ISColoringValue::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{ISColoring{Float64}},StridedArray{ISColoring{Float64}},Ptr{Void}))
    ccall((:MatColoringPatch,petsc1),PetscErrorCode,(Mat{Float64},Int32,Int32,Ptr{Cint},Ptr{ISColoring{Float64}}),arg1,arg2,arg3,ISColoringValue,arg4)
end

#= skipping function with undefined symbols: 
 function MatColoringSetWeightType(arg0::Type{Float64},arg1::MatColoring,arg2::MatColoringWeightType)
    ccall((:MatColoringSetWeightType,petsc1),PetscErrorCode,(MatColoring,MatColoringWeightType),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatColoringSetWeights(arg0::Type{Float64},arg1::MatColoring,arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatColoringSetWeights,petsc1),PetscErrorCode,(MatColoring,Ptr{Cint},Ptr{Int32}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function MatColoringCreateWeights(arg0::Type{Float64},arg1::MatColoring,arg2::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),lperm::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:MatColoringCreateWeights,petsc1),PetscErrorCode,(MatColoring,Ptr{Ptr{Cint}},Ptr{Ptr{Int32}}),arg1,arg2,lperm)
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
 function MatFDColoringSetParameters(arg0::Type{Float64},arg1::MatFDColoring,Float64::Integer,arg2::Integer)
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
 function MatFDColoringGetPerturbedColumns(arg0::Type{Float64},arg1::MatFDColoring,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:MatFDColoringGetPerturbedColumns,petsc1),PetscErrorCode,(MatFDColoring,Ptr{Int32},Ptr{Ptr{Int32}}),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function MatFDColoringSetUp(arg1::Mat{Float64},arg2::ISColoring{Float64},arg3::MatFDColoring)
    ccall((:MatFDColoringSetUp,petsc1),PetscErrorCode,(Mat{Float64},ISColoring{Float64},MatFDColoring),arg1,arg2,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function MatFDColoringSetBlockSize(arg0::Type{Float64},arg1::MatFDColoring,arg2::Integer,arg3::Integer)
    ccall((:MatFDColoringSetBlockSize,petsc1),PetscErrorCode,(MatFDColoring,Int32,Int32),arg1,arg2,arg3)
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
    ccall((:MatPartitioningSetNParts,petsc1),PetscErrorCode,(MatPartitioning,Int32),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningSetAdjacency(arg1::MatPartitioning,arg2::Mat{Float64})
    ccall((:MatPartitioningSetAdjacency,petsc1),PetscErrorCode,(MatPartitioning,Mat{Float64}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningSetVertexWeights(arg0::Type{Float64},arg1::MatPartitioning,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatPartitioningSetVertexWeights,petsc1),PetscErrorCode,(MatPartitioning,Ptr{Int32}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningSetPartitionWeights(arg0::Type{Float64},arg1::MatPartitioning,Float64::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
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
function MatPartitioningRegister(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:MatPartitioningRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
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
    ccall((:MatPartitioningGetType,petsc1),PetscErrorCode,(MatPartitioning,Ptr{Cstring}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningParmetisSetCoarseSequential(arg0::Type{Float64},arg1::MatPartitioning)
    ccall((:MatPartitioningParmetisSetCoarseSequential,petsc1),PetscErrorCode,(MatPartitioning,),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningParmetisGetEdgeCut(arg0::Type{Float64},arg1::MatPartitioning,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatPartitioningParmetisGetEdgeCut,petsc1),PetscErrorCode,(MatPartitioning,Ptr{Int32}),arg1,arg2)
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
 function MatPartitioningChacoSetCoarseLevel(arg0::Type{Float64},arg1::MatPartitioning,Float64::Integer)
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
 function MatPartitioningChacoSetEigenTol(arg0::Type{Float64},arg1::MatPartitioning,Float64::Integer)
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
    ccall((:MatPartitioningChacoSetEigenNumber,petsc1),PetscErrorCode,(MatPartitioning,Int32),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningChacoGetEigenNumber(arg0::Type{Float64},arg1::MatPartitioning,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatPartitioningChacoGetEigenNumber,petsc1),PetscErrorCode,(MatPartitioning,Ptr{Int32}),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningPartySetGlobal(arg1::MatPartitioning,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:MatPartitioningPartySetGlobal,petsc1),PetscErrorCode,(MatPartitioning,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningPartySetLocal(arg1::MatPartitioning,arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:MatPartitioningPartySetLocal,petsc1),PetscErrorCode,(MatPartitioning,Cstring),arg1,arg2)
end 
=#
#= skipping function with undefined symbols: 
 function MatPartitioningPartySetCoarseLevel(arg0::Type{Float64},arg1::MatPartitioning,Float64::Integer)
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
 function MatPartitioningPTScotchSetImbalance(arg0::Type{Float64},arg1::MatPartitioning,Float64::Integer)
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
function MatCoarsenRegister(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:MatCoarsenRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
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
    ccall((:MatCoarsenGetType,petsc1),PetscErrorCode,(MatCoarsen,Ptr{Cstring}),arg1,arg2)
end 
=#
function MatMeshToCellGraph(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatMeshToCellGraph,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Mat{Float64}}),arg1,arg2,arg3)
end

function MatHasOperation(arg1::Mat{Float64},arg2::MatOperation,arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:MatHasOperation,petsc1),PetscErrorCode,(Mat{Float64},MatOperation,Ptr{PetscBool}),arg1,arg2,arg3)
end

function MatShellSetOperation(arg1::Mat{Float64},arg2::MatOperation,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:MatShellSetOperation,petsc1),PetscErrorCode,(Mat{Float64},MatOperation,Ptr{Void}),arg1,arg2,arg3)
end

function MatShellGetOperation(arg1::Mat{Float64},arg2::MatOperation,arg3::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:MatShellGetOperation,petsc1),PetscErrorCode,(Mat{Float64},MatOperation,Ptr{Ptr{Void}}),arg1,arg2,arg3)
end

function MatShellSetContext(arg1::Mat{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:MatShellSetContext,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Void}),arg1,arg2)
end

function MatMPIBAIJSetHashTableFactor(arg1::Mat{Float64},Float64::Integer)
    ccall((:MatMPIBAIJSetHashTableFactor,petsc1),PetscErrorCode,(Mat{Float64},Cint),arg1,PetscReal)
end

function MatISGetLocalMat(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatISGetLocalMat,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2)
end

function MatISSetLocalMat(arg1::Mat{Float64},arg2::Mat{Float64})
    ccall((:MatISSetLocalMat,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64}),arg1,arg2)
end

function MatISGetMPIXAIJ(arg1::Mat{Float64},arg2::MatReuse,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatISGetMPIXAIJ,petsc1),PetscErrorCode,(Mat{Float64},MatReuse,Ptr{Mat{Float64}}),arg1,arg2,arg3)
end

#= skipping function with undefined symbols: 
 function MatNullSpaceCreate(arg1::MPI_Comm,arg2::PetscBool,arg3::Integer,arg4::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}),arg5::Union(Ptr{MatNullSpace},StridedArray{MatNullSpace},Ptr{Void}))
    ccall((:MatNullSpaceCreate,petsc1),PetscErrorCode,(comm_type,PetscBool,Int32,Ptr{Vec{Float64}},Ptr{MatNullSpace}),arg1.val,arg2,arg3,arg4,arg5)
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
 function MatNullSpaceGetVecs(arg1::MatNullSpace,arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Ptr{Vec{Float64}}},StridedArray{Ptr{Vec{Float64}}},Ptr{Void}))
    ccall((:MatNullSpaceGetVecs,petsc1),PetscErrorCode,(MatNullSpace,Ptr{PetscBool},Ptr{Int32},Ptr{Ptr{Vec{Float64}}}),arg1,arg2,arg3,arg4)
end 
=#
#= skipping function with undefined symbols: 
 function MatNullSpaceCreateRigidBody(arg1::Vec{Float64},arg2::Union(Ptr{MatNullSpace},StridedArray{MatNullSpace},Ptr{Void}))
    ccall((:MatNullSpaceCreateRigidBody,petsc1),PetscErrorCode,(Vec{Float64},Ptr{MatNullSpace}),arg1,arg2)
end 
=#
function MatReorderingSeqSBAIJ(arg1::Mat{Float64},arg2::IS{Float64})
    ccall((:MatReorderingSeqSBAIJ,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64}),arg1,arg2)
end

function MatMPISBAIJSetHashTableFactor(arg1::Mat{Float64},Float64::Integer)
    ccall((:MatMPISBAIJSetHashTableFactor,petsc1),PetscErrorCode,(Mat{Float64},Cint),arg1,PetscReal)
end

function MatSeqSBAIJSetColumnIndices(arg1::Mat{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatSeqSBAIJSetColumnIndices,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Int32}),arg1,arg2)
end

function MatSeqBAIJInvertBlockDiagonal(arg1::Mat{Float64})
    ccall((:MatSeqBAIJInvertBlockDiagonal,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
end

function MatCreateMAIJ(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateMAIJ,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Mat{Float64}}),arg1,arg2,arg3)
end

function MatMAIJRedimension(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatMAIJRedimension,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Mat{Float64}}),arg1,arg2,arg3)
end

function MatMAIJGetAIJ(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatMAIJGetAIJ,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2)
end

function MatComputeExplicitOperator(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatComputeExplicitOperator,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2)
end

function MatDiagonalScaleLocal(arg1::Mat{Float64},arg2::Vec{Float64})
    ccall((:MatDiagonalScaleLocal,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64}),arg1,arg2)
end

function MatCreateMFFD(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateMFFD,petsc1),PetscErrorCode,(comm_type,Int32,Int32,Int32,Int32,Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6)
end

function MatMFFDSetBase(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:MatMFFDSetBase,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function MatMFFDSetFunction(arg1::Mat{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:MatMFFDSetFunction,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function MatMFFDSetFunctioni(arg1::Mat{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:MatMFFDSetFunctioni,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Void}),arg1,arg2)
end

function MatMFFDSetFunctioniBase(arg1::Mat{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:MatMFFDSetFunctioniBase,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Void}),arg1,arg2)
end

#= skipping function with undefined symbols: 
 function MatMFFDAddNullSpace(arg1::Mat{Float64},arg2::MatNullSpace)
    ccall((:MatMFFDAddNullSpace,petsc1),PetscErrorCode,(Mat{Float64},MatNullSpace),arg1,arg2)
end 
=#
function MatMFFDSetHHistory(arg1::Mat{Float64},arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}),arg3::Integer)
    ccall((:MatMFFDSetHHistory,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Float64},Int32),arg1,arg2,arg3)
end

function MatMFFDResetHHistory(arg1::Mat{Float64})
    ccall((:MatMFFDResetHHistory,petsc1),PetscErrorCode,(Mat{Float64},),arg1)
end

function MatMFFDSetFunctionError(arg1::Mat{Float64},Float64::Integer)
    ccall((:MatMFFDSetFunctionError,petsc1),PetscErrorCode,(Mat{Float64},Cint),arg1,PetscReal)
end

function MatMFFDSetPeriod(arg1::Mat{Float64},arg2::Integer)
    ccall((:MatMFFDSetPeriod,petsc1),PetscErrorCode,(Mat{Float64},Int32),arg1,arg2)
end

function MatMFFDGetH(arg1::Mat{Float64},arg2::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:MatMFFDGetH,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Float64}),arg1,arg2)
end

function MatMFFDSetOptionsPrefix(arg1::Mat{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:MatMFFDSetOptionsPrefix,petsc1),PetscErrorCode,(Mat{Float64},Cstring),arg1,arg2)
end

function MatMFFDCheckPositivity(arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{Float64},StridedArray{Float64},Ptr{Void}))
    ccall((:MatMFFDCheckPositivity,petsc1),PetscErrorCode,(Ptr{Void},Vec{Float64},Vec{Float64},Ptr{Float64}),arg1,arg2,arg3,arg4)
end

function MatMFFDSetCheckh(arg1::Mat{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:MatMFFDSetCheckh,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function MatMFFDSetType(arg1::Mat{Float64},arg2::MatMFFDType)
    ccall((:MatMFFDSetType,petsc1),PetscErrorCode,(Mat{Float64},Cstring),arg1,arg2)
end

function MatMFFDRegister(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:MatMFFDRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
end

function MatMFFDDSSetUmin(arg1::Mat{Float64},Float64::Integer)
    ccall((:MatMFFDDSSetUmin,petsc1),PetscErrorCode,(Mat{Float64},Cint),arg1,PetscReal)
end

function MatMFFDWPSetComputeNormU(arg1::Mat{Float64},arg2::PetscBool)
    ccall((:MatMFFDWPSetComputeNormU,petsc1),PetscErrorCode,(Mat{Float64},PetscBool),arg1,arg2)
end

function PetscViewerMathematicaPutMatrix(arg1::PetscViewer{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscViewerMathematicaPutMatrix,petsc1),PetscErrorCode,(PetscViewer{Float64},Int32,Int32,Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function PetscViewerMathematicaPutCSRMatrix(arg1::PetscViewer{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PetscViewerMathematicaPutCSRMatrix,petsc1),PetscErrorCode,(PetscViewer{Float64},Int32,Int32,Ptr{Int32},Ptr{Int32},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatCreateNest(arg1::MPI_Comm,arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg6::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateNest,petsc1),PetscErrorCode,(comm_type,Int32,Ptr{IS{Float64}},Int32,Ptr{IS{Float64}},Ptr{Mat{Float64}},Ptr{Mat{Float64}}),arg1.val,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatNestGetSize(arg1::Mat{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatNestGetSize,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3)
end

function MatNestGetISs(arg1::Mat{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:MatNestGetISs,petsc1),PetscErrorCode,(Mat{Float64},Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3)
end

function MatNestGetLocalISs(arg1::Mat{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:MatNestGetLocalISs,petsc1),PetscErrorCode,(Mat{Float64},Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3)
end

function MatNestGetSubMats(arg1::Mat{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Ptr{Ptr{Mat{Float64}}}},StridedArray{Ptr{Ptr{Mat{Float64}}}},Ptr{Void}))
    ccall((:MatNestGetSubMats,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Int32},Ptr{Int32},Ptr{Ptr{Ptr{Mat{Float64}}}}),arg1,arg2,arg3,arg4)
end

function MatNestGetSubMat(arg1::Mat{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatNestGetSubMat,petsc1),PetscErrorCode,(Mat{Float64},Int32,Int32,Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4)
end

function MatNestSetVecType(arg1::Mat{Float64},arg2::VecType)
    ccall((:MatNestSetVecType,petsc1),PetscErrorCode,(Mat{Float64},Cstring),arg1,arg2)
end

function MatNestSetSubMats(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg6::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatNestSetSubMats,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{IS{Float64}},Int32,Ptr{IS{Float64}},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatNestSetSubMat(arg1::Mat{Float64},arg2::Integer,arg3::Integer,arg4::Mat{Float64})
    ccall((:MatNestSetSubMat,petsc1),PetscErrorCode,(Mat{Float64},Int32,Int32,Mat{Float64}),arg1,arg2,arg3,arg4)
end

function MatChop(arg1::Mat{Float64},Float64::Integer)
    ccall((:MatChop,petsc1),PetscErrorCode,(Mat{Float64},Cint),arg1,PetscReal)
end

function MatComputeBandwidth(arg1::Mat{Float64},Float64::Integer,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatComputeBandwidth,petsc1),PetscErrorCode,(Mat{Float64},Cint,Ptr{Int32}),arg1,PetscReal,arg2)
end

function MatSubdomainsCreateCoalesce(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    ccall((:MatSubdomainsCreateCoalesce,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Ptr{Ptr{IS{Float64}}}),arg1,arg2,arg3,arg4)
end

function PCExoticSetType(arg1::PC{Float64},arg2::PCExoticType)
    ccall((:PCExoticSetType,petsc1),PetscErrorCode,(PC{Float64},PCExoticType),arg1,arg2)
end

function PCInitializePackage(arg0::Type{Float64})
    ccall((:PCInitializePackage,petsc1),PetscErrorCode,())
end

function PCCreate(arg1::MPI_Comm,arg2::Union(Ptr{PC{Float64}},StridedArray{PC{Float64}},Ptr{Void}))
    ccall((:PCCreate,petsc1),PetscErrorCode,(comm_type,Ptr{PC{Float64}}),arg1.val,arg2)
end

function PCSetType(arg1::PC{Float64},arg2::PCType)
    ccall((:PCSetType,petsc1),PetscErrorCode,(PC{Float64},Cstring),arg1,arg2)
end

function PCGetType(arg1::PC{Float64},arg2::Union(Ptr{PCType},StridedArray{PCType},Ptr{Void}))
    ccall((:PCGetType,petsc1),PetscErrorCode,(PC{Float64},Ptr{Ptr{UInt8}}),arg1,arg2)
end

function PCSetUp(arg1::PC{Float64})
    ccall((:PCSetUp,petsc1),PetscErrorCode,(PC{Float64},),arg1)
end

function PCGetSetUpFailedReason(arg1::PC{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PCGetSetUpFailedReason,petsc1),PetscErrorCode,(PC{Float64},Ptr{Int32}),arg1,arg2)
end

function PCSetUpOnBlocks(arg1::PC{Float64})
    ccall((:PCSetUpOnBlocks,petsc1),PetscErrorCode,(PC{Float64},),arg1)
end

function PCApply(arg1::PC{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:PCApply,petsc1),PetscErrorCode,(PC{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function PCApplySymmetricLeft(arg1::PC{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:PCApplySymmetricLeft,petsc1),PetscErrorCode,(PC{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function PCApplySymmetricRight(arg1::PC{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:PCApplySymmetricRight,petsc1),PetscErrorCode,(PC{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function PCApplyBAorAB(arg1::PC{Float64},arg2::PCSide,arg3::Vec{Float64},arg4::Vec{Float64},arg5::Vec{Float64})
    ccall((:PCApplyBAorAB,petsc1),PetscErrorCode,(PC{Float64},PCSide,Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5)
end

function PCApplyTranspose(arg1::PC{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:PCApplyTranspose,petsc1),PetscErrorCode,(PC{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function PCApplyTransposeExists(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PCApplyTransposeExists,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
end

function PCApplyBAorABTranspose(arg1::PC{Float64},arg2::PCSide,arg3::Vec{Float64},arg4::Vec{Float64},arg5::Vec{Float64})
    ccall((:PCApplyBAorABTranspose,petsc1),PetscErrorCode,(PC{Float64},PCSide,Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5)
end

function PCSetReusePreconditioner(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCSetReusePreconditioner,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCGetReusePreconditioner(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PCGetReusePreconditioner,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
end

function PCSetErrorIfFailure(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCSetErrorIfFailure,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCApplyRichardson(arg1::PC{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64},Float64::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::PetscBool,arg9::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg10::Union(Ptr{PCRichardsonConvergedReason},StridedArray{PCRichardsonConvergedReason},Ptr{Void}))
    ccall((:PCApplyRichardson,petsc1),PetscErrorCode,(PC{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Cint,Cint,Cint,Int32,PetscBool,Ptr{Int32},Ptr{PCRichardsonConvergedReason}),arg1,arg2,arg3,arg4,PetscReal,arg5,arg6,arg7,arg8,arg9,arg10)
end

function PCApplyRichardsonExists(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PCApplyRichardsonExists,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
end

function PCSetInitialGuessNonzero(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCSetInitialGuessNonzero,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCGetInitialGuessNonzero(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PCGetInitialGuessNonzero,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
end

function PCSetUseAmat(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCSetUseAmat,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCGetUseAmat(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PCGetUseAmat,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
end

function PCRegister(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PCRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
end

function PCReset(arg1::PC{Float64})
    ccall((:PCReset,petsc1),PetscErrorCode,(PC{Float64},),arg1)
end

function PCDestroy(arg1::Union(Ptr{PC{Float64}},StridedArray{PC{Float64}},Ptr{Void}))
    ccall((:PCDestroy,petsc1),PetscErrorCode,(Ptr{PC{Float64}},),arg1)
end

function PCSetFromOptions(arg1::PC{Float64})
    ccall((:PCSetFromOptions,petsc1),PetscErrorCode,(PC{Float64},),arg1)
end

function PCFactorGetMatrix(arg1::PC{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:PCFactorGetMatrix,petsc1),PetscErrorCode,(PC{Float64},Ptr{Mat{Float64}}),arg1,arg2)
end

function PCSetModifySubMatrices(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PCSetModifySubMatrices,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function PCModifySubMatrices(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg5::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg6::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PCModifySubMatrices,petsc1),PetscErrorCode,(PC{Float64},Int32,Ptr{IS{Float64}},Ptr{IS{Float64}},Ptr{Mat{Float64}},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PCSetOperators(arg1::PC{Float64},arg2::Mat{Float64},arg3::Mat{Float64})
    ccall((:PCSetOperators,petsc1),PetscErrorCode,(PC{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3)
end

function PCGetOperators(arg1::PC{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:PCGetOperators,petsc1),PetscErrorCode,(PC{Float64},Ptr{Mat{Float64}},Ptr{Mat{Float64}}),arg1,arg2,arg3)
end

function PCGetOperatorsSet(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PCGetOperatorsSet,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3)
end

function PCView(arg1::PC{Float64},arg2::PetscViewer{Float64})
    ccall((:PCView,petsc1),PetscErrorCode,(PC{Float64},PetscViewer{Float64}),arg1,arg2)
end

function PCLoad(arg1::PC{Float64},arg2::PetscViewer{Float64})
    ccall((:PCLoad,petsc1),PetscErrorCode,(PC{Float64},PetscViewer{Float64}),arg1,arg2)
end

function PCAppendOptionsPrefix(arg1::PC{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PCAppendOptionsPrefix,petsc1),PetscErrorCode,(PC{Float64},Cstring),arg1,arg2)
end

function PCGetOptionsPrefix(arg1::PC{Float64},arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PCGetOptionsPrefix,petsc1),PetscErrorCode,(PC{Float64},Ptr{Ptr{UInt8}}),arg1,arg2)
end

function PCComputeExplicitOperator(arg1::PC{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:PCComputeExplicitOperator,petsc1),PetscErrorCode,(PC{Float64},Ptr{Mat{Float64}}),arg1,arg2)
end

function PCGetDiagonalScale(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PCGetDiagonalScale,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
end

function PCDiagonalScaleLeft(arg1::PC{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:PCDiagonalScaleLeft,petsc1),PetscErrorCode,(PC{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function PCDiagonalScaleRight(arg1::PC{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:PCDiagonalScaleRight,petsc1),PetscErrorCode,(PC{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function PCSetDiagonalScale(arg1::PC{Float64},arg2::Vec{Float64})
    ccall((:PCSetDiagonalScale,petsc1),PetscErrorCode,(PC{Float64},Vec{Float64}),arg1,arg2)
end

function PCJacobiSetType(arg1::PC{Float64},arg2::PCJacobiType)
    ccall((:PCJacobiSetType,petsc1),PetscErrorCode,(PC{Float64},PCJacobiType),arg1,arg2)
end

function PCJacobiGetType(arg1::PC{Float64},arg2::Union(Ptr{PCJacobiType},StridedArray{PCJacobiType},Ptr{Void}))
    ccall((:PCJacobiGetType,petsc1),PetscErrorCode,(PC{Float64},Ptr{PCJacobiType}),arg1,arg2)
end

function PCJacobiSetUseAbs(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCJacobiSetUseAbs,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCJacobiGetUseAbs(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PCJacobiGetUseAbs,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
end

function PCSORSetSymmetric(arg1::PC{Float64},arg2::MatSORType)
    ccall((:PCSORSetSymmetric,petsc1),PetscErrorCode,(PC{Float64},MatSORType),arg1,arg2)
end

function PCSORGetSymmetric(arg1::PC{Float64},arg2::Union(Ptr{MatSORType},StridedArray{MatSORType},Ptr{Void}))
    ccall((:PCSORGetSymmetric,petsc1),PetscErrorCode,(PC{Float64},Ptr{MatSORType}),arg1,arg2)
end

function PCSORSetOmega(arg1::PC{Float64},Float64::Integer)
    ccall((:PCSORSetOmega,petsc1),PetscErrorCode,(PC{Float64},Cint),arg1,PetscReal)
end

function PCSORGetOmega(arg1::PC{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PCSORGetOmega,petsc1),PetscErrorCode,(PC{Float64},Ptr{Cint}),arg1,arg2)
end

function PCSORSetIterations(arg1::PC{Float64},arg2::Integer,arg3::Integer)
    ccall((:PCSORSetIterations,petsc1),PetscErrorCode,(PC{Float64},Int32,Int32),arg1,arg2,arg3)
end

function PCSORGetIterations(arg1::PC{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PCSORGetIterations,petsc1),PetscErrorCode,(PC{Float64},Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3)
end

function PCEisenstatSetOmega(arg1::PC{Float64},Float64::Integer)
    ccall((:PCEisenstatSetOmega,petsc1),PetscErrorCode,(PC{Float64},Cint),arg1,PetscReal)
end

function PCEisenstatGetOmega(arg1::PC{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PCEisenstatGetOmega,petsc1),PetscErrorCode,(PC{Float64},Ptr{Cint}),arg1,arg2)
end

function PCEisenstatSetNoDiagonalScaling(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCEisenstatSetNoDiagonalScaling,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCEisenstatGetNoDiagonalScaling(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PCEisenstatGetNoDiagonalScaling,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
end

function PCBJacobiSetTotalBlocks(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PCBJacobiSetTotalBlocks,petsc1),PetscErrorCode,(PC{Float64},Int32,Ptr{Int32}),arg1,arg2,arg3)
end

function PCBJacobiGetTotalBlocks(arg1::PC{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:PCBJacobiGetTotalBlocks,petsc1),PetscErrorCode,(PC{Float64},Ptr{Int32},Ptr{Ptr{Int32}}),arg1,arg2,arg3)
end

function PCBJacobiSetLocalBlocks(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PCBJacobiSetLocalBlocks,petsc1),PetscErrorCode,(PC{Float64},Int32,Ptr{Int32}),arg1,arg2,arg3)
end

function PCBJacobiGetLocalBlocks(arg1::PC{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:PCBJacobiGetLocalBlocks,petsc1),PetscErrorCode,(PC{Float64},Ptr{Int32},Ptr{Ptr{Int32}}),arg1,arg2,arg3)
end

function PCShellSetApply(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PCShellSetApply,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
end

function PCShellSetApplyBA(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PCShellSetApplyBA,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
end

function PCShellSetApplyTranspose(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PCShellSetApplyTranspose,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
end

function PCShellSetSetUp(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PCShellSetSetUp,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
end

function PCShellSetApplyRichardson(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PCShellSetApplyRichardson,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
end

function PCShellSetView(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PCShellSetView,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
end

function PCShellSetDestroy(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PCShellSetDestroy,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
end

function PCShellSetContext(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PCShellSetContext,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
end

function PCShellGetContext(arg1::PC{Float64},arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:PCShellGetContext,petsc1),PetscErrorCode,(PC{Float64},Ptr{Ptr{Void}}),arg1,arg2)
end

function PCShellSetName(arg1::PC{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PCShellSetName,petsc1),PetscErrorCode,(PC{Float64},Cstring),arg1,arg2)
end

function PCShellGetName(arg1::PC{Float64},arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PCShellGetName,petsc1),PetscErrorCode,(PC{Float64},Ptr{Ptr{UInt8}}),arg1,arg2)
end

function PCFactorSetZeroPivot(arg1::PC{Float64},Float64::Integer)
    ccall((:PCFactorSetZeroPivot,petsc1),PetscErrorCode,(PC{Float64},Cint),arg1,PetscReal)
end

function PCFactorSetShiftType(arg1::PC{Float64},arg2::MatFactorShiftType)
    ccall((:PCFactorSetShiftType,petsc1),PetscErrorCode,(PC{Float64},MatFactorShiftType),arg1,arg2)
end

function PCFactorSetShiftAmount(arg1::PC{Float64},Float64::Integer)
    ccall((:PCFactorSetShiftAmount,petsc1),PetscErrorCode,(PC{Float64},Cint),arg1,PetscReal)
end

function PCFactorSetMatSolverPackage(arg1::PC{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PCFactorSetMatSolverPackage,petsc1),PetscErrorCode,(PC{Float64},Cstring),arg1,arg2)
end

function PCFactorGetMatSolverPackage(arg1::PC{Float64},arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PCFactorGetMatSolverPackage,petsc1),PetscErrorCode,(PC{Float64},Ptr{Ptr{UInt8}}),arg1,arg2)
end

function PCFactorSetUpMatSolverPackage(arg1::PC{Float64})
    ccall((:PCFactorSetUpMatSolverPackage,petsc1),PetscErrorCode,(PC{Float64},),arg1)
end

function PCFactorSetFill(arg1::PC{Float64},Float64::Integer)
    ccall((:PCFactorSetFill,petsc1),PetscErrorCode,(PC{Float64},Cint),arg1,PetscReal)
end

function PCFactorSetColumnPivot(arg1::PC{Float64},Float64::Integer)
    ccall((:PCFactorSetColumnPivot,petsc1),PetscErrorCode,(PC{Float64},Cint),arg1,PetscReal)
end

function PCFactorReorderForNonzeroDiagonal(arg1::PC{Float64},Float64::Integer)
    ccall((:PCFactorReorderForNonzeroDiagonal,petsc1),PetscErrorCode,(PC{Float64},Cint),arg1,PetscReal)
end

function PCFactorSetMatOrderingType(arg1::PC{Float64},arg2::MatOrderingType)
    ccall((:PCFactorSetMatOrderingType,petsc1),PetscErrorCode,(PC{Float64},Cstring),arg1,arg2)
end

function PCFactorSetReuseOrdering(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCFactorSetReuseOrdering,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCFactorSetReuseFill(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCFactorSetReuseFill,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCFactorSetUseInPlace(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCFactorSetUseInPlace,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCFactorGetUseInPlace(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PCFactorGetUseInPlace,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
end

function PCFactorSetAllowDiagonalFill(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCFactorSetAllowDiagonalFill,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCFactorGetAllowDiagonalFill(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PCFactorGetAllowDiagonalFill,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
end

function PCFactorSetPivotInBlocks(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCFactorSetPivotInBlocks,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCFactorSetLevels(arg1::PC{Float64},arg2::Integer)
    ccall((:PCFactorSetLevels,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCFactorGetLevels(arg1::PC{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PCFactorGetLevels,petsc1),PetscErrorCode,(PC{Float64},Ptr{Int32}),arg1,arg2)
end

function PCFactorSetDropTolerance(arg1::PC{Float64},Float64::Integer,arg2::Integer,arg3::Integer)
    ccall((:PCFactorSetDropTolerance,petsc1),PetscErrorCode,(PC{Float64},Cint,Cint,Int32),arg1,PetscReal,arg2,arg3)
end

function PCASMSetLocalSubdomains(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:PCASMSetLocalSubdomains,petsc1),PetscErrorCode,(PC{Float64},Int32,Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4)
end

function PCASMSetTotalSubdomains(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:PCASMSetTotalSubdomains,petsc1),PetscErrorCode,(PC{Float64},Int32,Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4)
end

function PCASMSetOverlap(arg1::PC{Float64},arg2::Integer)
    ccall((:PCASMSetOverlap,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCASMSetDMSubdomains(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCASMSetDMSubdomains,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCASMGetDMSubdomains(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PCASMGetDMSubdomains,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
end

function PCASMSetSortIndices(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCASMSetSortIndices,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCASMSetType(arg1::PC{Float64},arg2::PCASMType)
    ccall((:PCASMSetType,petsc1),PetscErrorCode,(PC{Float64},PCASMType),arg1,arg2)
end

function PCASMGetType(arg1::PC{Float64},arg2::Union(Ptr{PCASMType},StridedArray{PCASMType},Ptr{Void}))
    ccall((:PCASMGetType,petsc1),PetscErrorCode,(PC{Float64},Ptr{PCASMType}),arg1,arg2)
end

function PCASMSetLocalType(arg1::PC{Float64},arg2::PCCompositeType)
    ccall((:PCASMSetLocalType,petsc1),PetscErrorCode,(PC{Float64},PCCompositeType),arg1,arg2)
end

function PCASMGetLocalType(arg1::PC{Float64},arg2::Union(Ptr{PCCompositeType},StridedArray{PCCompositeType},Ptr{Void}))
    ccall((:PCASMGetLocalType,petsc1),PetscErrorCode,(PC{Float64},Ptr{PCCompositeType}),arg1,arg2)
end

function PCASMCreateSubdomains(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    ccall((:PCASMCreateSubdomains,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Ptr{IS{Float64}}}),arg1,arg2,arg3)
end

function PCASMDestroySubdomains(arg1::Integer,arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:PCASMDestroySubdomains,petsc1),PetscErrorCode,(Int32,Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3)
end

function PCASMCreateSubdomains2D(arg1::Integer,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg8::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}),arg9::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    ccall((:PCASMCreateSubdomains2D,petsc1),PetscErrorCode,(Int32,Int32,Int32,Int32,Int32,Int32,Ptr{Int32},Ptr{Ptr{IS{Float64}}},Ptr{Ptr{IS{Float64}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function PCASMGetLocalSubdomains(arg1::PC{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}),arg4::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    ccall((:PCASMGetLocalSubdomains,petsc1),PetscErrorCode,(PC{Float64},Ptr{Int32},Ptr{Ptr{IS{Float64}}},Ptr{Ptr{IS{Float64}}}),arg1,arg2,arg3,arg4)
end

function PCASMGetLocalSubmatrices(arg1::PC{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Ptr{Mat{Float64}}},StridedArray{Ptr{Mat{Float64}}},Ptr{Void}))
    ccall((:PCASMGetLocalSubmatrices,petsc1),PetscErrorCode,(PC{Float64},Ptr{Int32},Ptr{Ptr{Mat{Float64}}}),arg1,arg2,arg3)
end

function PCGASMSetTotalSubdomains(arg1::PC{Float64},arg2::Integer)
    ccall((:PCGASMSetTotalSubdomains,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCGASMSetSubdomains(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}),arg4::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:PCGASMSetSubdomains,petsc1),PetscErrorCode,(PC{Float64},Int32,Ptr{IS{Float64}},Ptr{IS{Float64}}),arg1,arg2,arg3,arg4)
end

function PCGASMSetOverlap(arg1::PC{Float64},arg2::Integer)
    ccall((:PCGASMSetOverlap,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCGASMSetUseDMSubdomains(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCGASMSetUseDMSubdomains,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCGASMGetUseDMSubdomains(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PCGASMGetUseDMSubdomains,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
end

function PCGASMSetSortIndices(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCGASMSetSortIndices,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCGASMSetType(arg1::PC{Float64},arg2::PCGASMType)
    ccall((:PCGASMSetType,petsc1),PetscErrorCode,(PC{Float64},PCGASMType),arg1,arg2)
end

function PCGASMCreateSubdomains(arg1::Mat{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    ccall((:PCGASMCreateSubdomains,petsc1),PetscErrorCode,(Mat{Float64},Int32,Ptr{Int32},Ptr{Ptr{IS{Float64}}}),arg1,arg2,arg3,arg4)
end

function PCGASMDestroySubdomains(arg1::Integer,arg2::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}),arg3::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    ccall((:PCGASMDestroySubdomains,petsc1),PetscErrorCode,(Int32,Ptr{Ptr{IS{Float64}}},Ptr{Ptr{IS{Float64}}}),arg1,arg2,arg3)
end

function PCGASMCreateSubdomains2D(arg1::PC{Float64},arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Integer,arg8::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg9::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}),arg10::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    ccall((:PCGASMCreateSubdomains2D,petsc1),PetscErrorCode,(PC{Float64},Int32,Int32,Int32,Int32,Int32,Int32,Ptr{Int32},Ptr{Ptr{IS{Float64}}},Ptr{Ptr{IS{Float64}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function PCGASMGetSubdomains(arg1::PC{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}),arg4::Union(Ptr{Ptr{IS{Float64}}},StridedArray{Ptr{IS{Float64}}},Ptr{Void}))
    ccall((:PCGASMGetSubdomains,petsc1),PetscErrorCode,(PC{Float64},Ptr{Int32},Ptr{Ptr{IS{Float64}}},Ptr{Ptr{IS{Float64}}}),arg1,arg2,arg3,arg4)
end

function PCGASMGetSubmatrices(arg1::PC{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Ptr{Mat{Float64}}},StridedArray{Ptr{Mat{Float64}}},Ptr{Void}))
    ccall((:PCGASMGetSubmatrices,petsc1),PetscErrorCode,(PC{Float64},Ptr{Int32},Ptr{Ptr{Mat{Float64}}}),arg1,arg2,arg3)
end

function PCCompositeSetType(arg1::PC{Float64},arg2::PCCompositeType)
    ccall((:PCCompositeSetType,petsc1),PetscErrorCode,(PC{Float64},PCCompositeType),arg1,arg2)
end

function PCCompositeGetType(arg1::PC{Float64},arg2::Union(Ptr{PCCompositeType},StridedArray{PCCompositeType},Ptr{Void}))
    ccall((:PCCompositeGetType,petsc1),PetscErrorCode,(PC{Float64},Ptr{PCCompositeType}),arg1,arg2)
end

function PCCompositeAddPC(arg1::PC{Float64},arg2::PCType)
    ccall((:PCCompositeAddPC,petsc1),PetscErrorCode,(PC{Float64},Cstring),arg1,arg2)
end

function PCCompositeGetPC(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{PC{Float64}},StridedArray{PC{Float64}},Ptr{Void}))
    ccall((:PCCompositeGetPC,petsc1),PetscErrorCode,(PC{Float64},Int32,Ptr{PC{Float64}}),arg1,arg2,arg3)
end

function PCCompositeSpecialSetAlpha(arg1::PC{Float64},arg2::Float64)
    ccall((:PCCompositeSpecialSetAlpha,petsc1),PetscErrorCode,(PC{Float64},Float64),arg1,arg2)
end

function PCRedundantSetNumber(arg1::PC{Float64},arg2::Integer)
    ccall((:PCRedundantSetNumber,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCRedundantSetScatter(arg1::PC{Float64},arg2::VecScatter{Float64},arg3::VecScatter{Float64})
    ccall((:PCRedundantSetScatter,petsc1),PetscErrorCode,(PC{Float64},VecScatter{Float64},VecScatter{Float64}),arg1,arg2,arg3)
end

function PCRedundantGetOperators(arg1::PC{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:PCRedundantGetOperators,petsc1),PetscErrorCode,(PC{Float64},Ptr{Mat{Float64}},Ptr{Mat{Float64}}),arg1,arg2,arg3)
end

function PCSPAISetEpsilon(arg1::PC{Float64},arg2::Cdouble)
    ccall((:PCSPAISetEpsilon,petsc1),PetscErrorCode,(PC{Float64},Cdouble),arg1,arg2)
end

function PCSPAISetNBSteps(arg1::PC{Float64},arg2::Integer)
    ccall((:PCSPAISetNBSteps,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCSPAISetMax(arg1::PC{Float64},arg2::Integer)
    ccall((:PCSPAISetMax,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCSPAISetMaxNew(arg1::PC{Float64},arg2::Integer)
    ccall((:PCSPAISetMaxNew,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCSPAISetBlockSize(arg1::PC{Float64},arg2::Integer)
    ccall((:PCSPAISetBlockSize,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCSPAISetCacheSize(arg1::PC{Float64},arg2::Integer)
    ccall((:PCSPAISetCacheSize,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCSPAISetVerbose(arg1::PC{Float64},arg2::Integer)
    ccall((:PCSPAISetVerbose,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCSPAISetSp(arg1::PC{Float64},arg2::Integer)
    ccall((:PCSPAISetSp,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCHYPRESetType(arg1::PC{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PCHYPRESetType,petsc1),PetscErrorCode,(PC{Float64},Cstring),arg1,arg2)
end

function PCHYPREGetType(arg1::PC{Float64},arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:PCHYPREGetType,petsc1),PetscErrorCode,(PC{Float64},Ptr{Ptr{UInt8}}),arg1,arg2)
end

function PCHYPRESetDiscreteGradient(arg1::PC{Float64},arg2::Mat{Float64})
    ccall((:PCHYPRESetDiscreteGradient,petsc1),PetscErrorCode,(PC{Float64},Mat{Float64}),arg1,arg2)
end

function PCHYPRESetDiscreteCurl(arg1::PC{Float64},arg2::Mat{Float64})
    ccall((:PCHYPRESetDiscreteCurl,petsc1),PetscErrorCode,(PC{Float64},Mat{Float64}),arg1,arg2)
end

function PCHYPRESetEdgeConstantVectors(arg1::PC{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64})
    ccall((:PCHYPRESetEdgeConstantVectors,petsc1),PetscErrorCode,(PC{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
end

function PCHYPRESetAlphaPoissonMatrix(arg1::PC{Float64},arg2::Mat{Float64})
    ccall((:PCHYPRESetAlphaPoissonMatrix,petsc1),PetscErrorCode,(PC{Float64},Mat{Float64}),arg1,arg2)
end

function PCHYPRESetBetaPoissonMatrix(arg1::PC{Float64},arg2::Mat{Float64})
    ccall((:PCHYPRESetBetaPoissonMatrix,petsc1),PetscErrorCode,(PC{Float64},Mat{Float64}),arg1,arg2)
end

function PCBJacobiGetLocalBlocks(arg1::PC{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:PCBJacobiGetLocalBlocks,petsc1),PetscErrorCode,(PC{Float64},Ptr{Int32},Ptr{Ptr{Int32}}),arg1,arg2,arg3)
end

function PCBJacobiGetTotalBlocks(arg1::PC{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Ptr{Int32}},StridedArray{Ptr{Int32}},Ptr{Void}))
    ccall((:PCBJacobiGetTotalBlocks,petsc1),PetscErrorCode,(PC{Float64},Ptr{Int32},Ptr{Ptr{Int32}}),arg1,arg2,arg3)
end

function PCFieldSplitSetFields(arg1::PC{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Integer,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PCFieldSplitSetFields,petsc1),PetscErrorCode,(PC{Float64},Cstring,Int32,Ptr{Int32},Ptr{Int32}),arg1,arg2,arg3,arg4,arg5)
end

function PCFieldSplitSetType(arg1::PC{Float64},arg2::PCCompositeType)
    ccall((:PCFieldSplitSetType,petsc1),PetscErrorCode,(PC{Float64},PCCompositeType),arg1,arg2)
end

function PCFieldSplitGetType(arg1::PC{Float64},arg2::Union(Ptr{PCCompositeType},StridedArray{PCCompositeType},Ptr{Void}))
    ccall((:PCFieldSplitGetType,petsc1),PetscErrorCode,(PC{Float64},Ptr{PCCompositeType}),arg1,arg2)
end

function PCFieldSplitSetBlockSize(arg1::PC{Float64},arg2::Integer)
    ccall((:PCFieldSplitSetBlockSize,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCFieldSplitSetIS(arg1::PC{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::IS{Float64})
    ccall((:PCFieldSplitSetIS,petsc1),PetscErrorCode,(PC{Float64},Cstring,IS{Float64}),arg1,arg2,arg3)
end

function PCFieldSplitGetIS(arg1::PC{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:PCFieldSplitGetIS,petsc1),PetscErrorCode,(PC{Float64},Cstring,Ptr{IS{Float64}}),arg1,arg2,arg3)
end

function PCFieldSplitSetDMSplits(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCFieldSplitSetDMSplits,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCFieldSplitGetDMSplits(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PCFieldSplitGetDMSplits,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
end

function PCFieldSplitSetDiagUseAmat(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCFieldSplitSetDiagUseAmat,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCFieldSplitGetDiagUseAmat(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PCFieldSplitGetDiagUseAmat,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
end

function PCFieldSplitSetOffDiagUseAmat(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCFieldSplitSetOffDiagUseAmat,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCFieldSplitGetOffDiagUseAmat(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PCFieldSplitGetOffDiagUseAmat,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
end

function PETSC_DEPRECATED(arg0::Type{Float64})
    ccall((:PETSC_DEPRECATED,petsc1),Cint,())
end

function PCFieldSplitSchurPrecondition(arg1::PC{Float64},arg2::PCFieldSplitSchurPreType,arg3::Mat{Float64})
    ccall((:PCFieldSplitSchurPrecondition,petsc1),PetscErrorCode,(PC{Float64},PCFieldSplitSchurPreType,Mat{Float64}),arg1,arg2,arg3)
end

function PCFieldSplitSetSchurPre(arg1::PC{Float64},arg2::PCFieldSplitSchurPreType,arg3::Mat{Float64})
    ccall((:PCFieldSplitSetSchurPre,petsc1),PetscErrorCode,(PC{Float64},PCFieldSplitSchurPreType,Mat{Float64}),arg1,arg2,arg3)
end

function PCFieldSplitGetSchurPre(arg1::PC{Float64},arg2::Union(Ptr{PCFieldSplitSchurPreType},StridedArray{PCFieldSplitSchurPreType},Ptr{Void}),arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:PCFieldSplitGetSchurPre,petsc1),PetscErrorCode,(PC{Float64},Ptr{PCFieldSplitSchurPreType},Ptr{Mat{Float64}}),arg1,arg2,arg3)
end

function PCFieldSplitSetSchurFactType(arg1::PC{Float64},arg2::PCFieldSplitSchurFactType)
    ccall((:PCFieldSplitSetSchurFactType,petsc1),PetscErrorCode,(PC{Float64},PCFieldSplitSchurFactType),arg1,arg2)
end

function PCFieldSplitGetSchurBlocks(arg1::PC{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg5::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:PCFieldSplitGetSchurBlocks,petsc1),PetscErrorCode,(PC{Float64},Ptr{Mat{Float64}},Ptr{Mat{Float64}},Ptr{Mat{Float64}},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4,arg5)
end

function PCFieldSplitSchurGetS(arg1::PC{Float64},S::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:PCFieldSplitSchurGetS,petsc1),PetscErrorCode,(PC{Float64},Ptr{Mat{Float64}}),arg1,S)
end

function PCFieldSplitSchurRestoreS(arg1::PC{Float64},S::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:PCFieldSplitSchurRestoreS,petsc1),PetscErrorCode,(PC{Float64},Ptr{Mat{Float64}}),arg1,S)
end

function PCGalerkinSetRestriction(arg1::PC{Float64},arg2::Mat{Float64})
    ccall((:PCGalerkinSetRestriction,petsc1),PetscErrorCode,(PC{Float64},Mat{Float64}),arg1,arg2)
end

function PCGalerkinSetInterpolation(arg1::PC{Float64},arg2::Mat{Float64})
    ccall((:PCGalerkinSetInterpolation,petsc1),PetscErrorCode,(PC{Float64},Mat{Float64}),arg1,arg2)
end

function PCSetCoordinates(arg1::PC{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:PCSetCoordinates,petsc1),PetscErrorCode,(PC{Float64},Int32,Int32,Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function PCPythonSetType(arg1::PC{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:PCPythonSetType,petsc1),PetscErrorCode,(PC{Float64},Cstring),arg1,arg2)
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
    ccall((:PCSetApplicationContext,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
end

function PCGetApplicationContext(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PCGetApplicationContext,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
end

function PCBiCGStabCUSPSetTolerance(arg1::PC{Float64},Float64::Integer)
    ccall((:PCBiCGStabCUSPSetTolerance,petsc1),PetscErrorCode,(PC{Float64},Cint),arg1,PetscReal)
end

function PCBiCGStabCUSPSetIterations(arg1::PC{Float64},arg2::Integer)
    ccall((:PCBiCGStabCUSPSetIterations,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCBiCGStabCUSPSetUseVerboseMonitor(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCBiCGStabCUSPSetUseVerboseMonitor,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCAINVCUSPSetDropTolerance(arg1::PC{Float64},Float64::Integer)
    ccall((:PCAINVCUSPSetDropTolerance,petsc1),PetscErrorCode,(PC{Float64},Cint),arg1,PetscReal)
end

function PCAINVCUSPUseScaling(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCAINVCUSPUseScaling,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCAINVCUSPSetNonzeros(arg1::PC{Float64},arg2::Integer)
    ccall((:PCAINVCUSPSetNonzeros,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCAINVCUSPSetLinParameter(arg1::PC{Float64},arg2::Integer)
    ccall((:PCAINVCUSPSetLinParameter,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCPARMSSetGlobal(arg1::PC{Float64},arg2::PCPARMSGlobalType)
    ccall((:PCPARMSSetGlobal,petsc1),PetscErrorCode,(PC{Float64},PCPARMSGlobalType),arg1,arg2)
end

function PCPARMSSetLocal(arg1::PC{Float64},arg2::PCPARMSLocalType)
    ccall((:PCPARMSSetLocal,petsc1),PetscErrorCode,(PC{Float64},PCPARMSLocalType),arg1,arg2)
end

function PCPARMSSetSolveTolerances(arg1::PC{Float64},Float64::Integer,arg2::Integer)
    ccall((:PCPARMSSetSolveTolerances,petsc1),PetscErrorCode,(PC{Float64},Cint,Int32),arg1,PetscReal,arg2)
end

function PCPARMSSetSolveRestart(arg1::PC{Float64},arg2::Integer)
    ccall((:PCPARMSSetSolveRestart,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCPARMSSetNonsymPerm(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCPARMSSetNonsymPerm,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCPARMSSetFill(arg1::PC{Float64},arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:PCPARMSSetFill,petsc1),PetscErrorCode,(PC{Float64},Int32,Int32,Int32),arg1,arg2,arg3,arg4)
end

function PCGAMGSetType(arg1::PC{Float64},arg2::PCGAMGType)
    ccall((:PCGAMGSetType,petsc1),PetscErrorCode,(PC{Float64},Cstring),arg1,arg2)
end

function PCGAMGGetType(arg1::PC{Float64},arg2::Union(Ptr{PCGAMGType},StridedArray{PCGAMGType},Ptr{Void}))
    ccall((:PCGAMGGetType,petsc1),PetscErrorCode,(PC{Float64},Ptr{Ptr{UInt8}}),arg1,arg2)
end

function PCGAMGSetProcEqLim(arg1::PC{Float64},arg2::Integer)
    ccall((:PCGAMGSetProcEqLim,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCGAMGSetRepartitioning(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCGAMGSetRepartitioning,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCGAMGSetUseASMAggs(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCGAMGSetUseASMAggs,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCGAMGSetSolverType(arg1::PC{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Integer)
    ccall((:PCGAMGSetSolverType,petsc1),PetscErrorCode,(PC{Float64},Cstring,Int32),arg1,arg2,arg3)
end

function PCGAMGSetThreshold(arg1::PC{Float64},Float64::Integer)
    ccall((:PCGAMGSetThreshold,petsc1),PetscErrorCode,(PC{Float64},Cint),arg1,PetscReal)
end

function PCGAMGSetCoarseEqLim(arg1::PC{Float64},arg2::Integer)
    ccall((:PCGAMGSetCoarseEqLim,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCGAMGSetNlevels(arg1::PC{Float64},arg2::Integer)
    ccall((:PCGAMGSetNlevels,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCGAMGSetNSmooths(arg1::PC{Float64},arg2::Integer)
    ccall((:PCGAMGSetNSmooths,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCGAMGSetSymGraph(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCGAMGSetSymGraph,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCGAMGSetSquareGraph(arg1::PC{Float64},arg2::Integer)
    ccall((:PCGAMGSetSquareGraph,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCGAMGSetReuseInterpolation(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCGAMGSetReuseInterpolation,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCGAMGFinalizePackage(arg0::Type{Float64})
    ccall((:PCGAMGFinalizePackage,petsc1),PetscErrorCode,())
end

function PCGAMGInitializePackage(arg0::Type{Float64})
    ccall((:PCGAMGInitializePackage,petsc1),PetscErrorCode,())
end

function PCGAMGRegister(arg0::Type{Float64},arg1::PCGAMGType,arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PCGAMGRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
end

function PCGAMGClassicalSetType(arg1::PC{Float64},arg2::PCGAMGClassicalType)
    ccall((:PCGAMGClassicalSetType,petsc1),PetscErrorCode,(PC{Float64},Cstring),arg1,arg2)
end

function PCGAMGClassicalGetType(arg1::PC{Float64},arg2::Union(Ptr{PCGAMGClassicalType},StridedArray{PCGAMGClassicalType},Ptr{Void}))
    ccall((:PCGAMGClassicalGetType,petsc1),PetscErrorCode,(PC{Float64},Ptr{Ptr{UInt8}}),arg1,arg2)
end

function PCBDDCSetChangeOfBasisMat(arg1::PC{Float64},arg2::Mat{Float64})
    ccall((:PCBDDCSetChangeOfBasisMat,petsc1),PetscErrorCode,(PC{Float64},Mat{Float64}),arg1,arg2)
end

function PCBDDCSetPrimalVerticesLocalIS(arg1::PC{Float64},arg2::IS{Float64})
    ccall((:PCBDDCSetPrimalVerticesLocalIS,petsc1),PetscErrorCode,(PC{Float64},IS{Float64}),arg1,arg2)
end

function PCBDDCSetCoarseningRatio(arg1::PC{Float64},arg2::Integer)
    ccall((:PCBDDCSetCoarseningRatio,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCBDDCSetLevels(arg1::PC{Float64},arg2::Integer)
    ccall((:PCBDDCSetLevels,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

#= skipping function with undefined symbols: 
 function PCBDDCSetNullSpace(arg1::PC{Float64},arg2::MatNullSpace)
    ccall((:PCBDDCSetNullSpace,petsc1),PetscErrorCode,(PC{Float64},MatNullSpace),arg1,arg2)
end 
=#
function PCBDDCSetDirichletBoundaries(arg1::PC{Float64},arg2::IS{Float64})
    ccall((:PCBDDCSetDirichletBoundaries,petsc1),PetscErrorCode,(PC{Float64},IS{Float64}),arg1,arg2)
end

function PCBDDCSetDirichletBoundariesLocal(arg1::PC{Float64},arg2::IS{Float64})
    ccall((:PCBDDCSetDirichletBoundariesLocal,petsc1),PetscErrorCode,(PC{Float64},IS{Float64}),arg1,arg2)
end

function PCBDDCGetDirichletBoundaries(arg1::PC{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:PCBDDCGetDirichletBoundaries,petsc1),PetscErrorCode,(PC{Float64},Ptr{IS{Float64}}),arg1,arg2)
end

function PCBDDCGetDirichletBoundariesLocal(arg1::PC{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:PCBDDCGetDirichletBoundariesLocal,petsc1),PetscErrorCode,(PC{Float64},Ptr{IS{Float64}}),arg1,arg2)
end

function PCBDDCSetNeumannBoundaries(arg1::PC{Float64},arg2::IS{Float64})
    ccall((:PCBDDCSetNeumannBoundaries,petsc1),PetscErrorCode,(PC{Float64},IS{Float64}),arg1,arg2)
end

function PCBDDCSetNeumannBoundariesLocal(arg1::PC{Float64},arg2::IS{Float64})
    ccall((:PCBDDCSetNeumannBoundariesLocal,petsc1),PetscErrorCode,(PC{Float64},IS{Float64}),arg1,arg2)
end

function PCBDDCGetNeumannBoundaries(arg1::PC{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:PCBDDCGetNeumannBoundaries,petsc1),PetscErrorCode,(PC{Float64},Ptr{IS{Float64}}),arg1,arg2)
end

function PCBDDCGetNeumannBoundariesLocal(arg1::PC{Float64},arg2::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:PCBDDCGetNeumannBoundariesLocal,petsc1),PetscErrorCode,(PC{Float64},Ptr{IS{Float64}}),arg1,arg2)
end

function PCBDDCSetDofsSplitting(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:PCBDDCSetDofsSplitting,petsc1),PetscErrorCode,(PC{Float64},Int32,Ptr{IS{Float64}}),arg1,arg2,arg3)
end

function PCBDDCSetDofsSplittingLocal(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{Void}))
    ccall((:PCBDDCSetDofsSplittingLocal,petsc1),PetscErrorCode,(PC{Float64},Int32,Ptr{IS{Float64}}),arg1,arg2,arg3)
end

function PCBDDCSetLocalAdjacencyGraph(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),PetscCopyMode::Integer)
    ccall((:PCBDDCSetLocalAdjacencyGraph,petsc1),PetscErrorCode,(PC{Float64},Int32,Ptr{Int32},Ptr{Int32},Cint),arg1,arg2,arg3,arg4,PetscCopyMode)
end

function PCBDDCCreateFETIDPOperators(arg1::PC{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg3::Union(Ptr{PC{Float64}},StridedArray{PC{Float64}},Ptr{Void}))
    ccall((:PCBDDCCreateFETIDPOperators,petsc1),PetscErrorCode,(PC{Float64},Ptr{Mat{Float64}},Ptr{PC{Float64}}),arg1,arg2,arg3)
end

function PCBDDCMatFETIDPGetRHS(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:PCBDDCMatFETIDPGetRHS,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function PCBDDCMatFETIDPGetSolution(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:PCBDDCMatFETIDPGetSolution,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function PCISSetUseStiffnessScaling(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCISSetUseStiffnessScaling,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCISSetSubdomainScalingFactor(arg1::PC{Float64},arg2::Float64)
    ccall((:PCISSetSubdomainScalingFactor,petsc1),PetscErrorCode,(PC{Float64},Float64),arg1,arg2)
end

function PCISSetSubdomainDiagonalScaling(arg1::PC{Float64},arg2::Vec{Float64})
    ccall((:PCISSetSubdomainDiagonalScaling,petsc1),PetscErrorCode,(PC{Float64},Vec{Float64}),arg1,arg2)
end

function PCMGSetType(arg1::PC{Float64},arg2::PCMGType)
    ccall((:PCMGSetType,petsc1),PetscErrorCode,(PC{Float64},PCMGType),arg1,arg2)
end

function PCMGGetType(arg1::PC{Float64},arg2::Union(Ptr{PCMGType},StridedArray{PCMGType},Ptr{Void}))
    ccall((:PCMGGetType,petsc1),PetscErrorCode,(PC{Float64},Ptr{PCMGType}),arg1,arg2)
end

function PCMGSetLevels(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{MPI_Comm},StridedArray{MPI_Comm},Ptr{Void}))
    ccall((:PCMGSetLevels,petsc1),PetscErrorCode,(PC{Float64},Int32,Ptr{comm_type}),arg1,arg2,arg3)
end

function PCMGGetLevels(arg1::PC{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PCMGGetLevels,petsc1),PetscErrorCode,(PC{Float64},Ptr{Int32}),arg1,arg2)
end

function PCMGSetNumberSmoothUp(arg1::PC{Float64},arg2::Integer)
    ccall((:PCMGSetNumberSmoothUp,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCMGSetNumberSmoothDown(arg1::PC{Float64},arg2::Integer)
    ccall((:PCMGSetNumberSmoothDown,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCMGSetCycleType(arg1::PC{Float64},arg2::PCMGCycleType)
    ccall((:PCMGSetCycleType,petsc1),PetscErrorCode,(PC{Float64},PCMGCycleType),arg1,arg2)
end

function PCMGSetCycleTypeOnLevel(arg1::PC{Float64},arg2::Integer,arg3::PCMGCycleType)
    ccall((:PCMGSetCycleTypeOnLevel,petsc1),PetscErrorCode,(PC{Float64},Int32,PCMGCycleType),arg1,arg2,arg3)
end

function PCMGSetCyclesOnLevel(arg1::PC{Float64},arg2::Integer,arg3::Integer)
    ccall((:PCMGSetCyclesOnLevel,petsc1),PetscErrorCode,(PC{Float64},Int32,Int32),arg1,arg2,arg3)
end

function PCMGMultiplicativeSetCycles(arg1::PC{Float64},arg2::Integer)
    ccall((:PCMGMultiplicativeSetCycles,petsc1),PetscErrorCode,(PC{Float64},Int32),arg1,arg2)
end

function PCMGSetGalerkin(arg1::PC{Float64},arg2::PetscBool)
    ccall((:PCMGSetGalerkin,petsc1),PetscErrorCode,(PC{Float64},PetscBool),arg1,arg2)
end

function PCMGGetGalerkin(arg1::PC{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:PCMGGetGalerkin,petsc1),PetscErrorCode,(PC{Float64},Ptr{PetscBool}),arg1,arg2)
end

function PCMGSetRhs(arg1::PC{Float64},arg2::Integer,arg3::Vec{Float64})
    ccall((:PCMGSetRhs,petsc1),PetscErrorCode,(PC{Float64},Int32,Vec{Float64}),arg1,arg2,arg3)
end

function PCMGSetX(arg1::PC{Float64},arg2::Integer,arg3::Vec{Float64})
    ccall((:PCMGSetX,petsc1),PetscErrorCode,(PC{Float64},Int32,Vec{Float64}),arg1,arg2,arg3)
end

function PCMGSetR(arg1::PC{Float64},arg2::Integer,arg3::Vec{Float64})
    ccall((:PCMGSetR,petsc1),PetscErrorCode,(PC{Float64},Int32,Vec{Float64}),arg1,arg2,arg3)
end

function PCMGSetRestriction(arg1::PC{Float64},arg2::Integer,arg3::Mat{Float64})
    ccall((:PCMGSetRestriction,petsc1),PetscErrorCode,(PC{Float64},Int32,Mat{Float64}),arg1,arg2,arg3)
end

function PCMGGetRestriction(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:PCMGGetRestriction,petsc1),PetscErrorCode,(PC{Float64},Int32,Ptr{Mat{Float64}}),arg1,arg2,arg3)
end

function PCMGSetInterpolation(arg1::PC{Float64},arg2::Integer,arg3::Mat{Float64})
    ccall((:PCMGSetInterpolation,petsc1),PetscErrorCode,(PC{Float64},Int32,Mat{Float64}),arg1,arg2,arg3)
end

function PCMGGetInterpolation(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:PCMGGetInterpolation,petsc1),PetscErrorCode,(PC{Float64},Int32,Ptr{Mat{Float64}}),arg1,arg2,arg3)
end

function PCMGSetRScale(arg1::PC{Float64},arg2::Integer,arg3::Vec{Float64})
    ccall((:PCMGSetRScale,petsc1),PetscErrorCode,(PC{Float64},Int32,Vec{Float64}),arg1,arg2,arg3)
end

function PCMGGetRScale(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:PCMGGetRScale,petsc1),PetscErrorCode,(PC{Float64},Int32,Ptr{Vec{Float64}}),arg1,arg2,arg3)
end

function PCMGSetResidual(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Mat{Float64})
    ccall((:PCMGSetResidual,petsc1),PetscErrorCode,(PC{Float64},Int32,Ptr{Void},Mat{Float64}),arg1,arg2,arg3,arg4)
end

function PCMGResidualDefault(arg1::Mat{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64})
    ccall((:PCMGResidualDefault,petsc1),PetscErrorCode,(Mat{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4)
end

function KSPInitializePackage(arg0::Type{Float64})
    ccall((:KSPInitializePackage,petsc1),PetscErrorCode,())
end

function KSPCreate(arg1::MPI_Comm,arg2::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    ccall((:KSPCreate,petsc1),PetscErrorCode,(comm_type,Ptr{KSP{Float64}}),arg1.val,arg2)
end

function KSPSetType(arg1::KSP{Float64},arg2::KSPType)
    ccall((:KSPSetType,petsc1),PetscErrorCode,(KSP{Float64},Cstring),arg1,arg2)
end

function KSPGetType(arg1::KSP{Float64},arg2::Union(Ptr{KSPType},StridedArray{KSPType},Ptr{Void}))
    ccall((:KSPGetType,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Ptr{UInt8}}),arg1,arg2)
end

function KSPSetUp(arg1::KSP{Float64})
    ccall((:KSPSetUp,petsc1),PetscErrorCode,(KSP{Float64},),arg1)
end

function KSPSetUpOnBlocks(arg1::KSP{Float64})
    ccall((:KSPSetUpOnBlocks,petsc1),PetscErrorCode,(KSP{Float64},),arg1)
end

function KSPSolve(arg1::KSP{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:KSPSolve,petsc1),PetscErrorCode,(KSP{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function KSPSolveTranspose(arg1::KSP{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:KSPSolveTranspose,petsc1),PetscErrorCode,(KSP{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function KSPReset(arg1::KSP{Float64})
    ccall((:KSPReset,petsc1),PetscErrorCode,(KSP{Float64},),arg1)
end

function KSPDestroy(arg1::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    ccall((:KSPDestroy,petsc1),PetscErrorCode,(Ptr{KSP{Float64}},),arg1)
end

function KSPSetReusePreconditioner(arg1::KSP{Float64},arg2::PetscBool)
    ccall((:KSPSetReusePreconditioner,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
end

function KSPSetSkipPCSetFromOptions(arg1::KSP{Float64},arg2::PetscBool)
    ccall((:KSPSetSkipPCSetFromOptions,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
end

function KSPRegister(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPRegister,petsc1),PetscErrorCode,(Cstring,Ptr{Void}),arg1,arg2)
end

function KSPSetPCSide(arg1::KSP{Float64},arg2::PCSide)
    ccall((:KSPSetPCSide,petsc1),PetscErrorCode,(KSP{Float64},PCSide),arg1,arg2)
end

function KSPGetPCSide(arg1::KSP{Float64},arg2::Union(Ptr{PCSide},StridedArray{PCSide},Ptr{Void}))
    ccall((:KSPGetPCSide,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PCSide}),arg1,arg2)
end

function KSPSetTolerances(arg1::KSP{Float64},Float64::Integer,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:KSPSetTolerances,petsc1),PetscErrorCode,(KSP{Float64},Cint,Cint,Cint,Int32),arg1,PetscReal,arg2,arg3,arg4)
end

function KSPGetTolerances(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:KSPGetTolerances,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Int32}),arg1,arg2,arg3,arg4,arg5)
end

function KSPSetInitialGuessNonzero(arg1::KSP{Float64},arg2::PetscBool)
    ccall((:KSPSetInitialGuessNonzero,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
end

function KSPGetInitialGuessNonzero(arg1::KSP{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:KSPGetInitialGuessNonzero,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PetscBool}),arg1,arg2)
end

function KSPSetInitialGuessKnoll(arg1::KSP{Float64},arg2::PetscBool)
    ccall((:KSPSetInitialGuessKnoll,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
end

function KSPGetInitialGuessKnoll(arg1::KSP{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:KSPGetInitialGuessKnoll,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PetscBool}),arg1,arg2)
end

function KSPSetErrorIfNotConverged(arg1::KSP{Float64},arg2::PetscBool)
    ccall((:KSPSetErrorIfNotConverged,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
end

function KSPGetErrorIfNotConverged(arg1::KSP{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:KSPGetErrorIfNotConverged,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PetscBool}),arg1,arg2)
end

function KSPSetComputeEigenvalues(arg1::KSP{Float64},arg2::PetscBool)
    ccall((:KSPSetComputeEigenvalues,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
end

function KSPGetComputeEigenvalues(arg1::KSP{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:KSPGetComputeEigenvalues,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PetscBool}),arg1,arg2)
end

function KSPSetComputeSingularValues(arg1::KSP{Float64},arg2::PetscBool)
    ccall((:KSPSetComputeSingularValues,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
end

function KSPGetComputeSingularValues(arg1::KSP{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:KSPGetComputeSingularValues,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PetscBool}),arg1,arg2)
end

function KSPGetRhs(arg1::KSP{Float64},arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:KSPGetRhs,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Vec{Float64}}),arg1,arg2)
end

function KSPGetSolution(arg1::KSP{Float64},arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:KSPGetSolution,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Vec{Float64}}),arg1,arg2)
end

function KSPGetResidualNorm(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:KSPGetResidualNorm,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint}),arg1,arg2)
end

function KSPGetIterationNumber(arg1::KSP{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:KSPGetIterationNumber,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Int32}),arg1,arg2)
end

function KSPGetTotalIterations(arg1::KSP{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:KSPGetTotalIterations,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Int32}),arg1,arg2)
end

function KSPCreateVecs(arg1::KSP{Float64},arg2::Integer,arg3::Union(Ptr{Ptr{Vec{Float64}}},StridedArray{Ptr{Vec{Float64}}},Ptr{Void}),arg4::Integer,arg5::Union(Ptr{Ptr{Vec{Float64}}},StridedArray{Ptr{Vec{Float64}}},Ptr{Void}))
    ccall((:KSPCreateVecs,petsc1),PetscErrorCode,(KSP{Float64},Int32,Ptr{Ptr{Vec{Float64}}},Int32,Ptr{Ptr{Vec{Float64}}}),arg1,arg2,arg3,arg4,arg5)
end

function KSPSetPostSolve(arg1::KSP{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPSetPostSolve,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function KSPSetPC(arg1::KSP{Float64},arg2::PC{Float64})
    ccall((:KSPSetPC,petsc1),PetscErrorCode,(KSP{Float64},PC{Float64}),arg1,arg2)
end

function KSPGetPC(arg1::KSP{Float64},arg2::Union(Ptr{PC{Float64}},StridedArray{PC{Float64}},Ptr{Void}))
    ccall((:KSPGetPC,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PC{Float64}}),arg1,arg2)
end

function KSPMonitor(arg1::KSP{Float64},arg2::Integer,Float64::Integer)
    ccall((:KSPMonitor,petsc1),PetscErrorCode,(KSP{Float64},Int32,Cint),arg1,arg2,PetscReal)
end

function KSPMonitorSet(arg1::KSP{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPMonitorSet,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function KSPMonitorCancel(arg1::KSP{Float64})
    ccall((:KSPMonitorCancel,petsc1),PetscErrorCode,(KSP{Float64},),arg1)
end

function KSPGetMonitorContext(arg1::KSP{Float64},arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:KSPGetMonitorContext,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Ptr{Void}}),arg1,arg2)
end

function KSPGetResidualHistory(arg1::KSP{Float64},arg2::Union(Ptr{Ptr{Cint}},StridedArray{Ptr{Cint}},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:KSPGetResidualHistory,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Ptr{Cint}},Ptr{Int32}),arg1,arg2,arg3)
end

function KSPSetResidualHistory(arg1::KSP{Float64},Float64::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg2::Integer,arg3::PetscBool)
    ccall((:KSPSetResidualHistory,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint},Int32,PetscBool),arg1,PetscReal,arg2,arg3)
end

function KSPBuildSolutionDefault(arg1::KSP{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:KSPBuildSolutionDefault,petsc1),PetscErrorCode,(KSP{Float64},Vec{Float64},Ptr{Vec{Float64}}),arg1,arg2,arg3)
end

function KSPBuildResidualDefault(arg1::KSP{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:KSPBuildResidualDefault,petsc1),PetscErrorCode,(KSP{Float64},Vec{Float64},Vec{Float64},Ptr{Vec{Float64}}),arg1,arg2,arg3,arg4)
end

function KSPDestroyDefault(arg1::KSP{Float64})
    ccall((:KSPDestroyDefault,petsc1),PetscErrorCode,(KSP{Float64},),arg1)
end

function KSPSetWorkVecs(arg1::KSP{Float64},arg2::Integer)
    ccall((:KSPSetWorkVecs,petsc1),PetscErrorCode,(KSP{Float64},Int32),arg1,arg2)
end

function PCKSPGetKSP(arg1::PC{Float64},arg2::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    ccall((:PCKSPGetKSP,petsc1),PetscErrorCode,(PC{Float64},Ptr{KSP{Float64}}),arg1,arg2)
end

function PCBJacobiGetSubKSP(arg1::PC{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Ptr{KSP{Float64}}},StridedArray{Ptr{KSP{Float64}}},Ptr{Void}))
    ccall((:PCBJacobiGetSubKSP,petsc1),PetscErrorCode,(PC{Float64},Ptr{Int32},Ptr{Int32},Ptr{Ptr{KSP{Float64}}}),arg1,arg2,arg3,arg4)
end

function PCASMGetSubKSP(arg1::PC{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Ptr{KSP{Float64}}},StridedArray{Ptr{KSP{Float64}}},Ptr{Void}))
    ccall((:PCASMGetSubKSP,petsc1),PetscErrorCode,(PC{Float64},Ptr{Int32},Ptr{Int32},Ptr{Ptr{KSP{Float64}}}),arg1,arg2,arg3,arg4)
end

function PCGASMGetSubKSP(arg1::PC{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Ptr{KSP{Float64}}},StridedArray{Ptr{KSP{Float64}}},Ptr{Void}))
    ccall((:PCGASMGetSubKSP,petsc1),PetscErrorCode,(PC{Float64},Ptr{Int32},Ptr{Int32},Ptr{Ptr{KSP{Float64}}}),arg1,arg2,arg3,arg4)
end

function PCFieldSplitGetSubKSP(arg1::PC{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Ptr{KSP{Float64}}},StridedArray{Ptr{KSP{Float64}}},Ptr{Void}))
    ccall((:PCFieldSplitGetSubKSP,petsc1),PetscErrorCode,(PC{Float64},Ptr{Int32},Ptr{Ptr{KSP{Float64}}}),arg1,arg2,arg3)
end

function PCMGGetSmoother(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    ccall((:PCMGGetSmoother,petsc1),PetscErrorCode,(PC{Float64},Int32,Ptr{KSP{Float64}}),arg1,arg2,arg3)
end

function PCMGGetSmootherDown(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    ccall((:PCMGGetSmootherDown,petsc1),PetscErrorCode,(PC{Float64},Int32,Ptr{KSP{Float64}}),arg1,arg2,arg3)
end

function PCMGGetSmootherUp(arg1::PC{Float64},arg2::Integer,arg3::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    ccall((:PCMGGetSmootherUp,petsc1),PetscErrorCode,(PC{Float64},Int32,Ptr{KSP{Float64}}),arg1,arg2,arg3)
end

function PCMGGetCoarseSolve(arg1::PC{Float64},arg2::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    ccall((:PCMGGetCoarseSolve,petsc1),PetscErrorCode,(PC{Float64},Ptr{KSP{Float64}}),arg1,arg2)
end

function PCGalerkinGetKSP(arg1::PC{Float64},arg2::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    ccall((:PCGalerkinGetKSP,petsc1),PetscErrorCode,(PC{Float64},Ptr{KSP{Float64}}),arg1,arg2)
end

function KSPBuildSolution(arg1::KSP{Float64},arg2::Vec{Float64},arg3::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:KSPBuildSolution,petsc1),PetscErrorCode,(KSP{Float64},Vec{Float64},Ptr{Vec{Float64}}),arg1,arg2,arg3)
end

function KSPBuildResidual(arg1::KSP{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:KSPBuildResidual,petsc1),PetscErrorCode,(KSP{Float64},Vec{Float64},Vec{Float64},Ptr{Vec{Float64}}),arg1,arg2,arg3,arg4)
end

function KSPRichardsonSetScale(arg1::KSP{Float64},Float64::Integer)
    ccall((:KSPRichardsonSetScale,petsc1),PetscErrorCode,(KSP{Float64},Cint),arg1,PetscReal)
end

function KSPRichardsonSetSelfScale(arg1::KSP{Float64},arg2::PetscBool)
    ccall((:KSPRichardsonSetSelfScale,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
end

function KSPChebyshevSetEigenvalues(arg1::KSP{Float64},Float64::Integer,arg2::Integer)
    ccall((:KSPChebyshevSetEigenvalues,petsc1),PetscErrorCode,(KSP{Float64},Cint,Cint),arg1,PetscReal,arg2)
end

function KSPChebyshevEstEigSet(arg1::KSP{Float64},Float64::Integer,arg2::Integer,arg3::Integer,arg4::Integer)
    ccall((:KSPChebyshevEstEigSet,petsc1),PetscErrorCode,(KSP{Float64},Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4)
end

#= skipping function with undefined symbols: 
 function KSPChebyshevEstEigSetRandom(arg1::KSP{Float64},arg2::PetscRandom)
    ccall((:KSPChebyshevEstEigSetRandom,petsc1),PetscErrorCode,(KSP{Float64},PetscRandom),arg1,arg2)
end 
=#
function KSPChebyshevEstEigGetKSP(arg1::KSP{Float64},arg2::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    ccall((:KSPChebyshevEstEigGetKSP,petsc1),PetscErrorCode,(KSP{Float64},Ptr{KSP{Float64}}),arg1,arg2)
end

function KSPComputeExtremeSingularValues(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:KSPComputeExtremeSingularValues,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3)
end

function KSPComputeEigenvalues(arg1::KSP{Float64},arg2::Integer,Float64::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:KSPComputeEigenvalues,petsc1),PetscErrorCode,(KSP{Float64},Int32,Ptr{Cint},Ptr{Cint},Ptr{Int32}),arg1,arg2,PetscReal,arg3,arg4)
end

function KSPComputeEigenvaluesExplicitly(arg1::KSP{Float64},arg2::Integer,Float64::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}),arg3::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:KSPComputeEigenvaluesExplicitly,petsc1),PetscErrorCode,(KSP{Float64},Int32,Ptr{Cint},Ptr{Cint}),arg1,arg2,PetscReal,arg3)
end

function KSPFCGSetMmax(arg1::KSP{Float64},arg2::Integer)
    ccall((:KSPFCGSetMmax,petsc1),PetscErrorCode,(KSP{Float64},Int32),arg1,arg2)
end

function KSPFCGGetMmax(arg1::KSP{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:KSPFCGGetMmax,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Int32}),arg1,arg2)
end

function KSPFCGSetNprealloc(arg1::KSP{Float64},arg2::Integer)
    ccall((:KSPFCGSetNprealloc,petsc1),PetscErrorCode,(KSP{Float64},Int32),arg1,arg2)
end

function KSPFCGGetNprealloc(arg1::KSP{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:KSPFCGGetNprealloc,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Int32}),arg1,arg2)
end

function KSPFCGSetTruncationType(arg1::KSP{Float64},arg2::KSPFCGTruncationType)
    ccall((:KSPFCGSetTruncationType,petsc1),PetscErrorCode,(KSP{Float64},KSPFCGTruncationType),arg1,arg2)
end

function KSPFCGGetTruncationType(arg1::KSP{Float64},arg2::Union(Ptr{KSPFCGTruncationType},StridedArray{KSPFCGTruncationType},Ptr{Void}))
    ccall((:KSPFCGGetTruncationType,petsc1),PetscErrorCode,(KSP{Float64},Ptr{KSPFCGTruncationType}),arg1,arg2)
end

function KSPGMRESSetRestart(arg1::KSP{Float64},arg2::Integer)
    ccall((:KSPGMRESSetRestart,petsc1),PetscErrorCode,(KSP{Float64},Int32),arg1,arg2)
end

function KSPGMRESGetRestart(arg1::KSP{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:KSPGMRESGetRestart,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Int32}),arg1,arg2)
end

function KSPGMRESSetHapTol(arg1::KSP{Float64},Float64::Integer)
    ccall((:KSPGMRESSetHapTol,petsc1),PetscErrorCode,(KSP{Float64},Cint),arg1,PetscReal)
end

function KSPGMRESSetPreAllocateVectors(arg1::KSP{Float64})
    ccall((:KSPGMRESSetPreAllocateVectors,petsc1),PetscErrorCode,(KSP{Float64},),arg1)
end

function KSPGMRESSetOrthogonalization(arg1::KSP{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPGMRESSetOrthogonalization,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Void}),arg1,arg2)
end

function KSPGMRESGetOrthogonalization(arg1::KSP{Float64},arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:KSPGMRESGetOrthogonalization,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Ptr{Void}}),arg1,arg2)
end

function KSPGMRESModifiedGramSchmidtOrthogonalization(arg1::KSP{Float64},arg2::Integer)
    ccall((:KSPGMRESModifiedGramSchmidtOrthogonalization,petsc1),PetscErrorCode,(KSP{Float64},Int32),arg1,arg2)
end

function KSPGMRESClassicalGramSchmidtOrthogonalization(arg1::KSP{Float64},arg2::Integer)
    ccall((:KSPGMRESClassicalGramSchmidtOrthogonalization,petsc1),PetscErrorCode,(KSP{Float64},Int32),arg1,arg2)
end

function KSPLGMRESSetAugDim(arg1::KSP{Float64},arg2::Integer)
    ccall((:KSPLGMRESSetAugDim,petsc1),PetscErrorCode,(KSP{Float64},Int32),arg1,arg2)
end

function KSPLGMRESSetConstant(arg1::KSP{Float64})
    ccall((:KSPLGMRESSetConstant,petsc1),PetscErrorCode,(KSP{Float64},),arg1)
end

function KSPGCRSetRestart(arg1::KSP{Float64},arg2::Integer)
    ccall((:KSPGCRSetRestart,petsc1),PetscErrorCode,(KSP{Float64},Int32),arg1,arg2)
end

function KSPGCRGetRestart(arg1::KSP{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:KSPGCRGetRestart,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Int32}),arg1,arg2)
end

function KSPGCRSetModifyPC(arg1::KSP{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPGCRSetModifyPC,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function KSPGMRESSetCGSRefinementType(arg1::KSP{Float64},arg2::KSPGMRESCGSRefinementType)
    ccall((:KSPGMRESSetCGSRefinementType,petsc1),PetscErrorCode,(KSP{Float64},KSPGMRESCGSRefinementType),arg1,arg2)
end

function KSPGMRESGetCGSRefinementType(arg1::KSP{Float64},arg2::Union(Ptr{KSPGMRESCGSRefinementType},StridedArray{KSPGMRESCGSRefinementType},Ptr{Void}))
    ccall((:KSPGMRESGetCGSRefinementType,petsc1),PetscErrorCode,(KSP{Float64},Ptr{KSPGMRESCGSRefinementType}),arg1,arg2)
end

function KSPFGMRESModifyPCNoChange(arg1::KSP{Float64},arg2::Integer,arg3::Integer,Float64::Integer,arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPFGMRESModifyPCNoChange,petsc1),PetscErrorCode,(KSP{Float64},Int32,Int32,Cint,Ptr{Void}),arg1,arg2,arg3,PetscReal,arg4)
end

function KSPFGMRESModifyPCKSP(arg1::KSP{Float64},arg2::Integer,arg3::Integer,Float64::Integer,arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPFGMRESModifyPCKSP,petsc1),PetscErrorCode,(KSP{Float64},Int32,Int32,Cint,Ptr{Void}),arg1,arg2,arg3,PetscReal,arg4)
end

function KSPFGMRESSetModifyPC(arg1::KSP{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPFGMRESSetModifyPC,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function KSPQCGSetTrustRegionRadius(arg1::KSP{Float64},Float64::Integer)
    ccall((:KSPQCGSetTrustRegionRadius,petsc1),PetscErrorCode,(KSP{Float64},Cint),arg1,PetscReal)
end

function KSPQCGGetQuadratic(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:KSPQCGGetQuadratic,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint}),arg1,arg2)
end

function KSPQCGGetTrialStepNorm(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:KSPQCGGetTrialStepNorm,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint}),arg1,arg2)
end

function KSPBCGSLSetXRes(arg1::KSP{Float64},Float64::Integer)
    ccall((:KSPBCGSLSetXRes,petsc1),PetscErrorCode,(KSP{Float64},Cint),arg1,PetscReal)
end

function KSPBCGSLSetPol(arg1::KSP{Float64},arg2::PetscBool)
    ccall((:KSPBCGSLSetPol,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
end

function KSPBCGSLSetEll(arg1::KSP{Float64},arg2::Integer)
    ccall((:KSPBCGSLSetEll,petsc1),PetscErrorCode,(KSP{Float64},Int32),arg1,arg2)
end

function KSPBCGSLSetUsePseudoinverse(arg1::KSP{Float64},arg2::PetscBool)
    ccall((:KSPBCGSLSetUsePseudoinverse,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
end

function KSPSetFromOptions(arg1::KSP{Float64})
    ccall((:KSPSetFromOptions,petsc1),PetscErrorCode,(KSP{Float64},),arg1)
end

function KSPAddOptionsChecker(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPAddOptionsChecker,petsc1),PetscErrorCode,(Ptr{Void},),arg1)
end

function KSPMonitorSingularValue(arg1::KSP{Float64},arg2::Integer,Float64::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPMonitorSingularValue,petsc1),PetscErrorCode,(KSP{Float64},Int32,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorDefault(arg1::KSP{Float64},arg2::Integer,Float64::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPMonitorDefault,petsc1),PetscErrorCode,(KSP{Float64},Int32,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPLSQRMonitorDefault(arg1::KSP{Float64},arg2::Integer,Float64::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPLSQRMonitorDefault,petsc1),PetscErrorCode,(KSP{Float64},Int32,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorRange(arg1::KSP{Float64},arg2::Integer,Float64::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPMonitorRange,petsc1),PetscErrorCode,(KSP{Float64},Int32,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorDynamicTolerance(ksp::KSP{Float64},its::Integer,fnorm::Integer,dummy::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPMonitorDynamicTolerance,petsc1),PetscErrorCode,(KSP{Float64},Int32,Cint,Ptr{Void}),ksp,its,fnorm,dummy)
end

function KSPMonitorDynamicToleranceDestroy(arg0::Type{Float64},dummy::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:KSPMonitorDynamicToleranceDestroy,petsc1),PetscErrorCode,(Ptr{Ptr{Void}},),dummy)
end

function KSPMonitorTrueResidualNorm(arg1::KSP{Float64},arg2::Integer,Float64::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPMonitorTrueResidualNorm,petsc1),PetscErrorCode,(KSP{Float64},Int32,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorTrueResidualMaxNorm(arg1::KSP{Float64},arg2::Integer,Float64::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPMonitorTrueResidualMaxNorm,petsc1),PetscErrorCode,(KSP{Float64},Int32,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorDefaultShort(arg1::KSP{Float64},arg2::Integer,Float64::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPMonitorDefaultShort,petsc1),PetscErrorCode,(KSP{Float64},Int32,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorSolution(arg1::KSP{Float64},arg2::Integer,Float64::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPMonitorSolution,petsc1),PetscErrorCode,(KSP{Float64},Int32,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorSAWs(arg1::KSP{Float64},arg2::Integer,Float64::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPMonitorSAWs,petsc1),PetscErrorCode,(KSP{Float64},Int32,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorSAWsCreate(arg1::KSP{Float64},arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:KSPMonitorSAWsCreate,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Ptr{Void}}),arg1,arg2)
end

function KSPMonitorSAWsDestroy(arg0::Type{Float64},arg1::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:KSPMonitorSAWsDestroy,petsc1),PetscErrorCode,(Ptr{Ptr{Void}},),arg1)
end

function KSPGMRESMonitorKrylov(arg1::KSP{Float64},arg2::Integer,Float64::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPGMRESMonitorKrylov,petsc1),PetscErrorCode,(KSP{Float64},Int32,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPUnwindPreconditioner(arg1::KSP{Float64},arg2::Vec{Float64},arg3::Vec{Float64})
    ccall((:KSPUnwindPreconditioner,petsc1),PetscErrorCode,(KSP{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3)
end

function KSPInitialResidual(arg1::KSP{Float64},arg2::Vec{Float64},arg3::Vec{Float64},arg4::Vec{Float64},arg5::Vec{Float64},arg6::Vec{Float64})
    ccall((:KSPInitialResidual,petsc1),PetscErrorCode,(KSP{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64},Vec{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function KSPSetOperators(arg1::KSP{Float64},arg2::Mat{Float64},arg3::Mat{Float64})
    ccall((:KSPSetOperators,petsc1),PetscErrorCode,(KSP{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3)
end

function KSPGetOperators(arg1::KSP{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:KSPGetOperators,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Mat{Float64}},Ptr{Mat{Float64}}),arg1,arg2,arg3)
end

function KSPGetOperatorsSet(arg1::KSP{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}),arg3::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:KSPGetOperatorsSet,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3)
end

function KSPSetOptionsPrefix(arg1::KSP{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:KSPSetOptionsPrefix,petsc1),PetscErrorCode,(KSP{Float64},Cstring),arg1,arg2)
end

function KSPAppendOptionsPrefix(arg1::KSP{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:KSPAppendOptionsPrefix,petsc1),PetscErrorCode,(KSP{Float64},Cstring),arg1,arg2)
end

function KSPGetOptionsPrefix(arg1::KSP{Float64},arg2::Union(Ptr{ASCIIString},StridedArray{ASCIIString},Ptr{Void}))
    ccall((:KSPGetOptionsPrefix,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Ptr{UInt8}}),arg1,arg2)
end

function KSPSetTabLevel(arg1::KSP{Float64},arg2::Integer)
    ccall((:KSPSetTabLevel,petsc1),PetscErrorCode,(KSP{Float64},Int32),arg1,arg2)
end

function KSPGetTabLevel(arg1::KSP{Float64},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:KSPGetTabLevel,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Int32}),arg1,arg2)
end

function KSPSetDiagonalScale(arg1::KSP{Float64},arg2::PetscBool)
    ccall((:KSPSetDiagonalScale,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
end

function KSPGetDiagonalScale(arg1::KSP{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:KSPGetDiagonalScale,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PetscBool}),arg1,arg2)
end

function KSPSetDiagonalScaleFix(arg1::KSP{Float64},arg2::PetscBool)
    ccall((:KSPSetDiagonalScaleFix,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
end

function KSPGetDiagonalScaleFix(arg1::KSP{Float64},arg2::Union(Ptr{PetscBool},StridedArray{PetscBool},Ptr{Void}))
    ccall((:KSPGetDiagonalScaleFix,petsc1),PetscErrorCode,(KSP{Float64},Ptr{PetscBool}),arg1,arg2)
end

function KSPView(arg1::KSP{Float64},arg2::PetscViewer{Float64})
    ccall((:KSPView,petsc1),PetscErrorCode,(KSP{Float64},PetscViewer{Float64}),arg1,arg2)
end

function KSPLoad(arg1::KSP{Float64},arg2::PetscViewer{Float64})
    ccall((:KSPLoad,petsc1),PetscErrorCode,(KSP{Float64},PetscViewer{Float64}),arg1,arg2)
end

function KSPReasonViewFromOptions(arg1::KSP{Float64})
    ccall((:KSPReasonViewFromOptions,petsc1),PetscErrorCode,(KSP{Float64},),arg1)
end

function KSPLSQRSetStandardErrorVec(arg1::KSP{Float64},arg2::Vec{Float64})
    ccall((:KSPLSQRSetStandardErrorVec,petsc1),PetscErrorCode,(KSP{Float64},Vec{Float64}),arg1,arg2)
end

function KSPLSQRGetStandardErrorVec(arg1::KSP{Float64},arg2::Union(Ptr{Vec{Float64}},StridedArray{Vec{Float64}},Ptr{Void}))
    ccall((:KSPLSQRGetStandardErrorVec,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Vec{Float64}}),arg1,arg2)
end

function PCRedundantGetKSP(arg1::PC{Float64},arg2::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    ccall((:PCRedundantGetKSP,petsc1),PetscErrorCode,(PC{Float64},Ptr{KSP{Float64}}),arg1,arg2)
end

function PCRedistributeGetKSP(arg1::PC{Float64},arg2::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    ccall((:PCRedistributeGetKSP,petsc1),PetscErrorCode,(PC{Float64},Ptr{KSP{Float64}}),arg1,arg2)
end

function KSPSetNormType(arg1::KSP{Float64},arg2::KSPNormType)
    ccall((:KSPSetNormType,petsc1),PetscErrorCode,(KSP{Float64},KSPNormType),arg1,arg2)
end

function KSPGetNormType(arg1::KSP{Float64},arg2::Union(Ptr{KSPNormType},StridedArray{KSPNormType},Ptr{Void}))
    ccall((:KSPGetNormType,petsc1),PetscErrorCode,(KSP{Float64},Ptr{KSPNormType}),arg1,arg2)
end

function KSPSetSupportedNorm(ksp::KSP{Float64},arg1::KSPNormType,arg2::PCSide,arg3::Integer)
    ccall((:KSPSetSupportedNorm,petsc1),PetscErrorCode,(KSP{Float64},KSPNormType,PCSide,Int32),ksp,arg1,arg2,arg3)
end

function KSPSetCheckNormIteration(arg1::KSP{Float64},arg2::Integer)
    ccall((:KSPSetCheckNormIteration,petsc1),PetscErrorCode,(KSP{Float64},Int32),arg1,arg2)
end

function KSPSetLagNorm(arg1::KSP{Float64},arg2::PetscBool)
    ccall((:KSPSetLagNorm,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
end

function KSPSetConvergenceTest(arg1::KSP{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPSetConvergenceTest,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function KSPGetConvergenceContext(arg1::KSP{Float64},arg2::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:KSPGetConvergenceContext,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Ptr{Void}}),arg1,arg2)
end

function KSPConvergedDefault(arg1::KSP{Float64},arg2::Integer,Float64::Integer,arg3::Union(Ptr{KSPConvergedReason},StridedArray{KSPConvergedReason},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPConvergedDefault,petsc1),PetscErrorCode,(KSP{Float64},Int32,Cint,Ptr{KSPConvergedReason},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end

function KSPConvergedLSQR(arg1::KSP{Float64},arg2::Integer,Float64::Integer,arg3::Union(Ptr{KSPConvergedReason},StridedArray{KSPConvergedReason},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPConvergedLSQR,petsc1),PetscErrorCode,(KSP{Float64},Int32,Cint,Ptr{KSPConvergedReason},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end

function KSPConvergedDefaultDestroy(arg0::Type{Float64},arg1::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPConvergedDefaultDestroy,petsc1),PetscErrorCode,(Ptr{Void},),arg1)
end

function KSPConvergedDefaultCreate(arg0::Type{Float64},arg1::Union(Ptr{Ptr{Void}},StridedArray{Ptr{Void}},Ptr{Void}))
    ccall((:KSPConvergedDefaultCreate,petsc1),PetscErrorCode,(Ptr{Ptr{Void}},),arg1)
end

function KSPConvergedDefaultSetUIRNorm(arg1::KSP{Float64})
    ccall((:KSPConvergedDefaultSetUIRNorm,petsc1),PetscErrorCode,(KSP{Float64},),arg1)
end

function KSPConvergedDefaultSetUMIRNorm(arg1::KSP{Float64})
    ccall((:KSPConvergedDefaultSetUMIRNorm,petsc1),PetscErrorCode,(KSP{Float64},),arg1)
end

function KSPConvergedSkip(arg1::KSP{Float64},arg2::Integer,Float64::Integer,arg3::Union(Ptr{KSPConvergedReason},StridedArray{KSPConvergedReason},Ptr{Void}),arg4::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPConvergedSkip,petsc1),PetscErrorCode,(KSP{Float64},Int32,Cint,Ptr{KSPConvergedReason},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end

function KSPGetConvergedReason(arg1::KSP{Float64},arg2::Union(Ptr{KSPConvergedReason},StridedArray{KSPConvergedReason},Ptr{Void}))
    ccall((:KSPGetConvergedReason,petsc1),PetscErrorCode,(KSP{Float64},Ptr{KSPConvergedReason}),arg1,arg2)
end

function KSPCGSetType(arg1::KSP{Float64},arg2::KSPCGType)
    ccall((:KSPCGSetType,petsc1),PetscErrorCode,(KSP{Float64},KSPCGType),arg1,arg2)
end

function KSPCGUseSingleReduction(arg1::KSP{Float64},arg2::PetscBool)
    ccall((:KSPCGUseSingleReduction,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
end

function KSPNASHSetRadius(arg1::KSP{Float64},Float64::Integer)
    ccall((:KSPNASHSetRadius,petsc1),PetscErrorCode,(KSP{Float64},Cint),arg1,PetscReal)
end

function KSPNASHGetNormD(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:KSPNASHGetNormD,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint}),arg1,arg2)
end

function KSPNASHGetObjFcn(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:KSPNASHGetObjFcn,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint}),arg1,arg2)
end

function KSPSTCGSetRadius(arg1::KSP{Float64},Float64::Integer)
    ccall((:KSPSTCGSetRadius,petsc1),PetscErrorCode,(KSP{Float64},Cint),arg1,PetscReal)
end

function KSPSTCGGetNormD(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:KSPSTCGGetNormD,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint}),arg1,arg2)
end

function KSPSTCGGetObjFcn(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:KSPSTCGGetObjFcn,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint}),arg1,arg2)
end

function KSPGLTRSetRadius(arg1::KSP{Float64},Float64::Integer)
    ccall((:KSPGLTRSetRadius,petsc1),PetscErrorCode,(KSP{Float64},Cint),arg1,PetscReal)
end

function KSPGLTRGetNormD(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:KSPGLTRGetNormD,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint}),arg1,arg2)
end

function KSPGLTRGetObjFcn(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:KSPGLTRGetObjFcn,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint}),arg1,arg2)
end

function KSPGLTRGetMinEig(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:KSPGLTRGetMinEig,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint}),arg1,arg2)
end

function KSPGLTRGetLambda(arg1::KSP{Float64},arg2::Union(Ptr{Cint},StridedArray{Cint},Ptr{Void}))
    ccall((:KSPGLTRGetLambda,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Cint}),arg1,arg2)
end

function KSPPythonSetType(arg1::KSP{Float64},arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}))
    ccall((:KSPPythonSetType,petsc1),PetscErrorCode,(KSP{Float64},Cstring),arg1,arg2)
end

function PCPreSolve(arg1::PC{Float64},arg2::KSP{Float64})
    ccall((:PCPreSolve,petsc1),PetscErrorCode,(PC{Float64},KSP{Float64}),arg1,arg2)
end

function PCPostSolve(arg1::PC{Float64},arg2::KSP{Float64})
    ccall((:PCPostSolve,petsc1),PetscErrorCode,(PC{Float64},KSP{Float64}),arg1,arg2)
end

#= skipping function with undefined symbols: 
 function KSPMonitorLGResidualNormCreate(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Union(Ptr{Ptr{PetscObject}},StridedArray{Ptr{PetscObject}},Ptr{Void}))
    ccall((:KSPMonitorLGResidualNormCreate,petsc1),PetscErrorCode,(Cstring,Cstring,Cint,Cint,Cint,Cint,Ptr{Ptr{PetscObject}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function KSPMonitorLGResidualNorm(arg1::KSP{Float64},arg2::Integer,Float64::Integer,arg3::Union(Ptr{PetscObject},StridedArray{PetscObject},Ptr{Void}))
    ccall((:KSPMonitorLGResidualNorm,petsc1),PetscErrorCode,(KSP{Float64},Int32,Cint,Ptr{PetscObject}),arg1,arg2,PetscReal,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function KSPMonitorLGResidualNormDestroy(arg0::Type{Float64},arg1::Union(Ptr{Ptr{PetscObject}},StridedArray{Ptr{PetscObject}},Ptr{Void}))
    ccall((:KSPMonitorLGResidualNormDestroy,petsc1),PetscErrorCode,(Ptr{Ptr{PetscObject}},),arg1)
end 
=#
#= skipping function with undefined symbols: 
 function KSPMonitorLGTrueResidualNormCreate(arg1::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg2::Union(ASCIIString,StridedArray{Uint8},Ptr{Void}),arg3::Integer,arg4::Integer,arg5::Integer,arg6::Integer,arg7::Union(Ptr{Ptr{PetscObject}},StridedArray{Ptr{PetscObject}},Ptr{Void}))
    ccall((:KSPMonitorLGTrueResidualNormCreate,petsc1),PetscErrorCode,(Cstring,Cstring,Cint,Cint,Cint,Cint,Ptr{Ptr{PetscObject}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end 
=#
#= skipping function with undefined symbols: 
 function KSPMonitorLGTrueResidualNorm(arg1::KSP{Float64},arg2::Integer,Float64::Integer,arg3::Union(Ptr{PetscObject},StridedArray{PetscObject},Ptr{Void}))
    ccall((:KSPMonitorLGTrueResidualNorm,petsc1),PetscErrorCode,(KSP{Float64},Int32,Cint,Ptr{PetscObject}),arg1,arg2,PetscReal,arg3)
end 
=#
#= skipping function with undefined symbols: 
 function KSPMonitorLGTrueResidualNormDestroy(arg0::Type{Float64},arg1::Union(Ptr{Ptr{PetscObject}},StridedArray{Ptr{PetscObject}},Ptr{Void}))
    ccall((:KSPMonitorLGTrueResidualNormDestroy,petsc1),PetscErrorCode,(Ptr{Ptr{PetscObject}},),arg1)
end 
=#
function KSPMonitorLGRange(arg1::KSP{Float64},arg2::Integer,Float64::Integer,arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPMonitorLGRange,petsc1),PetscErrorCode,(KSP{Float64},Int32,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function PCShellSetPreSolve(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PCShellSetPreSolve,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
end

function PCShellSetPostSolve(arg1::PC{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:PCShellSetPostSolve,petsc1),PetscErrorCode,(PC{Float64},Ptr{Void}),arg1,arg2)
end

#= skipping function with undefined symbols: 
 function KSPFischerGuessCreate(arg1::KSP{Float64},arg2::Integer,arg3::Integer,arg4::Union(Ptr{KSPFischerGuess},StridedArray{KSPFischerGuess},Ptr{Void}))
    ccall((:KSPFischerGuessCreate,petsc1),PetscErrorCode,(KSP{Float64},Int32,Int32,Ptr{KSPFischerGuess}),arg1,arg2,arg3,arg4)
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
    ccall((:KSPSetUseFischerGuess,petsc1),PetscErrorCode,(KSP{Float64},Int32,Int32),arg1,arg2,arg3)
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
    ccall((:MatCreateSchurComplement,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64},Mat{Float64},Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatSchurComplementGetKSP(arg1::Mat{Float64},arg2::Union(Ptr{KSP{Float64}},StridedArray{KSP{Float64}},Ptr{Void}))
    ccall((:MatSchurComplementGetKSP,petsc1),PetscErrorCode,(Mat{Float64},Ptr{KSP{Float64}}),arg1,arg2)
end

function MatSchurComplementSetKSP(arg1::Mat{Float64},arg2::KSP{Float64})
    ccall((:MatSchurComplementSetKSP,petsc1),PetscErrorCode,(Mat{Float64},KSP{Float64}),arg1,arg2)
end

function MatSchurComplementSetSubMatrices(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64},arg4::Mat{Float64},arg5::Mat{Float64},arg6::Mat{Float64})
    ccall((:MatSchurComplementSetSubMatrices,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64},Mat{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatSchurComplementUpdateSubMatrices(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64},arg4::Mat{Float64},arg5::Mat{Float64},arg6::Mat{Float64})
    ccall((:MatSchurComplementUpdateSubMatrices,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64},Mat{Float64},Mat{Float64},Mat{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatSchurComplementGetSubMatrices(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg4::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg5::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg6::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatSchurComplementGetSubMatrices,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}},Ptr{Mat{Float64}},Ptr{Mat{Float64}},Ptr{Mat{Float64}},Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatSchurComplementSetAinvType(arg1::Mat{Float64},arg2::MatSchurComplementAinvType)
    ccall((:MatSchurComplementSetAinvType,petsc1),PetscErrorCode,(Mat{Float64},MatSchurComplementAinvType),arg1,arg2)
end

function MatSchurComplementGetAinvType(arg1::Mat{Float64},arg2::Union(Ptr{MatSchurComplementAinvType},StridedArray{MatSchurComplementAinvType},Ptr{Void}))
    ccall((:MatSchurComplementGetAinvType,petsc1),PetscErrorCode,(Mat{Float64},Ptr{MatSchurComplementAinvType}),arg1,arg2)
end

function MatSchurComplementGetPmat(arg1::Mat{Float64},arg2::MatReuse,arg3::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatSchurComplementGetPmat,petsc1),PetscErrorCode,(Mat{Float64},MatReuse,Ptr{Mat{Float64}}),arg1,arg2,arg3)
end

function MatSchurComplementComputeExplicitOperator(arg1::Mat{Float64},arg2::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatSchurComplementComputeExplicitOperator,petsc1),PetscErrorCode,(Mat{Float64},Ptr{Mat{Float64}}),arg1,arg2)
end

function MatGetSchurComplement(arg1::Mat{Float64},arg2::IS{Float64},arg3::IS{Float64},arg4::IS{Float64},arg5::IS{Float64},arg6::MatReuse,arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}),arg8::MatSchurComplementAinvType,arg9::MatReuse,arg10::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatGetSchurComplement,petsc1),PetscErrorCode,(Mat{Float64},IS{Float64},IS{Float64},IS{Float64},IS{Float64},MatReuse,Ptr{Mat{Float64}},MatSchurComplementAinvType,MatReuse,Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function MatCreateSchurComplementPmat(arg1::Mat{Float64},arg2::Mat{Float64},arg3::Mat{Float64},arg4::Mat{Float64},arg5::MatSchurComplementAinvType,arg6::MatReuse,arg7::Union(Ptr{Mat{Float64}},StridedArray{Mat{Float64}},Ptr{Void}))
    ccall((:MatCreateSchurComplementPmat,petsc1),PetscErrorCode,(Mat{Float64},Mat{Float64},Mat{Float64},Mat{Float64},MatSchurComplementAinvType,MatReuse,Ptr{Mat{Float64}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

#= skipping function with undefined symbols: 
 function KSPSetDM(arg1::KSP{Float64},arg2::DM)
    ccall((:KSPSetDM,petsc1),PetscErrorCode,(KSP{Float64},DM),arg1,arg2)
end 
=#
function KSPSetDMActive(arg1::KSP{Float64},arg2::PetscBool)
    ccall((:KSPSetDMActive,petsc1),PetscErrorCode,(KSP{Float64},PetscBool),arg1,arg2)
end

#= skipping function with undefined symbols: 
 function KSPGetDM(arg1::KSP{Float64},arg2::Union(Ptr{DM},StridedArray{DM},Ptr{Void}))
    ccall((:KSPGetDM,petsc1),PetscErrorCode,(KSP{Float64},Ptr{DM}),arg1,arg2)
end 
=#
function KSPSetApplicationContext(arg1::KSP{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPSetApplicationContext,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Void}),arg1,arg2)
end

function KSPGetApplicationContext(arg1::KSP{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPGetApplicationContext,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Void}),arg1,arg2)
end

function KSPSetComputeRHS(arg1::KSP{Float64},func::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPSetComputeRHS,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Void},Ptr{Void}),arg1,func,arg2)
end

function KSPSetComputeOperators(arg1::KSP{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPSetComputeOperators,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function KSPSetComputeInitialGuess(arg1::KSP{Float64},arg2::Union(Ptr{Void},StridedArray{Void},Ptr{Void}),arg3::Union(Ptr{Void},StridedArray{Void},Ptr{Void}))
    ccall((:KSPSetComputeInitialGuess,petsc1),PetscErrorCode,(KSP{Float64},Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
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
