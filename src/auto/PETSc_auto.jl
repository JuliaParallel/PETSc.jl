# Julia wrapper for header: /home/jared/build/petsc-3.6.0/include/petsc.h
# Automatically generated using Clang.jl wrap_c, version 0.0.0

module PETSc

include("libPETSc_common.jl")
function PetscIsInfOrNanReal()
    ccall((:PetscIsInfOrNanReal,petsc),PetscErrorCode,())
end

function PetscIsNormalReal()
    ccall((:PetscIsNormalReal,petsc),PetscBool,())
end

function PetscSetHelpVersionFunctions(arg1::Ptr{Void},arg2::Ptr{Void})
    ccall((:PetscSetHelpVersionFunctions,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void}),arg1,arg2)
end

function PetscCommDuplicate(arg1::MPI_Comm,arg2::Ptr{MPI_Comm},arg3::Ptr{Cint})
    ccall((:PetscCommDuplicate,petsc),PetscErrorCode,(MPI_Comm,Ptr{MPI_Comm},Ptr{Cint}),arg1,arg2,arg3)
end

function PetscCommDestroy(arg1::Ptr{MPI_Comm})
    ccall((:PetscCommDestroy,petsc),PetscErrorCode,(Ptr{MPI_Comm},),arg1)
end

function PetscMallocSet(arg1::Ptr{Void},arg2::Ptr{Void})
    ccall((:PetscMallocSet,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void}),arg1,arg2)
end

function PetscMallocClear()
    ccall((:PetscMallocClear,petsc),PetscErrorCode,())
end

function PetscMallocDump(arg1::Ptr{Void})
    ccall((:PetscMallocDump,petsc),PetscErrorCode,(Ptr{Void},),arg1)
end

function PetscMallocDumpLog(arg1::Ptr{Void})
    ccall((:PetscMallocDumpLog,petsc),PetscErrorCode,(Ptr{Void},),arg1)
end

function PetscMallocGetCurrentUsage(arg1::Ptr{PetscLogDouble})
    ccall((:PetscMallocGetCurrentUsage,petsc),PetscErrorCode,(Ptr{PetscLogDouble},),arg1)
end

function PetscMallocGetMaximumUsage(arg1::Ptr{PetscLogDouble})
    ccall((:PetscMallocGetMaximumUsage,petsc),PetscErrorCode,(Ptr{PetscLogDouble},),arg1)
end

function PetscMallocDebug(arg1::PetscBool)
    ccall((:PetscMallocDebug,petsc),PetscErrorCode,(PetscBool,),arg1)
end

function PetscMallocGetDebug(arg1::Ptr{PetscBool})
    ccall((:PetscMallocGetDebug,petsc),PetscErrorCode,(Ptr{PetscBool},),arg1)
end

function PetscMallocValidate(arg1::Cint,arg2::Ptr{Uint8},arg3::Ptr{Uint8})
    ccall((:PetscMallocValidate,petsc),PetscErrorCode,(Cint,Ptr{Uint8},Ptr{Uint8}),arg1,arg2,arg3)
end

function PetscMallocSetDumpLog()
    ccall((:PetscMallocSetDumpLog,petsc),PetscErrorCode,())
end

function PetscMallocSetDumpLogThreshold(arg1::PetscLogDouble)
    ccall((:PetscMallocSetDumpLogThreshold,petsc),PetscErrorCode,(PetscLogDouble,),arg1)
end

function PetscMallocGetDumpLog(arg1::Ptr{PetscBool})
    ccall((:PetscMallocGetDumpLog,petsc),PetscErrorCode,(Ptr{PetscBool},),arg1)
end

typealias MPI_Datatype MPI.MPIDatatype

function PetscDataTypeToMPIDataType(arg1::PetscDataType,arg2::Ptr{MPI_Datatype})
    ccall((:PetscDataTypeToMPIDataType,petsc),PetscErrorCode,(PetscDataType,Ptr{MPI_Datatype}),arg1,arg2)
end

function PetscMPIDataTypeToPetscDataType(arg1::MPI_Datatype,arg2::Ptr{PetscDataType})
    ccall((:PetscMPIDataTypeToPetscDataType,petsc),PetscErrorCode,(MPI_Datatype,Ptr{PetscDataType}),arg1,arg2)
end

function PetscDataTypeGetSize(arg1::PetscDataType,arg2::Ptr{Cint})
    ccall((:PetscDataTypeGetSize,petsc),PetscErrorCode,(PetscDataType,Ptr{Cint}),arg1,arg2)
end

function PetscDataTypeFromString(arg1::Ptr{Uint8},arg2::Ptr{PetscDataType},arg3::Ptr{PetscBool})
    ccall((:PetscDataTypeFromString,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{PetscDataType},Ptr{PetscBool}),arg1,arg2,arg3)
end

function PetscBitMemcpy(arg1::Ptr{Void},arg2::PetscInt,arg3::Ptr{Void},arg4::PetscInt,arg5::PetscInt,arg6::PetscDataType)
    ccall((:PetscBitMemcpy,petsc),PetscErrorCode,(Ptr{Void},PetscInt,Ptr{Void},PetscInt,PetscInt,PetscDataType),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscMemmove(arg1::Ptr{Void},arg2::Ptr{Void},size_t::Cint)
    ccall((:PetscMemmove,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Cint),arg1,arg2,size_t)
end

function PetscMemcmp(arg1::Ptr{Void},arg2::Ptr{Void},size_t::Cint,arg3::Ptr{PetscBool})
    ccall((:PetscMemcmp,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Cint,Ptr{PetscBool}),arg1,arg2,size_t,arg3)
end

function PetscStrlen(arg1::Ptr{Uint8},arg2::Ptr{Cint})
    ccall((:PetscStrlen,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Cint}),arg1,arg2)
end

function PetscStrToArray(arg1::Ptr{Uint8},arg2::Uint8,arg3::Ptr{Cint},arg4::Ptr{Ptr{Ptr{Uint8}}})
    ccall((:PetscStrToArray,petsc),PetscErrorCode,(Ptr{Uint8},Uint8,Ptr{Cint},Ptr{Ptr{Ptr{Uint8}}}),arg1,arg2,arg3,arg4)
end

function PetscStrToArrayDestroy(arg1::Cint,arg2::Ptr{Ptr{Uint8}})
    ccall((:PetscStrToArrayDestroy,petsc),PetscErrorCode,(Cint,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function PetscStrcmp(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{PetscBool})
    ccall((:PetscStrcmp,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{PetscBool}),arg1,arg2,arg3)
end

function PetscStrgrt(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{PetscBool})
    ccall((:PetscStrgrt,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{PetscBool}),arg1,arg2,arg3)
end

function PetscStrcasecmp(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{PetscBool})
    ccall((:PetscStrcasecmp,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{PetscBool}),arg1,arg2,arg3)
end

function PetscStrncmp(arg1::Ptr{Uint8},arg2::Ptr{Uint8},size_t::Cint,arg3::Ptr{PetscBool})
    ccall((:PetscStrncmp,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Cint,Ptr{PetscBool}),arg1,arg2,size_t,arg3)
end

function PetscStrcpy(arg1::Ptr{Uint8},arg2::Ptr{Uint8})
    ccall((:PetscStrcpy,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8}),arg1,arg2)
end

function PetscStrcat(arg1::Ptr{Uint8},arg2::Ptr{Uint8})
    ccall((:PetscStrcat,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8}),arg1,arg2)
end

function PetscStrncat(arg1::Ptr{Uint8},arg2::Ptr{Uint8},size_t::Cint)
    ccall((:PetscStrncat,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Cint),arg1,arg2,size_t)
end

function PetscStrncpy(arg1::Ptr{Uint8},arg2::Ptr{Uint8},size_t::Cint)
    ccall((:PetscStrncpy,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Cint),arg1,arg2,size_t)
end

function PetscStrchr(arg1::Ptr{Uint8},arg2::Uint8,arg3::Ptr{Ptr{Uint8}})
    ccall((:PetscStrchr,petsc),PetscErrorCode,(Ptr{Uint8},Uint8,Ptr{Ptr{Uint8}}),arg1,arg2,arg3)
end

function PetscStrtolower(arg1::Ptr{Uint8})
    ccall((:PetscStrtolower,petsc),PetscErrorCode,(Ptr{Uint8},),arg1)
end

function PetscStrtoupper(arg1::Ptr{Uint8})
    ccall((:PetscStrtoupper,petsc),PetscErrorCode,(Ptr{Uint8},),arg1)
end

function PetscStrrchr(arg1::Ptr{Uint8},arg2::Uint8,arg3::Ptr{Ptr{Uint8}})
    ccall((:PetscStrrchr,petsc),PetscErrorCode,(Ptr{Uint8},Uint8,Ptr{Ptr{Uint8}}),arg1,arg2,arg3)
end

function PetscStrstr(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{Ptr{Uint8}})
    ccall((:PetscStrstr,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{Ptr{Uint8}}),arg1,arg2,arg3)
end

function PetscStrrstr(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{Ptr{Uint8}})
    ccall((:PetscStrrstr,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{Ptr{Uint8}}),arg1,arg2,arg3)
end

function PetscStrendswith(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{PetscBool})
    ccall((:PetscStrendswith,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{PetscBool}),arg1,arg2,arg3)
end

function PetscStrbeginswith(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{PetscBool})
    ccall((:PetscStrbeginswith,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{PetscBool}),arg1,arg2,arg3)
end

function PetscStrendswithwhich(arg1::Ptr{Uint8},arg2::Ptr{Ptr{Uint8}},arg3::Ptr{PetscInt})
    ccall((:PetscStrendswithwhich,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Ptr{Uint8}},Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscStrallocpy(arg1::Ptr{Uint8},arg2::Ptr{Ptr{Uint8}})
    ccall((:PetscStrallocpy,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Ptr{Uint8}}),arg1,arg2)
end

function PetscStrArrayallocpy(arg1::Ptr{Ptr{Uint8}},arg2::Ptr{Ptr{Ptr{Uint8}}})
    ccall((:PetscStrArrayallocpy,petsc),PetscErrorCode,(Ptr{Ptr{Uint8}},Ptr{Ptr{Ptr{Uint8}}}),arg1,arg2)
end

function PetscStrArrayDestroy(arg1::Ptr{Ptr{Ptr{Uint8}}})
    ccall((:PetscStrArrayDestroy,petsc),PetscErrorCode,(Ptr{Ptr{Ptr{Uint8}}},),arg1)
end

function PetscStrNArrayallocpy(arg1::PetscInt,arg2::Ptr{Ptr{Uint8}},arg3::Ptr{Ptr{Ptr{Uint8}}})
    ccall((:PetscStrNArrayallocpy,petsc),PetscErrorCode,(PetscInt,Ptr{Ptr{Uint8}},Ptr{Ptr{Ptr{Uint8}}}),arg1,arg2,arg3)
end

function PetscStrNArrayDestroy(arg1::PetscInt,arg2::Ptr{Ptr{Ptr{Uint8}}})
    ccall((:PetscStrNArrayDestroy,petsc),PetscErrorCode,(PetscInt,Ptr{Ptr{Ptr{Uint8}}}),arg1,arg2)
end

function PetscStrreplace(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},size_t::Cint)
    ccall((:PetscStrreplace,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Cint),arg1,arg2,arg3,size_t)
end

function PetscStrcmpNoError(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{PetscBool})
    ccall((:PetscStrcmpNoError,petsc),Void,(Ptr{Uint8},Ptr{Uint8},Ptr{PetscBool}),arg1,arg2,arg3)
end

function PetscTokenCreate(arg1::Ptr{Uint8},arg2::Uint8,arg3::Ptr{PetscToken})
    ccall((:PetscTokenCreate,petsc),PetscErrorCode,(Ptr{Uint8},Uint8,Ptr{PetscToken}),arg1,arg2,arg3)
end

function PetscTokenFind(arg1::PetscToken,arg2::Ptr{Ptr{Uint8}})
    ccall((:PetscTokenFind,petsc),PetscErrorCode,(PetscToken,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function PetscTokenDestroy(arg1::Ptr{PetscToken})
    ccall((:PetscTokenDestroy,petsc),PetscErrorCode,(Ptr{PetscToken},),arg1)
end

function PetscEListFind(arg1::PetscInt,arg2::Ptr{Ptr{Uint8}},arg3::Ptr{Uint8},arg4::Ptr{PetscInt},arg5::Ptr{PetscBool})
    ccall((:PetscEListFind,petsc),PetscErrorCode,(PetscInt,Ptr{Ptr{Uint8}},Ptr{Uint8},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function PetscEnumFind(arg1::Ptr{Ptr{Uint8}},arg2::Ptr{Uint8},arg3::Ptr{PetscEnum},arg4::Ptr{PetscBool})
    ccall((:PetscEnumFind,petsc),PetscErrorCode,(Ptr{Ptr{Uint8}},Ptr{Uint8},Ptr{PetscEnum},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function PetscMaxSum(arg1::MPI_Comm,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:PetscMaxSum,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function MPIULong_Send(arg1::Ptr{Void},arg2::PetscInt,arg3::MPI_Datatype,arg4::PetscMPIInt,arg5::PetscMPIInt,arg6::MPI_Comm)
    ccall((:MPIULong_Send,petsc),PetscErrorCode,(Ptr{Void},PetscInt,MPI_Datatype,PetscMPIInt,PetscMPIInt,MPI_Comm),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MPIULong_Recv(arg1::Ptr{Void},arg2::PetscInt,arg3::MPI_Datatype,arg4::PetscMPIInt,arg5::PetscMPIInt,arg6::MPI_Comm)
    ccall((:MPIULong_Recv,petsc),PetscErrorCode,(Ptr{Void},PetscInt,MPI_Datatype,PetscMPIInt,PetscMPIInt,MPI_Comm),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscErrorPrintfInitialize()
    ccall((:PetscErrorPrintfInitialize,petsc),PetscErrorCode,())
end

function PetscErrorMessage(arg1::Cint,arg2::Ptr{Ptr{Uint8}},arg3::Ptr{Ptr{Uint8}})
    ccall((:PetscErrorMessage,petsc),PetscErrorCode,(Cint,Ptr{Ptr{Uint8}},Ptr{Ptr{Uint8}}),arg1,arg2,arg3)
end

function PetscTraceBackErrorHandler(arg1::MPI_Comm,arg2::Cint,arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::PetscErrorCode,arg6::PetscErrorType,arg7::Ptr{Uint8},arg8::Ptr{Void})
    ccall((:PetscTraceBackErrorHandler,petsc),PetscErrorCode,(MPI_Comm,Cint,Ptr{Uint8},Ptr{Uint8},PetscErrorCode,PetscErrorType,Ptr{Uint8},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscIgnoreErrorHandler(arg1::MPI_Comm,arg2::Cint,arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::PetscErrorCode,arg6::PetscErrorType,arg7::Ptr{Uint8},arg8::Ptr{Void})
    ccall((:PetscIgnoreErrorHandler,petsc),PetscErrorCode,(MPI_Comm,Cint,Ptr{Uint8},Ptr{Uint8},PetscErrorCode,PetscErrorType,Ptr{Uint8},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscEmacsClientErrorHandler(arg1::MPI_Comm,arg2::Cint,arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::PetscErrorCode,arg6::PetscErrorType,arg7::Ptr{Uint8},arg8::Ptr{Void})
    ccall((:PetscEmacsClientErrorHandler,petsc),PetscErrorCode,(MPI_Comm,Cint,Ptr{Uint8},Ptr{Uint8},PetscErrorCode,PetscErrorType,Ptr{Uint8},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscMPIAbortErrorHandler(arg1::MPI_Comm,arg2::Cint,arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::PetscErrorCode,arg6::PetscErrorType,arg7::Ptr{Uint8},arg8::Ptr{Void})
    ccall((:PetscMPIAbortErrorHandler,petsc),PetscErrorCode,(MPI_Comm,Cint,Ptr{Uint8},Ptr{Uint8},PetscErrorCode,PetscErrorType,Ptr{Uint8},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscAbortErrorHandler(arg1::MPI_Comm,arg2::Cint,arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::PetscErrorCode,arg6::PetscErrorType,arg7::Ptr{Uint8},arg8::Ptr{Void})
    ccall((:PetscAbortErrorHandler,petsc),PetscErrorCode,(MPI_Comm,Cint,Ptr{Uint8},Ptr{Uint8},PetscErrorCode,PetscErrorType,Ptr{Uint8},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscAttachDebuggerErrorHandler(arg1::MPI_Comm,arg2::Cint,arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::PetscErrorCode,arg6::PetscErrorType,arg7::Ptr{Uint8},arg8::Ptr{Void})
    ccall((:PetscAttachDebuggerErrorHandler,petsc),PetscErrorCode,(MPI_Comm,Cint,Ptr{Uint8},Ptr{Uint8},PetscErrorCode,PetscErrorType,Ptr{Uint8},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscReturnErrorHandler(arg1::MPI_Comm,arg2::Cint,arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::PetscErrorCode,arg6::PetscErrorType,arg7::Ptr{Uint8},arg8::Ptr{Void})
    ccall((:PetscReturnErrorHandler,petsc),PetscErrorCode,(MPI_Comm,Cint,Ptr{Uint8},Ptr{Uint8},PetscErrorCode,PetscErrorType,Ptr{Uint8},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscPushErrorHandler(handler::Ptr{Void},arg1::Ptr{Void})
    ccall((:PetscPushErrorHandler,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void}),handler,arg1)
end

function PetscPopErrorHandler()
    ccall((:PetscPopErrorHandler,petsc),PetscErrorCode,())
end

function PetscSignalHandlerDefault(arg1::Cint,arg2::Ptr{Void})
    ccall((:PetscSignalHandlerDefault,petsc),PetscErrorCode,(Cint,Ptr{Void}),arg1,arg2)
end

function PetscPushSignalHandler(arg1::Ptr{Void},arg2::Ptr{Void})
    ccall((:PetscPushSignalHandler,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void}),arg1,arg2)
end

function PetscPopSignalHandler()
    ccall((:PetscPopSignalHandler,petsc),PetscErrorCode,())
end

function PetscCheckPointerSetIntensity(arg1::PetscInt)
    ccall((:PetscCheckPointerSetIntensity,petsc),PetscErrorCode,(PetscInt,),arg1)
end

function PetscSetFPTrap(arg1::PetscFPTrap)
    ccall((:PetscSetFPTrap,petsc),PetscErrorCode,(PetscFPTrap,),arg1)
end

function PetscFPTrapPush(arg1::PetscFPTrap)
    ccall((:PetscFPTrapPush,petsc),PetscErrorCode,(PetscFPTrap,),arg1)
end

function PetscFPTrapPop()
    ccall((:PetscFPTrapPop,petsc),PetscErrorCode,())
end

function PetscStackCopy(arg1::Ptr{PetscStack},arg2::Ptr{PetscStack})
    ccall((:PetscStackCopy,petsc),PetscErrorCode,(Ptr{PetscStack},Ptr{PetscStack}),arg1,arg2)
end

function PetscStackPrint(arg1::Ptr{PetscStack},arg2::Ptr{Void})
    ccall((:PetscStackPrint,petsc),PetscErrorCode,(Ptr{PetscStack},Ptr{Void}),arg1,arg2)
end

function PetscStackView(arg1::Ptr{Void})
    ccall((:PetscStackView,petsc),PetscErrorCode,(Ptr{Void},),arg1)
end

function PetscStackDestroy()
    ccall((:PetscStackDestroy,petsc),PetscErrorCode,())
end

function PetscClassIdRegister(arg1::Ptr{Uint8},arg2::Ptr{PetscClassId})
    ccall((:PetscClassIdRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{PetscClassId}),arg1,arg2)
end

function PetscMemoryGetCurrentUsage(arg1::Ptr{PetscLogDouble})
    ccall((:PetscMemoryGetCurrentUsage,petsc),PetscErrorCode,(Ptr{PetscLogDouble},),arg1)
end

function PetscMemoryGetMaximumUsage(arg1::Ptr{PetscLogDouble})
    ccall((:PetscMemoryGetMaximumUsage,petsc),PetscErrorCode,(Ptr{PetscLogDouble},),arg1)
end

function PetscMemorySetGetMaximumUsage()
    ccall((:PetscMemorySetGetMaximumUsage,petsc),PetscErrorCode,())
end

function PetscMemoryTrace(arg1::Ptr{Uint8})
    ccall((:PetscMemoryTrace,petsc),PetscErrorCode,(Ptr{Uint8},),arg1)
end

function PetscInfoAllow(arg1::PetscBool,arg2::Ptr{Uint8})
    ccall((:PetscInfoAllow,petsc),PetscErrorCode,(PetscBool,Ptr{Uint8}),arg1,arg2)
end

function PetscSleep()
    ccall((:PetscSleep,petsc),PetscErrorCode,())
end

function PetscInitialize(arg1::Ptr{Cint},arg2::Ptr{Ptr{Ptr{Uint8}}},arg3::Ptr{Uint8},arg4::Ptr{Uint8})
    ccall((:PetscInitialize,petsc),PetscErrorCode,(Ptr{Cint},Ptr{Ptr{Ptr{Uint8}}},Ptr{Uint8},Ptr{Uint8}),arg1,arg2,arg3,arg4)
end

function PetscInitializeNoPointers(arg1::Cint,arg2::Ptr{Ptr{Uint8}},arg3::Ptr{Uint8},arg4::Ptr{Uint8})
    ccall((:PetscInitializeNoPointers,petsc),PetscErrorCode,(Cint,Ptr{Ptr{Uint8}},Ptr{Uint8},Ptr{Uint8}),arg1,arg2,arg3,arg4)
end

function PetscInitializeNoArguments()
    ccall((:PetscInitializeNoArguments,petsc),PetscErrorCode,())
end

function PetscInitialized(arg1::Ptr{PetscBool})
    ccall((:PetscInitialized,petsc),PetscErrorCode,(Ptr{PetscBool},),arg1)
end

function PetscFinalized(arg1::Ptr{PetscBool})
    ccall((:PetscFinalized,petsc),PetscErrorCode,(Ptr{PetscBool},),arg1)
end

function PetscFinalize()
    ccall((:PetscFinalize,petsc),PetscErrorCode,())
end

function PetscInitializeFortran()
    ccall((:PetscInitializeFortran,petsc),PetscErrorCode,())
end

function PetscGetArgs(arg1::Ptr{Cint},arg2::Ptr{Ptr{Ptr{Uint8}}})
    ccall((:PetscGetArgs,petsc),PetscErrorCode,(Ptr{Cint},Ptr{Ptr{Ptr{Uint8}}}),arg1,arg2)
end

function PetscGetArguments(arg1::Ptr{Ptr{Ptr{Uint8}}})
    ccall((:PetscGetArguments,petsc),PetscErrorCode,(Ptr{Ptr{Ptr{Uint8}}},),arg1)
end

function PetscFreeArguments(arg1::Ptr{Ptr{Uint8}})
    ccall((:PetscFreeArguments,petsc),PetscErrorCode,(Ptr{Ptr{Uint8}},),arg1)
end

function PetscEnd()
    ccall((:PetscEnd,petsc),PetscErrorCode,())
end

function PetscSysInitializePackage()
    ccall((:PetscSysInitializePackage,petsc),PetscErrorCode,())
end

function PetscPythonInitialize(arg1::Ptr{Uint8},arg2::Ptr{Uint8})
    ccall((:PetscPythonInitialize,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8}),arg1,arg2)
end

function PetscPythonFinalize()
    ccall((:PetscPythonFinalize,petsc),PetscErrorCode,())
end

function PetscPythonPrintError()
    ccall((:PetscPythonPrintError,petsc),PetscErrorCode,())
end

function PetscPythonMonitorSet(arg1::PetscObject,arg2::Ptr{Uint8})
    ccall((:PetscPythonMonitorSet,petsc),PetscErrorCode,(PetscObject,Ptr{Uint8}),arg1,arg2)
end

function PetscObjectDestroy(arg1::Ptr{PetscObject})
    ccall((:PetscObjectDestroy,petsc),PetscErrorCode,(Ptr{PetscObject},),arg1)
end

function PetscObjectGetComm(arg1::PetscObject,arg2::Ptr{MPI_Comm})
    ccall((:PetscObjectGetComm,petsc),PetscErrorCode,(PetscObject,Ptr{MPI_Comm}),arg1,arg2)
end

function PetscObjectGetClassId(arg1::PetscObject,arg2::Ptr{PetscClassId})
    ccall((:PetscObjectGetClassId,petsc),PetscErrorCode,(PetscObject,Ptr{PetscClassId}),arg1,arg2)
end

function PetscObjectGetClassName(arg1::PetscObject,arg2::Ptr{Ptr{Uint8}})
    ccall((:PetscObjectGetClassName,petsc),PetscErrorCode,(PetscObject,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function PetscObjectSetType(arg1::PetscObject,arg2::Ptr{Uint8})
    ccall((:PetscObjectSetType,petsc),PetscErrorCode,(PetscObject,Ptr{Uint8}),arg1,arg2)
end

function PetscObjectSetPrecision(arg1::PetscObject,arg2::PetscPrecision)
    ccall((:PetscObjectSetPrecision,petsc),PetscErrorCode,(PetscObject,PetscPrecision),arg1,arg2)
end

function PetscObjectGetType(arg1::PetscObject,arg2::Ptr{Ptr{Uint8}})
    ccall((:PetscObjectGetType,petsc),PetscErrorCode,(PetscObject,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function PetscObjectSetName(arg1::PetscObject,arg2::Ptr{Uint8})
    ccall((:PetscObjectSetName,petsc),PetscErrorCode,(PetscObject,Ptr{Uint8}),arg1,arg2)
end

function PetscObjectGetName(arg1::PetscObject,arg2::Ptr{Ptr{Uint8}})
    ccall((:PetscObjectGetName,petsc),PetscErrorCode,(PetscObject,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function PetscObjectSetTabLevel(arg1::PetscObject,arg2::PetscInt)
    ccall((:PetscObjectSetTabLevel,petsc),PetscErrorCode,(PetscObject,PetscInt),arg1,arg2)
end

function PetscObjectGetTabLevel(arg1::PetscObject,arg2::Ptr{PetscInt})
    ccall((:PetscObjectGetTabLevel,petsc),PetscErrorCode,(PetscObject,Ptr{PetscInt}),arg1,arg2)
end

function PetscObjectIncrementTabLevel(arg1::PetscObject,arg2::PetscObject,arg3::PetscInt)
    ccall((:PetscObjectIncrementTabLevel,petsc),PetscErrorCode,(PetscObject,PetscObject,PetscInt),arg1,arg2,arg3)
end

function PetscObjectReference(arg1::PetscObject)
    ccall((:PetscObjectReference,petsc),PetscErrorCode,(PetscObject,),arg1)
end

function PetscObjectGetReference(arg1::PetscObject,arg2::Ptr{PetscInt})
    ccall((:PetscObjectGetReference,petsc),PetscErrorCode,(PetscObject,Ptr{PetscInt}),arg1,arg2)
end

function PetscObjectDereference(arg1::PetscObject)
    ccall((:PetscObjectDereference,petsc),PetscErrorCode,(PetscObject,),arg1)
end

function PetscObjectGetNewTag(arg1::PetscObject,arg2::Ptr{PetscMPIInt})
    ccall((:PetscObjectGetNewTag,petsc),PetscErrorCode,(PetscObject,Ptr{PetscMPIInt}),arg1,arg2)
end

function PetscObjectCompose(arg1::PetscObject,arg2::Ptr{Uint8},arg3::PetscObject)
    ccall((:PetscObjectCompose,petsc),PetscErrorCode,(PetscObject,Ptr{Uint8},PetscObject),arg1,arg2,arg3)
end

function PetscObjectRemoveReference(arg1::PetscObject,arg2::Ptr{Uint8})
    ccall((:PetscObjectRemoveReference,petsc),PetscErrorCode,(PetscObject,Ptr{Uint8}),arg1,arg2)
end

function PetscObjectQuery(arg1::PetscObject,arg2::Ptr{Uint8},arg3::Ptr{PetscObject})
    ccall((:PetscObjectQuery,petsc),PetscErrorCode,(PetscObject,Ptr{Uint8},Ptr{PetscObject}),arg1,arg2,arg3)
end

function PetscObjectComposeFunction_Private(arg1::PetscObject,arg2::Ptr{Uint8},arg3::Ptr{Void})
    ccall((:PetscObjectComposeFunction_Private,petsc),PetscErrorCode,(PetscObject,Ptr{Uint8},Ptr{Void}),arg1,arg2,arg3)
end

function PetscObjectSetFromOptions(arg1::PetscObject)
    ccall((:PetscObjectSetFromOptions,petsc),PetscErrorCode,(PetscObject,),arg1)
end

function PetscObjectSetUp(arg1::PetscObject)
    ccall((:PetscObjectSetUp,petsc),PetscErrorCode,(PetscObject,),arg1)
end

function PetscCommGetNewTag(arg1::MPI_Comm,arg2::Ptr{PetscMPIInt})
    ccall((:PetscCommGetNewTag,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscMPIInt}),arg1,arg2)
end

function PetscObjectAddOptionsHandler(arg1::PetscObject,arg2::Ptr{Void},arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:PetscObjectAddOptionsHandler,petsc),PetscErrorCode,(PetscObject,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function PetscObjectProcessOptionsHandlers(arg1::PetscObject)
    ccall((:PetscObjectProcessOptionsHandlers,petsc),PetscErrorCode,(PetscObject,),arg1)
end

function PetscObjectDestroyOptionsHandlers(arg1::PetscObject)
    ccall((:PetscObjectDestroyOptionsHandlers,petsc),PetscErrorCode,(PetscObject,),arg1)
end

function PetscObjectsListGetGlobalNumbering(arg1::MPI_Comm,arg2::PetscInt,arg3::Ptr{PetscObject},arg4::Ptr{PetscInt},arg5::Ptr{PetscInt})
    ccall((:PetscObjectsListGetGlobalNumbering,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{PetscObject},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end

function PetscOptionsHasName(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{PetscBool})
    ccall((:PetscOptionsHasName,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{PetscBool}),arg1,arg2,arg3)
end

function PetscOptionsGetInt(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{PetscInt},arg4::Ptr{PetscBool})
    ccall((:PetscOptionsGetInt,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function PetscOptionsGetBool(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{PetscBool},arg4::Ptr{PetscBool})
    ccall((:PetscOptionsGetBool,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function PetscOptionsGetReal(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{Cint},arg4::Ptr{PetscBool})
    ccall((:PetscOptionsGetReal,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{Cint},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function PetscOptionsGetScalar(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{PetscScalar},arg4::Ptr{PetscBool})
    ccall((:PetscOptionsGetScalar,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{PetscScalar},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function PetscOptionsGetIntArray(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{PetscBool})
    ccall((:PetscOptionsGetIntArray,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function PetscOptionsGetRealArray(arg1::Ptr{Uint8},arg2::Ptr{Uint8},PetscReal::Ptr{Cint},arg3::Ptr{PetscInt},arg4::Ptr{PetscBool})
    ccall((:PetscOptionsGetRealArray,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{Cint},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,PetscReal,arg3,arg4)
end

function PetscOptionsGetScalarArray(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{PetscScalar},arg4::Ptr{PetscInt},arg5::Ptr{PetscBool})
    ccall((:PetscOptionsGetScalarArray,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{PetscScalar},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function PetscOptionsGetBoolArray(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{PetscBool},arg4::Ptr{PetscInt},arg5::Ptr{PetscBool})
    ccall((:PetscOptionsGetBoolArray,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{PetscBool},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function PetscOptionsGetString(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{Uint8},size_t::Cint,arg4::Ptr{PetscBool})
    ccall((:PetscOptionsGetString,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Cint,Ptr{PetscBool}),arg1,arg2,arg3,size_t,arg4)
end

function PetscOptionsGetStringArray(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{Ptr{Uint8}},arg4::Ptr{PetscInt},arg5::Ptr{PetscBool})
    ccall((:PetscOptionsGetStringArray,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{Ptr{Uint8}},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function PetscOptionsGetEList(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{Ptr{Uint8}},arg4::PetscInt,arg5::Ptr{PetscInt},arg6::Ptr{PetscBool})
    ccall((:PetscOptionsGetEList,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{Ptr{Uint8}},PetscInt,Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscOptionsGetEnum(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{Ptr{Uint8}},arg4::Ptr{PetscEnum},arg5::Ptr{PetscBool})
    ccall((:PetscOptionsGetEnum,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{Ptr{Uint8}},Ptr{PetscEnum},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function PetscOptionsGetEnumArray(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{Ptr{Uint8}},arg4::Ptr{PetscEnum},arg5::Ptr{PetscInt},arg6::Ptr{PetscBool})
    ccall((:PetscOptionsGetEnumArray,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{Ptr{Uint8}},Ptr{PetscEnum},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscOptionsValidKey(arg1::Ptr{Uint8},arg2::Ptr{PetscBool})
    ccall((:PetscOptionsValidKey,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{PetscBool}),arg1,arg2)
end

function PetscOptionsSetAlias(arg1::Ptr{Uint8},arg2::Ptr{Uint8})
    ccall((:PetscOptionsSetAlias,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8}),arg1,arg2)
end

function PetscOptionsSetValue(arg1::Ptr{Uint8},arg2::Ptr{Uint8})
    ccall((:PetscOptionsSetValue,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8}),arg1,arg2)
end

function PetscOptionsClearValue(arg1::Ptr{Uint8})
    ccall((:PetscOptionsClearValue,petsc),PetscErrorCode,(Ptr{Uint8},),arg1)
end

function PetscOptionsAllUsed(arg1::Ptr{PetscInt})
    ccall((:PetscOptionsAllUsed,petsc),PetscErrorCode,(Ptr{PetscInt},),arg1)
end

function PetscOptionsUsed(arg1::Ptr{Uint8},arg2::Ptr{PetscBool})
    ccall((:PetscOptionsUsed,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{PetscBool}),arg1,arg2)
end

function PetscOptionsLeft()
    ccall((:PetscOptionsLeft,petsc),PetscErrorCode,())
end

function PetscOptionsView(arg1::PetscViewer)
    ccall((:PetscOptionsView,petsc),PetscErrorCode,(PetscViewer,),arg1)
end

function PetscOptionsCreate()
    ccall((:PetscOptionsCreate,petsc),PetscErrorCode,())
end

function PetscOptionsInsert(arg1::Ptr{Cint},arg2::Ptr{Ptr{Ptr{Uint8}}},arg3::Ptr{Uint8})
    ccall((:PetscOptionsInsert,petsc),PetscErrorCode,(Ptr{Cint},Ptr{Ptr{Ptr{Uint8}}},Ptr{Uint8}),arg1,arg2,arg3)
end

function PetscOptionsInsertFile(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::PetscBool)
    ccall((:PetscOptionsInsertFile,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},PetscBool),arg1,arg2,arg3)
end

function PetscOptionsInsertString(arg1::Ptr{Uint8})
    ccall((:PetscOptionsInsertString,petsc),PetscErrorCode,(Ptr{Uint8},),arg1)
end

function PetscOptionsDestroy()
    ccall((:PetscOptionsDestroy,petsc),PetscErrorCode,())
end

function PetscOptionsClear()
    ccall((:PetscOptionsClear,petsc),PetscErrorCode,())
end

function PetscOptionsPrefixPush(arg1::Ptr{Uint8})
    ccall((:PetscOptionsPrefixPush,petsc),PetscErrorCode,(Ptr{Uint8},),arg1)
end

function PetscOptionsPrefixPop()
    ccall((:PetscOptionsPrefixPop,petsc),PetscErrorCode,())
end

function PetscOptionsReject(arg1::Ptr{Uint8},arg2::Ptr{Uint8})
    ccall((:PetscOptionsReject,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8}),arg1,arg2)
end

function PetscOptionsGetAll(arg1::Ptr{Ptr{Uint8}})
    ccall((:PetscOptionsGetAll,petsc),PetscErrorCode,(Ptr{Ptr{Uint8}},),arg1)
end

function PetscOptionsGetenv(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},size_t::Cint,arg4::Ptr{PetscBool})
    ccall((:PetscOptionsGetenv,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Cint,Ptr{PetscBool}),arg1,arg2,arg3,size_t,arg4)
end

function PetscOptionsStringToInt(arg1::Ptr{Uint8},arg2::Ptr{PetscInt})
    ccall((:PetscOptionsStringToInt,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{PetscInt}),arg1,arg2)
end

function PetscOptionsStringToReal(arg1::Ptr{Uint8},arg2::Ptr{Cint})
    ccall((:PetscOptionsStringToReal,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Cint}),arg1,arg2)
end

function PetscOptionsStringToBool(arg1::Ptr{Uint8},arg2::Ptr{PetscBool})
    ccall((:PetscOptionsStringToBool,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{PetscBool}),arg1,arg2)
end

function PetscOptionsMonitorSet(arg1::Ptr{Void},arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:PetscOptionsMonitorSet,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function PetscOptionsMonitorCancel()
    ccall((:PetscOptionsMonitorCancel,petsc),PetscErrorCode,())
end

function PetscOptionsMonitorDefault(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{Void})
    ccall((:PetscOptionsMonitorDefault,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{Void}),arg1,arg2,arg3)
end

function PetscOptionsBegin_Private(arg1::Ptr{PetscOptions},arg2::MPI_Comm,arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Ptr{Uint8})
    ccall((:PetscOptionsBegin_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},MPI_Comm,Ptr{Uint8},Ptr{Uint8},Ptr{Uint8}),arg1,arg2,arg3,arg4,arg5)
end

function PetscObjectOptionsBegin_Private(arg1::Ptr{PetscOptions},arg2::PetscObject)
    ccall((:PetscObjectOptionsBegin_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},PetscObject),arg1,arg2)
end

function PetscOptionsEnd_Private(arg1::Ptr{PetscOptions})
    ccall((:PetscOptionsEnd_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},),arg1)
end

function PetscOptionsHead(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8})
    ccall((:PetscOptionsHead,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8}),arg1,arg2)
end

function PetscOptionsEnum_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Ptr{Ptr{Uint8}},arg6::PetscEnum,arg7::Ptr{PetscEnum},arg8::Ptr{PetscBool})
    ccall((:PetscOptionsEnum_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Ptr{Uint8}},PetscEnum,Ptr{PetscEnum},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscOptionsInt_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::PetscInt,arg6::Ptr{PetscInt},arg7::Ptr{PetscBool})
    ccall((:PetscOptionsInt_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},PetscInt,Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscOptionsReal_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},PetscReal::Cint,arg5::Ptr{Cint},arg6::Ptr{PetscBool})
    ccall((:PetscOptionsReal_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Cint,Ptr{Cint},Ptr{PetscBool}),arg1,arg2,arg3,arg4,PetscReal,arg5,arg6)
end

function PetscOptionsScalar_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::PetscScalar,arg6::Ptr{PetscScalar},arg7::Ptr{PetscBool})
    ccall((:PetscOptionsScalar_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},PetscScalar,Ptr{PetscScalar},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscOptionsName_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Ptr{PetscBool})
    ccall((:PetscOptionsName_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function PetscOptionsString_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Ptr{Uint8},arg6::Ptr{Uint8},size_t::Cint,arg7::Ptr{PetscBool})
    ccall((:PetscOptionsString_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Cint,Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,size_t,arg7)
end

function PetscOptionsBool_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::PetscBool,arg6::Ptr{PetscBool},arg7::Ptr{PetscBool})
    ccall((:PetscOptionsBool_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},PetscBool,Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscOptionsBoolGroupBegin_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Ptr{PetscBool})
    ccall((:PetscOptionsBoolGroupBegin_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function PetscOptionsBoolGroup_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Ptr{PetscBool})
    ccall((:PetscOptionsBoolGroup_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function PetscOptionsBoolGroupEnd_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Ptr{PetscBool})
    ccall((:PetscOptionsBoolGroupEnd_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function PetscOptionsFList_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::PetscFunctionList,arg6::Ptr{Uint8},arg7::Ptr{Uint8},size_t::Cint,arg8::Ptr{PetscBool})
    ccall((:PetscOptionsFList_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},PetscFunctionList,Ptr{Uint8},Ptr{Uint8},Cint,Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,size_t,arg8)
end

function PetscOptionsEList_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Ptr{Ptr{Uint8}},arg6::PetscInt,arg7::Ptr{Uint8},arg8::Ptr{PetscInt},arg9::Ptr{PetscBool})
    ccall((:PetscOptionsEList_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Ptr{Uint8}},PetscInt,Ptr{Uint8},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function PetscOptionsRealArray_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},PetscReal::Ptr{Cint},arg5::Ptr{PetscInt},arg6::Ptr{PetscBool})
    ccall((:PetscOptionsRealArray_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Cint},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,PetscReal,arg5,arg6)
end

function PetscOptionsScalarArray_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Ptr{PetscScalar},arg6::Ptr{PetscInt},arg7::Ptr{PetscBool})
    ccall((:PetscOptionsScalarArray_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{PetscScalar},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscOptionsIntArray_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Ptr{PetscInt},arg6::Ptr{PetscInt},arg7::Ptr{PetscBool})
    ccall((:PetscOptionsIntArray_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscOptionsStringArray_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Ptr{Ptr{Uint8}},arg6::Ptr{PetscInt},arg7::Ptr{PetscBool})
    ccall((:PetscOptionsStringArray_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Ptr{Uint8}},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscOptionsBoolArray_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Ptr{PetscBool},arg6::Ptr{PetscInt},arg7::Ptr{PetscBool})
    ccall((:PetscOptionsBoolArray_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{PetscBool},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscOptionsEnumArray_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Ptr{Ptr{Uint8}},arg6::Ptr{PetscEnum},arg7::Ptr{PetscInt},arg8::Ptr{PetscBool})
    ccall((:PetscOptionsEnumArray_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Ptr{Uint8}},Ptr{PetscEnum},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscOptionsSetFromOptions()
    ccall((:PetscOptionsSetFromOptions,petsc),PetscErrorCode,())
end

function PetscOptionsSAWsDestroy()
    ccall((:PetscOptionsSAWsDestroy,petsc),PetscErrorCode,())
end

function PetscMemoryShowUsage(arg1::PetscViewer,arg2::Ptr{Uint8})
    ccall((:PetscMemoryShowUsage,petsc),PetscErrorCode,(PetscViewer,Ptr{Uint8}),arg1,arg2)
end

function PetscObjectPrintClassNamePrefixType(arg1::PetscObject,arg2::PetscViewer)
    ccall((:PetscObjectPrintClassNamePrefixType,petsc),PetscErrorCode,(PetscObject,PetscViewer),arg1,arg2)
end

function PetscObjectView(arg1::PetscObject,arg2::PetscViewer)
    ccall((:PetscObjectView,petsc),PetscErrorCode,(PetscObject,PetscViewer),arg1,arg2)
end

function PetscObjectQueryFunction_Private(arg1::PetscObject,arg2::Ptr{Uint8},arg3::Ptr{Ptr{Void}})
    ccall((:PetscObjectQueryFunction_Private,petsc),PetscErrorCode,(PetscObject,Ptr{Uint8},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end

function PetscObjectSetOptionsPrefix(arg1::PetscObject,arg2::Ptr{Uint8})
    ccall((:PetscObjectSetOptionsPrefix,petsc),PetscErrorCode,(PetscObject,Ptr{Uint8}),arg1,arg2)
end

function PetscObjectAppendOptionsPrefix(arg1::PetscObject,arg2::Ptr{Uint8})
    ccall((:PetscObjectAppendOptionsPrefix,petsc),PetscErrorCode,(PetscObject,Ptr{Uint8}),arg1,arg2)
end

function PetscObjectPrependOptionsPrefix(arg1::PetscObject,arg2::Ptr{Uint8})
    ccall((:PetscObjectPrependOptionsPrefix,petsc),PetscErrorCode,(PetscObject,Ptr{Uint8}),arg1,arg2)
end

function PetscObjectGetOptionsPrefix(arg1::PetscObject,arg2::Ptr{Ptr{Uint8}})
    ccall((:PetscObjectGetOptionsPrefix,petsc),PetscErrorCode,(PetscObject,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function PetscObjectChangeTypeName(arg1::PetscObject,arg2::Ptr{Uint8})
    ccall((:PetscObjectChangeTypeName,petsc),PetscErrorCode,(PetscObject,Ptr{Uint8}),arg1,arg2)
end

function PetscObjectRegisterDestroy(arg1::PetscObject)
    ccall((:PetscObjectRegisterDestroy,petsc),PetscErrorCode,(PetscObject,),arg1)
end

function PetscObjectRegisterDestroyAll()
    ccall((:PetscObjectRegisterDestroyAll,petsc),PetscErrorCode,())
end

function PetscObjectViewFromOptions(arg1::PetscObject,arg2::PetscObject,arg3::Ptr{Uint8})
    ccall((:PetscObjectViewFromOptions,petsc),PetscErrorCode,(PetscObject,PetscObject,Ptr{Uint8}),arg1,arg2,arg3)
end

function PetscObjectName(arg1::PetscObject)
    ccall((:PetscObjectName,petsc),PetscErrorCode,(PetscObject,),arg1)
end

function PetscObjectTypeCompare(arg1::PetscObject,arg2::Ptr{Uint8},arg3::Ptr{PetscBool})
    ccall((:PetscObjectTypeCompare,petsc),PetscErrorCode,(PetscObject,Ptr{Uint8},Ptr{PetscBool}),arg1,arg2,arg3)
end

function PetscRegisterFinalize(arg1::Ptr{Void})
    ccall((:PetscRegisterFinalize,petsc),PetscErrorCode,(Ptr{Void},),arg1)
end

function PetscRegisterFinalizeAll()
    ccall((:PetscRegisterFinalizeAll,petsc),PetscErrorCode,())
end

function PetscDLOpen(arg1::Ptr{Uint8},arg2::PetscDLMode,arg3::Ptr{PetscDLHandle})
    ccall((:PetscDLOpen,petsc),PetscErrorCode,(Ptr{Uint8},PetscDLMode,Ptr{PetscDLHandle}),arg1,arg2,arg3)
end

function PetscDLClose(arg1::Ptr{PetscDLHandle})
    ccall((:PetscDLClose,petsc),PetscErrorCode,(Ptr{PetscDLHandle},),arg1)
end

function PetscDLSym(arg1::PetscDLHandle,arg2::Ptr{Uint8},arg3::Ptr{Ptr{Void}})
    ccall((:PetscDLSym,petsc),PetscErrorCode,(PetscDLHandle,Ptr{Uint8},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end

function PetscObjectsDump(arg1::Ptr{Void},arg2::PetscBool)
    ccall((:PetscObjectsDump,petsc),PetscErrorCode,(Ptr{Void},PetscBool),arg1,arg2)
end

function PetscObjectListDestroy(arg1::Ptr{PetscObjectList})
    ccall((:PetscObjectListDestroy,petsc),PetscErrorCode,(Ptr{PetscObjectList},),arg1)
end

function PetscObjectListFind(arg1::PetscObjectList,arg2::Ptr{Uint8},arg3::Ptr{PetscObject})
    ccall((:PetscObjectListFind,petsc),PetscErrorCode,(PetscObjectList,Ptr{Uint8},Ptr{PetscObject}),arg1,arg2,arg3)
end

function PetscObjectListReverseFind(arg1::PetscObjectList,arg2::PetscObject,arg3::Ptr{Ptr{Uint8}},arg4::Ptr{PetscBool})
    ccall((:PetscObjectListReverseFind,petsc),PetscErrorCode,(PetscObjectList,PetscObject,Ptr{Ptr{Uint8}},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function PetscObjectListAdd(arg1::Ptr{PetscObjectList},arg2::Ptr{Uint8},arg3::PetscObject)
    ccall((:PetscObjectListAdd,petsc),PetscErrorCode,(Ptr{PetscObjectList},Ptr{Uint8},PetscObject),arg1,arg2,arg3)
end

function PetscObjectListRemoveReference(arg1::Ptr{PetscObjectList},arg2::Ptr{Uint8})
    ccall((:PetscObjectListRemoveReference,petsc),PetscErrorCode,(Ptr{PetscObjectList},Ptr{Uint8}),arg1,arg2)
end

function PetscObjectListDuplicate(arg1::PetscObjectList,arg2::Ptr{PetscObjectList})
    ccall((:PetscObjectListDuplicate,petsc),PetscErrorCode,(PetscObjectList,Ptr{PetscObjectList}),arg1,arg2)
end

function PetscFunctionListAdd_Private(arg1::Ptr{PetscFunctionList},arg2::Ptr{Uint8},arg3::Ptr{Void})
    ccall((:PetscFunctionListAdd_Private,petsc),PetscErrorCode,(Ptr{PetscFunctionList},Ptr{Uint8},Ptr{Void}),arg1,arg2,arg3)
end

function PetscFunctionListDestroy(arg1::Ptr{PetscFunctionList})
    ccall((:PetscFunctionListDestroy,petsc),PetscErrorCode,(Ptr{PetscFunctionList},),arg1)
end

function PetscFunctionListFind_Private(arg1::PetscFunctionList,arg2::Ptr{Uint8},arg3::Ptr{Ptr{Void}})
    ccall((:PetscFunctionListFind_Private,petsc),PetscErrorCode,(PetscFunctionList,Ptr{Uint8},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end

function PetscFunctionListPrintTypes(arg1::MPI_Comm,arg2::Ptr{Void},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Ptr{Uint8},arg6::Ptr{Uint8},arg7::PetscFunctionList,arg8::Ptr{Uint8})
    ccall((:PetscFunctionListPrintTypes,petsc),PetscErrorCode,(MPI_Comm,Ptr{Void},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},PetscFunctionList,Ptr{Uint8}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscFunctionListDuplicate(arg1::PetscFunctionList,arg2::Ptr{PetscFunctionList})
    ccall((:PetscFunctionListDuplicate,petsc),PetscErrorCode,(PetscFunctionList,Ptr{PetscFunctionList}),arg1,arg2)
end

function PetscFunctionListView(arg1::PetscFunctionList,arg2::PetscViewer)
    ccall((:PetscFunctionListView,petsc),PetscErrorCode,(PetscFunctionList,PetscViewer),arg1,arg2)
end

function PetscFunctionListGet(arg1::PetscFunctionList,arg2::Ptr{Ptr{Ptr{Uint8}}},arg3::Ptr{Cint})
    ccall((:PetscFunctionListGet,petsc),PetscErrorCode,(PetscFunctionList,Ptr{Ptr{Ptr{Uint8}}},Ptr{Cint}),arg1,arg2,arg3)
end

function PetscDLLibraryAppend(arg1::MPI_Comm,arg2::Ptr{PetscDLLibrary},arg3::Ptr{Uint8})
    ccall((:PetscDLLibraryAppend,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscDLLibrary},Ptr{Uint8}),arg1,arg2,arg3)
end

function PetscDLLibraryPrepend(arg1::MPI_Comm,arg2::Ptr{PetscDLLibrary},arg3::Ptr{Uint8})
    ccall((:PetscDLLibraryPrepend,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscDLLibrary},Ptr{Uint8}),arg1,arg2,arg3)
end

function PetscDLLibrarySym(arg1::MPI_Comm,arg2::Ptr{PetscDLLibrary},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Ptr{Ptr{Void}})
    ccall((:PetscDLLibrarySym,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscDLLibrary},Ptr{Uint8},Ptr{Uint8},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4,arg5)
end

function PetscDLLibraryPrintPath(arg1::PetscDLLibrary)
    ccall((:PetscDLLibraryPrintPath,petsc),PetscErrorCode,(PetscDLLibrary,),arg1)
end

function PetscDLLibraryRetrieve(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},size_t::Cint,arg4::Ptr{PetscBool})
    ccall((:PetscDLLibraryRetrieve,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Cint,Ptr{PetscBool}),arg1,arg2,arg3,size_t,arg4)
end

function PetscDLLibraryOpen(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{PetscDLLibrary})
    ccall((:PetscDLLibraryOpen,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{PetscDLLibrary}),arg1,arg2,arg3)
end

function PetscDLLibraryClose(arg1::PetscDLLibrary)
    ccall((:PetscDLLibraryClose,petsc),PetscErrorCode,(PetscDLLibrary,),arg1)
end

function PetscSplitOwnership(arg1::MPI_Comm,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt})
    ccall((:PetscSplitOwnership,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscSplitOwnershipBlock(arg1::MPI_Comm,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:PetscSplitOwnershipBlock,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function PetscSequentialPhaseBegin(arg1::MPI_Comm,arg2::PetscMPIInt)
    ccall((:PetscSequentialPhaseBegin,petsc),PetscErrorCode,(MPI_Comm,PetscMPIInt),arg1,arg2)
end

function PetscSequentialPhaseEnd(arg1::MPI_Comm,arg2::PetscMPIInt)
    ccall((:PetscSequentialPhaseEnd,petsc),PetscErrorCode,(MPI_Comm,PetscMPIInt),arg1,arg2)
end

function PetscBarrier(arg1::PetscObject)
    ccall((:PetscBarrier,petsc),PetscErrorCode,(PetscObject,),arg1)
end

function PetscMPIDump(arg1::Ptr{Void})
    ccall((:PetscMPIDump,petsc),PetscErrorCode,(Ptr{Void},),arg1)
end

function PetscInfoDeactivateClass(arg1::PetscClassId)
    ccall((:PetscInfoDeactivateClass,petsc),PetscErrorCode,(PetscClassId,),arg1)
end

function PetscInfoActivateClass(arg1::PetscClassId)
    ccall((:PetscInfoActivateClass,petsc),PetscErrorCode,(PetscClassId,),arg1)
end

function PetscLogGetStageLog(arg1::Ptr{PetscStageLog})
    ccall((:PetscLogGetStageLog,petsc),PetscErrorCode,(Ptr{PetscStageLog},),arg1)
end

function PetscStageLogGetCurrent(arg1::PetscStageLog,arg2::Ptr{Cint})
    ccall((:PetscStageLogGetCurrent,petsc),PetscErrorCode,(PetscStageLog,Ptr{Cint}),arg1,arg2)
end

function PetscStageLogGetEventPerfLog(arg1::PetscStageLog,arg2::Cint,arg3::Ptr{PetscEventPerfLog})
    ccall((:PetscStageLogGetEventPerfLog,petsc),PetscErrorCode,(PetscStageLog,Cint,Ptr{PetscEventPerfLog}),arg1,arg2,arg3)
end

function PetscLogObjectParent(arg1::PetscObject,arg2::PetscObject)
    ccall((:PetscLogObjectParent,petsc),PetscErrorCode,(PetscObject,PetscObject),arg1,arg2)
end

function PetscLogObjectMemory(arg1::PetscObject,arg2::PetscLogDouble)
    ccall((:PetscLogObjectMemory,petsc),PetscErrorCode,(PetscObject,PetscLogDouble),arg1,arg2)
end

function PetscIntStackCreate(arg1::Ptr{PetscIntStack})
    ccall((:PetscIntStackCreate,petsc),PetscErrorCode,(Ptr{PetscIntStack},),arg1)
end

function PetscIntStackDestroy(arg1::PetscIntStack)
    ccall((:PetscIntStackDestroy,petsc),PetscErrorCode,(PetscIntStack,),arg1)
end

function PetscIntStackPush(arg1::PetscIntStack,arg2::Cint)
    ccall((:PetscIntStackPush,petsc),PetscErrorCode,(PetscIntStack,Cint),arg1,arg2)
end

function PetscIntStackPop(arg1::PetscIntStack,arg2::Ptr{Cint})
    ccall((:PetscIntStackPop,petsc),PetscErrorCode,(PetscIntStack,Ptr{Cint}),arg1,arg2)
end

function PetscIntStackTop(arg1::PetscIntStack,arg2::Ptr{Cint})
    ccall((:PetscIntStackTop,petsc),PetscErrorCode,(PetscIntStack,Ptr{Cint}),arg1,arg2)
end

function PetscIntStackEmpty(arg1::PetscIntStack,arg2::Ptr{PetscBool})
    ccall((:PetscIntStackEmpty,petsc),PetscErrorCode,(PetscIntStack,Ptr{PetscBool}),arg1,arg2)
end

function PetscFixFilename(arg1::Ptr{Uint8},arg2::Ptr{Uint8})
    ccall((:PetscFixFilename,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8}),arg1,arg2)
end

function PetscFOpen(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Ptr{Void}})
    ccall((:PetscFOpen,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4)
end

function PetscFClose(arg1::MPI_Comm,arg2::Ptr{Void})
    ccall((:PetscFClose,petsc),PetscErrorCode,(MPI_Comm,Ptr{Void}),arg1,arg2)
end

#=
function PetscVSNPrintf(arg1::Ptr{Uint8},size_t::Cint,arg2::Ptr{Uint8},arg3::Ptr{Cint},arg4::va_list)
    ccall((:PetscVSNPrintf,petsc),PetscErrorCode,(Ptr{Uint8},Cint,Ptr{Uint8},Ptr{Cint},va_list),arg1,size_t,arg2,arg3,arg4)
end


function PetscVFPrintfDefault(arg1::Ptr{Void},arg2::Ptr{Uint8},arg3::va_list)
    ccall((:PetscVFPrintfDefault,petsc),PetscErrorCode,(Ptr{Void},Ptr{Uint8},va_list),arg1,arg2,arg3)
end

=#

function PetscSynchronizedFlush(arg1::MPI_Comm,arg2::Ptr{Void})
    ccall((:PetscSynchronizedFlush,petsc),PetscErrorCode,(MPI_Comm,Ptr{Void}),arg1,arg2)
end

function PetscSynchronizedFGets(arg1::MPI_Comm,arg2::Ptr{Void},size_t::Cint,arg3::Ptr{Uint8})
    ccall((:PetscSynchronizedFGets,petsc),PetscErrorCode,(MPI_Comm,Ptr{Void},Cint,Ptr{Uint8}),arg1,arg2,size_t,arg3)
end

function PetscStartMatlab(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Ptr{Void}})
    ccall((:PetscStartMatlab,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4)
end

function PetscStartJava(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Ptr{Void}})
    ccall((:PetscStartJava,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4)
end

function PetscGetPetscDir(arg1::Ptr{Ptr{Uint8}})
    ccall((:PetscGetPetscDir,petsc),PetscErrorCode,(Ptr{Ptr{Uint8}},),arg1)
end

function PetscPopUpSelect(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Cint,arg5::Ptr{Ptr{Uint8}},arg6::Ptr{Cint})
    ccall((:PetscPopUpSelect,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Cint,Ptr{Ptr{Uint8}},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscContainerGetPointer(arg1::PetscContainer,arg2::Ptr{Ptr{Void}})
    ccall((:PetscContainerGetPointer,petsc),PetscErrorCode,(PetscContainer,Ptr{Ptr{Void}}),arg1,arg2)
end

function PetscContainerSetPointer(arg1::PetscContainer,arg2::Ptr{Void})
    ccall((:PetscContainerSetPointer,petsc),PetscErrorCode,(PetscContainer,Ptr{Void}),arg1,arg2)
end

function PetscContainerDestroy(arg1::Ptr{PetscContainer})
    ccall((:PetscContainerDestroy,petsc),PetscErrorCode,(Ptr{PetscContainer},),arg1)
end

function PetscContainerCreate(arg1::MPI_Comm,arg2::Ptr{PetscContainer})
    ccall((:PetscContainerCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscContainer}),arg1,arg2)
end

function PetscContainerSetUserDestroy(arg1::PetscContainer,arg2::Ptr{Void})
    ccall((:PetscContainerSetUserDestroy,petsc),PetscErrorCode,(PetscContainer,Ptr{Void}),arg1,arg2)
end

function PetscIntView(arg1::PetscInt,arg2::Ptr{PetscInt},arg3::PetscViewer)
    ccall((:PetscIntView,petsc),PetscErrorCode,(PetscInt,Ptr{PetscInt},PetscViewer),arg1,arg2,arg3)
end

function PetscRealView(arg1::PetscInt,PetscReal::Ptr{Cint},arg2::PetscViewer)
    ccall((:PetscRealView,petsc),PetscErrorCode,(PetscInt,Ptr{Cint},PetscViewer),arg1,PetscReal,arg2)
end

function PetscScalarView(arg1::PetscInt,arg2::Ptr{PetscScalar},arg3::PetscViewer)
    ccall((:PetscScalarView,petsc),PetscErrorCode,(PetscInt,Ptr{PetscScalar},PetscViewer),arg1,arg2,arg3)
end

function PetscGetHostName(arg1::Ptr{Uint8},size_t::Cint)
    ccall((:PetscGetHostName,petsc),PetscErrorCode,(Ptr{Uint8},Cint),arg1,size_t)
end

function PetscGetUserName(arg1::Ptr{Uint8},size_t::Cint)
    ccall((:PetscGetUserName,petsc),PetscErrorCode,(Ptr{Uint8},Cint),arg1,size_t)
end

function PetscGetProgramName(arg1::Ptr{Uint8},size_t::Cint)
    ccall((:PetscGetProgramName,petsc),PetscErrorCode,(Ptr{Uint8},Cint),arg1,size_t)
end

function PetscSetProgramName(arg1::Ptr{Uint8})
    ccall((:PetscSetProgramName,petsc),PetscErrorCode,(Ptr{Uint8},),arg1)
end

function PetscGetDate(arg1::Ptr{Uint8},size_t::Cint)
    ccall((:PetscGetDate,petsc),PetscErrorCode,(Ptr{Uint8},Cint),arg1,size_t)
end

function PetscGetVersion(arg1::Ptr{Uint8},size_t::Cint)
    ccall((:PetscGetVersion,petsc),PetscErrorCode,(Ptr{Uint8},Cint),arg1,size_t)
end

function PetscSortInt(arg1::PetscInt,arg2::Ptr{PetscInt})
    ccall((:PetscSortInt,petsc),PetscErrorCode,(PetscInt,Ptr{PetscInt}),arg1,arg2)
end

function PetscSortRemoveDupsInt(arg1::Ptr{PetscInt},arg2::Ptr{PetscInt})
    ccall((:PetscSortRemoveDupsInt,petsc),PetscErrorCode,(Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2)
end

function PetscFindInt(arg1::PetscInt,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:PetscFindInt,petsc),PetscErrorCode,(PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function PetscSortIntWithPermutation(arg1::PetscInt,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt})
    ccall((:PetscSortIntWithPermutation,petsc),PetscErrorCode,(PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscSortStrWithPermutation(arg1::PetscInt,arg2::Ptr{Ptr{Uint8}},arg3::Ptr{PetscInt})
    ccall((:PetscSortStrWithPermutation,petsc),PetscErrorCode,(PetscInt,Ptr{Ptr{Uint8}},Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscSortIntWithArray(arg1::PetscInt,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt})
    ccall((:PetscSortIntWithArray,petsc),PetscErrorCode,(PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscSortIntWithArrayPair(arg1::PetscInt,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:PetscSortIntWithArrayPair,petsc),PetscErrorCode,(PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function PetscSortMPIInt(arg1::PetscInt,arg2::Ptr{PetscMPIInt})
    ccall((:PetscSortMPIInt,petsc),PetscErrorCode,(PetscInt,Ptr{PetscMPIInt}),arg1,arg2)
end

function PetscSortRemoveDupsMPIInt(arg1::Ptr{PetscInt},arg2::Ptr{PetscMPIInt})
    ccall((:PetscSortRemoveDupsMPIInt,petsc),PetscErrorCode,(Ptr{PetscInt},Ptr{PetscMPIInt}),arg1,arg2)
end

function PetscSortMPIIntWithArray(arg1::PetscMPIInt,arg2::Ptr{PetscMPIInt},arg3::Ptr{PetscMPIInt})
    ccall((:PetscSortMPIIntWithArray,petsc),PetscErrorCode,(PetscMPIInt,Ptr{PetscMPIInt},Ptr{PetscMPIInt}),arg1,arg2,arg3)
end

function PetscSortIntWithScalarArray(arg1::PetscInt,arg2::Ptr{PetscInt},arg3::Ptr{PetscScalar})
    ccall((:PetscSortIntWithScalarArray,petsc),PetscErrorCode,(PetscInt,Ptr{PetscInt},Ptr{PetscScalar}),arg1,arg2,arg3)
end

function PetscSortIntWithDataArray(arg1::PetscInt,arg2::Ptr{PetscInt},arg3::Ptr{Void},size_t::Cint,arg4::Ptr{Void})
    ccall((:PetscSortIntWithDataArray,petsc),PetscErrorCode,(PetscInt,Ptr{PetscInt},Ptr{Void},Cint,Ptr{Void}),arg1,arg2,arg3,size_t,arg4)
end

function PetscSortReal(arg1::PetscInt,PetscReal::Ptr{Cint})
    ccall((:PetscSortReal,petsc),PetscErrorCode,(PetscInt,Ptr{Cint}),arg1,PetscReal)
end

function PetscSortRealWithPermutation(arg1::PetscInt,PetscReal::Ptr{Cint},arg2::Ptr{PetscInt})
    ccall((:PetscSortRealWithPermutation,petsc),PetscErrorCode,(PetscInt,Ptr{Cint},Ptr{PetscInt}),arg1,PetscReal,arg2)
end

function PetscSortRemoveDupsReal(arg1::Ptr{PetscInt},PetscReal::Ptr{Cint})
    ccall((:PetscSortRemoveDupsReal,petsc),PetscErrorCode,(Ptr{PetscInt},Ptr{Cint}),arg1,PetscReal)
end

function PetscSortSplit(arg1::PetscInt,arg2::PetscInt,arg3::Ptr{PetscScalar},arg4::Ptr{PetscInt})
    ccall((:PetscSortSplit,petsc),PetscErrorCode,(PetscInt,PetscInt,Ptr{PetscScalar},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function PetscSortSplitReal(arg1::PetscInt,arg2::PetscInt,PetscReal::Ptr{Cint},arg3::Ptr{PetscInt})
    ccall((:PetscSortSplitReal,petsc),PetscErrorCode,(PetscInt,PetscInt,Ptr{Cint},Ptr{PetscInt}),arg1,arg2,PetscReal,arg3)
end

function PetscProcessTree(arg1::PetscInt,arg2::Ptr{PetscBool},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{Ptr{PetscInt}},arg6::Ptr{Ptr{PetscInt}},arg7::Ptr{Ptr{PetscInt}},arg8::Ptr{Ptr{PetscInt}})
    ccall((:PetscProcessTree,petsc),PetscErrorCode,(PetscInt,Ptr{PetscBool},Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscMergeIntArrayPair(arg1::PetscInt,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::PetscInt,arg5::Ptr{PetscInt},arg6::Ptr{PetscInt},arg7::Ptr{PetscInt},arg8::Ptr{Ptr{PetscInt}},arg9::Ptr{Ptr{PetscInt}})
    ccall((:PetscMergeIntArrayPair,petsc),PetscErrorCode,(PetscInt,Ptr{PetscInt},Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function PetscMergeIntArray(arg1::PetscInt,arg2::Ptr{PetscInt},arg3::PetscInt,arg4::Ptr{PetscInt},arg5::Ptr{PetscInt},arg6::Ptr{Ptr{PetscInt}})
    ccall((:PetscMergeIntArray,petsc),PetscErrorCode,(PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscSetDisplay()
    ccall((:PetscSetDisplay,petsc),PetscErrorCode,())
end

function PetscGetDisplay(arg1::Ptr{Uint8},size_t::Cint)
    ccall((:PetscGetDisplay,petsc),PetscErrorCode,(Ptr{Uint8},Cint),arg1,size_t)
end

function PetscRandomInitializePackage()
    ccall((:PetscRandomInitializePackage,petsc),PetscErrorCode,())
end

function PetscRandomRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:PetscRandomRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function PetscRandomSetType(arg1::PetscRandom,arg2::PetscRandomType)
    ccall((:PetscRandomSetType,petsc),PetscErrorCode,(PetscRandom,PetscRandomType),arg1,arg2)
end

function PetscRandomSetFromOptions(arg1::PetscRandom)
    ccall((:PetscRandomSetFromOptions,petsc),PetscErrorCode,(PetscRandom,),arg1)
end

function PetscRandomGetType(arg1::PetscRandom,arg2::Ptr{PetscRandomType})
    ccall((:PetscRandomGetType,petsc),PetscErrorCode,(PetscRandom,Ptr{PetscRandomType}),arg1,arg2)
end

function PetscRandomCreate(arg1::MPI_Comm,arg2::Ptr{PetscRandom})
    ccall((:PetscRandomCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscRandom}),arg1,arg2)
end

function PetscRandomGetValue(arg1::PetscRandom,arg2::Ptr{PetscScalar})
    ccall((:PetscRandomGetValue,petsc),PetscErrorCode,(PetscRandom,Ptr{PetscScalar}),arg1,arg2)
end

function PetscRandomGetValueReal(arg1::PetscRandom,arg2::Ptr{Cint})
    ccall((:PetscRandomGetValueReal,petsc),PetscErrorCode,(PetscRandom,Ptr{Cint}),arg1,arg2)
end

function PetscRandomGetInterval(arg1::PetscRandom,arg2::Ptr{PetscScalar},arg3::Ptr{PetscScalar})
    ccall((:PetscRandomGetInterval,petsc),PetscErrorCode,(PetscRandom,Ptr{PetscScalar},Ptr{PetscScalar}),arg1,arg2,arg3)
end

function PetscRandomSetInterval(arg1::PetscRandom,arg2::PetscScalar,arg3::PetscScalar)
    ccall((:PetscRandomSetInterval,petsc),PetscErrorCode,(PetscRandom,PetscScalar,PetscScalar),arg1,arg2,arg3)
end

function PetscRandomSetSeed(arg1::PetscRandom,arg2::Culong)
    ccall((:PetscRandomSetSeed,petsc),PetscErrorCode,(PetscRandom,Culong),arg1,arg2)
end

function PetscRandomGetSeed(arg1::PetscRandom,arg2::Ptr{Culong})
    ccall((:PetscRandomGetSeed,petsc),PetscErrorCode,(PetscRandom,Ptr{Culong}),arg1,arg2)
end

function PetscRandomSeed(arg1::PetscRandom)
    ccall((:PetscRandomSeed,petsc),PetscErrorCode,(PetscRandom,),arg1)
end

function PetscRandomDestroy(arg1::Ptr{PetscRandom})
    ccall((:PetscRandomDestroy,petsc),PetscErrorCode,(Ptr{PetscRandom},),arg1)
end

function PetscGetFullPath(arg1::Ptr{Uint8},arg2::Ptr{Uint8},size_t::Cint)
    ccall((:PetscGetFullPath,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Cint),arg1,arg2,size_t)
end

function PetscGetRelativePath(arg1::Ptr{Uint8},arg2::Ptr{Uint8},size_t::Cint)
    ccall((:PetscGetRelativePath,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Cint),arg1,arg2,size_t)
end

function PetscGetWorkingDirectory(arg1::Ptr{Uint8},size_t::Cint)
    ccall((:PetscGetWorkingDirectory,petsc),PetscErrorCode,(Ptr{Uint8},Cint),arg1,size_t)
end

function PetscGetRealPath(arg1::Ptr{Uint8},arg2::Ptr{Uint8})
    ccall((:PetscGetRealPath,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8}),arg1,arg2)
end

function PetscGetHomeDirectory(arg1::Ptr{Uint8},size_t::Cint)
    ccall((:PetscGetHomeDirectory,petsc),PetscErrorCode,(Ptr{Uint8},Cint),arg1,size_t)
end

function PetscTestFile(arg1::Ptr{Uint8},arg2::Uint8,arg3::Ptr{PetscBool})
    ccall((:PetscTestFile,petsc),PetscErrorCode,(Ptr{Uint8},Uint8,Ptr{PetscBool}),arg1,arg2,arg3)
end

function PetscTestDirectory(arg1::Ptr{Uint8},arg2::Uint8,arg3::Ptr{PetscBool})
    ccall((:PetscTestDirectory,petsc),PetscErrorCode,(Ptr{Uint8},Uint8,Ptr{PetscBool}),arg1,arg2,arg3)
end

function PetscBinaryRead(arg1::Cint,arg2::Ptr{Void},arg3::PetscInt,arg4::PetscDataType)
    ccall((:PetscBinaryRead,petsc),PetscErrorCode,(Cint,Ptr{Void},PetscInt,PetscDataType),arg1,arg2,arg3,arg4)
end

function PetscBinarySynchronizedRead(arg1::MPI_Comm,arg2::Cint,arg3::Ptr{Void},arg4::PetscInt,arg5::PetscDataType)
    ccall((:PetscBinarySynchronizedRead,petsc),PetscErrorCode,(MPI_Comm,Cint,Ptr{Void},PetscInt,PetscDataType),arg1,arg2,arg3,arg4,arg5)
end

function PetscBinarySynchronizedWrite(arg1::MPI_Comm,arg2::Cint,arg3::Ptr{Void},arg4::PetscInt,arg5::PetscDataType,arg6::PetscBool)
    ccall((:PetscBinarySynchronizedWrite,petsc),PetscErrorCode,(MPI_Comm,Cint,Ptr{Void},PetscInt,PetscDataType,PetscBool),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscBinaryWrite(arg1::Cint,arg2::Ptr{Void},arg3::PetscInt,arg4::PetscDataType,arg5::PetscBool)
    ccall((:PetscBinaryWrite,petsc),PetscErrorCode,(Cint,Ptr{Void},PetscInt,PetscDataType,PetscBool),arg1,arg2,arg3,arg4,arg5)
end

function PetscBinaryOpen(arg1::Ptr{Uint8},arg2::PetscFileMode,arg3::Ptr{Cint})
    ccall((:PetscBinaryOpen,petsc),PetscErrorCode,(Ptr{Uint8},PetscFileMode,Ptr{Cint}),arg1,arg2,arg3)
end

function PetscBinaryClose(arg1::Cint)
    ccall((:PetscBinaryClose,petsc),PetscErrorCode,(Cint,),arg1)
end

function PetscSharedTmp(arg1::MPI_Comm,arg2::Ptr{PetscBool})
    ccall((:PetscSharedTmp,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscBool}),arg1,arg2)
end

function PetscSharedWorkingDirectory(arg1::MPI_Comm,arg2::Ptr{PetscBool})
    ccall((:PetscSharedWorkingDirectory,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscBool}),arg1,arg2)
end

function PetscGetTmp(arg1::MPI_Comm,arg2::Ptr{Uint8},size_t::Cint)
    ccall((:PetscGetTmp,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Cint),arg1,arg2,size_t)
end

function PetscFileRetrieve(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},size_t::Cint,arg4::Ptr{PetscBool})
    ccall((:PetscFileRetrieve,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Cint,Ptr{PetscBool}),arg1,arg2,arg3,size_t,arg4)
end

function PetscLs(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},size_t::Cint,arg4::Ptr{PetscBool})
    ccall((:PetscLs,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Cint,Ptr{PetscBool}),arg1,arg2,arg3,size_t,arg4)
end

function PetscOpenSocket(arg1::Ptr{Uint8},arg2::Cint,arg3::Ptr{Cint})
    ccall((:PetscOpenSocket,petsc),PetscErrorCode,(Ptr{Uint8},Cint,Ptr{Cint}),arg1,arg2,arg3)
end

#=
function PetscBinarySeek(arg1::Cint,arg2::off_t,arg3::PetscBinarySeekType,arg4::Ptr{off_t})
    ccall((:PetscBinarySeek,petsc),PetscErrorCode,(Cint,off_t,PetscBinarySeekType,Ptr{off_t}),arg1,arg2,arg3,arg4)
end

function PetscBinarySynchronizedSeek(arg1::MPI_Comm,arg2::Cint,arg3::off_t,arg4::PetscBinarySeekType,arg5::Ptr{off_t})
    ccall((:PetscBinarySynchronizedSeek,petsc),PetscErrorCode,(MPI_Comm,Cint,off_t,PetscBinarySeekType,Ptr{off_t}),arg1,arg2,arg3,arg4,arg5)
end

=#

function PetscByteSwap(arg1::Ptr{Void},arg2::PetscDataType,arg3::PetscInt)
    ccall((:PetscByteSwap,petsc),PetscErrorCode,(Ptr{Void},PetscDataType,PetscInt),arg1,arg2,arg3)
end

function PetscSetDebugTerminal(arg1::Ptr{Uint8})
    ccall((:PetscSetDebugTerminal,petsc),PetscErrorCode,(Ptr{Uint8},),arg1)
end

function PetscSetDebugger(arg1::Ptr{Uint8},arg2::PetscBool)
    ccall((:PetscSetDebugger,petsc),PetscErrorCode,(Ptr{Uint8},PetscBool),arg1,arg2)
end

function PetscSetDefaultDebugger()
    ccall((:PetscSetDefaultDebugger,petsc),PetscErrorCode,())
end

function PetscSetDebuggerFromString(arg1::Ptr{Uint8})
    ccall((:PetscSetDebuggerFromString,petsc),PetscErrorCode,(Ptr{Uint8},),arg1)
end

function PetscAttachDebugger()
    ccall((:PetscAttachDebugger,petsc),PetscErrorCode,())
end

function PetscStopForDebugger()
    ccall((:PetscStopForDebugger,petsc),PetscErrorCode,())
end

function PetscGatherNumberOfMessages(arg1::MPI_Comm,arg2::Ptr{PetscMPIInt},arg3::Ptr{PetscMPIInt},arg4::Ptr{PetscMPIInt})
    ccall((:PetscGatherNumberOfMessages,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscMPIInt},Ptr{PetscMPIInt},Ptr{PetscMPIInt}),arg1,arg2,arg3,arg4)
end

function PetscGatherMessageLengths(arg1::MPI_Comm,arg2::PetscMPIInt,arg3::PetscMPIInt,arg4::Ptr{PetscMPIInt},arg5::Ptr{Ptr{PetscMPIInt}},arg6::Ptr{Ptr{PetscMPIInt}})
    ccall((:PetscGatherMessageLengths,petsc),PetscErrorCode,(MPI_Comm,PetscMPIInt,PetscMPIInt,Ptr{PetscMPIInt},Ptr{Ptr{PetscMPIInt}},Ptr{Ptr{PetscMPIInt}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscGatherMessageLengths2(arg1::MPI_Comm,arg2::PetscMPIInt,arg3::PetscMPIInt,arg4::Ptr{PetscMPIInt},arg5::Ptr{PetscMPIInt},arg6::Ptr{Ptr{PetscMPIInt}},arg7::Ptr{Ptr{PetscMPIInt}},arg8::Ptr{Ptr{PetscMPIInt}})
    ccall((:PetscGatherMessageLengths2,petsc),PetscErrorCode,(MPI_Comm,PetscMPIInt,PetscMPIInt,Ptr{PetscMPIInt},Ptr{PetscMPIInt},Ptr{Ptr{PetscMPIInt}},Ptr{Ptr{PetscMPIInt}},Ptr{Ptr{PetscMPIInt}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

typealias MPI_Request Int32

function PetscPostIrecvInt(arg1::MPI_Comm,arg2::PetscMPIInt,arg3::PetscMPIInt,arg4::Ptr{PetscMPIInt},arg5::Ptr{PetscMPIInt},arg6::Ptr{Ptr{Ptr{PetscInt}}},arg7::Ptr{Ptr{MPI_Request}})
    ccall((:PetscPostIrecvInt,petsc),PetscErrorCode,(MPI_Comm,PetscMPIInt,PetscMPIInt,Ptr{PetscMPIInt},Ptr{PetscMPIInt},Ptr{Ptr{Ptr{PetscInt}}},Ptr{Ptr{MPI_Request}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscPostIrecvScalar(arg1::MPI_Comm,arg2::PetscMPIInt,arg3::PetscMPIInt,arg4::Ptr{PetscMPIInt},arg5::Ptr{PetscMPIInt},arg6::Ptr{Ptr{Ptr{PetscScalar}}},arg7::Ptr{Ptr{MPI_Request}})
    ccall((:PetscPostIrecvScalar,petsc),PetscErrorCode,(MPI_Comm,PetscMPIInt,PetscMPIInt,Ptr{PetscMPIInt},Ptr{PetscMPIInt},Ptr{Ptr{Ptr{PetscScalar}}},Ptr{Ptr{MPI_Request}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscCommBuildTwoSided(arg1::MPI_Comm,arg2::PetscMPIInt,arg3::MPI_Datatype,arg4::PetscInt,arg5::Ptr{PetscMPIInt},arg6::Ptr{Void},arg7::Ptr{PetscInt},arg8::Ptr{Ptr{PetscMPIInt}},arg9::Ptr{Void})
    ccall((:PetscCommBuildTwoSided,petsc),PetscErrorCode,(MPI_Comm,PetscMPIInt,MPI_Datatype,PetscInt,Ptr{PetscMPIInt},Ptr{Void},Ptr{PetscInt},Ptr{Ptr{PetscMPIInt}},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function PetscCommBuildTwoSidedSetType(arg1::MPI_Comm,arg2::PetscBuildTwoSidedType)
    ccall((:PetscCommBuildTwoSidedSetType,petsc),PetscErrorCode,(MPI_Comm,PetscBuildTwoSidedType),arg1,arg2)
end

function PetscCommBuildTwoSidedGetType(arg1::MPI_Comm,arg2::Ptr{PetscBuildTwoSidedType})
    ccall((:PetscCommBuildTwoSidedGetType,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscBuildTwoSidedType}),arg1,arg2)
end

function PetscSSEIsEnabled(arg1::MPI_Comm,arg2::Ptr{PetscBool},arg3::Ptr{PetscBool})
    ccall((:PetscSSEIsEnabled,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3)
end

function PetscObjectComm(arg1::PetscObject)
    ccall((:PetscObjectComm,petsc),MPI_Comm,(PetscObject,),arg1)
end

function PetscSubcommDestroy(arg1::Ptr{PetscSubcomm})
    ccall((:PetscSubcommDestroy,petsc),PetscErrorCode,(Ptr{PetscSubcomm},),arg1)
end

function PetscSubcommSetNumber(arg1::PetscSubcomm,arg2::PetscInt)
    ccall((:PetscSubcommSetNumber,petsc),PetscErrorCode,(PetscSubcomm,PetscInt),arg1,arg2)
end

function PetscSubcommSetType(arg1::PetscSubcomm,arg2::PetscSubcommType)
    ccall((:PetscSubcommSetType,petsc),PetscErrorCode,(PetscSubcomm,PetscSubcommType),arg1,arg2)
end

function PetscSubcommSetTypeGeneral(arg1::PetscSubcomm,arg2::PetscMPIInt,arg3::PetscMPIInt)
    ccall((:PetscSubcommSetTypeGeneral,petsc),PetscErrorCode,(PetscSubcomm,PetscMPIInt,PetscMPIInt),arg1,arg2,arg3)
end

function PetscSubcommView(arg1::PetscSubcomm,arg2::PetscViewer)
    ccall((:PetscSubcommView,petsc),PetscErrorCode,(PetscSubcomm,PetscViewer),arg1,arg2)
end

function PetscSubcommSetFromOptions(arg1::PetscSubcomm)
    ccall((:PetscSubcommSetFromOptions,petsc),PetscErrorCode,(PetscSubcomm,),arg1)
end

function PetscSegBufferCreate()
    ccall((:PetscSegBufferCreate,petsc),PetscErrorCode,())
end

function PetscSegBufferDestroy(arg1::Ptr{PetscSegBuffer})
    ccall((:PetscSegBufferDestroy,petsc),PetscErrorCode,(Ptr{PetscSegBuffer},),arg1)
end

function PetscSegBufferGet(arg1::PetscSegBuffer,size_t::Cint,arg2::Ptr{Void})
    ccall((:PetscSegBufferGet,petsc),PetscErrorCode,(PetscSegBuffer,Cint,Ptr{Void}),arg1,size_t,arg2)
end

function PetscSegBufferExtractAlloc(arg1::PetscSegBuffer,arg2::Ptr{Void})
    ccall((:PetscSegBufferExtractAlloc,petsc),PetscErrorCode,(PetscSegBuffer,Ptr{Void}),arg1,arg2)
end

function PetscSegBufferExtractTo(arg1::PetscSegBuffer,arg2::Ptr{Void})
    ccall((:PetscSegBufferExtractTo,petsc),PetscErrorCode,(PetscSegBuffer,Ptr{Void}),arg1,arg2)
end

function PetscSegBufferExtractInPlace(arg1::PetscSegBuffer,arg2::Ptr{Void})
    ccall((:PetscSegBufferExtractInPlace,petsc),PetscErrorCode,(PetscSegBuffer,Ptr{Void}),arg1,arg2)
end

function PetscSegBufferGetSize(arg1::PetscSegBuffer,arg2::Ptr{Cint})
    ccall((:PetscSegBufferGetSize,petsc),PetscErrorCode,(PetscSegBuffer,Ptr{Cint}),arg1,arg2)
end

function PetscSegBufferUnuse(arg1::PetscSegBuffer,size_t::Cint)
    ccall((:PetscSegBufferUnuse,petsc),PetscErrorCode,(PetscSegBuffer,Cint),arg1,size_t)
end

function PetscGoogleDriveAuthorize(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},size_t::Cint)
    ccall((:PetscGoogleDriveAuthorize,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Cint),arg1,arg2,arg3,size_t)
end

function PetscGoogleDriveRefresh(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},size_t::Cint)
    ccall((:PetscGoogleDriveRefresh,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Cint),arg1,arg2,arg3,size_t)
end

function PetscGoogleDriveUpload(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8})
    ccall((:PetscGoogleDriveUpload,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8}),arg1,arg2,arg3)
end

function PetscBoxAuthorize(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},size_t::Cint)
    ccall((:PetscBoxAuthorize,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Cint),arg1,arg2,arg3,size_t)
end

function PetscBoxRefresh(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},size_t::Cint)
    ccall((:PetscBoxRefresh,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Cint),arg1,arg2,arg3,arg4,size_t)
end

function PetscTextBelt(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{PetscBool})
    ccall((:PetscTextBelt,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function PetscPullJSONValue(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{Uint8},size_t::Cint,arg4::Ptr{PetscBool})
    ccall((:PetscPullJSONValue,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Cint,Ptr{PetscBool}),arg1,arg2,arg3,size_t,arg4)
end

function PetscPushJSONValue(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{Uint8},size_t::Cint)
    ccall((:PetscPushJSONValue,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Cint),arg1,arg2,arg3,size_t)
end

function PetscBagCreate(arg1::MPI_Comm,size_t::Cint,arg2::Ptr{PetscBag})
    ccall((:PetscBagCreate,petsc),PetscErrorCode,(MPI_Comm,Cint,Ptr{PetscBag}),arg1,size_t,arg2)
end

function PetscBagDestroy(arg1::Ptr{PetscBag})
    ccall((:PetscBagDestroy,petsc),PetscErrorCode,(Ptr{PetscBag},),arg1)
end

function PetscBagGetData(arg1::PetscBag,arg2::Ptr{Ptr{Void}})
    ccall((:PetscBagGetData,petsc),PetscErrorCode,(PetscBag,Ptr{Ptr{Void}}),arg1,arg2)
end

function PetscBagRegisterReal(arg1::PetscBag,arg2::Ptr{Void},PetscReal::Cint,arg3::Ptr{Uint8},arg4::Ptr{Uint8})
    ccall((:PetscBagRegisterReal,petsc),PetscErrorCode,(PetscBag,Ptr{Void},Cint,Ptr{Uint8},Ptr{Uint8}),arg1,arg2,PetscReal,arg3,arg4)
end

function PetscBagRegisterRealArray(arg1::PetscBag,arg2::Ptr{Void},arg3::PetscInt,arg4::Ptr{Uint8},arg5::Ptr{Uint8})
    ccall((:PetscBagRegisterRealArray,petsc),PetscErrorCode,(PetscBag,Ptr{Void},PetscInt,Ptr{Uint8},Ptr{Uint8}),arg1,arg2,arg3,arg4,arg5)
end

function PetscBagRegisterString(arg1::PetscBag,arg2::Ptr{Void},arg3::PetscInt,arg4::Ptr{Uint8},arg5::Ptr{Uint8},arg6::Ptr{Uint8})
    ccall((:PetscBagRegisterString,petsc),PetscErrorCode,(PetscBag,Ptr{Void},PetscInt,Ptr{Uint8},Ptr{Uint8},Ptr{Uint8}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscBagRegisterScalar(arg1::PetscBag,arg2::Ptr{Void},arg3::PetscScalar,arg4::Ptr{Uint8},arg5::Ptr{Uint8})
    ccall((:PetscBagRegisterScalar,petsc),PetscErrorCode,(PetscBag,Ptr{Void},PetscScalar,Ptr{Uint8},Ptr{Uint8}),arg1,arg2,arg3,arg4,arg5)
end

function PetscBagRegisterInt(arg1::PetscBag,arg2::Ptr{Void},arg3::PetscInt,arg4::Ptr{Uint8},arg5::Ptr{Uint8})
    ccall((:PetscBagRegisterInt,petsc),PetscErrorCode,(PetscBag,Ptr{Void},PetscInt,Ptr{Uint8},Ptr{Uint8}),arg1,arg2,arg3,arg4,arg5)
end

function PetscBagRegister64bitInt(arg1::PetscBag,arg2::Ptr{Void},arg3::Petsc64bitInt,arg4::Ptr{Uint8},arg5::Ptr{Uint8})
    ccall((:PetscBagRegister64bitInt,petsc),PetscErrorCode,(PetscBag,Ptr{Void},Petsc64bitInt,Ptr{Uint8},Ptr{Uint8}),arg1,arg2,arg3,arg4,arg5)
end

function PetscBagRegisterIntArray(arg1::PetscBag,arg2::Ptr{Void},arg3::PetscInt,arg4::Ptr{Uint8},arg5::Ptr{Uint8})
    ccall((:PetscBagRegisterIntArray,petsc),PetscErrorCode,(PetscBag,Ptr{Void},PetscInt,Ptr{Uint8},Ptr{Uint8}),arg1,arg2,arg3,arg4,arg5)
end

function PetscBagRegisterEnum(arg1::PetscBag,arg2::Ptr{Void},arg3::Ptr{Ptr{Uint8}},arg4::PetscEnum,arg5::Ptr{Uint8},arg6::Ptr{Uint8})
    ccall((:PetscBagRegisterEnum,petsc),PetscErrorCode,(PetscBag,Ptr{Void},Ptr{Ptr{Uint8}},PetscEnum,Ptr{Uint8},Ptr{Uint8}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscBagRegisterBool(arg1::PetscBag,arg2::Ptr{Void},arg3::PetscBool,arg4::Ptr{Uint8},arg5::Ptr{Uint8})
    ccall((:PetscBagRegisterBool,petsc),PetscErrorCode,(PetscBag,Ptr{Void},PetscBool,Ptr{Uint8},Ptr{Uint8}),arg1,arg2,arg3,arg4,arg5)
end

function PetscBagRegisterBoolArray(arg1::PetscBag,arg2::Ptr{Void},arg3::PetscInt,arg4::Ptr{Uint8},arg5::Ptr{Uint8})
    ccall((:PetscBagRegisterBoolArray,petsc),PetscErrorCode,(PetscBag,Ptr{Void},PetscInt,Ptr{Uint8},Ptr{Uint8}),arg1,arg2,arg3,arg4,arg5)
end

function PetscBagGetNames(arg1::PetscBag,arg2::Ptr{Ptr{Uint8}})
    ccall((:PetscBagGetNames,petsc),PetscErrorCode,(PetscBag,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function PetscBagSetFromOptions(arg1::PetscBag)
    ccall((:PetscBagSetFromOptions,petsc),PetscErrorCode,(PetscBag,),arg1)
end

function PetscBagGetName(arg1::PetscBag,arg2::Ptr{Ptr{Uint8}})
    ccall((:PetscBagGetName,petsc),PetscErrorCode,(PetscBag,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function PetscBagSetName(arg1::PetscBag,arg2::Ptr{Uint8},arg3::Ptr{Uint8})
    ccall((:PetscBagSetName,petsc),PetscErrorCode,(PetscBag,Ptr{Uint8},Ptr{Uint8}),arg1,arg2,arg3)
end

function PetscBagSetOptionsPrefix(arg1::PetscBag,arg2::Ptr{Uint8})
    ccall((:PetscBagSetOptionsPrefix,petsc),PetscErrorCode,(PetscBag,Ptr{Uint8}),arg1,arg2)
end

function PetscBagView(arg1::PetscBag,arg2::PetscViewer)
    ccall((:PetscBagView,petsc),PetscErrorCode,(PetscBag,PetscViewer),arg1,arg2)
end

function PetscBagLoad(arg1::PetscViewer,arg2::PetscBag)
    ccall((:PetscBagLoad,petsc),PetscErrorCode,(PetscViewer,PetscBag),arg1,arg2)
end

function PetscBagSetViewer(arg1::PetscBag,arg2::Ptr{Void})
    ccall((:PetscBagSetViewer,petsc),PetscErrorCode,(PetscBag,Ptr{Void}),arg1,arg2)
end

function PetscBagSetLoader(arg1::PetscBag,arg2::Ptr{Void})
    ccall((:PetscBagSetLoader,petsc),PetscErrorCode,(PetscBag,Ptr{Void}),arg1,arg2)
end

function PetscBagSetDestroy(arg1::PetscBag,arg2::Ptr{Void})
    ccall((:PetscBagSetDestroy,petsc),PetscErrorCode,(PetscBag,Ptr{Void}),arg1,arg2)
end

function PetscGetCPUTime(arg1::Ptr{PetscLogDouble})
    ccall((:PetscGetCPUTime,petsc),PetscErrorCode,(Ptr{PetscLogDouble},),arg1)
end

function PetscViewerInitializePackage()
    ccall((:PetscViewerInitializePackage,petsc),PetscErrorCode,())
end

function PetscViewerRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:PetscViewerRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function PetscViewerCreate(arg1::MPI_Comm,arg2::Ptr{PetscViewer})
    ccall((:PetscViewerCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscViewer}),arg1,arg2)
end

function PetscViewerSetFromOptions(arg1::PetscViewer)
    ccall((:PetscViewerSetFromOptions,petsc),PetscErrorCode,(PetscViewer,),arg1)
end

function PetscViewerASCIIOpenWithVoid(arg1::MPI_Comm,arg2::Ptr{Void},arg3::Ptr{PetscViewer})
    ccall((:PetscViewerASCIIOpenWithVoid,petsc),PetscErrorCode,(MPI_Comm,Ptr{Void},Ptr{PetscViewer}),arg1,arg2,arg3)
end

function PetscViewerASCIIOpen(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{PetscViewer})
    ccall((:PetscViewerASCIIOpen,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{PetscViewer}),arg1,arg2,arg3)
end

function PetscViewerASCIISetVoid(arg1::PetscViewer,arg2::Ptr{Void})
    ccall((:PetscViewerASCIISetVoid,petsc),PetscErrorCode,(PetscViewer,Ptr{Void}),arg1,arg2)
end

function PetscViewerBinaryOpen(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::PetscFileMode,arg4::Ptr{PetscViewer})
    ccall((:PetscViewerBinaryOpen,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},PetscFileMode,Ptr{PetscViewer}),arg1,arg2,arg3,arg4)
end

function PetscViewerBinaryGetFlowControl(arg1::PetscViewer,arg2::Ptr{PetscInt})
    ccall((:PetscViewerBinaryGetFlowControl,petsc),PetscErrorCode,(PetscViewer,Ptr{PetscInt}),arg1,arg2)
end

function PetscViewerBinarySetFlowControl(arg1::PetscViewer,arg2::PetscInt)
    ccall((:PetscViewerBinarySetFlowControl,petsc),PetscErrorCode,(PetscViewer,PetscInt),arg1,arg2)
end

function PetscViewerBinarySetUseMPIIO(arg1::PetscViewer,arg2::PetscBool)
    ccall((:PetscViewerBinarySetUseMPIIO,petsc),PetscErrorCode,(PetscViewer,PetscBool),arg1,arg2)
end

function PetscViewerBinaryGetUseMPIIO(arg1::PetscViewer,arg2::Ptr{PetscBool})
    ccall((:PetscViewerBinaryGetUseMPIIO,petsc),PetscErrorCode,(PetscViewer,Ptr{PetscBool}),arg1,arg2)
end

function PetscViewerSocketOpen(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Cint,arg4::Ptr{PetscViewer})
    ccall((:PetscViewerSocketOpen,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Cint,Ptr{PetscViewer}),arg1,arg2,arg3,arg4)
end

function PetscViewerStringOpen(arg1::MPI_Comm,arg2::Ptr{Uint8},size_t::Cint,arg3::Ptr{PetscViewer})
    ccall((:PetscViewerStringOpen,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Cint,Ptr{PetscViewer}),arg1,arg2,size_t,arg3)
end

function PetscViewerDrawOpen(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Cint,arg5::Cint,arg6::Cint,arg7::Cint,arg8::Ptr{PetscViewer})
    ccall((:PetscViewerDrawOpen,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Cint,Cint,Cint,Cint,Ptr{PetscViewer}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscViewerDrawSetDrawType(arg1::PetscViewer,arg2::PetscDrawType)
    ccall((:PetscViewerDrawSetDrawType,petsc),PetscErrorCode,(PetscViewer,PetscDrawType),arg1,arg2)
end

function PetscViewerMathematicaOpen(arg1::MPI_Comm,arg2::Cint,arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Ptr{PetscViewer})
    ccall((:PetscViewerMathematicaOpen,petsc),PetscErrorCode,(MPI_Comm,Cint,Ptr{Uint8},Ptr{Uint8},Ptr{PetscViewer}),arg1,arg2,arg3,arg4,arg5)
end

function PetscViewerSiloOpen(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{PetscViewer})
    ccall((:PetscViewerSiloOpen,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{PetscViewer}),arg1,arg2,arg3)
end

function PetscViewerMatlabOpen(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::PetscFileMode,arg4::Ptr{PetscViewer})
    ccall((:PetscViewerMatlabOpen,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},PetscFileMode,Ptr{PetscViewer}),arg1,arg2,arg3,arg4)
end

function PetscViewerGetType(arg1::PetscViewer,arg2::Ptr{PetscViewerType})
    ccall((:PetscViewerGetType,petsc),PetscErrorCode,(PetscViewer,Ptr{PetscViewerType}),arg1,arg2)
end

function PetscViewerSetType(arg1::PetscViewer,arg2::PetscViewerType)
    ccall((:PetscViewerSetType,petsc),PetscErrorCode,(PetscViewer,PetscViewerType),arg1,arg2)
end

function PetscViewerDestroy(arg1::Ptr{PetscViewer})
    ccall((:PetscViewerDestroy,petsc),PetscErrorCode,(Ptr{PetscViewer},),arg1)
end

function PetscViewerGetSingleton(arg1::PetscViewer,arg2::Ptr{PetscViewer})
    ccall((:PetscViewerGetSingleton,petsc),PetscErrorCode,(PetscViewer,Ptr{PetscViewer}),arg1,arg2)
end

function PetscViewerRestoreSingleton(arg1::PetscViewer,arg2::Ptr{PetscViewer})
    ccall((:PetscViewerRestoreSingleton,petsc),PetscErrorCode,(PetscViewer,Ptr{PetscViewer}),arg1,arg2)
end

function PetscViewerGetSubcomm(arg1::PetscViewer,arg2::MPI_Comm,arg3::Ptr{PetscViewer})
    ccall((:PetscViewerGetSubcomm,petsc),PetscErrorCode,(PetscViewer,MPI_Comm,Ptr{PetscViewer}),arg1,arg2,arg3)
end

function PetscViewerRestoreSubcomm(arg1::PetscViewer,arg2::MPI_Comm,arg3::Ptr{PetscViewer})
    ccall((:PetscViewerRestoreSubcomm,petsc),PetscErrorCode,(PetscViewer,MPI_Comm,Ptr{PetscViewer}),arg1,arg2,arg3)
end

function PetscViewerSetUp(arg1::PetscViewer)
    ccall((:PetscViewerSetUp,petsc),PetscErrorCode,(PetscViewer,),arg1)
end

function PetscViewerView(arg1::PetscViewer,arg2::PetscViewer)
    ccall((:PetscViewerView,petsc),PetscErrorCode,(PetscViewer,PetscViewer),arg1,arg2)
end

function PetscViewerAppendOptionsPrefix(arg1::PetscViewer,arg2::Ptr{Uint8})
    ccall((:PetscViewerAppendOptionsPrefix,petsc),PetscErrorCode,(PetscViewer,Ptr{Uint8}),arg1,arg2)
end

function PetscViewerGetOptionsPrefix(arg1::PetscViewer,arg2::Ptr{Ptr{Uint8}})
    ccall((:PetscViewerGetOptionsPrefix,petsc),PetscErrorCode,(PetscViewer,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function PetscViewerSetFormat(arg1::PetscViewer,arg2::PetscViewerFormat)
    ccall((:PetscViewerSetFormat,petsc),PetscErrorCode,(PetscViewer,PetscViewerFormat),arg1,arg2)
end

function PetscViewerPushFormat(arg1::PetscViewer,arg2::PetscViewerFormat)
    ccall((:PetscViewerPushFormat,petsc),PetscErrorCode,(PetscViewer,PetscViewerFormat),arg1,arg2)
end

function PetscViewerPopFormat(arg1::PetscViewer)
    ccall((:PetscViewerPopFormat,petsc),PetscErrorCode,(PetscViewer,),arg1)
end

function PetscViewerGetFormat(arg1::PetscViewer,arg2::Ptr{PetscViewerFormat})
    ccall((:PetscViewerGetFormat,petsc),PetscErrorCode,(PetscViewer,Ptr{PetscViewerFormat}),arg1,arg2)
end

function PetscViewerFlush(arg1::PetscViewer)
    ccall((:PetscViewerFlush,petsc),PetscErrorCode,(PetscViewer,),arg1)
end

function PetscOptionsGetViewer(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{PetscViewer},arg5::Ptr{PetscViewerFormat},arg6::Ptr{PetscBool})
    ccall((:PetscOptionsGetViewer,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Ptr{PetscViewer},Ptr{PetscViewerFormat},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscOptionsViewer_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Ptr{PetscViewer},arg6::Ptr{PetscViewerFormat},arg7::Ptr{PetscBool})
    ccall((:PetscOptionsViewer_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{PetscViewer},Ptr{PetscViewerFormat},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscViewerASCIIGetPointer(arg1::PetscViewer,arg2::Ptr{Ptr{Void}})
    ccall((:PetscViewerASCIIGetPointer,petsc),PetscErrorCode,(PetscViewer,Ptr{Ptr{Void}}),arg1,arg2)
end

function PetscViewerFileGetMode(arg1::PetscViewer,arg2::Ptr{PetscFileMode})
    ccall((:PetscViewerFileGetMode,petsc),PetscErrorCode,(PetscViewer,Ptr{PetscFileMode}),arg1,arg2)
end

function PetscViewerFileSetMode(arg1::PetscViewer,arg2::PetscFileMode)
    ccall((:PetscViewerFileSetMode,petsc),PetscErrorCode,(PetscViewer,PetscFileMode),arg1,arg2)
end

function PetscViewerRead(arg1::PetscViewer,arg2::Ptr{Void},arg3::PetscInt,arg4::Ptr{PetscInt},arg5::PetscDataType)
    ccall((:PetscViewerRead,petsc),PetscErrorCode,(PetscViewer,Ptr{Void},PetscInt,Ptr{PetscInt},PetscDataType),arg1,arg2,arg3,arg4,arg5)
end

function PetscViewerASCIISynchronizedAllow(arg1::PetscViewer,arg2::PetscBool)
    ccall((:PetscViewerASCIISynchronizedAllow,petsc),PetscErrorCode,(PetscViewer,PetscBool),arg1,arg2)
end

function PetscViewerASCIIPushTab(arg1::PetscViewer)
    ccall((:PetscViewerASCIIPushTab,petsc),PetscErrorCode,(PetscViewer,),arg1)
end

function PetscViewerASCIIPopTab(arg1::PetscViewer)
    ccall((:PetscViewerASCIIPopTab,petsc),PetscErrorCode,(PetscViewer,),arg1)
end

function PetscViewerASCIIUseTabs(arg1::PetscViewer,arg2::PetscBool)
    ccall((:PetscViewerASCIIUseTabs,petsc),PetscErrorCode,(PetscViewer,PetscBool),arg1,arg2)
end

function PetscViewerASCIISetTab(arg1::PetscViewer,arg2::PetscInt)
    ccall((:PetscViewerASCIISetTab,petsc),PetscErrorCode,(PetscViewer,PetscInt),arg1,arg2)
end

function PetscViewerASCIIGetTab(arg1::PetscViewer,arg2::Ptr{PetscInt})
    ccall((:PetscViewerASCIIGetTab,petsc),PetscErrorCode,(PetscViewer,Ptr{PetscInt}),arg1,arg2)
end

function PetscViewerASCIIAddTab(arg1::PetscViewer,arg2::PetscInt)
    ccall((:PetscViewerASCIIAddTab,petsc),PetscErrorCode,(PetscViewer,PetscInt),arg1,arg2)
end

function PetscViewerASCIISubtractTab(arg1::PetscViewer,arg2::PetscInt)
    ccall((:PetscViewerASCIISubtractTab,petsc),PetscErrorCode,(PetscViewer,PetscInt),arg1,arg2)
end

function PetscViewerASCIIRead(arg1::PetscViewer,arg2::Ptr{Void},arg3::PetscInt,arg4::Ptr{PetscInt},arg5::PetscDataType)
    ccall((:PetscViewerASCIIRead,petsc),PetscErrorCode,(PetscViewer,Ptr{Void},PetscInt,Ptr{PetscInt},PetscDataType),arg1,arg2,arg3,arg4,arg5)
end

function PetscViewerBinaryGetDescriptor(arg1::PetscViewer,arg2::Ptr{Cint})
    ccall((:PetscViewerBinaryGetDescriptor,petsc),PetscErrorCode,(PetscViewer,Ptr{Cint}),arg1,arg2)
end

function PetscViewerBinaryGetInfoPointer(arg1::PetscViewer,arg2::Ptr{Ptr{Void}})
    ccall((:PetscViewerBinaryGetInfoPointer,petsc),PetscErrorCode,(PetscViewer,Ptr{Ptr{Void}}),arg1,arg2)
end

function PetscViewerBinaryRead(arg1::PetscViewer,arg2::Ptr{Void},arg3::PetscInt,arg4::Ptr{PetscInt},arg5::PetscDataType)
    ccall((:PetscViewerBinaryRead,petsc),PetscErrorCode,(PetscViewer,Ptr{Void},PetscInt,Ptr{PetscInt},PetscDataType),arg1,arg2,arg3,arg4,arg5)
end

function PetscViewerBinaryWrite(arg1::PetscViewer,arg2::Ptr{Void},arg3::PetscInt,arg4::PetscDataType,arg5::PetscBool)
    ccall((:PetscViewerBinaryWrite,petsc),PetscErrorCode,(PetscViewer,Ptr{Void},PetscInt,PetscDataType,PetscBool),arg1,arg2,arg3,arg4,arg5)
end

function PetscViewerStringSetString(arg1::PetscViewer,arg2::Ptr{Uint8},arg3::PetscInt)
    ccall((:PetscViewerStringSetString,petsc),PetscErrorCode,(PetscViewer,Ptr{Uint8},PetscInt),arg1,arg2,arg3)
end

function PetscViewerDrawClear(arg1::PetscViewer)
    ccall((:PetscViewerDrawClear,petsc),PetscErrorCode,(PetscViewer,),arg1)
end

function PetscViewerDrawSetHold(arg1::PetscViewer,arg2::PetscBool)
    ccall((:PetscViewerDrawSetHold,petsc),PetscErrorCode,(PetscViewer,PetscBool),arg1,arg2)
end

function PetscViewerDrawGetHold(arg1::PetscViewer,arg2::Ptr{PetscBool})
    ccall((:PetscViewerDrawGetHold,petsc),PetscErrorCode,(PetscViewer,Ptr{PetscBool}),arg1,arg2)
end

function PetscViewerDrawSetPause(arg1::PetscViewer,PetscReal::Cint)
    ccall((:PetscViewerDrawSetPause,petsc),PetscErrorCode,(PetscViewer,Cint),arg1,PetscReal)
end

function PetscViewerDrawGetPause(arg1::PetscViewer,arg2::Ptr{Cint})
    ccall((:PetscViewerDrawGetPause,petsc),PetscErrorCode,(PetscViewer,Ptr{Cint}),arg1,arg2)
end

function PetscViewerDrawSetInfo(arg1::PetscViewer,arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Cint,arg5::Cint,arg6::Cint,arg7::Cint)
    ccall((:PetscViewerDrawSetInfo,petsc),PetscErrorCode,(PetscViewer,Ptr{Uint8},Ptr{Uint8},Cint,Cint,Cint,Cint),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscViewerDrawResize(arg1::PetscViewer,arg2::Cint,arg3::Cint)
    ccall((:PetscViewerDrawResize,petsc),PetscErrorCode,(PetscViewer,Cint,Cint),arg1,arg2,arg3)
end

function PetscViewerDrawSetBounds(arg1::PetscViewer,arg2::PetscInt,arg3::Ptr{Cint})
    ccall((:PetscViewerDrawSetBounds,petsc),PetscErrorCode,(PetscViewer,PetscInt,Ptr{Cint}),arg1,arg2,arg3)
end

function PetscViewerDrawGetBounds(arg1::PetscViewer,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{Cint}})
    ccall((:PetscViewerDrawGetBounds,petsc),PetscErrorCode,(PetscViewer,Ptr{PetscInt},Ptr{Ptr{Cint}}),arg1,arg2,arg3)
end

function PetscViewerSocketSetConnection(arg1::PetscViewer,arg2::Ptr{Uint8},arg3::Cint)
    ccall((:PetscViewerSocketSetConnection,petsc),PetscErrorCode,(PetscViewer,Ptr{Uint8},Cint),arg1,arg2,arg3)
end

function PetscViewerBinarySkipInfo(arg1::PetscViewer)
    ccall((:PetscViewerBinarySkipInfo,petsc),PetscErrorCode,(PetscViewer,),arg1)
end

function PetscViewerBinarySetSkipInfo(arg1::PetscViewer,arg2::PetscBool)
    ccall((:PetscViewerBinarySetSkipInfo,petsc),PetscErrorCode,(PetscViewer,PetscBool),arg1,arg2)
end

function PetscViewerBinaryGetSkipInfo(arg1::PetscViewer,arg2::Ptr{PetscBool})
    ccall((:PetscViewerBinaryGetSkipInfo,petsc),PetscErrorCode,(PetscViewer,Ptr{PetscBool}),arg1,arg2)
end

function PetscViewerBinarySetSkipOptions(arg1::PetscViewer,arg2::PetscBool)
    ccall((:PetscViewerBinarySetSkipOptions,petsc),PetscErrorCode,(PetscViewer,PetscBool),arg1,arg2)
end

function PetscViewerBinaryGetSkipOptions(arg1::PetscViewer,arg2::Ptr{PetscBool})
    ccall((:PetscViewerBinaryGetSkipOptions,petsc),PetscErrorCode,(PetscViewer,Ptr{PetscBool}),arg1,arg2)
end

function PetscViewerBinarySetSkipHeader(arg1::PetscViewer,arg2::PetscBool)
    ccall((:PetscViewerBinarySetSkipHeader,petsc),PetscErrorCode,(PetscViewer,PetscBool),arg1,arg2)
end

function PetscViewerBinaryGetSkipHeader(arg1::PetscViewer,arg2::Ptr{PetscBool})
    ccall((:PetscViewerBinaryGetSkipHeader,petsc),PetscErrorCode,(PetscViewer,Ptr{PetscBool}),arg1,arg2)
end

function PetscViewerBinaryReadStringArray(arg1::PetscViewer,arg2::Ptr{Ptr{Ptr{Uint8}}})
    ccall((:PetscViewerBinaryReadStringArray,petsc),PetscErrorCode,(PetscViewer,Ptr{Ptr{Ptr{Uint8}}}),arg1,arg2)
end

function PetscViewerBinaryWriteStringArray(arg1::PetscViewer,arg2::Ptr{Ptr{Uint8}})
    ccall((:PetscViewerBinaryWriteStringArray,petsc),PetscErrorCode,(PetscViewer,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function PetscViewerFileSetName(arg1::PetscViewer,arg2::Ptr{Uint8})
    ccall((:PetscViewerFileSetName,petsc),PetscErrorCode,(PetscViewer,Ptr{Uint8}),arg1,arg2)
end

function PetscViewerFileGetName(arg1::PetscViewer,arg2::Ptr{Ptr{Uint8}})
    ccall((:PetscViewerFileGetName,petsc),PetscErrorCode,(PetscViewer,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function PetscViewerVUGetPointer(arg1::PetscViewer,arg2::Ptr{Ptr{Void}})
    ccall((:PetscViewerVUGetPointer,petsc),PetscErrorCode,(PetscViewer,Ptr{Ptr{Void}}),arg1,arg2)
end

function PetscViewerVUSetVecSeen(arg1::PetscViewer,arg2::PetscBool)
    ccall((:PetscViewerVUSetVecSeen,petsc),PetscErrorCode,(PetscViewer,PetscBool),arg1,arg2)
end

function PetscViewerVUGetVecSeen(arg1::PetscViewer,arg2::Ptr{PetscBool})
    ccall((:PetscViewerVUGetVecSeen,petsc),PetscErrorCode,(PetscViewer,Ptr{PetscBool}),arg1,arg2)
end

function PetscViewerVUFlushDeferred(arg1::PetscViewer)
    ccall((:PetscViewerVUFlushDeferred,petsc),PetscErrorCode,(PetscViewer,),arg1)
end

function PetscViewerMathematicaInitializePackage()
    ccall((:PetscViewerMathematicaInitializePackage,petsc),PetscErrorCode,())
end

function PetscViewerMathematicaFinalizePackage()
    ccall((:PetscViewerMathematicaFinalizePackage,petsc),PetscErrorCode,())
end

function PetscViewerMathematicaGetName(arg1::PetscViewer,arg2::Ptr{Ptr{Uint8}})
    ccall((:PetscViewerMathematicaGetName,petsc),PetscErrorCode,(PetscViewer,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function PetscViewerMathematicaSetName(arg1::PetscViewer,arg2::Ptr{Uint8})
    ccall((:PetscViewerMathematicaSetName,petsc),PetscErrorCode,(PetscViewer,Ptr{Uint8}),arg1,arg2)
end

function PetscViewerMathematicaClearName(arg1::PetscViewer)
    ccall((:PetscViewerMathematicaClearName,petsc),PetscErrorCode,(PetscViewer,),arg1)
end

function PetscViewerMathematicaSkipPackets(arg1::PetscViewer,arg2::Cint)
    ccall((:PetscViewerMathematicaSkipPackets,petsc),PetscErrorCode,(PetscViewer,Cint),arg1,arg2)
end

function PetscViewerSiloGetName(arg1::PetscViewer,arg2::Ptr{Ptr{Uint8}})
    ccall((:PetscViewerSiloGetName,petsc),PetscErrorCode,(PetscViewer,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function PetscViewerSiloSetName(arg1::PetscViewer,arg2::Ptr{Uint8})
    ccall((:PetscViewerSiloSetName,petsc),PetscErrorCode,(PetscViewer,Ptr{Uint8}),arg1,arg2)
end

function PetscViewerSiloClearName(arg1::PetscViewer)
    ccall((:PetscViewerSiloClearName,petsc),PetscErrorCode,(PetscViewer,),arg1)
end

function PetscViewerSiloGetMeshName(arg1::PetscViewer,arg2::Ptr{Ptr{Uint8}})
    ccall((:PetscViewerSiloGetMeshName,petsc),PetscErrorCode,(PetscViewer,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function PetscViewerSiloSetMeshName(arg1::PetscViewer,arg2::Ptr{Uint8})
    ccall((:PetscViewerSiloSetMeshName,petsc),PetscErrorCode,(PetscViewer,Ptr{Uint8}),arg1,arg2)
end

function PetscViewerSiloClearMeshName(arg1::PetscViewer)
    ccall((:PetscViewerSiloClearMeshName,petsc),PetscErrorCode,(PetscViewer,),arg1)
end

function PetscViewerNetcdfOpen(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::PetscFileMode,arg4::Ptr{PetscViewer})
    ccall((:PetscViewerNetcdfOpen,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},PetscFileMode,Ptr{PetscViewer}),arg1,arg2,arg3,arg4)
end

function PetscViewerNetcdfGetID(arg1::PetscViewer,arg2::Ptr{Cint})
    ccall((:PetscViewerNetcdfGetID,petsc),PetscErrorCode,(PetscViewer,Ptr{Cint}),arg1,arg2)
end

function PetscViewerVTKAddField(arg1::PetscViewer,arg2::PetscObject,PetscViewerVTKWriteFunction::Ptr{Void},arg3::PetscViewerVTKFieldType,arg4::PetscObject)
    ccall((:PetscViewerVTKAddField,petsc),PetscErrorCode,(PetscViewer,PetscObject,Ptr{Void},PetscViewerVTKFieldType,PetscObject),arg1,arg2,PetscViewerVTKWriteFunction,arg3,arg4)
end

function PetscViewerVTKOpen(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::PetscFileMode,arg4::Ptr{PetscViewer})
    ccall((:PetscViewerVTKOpen,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},PetscFileMode,Ptr{PetscViewer}),arg1,arg2,arg3,arg4)
end

function PETSC_VIEWER_STDOUT_(arg1::MPI_Comm)
    ccall((:PETSC_VIEWER_STDOUT_,petsc),PetscViewer,(MPI_Comm,),arg1)
end

function PetscViewerASCIIGetStdout(arg1::MPI_Comm,arg2::Ptr{PetscViewer})
    ccall((:PetscViewerASCIIGetStdout,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscViewer}),arg1,arg2)
end

function PETSC_VIEWER_STDERR_(arg1::MPI_Comm)
    ccall((:PETSC_VIEWER_STDERR_,petsc),PetscViewer,(MPI_Comm,),arg1)
end

function PetscViewerASCIIGetStderr(arg1::MPI_Comm,arg2::Ptr{PetscViewer})
    ccall((:PetscViewerASCIIGetStderr,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscViewer}),arg1,arg2)
end

function PETSC_VIEWER_DRAW_(arg1::MPI_Comm)
    ccall((:PETSC_VIEWER_DRAW_,petsc),PetscViewer,(MPI_Comm,),arg1)
end

function PETSC_VIEWER_SOCKET_(arg1::MPI_Comm)
    ccall((:PETSC_VIEWER_SOCKET_,petsc),PetscViewer,(MPI_Comm,),arg1)
end

function PETSC_VIEWER_BINARY_(arg1::MPI_Comm)
    ccall((:PETSC_VIEWER_BINARY_,petsc),PetscViewer,(MPI_Comm,),arg1)
end

function PETSC_VIEWER_MATLAB_(arg1::MPI_Comm)
    ccall((:PETSC_VIEWER_MATLAB_,petsc),PetscViewer,(MPI_Comm,),arg1)
end

function PETSC_VIEWER_HDF5_(arg1::MPI_Comm)
    ccall((:PETSC_VIEWER_HDF5_,petsc),PetscViewer,(MPI_Comm,),arg1)
end

function PetscViewerMatlabGetArray(arg1::PetscViewer,arg2::Cint,arg3::Cint,arg4::Ptr{PetscScalar},arg5::Ptr{Uint8})
    ccall((:PetscViewerMatlabGetArray,petsc),PetscErrorCode,(PetscViewer,Cint,Cint,Ptr{PetscScalar},Ptr{Uint8}),arg1,arg2,arg3,arg4,arg5)
end

function PetscViewerMatlabPutVariable(arg1::PetscViewer,arg2::Ptr{Uint8},arg3::Ptr{Void})
    ccall((:PetscViewerMatlabPutVariable,petsc),PetscErrorCode,(PetscViewer,Ptr{Uint8},Ptr{Void}),arg1,arg2,arg3)
end

function PetscViewersCreate(arg1::MPI_Comm,arg2::Ptr{PetscViewers})
    ccall((:PetscViewersCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscViewers}),arg1,arg2)
end

function PetscViewersDestroy(arg1::Ptr{PetscViewers})
    ccall((:PetscViewersDestroy,petsc),PetscErrorCode,(Ptr{PetscViewers},),arg1)
end

function PetscViewersGetViewer(arg1::PetscViewers,arg2::PetscInt,arg3::Ptr{PetscViewer})
    ccall((:PetscViewersGetViewer,petsc),PetscErrorCode,(PetscViewers,PetscInt,Ptr{PetscViewer}),arg1,arg2,arg3)
end

function PetscTableCreate(arg1::PetscInt,arg2::PetscInt,arg3::Ptr{PetscTable})
    ccall((:PetscTableCreate,petsc),PetscErrorCode,(PetscInt,PetscInt,Ptr{PetscTable}),arg1,arg2,arg3)
end

function PetscTableCreateCopy(arg1::PetscTable,arg2::Ptr{PetscTable})
    ccall((:PetscTableCreateCopy,petsc),PetscErrorCode,(PetscTable,Ptr{PetscTable}),arg1,arg2)
end

function PetscTableDestroy(arg1::Ptr{PetscTable})
    ccall((:PetscTableDestroy,petsc),PetscErrorCode,(Ptr{PetscTable},),arg1)
end

function PetscTableGetCount(arg1::PetscTable,arg2::Ptr{PetscInt})
    ccall((:PetscTableGetCount,petsc),PetscErrorCode,(PetscTable,Ptr{PetscInt}),arg1,arg2)
end

function PetscTableIsEmpty(arg1::PetscTable,arg2::Ptr{PetscInt})
    ccall((:PetscTableIsEmpty,petsc),PetscErrorCode,(PetscTable,Ptr{PetscInt}),arg1,arg2)
end

function PetscTableAddExpand(arg1::PetscTable,arg2::PetscInt,arg3::PetscInt,arg4::InsertMode)
    ccall((:PetscTableAddExpand,petsc),PetscErrorCode,(PetscTable,PetscInt,PetscInt,InsertMode),arg1,arg2,arg3,arg4)
end

function PetscTableAddCountExpand(arg1::PetscTable,arg2::PetscInt)
    ccall((:PetscTableAddCountExpand,petsc),PetscErrorCode,(PetscTable,PetscInt),arg1,arg2)
end

function PetscTableGetHeadPosition(arg1::PetscTable,arg2::Ptr{PetscTablePosition})
    ccall((:PetscTableGetHeadPosition,petsc),PetscErrorCode,(PetscTable,Ptr{PetscTablePosition}),arg1,arg2)
end

function PetscTableGetNext(arg1::PetscTable,arg2::Ptr{PetscTablePosition},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:PetscTableGetNext,petsc),PetscErrorCode,(PetscTable,Ptr{PetscTablePosition},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function PetscTableRemoveAll(arg1::PetscTable)
    ccall((:PetscTableRemoveAll,petsc),PetscErrorCode,(PetscTable,),arg1)
end

function PetscMatlabEngineCreate(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{PetscMatlabEngine})
    ccall((:PetscMatlabEngineCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{PetscMatlabEngine}),arg1,arg2,arg3)
end

function PetscMatlabEngineDestroy(arg1::Ptr{PetscMatlabEngine})
    ccall((:PetscMatlabEngineDestroy,petsc),PetscErrorCode,(Ptr{PetscMatlabEngine},),arg1)
end

function PetscMatlabEngineGetOutput(arg1::PetscMatlabEngine,arg2::Ptr{Ptr{Uint8}})
    ccall((:PetscMatlabEngineGetOutput,petsc),PetscErrorCode,(PetscMatlabEngine,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function PetscMatlabEnginePrintOutput(arg1::PetscMatlabEngine,arg2::Ptr{Void})
    ccall((:PetscMatlabEnginePrintOutput,petsc),PetscErrorCode,(PetscMatlabEngine,Ptr{Void}),arg1,arg2)
end

function PetscMatlabEnginePut(arg1::PetscMatlabEngine,arg2::PetscObject)
    ccall((:PetscMatlabEnginePut,petsc),PetscErrorCode,(PetscMatlabEngine,PetscObject),arg1,arg2)
end

function PetscMatlabEngineGet(arg1::PetscMatlabEngine,arg2::PetscObject)
    ccall((:PetscMatlabEngineGet,petsc),PetscErrorCode,(PetscMatlabEngine,PetscObject),arg1,arg2)
end

function PetscMatlabEnginePutArray(arg1::PetscMatlabEngine,arg2::Cint,arg3::Cint,arg4::Ptr{PetscScalar},arg5::Ptr{Uint8})
    ccall((:PetscMatlabEnginePutArray,petsc),PetscErrorCode,(PetscMatlabEngine,Cint,Cint,Ptr{PetscScalar},Ptr{Uint8}),arg1,arg2,arg3,arg4,arg5)
end

function PetscMatlabEngineGetArray(arg1::PetscMatlabEngine,arg2::Cint,arg3::Cint,arg4::Ptr{PetscScalar},arg5::Ptr{Uint8})
    ccall((:PetscMatlabEngineGetArray,petsc),PetscErrorCode,(PetscMatlabEngine,Cint,Cint,Ptr{PetscScalar},Ptr{Uint8}),arg1,arg2,arg3,arg4,arg5)
end

function PETSC_MATLAB_ENGINE_(arg1::MPI_Comm)
    ccall((:PETSC_MATLAB_ENGINE_,petsc),PetscMatlabEngine,(MPI_Comm,),arg1)
end

function PetscDrawInitializePackage()
    ccall((:PetscDrawInitializePackage,petsc),PetscErrorCode,())
end

function PetscDrawRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:PetscDrawRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function PetscDrawGetType(arg1::PetscDraw,arg2::Ptr{PetscDrawType})
    ccall((:PetscDrawGetType,petsc),PetscErrorCode,(PetscDraw,Ptr{PetscDrawType}),arg1,arg2)
end

function PetscDrawSetType(arg1::PetscDraw,arg2::PetscDrawType)
    ccall((:PetscDrawSetType,petsc),PetscErrorCode,(PetscDraw,PetscDrawType),arg1,arg2)
end

function PetscDrawCreate(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Cint,arg5::Cint,arg6::Cint,arg7::Cint,arg8::Ptr{PetscDraw})
    ccall((:PetscDrawCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Cint,Cint,Cint,Cint,Ptr{PetscDraw}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscDrawSetFromOptions(arg1::PetscDraw)
    ccall((:PetscDrawSetFromOptions,petsc),PetscErrorCode,(PetscDraw,),arg1)
end

function PetscDrawSetSave(arg1::PetscDraw,arg2::Ptr{Uint8},arg3::PetscBool)
    ccall((:PetscDrawSetSave,petsc),PetscErrorCode,(PetscDraw,Ptr{Uint8},PetscBool),arg1,arg2,arg3)
end

function PetscDrawSetSaveFinalImage(arg1::PetscDraw,arg2::Ptr{Uint8})
    ccall((:PetscDrawSetSaveFinalImage,petsc),PetscErrorCode,(PetscDraw,Ptr{Uint8}),arg1,arg2)
end

function PetscDrawView(arg1::PetscDraw,arg2::PetscViewer)
    ccall((:PetscDrawView,petsc),PetscErrorCode,(PetscDraw,PetscViewer),arg1,arg2)
end

function PetscDrawOpenGLUT(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Cint,arg5::Cint,arg6::Cint,arg7::Cint,arg8::Ptr{PetscDraw})
    ccall((:PetscDrawOpenGLUT,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Cint,Cint,Cint,Cint,Ptr{PetscDraw}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscDrawOpenNull(arg1::MPI_Comm,arg2::Ptr{PetscDraw})
    ccall((:PetscDrawOpenNull,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscDraw}),arg1,arg2)
end

function PetscDrawDestroy(arg1::Ptr{PetscDraw})
    ccall((:PetscDrawDestroy,petsc),PetscErrorCode,(Ptr{PetscDraw},),arg1)
end

function PetscDrawIsNull(arg1::PetscDraw,arg2::Ptr{PetscBool})
    ccall((:PetscDrawIsNull,petsc),PetscErrorCode,(PetscDraw,Ptr{PetscBool}),arg1,arg2)
end

function PetscDrawGetPopup(arg1::PetscDraw,arg2::Ptr{PetscDraw})
    ccall((:PetscDrawGetPopup,petsc),PetscErrorCode,(PetscDraw,Ptr{PetscDraw}),arg1,arg2)
end

function PetscDrawCheckResizedWindow(arg1::PetscDraw)
    ccall((:PetscDrawCheckResizedWindow,petsc),PetscErrorCode,(PetscDraw,),arg1)
end

function PetscDrawResizeWindow(arg1::PetscDraw,arg2::Cint,arg3::Cint)
    ccall((:PetscDrawResizeWindow,petsc),PetscErrorCode,(PetscDraw,Cint,Cint),arg1,arg2,arg3)
end

function PetscDrawScalePopup(arg1::PetscDraw,PetscReal::Cint,arg2::Cint)
    ccall((:PetscDrawScalePopup,petsc),PetscErrorCode,(PetscDraw,Cint,Cint),arg1,PetscReal,arg2)
end

function PetscDrawPixelToCoordinate(arg1::PetscDraw,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{Cint},arg5::Ptr{Cint})
    ccall((:PetscDrawPixelToCoordinate,petsc),PetscErrorCode,(PetscDraw,PetscInt,PetscInt,Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
end

function PetscDrawCoordinateToPixel(arg1::PetscDraw,PetscReal::Cint,arg2::Cint,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:PetscDrawCoordinateToPixel,petsc),PetscErrorCode,(PetscDraw,Cint,Cint,Ptr{PetscInt},Ptr{PetscInt}),arg1,PetscReal,arg2,arg3,arg4)
end

function PetscDrawIndicatorFunction(arg1::PetscDraw,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::Cint,arg5::Cint,arg6::Ptr{Void},arg7::Ptr{Void})
    ccall((:PetscDrawIndicatorFunction,petsc),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Cint,Cint,Ptr{Void},Ptr{Void}),arg1,PetscReal,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscDrawLine(arg1::PetscDraw,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::Cint,arg5::Cint)
    ccall((:PetscDrawLine,petsc),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4,arg5)
end

function PetscDrawArrow(arg1::PetscDraw,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::Cint,arg5::Cint)
    ccall((:PetscDrawArrow,petsc),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4,arg5)
end

function PetscDrawLineSetWidth(arg1::PetscDraw,PetscReal::Cint)
    ccall((:PetscDrawLineSetWidth,petsc),PetscErrorCode,(PetscDraw,Cint),arg1,PetscReal)
end

function PetscDrawLineGetWidth(arg1::PetscDraw,arg2::Ptr{Cint})
    ccall((:PetscDrawLineGetWidth,petsc),PetscErrorCode,(PetscDraw,Ptr{Cint}),arg1,arg2)
end

function PetscDrawMarker(arg1::PetscDraw,PetscReal::Cint,arg2::Cint,arg3::Cint)
    ccall((:PetscDrawMarker,petsc),PetscErrorCode,(PetscDraw,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3)
end

function PetscDrawSetMarkerType(arg1::PetscDraw,arg2::PetscDrawMarkerType)
    ccall((:PetscDrawSetMarkerType,petsc),PetscErrorCode,(PetscDraw,PetscDrawMarkerType),arg1,arg2)
end

function PetscDrawGetMarkerType(arg1::PetscDraw,arg2::Ptr{PetscDrawMarkerType})
    ccall((:PetscDrawGetMarkerType,petsc),PetscErrorCode,(PetscDraw,Ptr{PetscDrawMarkerType}),arg1,arg2)
end

function PetscDrawPoint(arg1::PetscDraw,PetscReal::Cint,arg2::Cint,arg3::Cint)
    ccall((:PetscDrawPoint,petsc),PetscErrorCode,(PetscDraw,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3)
end

function PetscDrawPointPixel(arg1::PetscDraw,arg2::PetscInt,arg3::PetscInt,arg4::Cint)
    ccall((:PetscDrawPointPixel,petsc),PetscErrorCode,(PetscDraw,PetscInt,PetscInt,Cint),arg1,arg2,arg3,arg4)
end

function PetscDrawPointSetSize(arg1::PetscDraw,PetscReal::Cint)
    ccall((:PetscDrawPointSetSize,petsc),PetscErrorCode,(PetscDraw,Cint),arg1,PetscReal)
end

function PetscDrawRectangle(arg1::PetscDraw,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::Cint,arg5::Cint,arg6::Cint,arg7::Cint,arg8::Cint)
    ccall((:PetscDrawRectangle,petsc),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Cint,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscDrawTriangle(arg1::PetscDraw,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::Cint,arg5::Cint,arg6::Cint,arg7::Cint,arg8::Cint,arg9::Cint)
    ccall((:PetscDrawTriangle,petsc),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Cint,Cint,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function PetscDrawEllipse(arg1::PetscDraw,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::Cint,arg5::Cint)
    ccall((:PetscDrawEllipse,petsc),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4,arg5)
end

function PetscDrawTensorContourPatch(arg1::PetscDraw,arg2::Cint,arg3::Cint,arg4::Ptr{Cint},arg5::Ptr{Cint},PetscReal::Cint,arg6::Cint,arg7::Ptr{Cint})
    ccall((:PetscDrawTensorContourPatch,petsc),PetscErrorCode,(PetscDraw,Cint,Cint,Ptr{Cint},Ptr{Cint},Cint,Cint,Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,PetscReal,arg6,arg7)
end

function PetscDrawTensorContour(arg1::PetscDraw,arg2::Cint,arg3::Cint,PetscReal::Ptr{Cint},arg4::Ptr{Cint},arg5::Ptr{Cint})
    ccall((:PetscDrawTensorContour,petsc),PetscErrorCode,(PetscDraw,Cint,Cint,Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,PetscReal,arg4,arg5)
end

function PetscDrawString(arg1::PetscDraw,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::Ptr{Uint8})
    ccall((:PetscDrawString,petsc),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Ptr{Uint8}),arg1,PetscReal,arg2,arg3,arg4)
end

function PetscDrawStringCentered(arg1::PetscDraw,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::Ptr{Uint8})
    ccall((:PetscDrawStringCentered,petsc),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Ptr{Uint8}),arg1,PetscReal,arg2,arg3,arg4)
end

function PetscDrawStringBoxed(arg1::PetscDraw,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::Cint,arg5::Ptr{Uint8},arg6::Ptr{Cint},arg7::Ptr{Cint})
    ccall((:PetscDrawStringBoxed,petsc),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Cint,Ptr{Uint8},Ptr{Cint},Ptr{Cint}),arg1,PetscReal,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscDrawStringBoxedSize(arg1::PetscDraw,arg2::Ptr{Uint8},arg3::Ptr{Cint},arg4::Ptr{Cint})
    ccall((:PetscDrawStringBoxedSize,petsc),PetscErrorCode,(PetscDraw,Ptr{Uint8},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function PetscDrawStringVertical(arg1::PetscDraw,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::Ptr{Uint8})
    ccall((:PetscDrawStringVertical,petsc),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Ptr{Uint8}),arg1,PetscReal,arg2,arg3,arg4)
end

function PetscDrawStringSetSize(arg1::PetscDraw,PetscReal::Cint,arg2::Cint)
    ccall((:PetscDrawStringSetSize,petsc),PetscErrorCode,(PetscDraw,Cint,Cint),arg1,PetscReal,arg2)
end

function PetscDrawStringGetSize(arg1::PetscDraw,arg2::Ptr{Cint},arg3::Ptr{Cint})
    ccall((:PetscDrawStringGetSize,petsc),PetscErrorCode,(PetscDraw,Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3)
end

function PetscDrawSetViewPort(arg1::PetscDraw,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::Cint)
    ccall((:PetscDrawSetViewPort,petsc),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4)
end

function PetscDrawGetViewPort(arg1::PetscDraw,arg2::Ptr{Cint},arg3::Ptr{Cint},arg4::Ptr{Cint},arg5::Ptr{Cint})
    ccall((:PetscDrawGetViewPort,petsc),PetscErrorCode,(PetscDraw,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
end

function PetscDrawSplitViewPort(arg1::PetscDraw)
    ccall((:PetscDrawSplitViewPort,petsc),PetscErrorCode,(PetscDraw,),arg1)
end

function PetscDrawSetCoordinates(arg1::PetscDraw,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::Cint)
    ccall((:PetscDrawSetCoordinates,petsc),PetscErrorCode,(PetscDraw,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4)
end

function PetscDrawGetCoordinates(arg1::PetscDraw,arg2::Ptr{Cint},arg3::Ptr{Cint},arg4::Ptr{Cint},arg5::Ptr{Cint})
    ccall((:PetscDrawGetCoordinates,petsc),PetscErrorCode,(PetscDraw,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
end

function PetscDrawSetTitle(arg1::PetscDraw,arg2::Ptr{Uint8})
    ccall((:PetscDrawSetTitle,petsc),PetscErrorCode,(PetscDraw,Ptr{Uint8}),arg1,arg2)
end

function PetscDrawAppendTitle(arg1::PetscDraw,arg2::Ptr{Uint8})
    ccall((:PetscDrawAppendTitle,petsc),PetscErrorCode,(PetscDraw,Ptr{Uint8}),arg1,arg2)
end

function PetscDrawGetTitle(arg1::PetscDraw,arg2::Ptr{Ptr{Uint8}})
    ccall((:PetscDrawGetTitle,petsc),PetscErrorCode,(PetscDraw,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function PetscDrawSetPause(arg1::PetscDraw,PetscReal::Cint)
    ccall((:PetscDrawSetPause,petsc),PetscErrorCode,(PetscDraw,Cint),arg1,PetscReal)
end

function PetscDrawGetPause(arg1::PetscDraw,arg2::Ptr{Cint})
    ccall((:PetscDrawGetPause,petsc),PetscErrorCode,(PetscDraw,Ptr{Cint}),arg1,arg2)
end

function PetscDrawPause(arg1::PetscDraw)
    ccall((:PetscDrawPause,petsc),PetscErrorCode,(PetscDraw,),arg1)
end

function PetscDrawSetDoubleBuffer(arg1::PetscDraw)
    ccall((:PetscDrawSetDoubleBuffer,petsc),PetscErrorCode,(PetscDraw,),arg1)
end

function PetscDrawFlush(arg1::PetscDraw)
    ccall((:PetscDrawFlush,petsc),PetscErrorCode,(PetscDraw,),arg1)
end

function PetscDrawSynchronizedFlush(arg1::PetscDraw)
    ccall((:PetscDrawSynchronizedFlush,petsc),PetscErrorCode,(PetscDraw,),arg1)
end

function PetscDrawClear(arg1::PetscDraw)
    ccall((:PetscDrawClear,petsc),PetscErrorCode,(PetscDraw,),arg1)
end

function PetscDrawSave(arg1::PetscDraw)
    ccall((:PetscDrawSave,petsc),PetscErrorCode,(PetscDraw,),arg1)
end

function PetscDrawSynchronizedClear(arg1::PetscDraw)
    ccall((:PetscDrawSynchronizedClear,petsc),PetscErrorCode,(PetscDraw,),arg1)
end

function PetscDrawBOP(arg1::PetscDraw)
    ccall((:PetscDrawBOP,petsc),PetscErrorCode,(PetscDraw,),arg1)
end

function PetscDrawEOP(arg1::PetscDraw)
    ccall((:PetscDrawEOP,petsc),PetscErrorCode,(PetscDraw,),arg1)
end

function PetscDrawSetDisplay(arg1::PetscDraw,arg2::Ptr{Uint8})
    ccall((:PetscDrawSetDisplay,petsc),PetscErrorCode,(PetscDraw,Ptr{Uint8}),arg1,arg2)
end

function PetscDrawGetSingleton(arg1::PetscDraw,arg2::Ptr{PetscDraw})
    ccall((:PetscDrawGetSingleton,petsc),PetscErrorCode,(PetscDraw,Ptr{PetscDraw}),arg1,arg2)
end

function PetscDrawRestoreSingleton(arg1::PetscDraw,arg2::Ptr{PetscDraw})
    ccall((:PetscDrawRestoreSingleton,petsc),PetscErrorCode,(PetscDraw,Ptr{PetscDraw}),arg1,arg2)
end

function PetscDrawGetCurrentPoint(arg1::PetscDraw,arg2::Ptr{Cint},arg3::Ptr{Cint})
    ccall((:PetscDrawGetCurrentPoint,petsc),PetscErrorCode,(PetscDraw,Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3)
end

function PetscDrawSetCurrentPoint(arg1::PetscDraw,PetscReal::Cint,arg2::Cint)
    ccall((:PetscDrawSetCurrentPoint,petsc),PetscErrorCode,(PetscDraw,Cint,Cint),arg1,PetscReal,arg2)
end

function PetscDrawPushCurrentPoint(arg1::PetscDraw,PetscReal::Cint,arg2::Cint)
    ccall((:PetscDrawPushCurrentPoint,petsc),PetscErrorCode,(PetscDraw,Cint,Cint),arg1,PetscReal,arg2)
end

function PetscDrawPopCurrentPoint(arg1::PetscDraw)
    ccall((:PetscDrawPopCurrentPoint,petsc),PetscErrorCode,(PetscDraw,),arg1)
end

function PetscDrawGetBoundingBox(arg1::PetscDraw,arg2::Ptr{Cint},arg3::Ptr{Cint},arg4::Ptr{Cint},arg5::Ptr{Cint})
    ccall((:PetscDrawGetBoundingBox,petsc),PetscErrorCode,(PetscDraw,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
end

function PetscDrawGetMouseButton(arg1::PetscDraw,arg2::Ptr{PetscDrawButton},arg3::Ptr{Cint},arg4::Ptr{Cint},arg5::Ptr{Cint},arg6::Ptr{Cint})
    ccall((:PetscDrawGetMouseButton,petsc),PetscErrorCode,(PetscDraw,Ptr{PetscDrawButton},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscDrawSynchronizedGetMouseButton(arg1::PetscDraw,arg2::Ptr{PetscDrawButton},arg3::Ptr{Cint},arg4::Ptr{Cint},arg5::Ptr{Cint},arg6::Ptr{Cint})
    ccall((:PetscDrawSynchronizedGetMouseButton,petsc),PetscErrorCode,(PetscDraw,Ptr{PetscDrawButton},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscDrawZoom(arg1::PetscDraw,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:PetscDrawZoom,petsc),PetscErrorCode,(PetscDraw,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function PetscDrawViewPortsCreate(arg1::PetscDraw,arg2::PetscInt,arg3::Ptr{Ptr{PetscDrawViewPorts}})
    ccall((:PetscDrawViewPortsCreate,petsc),PetscErrorCode,(PetscDraw,PetscInt,Ptr{Ptr{PetscDrawViewPorts}}),arg1,arg2,arg3)
end

function PetscDrawViewPortsCreateRect(arg1::PetscDraw,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{Ptr{PetscDrawViewPorts}})
    ccall((:PetscDrawViewPortsCreateRect,petsc),PetscErrorCode,(PetscDraw,PetscInt,PetscInt,Ptr{Ptr{PetscDrawViewPorts}}),arg1,arg2,arg3,arg4)
end

function PetscDrawViewPortsDestroy(arg1::Ptr{PetscDrawViewPorts})
    ccall((:PetscDrawViewPortsDestroy,petsc),PetscErrorCode,(Ptr{PetscDrawViewPorts},),arg1)
end

function PetscDrawViewPortsSet(arg1::Ptr{PetscDrawViewPorts},arg2::PetscInt)
    ccall((:PetscDrawViewPortsSet,petsc),PetscErrorCode,(Ptr{PetscDrawViewPorts},PetscInt),arg1,arg2)
end

function PetscDrawAxisCreate(arg1::PetscDraw,arg2::Ptr{PetscDrawAxis})
    ccall((:PetscDrawAxisCreate,petsc),PetscErrorCode,(PetscDraw,Ptr{PetscDrawAxis}),arg1,arg2)
end

function PetscDrawAxisDestroy(arg1::Ptr{PetscDrawAxis})
    ccall((:PetscDrawAxisDestroy,petsc),PetscErrorCode,(Ptr{PetscDrawAxis},),arg1)
end

function PetscDrawAxisDraw(arg1::PetscDrawAxis)
    ccall((:PetscDrawAxisDraw,petsc),PetscErrorCode,(PetscDrawAxis,),arg1)
end

function PetscDrawAxisSetLimits(arg1::PetscDrawAxis,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::Cint)
    ccall((:PetscDrawAxisSetLimits,petsc),PetscErrorCode,(PetscDrawAxis,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4)
end

function PetscDrawAxisGetLimits(arg1::PetscDrawAxis,arg2::Ptr{Cint},arg3::Ptr{Cint},arg4::Ptr{Cint},arg5::Ptr{Cint})
    ccall((:PetscDrawAxisGetLimits,petsc),PetscErrorCode,(PetscDrawAxis,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
end

function PetscDrawAxisSetHoldLimits(arg1::PetscDrawAxis,arg2::PetscBool)
    ccall((:PetscDrawAxisSetHoldLimits,petsc),PetscErrorCode,(PetscDrawAxis,PetscBool),arg1,arg2)
end

function PetscDrawAxisSetColors(arg1::PetscDrawAxis,arg2::Cint,arg3::Cint,arg4::Cint)
    ccall((:PetscDrawAxisSetColors,petsc),PetscErrorCode,(PetscDrawAxis,Cint,Cint,Cint),arg1,arg2,arg3,arg4)
end

function PetscDrawAxisSetLabels(arg1::PetscDrawAxis,arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8})
    ccall((:PetscDrawAxisSetLabels,petsc),PetscErrorCode,(PetscDrawAxis,Ptr{Uint8},Ptr{Uint8},Ptr{Uint8}),arg1,arg2,arg3,arg4)
end

function PetscDrawLGCreate(arg1::PetscDraw,arg2::PetscInt,arg3::Ptr{PetscDrawLG})
    ccall((:PetscDrawLGCreate,petsc),PetscErrorCode,(PetscDraw,PetscInt,Ptr{PetscDrawLG}),arg1,arg2,arg3)
end

function PetscDrawLGDestroy(arg1::Ptr{PetscDrawLG})
    ccall((:PetscDrawLGDestroy,petsc),PetscErrorCode,(Ptr{PetscDrawLG},),arg1)
end

function PetscDrawLGAddPoint(arg1::PetscDrawLG,arg2::Ptr{Cint},arg3::Ptr{Cint})
    ccall((:PetscDrawLGAddPoint,petsc),PetscErrorCode,(PetscDrawLG,Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3)
end

function PetscDrawLGAddCommonPoint(arg1::PetscDrawLG,PetscReal::Cint,arg2::Ptr{Cint})
    ccall((:PetscDrawLGAddCommonPoint,petsc),PetscErrorCode,(PetscDrawLG,Cint,Ptr{Cint}),arg1,PetscReal,arg2)
end

function PetscDrawLGAddPoints(arg1::PetscDrawLG,arg2::PetscInt,arg3::Ptr{Ptr{Cint}},arg4::Ptr{Ptr{Cint}})
    ccall((:PetscDrawLGAddPoints,petsc),PetscErrorCode,(PetscDrawLG,PetscInt,Ptr{Ptr{Cint}},Ptr{Ptr{Cint}}),arg1,arg2,arg3,arg4)
end

function PetscDrawLGDraw(arg1::PetscDrawLG)
    ccall((:PetscDrawLGDraw,petsc),PetscErrorCode,(PetscDrawLG,),arg1)
end

function PetscDrawLGView(arg1::PetscDrawLG,arg2::PetscViewer)
    ccall((:PetscDrawLGView,petsc),PetscErrorCode,(PetscDrawLG,PetscViewer),arg1,arg2)
end

function PetscDrawLGReset(arg1::PetscDrawLG)
    ccall((:PetscDrawLGReset,petsc),PetscErrorCode,(PetscDrawLG,),arg1)
end

function PetscDrawLGSetDimension(arg1::PetscDrawLG,arg2::PetscInt)
    ccall((:PetscDrawLGSetDimension,petsc),PetscErrorCode,(PetscDrawLG,PetscInt),arg1,arg2)
end

function PetscDrawLGGetDimension(arg1::PetscDrawLG,arg2::Ptr{PetscInt})
    ccall((:PetscDrawLGGetDimension,petsc),PetscErrorCode,(PetscDrawLG,Ptr{PetscInt}),arg1,arg2)
end

function PetscDrawLGSetLegend(arg1::PetscDrawLG,arg2::Ptr{Ptr{Uint8}})
    ccall((:PetscDrawLGSetLegend,petsc),PetscErrorCode,(PetscDrawLG,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function PetscDrawLGGetAxis(arg1::PetscDrawLG,arg2::Ptr{PetscDrawAxis})
    ccall((:PetscDrawLGGetAxis,petsc),PetscErrorCode,(PetscDrawLG,Ptr{PetscDrawAxis}),arg1,arg2)
end

function PetscDrawLGGetDraw(arg1::PetscDrawLG,arg2::Ptr{PetscDraw})
    ccall((:PetscDrawLGGetDraw,petsc),PetscErrorCode,(PetscDrawLG,Ptr{PetscDraw}),arg1,arg2)
end

function PetscDrawLGSetUseMarkers(arg1::PetscDrawLG,arg2::PetscBool)
    ccall((:PetscDrawLGSetUseMarkers,petsc),PetscErrorCode,(PetscDrawLG,PetscBool),arg1,arg2)
end

function PetscDrawLGSetLimits(arg1::PetscDrawLG,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::Cint)
    ccall((:PetscDrawLGSetLimits,petsc),PetscErrorCode,(PetscDrawLG,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4)
end

function PetscDrawLGSetColors(arg1::PetscDrawLG,arg2::Ptr{Cint})
    ccall((:PetscDrawLGSetColors,petsc),PetscErrorCode,(PetscDrawLG,Ptr{Cint}),arg1,arg2)
end

function PetscDrawLGSetFromOptions(arg1::PetscDrawLG)
    ccall((:PetscDrawLGSetFromOptions,petsc),PetscErrorCode,(PetscDrawLG,),arg1)
end

function PetscDrawSPCreate(arg1::PetscDraw,arg2::Cint,arg3::Ptr{PetscDrawSP})
    ccall((:PetscDrawSPCreate,petsc),PetscErrorCode,(PetscDraw,Cint,Ptr{PetscDrawSP}),arg1,arg2,arg3)
end

function PetscDrawSPDestroy(arg1::Ptr{PetscDrawSP})
    ccall((:PetscDrawSPDestroy,petsc),PetscErrorCode,(Ptr{PetscDrawSP},),arg1)
end

function PetscDrawSPAddPoint(arg1::PetscDrawSP,arg2::Ptr{Cint},arg3::Ptr{Cint})
    ccall((:PetscDrawSPAddPoint,petsc),PetscErrorCode,(PetscDrawSP,Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3)
end

function PetscDrawSPAddPoints(arg1::PetscDrawSP,arg2::Cint,arg3::Ptr{Ptr{Cint}},arg4::Ptr{Ptr{Cint}})
    ccall((:PetscDrawSPAddPoints,petsc),PetscErrorCode,(PetscDrawSP,Cint,Ptr{Ptr{Cint}},Ptr{Ptr{Cint}}),arg1,arg2,arg3,arg4)
end

function PetscDrawSPDraw(arg1::PetscDrawSP,arg2::PetscBool)
    ccall((:PetscDrawSPDraw,petsc),PetscErrorCode,(PetscDrawSP,PetscBool),arg1,arg2)
end

function PetscDrawSPReset(arg1::PetscDrawSP)
    ccall((:PetscDrawSPReset,petsc),PetscErrorCode,(PetscDrawSP,),arg1)
end

function PetscDrawSPSetDimension(arg1::PetscDrawSP,arg2::Cint)
    ccall((:PetscDrawSPSetDimension,petsc),PetscErrorCode,(PetscDrawSP,Cint),arg1,arg2)
end

function PetscDrawSPGetAxis(arg1::PetscDrawSP,arg2::Ptr{PetscDrawAxis})
    ccall((:PetscDrawSPGetAxis,petsc),PetscErrorCode,(PetscDrawSP,Ptr{PetscDrawAxis}),arg1,arg2)
end

function PetscDrawSPGetDraw(arg1::PetscDrawSP,arg2::Ptr{PetscDraw})
    ccall((:PetscDrawSPGetDraw,petsc),PetscErrorCode,(PetscDrawSP,Ptr{PetscDraw}),arg1,arg2)
end

function PetscDrawSPSetLimits(arg1::PetscDrawSP,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::Cint)
    ccall((:PetscDrawSPSetLimits,petsc),PetscErrorCode,(PetscDrawSP,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4)
end

function PetscDrawLGSPDraw(arg1::PetscDrawLG,arg2::PetscDrawSP)
    ccall((:PetscDrawLGSPDraw,petsc),PetscErrorCode,(PetscDrawLG,PetscDrawSP),arg1,arg2)
end

function PetscDrawHGCreate(arg1::PetscDraw,arg2::Cint,arg3::Ptr{PetscDrawHG})
    ccall((:PetscDrawHGCreate,petsc),PetscErrorCode,(PetscDraw,Cint,Ptr{PetscDrawHG}),arg1,arg2,arg3)
end

function PetscDrawHGDestroy(arg1::Ptr{PetscDrawHG})
    ccall((:PetscDrawHGDestroy,petsc),PetscErrorCode,(Ptr{PetscDrawHG},),arg1)
end

function PetscDrawHGAddValue(arg1::PetscDrawHG,PetscReal::Cint)
    ccall((:PetscDrawHGAddValue,petsc),PetscErrorCode,(PetscDrawHG,Cint),arg1,PetscReal)
end

function PetscDrawHGDraw(arg1::PetscDrawHG)
    ccall((:PetscDrawHGDraw,petsc),PetscErrorCode,(PetscDrawHG,),arg1)
end

function PetscDrawHGView(arg1::PetscDrawHG,arg2::PetscViewer)
    ccall((:PetscDrawHGView,petsc),PetscErrorCode,(PetscDrawHG,PetscViewer),arg1,arg2)
end

function PetscDrawHGReset(arg1::PetscDrawHG)
    ccall((:PetscDrawHGReset,petsc),PetscErrorCode,(PetscDrawHG,),arg1)
end

function PetscDrawHGGetAxis(arg1::PetscDrawHG,arg2::Ptr{PetscDrawAxis})
    ccall((:PetscDrawHGGetAxis,petsc),PetscErrorCode,(PetscDrawHG,Ptr{PetscDrawAxis}),arg1,arg2)
end

function PetscDrawHGGetDraw(arg1::PetscDrawHG,arg2::Ptr{PetscDraw})
    ccall((:PetscDrawHGGetDraw,petsc),PetscErrorCode,(PetscDrawHG,Ptr{PetscDraw}),arg1,arg2)
end

function PetscDrawHGSetLimits(arg1::PetscDrawHG,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::Cint)
    ccall((:PetscDrawHGSetLimits,petsc),PetscErrorCode,(PetscDrawHG,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4)
end

function PetscDrawHGSetNumberBins(arg1::PetscDrawHG,arg2::Cint)
    ccall((:PetscDrawHGSetNumberBins,petsc),PetscErrorCode,(PetscDrawHG,Cint),arg1,arg2)
end

function PetscDrawHGSetColor(arg1::PetscDrawHG,arg2::Cint)
    ccall((:PetscDrawHGSetColor,petsc),PetscErrorCode,(PetscDrawHG,Cint),arg1,arg2)
end

function PetscDrawHGCalcStats(arg1::PetscDrawHG,arg2::PetscBool)
    ccall((:PetscDrawHGCalcStats,petsc),PetscErrorCode,(PetscDrawHG,PetscBool),arg1,arg2)
end

function PetscDrawHGIntegerBins(arg1::PetscDrawHG,arg2::PetscBool)
    ccall((:PetscDrawHGIntegerBins,petsc),PetscErrorCode,(PetscDrawHG,PetscBool),arg1,arg2)
end

function PetscDrawBarCreate(arg1::PetscDraw,arg2::Ptr{PetscDrawBar})
    ccall((:PetscDrawBarCreate,petsc),PetscErrorCode,(PetscDraw,Ptr{PetscDrawBar}),arg1,arg2)
end

function PetscDrawBarSetData(arg1::PetscDrawBar,arg2::PetscInt,PetscReal::Ptr{Cint},arg3::Ptr{Ptr{Uint8}})
    ccall((:PetscDrawBarSetData,petsc),PetscErrorCode,(PetscDrawBar,PetscInt,Ptr{Cint},Ptr{Ptr{Uint8}}),arg1,arg2,PetscReal,arg3)
end

function PetscDrawBarDestroy(arg1::Ptr{PetscDrawBar})
    ccall((:PetscDrawBarDestroy,petsc),PetscErrorCode,(Ptr{PetscDrawBar},),arg1)
end

function PetscDrawBarDraw(arg1::PetscDrawBar)
    ccall((:PetscDrawBarDraw,petsc),PetscErrorCode,(PetscDrawBar,),arg1)
end

function PetscDrawBarSetColor(arg1::PetscDrawBar,arg2::Cint)
    ccall((:PetscDrawBarSetColor,petsc),PetscErrorCode,(PetscDrawBar,Cint),arg1,arg2)
end

function PetscDrawBarSetLimits(arg1::PetscDrawBar,PetscReal::Cint,arg2::Cint)
    ccall((:PetscDrawBarSetLimits,petsc),PetscErrorCode,(PetscDrawBar,Cint,Cint),arg1,PetscReal,arg2)
end

function PetscDrawBarSort(arg1::PetscDrawBar,arg2::PetscBool,PetscReal::Cint)
    ccall((:PetscDrawBarSort,petsc),PetscErrorCode,(PetscDrawBar,PetscBool,Cint),arg1,arg2,PetscReal)
end

function PetscDrawBarSetFromOptions(arg1::PetscDrawBar)
    ccall((:PetscDrawBarSetFromOptions,petsc),PetscErrorCode,(PetscDrawBar,),arg1)
end

function PetscDrawBarGetAxis(arg1::PetscDrawBar,arg2::Ptr{PetscDrawAxis})
    ccall((:PetscDrawBarGetAxis,petsc),PetscErrorCode,(PetscDrawBar,Ptr{PetscDrawAxis}),arg1,arg2)
end

function PetscDrawBarGetDraw(arg1::PetscDrawBar,arg2::Ptr{PetscDraw})
    ccall((:PetscDrawBarGetDraw,petsc),PetscErrorCode,(PetscDrawBar,Ptr{PetscDraw}),arg1,arg2)
end

function PetscViewerDrawGetDraw(arg1::PetscViewer,arg2::PetscInt,arg3::Ptr{PetscDraw})
    ccall((:PetscViewerDrawGetDraw,petsc),PetscErrorCode,(PetscViewer,PetscInt,Ptr{PetscDraw}),arg1,arg2,arg3)
end

function PetscViewerDrawBaseAdd(arg1::PetscViewer,arg2::PetscInt)
    ccall((:PetscViewerDrawBaseAdd,petsc),PetscErrorCode,(PetscViewer,PetscInt),arg1,arg2)
end

function PetscViewerDrawBaseSet(arg1::PetscViewer,arg2::PetscInt)
    ccall((:PetscViewerDrawBaseSet,petsc),PetscErrorCode,(PetscViewer,PetscInt),arg1,arg2)
end

function PetscViewerDrawGetDrawLG(arg1::PetscViewer,arg2::PetscInt,arg3::Ptr{PetscDrawLG})
    ccall((:PetscViewerDrawGetDrawLG,petsc),PetscErrorCode,(PetscViewer,PetscInt,Ptr{PetscDrawLG}),arg1,arg2,arg3)
end

function PetscViewerDrawGetDrawAxis(arg1::PetscViewer,arg2::PetscInt,arg3::Ptr{PetscDrawAxis})
    ccall((:PetscViewerDrawGetDrawAxis,petsc),PetscErrorCode,(PetscViewer,PetscInt,Ptr{PetscDrawAxis}),arg1,arg2,arg3)
end

function PetscDrawUtilitySetCmapHue(arg1::Ptr{Cuchar},arg2::Ptr{Cuchar},arg3::Ptr{Cuchar},arg4::Cint)
    ccall((:PetscDrawUtilitySetCmapHue,petsc),PetscErrorCode,(Ptr{Cuchar},Ptr{Cuchar},Ptr{Cuchar},Cint),arg1,arg2,arg3,arg4)
end

function PetscDrawUtilitySetGamma()
    ccall((:PetscDrawUtilitySetGamma,petsc),PetscErrorCode,())
end

function ISInitializePackage()
    ccall((:ISInitializePackage,petsc),PetscErrorCode,())
end

function ISSetType(arg1::IS,arg2::ISType)
    ccall((:ISSetType,petsc),PetscErrorCode,(IS,ISType),arg1,arg2)
end

function ISGetType(arg1::IS,arg2::Ptr{ISType})
    ccall((:ISGetType,petsc),PetscErrorCode,(IS,Ptr{ISType}),arg1,arg2)
end

function ISRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:ISRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function ISCreate(arg1::MPI_Comm,arg2::Ptr{IS})
    ccall((:ISCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{IS}),arg1,arg2)
end

function ISCreateGeneral(arg1::MPI_Comm,arg2::PetscInt,arg3::Ptr{PetscInt},PetscCopyMode::Cint,arg4::Ptr{IS})
    ccall((:ISCreateGeneral,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{PetscInt},Cint,Ptr{IS}),arg1,arg2,arg3,PetscCopyMode,arg4)
end

function ISGeneralSetIndices(arg1::IS,arg2::PetscInt,arg3::Ptr{PetscInt},PetscCopyMode::Cint)
    ccall((:ISGeneralSetIndices,petsc),PetscErrorCode,(IS,PetscInt,Ptr{PetscInt},Cint),arg1,arg2,arg3,PetscCopyMode)
end

function ISCreateBlock(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscInt},PetscCopyMode::Cint,arg5::Ptr{IS})
    ccall((:ISCreateBlock,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,Ptr{PetscInt},Cint,Ptr{IS}),arg1,arg2,arg3,arg4,PetscCopyMode,arg5)
end

function ISBlockSetIndices(arg1::IS,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscInt},PetscCopyMode::Cint)
    ccall((:ISBlockSetIndices,petsc),PetscErrorCode,(IS,PetscInt,PetscInt,Ptr{PetscInt},Cint),arg1,arg2,arg3,arg4,PetscCopyMode)
end

function ISCreateStride(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::Ptr{IS})
    ccall((:ISCreateStride,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,Ptr{IS}),arg1,arg2,arg3,arg4,arg5)
end

function ISStrideSetStride(arg1::IS,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt)
    ccall((:ISStrideSetStride,petsc),PetscErrorCode,(IS,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end

function ISDestroy(arg1::Ptr{IS})
    ccall((:ISDestroy,petsc),PetscErrorCode,(Ptr{IS},),arg1)
end

function ISSetPermutation(arg1::IS)
    ccall((:ISSetPermutation,petsc),PetscErrorCode,(IS,),arg1)
end

function ISPermutation(arg1::IS,arg2::Ptr{PetscBool})
    ccall((:ISPermutation,petsc),PetscErrorCode,(IS,Ptr{PetscBool}),arg1,arg2)
end

function ISSetIdentity(arg1::IS)
    ccall((:ISSetIdentity,petsc),PetscErrorCode,(IS,),arg1)
end

function ISIdentity(arg1::IS,arg2::Ptr{PetscBool})
    ccall((:ISIdentity,petsc),PetscErrorCode,(IS,Ptr{PetscBool}),arg1,arg2)
end

function ISContiguousLocal(arg1::IS,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscInt},arg5::Ptr{PetscBool})
    ccall((:ISContiguousLocal,petsc),PetscErrorCode,(IS,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function ISGetIndices(arg1::IS,arg2::Ptr{Ptr{PetscInt}})
    ccall((:ISGetIndices,petsc),PetscErrorCode,(IS,Ptr{Ptr{PetscInt}}),arg1,arg2)
end

function ISRestoreIndices(arg1::IS,arg2::Ptr{Ptr{PetscInt}})
    ccall((:ISRestoreIndices,petsc),PetscErrorCode,(IS,Ptr{Ptr{PetscInt}}),arg1,arg2)
end

function ISGetTotalIndices(arg1::IS,arg2::Ptr{Ptr{PetscInt}})
    ccall((:ISGetTotalIndices,petsc),PetscErrorCode,(IS,Ptr{Ptr{PetscInt}}),arg1,arg2)
end

function ISRestoreTotalIndices(arg1::IS,arg2::Ptr{Ptr{PetscInt}})
    ccall((:ISRestoreTotalIndices,petsc),PetscErrorCode,(IS,Ptr{Ptr{PetscInt}}),arg1,arg2)
end

function ISGetNonlocalIndices(arg1::IS,arg2::Ptr{Ptr{PetscInt}})
    ccall((:ISGetNonlocalIndices,petsc),PetscErrorCode,(IS,Ptr{Ptr{PetscInt}}),arg1,arg2)
end

function ISRestoreNonlocalIndices(arg1::IS,arg2::Ptr{Ptr{PetscInt}})
    ccall((:ISRestoreNonlocalIndices,petsc),PetscErrorCode,(IS,Ptr{Ptr{PetscInt}}),arg1,arg2)
end

function ISGetNonlocalIS(arg1::IS,is::Ptr{IS})
    ccall((:ISGetNonlocalIS,petsc),PetscErrorCode,(IS,Ptr{IS}),arg1,is)
end

function ISRestoreNonlocalIS(arg1::IS,is::Ptr{IS})
    ccall((:ISRestoreNonlocalIS,petsc),PetscErrorCode,(IS,Ptr{IS}),arg1,is)
end

function ISGetSize(arg1::IS,arg2::Ptr{PetscInt})
    ccall((:ISGetSize,petsc),PetscErrorCode,(IS,Ptr{PetscInt}),arg1,arg2)
end

function ISGetLocalSize(arg1::IS,arg2::Ptr{PetscInt})
    ccall((:ISGetLocalSize,petsc),PetscErrorCode,(IS,Ptr{PetscInt}),arg1,arg2)
end

function ISInvertPermutation(arg1::IS,arg2::PetscInt,arg3::Ptr{IS})
    ccall((:ISInvertPermutation,petsc),PetscErrorCode,(IS,PetscInt,Ptr{IS}),arg1,arg2,arg3)
end

function ISView(arg1::IS,arg2::PetscViewer)
    ccall((:ISView,petsc),PetscErrorCode,(IS,PetscViewer),arg1,arg2)
end

function ISEqual(arg1::IS,arg2::IS,arg3::Ptr{PetscBool})
    ccall((:ISEqual,petsc),PetscErrorCode,(IS,IS,Ptr{PetscBool}),arg1,arg2,arg3)
end

function ISSort(arg1::IS)
    ccall((:ISSort,petsc),PetscErrorCode,(IS,),arg1)
end

function ISSortRemoveDups(arg1::IS)
    ccall((:ISSortRemoveDups,petsc),PetscErrorCode,(IS,),arg1)
end

function ISSorted(arg1::IS,arg2::Ptr{PetscBool})
    ccall((:ISSorted,petsc),PetscErrorCode,(IS,Ptr{PetscBool}),arg1,arg2)
end

function ISDifference(arg1::IS,arg2::IS,arg3::Ptr{IS})
    ccall((:ISDifference,petsc),PetscErrorCode,(IS,IS,Ptr{IS}),arg1,arg2,arg3)
end

function ISSum(arg1::IS,arg2::IS,arg3::Ptr{IS})
    ccall((:ISSum,petsc),PetscErrorCode,(IS,IS,Ptr{IS}),arg1,arg2,arg3)
end

function ISExpand(arg1::IS,arg2::IS,arg3::Ptr{IS})
    ccall((:ISExpand,petsc),PetscErrorCode,(IS,IS,Ptr{IS}),arg1,arg2,arg3)
end

function ISGetMinMax(arg1::IS,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt})
    ccall((:ISGetMinMax,petsc),PetscErrorCode,(IS,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function ISBlockGetIndices(arg1::IS,arg2::Ptr{Ptr{PetscInt}})
    ccall((:ISBlockGetIndices,petsc),PetscErrorCode,(IS,Ptr{Ptr{PetscInt}}),arg1,arg2)
end

function ISBlockRestoreIndices(arg1::IS,arg2::Ptr{Ptr{PetscInt}})
    ccall((:ISBlockRestoreIndices,petsc),PetscErrorCode,(IS,Ptr{Ptr{PetscInt}}),arg1,arg2)
end

function ISBlockGetLocalSize(arg1::IS,arg2::Ptr{PetscInt})
    ccall((:ISBlockGetLocalSize,petsc),PetscErrorCode,(IS,Ptr{PetscInt}),arg1,arg2)
end

function ISBlockGetSize(arg1::IS,arg2::Ptr{PetscInt})
    ccall((:ISBlockGetSize,petsc),PetscErrorCode,(IS,Ptr{PetscInt}),arg1,arg2)
end

function ISGetBlockSize(arg1::IS,arg2::Ptr{PetscInt})
    ccall((:ISGetBlockSize,petsc),PetscErrorCode,(IS,Ptr{PetscInt}),arg1,arg2)
end

function ISSetBlockSize(arg1::IS,arg2::PetscInt)
    ccall((:ISSetBlockSize,petsc),PetscErrorCode,(IS,PetscInt),arg1,arg2)
end

function ISStrideGetInfo(arg1::IS,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt})
    ccall((:ISStrideGetInfo,petsc),PetscErrorCode,(IS,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function ISToGeneral(arg1::IS)
    ccall((:ISToGeneral,petsc),PetscErrorCode,(IS,),arg1)
end

function ISDuplicate(arg1::IS,arg2::Ptr{IS})
    ccall((:ISDuplicate,petsc),PetscErrorCode,(IS,Ptr{IS}),arg1,arg2)
end

function ISCopy(arg1::IS,arg2::IS)
    ccall((:ISCopy,petsc),PetscErrorCode,(IS,IS),arg1,arg2)
end

function ISAllGather(arg1::IS,arg2::Ptr{IS})
    ccall((:ISAllGather,petsc),PetscErrorCode,(IS,Ptr{IS}),arg1,arg2)
end

function ISComplement(arg1::IS,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{IS})
    ccall((:ISComplement,petsc),PetscErrorCode,(IS,PetscInt,PetscInt,Ptr{IS}),arg1,arg2,arg3,arg4)
end

function ISConcatenate(arg1::MPI_Comm,arg2::PetscInt,arg3::Ptr{IS},arg4::Ptr{IS})
    ccall((:ISConcatenate,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{IS},Ptr{IS}),arg1,arg2,arg3,arg4)
end

function ISListToPair(arg1::MPI_Comm,arg2::PetscInt,arg3::Ptr{IS},arg4::Ptr{IS},arg5::Ptr{IS})
    ccall((:ISListToPair,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{IS},Ptr{IS},Ptr{IS}),arg1,arg2,arg3,arg4,arg5)
end

function ISPairToList(arg1::IS,arg2::IS,arg3::Ptr{PetscInt},arg4::Ptr{Ptr{IS}})
    ccall((:ISPairToList,petsc),PetscErrorCode,(IS,IS,Ptr{PetscInt},Ptr{Ptr{IS}}),arg1,arg2,arg3,arg4)
end

function ISEmbed(arg1::IS,arg2::IS,arg3::PetscBool,arg4::Ptr{IS})
    ccall((:ISEmbed,petsc),PetscErrorCode,(IS,IS,PetscBool,Ptr{IS}),arg1,arg2,arg3,arg4)
end

function ISSortPermutation(arg1::IS,arg2::PetscBool,arg3::Ptr{IS})
    ccall((:ISSortPermutation,petsc),PetscErrorCode,(IS,PetscBool,Ptr{IS}),arg1,arg2,arg3)
end

function ISOnComm(arg1::IS,arg2::MPI_Comm,PetscCopyMode::Cint,arg3::Ptr{IS})
    ccall((:ISOnComm,petsc),PetscErrorCode,(IS,MPI_Comm,Cint,Ptr{IS}),arg1,arg2,PetscCopyMode,arg3)
end

function ISLocalToGlobalMappingCreate(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscInt},PetscCopyMode::Cint,arg5::Ptr{ISLocalToGlobalMapping})
    ccall((:ISLocalToGlobalMappingCreate,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,Ptr{PetscInt},Cint,Ptr{ISLocalToGlobalMapping}),arg1,arg2,arg3,arg4,PetscCopyMode,arg5)
end

function ISLocalToGlobalMappingCreateIS(arg1::IS,arg2::Ptr{ISLocalToGlobalMapping})
    ccall((:ISLocalToGlobalMappingCreateIS,petsc),PetscErrorCode,(IS,Ptr{ISLocalToGlobalMapping}),arg1,arg2)
end

function ISLocalToGlobalMappingCreateSF(arg1::PetscSF,arg2::PetscInt,arg3::Ptr{ISLocalToGlobalMapping})
    ccall((:ISLocalToGlobalMappingCreateSF,petsc),PetscErrorCode,(PetscSF,PetscInt,Ptr{ISLocalToGlobalMapping}),arg1,arg2,arg3)
end

function ISLocalToGlobalMappingView(arg1::ISLocalToGlobalMapping,arg2::PetscViewer)
    ccall((:ISLocalToGlobalMappingView,petsc),PetscErrorCode,(ISLocalToGlobalMapping,PetscViewer),arg1,arg2)
end

function ISLocalToGlobalMappingDestroy(arg1::Ptr{ISLocalToGlobalMapping})
    ccall((:ISLocalToGlobalMappingDestroy,petsc),PetscErrorCode,(Ptr{ISLocalToGlobalMapping},),arg1)
end

function ISLocalToGlobalMappingApply(arg1::ISLocalToGlobalMapping,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:ISLocalToGlobalMappingApply,petsc),PetscErrorCode,(ISLocalToGlobalMapping,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function ISLocalToGlobalMappingApplyBlock(arg1::ISLocalToGlobalMapping,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:ISLocalToGlobalMappingApplyBlock,petsc),PetscErrorCode,(ISLocalToGlobalMapping,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function ISLocalToGlobalMappingApplyIS(arg1::ISLocalToGlobalMapping,arg2::IS,arg3::Ptr{IS})
    ccall((:ISLocalToGlobalMappingApplyIS,petsc),PetscErrorCode,(ISLocalToGlobalMapping,IS,Ptr{IS}),arg1,arg2,arg3)
end

function ISGlobalToLocalMappingApply(arg1::ISLocalToGlobalMapping,arg2::ISGlobalToLocalMappingType,arg3::PetscInt,arg4::Ptr{PetscInt},arg5::Ptr{PetscInt},arg6::Ptr{PetscInt})
    ccall((:ISGlobalToLocalMappingApply,petsc),PetscErrorCode,(ISLocalToGlobalMapping,ISGlobalToLocalMappingType,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function ISGlobalToLocalMappingApplyBlock(arg1::ISLocalToGlobalMapping,arg2::ISGlobalToLocalMappingType,arg3::PetscInt,arg4::Ptr{PetscInt},arg5::Ptr{PetscInt},arg6::Ptr{PetscInt})
    ccall((:ISGlobalToLocalMappingApplyBlock,petsc),PetscErrorCode,(ISLocalToGlobalMapping,ISGlobalToLocalMappingType,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function ISGlobalToLocalMappingApplyIS(arg1::ISLocalToGlobalMapping,arg2::ISGlobalToLocalMappingType,arg3::IS,arg4::Ptr{IS})
    ccall((:ISGlobalToLocalMappingApplyIS,petsc),PetscErrorCode,(ISLocalToGlobalMapping,ISGlobalToLocalMappingType,IS,Ptr{IS}),arg1,arg2,arg3,arg4)
end

function ISLocalToGlobalMappingGetSize(arg1::ISLocalToGlobalMapping,arg2::Ptr{PetscInt})
    ccall((:ISLocalToGlobalMappingGetSize,petsc),PetscErrorCode,(ISLocalToGlobalMapping,Ptr{PetscInt}),arg1,arg2)
end

function ISLocalToGlobalMappingGetInfo(arg1::ISLocalToGlobalMapping,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{PetscInt}},arg4::Ptr{Ptr{PetscInt}},arg5::Ptr{Ptr{Ptr{PetscInt}}})
    ccall((:ISLocalToGlobalMappingGetInfo,petsc),PetscErrorCode,(ISLocalToGlobalMapping,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{Ptr{Ptr{PetscInt}}}),arg1,arg2,arg3,arg4,arg5)
end

function ISLocalToGlobalMappingRestoreInfo(arg1::ISLocalToGlobalMapping,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{PetscInt}},arg4::Ptr{Ptr{PetscInt}},arg5::Ptr{Ptr{Ptr{PetscInt}}})
    ccall((:ISLocalToGlobalMappingRestoreInfo,petsc),PetscErrorCode,(ISLocalToGlobalMapping,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{Ptr{Ptr{PetscInt}}}),arg1,arg2,arg3,arg4,arg5)
end

function ISLocalToGlobalMappingGetBlockInfo(arg1::ISLocalToGlobalMapping,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{PetscInt}},arg4::Ptr{Ptr{PetscInt}},arg5::Ptr{Ptr{Ptr{PetscInt}}})
    ccall((:ISLocalToGlobalMappingGetBlockInfo,petsc),PetscErrorCode,(ISLocalToGlobalMapping,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{Ptr{Ptr{PetscInt}}}),arg1,arg2,arg3,arg4,arg5)
end

function ISLocalToGlobalMappingRestoreBlockInfo(arg1::ISLocalToGlobalMapping,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{PetscInt}},arg4::Ptr{Ptr{PetscInt}},arg5::Ptr{Ptr{Ptr{PetscInt}}})
    ccall((:ISLocalToGlobalMappingRestoreBlockInfo,petsc),PetscErrorCode,(ISLocalToGlobalMapping,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{Ptr{Ptr{PetscInt}}}),arg1,arg2,arg3,arg4,arg5)
end

function ISLocalToGlobalMappingGetIndices(arg1::ISLocalToGlobalMapping,arg2::Ptr{Ptr{PetscInt}})
    ccall((:ISLocalToGlobalMappingGetIndices,petsc),PetscErrorCode,(ISLocalToGlobalMapping,Ptr{Ptr{PetscInt}}),arg1,arg2)
end

function ISLocalToGlobalMappingRestoreIndices(arg1::ISLocalToGlobalMapping,arg2::Ptr{Ptr{PetscInt}})
    ccall((:ISLocalToGlobalMappingRestoreIndices,petsc),PetscErrorCode,(ISLocalToGlobalMapping,Ptr{Ptr{PetscInt}}),arg1,arg2)
end

function ISLocalToGlobalMappingGetBlockIndices(arg1::ISLocalToGlobalMapping,arg2::Ptr{Ptr{PetscInt}})
    ccall((:ISLocalToGlobalMappingGetBlockIndices,petsc),PetscErrorCode,(ISLocalToGlobalMapping,Ptr{Ptr{PetscInt}}),arg1,arg2)
end

function ISLocalToGlobalMappingRestoreBlockIndices(arg1::ISLocalToGlobalMapping,arg2::Ptr{Ptr{PetscInt}})
    ccall((:ISLocalToGlobalMappingRestoreBlockIndices,petsc),PetscErrorCode,(ISLocalToGlobalMapping,Ptr{Ptr{PetscInt}}),arg1,arg2)
end

function ISLocalToGlobalMappingConcatenate(arg1::MPI_Comm,arg2::PetscInt,arg3::Ptr{ISLocalToGlobalMapping},arg4::Ptr{ISLocalToGlobalMapping})
    ccall((:ISLocalToGlobalMappingConcatenate,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{ISLocalToGlobalMapping},Ptr{ISLocalToGlobalMapping}),arg1,arg2,arg3,arg4)
end

function ISG2LMapApply(arg1::ISLocalToGlobalMapping,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:ISG2LMapApply,petsc),PetscErrorCode,(ISLocalToGlobalMapping,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function ISLocalToGlobalMappingGetBlockSize(arg1::ISLocalToGlobalMapping,arg2::Ptr{PetscInt})
    ccall((:ISLocalToGlobalMappingGetBlockSize,petsc),PetscErrorCode,(ISLocalToGlobalMapping,Ptr{PetscInt}),arg1,arg2)
end

function ISAllGatherColors(arg1::MPI_Comm,arg2::PetscInt,arg3::Ptr{Cint},arg4::Ptr{PetscInt},arg5::Ptr{Ptr{Cint}})
    ccall((:ISAllGatherColors,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{Cint},Ptr{PetscInt},Ptr{Ptr{Cint}}),arg1,arg2,arg3,arg4,arg5)
end

function ISColoringCreate(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,ISColoringValue::Ptr{Cint},PetscCopyMode::Cint,arg4::Ptr{ISColoring})
    ccall((:ISColoringCreate,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,Ptr{Cint},Cint,Ptr{ISColoring}),arg1,arg2,arg3,ISColoringValue,PetscCopyMode,arg4)
end

function ISColoringDestroy(arg1::Ptr{ISColoring})
    ccall((:ISColoringDestroy,petsc),PetscErrorCode,(Ptr{ISColoring},),arg1)
end

function ISColoringView(arg1::ISColoring,arg2::PetscViewer)
    ccall((:ISColoringView,petsc),PetscErrorCode,(ISColoring,PetscViewer),arg1,arg2)
end

function ISColoringViewFromOptions(arg1::ISColoring,arg2::PetscObject,arg3::Ptr{Uint8})
    ccall((:ISColoringViewFromOptions,petsc),PetscErrorCode,(ISColoring,PetscObject,Ptr{Uint8}),arg1,arg2,arg3)
end

function ISColoringGetIS(arg1::ISColoring,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{IS}})
    ccall((:ISColoringGetIS,petsc),PetscErrorCode,(ISColoring,Ptr{PetscInt},Ptr{Ptr{IS}}),arg1,arg2,arg3)
end

function ISColoringRestoreIS(arg1::ISColoring,arg2::Ptr{Ptr{IS}})
    ccall((:ISColoringRestoreIS,petsc),PetscErrorCode,(ISColoring,Ptr{Ptr{IS}}),arg1,arg2)
end

function ISColoringReference(arg1::ISColoring)
    ccall((:ISColoringReference,petsc),PetscErrorCode,(ISColoring,),arg1)
end

function ISColoringSetType(arg1::ISColoring,arg2::ISColoringType)
    ccall((:ISColoringSetType,petsc),PetscErrorCode,(ISColoring,ISColoringType),arg1,arg2)
end

function ISPartitioningToNumbering(arg1::IS,arg2::Ptr{IS})
    ccall((:ISPartitioningToNumbering,petsc),PetscErrorCode,(IS,Ptr{IS}),arg1,arg2)
end

function ISPartitioningCount(arg1::IS,arg2::PetscInt,arg3::Ptr{PetscInt})
    ccall((:ISPartitioningCount,petsc),PetscErrorCode,(IS,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function ISCompressIndicesGeneral(arg1::PetscInt,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::Ptr{IS},arg6::Ptr{IS})
    ccall((:ISCompressIndicesGeneral,petsc),PetscErrorCode,(PetscInt,PetscInt,PetscInt,PetscInt,Ptr{IS},Ptr{IS}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function ISCompressIndicesSorted(arg1::PetscInt,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{IS},arg5::Ptr{IS})
    ccall((:ISCompressIndicesSorted,petsc),PetscErrorCode,(PetscInt,PetscInt,PetscInt,Ptr{IS},Ptr{IS}),arg1,arg2,arg3,arg4,arg5)
end

function ISExpandIndicesGeneral(arg1::PetscInt,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::Ptr{IS},arg6::Ptr{IS})
    ccall((:ISExpandIndicesGeneral,petsc),PetscErrorCode,(PetscInt,PetscInt,PetscInt,PetscInt,Ptr{IS},Ptr{IS}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscLayoutSetUp(arg1::PetscLayout)
    ccall((:PetscLayoutSetUp,petsc),PetscErrorCode,(PetscLayout,),arg1)
end

function PetscLayoutDestroy(arg1::Ptr{PetscLayout})
    ccall((:PetscLayoutDestroy,petsc),PetscErrorCode,(Ptr{PetscLayout},),arg1)
end

function PetscLayoutDuplicate(arg1::PetscLayout,arg2::Ptr{PetscLayout})
    ccall((:PetscLayoutDuplicate,petsc),PetscErrorCode,(PetscLayout,Ptr{PetscLayout}),arg1,arg2)
end

function PetscLayoutReference(arg1::PetscLayout,arg2::Ptr{PetscLayout})
    ccall((:PetscLayoutReference,petsc),PetscErrorCode,(PetscLayout,Ptr{PetscLayout}),arg1,arg2)
end

function PetscLayoutSetLocalSize(arg1::PetscLayout,arg2::PetscInt)
    ccall((:PetscLayoutSetLocalSize,petsc),PetscErrorCode,(PetscLayout,PetscInt),arg1,arg2)
end

function PetscLayoutGetLocalSize(arg1::PetscLayout,arg2::Ptr{PetscInt})
    ccall((:PetscLayoutGetLocalSize,petsc),PetscErrorCode,(PetscLayout,Ptr{PetscInt}),arg1,arg2)
end

function PetscLayoutSetSize(arg1::PetscLayout,arg2::PetscInt)
    ccall((:PetscLayoutSetSize,petsc),PetscErrorCode,(PetscLayout,PetscInt),arg1,arg2)
end

function PetscLayoutGetSize(arg1::PetscLayout,arg2::Ptr{PetscInt})
    ccall((:PetscLayoutGetSize,petsc),PetscErrorCode,(PetscLayout,Ptr{PetscInt}),arg1,arg2)
end

function PetscLayoutSetBlockSize(arg1::PetscLayout,arg2::PetscInt)
    ccall((:PetscLayoutSetBlockSize,petsc),PetscErrorCode,(PetscLayout,PetscInt),arg1,arg2)
end

function PetscLayoutGetBlockSize(arg1::PetscLayout,arg2::Ptr{PetscInt})
    ccall((:PetscLayoutGetBlockSize,petsc),PetscErrorCode,(PetscLayout,Ptr{PetscInt}),arg1,arg2)
end

function PetscLayoutGetRange(arg1::PetscLayout,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt})
    ccall((:PetscLayoutGetRange,petsc),PetscErrorCode,(PetscLayout,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscLayoutGetRanges(arg1::PetscLayout,arg2::Ptr{Ptr{PetscInt}})
    ccall((:PetscLayoutGetRanges,petsc),PetscErrorCode,(PetscLayout,Ptr{Ptr{PetscInt}}),arg1,arg2)
end

function PetscLayoutSetISLocalToGlobalMapping(arg1::PetscLayout,arg2::ISLocalToGlobalMapping)
    ccall((:PetscLayoutSetISLocalToGlobalMapping,petsc),PetscErrorCode,(PetscLayout,ISLocalToGlobalMapping),arg1,arg2)
end

function PetscSFSetGraphLayout(arg1::PetscSF,arg2::PetscLayout,arg3::PetscInt,arg4::Ptr{PetscInt},PetscCopyMode::Cint,arg5::Ptr{PetscInt})
    ccall((:PetscSFSetGraphLayout,petsc),PetscErrorCode,(PetscSF,PetscLayout,PetscInt,Ptr{PetscInt},Cint,Ptr{PetscInt}),arg1,arg2,arg3,arg4,PetscCopyMode,arg5)
end

function PetscSectionCreate(arg1::MPI_Comm,arg2::Ptr{PetscSection})
    ccall((:PetscSectionCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscSection}),arg1,arg2)
end

function PetscSectionClone(arg1::PetscSection,arg2::Ptr{PetscSection})
    ccall((:PetscSectionClone,petsc),PetscErrorCode,(PetscSection,Ptr{PetscSection}),arg1,arg2)
end

function PetscSectionCopy(arg1::PetscSection,arg2::PetscSection)
    ccall((:PetscSectionCopy,petsc),PetscErrorCode,(PetscSection,PetscSection),arg1,arg2)
end

function PetscSectionGetNumFields(arg1::PetscSection,arg2::Ptr{PetscInt})
    ccall((:PetscSectionGetNumFields,petsc),PetscErrorCode,(PetscSection,Ptr{PetscInt}),arg1,arg2)
end

function PetscSectionSetNumFields(arg1::PetscSection,arg2::PetscInt)
    ccall((:PetscSectionSetNumFields,petsc),PetscErrorCode,(PetscSection,PetscInt),arg1,arg2)
end

function PetscSectionGetFieldName(arg1::PetscSection,arg2::PetscInt,arg3::Ptr{Ptr{Uint8}})
    ccall((:PetscSectionGetFieldName,petsc),PetscErrorCode,(PetscSection,PetscInt,Ptr{Ptr{Uint8}}),arg1,arg2,arg3)
end

function PetscSectionSetFieldName(arg1::PetscSection,arg2::PetscInt,arg3::Ptr{Uint8})
    ccall((:PetscSectionSetFieldName,petsc),PetscErrorCode,(PetscSection,PetscInt,Ptr{Uint8}),arg1,arg2,arg3)
end

function PetscSectionGetFieldComponents(arg1::PetscSection,arg2::PetscInt,arg3::Ptr{PetscInt})
    ccall((:PetscSectionGetFieldComponents,petsc),PetscErrorCode,(PetscSection,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscSectionSetFieldComponents(arg1::PetscSection,arg2::PetscInt,arg3::PetscInt)
    ccall((:PetscSectionSetFieldComponents,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt),arg1,arg2,arg3)
end

function PetscSectionGetChart(arg1::PetscSection,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt})
    ccall((:PetscSectionGetChart,petsc),PetscErrorCode,(PetscSection,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscSectionSetChart(arg1::PetscSection,arg2::PetscInt,arg3::PetscInt)
    ccall((:PetscSectionSetChart,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt),arg1,arg2,arg3)
end

function PetscSectionGetPermutation(arg1::PetscSection,arg2::Ptr{IS})
    ccall((:PetscSectionGetPermutation,petsc),PetscErrorCode,(PetscSection,Ptr{IS}),arg1,arg2)
end

function PetscSectionSetPermutation(arg1::PetscSection,arg2::IS)
    ccall((:PetscSectionSetPermutation,petsc),PetscErrorCode,(PetscSection,IS),arg1,arg2)
end

function PetscSectionGetDof(arg1::PetscSection,arg2::PetscInt,arg3::Ptr{PetscInt})
    ccall((:PetscSectionGetDof,petsc),PetscErrorCode,(PetscSection,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscSectionSetDof(arg1::PetscSection,arg2::PetscInt,arg3::PetscInt)
    ccall((:PetscSectionSetDof,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt),arg1,arg2,arg3)
end

function PetscSectionAddDof(arg1::PetscSection,arg2::PetscInt,arg3::PetscInt)
    ccall((:PetscSectionAddDof,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt),arg1,arg2,arg3)
end

function PetscSectionGetFieldDof(arg1::PetscSection,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscInt})
    ccall((:PetscSectionGetFieldDof,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function PetscSectionSetFieldDof(arg1::PetscSection,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt)
    ccall((:PetscSectionSetFieldDof,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end

function PetscSectionAddFieldDof(arg1::PetscSection,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt)
    ccall((:PetscSectionAddFieldDof,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end

function PetscSectionHasConstraints(arg1::PetscSection,arg2::Ptr{PetscBool})
    ccall((:PetscSectionHasConstraints,petsc),PetscErrorCode,(PetscSection,Ptr{PetscBool}),arg1,arg2)
end

function PetscSectionGetConstraintDof(arg1::PetscSection,arg2::PetscInt,arg3::Ptr{PetscInt})
    ccall((:PetscSectionGetConstraintDof,petsc),PetscErrorCode,(PetscSection,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscSectionSetConstraintDof(arg1::PetscSection,arg2::PetscInt,arg3::PetscInt)
    ccall((:PetscSectionSetConstraintDof,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt),arg1,arg2,arg3)
end

function PetscSectionAddConstraintDof(arg1::PetscSection,arg2::PetscInt,arg3::PetscInt)
    ccall((:PetscSectionAddConstraintDof,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt),arg1,arg2,arg3)
end

function PetscSectionGetFieldConstraintDof(arg1::PetscSection,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscInt})
    ccall((:PetscSectionGetFieldConstraintDof,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function PetscSectionSetFieldConstraintDof(arg1::PetscSection,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt)
    ccall((:PetscSectionSetFieldConstraintDof,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end

function PetscSectionAddFieldConstraintDof(arg1::PetscSection,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt)
    ccall((:PetscSectionAddFieldConstraintDof,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end

function PetscSectionGetConstraintIndices(arg1::PetscSection,arg2::PetscInt,arg3::Ptr{Ptr{PetscInt}})
    ccall((:PetscSectionGetConstraintIndices,petsc),PetscErrorCode,(PetscSection,PetscInt,Ptr{Ptr{PetscInt}}),arg1,arg2,arg3)
end

function PetscSectionSetConstraintIndices(arg1::PetscSection,arg2::PetscInt,arg3::Ptr{PetscInt})
    ccall((:PetscSectionSetConstraintIndices,petsc),PetscErrorCode,(PetscSection,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscSectionGetFieldConstraintIndices(arg1::PetscSection,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{Ptr{PetscInt}})
    ccall((:PetscSectionGetFieldConstraintIndices,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt,Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4)
end

function PetscSectionSetFieldConstraintIndices(arg1::PetscSection,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscInt})
    ccall((:PetscSectionSetFieldConstraintIndices,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function PetscSectionSetUpBC(arg1::PetscSection)
    ccall((:PetscSectionSetUpBC,petsc),PetscErrorCode,(PetscSection,),arg1)
end

function PetscSectionSetUp(arg1::PetscSection)
    ccall((:PetscSectionSetUp,petsc),PetscErrorCode,(PetscSection,),arg1)
end

function PetscSectionGetMaxDof(arg1::PetscSection,arg2::Ptr{PetscInt})
    ccall((:PetscSectionGetMaxDof,petsc),PetscErrorCode,(PetscSection,Ptr{PetscInt}),arg1,arg2)
end

function PetscSectionGetStorageSize(arg1::PetscSection,arg2::Ptr{PetscInt})
    ccall((:PetscSectionGetStorageSize,petsc),PetscErrorCode,(PetscSection,Ptr{PetscInt}),arg1,arg2)
end

function PetscSectionGetConstrainedStorageSize(arg1::PetscSection,arg2::Ptr{PetscInt})
    ccall((:PetscSectionGetConstrainedStorageSize,petsc),PetscErrorCode,(PetscSection,Ptr{PetscInt}),arg1,arg2)
end

function PetscSectionGetOffset(arg1::PetscSection,arg2::PetscInt,arg3::Ptr{PetscInt})
    ccall((:PetscSectionGetOffset,petsc),PetscErrorCode,(PetscSection,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscSectionSetOffset(arg1::PetscSection,arg2::PetscInt,arg3::PetscInt)
    ccall((:PetscSectionSetOffset,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt),arg1,arg2,arg3)
end

function PetscSectionGetFieldOffset(arg1::PetscSection,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscInt})
    ccall((:PetscSectionGetFieldOffset,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function PetscSectionSetFieldOffset(arg1::PetscSection,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt)
    ccall((:PetscSectionSetFieldOffset,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end

function PetscSectionGetOffsetRange(arg1::PetscSection,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt})
    ccall((:PetscSectionGetOffsetRange,petsc),PetscErrorCode,(PetscSection,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscSectionView(arg1::PetscSection,arg2::PetscViewer)
    ccall((:PetscSectionView,petsc),PetscErrorCode,(PetscSection,PetscViewer),arg1,arg2)
end

function PetscSectionDestroy(arg1::Ptr{PetscSection})
    ccall((:PetscSectionDestroy,petsc),PetscErrorCode,(Ptr{PetscSection},),arg1)
end

function PetscSectionCreateGlobalSection(arg1::PetscSection,arg2::PetscSF,arg3::PetscBool,arg4::Ptr{PetscSection})
    ccall((:PetscSectionCreateGlobalSection,petsc),PetscErrorCode,(PetscSection,PetscSF,PetscBool,Ptr{PetscSection}),arg1,arg2,arg3,arg4)
end

function PetscSectionCreateGlobalSectionCensored(arg1::PetscSection,arg2::PetscSF,arg3::PetscBool,arg4::PetscInt,arg5::Ptr{PetscInt},arg6::Ptr{PetscSection})
    ccall((:PetscSectionCreateGlobalSectionCensored,petsc),PetscErrorCode,(PetscSection,PetscSF,PetscBool,PetscInt,Ptr{PetscInt},Ptr{PetscSection}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscSectionCreateSubsection(arg1::PetscSection,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscSection})
    ccall((:PetscSectionCreateSubsection,petsc),PetscErrorCode,(PetscSection,PetscInt,Ptr{PetscInt},Ptr{PetscSection}),arg1,arg2,arg3,arg4)
end

function PetscSectionCreateSubmeshSection(arg1::PetscSection,arg2::IS,arg3::Ptr{PetscSection})
    ccall((:PetscSectionCreateSubmeshSection,petsc),PetscErrorCode,(PetscSection,IS,Ptr{PetscSection}),arg1,arg2,arg3)
end

function PetscSectionGetPointLayout(arg1::MPI_Comm,arg2::PetscSection,arg3::Ptr{PetscLayout})
    ccall((:PetscSectionGetPointLayout,petsc),PetscErrorCode,(MPI_Comm,PetscSection,Ptr{PetscLayout}),arg1,arg2,arg3)
end

function PetscSectionGetValueLayout(arg1::MPI_Comm,arg2::PetscSection,arg3::Ptr{PetscLayout})
    ccall((:PetscSectionGetValueLayout,petsc),PetscErrorCode,(MPI_Comm,PetscSection,Ptr{PetscLayout}),arg1,arg2,arg3)
end

function PetscSectionPermute(arg1::PetscSection,arg2::IS,arg3::Ptr{PetscSection})
    ccall((:PetscSectionPermute,petsc),PetscErrorCode,(PetscSection,IS,Ptr{PetscSection}),arg1,arg2,arg3)
end

function PetscSectionGetField(arg1::PetscSection,arg2::PetscInt,arg3::Ptr{PetscSection})
    ccall((:PetscSectionGetField,petsc),PetscErrorCode,(PetscSection,PetscInt,Ptr{PetscSection}),arg1,arg2,arg3)
end

function PetscSectionSetClosureIndex(arg1::PetscSection,arg2::PetscObject,arg3::PetscSection,arg4::IS)
    ccall((:PetscSectionSetClosureIndex,petsc),PetscErrorCode,(PetscSection,PetscObject,PetscSection,IS),arg1,arg2,arg3,arg4)
end

function PetscSectionGetClosureIndex(arg1::PetscSection,arg2::PetscObject,arg3::Ptr{PetscSection},arg4::Ptr{IS})
    ccall((:PetscSectionGetClosureIndex,petsc),PetscErrorCode,(PetscSection,PetscObject,Ptr{PetscSection},Ptr{IS}),arg1,arg2,arg3,arg4)
end

function PetscSFConvertPartition(arg1::PetscSF,arg2::PetscSection,arg3::IS,arg4::Ptr{ISLocalToGlobalMapping},arg5::Ptr{PetscSF})
    ccall((:PetscSFConvertPartition,petsc),PetscErrorCode,(PetscSF,PetscSection,IS,Ptr{ISLocalToGlobalMapping},Ptr{PetscSF}),arg1,arg2,arg3,arg4,arg5)
end

function PetscSFCreateRemoteOffsets(arg1::PetscSF,arg2::PetscSection,arg3::PetscSection,arg4::Ptr{Ptr{PetscInt}})
    ccall((:PetscSFCreateRemoteOffsets,petsc),PetscErrorCode,(PetscSF,PetscSection,PetscSection,Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4)
end

function PetscSFDistributeSection(arg1::PetscSF,arg2::PetscSection,arg3::Ptr{Ptr{PetscInt}},arg4::PetscSection)
    ccall((:PetscSFDistributeSection,petsc),PetscErrorCode,(PetscSF,PetscSection,Ptr{Ptr{PetscInt}},PetscSection),arg1,arg2,arg3,arg4)
end

function PetscSFCreateSectionSF(arg1::PetscSF,arg2::PetscSection,arg3::Ptr{PetscInt},arg4::PetscSection,arg5::Ptr{PetscSF})
    ccall((:PetscSFCreateSectionSF,petsc),PetscErrorCode,(PetscSF,PetscSection,Ptr{PetscInt},PetscSection,Ptr{PetscSF}),arg1,arg2,arg3,arg4,arg5)
end

function VecInitializePackage()
    ccall((:VecInitializePackage,petsc),PetscErrorCode,())
end

function VecFinalizePackage()
    ccall((:VecFinalizePackage,petsc),PetscErrorCode,())
end

function VecCreate(arg1::MPI_Comm,arg2::Ptr{Vec})
    ccall((:VecCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{Vec}),arg1,arg2)
end

function VecCreateSeq(arg1::MPI_Comm,arg2::PetscInt,arg3::Ptr{Vec})
    ccall((:VecCreateSeq,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{Vec}),arg1,arg2,arg3)
end

function VecCreateMPI(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{Vec})
    ccall((:VecCreateMPI,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,Ptr{Vec}),arg1,arg2,arg3,arg4)
end

function VecCreateSeqWithArray(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscScalar},arg5::Ptr{Vec})
    ccall((:VecCreateSeqWithArray,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,Ptr{PetscScalar},Ptr{Vec}),arg1,arg2,arg3,arg4,arg5)
end

function VecCreateMPIWithArray(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::Ptr{PetscScalar},arg6::Ptr{Vec})
    ccall((:VecCreateMPIWithArray,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,Ptr{PetscScalar},Ptr{Vec}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecCreateShared(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{Vec})
    ccall((:VecCreateShared,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,Ptr{Vec}),arg1,arg2,arg3,arg4)
end

function VecSetFromOptions(arg1::Vec)
    ccall((:VecSetFromOptions,petsc),PetscErrorCode,(Vec,),arg1)
end

function VecDestroy(arg1::Ptr{Vec})
    ccall((:VecDestroy,petsc),PetscErrorCode,(Ptr{Vec},),arg1)
end

function VecZeroEntries(arg1::Vec)
    ccall((:VecZeroEntries,petsc),PetscErrorCode,(Vec,),arg1)
end

function VecSetOptionsPrefix(arg1::Vec,arg2::Ptr{Uint8})
    ccall((:VecSetOptionsPrefix,petsc),PetscErrorCode,(Vec,Ptr{Uint8}),arg1,arg2)
end

function VecAppendOptionsPrefix(arg1::Vec,arg2::Ptr{Uint8})
    ccall((:VecAppendOptionsPrefix,petsc),PetscErrorCode,(Vec,Ptr{Uint8}),arg1,arg2)
end

function VecGetOptionsPrefix(arg1::Vec,arg2::Ptr{Ptr{Uint8}})
    ccall((:VecGetOptionsPrefix,petsc),PetscErrorCode,(Vec,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function VecSetSizes(arg1::Vec,arg2::PetscInt,arg3::PetscInt)
    ccall((:VecSetSizes,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt),arg1,arg2,arg3)
end

function VecDotNorm2(arg1::Vec,arg2::Vec,arg3::Ptr{PetscScalar},arg4::Ptr{Cint})
    ccall((:VecDotNorm2,petsc),PetscErrorCode,(Vec,Vec,Ptr{PetscScalar},Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function VecDot(arg1::Vec,arg2::Vec,arg3::Ptr{PetscScalar})
    ccall((:VecDot,petsc),PetscErrorCode,(Vec,Vec,Ptr{PetscScalar}),arg1,arg2,arg3)
end

function VecDotRealPart(arg1::Vec,arg2::Vec,arg3::Ptr{Cint})
    ccall((:VecDotRealPart,petsc),PetscErrorCode,(Vec,Vec,Ptr{Cint}),arg1,arg2,arg3)
end

function VecTDot(arg1::Vec,arg2::Vec,arg3::Ptr{PetscScalar})
    ccall((:VecTDot,petsc),PetscErrorCode,(Vec,Vec,Ptr{PetscScalar}),arg1,arg2,arg3)
end

function VecMDot(arg1::Vec,arg2::PetscInt,arg3::Ptr{Vec},arg4::Ptr{PetscScalar})
    ccall((:VecMDot,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{Vec},Ptr{PetscScalar}),arg1,arg2,arg3,arg4)
end

function VecMTDot(arg1::Vec,arg2::PetscInt,arg3::Ptr{Vec},arg4::Ptr{PetscScalar})
    ccall((:VecMTDot,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{Vec},Ptr{PetscScalar}),arg1,arg2,arg3,arg4)
end

function VecGetSubVector(arg1::Vec,arg2::IS,arg3::Ptr{Vec})
    ccall((:VecGetSubVector,petsc),PetscErrorCode,(Vec,IS,Ptr{Vec}),arg1,arg2,arg3)
end

function VecRestoreSubVector(arg1::Vec,arg2::IS,arg3::Ptr{Vec})
    ccall((:VecRestoreSubVector,petsc),PetscErrorCode,(Vec,IS,Ptr{Vec}),arg1,arg2,arg3)
end

function VecNorm(arg1::Vec,arg2::NormType,arg3::Ptr{Cint})
    ccall((:VecNorm,petsc),PetscErrorCode,(Vec,NormType,Ptr{Cint}),arg1,arg2,arg3)
end

function VecNormAvailable(arg1::Vec,arg2::NormType,arg3::Ptr{PetscBool},arg4::Ptr{Cint})
    ccall((:VecNormAvailable,petsc),PetscErrorCode,(Vec,NormType,Ptr{PetscBool},Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function VecNormalize(arg1::Vec,arg2::Ptr{Cint})
    ccall((:VecNormalize,petsc),PetscErrorCode,(Vec,Ptr{Cint}),arg1,arg2)
end

function VecSum(arg1::Vec,arg2::Ptr{PetscScalar})
    ccall((:VecSum,petsc),PetscErrorCode,(Vec,Ptr{PetscScalar}),arg1,arg2)
end

function VecMax(arg1::Vec,arg2::Ptr{PetscInt},arg3::Ptr{Cint})
    ccall((:VecMax,petsc),PetscErrorCode,(Vec,Ptr{PetscInt},Ptr{Cint}),arg1,arg2,arg3)
end

function VecMin(arg1::Vec,arg2::Ptr{PetscInt},arg3::Ptr{Cint})
    ccall((:VecMin,petsc),PetscErrorCode,(Vec,Ptr{PetscInt},Ptr{Cint}),arg1,arg2,arg3)
end

function VecScale(arg1::Vec,arg2::PetscScalar)
    ccall((:VecScale,petsc),PetscErrorCode,(Vec,PetscScalar),arg1,arg2)
end

function VecCopy(arg1::Vec,arg2::Vec)
    ccall((:VecCopy,petsc),PetscErrorCode,(Vec,Vec),arg1,arg2)
end

function VecSetRandom(arg1::Vec,arg2::PetscRandom)
    ccall((:VecSetRandom,petsc),PetscErrorCode,(Vec,PetscRandom),arg1,arg2)
end

function VecSet(arg1::Vec,arg2::PetscScalar)
    ccall((:VecSet,petsc),PetscErrorCode,(Vec,PetscScalar),arg1,arg2)
end

function VecSetInf(arg1::Vec)
    ccall((:VecSetInf,petsc),PetscErrorCode,(Vec,),arg1)
end

function VecSwap(arg1::Vec,arg2::Vec)
    ccall((:VecSwap,petsc),PetscErrorCode,(Vec,Vec),arg1,arg2)
end

function VecAXPY(arg1::Vec,arg2::PetscScalar,arg3::Vec)
    ccall((:VecAXPY,petsc),PetscErrorCode,(Vec,PetscScalar,Vec),arg1,arg2,arg3)
end

function VecAXPBY(arg1::Vec,arg2::PetscScalar,arg3::PetscScalar,arg4::Vec)
    ccall((:VecAXPBY,petsc),PetscErrorCode,(Vec,PetscScalar,PetscScalar,Vec),arg1,arg2,arg3,arg4)
end

function VecMAXPY(arg1::Vec,arg2::PetscInt,arg3::Ptr{PetscScalar},arg4::Ptr{Vec})
    ccall((:VecMAXPY,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{PetscScalar},Ptr{Vec}),arg1,arg2,arg3,arg4)
end

function VecAYPX(arg1::Vec,arg2::PetscScalar,arg3::Vec)
    ccall((:VecAYPX,petsc),PetscErrorCode,(Vec,PetscScalar,Vec),arg1,arg2,arg3)
end

function VecWAXPY(arg1::Vec,arg2::PetscScalar,arg3::Vec,arg4::Vec)
    ccall((:VecWAXPY,petsc),PetscErrorCode,(Vec,PetscScalar,Vec,Vec),arg1,arg2,arg3,arg4)
end

function VecAXPBYPCZ(arg1::Vec,arg2::PetscScalar,arg3::PetscScalar,arg4::PetscScalar,arg5::Vec,arg6::Vec)
    ccall((:VecAXPBYPCZ,petsc),PetscErrorCode,(Vec,PetscScalar,PetscScalar,PetscScalar,Vec,Vec),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecPointwiseMax(arg1::Vec,arg2::Vec,arg3::Vec)
    ccall((:VecPointwiseMax,petsc),PetscErrorCode,(Vec,Vec,Vec),arg1,arg2,arg3)
end

function VecPointwiseMaxAbs(arg1::Vec,arg2::Vec,arg3::Vec)
    ccall((:VecPointwiseMaxAbs,petsc),PetscErrorCode,(Vec,Vec,Vec),arg1,arg2,arg3)
end

function VecPointwiseMin(arg1::Vec,arg2::Vec,arg3::Vec)
    ccall((:VecPointwiseMin,petsc),PetscErrorCode,(Vec,Vec,Vec),arg1,arg2,arg3)
end

function VecPointwiseMult(arg1::Vec,arg2::Vec,arg3::Vec)
    ccall((:VecPointwiseMult,petsc),PetscErrorCode,(Vec,Vec,Vec),arg1,arg2,arg3)
end

function VecPointwiseDivide(arg1::Vec,arg2::Vec,arg3::Vec)
    ccall((:VecPointwiseDivide,petsc),PetscErrorCode,(Vec,Vec,Vec),arg1,arg2,arg3)
end

function VecMaxPointwiseDivide(arg1::Vec,arg2::Vec,arg3::Ptr{Cint})
    ccall((:VecMaxPointwiseDivide,petsc),PetscErrorCode,(Vec,Vec,Ptr{Cint}),arg1,arg2,arg3)
end

function VecShift(arg1::Vec,arg2::PetscScalar)
    ccall((:VecShift,petsc),PetscErrorCode,(Vec,PetscScalar),arg1,arg2)
end

function VecReciprocal(arg1::Vec)
    ccall((:VecReciprocal,petsc),PetscErrorCode,(Vec,),arg1)
end

function VecPermute(arg1::Vec,arg2::IS,arg3::PetscBool)
    ccall((:VecPermute,petsc),PetscErrorCode,(Vec,IS,PetscBool),arg1,arg2,arg3)
end

function VecSqrtAbs(arg1::Vec)
    ccall((:VecSqrtAbs,petsc),PetscErrorCode,(Vec,),arg1)
end

function VecLog(arg1::Vec)
    ccall((:VecLog,petsc),PetscErrorCode,(Vec,),arg1)
end

function VecExp(arg1::Vec)
    ccall((:VecExp,petsc),PetscErrorCode,(Vec,),arg1)
end

function VecAbs(arg1::Vec)
    ccall((:VecAbs,petsc),PetscErrorCode,(Vec,),arg1)
end

function VecDuplicate(arg1::Vec,arg2::Ptr{Vec})
    ccall((:VecDuplicate,petsc),PetscErrorCode,(Vec,Ptr{Vec}),arg1,arg2)
end

function VecDuplicateVecs(arg1::Vec,arg2::PetscInt,arg3::Ptr{Ptr{Vec}})
    ccall((:VecDuplicateVecs,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{Ptr{Vec}}),arg1,arg2,arg3)
end

function VecDestroyVecs(arg1::PetscInt,arg2::Ptr{Ptr{Vec}})
    ccall((:VecDestroyVecs,petsc),PetscErrorCode,(PetscInt,Ptr{Ptr{Vec}}),arg1,arg2)
end

function VecStrideNormAll(arg1::Vec,arg2::NormType,PetscReal::Ptr{Cint})
    ccall((:VecStrideNormAll,petsc),PetscErrorCode,(Vec,NormType,Ptr{Cint}),arg1,arg2,PetscReal)
end

function VecStrideMaxAll(arg1::Vec,arg2::Ptr{PetscInt},PetscReal::Ptr{Cint})
    ccall((:VecStrideMaxAll,petsc),PetscErrorCode,(Vec,Ptr{PetscInt},Ptr{Cint}),arg1,arg2,PetscReal)
end

function VecStrideMinAll(arg1::Vec,arg2::Ptr{PetscInt},PetscReal::Ptr{Cint})
    ccall((:VecStrideMinAll,petsc),PetscErrorCode,(Vec,Ptr{PetscInt},Ptr{Cint}),arg1,arg2,PetscReal)
end

function VecStrideScaleAll(arg1::Vec,arg2::Ptr{PetscScalar})
    ccall((:VecStrideScaleAll,petsc),PetscErrorCode,(Vec,Ptr{PetscScalar}),arg1,arg2)
end

function VecUniqueEntries(arg1::Vec,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{PetscScalar}})
    ccall((:VecUniqueEntries,petsc),PetscErrorCode,(Vec,Ptr{PetscInt},Ptr{Ptr{PetscScalar}}),arg1,arg2,arg3)
end

function VecStrideNorm(arg1::Vec,arg2::PetscInt,arg3::NormType,arg4::Ptr{Cint})
    ccall((:VecStrideNorm,petsc),PetscErrorCode,(Vec,PetscInt,NormType,Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function VecStrideMax(arg1::Vec,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{Cint})
    ccall((:VecStrideMax,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{PetscInt},Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function VecStrideMin(arg1::Vec,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{Cint})
    ccall((:VecStrideMin,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{PetscInt},Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function VecStrideScale(arg1::Vec,arg2::PetscInt,arg3::PetscScalar)
    ccall((:VecStrideScale,petsc),PetscErrorCode,(Vec,PetscInt,PetscScalar),arg1,arg2,arg3)
end

function VecStrideSet(arg1::Vec,arg2::PetscInt,arg3::PetscScalar)
    ccall((:VecStrideSet,petsc),PetscErrorCode,(Vec,PetscInt,PetscScalar),arg1,arg2,arg3)
end

function VecStrideGather(arg1::Vec,arg2::PetscInt,arg3::Vec,arg4::InsertMode)
    ccall((:VecStrideGather,petsc),PetscErrorCode,(Vec,PetscInt,Vec,InsertMode),arg1,arg2,arg3,arg4)
end

function VecStrideScatter(arg1::Vec,arg2::PetscInt,arg3::Vec,arg4::InsertMode)
    ccall((:VecStrideScatter,petsc),PetscErrorCode,(Vec,PetscInt,Vec,InsertMode),arg1,arg2,arg3,arg4)
end

function VecStrideGatherAll(arg1::Vec,arg2::Ptr{Vec},arg3::InsertMode)
    ccall((:VecStrideGatherAll,petsc),PetscErrorCode,(Vec,Ptr{Vec},InsertMode),arg1,arg2,arg3)
end

function VecStrideScatterAll(arg1::Ptr{Vec},arg2::Vec,arg3::InsertMode)
    ccall((:VecStrideScatterAll,petsc),PetscErrorCode,(Ptr{Vec},Vec,InsertMode),arg1,arg2,arg3)
end

function VecStrideSubSetScatter(arg1::Vec,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Vec,arg6::InsertMode)
    ccall((:VecStrideSubSetScatter,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Vec,InsertMode),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecStrideSubSetGather(arg1::Vec,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Vec,arg6::InsertMode)
    ccall((:VecStrideSubSetGather,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Vec,InsertMode),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecSetValues(arg1::Vec,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscScalar},arg5::InsertMode)
    ccall((:VecSetValues,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{PetscInt},Ptr{PetscScalar},InsertMode),arg1,arg2,arg3,arg4,arg5)
end

function VecGetValues(arg1::Vec,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscScalar})
    ccall((:VecGetValues,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{PetscInt},Ptr{PetscScalar}),arg1,arg2,arg3,arg4)
end

function VecAssemblyBegin(arg1::Vec)
    ccall((:VecAssemblyBegin,petsc),PetscErrorCode,(Vec,),arg1)
end

function VecAssemblyEnd(arg1::Vec)
    ccall((:VecAssemblyEnd,petsc),PetscErrorCode,(Vec,),arg1)
end

function VecStashSetInitialSize(arg1::Vec,arg2::PetscInt,arg3::PetscInt)
    ccall((:VecStashSetInitialSize,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt),arg1,arg2,arg3)
end

function VecStashView(arg1::Vec,arg2::PetscViewer)
    ccall((:VecStashView,petsc),PetscErrorCode,(Vec,PetscViewer),arg1,arg2)
end

function VecStashViewFromOptions(arg1::Vec,arg2::PetscObject,arg3::Ptr{Uint8})
    ccall((:VecStashViewFromOptions,petsc),PetscErrorCode,(Vec,PetscObject,Ptr{Uint8}),arg1,arg2,arg3)
end

function VecStashGetInfo(arg1::Vec,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{PetscInt})
    ccall((:VecStashGetInfo,petsc),PetscErrorCode,(Vec,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end

function VecGetBlockSize(arg1::Vec,arg2::Ptr{PetscInt})
    ccall((:VecGetBlockSize,petsc),PetscErrorCode,(Vec,Ptr{PetscInt}),arg1,arg2)
end

function VecSetValuesBlocked(arg1::Vec,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscScalar},arg5::InsertMode)
    ccall((:VecSetValuesBlocked,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{PetscInt},Ptr{PetscScalar},InsertMode),arg1,arg2,arg3,arg4,arg5)
end

function VecSetType(arg1::Vec,arg2::VecType)
    ccall((:VecSetType,petsc),PetscErrorCode,(Vec,VecType),arg1,arg2)
end

function VecGetType(arg1::Vec,arg2::Ptr{VecType})
    ccall((:VecGetType,petsc),PetscErrorCode,(Vec,Ptr{VecType}),arg1,arg2)
end

function VecRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:VecRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function VecScatterCreate(arg1::Vec,arg2::IS,arg3::Vec,arg4::IS,arg5::Ptr{VecScatter})
    ccall((:VecScatterCreate,petsc),PetscErrorCode,(Vec,IS,Vec,IS,Ptr{VecScatter}),arg1,arg2,arg3,arg4,arg5)
end

function VecScatterCreateEmpty(arg1::MPI_Comm,arg2::Ptr{VecScatter})
    ccall((:VecScatterCreateEmpty,petsc),PetscErrorCode,(MPI_Comm,Ptr{VecScatter}),arg1,arg2)
end

function VecScatterCreateLocal(arg1::VecScatter,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{PetscInt},arg6::PetscInt,arg7::Ptr{PetscInt},arg8::Ptr{PetscInt},arg9::Ptr{PetscInt},arg10::PetscInt)
    ccall((:VecScatterCreateLocal,petsc),PetscErrorCode,(VecScatter,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},PetscInt),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function VecScatterBegin(arg1::VecScatter,arg2::Vec,arg3::Vec,arg4::InsertMode,arg5::ScatterMode)
    ccall((:VecScatterBegin,petsc),PetscErrorCode,(VecScatter,Vec,Vec,InsertMode,ScatterMode),arg1,arg2,arg3,arg4,arg5)
end

function VecScatterEnd(arg1::VecScatter,arg2::Vec,arg3::Vec,arg4::InsertMode,arg5::ScatterMode)
    ccall((:VecScatterEnd,petsc),PetscErrorCode,(VecScatter,Vec,Vec,InsertMode,ScatterMode),arg1,arg2,arg3,arg4,arg5)
end

function VecScatterDestroy(arg1::Ptr{VecScatter})
    ccall((:VecScatterDestroy,petsc),PetscErrorCode,(Ptr{VecScatter},),arg1)
end

function VecScatterCopy(arg1::VecScatter,arg2::Ptr{VecScatter})
    ccall((:VecScatterCopy,petsc),PetscErrorCode,(VecScatter,Ptr{VecScatter}),arg1,arg2)
end

function VecScatterView(arg1::VecScatter,arg2::PetscViewer)
    ccall((:VecScatterView,petsc),PetscErrorCode,(VecScatter,PetscViewer),arg1,arg2)
end

function VecScatterGetMerged(arg1::VecScatter,arg2::Ptr{PetscBool})
    ccall((:VecScatterGetMerged,petsc),PetscErrorCode,(VecScatter,Ptr{PetscBool}),arg1,arg2)
end

function VecGetArray4d(arg1::Vec,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::PetscInt,arg7::PetscInt,arg8::PetscInt,arg9::PetscInt,arg10::Ptr{Ptr{Ptr{Ptr{Ptr{PetscScalar}}}}})
    ccall((:VecGetArray4d,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{Ptr{PetscScalar}}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function VecRestoreArray4d(arg1::Vec,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::PetscInt,arg7::PetscInt,arg8::PetscInt,arg9::PetscInt,arg10::Ptr{Ptr{Ptr{Ptr{Ptr{PetscScalar}}}}})
    ccall((:VecRestoreArray4d,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{Ptr{PetscScalar}}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function VecGetArray3d(arg1::Vec,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::PetscInt,arg7::PetscInt,arg8::Ptr{Ptr{Ptr{Ptr{PetscScalar}}}})
    ccall((:VecGetArray3d,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{PetscScalar}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function VecRestoreArray3d(arg1::Vec,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::PetscInt,arg7::PetscInt,arg8::Ptr{Ptr{Ptr{Ptr{PetscScalar}}}})
    ccall((:VecRestoreArray3d,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{PetscScalar}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function VecGetArray2d(arg1::Vec,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::Ptr{Ptr{Ptr{PetscScalar}}})
    ccall((:VecGetArray2d,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{PetscScalar}}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecRestoreArray2d(arg1::Vec,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::Ptr{Ptr{Ptr{PetscScalar}}})
    ccall((:VecRestoreArray2d,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{PetscScalar}}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecGetArray1d(arg1::Vec,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{Ptr{PetscScalar}})
    ccall((:VecGetArray1d,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,Ptr{Ptr{PetscScalar}}),arg1,arg2,arg3,arg4)
end

function VecRestoreArray1d(arg1::Vec,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{Ptr{PetscScalar}})
    ccall((:VecRestoreArray1d,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,Ptr{Ptr{PetscScalar}}),arg1,arg2,arg3,arg4)
end

function VecGetArray4dRead(arg1::Vec,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::PetscInt,arg7::PetscInt,arg8::PetscInt,arg9::PetscInt,arg10::Ptr{Ptr{Ptr{Ptr{Ptr{PetscScalar}}}}})
    ccall((:VecGetArray4dRead,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{Ptr{PetscScalar}}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function VecRestoreArray4dRead(arg1::Vec,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::PetscInt,arg7::PetscInt,arg8::PetscInt,arg9::PetscInt,arg10::Ptr{Ptr{Ptr{Ptr{Ptr{PetscScalar}}}}})
    ccall((:VecRestoreArray4dRead,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{Ptr{PetscScalar}}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function VecGetArray3dRead(arg1::Vec,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::PetscInt,arg7::PetscInt,arg8::Ptr{Ptr{Ptr{Ptr{PetscScalar}}}})
    ccall((:VecGetArray3dRead,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{PetscScalar}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function VecRestoreArray3dRead(arg1::Vec,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::PetscInt,arg7::PetscInt,arg8::Ptr{Ptr{Ptr{Ptr{PetscScalar}}}})
    ccall((:VecRestoreArray3dRead,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{PetscScalar}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function VecGetArray2dRead(arg1::Vec,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::Ptr{Ptr{Ptr{PetscScalar}}})
    ccall((:VecGetArray2dRead,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{PetscScalar}}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecRestoreArray2dRead(arg1::Vec,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::Ptr{Ptr{Ptr{PetscScalar}}})
    ccall((:VecRestoreArray2dRead,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{PetscScalar}}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecGetArray1dRead(arg1::Vec,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{Ptr{PetscScalar}})
    ccall((:VecGetArray1dRead,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,Ptr{Ptr{PetscScalar}}),arg1,arg2,arg3,arg4)
end

function VecRestoreArray1dRead(arg1::Vec,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{Ptr{PetscScalar}})
    ccall((:VecRestoreArray1dRead,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,Ptr{Ptr{PetscScalar}}),arg1,arg2,arg3,arg4)
end

function VecPlaceArray(arg1::Vec,arg2::Ptr{PetscScalar})
    ccall((:VecPlaceArray,petsc),PetscErrorCode,(Vec,Ptr{PetscScalar}),arg1,arg2)
end

function VecResetArray(arg1::Vec)
    ccall((:VecResetArray,petsc),PetscErrorCode,(Vec,),arg1)
end

function VecReplaceArray(arg1::Vec,arg2::Ptr{PetscScalar})
    ccall((:VecReplaceArray,petsc),PetscErrorCode,(Vec,Ptr{PetscScalar}),arg1,arg2)
end

function VecGetArrays(arg1::Ptr{Vec},arg2::PetscInt,arg3::Ptr{Ptr{Ptr{PetscScalar}}})
    ccall((:VecGetArrays,petsc),PetscErrorCode,(Ptr{Vec},PetscInt,Ptr{Ptr{Ptr{PetscScalar}}}),arg1,arg2,arg3)
end

function VecRestoreArrays(arg1::Ptr{Vec},arg2::PetscInt,arg3::Ptr{Ptr{Ptr{PetscScalar}}})
    ccall((:VecRestoreArrays,petsc),PetscErrorCode,(Ptr{Vec},PetscInt,Ptr{Ptr{Ptr{PetscScalar}}}),arg1,arg2,arg3)
end

function VecView(arg1::Vec,arg2::PetscViewer)
    ccall((:VecView,petsc),PetscErrorCode,(Vec,PetscViewer),arg1,arg2)
end

function VecEqual(arg1::Vec,arg2::Vec,arg3::Ptr{PetscBool})
    ccall((:VecEqual,petsc),PetscErrorCode,(Vec,Vec,Ptr{PetscBool}),arg1,arg2,arg3)
end

function VecLoad(arg1::Vec,arg2::PetscViewer)
    ccall((:VecLoad,petsc),PetscErrorCode,(Vec,PetscViewer),arg1,arg2)
end

function VecGetSize(arg1::Vec,arg2::Ptr{PetscInt})
    ccall((:VecGetSize,petsc),PetscErrorCode,(Vec,Ptr{PetscInt}),arg1,arg2)
end

function VecGetLocalSize(arg1::Vec,arg2::Ptr{PetscInt})
    ccall((:VecGetLocalSize,petsc),PetscErrorCode,(Vec,Ptr{PetscInt}),arg1,arg2)
end

function VecGetOwnershipRange(arg1::Vec,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt})
    ccall((:VecGetOwnershipRange,petsc),PetscErrorCode,(Vec,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function VecGetOwnershipRanges(arg1::Vec,arg2::Ptr{Ptr{PetscInt}})
    ccall((:VecGetOwnershipRanges,petsc),PetscErrorCode,(Vec,Ptr{Ptr{PetscInt}}),arg1,arg2)
end

function VecSetLocalToGlobalMapping(arg1::Vec,arg2::ISLocalToGlobalMapping)
    ccall((:VecSetLocalToGlobalMapping,petsc),PetscErrorCode,(Vec,ISLocalToGlobalMapping),arg1,arg2)
end

function VecSetValuesLocal(arg1::Vec,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscScalar},arg5::InsertMode)
    ccall((:VecSetValuesLocal,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{PetscInt},Ptr{PetscScalar},InsertMode),arg1,arg2,arg3,arg4,arg5)
end

function VecGetLocalToGlobalMapping(arg1::Vec,arg2::Ptr{ISLocalToGlobalMapping})
    ccall((:VecGetLocalToGlobalMapping,petsc),PetscErrorCode,(Vec,Ptr{ISLocalToGlobalMapping}),arg1,arg2)
end

function VecDotBegin(arg1::Vec,arg2::Vec,arg3::Ptr{PetscScalar})
    ccall((:VecDotBegin,petsc),PetscErrorCode,(Vec,Vec,Ptr{PetscScalar}),arg1,arg2,arg3)
end

function VecDotEnd(arg1::Vec,arg2::Vec,arg3::Ptr{PetscScalar})
    ccall((:VecDotEnd,petsc),PetscErrorCode,(Vec,Vec,Ptr{PetscScalar}),arg1,arg2,arg3)
end

function VecTDotBegin(arg1::Vec,arg2::Vec,arg3::Ptr{PetscScalar})
    ccall((:VecTDotBegin,petsc),PetscErrorCode,(Vec,Vec,Ptr{PetscScalar}),arg1,arg2,arg3)
end

function VecTDotEnd(arg1::Vec,arg2::Vec,arg3::Ptr{PetscScalar})
    ccall((:VecTDotEnd,petsc),PetscErrorCode,(Vec,Vec,Ptr{PetscScalar}),arg1,arg2,arg3)
end

function VecNormBegin(arg1::Vec,arg2::NormType,arg3::Ptr{Cint})
    ccall((:VecNormBegin,petsc),PetscErrorCode,(Vec,NormType,Ptr{Cint}),arg1,arg2,arg3)
end

function VecNormEnd(arg1::Vec,arg2::NormType,arg3::Ptr{Cint})
    ccall((:VecNormEnd,petsc),PetscErrorCode,(Vec,NormType,Ptr{Cint}),arg1,arg2,arg3)
end

function VecMDotBegin(arg1::Vec,arg2::PetscInt,arg3::Ptr{Vec},arg4::Ptr{PetscScalar})
    ccall((:VecMDotBegin,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{Vec},Ptr{PetscScalar}),arg1,arg2,arg3,arg4)
end

function VecMDotEnd(arg1::Vec,arg2::PetscInt,arg3::Ptr{Vec},arg4::Ptr{PetscScalar})
    ccall((:VecMDotEnd,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{Vec},Ptr{PetscScalar}),arg1,arg2,arg3,arg4)
end

function VecMTDotBegin(arg1::Vec,arg2::PetscInt,arg3::Ptr{Vec},arg4::Ptr{PetscScalar})
    ccall((:VecMTDotBegin,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{Vec},Ptr{PetscScalar}),arg1,arg2,arg3,arg4)
end

function VecMTDotEnd(arg1::Vec,arg2::PetscInt,arg3::Ptr{Vec},arg4::Ptr{PetscScalar})
    ccall((:VecMTDotEnd,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{Vec},Ptr{PetscScalar}),arg1,arg2,arg3,arg4)
end

function PetscCommSplitReductionBegin(arg1::MPI_Comm)
    ccall((:PetscCommSplitReductionBegin,petsc),PetscErrorCode,(MPI_Comm,),arg1)
end

function VecSetOption(arg1::Vec,arg2::VecOption,arg3::PetscBool)
    ccall((:VecSetOption,petsc),PetscErrorCode,(Vec,VecOption,PetscBool),arg1,arg2,arg3)
end

function VecGetArray(arg1::Vec,arg2::Ptr{Ptr{PetscScalar}})
    ccall((:VecGetArray,petsc),PetscErrorCode,(Vec,Ptr{Ptr{PetscScalar}}),arg1,arg2)
end

function VecGetArrayRead(arg1::Vec,arg2::Ptr{Ptr{PetscScalar}})
    ccall((:VecGetArrayRead,petsc),PetscErrorCode,(Vec,Ptr{Ptr{PetscScalar}}),arg1,arg2)
end

function VecRestoreArray(arg1::Vec,arg2::Ptr{Ptr{PetscScalar}})
    ccall((:VecRestoreArray,petsc),PetscErrorCode,(Vec,Ptr{Ptr{PetscScalar}}),arg1,arg2)
end

function VecRestoreArrayRead(arg1::Vec,arg2::Ptr{Ptr{PetscScalar}})
    ccall((:VecRestoreArrayRead,petsc),PetscErrorCode,(Vec,Ptr{Ptr{PetscScalar}}),arg1,arg2)
end

function VecGetLocalVector(arg1::Vec,arg2::Vec)
    ccall((:VecGetLocalVector,petsc),PetscErrorCode,(Vec,Vec),arg1,arg2)
end

function VecRestoreLocalVector(arg1::Vec,arg2::Vec)
    ccall((:VecRestoreLocalVector,petsc),PetscErrorCode,(Vec,Vec),arg1,arg2)
end

function VecGetLocalVectorRead(arg1::Vec,arg2::Vec)
    ccall((:VecGetLocalVectorRead,petsc),PetscErrorCode,(Vec,Vec),arg1,arg2)
end

function VecRestoreLocalVectorRead(arg1::Vec,arg2::Vec)
    ccall((:VecRestoreLocalVectorRead,petsc),PetscErrorCode,(Vec,Vec),arg1,arg2)
end

function VecContourScale(arg1::Vec,PetscReal::Cint,arg2::Cint)
    ccall((:VecContourScale,petsc),PetscErrorCode,(Vec,Cint,Cint),arg1,PetscReal,arg2)
end

function VecSetOperation(arg1::Vec,arg2::VecOperation,arg3::Ptr{Void})
    ccall((:VecSetOperation,petsc),PetscErrorCode,(Vec,VecOperation,Ptr{Void}),arg1,arg2,arg3)
end

function VecMPISetGhost(arg1::Vec,arg2::PetscInt,arg3::Ptr{PetscInt})
    ccall((:VecMPISetGhost,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function VecCreateGhost(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::Ptr{PetscInt},arg6::Ptr{Vec})
    ccall((:VecCreateGhost,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Vec}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecCreateGhostWithArray(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::Ptr{PetscInt},arg6::Ptr{PetscScalar},arg7::Ptr{Vec})
    ccall((:VecCreateGhostWithArray,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscScalar},Ptr{Vec}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function VecCreateGhostBlock(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::Ptr{PetscInt},arg7::Ptr{Vec})
    ccall((:VecCreateGhostBlock,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Vec}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function VecCreateGhostBlockWithArray(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::Ptr{PetscInt},arg7::Ptr{PetscScalar},arg8::Ptr{Vec})
    ccall((:VecCreateGhostBlockWithArray,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscScalar},Ptr{Vec}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function VecGhostGetLocalForm(arg1::Vec,arg2::Ptr{Vec})
    ccall((:VecGhostGetLocalForm,petsc),PetscErrorCode,(Vec,Ptr{Vec}),arg1,arg2)
end

function VecGhostRestoreLocalForm(arg1::Vec,arg2::Ptr{Vec})
    ccall((:VecGhostRestoreLocalForm,petsc),PetscErrorCode,(Vec,Ptr{Vec}),arg1,arg2)
end

function VecGhostIsLocalForm(arg1::Vec,arg2::Vec,arg3::Ptr{PetscBool})
    ccall((:VecGhostIsLocalForm,petsc),PetscErrorCode,(Vec,Vec,Ptr{PetscBool}),arg1,arg2,arg3)
end

function VecGhostUpdateBegin(arg1::Vec,arg2::InsertMode,arg3::ScatterMode)
    ccall((:VecGhostUpdateBegin,petsc),PetscErrorCode,(Vec,InsertMode,ScatterMode),arg1,arg2,arg3)
end

function VecGhostUpdateEnd(arg1::Vec,arg2::InsertMode,arg3::ScatterMode)
    ccall((:VecGhostUpdateEnd,petsc),PetscErrorCode,(Vec,InsertMode,ScatterMode),arg1,arg2,arg3)
end

function VecConjugate(arg1::Vec)
    ccall((:VecConjugate,petsc),PetscErrorCode,(Vec,),arg1)
end

function VecScatterCreateToAll(arg1::Vec,arg2::Ptr{VecScatter},arg3::Ptr{Vec})
    ccall((:VecScatterCreateToAll,petsc),PetscErrorCode,(Vec,Ptr{VecScatter},Ptr{Vec}),arg1,arg2,arg3)
end

function VecScatterCreateToZero(arg1::Vec,arg2::Ptr{VecScatter},arg3::Ptr{Vec})
    ccall((:VecScatterCreateToZero,petsc),PetscErrorCode,(Vec,Ptr{VecScatter},Ptr{Vec}),arg1,arg2,arg3)
end

function ISComplementVec(arg1::IS,arg2::Vec,arg3::Ptr{IS})
    ccall((:ISComplementVec,petsc),PetscErrorCode,(IS,Vec,Ptr{IS}),arg1,arg2,arg3)
end

function VecPow(arg1::Vec,arg2::PetscScalar)
    ccall((:VecPow,petsc),PetscErrorCode,(Vec,PetscScalar),arg1,arg2)
end

function VecMedian(arg1::Vec,arg2::Vec,arg3::Vec,arg4::Vec)
    ccall((:VecMedian,petsc),PetscErrorCode,(Vec,Vec,Vec,Vec),arg1,arg2,arg3,arg4)
end

function VecWhichBetween(arg1::Vec,arg2::Vec,arg3::Vec,arg4::Ptr{IS})
    ccall((:VecWhichBetween,petsc),PetscErrorCode,(Vec,Vec,Vec,Ptr{IS}),arg1,arg2,arg3,arg4)
end

function VecWhichBetweenOrEqual(arg1::Vec,arg2::Vec,arg3::Vec,arg4::Ptr{IS})
    ccall((:VecWhichBetweenOrEqual,petsc),PetscErrorCode,(Vec,Vec,Vec,Ptr{IS}),arg1,arg2,arg3,arg4)
end

function VecWhichGreaterThan(arg1::Vec,arg2::Vec,arg3::Ptr{IS})
    ccall((:VecWhichGreaterThan,petsc),PetscErrorCode,(Vec,Vec,Ptr{IS}),arg1,arg2,arg3)
end

function VecWhichLessThan(arg1::Vec,arg2::Vec,arg3::Ptr{IS})
    ccall((:VecWhichLessThan,petsc),PetscErrorCode,(Vec,Vec,Ptr{IS}),arg1,arg2,arg3)
end

function VecWhichEqual(arg1::Vec,arg2::Vec,arg3::Ptr{IS})
    ccall((:VecWhichEqual,petsc),PetscErrorCode,(Vec,Vec,Ptr{IS}),arg1,arg2,arg3)
end

function VecISAXPY(arg1::Vec,arg2::IS,arg3::PetscScalar,arg4::Vec)
    ccall((:VecISAXPY,petsc),PetscErrorCode,(Vec,IS,PetscScalar,Vec),arg1,arg2,arg3,arg4)
end

function VecISSet(arg1::Vec,arg2::IS,arg3::PetscScalar)
    ccall((:VecISSet,petsc),PetscErrorCode,(Vec,IS,PetscScalar),arg1,arg2,arg3)
end

function VecBoundGradientProjection(arg1::Vec,arg2::Vec,arg3::Vec,arg4::Vec,arg5::Vec)
    ccall((:VecBoundGradientProjection,petsc),PetscErrorCode,(Vec,Vec,Vec,Vec,Vec),arg1,arg2,arg3,arg4,arg5)
end

function VecStepBoundInfo(arg1::Vec,arg2::Vec,arg3::Vec,arg4::Vec,arg5::Ptr{Cint},arg6::Ptr{Cint},arg7::Ptr{Cint})
    ccall((:VecStepBoundInfo,petsc),PetscErrorCode,(Vec,Vec,Vec,Vec,Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function VecStepMax(arg1::Vec,arg2::Vec,arg3::Ptr{Cint})
    ccall((:VecStepMax,petsc),PetscErrorCode,(Vec,Vec,Ptr{Cint}),arg1,arg2,arg3)
end

function VecStepMaxBounded(arg1::Vec,arg2::Vec,arg3::Vec,arg4::Vec,arg5::Ptr{Cint})
    ccall((:VecStepMaxBounded,petsc),PetscErrorCode,(Vec,Vec,Vec,Vec,Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
end

function PetscViewerMathematicaGetVector(arg1::PetscViewer,arg2::Vec)
    ccall((:PetscViewerMathematicaGetVector,petsc),PetscErrorCode,(PetscViewer,Vec),arg1,arg2)
end

function PetscViewerMathematicaPutVector(arg1::PetscViewer,arg2::Vec)
    ccall((:PetscViewerMathematicaPutVector,petsc),PetscErrorCode,(PetscViewer,Vec),arg1,arg2)
end

function VecsDestroy(arg1::Vecs)
    ccall((:VecsDestroy,petsc),PetscErrorCode,(Vecs,),arg1)
end

function VecsCreateSeq(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{Vecs})
    ccall((:VecsCreateSeq,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,Ptr{Vecs}),arg1,arg2,arg3,arg4)
end

function VecsCreateSeqWithArray(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscScalar},arg5::Ptr{Vecs})
    ccall((:VecsCreateSeqWithArray,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,Ptr{PetscScalar},Ptr{Vecs}),arg1,arg2,arg3,arg4,arg5)
end

function VecsDuplicate(arg1::Vecs,arg2::Ptr{Vecs})
    ccall((:VecsDuplicate,petsc),PetscErrorCode,(Vecs,Ptr{Vecs}),arg1,arg2)
end

function VecNestGetSubVecs(arg1::Vec,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{Vec}})
    ccall((:VecNestGetSubVecs,petsc),PetscErrorCode,(Vec,Ptr{PetscInt},Ptr{Ptr{Vec}}),arg1,arg2,arg3)
end

function VecNestGetSubVec(arg1::Vec,arg2::PetscInt,arg3::Ptr{Vec})
    ccall((:VecNestGetSubVec,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{Vec}),arg1,arg2,arg3)
end

function VecNestSetSubVecs(arg1::Vec,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{Vec})
    ccall((:VecNestSetSubVecs,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{PetscInt},Ptr{Vec}),arg1,arg2,arg3,arg4)
end

function VecNestSetSubVec(arg1::Vec,arg2::PetscInt,arg3::Vec)
    ccall((:VecNestSetSubVec,petsc),PetscErrorCode,(Vec,PetscInt,Vec),arg1,arg2,arg3)
end

function VecCreateNest(arg1::MPI_Comm,arg2::PetscInt,arg3::Ptr{IS},arg4::Ptr{Vec},arg5::Ptr{Vec})
    ccall((:VecCreateNest,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{IS},Ptr{Vec},Ptr{Vec}),arg1,arg2,arg3,arg4,arg5)
end

function VecNestGetSize(arg1::Vec,arg2::Ptr{PetscInt})
    ccall((:VecNestGetSize,petsc),PetscErrorCode,(Vec,Ptr{PetscInt}),arg1,arg2)
end

function PetscOptionsGetVec(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Vec,arg4::Ptr{PetscBool})
    ccall((:PetscOptionsGetVec,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Vec,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function VecChop(arg1::Vec,PetscReal::Cint)
    ccall((:VecChop,petsc),PetscErrorCode,(Vec,Cint),arg1,PetscReal)
end

function VecGetLayout(arg1::Vec,arg2::Ptr{PetscLayout})
    ccall((:VecGetLayout,petsc),PetscErrorCode,(Vec,Ptr{PetscLayout}),arg1,arg2)
end

function VecSetLayout(arg1::Vec,arg2::PetscLayout)
    ccall((:VecSetLayout,petsc),PetscErrorCode,(Vec,PetscLayout),arg1,arg2)
end

function PetscSectionVecView(arg1::PetscSection,arg2::Vec,arg3::PetscViewer)
    ccall((:PetscSectionVecView,petsc),PetscErrorCode,(PetscSection,Vec,PetscViewer),arg1,arg2,arg3)
end

function VecGetValuesSection(arg1::Vec,arg2::PetscSection,arg3::PetscInt,arg4::Ptr{Ptr{PetscScalar}})
    ccall((:VecGetValuesSection,petsc),PetscErrorCode,(Vec,PetscSection,PetscInt,Ptr{Ptr{PetscScalar}}),arg1,arg2,arg3,arg4)
end

function VecSetValuesSection(arg1::Vec,arg2::PetscSection,arg3::PetscInt,arg4::Ptr{PetscScalar},arg5::InsertMode)
    ccall((:VecSetValuesSection,petsc),PetscErrorCode,(Vec,PetscSection,PetscInt,Ptr{PetscScalar},InsertMode),arg1,arg2,arg3,arg4,arg5)
end

function PetscSectionVecNorm(arg1::PetscSection,arg2::PetscSection,arg3::Vec,arg4::NormType,PetscReal::Ptr{Cint})
    ccall((:PetscSectionVecNorm,petsc),PetscErrorCode,(PetscSection,PetscSection,Vec,NormType,Ptr{Cint}),arg1,arg2,arg3,arg4,PetscReal)
end

function MatGetFactor(arg1::Mat,arg2::Ptr{Uint8},arg3::MatFactorType,arg4::Ptr{Mat})
    ccall((:MatGetFactor,petsc),PetscErrorCode,(Mat,Ptr{Uint8},MatFactorType,Ptr{Mat}),arg1,arg2,arg3,arg4)
end

function MatGetFactorAvailable(arg1::Mat,arg2::Ptr{Uint8},arg3::MatFactorType,arg4::Ptr{PetscBool})
    ccall((:MatGetFactorAvailable,petsc),PetscErrorCode,(Mat,Ptr{Uint8},MatFactorType,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function MatFactorGetSolverPackage(arg1::Mat,arg2::Ptr{Ptr{Uint8}})
    ccall((:MatFactorGetSolverPackage,petsc),PetscErrorCode,(Mat,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function MatGetFactorType(arg1::Mat,arg2::Ptr{MatFactorType})
    ccall((:MatGetFactorType,petsc),PetscErrorCode,(Mat,Ptr{MatFactorType}),arg1,arg2)
end

function MatSolverPackageRegister(arg1::Ptr{Uint8},arg2::MatType,arg3::MatFactorType,arg4::Ptr{Void})
    ccall((:MatSolverPackageRegister,petsc),PetscErrorCode,(Ptr{Uint8},MatType,MatFactorType,Ptr{Void}),arg1,arg2,arg3,arg4)
end

function MatSolverPackageGet(arg1::Ptr{Uint8},arg2::MatType,arg3::MatFactorType,arg4::Ptr{PetscBool},arg5::Ptr{PetscBool},arg6::Ptr{Ptr{Void}})
    ccall((:MatSolverPackageGet,petsc),PetscErrorCode,(Ptr{Uint8},MatType,MatFactorType,Ptr{PetscBool},Ptr{PetscBool},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatInitializePackage()
    ccall((:MatInitializePackage,petsc),PetscErrorCode,())
end

function MatCreate(arg1::MPI_Comm,arg2::Ptr{Mat})
    ccall((:MatCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{Mat}),arg1,arg2)
end

function MatSetSizes(arg1::Mat,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt)
    ccall((:MatSetSizes,petsc),PetscErrorCode,(Mat,PetscInt,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4,arg5)
end

function MatSetType(arg1::Mat,arg2::MatType)
    ccall((:MatSetType,petsc),PetscErrorCode,(Mat,MatType),arg1,arg2)
end

function MatSetFromOptions(arg1::Mat)
    ccall((:MatSetFromOptions,petsc),PetscErrorCode,(Mat,),arg1)
end

function MatRegisterBaseName(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{Uint8})
    ccall((:MatRegisterBaseName,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{Uint8}),arg1,arg2,arg3)
end

function MatSetOptionsPrefix(arg1::Mat,arg2::Ptr{Uint8})
    ccall((:MatSetOptionsPrefix,petsc),PetscErrorCode,(Mat,Ptr{Uint8}),arg1,arg2)
end

function MatAppendOptionsPrefix(arg1::Mat,arg2::Ptr{Uint8})
    ccall((:MatAppendOptionsPrefix,petsc),PetscErrorCode,(Mat,Ptr{Uint8}),arg1,arg2)
end

function MatGetOptionsPrefix(arg1::Mat,arg2::Ptr{Ptr{Uint8}})
    ccall((:MatGetOptionsPrefix,petsc),PetscErrorCode,(Mat,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function MatSetErrorIfFPE(arg1::Mat,arg2::PetscBool)
    ccall((:MatSetErrorIfFPE,petsc),PetscErrorCode,(Mat,PetscBool),arg1,arg2)
end

function MatCreateSeqDense(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscScalar},arg5::Ptr{Mat})
    ccall((:MatCreateSeqDense,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,Ptr{PetscScalar},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5)
end

function MatCreateDense(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::Ptr{PetscScalar},arg7::Ptr{Mat})
    ccall((:MatCreateDense,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscScalar},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateSeqAIJ(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::Ptr{PetscInt},arg6::Ptr{Mat})
    ccall((:MatCreateSeqAIJ,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatCreateAIJ(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::PetscInt,arg7::Ptr{PetscInt},arg8::PetscInt,arg9::Ptr{PetscInt},arg10::Ptr{Mat})
    ccall((:MatCreateAIJ,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function MatCreateMPIAIJWithArrays(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::Ptr{PetscInt},arg7::Ptr{PetscInt},arg8::Ptr{PetscScalar},arg9::Ptr{Mat})
    ccall((:MatCreateMPIAIJWithArrays,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscScalar},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function MatCreateMPIAIJWithSplitArrays(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::Ptr{PetscInt},arg7::Ptr{PetscInt},arg8::Ptr{PetscScalar},arg9::Ptr{PetscInt},arg10::Ptr{PetscInt},arg11::Ptr{PetscScalar},arg12::Ptr{Mat})
    ccall((:MatCreateMPIAIJWithSplitArrays,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscScalar},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscScalar},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12)
end

function MatCreateSeqBAIJ(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::Ptr{PetscInt},arg7::Ptr{Mat})
    ccall((:MatCreateSeqBAIJ,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateBAIJ(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::PetscInt,arg7::PetscInt,arg8::Ptr{PetscInt},arg9::PetscInt,arg10::Ptr{PetscInt},arg11::Ptr{Mat})
    ccall((:MatCreateBAIJ,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
end

function MatCreateMPIBAIJWithArrays(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::PetscInt,arg7::Ptr{PetscInt},arg8::Ptr{PetscInt},arg9::Ptr{PetscScalar},arg10::Ptr{Mat})
    ccall((:MatCreateMPIBAIJWithArrays,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscScalar},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function MatCreateMPIAdj(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscInt},arg5::Ptr{PetscInt},arg6::Ptr{PetscInt},arg7::Ptr{Mat})
    ccall((:MatCreateMPIAdj,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateSeqSBAIJ(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::Ptr{PetscInt},arg7::Ptr{Mat})
    ccall((:MatCreateSeqSBAIJ,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateSBAIJ(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::PetscInt,arg7::PetscInt,arg8::Ptr{PetscInt},arg9::PetscInt,arg10::Ptr{PetscInt},arg11::Ptr{Mat})
    ccall((:MatCreateSBAIJ,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
end

function MatCreateMPISBAIJWithArrays(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::PetscInt,arg7::Ptr{PetscInt},arg8::Ptr{PetscInt},arg9::Ptr{PetscScalar},arg10::Ptr{Mat})
    ccall((:MatCreateMPISBAIJWithArrays,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscScalar},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function MatSeqSBAIJSetPreallocationCSR(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{PetscScalar})
    ccall((:MatSeqSBAIJSetPreallocationCSR,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscScalar}),arg1,arg2,arg3,arg4,arg5)
end

function MatMPISBAIJSetPreallocationCSR(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{PetscScalar})
    ccall((:MatMPISBAIJSetPreallocationCSR,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscScalar}),arg1,arg2,arg3,arg4,arg5)
end

function MatXAIJSetPreallocation(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{PetscInt},arg6::Ptr{PetscInt})
    ccall((:MatXAIJSetPreallocation,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatCreateShell(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::Ptr{Void},arg7::Ptr{Mat})
    ccall((:MatCreateShell,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Void},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateNormal(arg1::Mat,arg2::Ptr{Mat})
    ccall((:MatCreateNormal,petsc),PetscErrorCode,(Mat,Ptr{Mat}),arg1,arg2)
end

function MatCreateLRC(arg1::Mat,arg2::Mat,arg3::Mat,arg4::Ptr{Mat})
    ccall((:MatCreateLRC,petsc),PetscErrorCode,(Mat,Mat,Mat,Ptr{Mat}),arg1,arg2,arg3,arg4)
end

function MatCreateIS(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::PetscInt,arg7::ISLocalToGlobalMapping,arg8::Ptr{Mat})
    ccall((:MatCreateIS,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,ISLocalToGlobalMapping,Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatCreateSeqAIJCRL(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::Ptr{PetscInt},arg6::Ptr{Mat})
    ccall((:MatCreateSeqAIJCRL,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatCreateMPIAIJCRL(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::Ptr{PetscInt},arg6::PetscInt,arg7::Ptr{PetscInt},arg8::Ptr{Mat})
    ccall((:MatCreateMPIAIJCRL,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatCreateSeqBSTRM(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::Ptr{PetscInt},arg7::Ptr{Mat})
    ccall((:MatCreateSeqBSTRM,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateMPIBSTRM(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::PetscInt,arg7::PetscInt,arg8::Ptr{PetscInt},arg9::PetscInt,arg10::Ptr{PetscInt},arg11::Ptr{Mat})
    ccall((:MatCreateMPIBSTRM,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
end

function MatCreateSeqSBSTRM(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::Ptr{PetscInt},arg7::Ptr{Mat})
    ccall((:MatCreateSeqSBSTRM,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateMPISBSTRM(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::PetscInt,arg7::PetscInt,arg8::Ptr{PetscInt},arg9::PetscInt,arg10::Ptr{PetscInt},arg11::Ptr{Mat})
    ccall((:MatCreateMPISBSTRM,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
end

function MatCreateScatter(arg1::MPI_Comm,arg2::VecScatter,arg3::Ptr{Mat})
    ccall((:MatCreateScatter,petsc),PetscErrorCode,(MPI_Comm,VecScatter,Ptr{Mat}),arg1,arg2,arg3)
end

function MatScatterSetVecScatter(arg1::Mat,arg2::VecScatter)
    ccall((:MatScatterSetVecScatter,petsc),PetscErrorCode,(Mat,VecScatter),arg1,arg2)
end

function MatScatterGetVecScatter(arg1::Mat,arg2::Ptr{VecScatter})
    ccall((:MatScatterGetVecScatter,petsc),PetscErrorCode,(Mat,Ptr{VecScatter}),arg1,arg2)
end

function MatCreateBlockMat(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::Ptr{PetscInt},arg7::Ptr{Mat})
    ccall((:MatCreateBlockMat,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCompositeAddMat(arg1::Mat,arg2::Mat)
    ccall((:MatCompositeAddMat,petsc),PetscErrorCode,(Mat,Mat),arg1,arg2)
end

function MatCompositeMerge(arg1::Mat)
    ccall((:MatCompositeMerge,petsc),PetscErrorCode,(Mat,),arg1)
end

function MatCreateComposite(arg1::MPI_Comm,arg2::PetscInt,arg3::Ptr{Mat},arg4::Ptr{Mat})
    ccall((:MatCreateComposite,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{Mat},Ptr{Mat}),arg1,arg2,arg3,arg4)
end

function MatCompositeSetType(arg1::Mat,arg2::MatCompositeType)
    ccall((:MatCompositeSetType,petsc),PetscErrorCode,(Mat,MatCompositeType),arg1,arg2)
end

function MatCreateFFT(arg1::MPI_Comm,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::MatType,arg5::Ptr{Mat})
    ccall((:MatCreateFFT,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{PetscInt},MatType,Ptr{Mat}),arg1,arg2,arg3,arg4,arg5)
end

function MatCreateSeqCUFFT(arg1::MPI_Comm,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{Mat})
    ccall((:MatCreateSeqCUFFT,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4)
end

function MatCreateTranspose(arg1::Mat,arg2::Ptr{Mat})
    ccall((:MatCreateTranspose,petsc),PetscErrorCode,(Mat,Ptr{Mat}),arg1,arg2)
end

function MatCreateHermitianTranspose(arg1::Mat,arg2::Ptr{Mat})
    ccall((:MatCreateHermitianTranspose,petsc),PetscErrorCode,(Mat,Ptr{Mat}),arg1,arg2)
end

function MatCreateSubMatrix(arg1::Mat,arg2::IS,arg3::IS,arg4::Ptr{Mat})
    ccall((:MatCreateSubMatrix,petsc),PetscErrorCode,(Mat,IS,IS,Ptr{Mat}),arg1,arg2,arg3,arg4)
end

function MatSubMatrixUpdate(arg1::Mat,arg2::Mat,arg3::IS,arg4::IS)
    ccall((:MatSubMatrixUpdate,petsc),PetscErrorCode,(Mat,Mat,IS,IS),arg1,arg2,arg3,arg4)
end

function MatCreateLocalRef(arg1::Mat,arg2::IS,arg3::IS,arg4::Ptr{Mat})
    ccall((:MatCreateLocalRef,petsc),PetscErrorCode,(Mat,IS,IS,Ptr{Mat}),arg1,arg2,arg3,arg4)
end

function MatPythonSetType(arg1::Mat,arg2::Ptr{Uint8})
    ccall((:MatPythonSetType,petsc),PetscErrorCode,(Mat,Ptr{Uint8}),arg1,arg2)
end

function MatSetUp(arg1::Mat)
    ccall((:MatSetUp,petsc),PetscErrorCode,(Mat,),arg1)
end

function MatDestroy(arg1::Ptr{Mat})
    ccall((:MatDestroy,petsc),PetscErrorCode,(Ptr{Mat},),arg1)
end

function MatGetNonzeroState(arg1::Mat,arg2::Ptr{PetscObjectState})
    ccall((:MatGetNonzeroState,petsc),PetscErrorCode,(Mat,Ptr{PetscObjectState}),arg1,arg2)
end

function MatConjugate(arg1::Mat)
    ccall((:MatConjugate,petsc),PetscErrorCode,(Mat,),arg1)
end

function MatRealPart(arg1::Mat)
    ccall((:MatRealPart,petsc),PetscErrorCode,(Mat,),arg1)
end

function MatImaginaryPart(arg1::Mat)
    ccall((:MatImaginaryPart,petsc),PetscErrorCode,(Mat,),arg1)
end

function MatGetDiagonalBlock(arg1::Mat,arg2::Ptr{Mat})
    ccall((:MatGetDiagonalBlock,petsc),PetscErrorCode,(Mat,Ptr{Mat}),arg1,arg2)
end

function MatGetTrace(arg1::Mat,arg2::Ptr{PetscScalar})
    ccall((:MatGetTrace,petsc),PetscErrorCode,(Mat,Ptr{PetscScalar}),arg1,arg2)
end

function MatInvertBlockDiagonal(arg1::Mat,arg2::Ptr{Ptr{PetscScalar}})
    ccall((:MatInvertBlockDiagonal,petsc),PetscErrorCode,(Mat,Ptr{Ptr{PetscScalar}}),arg1,arg2)
end

function MatSetValues(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::PetscInt,arg5::Ptr{PetscInt},arg6::Ptr{PetscScalar},arg7::InsertMode)
    ccall((:MatSetValues,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{PetscScalar},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatSetValuesBlocked(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::PetscInt,arg5::Ptr{PetscInt},arg6::Ptr{PetscScalar},arg7::InsertMode)
    ccall((:MatSetValuesBlocked,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{PetscScalar},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatSetValuesRow(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscScalar})
    ccall((:MatSetValuesRow,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscScalar}),arg1,arg2,arg3)
end

function MatSetValuesRowLocal(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscScalar})
    ccall((:MatSetValuesRowLocal,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscScalar}),arg1,arg2,arg3)
end

function MatSetValuesBatch(arg1::Mat,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscInt},arg5::Ptr{PetscScalar})
    ccall((:MatSetValuesBatch,petsc),PetscErrorCode,(Mat,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscScalar}),arg1,arg2,arg3,arg4,arg5)
end

function MatSetRandom(arg1::Mat,arg2::PetscRandom)
    ccall((:MatSetRandom,petsc),PetscErrorCode,(Mat,PetscRandom),arg1,arg2)
end

function MatSetValuesStencil(arg1::Mat,arg2::PetscInt,arg3::Ptr{MatStencil},arg4::PetscInt,arg5::Ptr{MatStencil},arg6::Ptr{PetscScalar},arg7::InsertMode)
    ccall((:MatSetValuesStencil,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{MatStencil},PetscInt,Ptr{MatStencil},Ptr{PetscScalar},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatSetValuesBlockedStencil(arg1::Mat,arg2::PetscInt,arg3::Ptr{MatStencil},arg4::PetscInt,arg5::Ptr{MatStencil},arg6::Ptr{PetscScalar},arg7::InsertMode)
    ccall((:MatSetValuesBlockedStencil,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{MatStencil},PetscInt,Ptr{MatStencil},Ptr{PetscScalar},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatSetStencil(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::PetscInt)
    ccall((:MatSetStencil,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{PetscInt},PetscInt),arg1,arg2,arg3,arg4,arg5)
end

function MatSetColoring(arg1::Mat,arg2::ISColoring)
    ccall((:MatSetColoring,petsc),PetscErrorCode,(Mat,ISColoring),arg1,arg2)
end

function MatSetValuesAdifor(arg1::Mat,arg2::PetscInt,arg3::Ptr{Void})
    ccall((:MatSetValuesAdifor,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{Void}),arg1,arg2,arg3)
end

function MatAssemblyBegin(arg1::Mat,arg2::MatAssemblyType)
    ccall((:MatAssemblyBegin,petsc),PetscErrorCode,(Mat,MatAssemblyType),arg1,arg2)
end

function MatAssemblyEnd(arg1::Mat,arg2::MatAssemblyType)
    ccall((:MatAssemblyEnd,petsc),PetscErrorCode,(Mat,MatAssemblyType),arg1,arg2)
end

function MatAssembled(arg1::Mat,arg2::Ptr{PetscBool})
    ccall((:MatAssembled,petsc),PetscErrorCode,(Mat,Ptr{PetscBool}),arg1,arg2)
end

function MatSetOption(arg1::Mat,arg2::MatOption,arg3::PetscBool)
    ccall((:MatSetOption,petsc),PetscErrorCode,(Mat,MatOption,PetscBool),arg1,arg2,arg3)
end

function MatGetOption(arg1::Mat,arg2::MatOption,arg3::Ptr{PetscBool})
    ccall((:MatGetOption,petsc),PetscErrorCode,(Mat,MatOption,Ptr{PetscBool}),arg1,arg2,arg3)
end

function MatGetType(arg1::Mat,arg2::Ptr{MatType})
    ccall((:MatGetType,petsc),PetscErrorCode,(Mat,Ptr{MatType}),arg1,arg2)
end

function MatGetValues(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::PetscInt,arg5::Ptr{PetscInt},arg6::Ptr{PetscScalar})
    ccall((:MatGetValues,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{PetscScalar}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatGetRow(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{Ptr{PetscInt}},arg5::Ptr{Ptr{PetscScalar}})
    ccall((:MatGetRow,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscScalar}}),arg1,arg2,arg3,arg4,arg5)
end

function MatRestoreRow(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{Ptr{PetscInt}},arg5::Ptr{Ptr{PetscScalar}})
    ccall((:MatRestoreRow,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscScalar}}),arg1,arg2,arg3,arg4,arg5)
end

function MatGetRowUpperTriangular(arg1::Mat)
    ccall((:MatGetRowUpperTriangular,petsc),PetscErrorCode,(Mat,),arg1)
end

function MatRestoreRowUpperTriangular(arg1::Mat)
    ccall((:MatRestoreRowUpperTriangular,petsc),PetscErrorCode,(Mat,),arg1)
end

function MatGetColumn(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{Ptr{PetscInt}},arg5::Ptr{Ptr{PetscScalar}})
    ccall((:MatGetColumn,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscScalar}}),arg1,arg2,arg3,arg4,arg5)
end

function MatRestoreColumn(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{Ptr{PetscInt}},arg5::Ptr{Ptr{PetscScalar}})
    ccall((:MatRestoreColumn,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscScalar}}),arg1,arg2,arg3,arg4,arg5)
end

function MatGetColumnVector(arg1::Mat,arg2::Vec,arg3::PetscInt)
    ccall((:MatGetColumnVector,petsc),PetscErrorCode,(Mat,Vec,PetscInt),arg1,arg2,arg3)
end

function MatSeqAIJGetArray(arg1::Mat,arg2::Ptr{Ptr{PetscScalar}})
    ccall((:MatSeqAIJGetArray,petsc),PetscErrorCode,(Mat,Ptr{Ptr{PetscScalar}}),arg1,arg2)
end

function MatSeqAIJRestoreArray(arg1::Mat,arg2::Ptr{Ptr{PetscScalar}})
    ccall((:MatSeqAIJRestoreArray,petsc),PetscErrorCode,(Mat,Ptr{Ptr{PetscScalar}}),arg1,arg2)
end

function MatSeqAIJGetMaxRowNonzeros(arg1::Mat,arg2::Ptr{PetscInt})
    ccall((:MatSeqAIJGetMaxRowNonzeros,petsc),PetscErrorCode,(Mat,Ptr{PetscInt}),arg1,arg2)
end

function MatSeqAIJSetValuesLocalFast(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::PetscInt,arg5::Ptr{PetscInt},arg6::Ptr{PetscScalar},arg7::InsertMode)
    ccall((:MatSeqAIJSetValuesLocalFast,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{PetscScalar},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatDenseGetArray(arg1::Mat,arg2::Ptr{Ptr{PetscScalar}})
    ccall((:MatDenseGetArray,petsc),PetscErrorCode,(Mat,Ptr{Ptr{PetscScalar}}),arg1,arg2)
end

function MatDenseRestoreArray(arg1::Mat,arg2::Ptr{Ptr{PetscScalar}})
    ccall((:MatDenseRestoreArray,petsc),PetscErrorCode,(Mat,Ptr{Ptr{PetscScalar}}),arg1,arg2)
end

function MatGetBlockSize(arg1::Mat,arg2::Ptr{PetscInt})
    ccall((:MatGetBlockSize,petsc),PetscErrorCode,(Mat,Ptr{PetscInt}),arg1,arg2)
end

function MatSetBlockSize(arg1::Mat,arg2::PetscInt)
    ccall((:MatSetBlockSize,petsc),PetscErrorCode,(Mat,PetscInt),arg1,arg2)
end

function MatGetBlockSizes(arg1::Mat,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt})
    ccall((:MatGetBlockSizes,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatSetBlockSizes(arg1::Mat,arg2::PetscInt,arg3::PetscInt)
    ccall((:MatSetBlockSizes,petsc),PetscErrorCode,(Mat,PetscInt,PetscInt),arg1,arg2,arg3)
end

function MatSetBlockSizesFromMats(arg1::Mat,arg2::Mat,arg3::Mat)
    ccall((:MatSetBlockSizesFromMats,petsc),PetscErrorCode,(Mat,Mat,Mat),arg1,arg2,arg3)
end

function MatSetNThreads(arg1::Mat,arg2::PetscInt)
    ccall((:MatSetNThreads,petsc),PetscErrorCode,(Mat,PetscInt),arg1,arg2)
end

function MatGetNThreads(arg1::Mat,arg2::Ptr{PetscInt})
    ccall((:MatGetNThreads,petsc),PetscErrorCode,(Mat,Ptr{PetscInt}),arg1,arg2)
end

function MatMult(arg1::Mat,arg2::Vec,arg3::Vec)
    ccall((:MatMult,petsc),PetscErrorCode,(Mat,Vec,Vec),arg1,arg2,arg3)
end

function MatMultDiagonalBlock(arg1::Mat,arg2::Vec,arg3::Vec)
    ccall((:MatMultDiagonalBlock,petsc),PetscErrorCode,(Mat,Vec,Vec),arg1,arg2,arg3)
end

function MatMultAdd(arg1::Mat,arg2::Vec,arg3::Vec,arg4::Vec)
    ccall((:MatMultAdd,petsc),PetscErrorCode,(Mat,Vec,Vec,Vec),arg1,arg2,arg3,arg4)
end

function MatMultTranspose(arg1::Mat,arg2::Vec,arg3::Vec)
    ccall((:MatMultTranspose,petsc),PetscErrorCode,(Mat,Vec,Vec),arg1,arg2,arg3)
end

function MatMultHermitianTranspose(arg1::Mat,arg2::Vec,arg3::Vec)
    ccall((:MatMultHermitianTranspose,petsc),PetscErrorCode,(Mat,Vec,Vec),arg1,arg2,arg3)
end

function MatIsTranspose(arg1::Mat,arg2::Mat,PetscReal::Cint,arg3::Ptr{PetscBool})
    ccall((:MatIsTranspose,petsc),PetscErrorCode,(Mat,Mat,Cint,Ptr{PetscBool}),arg1,arg2,PetscReal,arg3)
end

function MatIsHermitianTranspose(arg1::Mat,arg2::Mat,PetscReal::Cint,arg3::Ptr{PetscBool})
    ccall((:MatIsHermitianTranspose,petsc),PetscErrorCode,(Mat,Mat,Cint,Ptr{PetscBool}),arg1,arg2,PetscReal,arg3)
end

function MatMultTransposeAdd(arg1::Mat,arg2::Vec,arg3::Vec,arg4::Vec)
    ccall((:MatMultTransposeAdd,petsc),PetscErrorCode,(Mat,Vec,Vec,Vec),arg1,arg2,arg3,arg4)
end

function MatMultHermitianTransposeAdd(arg1::Mat,arg2::Vec,arg3::Vec,arg4::Vec)
    ccall((:MatMultHermitianTransposeAdd,petsc),PetscErrorCode,(Mat,Vec,Vec,Vec),arg1,arg2,arg3,arg4)
end

function MatMultConstrained(arg1::Mat,arg2::Vec,arg3::Vec)
    ccall((:MatMultConstrained,petsc),PetscErrorCode,(Mat,Vec,Vec),arg1,arg2,arg3)
end

function MatMultTransposeConstrained(arg1::Mat,arg2::Vec,arg3::Vec)
    ccall((:MatMultTransposeConstrained,petsc),PetscErrorCode,(Mat,Vec,Vec),arg1,arg2,arg3)
end

function MatMatSolve(arg1::Mat,arg2::Mat,arg3::Mat)
    ccall((:MatMatSolve,petsc),PetscErrorCode,(Mat,Mat,Mat),arg1,arg2,arg3)
end

function MatResidual(arg1::Mat,arg2::Vec,arg3::Vec,arg4::Vec)
    ccall((:MatResidual,petsc),PetscErrorCode,(Mat,Vec,Vec,Vec),arg1,arg2,arg3,arg4)
end

function MatConvert(arg1::Mat,arg2::MatType,arg3::MatReuse,arg4::Ptr{Mat})
    ccall((:MatConvert,petsc),PetscErrorCode,(Mat,MatType,MatReuse,Ptr{Mat}),arg1,arg2,arg3,arg4)
end

function MatDuplicate(arg1::Mat,arg2::MatDuplicateOption,arg3::Ptr{Mat})
    ccall((:MatDuplicate,petsc),PetscErrorCode,(Mat,MatDuplicateOption,Ptr{Mat}),arg1,arg2,arg3)
end

function MatCopy(arg1::Mat,arg2::Mat,arg3::MatStructure)
    ccall((:MatCopy,petsc),PetscErrorCode,(Mat,Mat,MatStructure),arg1,arg2,arg3)
end

function MatView(arg1::Mat,arg2::PetscViewer)
    ccall((:MatView,petsc),PetscErrorCode,(Mat,PetscViewer),arg1,arg2)
end

function MatIsSymmetric(arg1::Mat,PetscReal::Cint,arg2::Ptr{PetscBool})
    ccall((:MatIsSymmetric,petsc),PetscErrorCode,(Mat,Cint,Ptr{PetscBool}),arg1,PetscReal,arg2)
end

function MatIsStructurallySymmetric(arg1::Mat,arg2::Ptr{PetscBool})
    ccall((:MatIsStructurallySymmetric,petsc),PetscErrorCode,(Mat,Ptr{PetscBool}),arg1,arg2)
end

function MatIsHermitian(arg1::Mat,PetscReal::Cint,arg2::Ptr{PetscBool})
    ccall((:MatIsHermitian,petsc),PetscErrorCode,(Mat,Cint,Ptr{PetscBool}),arg1,PetscReal,arg2)
end

function MatIsSymmetricKnown(arg1::Mat,arg2::Ptr{PetscBool},arg3::Ptr{PetscBool})
    ccall((:MatIsSymmetricKnown,petsc),PetscErrorCode,(Mat,Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3)
end

function MatIsHermitianKnown(arg1::Mat,arg2::Ptr{PetscBool},arg3::Ptr{PetscBool})
    ccall((:MatIsHermitianKnown,petsc),PetscErrorCode,(Mat,Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3)
end

function MatMissingDiagonal(arg1::Mat,arg2::Ptr{PetscBool},arg3::Ptr{PetscInt})
    ccall((:MatMissingDiagonal,petsc),PetscErrorCode,(Mat,Ptr{PetscBool},Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatLoad(arg1::Mat,arg2::PetscViewer)
    ccall((:MatLoad,petsc),PetscErrorCode,(Mat,PetscViewer),arg1,arg2)
end

function MatGetRowIJ(arg1::Mat,arg2::PetscInt,arg3::PetscBool,arg4::PetscBool,arg5::Ptr{PetscInt},arg6::Ptr{Ptr{PetscInt}},arg7::Ptr{Ptr{PetscInt}},arg8::Ptr{PetscBool})
    ccall((:MatGetRowIJ,petsc),PetscErrorCode,(Mat,PetscInt,PetscBool,PetscBool,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatRestoreRowIJ(arg1::Mat,arg2::PetscInt,arg3::PetscBool,arg4::PetscBool,arg5::Ptr{PetscInt},arg6::Ptr{Ptr{PetscInt}},arg7::Ptr{Ptr{PetscInt}},arg8::Ptr{PetscBool})
    ccall((:MatRestoreRowIJ,petsc),PetscErrorCode,(Mat,PetscInt,PetscBool,PetscBool,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatGetColumnIJ(arg1::Mat,arg2::PetscInt,arg3::PetscBool,arg4::PetscBool,arg5::Ptr{PetscInt},arg6::Ptr{Ptr{PetscInt}},arg7::Ptr{Ptr{PetscInt}},arg8::Ptr{PetscBool})
    ccall((:MatGetColumnIJ,petsc),PetscErrorCode,(Mat,PetscInt,PetscBool,PetscBool,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatRestoreColumnIJ(arg1::Mat,arg2::PetscInt,arg3::PetscBool,arg4::PetscBool,arg5::Ptr{PetscInt},arg6::Ptr{Ptr{PetscInt}},arg7::Ptr{Ptr{PetscInt}},arg8::Ptr{PetscBool})
    ccall((:MatRestoreColumnIJ,petsc),PetscErrorCode,(Mat,PetscInt,PetscBool,PetscBool,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatGetInfo(arg1::Mat,arg2::MatInfoType,arg3::Ptr{MatInfo})
    ccall((:MatGetInfo,petsc),PetscErrorCode,(Mat,MatInfoType,Ptr{MatInfo}),arg1,arg2,arg3)
end

function MatGetDiagonal(arg1::Mat,arg2::Vec)
    ccall((:MatGetDiagonal,petsc),PetscErrorCode,(Mat,Vec),arg1,arg2)
end

function MatGetRowMax(arg1::Mat,arg2::Vec,arg3::Ptr{PetscInt})
    ccall((:MatGetRowMax,petsc),PetscErrorCode,(Mat,Vec,Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatGetRowMin(arg1::Mat,arg2::Vec,arg3::Ptr{PetscInt})
    ccall((:MatGetRowMin,petsc),PetscErrorCode,(Mat,Vec,Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatGetRowMaxAbs(arg1::Mat,arg2::Vec,arg3::Ptr{PetscInt})
    ccall((:MatGetRowMaxAbs,petsc),PetscErrorCode,(Mat,Vec,Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatGetRowMinAbs(arg1::Mat,arg2::Vec,arg3::Ptr{PetscInt})
    ccall((:MatGetRowMinAbs,petsc),PetscErrorCode,(Mat,Vec,Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatGetRowSum(arg1::Mat,arg2::Vec)
    ccall((:MatGetRowSum,petsc),PetscErrorCode,(Mat,Vec),arg1,arg2)
end

function MatTranspose(arg1::Mat,arg2::MatReuse,arg3::Ptr{Mat})
    ccall((:MatTranspose,petsc),PetscErrorCode,(Mat,MatReuse,Ptr{Mat}),arg1,arg2,arg3)
end

function MatHermitianTranspose(arg1::Mat,arg2::MatReuse,arg3::Ptr{Mat})
    ccall((:MatHermitianTranspose,petsc),PetscErrorCode,(Mat,MatReuse,Ptr{Mat}),arg1,arg2,arg3)
end

function MatPermute(arg1::Mat,arg2::IS,arg3::IS,arg4::Ptr{Mat})
    ccall((:MatPermute,petsc),PetscErrorCode,(Mat,IS,IS,Ptr{Mat}),arg1,arg2,arg3,arg4)
end

function MatDiagonalScale(arg1::Mat,arg2::Vec,arg3::Vec)
    ccall((:MatDiagonalScale,petsc),PetscErrorCode,(Mat,Vec,Vec),arg1,arg2,arg3)
end

function MatDiagonalSet(arg1::Mat,arg2::Vec,arg3::InsertMode)
    ccall((:MatDiagonalSet,petsc),PetscErrorCode,(Mat,Vec,InsertMode),arg1,arg2,arg3)
end

function MatEqual(arg1::Mat,arg2::Mat,arg3::Ptr{PetscBool})
    ccall((:MatEqual,petsc),PetscErrorCode,(Mat,Mat,Ptr{PetscBool}),arg1,arg2,arg3)
end

function MatMultEqual(arg1::Mat,arg2::Mat,arg3::PetscInt,arg4::Ptr{PetscBool})
    ccall((:MatMultEqual,petsc),PetscErrorCode,(Mat,Mat,PetscInt,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function MatMultAddEqual(arg1::Mat,arg2::Mat,arg3::PetscInt,arg4::Ptr{PetscBool})
    ccall((:MatMultAddEqual,petsc),PetscErrorCode,(Mat,Mat,PetscInt,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function MatMultTransposeEqual(arg1::Mat,arg2::Mat,arg3::PetscInt,arg4::Ptr{PetscBool})
    ccall((:MatMultTransposeEqual,petsc),PetscErrorCode,(Mat,Mat,PetscInt,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function MatMultTransposeAddEqual(arg1::Mat,arg2::Mat,arg3::PetscInt,arg4::Ptr{PetscBool})
    ccall((:MatMultTransposeAddEqual,petsc),PetscErrorCode,(Mat,Mat,PetscInt,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function MatNorm(arg1::Mat,arg2::NormType,arg3::Ptr{Cint})
    ccall((:MatNorm,petsc),PetscErrorCode,(Mat,NormType,Ptr{Cint}),arg1,arg2,arg3)
end

function MatGetColumnNorms(arg1::Mat,arg2::NormType,arg3::Ptr{Cint})
    ccall((:MatGetColumnNorms,petsc),PetscErrorCode,(Mat,NormType,Ptr{Cint}),arg1,arg2,arg3)
end

function MatZeroEntries(arg1::Mat)
    ccall((:MatZeroEntries,petsc),PetscErrorCode,(Mat,),arg1)
end

function MatZeroRows(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::PetscScalar,arg5::Vec,arg6::Vec)
    ccall((:MatZeroRows,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscScalar,Vec,Vec),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatZeroRowsIS(arg1::Mat,arg2::IS,arg3::PetscScalar,arg4::Vec,arg5::Vec)
    ccall((:MatZeroRowsIS,petsc),PetscErrorCode,(Mat,IS,PetscScalar,Vec,Vec),arg1,arg2,arg3,arg4,arg5)
end

function MatZeroRowsStencil(arg1::Mat,arg2::PetscInt,arg3::Ptr{MatStencil},arg4::PetscScalar,arg5::Vec,arg6::Vec)
    ccall((:MatZeroRowsStencil,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{MatStencil},PetscScalar,Vec,Vec),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatZeroRowsColumnsStencil(arg1::Mat,arg2::PetscInt,arg3::Ptr{MatStencil},arg4::PetscScalar,arg5::Vec,arg6::Vec)
    ccall((:MatZeroRowsColumnsStencil,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{MatStencil},PetscScalar,Vec,Vec),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatZeroRowsColumns(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::PetscScalar,arg5::Vec,arg6::Vec)
    ccall((:MatZeroRowsColumns,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscScalar,Vec,Vec),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatZeroRowsColumnsIS(arg1::Mat,arg2::IS,arg3::PetscScalar,arg4::Vec,arg5::Vec)
    ccall((:MatZeroRowsColumnsIS,petsc),PetscErrorCode,(Mat,IS,PetscScalar,Vec,Vec),arg1,arg2,arg3,arg4,arg5)
end

function MatGetSize(arg1::Mat,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt})
    ccall((:MatGetSize,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatGetLocalSize(arg1::Mat,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt})
    ccall((:MatGetLocalSize,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatGetOwnershipRange(arg1::Mat,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt})
    ccall((:MatGetOwnershipRange,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatGetOwnershipRanges(arg1::Mat,arg2::Ptr{Ptr{PetscInt}})
    ccall((:MatGetOwnershipRanges,petsc),PetscErrorCode,(Mat,Ptr{Ptr{PetscInt}}),arg1,arg2)
end

function MatGetOwnershipRangeColumn(arg1::Mat,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt})
    ccall((:MatGetOwnershipRangeColumn,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatGetOwnershipRangesColumn(arg1::Mat,arg2::Ptr{Ptr{PetscInt}})
    ccall((:MatGetOwnershipRangesColumn,petsc),PetscErrorCode,(Mat,Ptr{Ptr{PetscInt}}),arg1,arg2)
end

function MatGetOwnershipIS(arg1::Mat,arg2::Ptr{IS},arg3::Ptr{IS})
    ccall((:MatGetOwnershipIS,petsc),PetscErrorCode,(Mat,Ptr{IS},Ptr{IS}),arg1,arg2,arg3)
end

function MatGetSubMatrices(arg1::Mat,arg2::PetscInt,arg3::Ptr{IS},arg4::Ptr{IS},arg5::MatReuse,arg6::Ptr{Ptr{Mat}})
    ccall((:MatGetSubMatrices,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{IS},Ptr{IS},MatReuse,Ptr{Ptr{Mat}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatGetSubMatricesMPI(arg1::Mat,arg2::PetscInt,arg3::Ptr{IS},arg4::Ptr{IS},arg5::MatReuse,arg6::Ptr{Ptr{Mat}})
    ccall((:MatGetSubMatricesMPI,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{IS},Ptr{IS},MatReuse,Ptr{Ptr{Mat}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatDestroyMatrices(arg1::PetscInt,arg2::Ptr{Ptr{Mat}})
    ccall((:MatDestroyMatrices,petsc),PetscErrorCode,(PetscInt,Ptr{Ptr{Mat}}),arg1,arg2)
end

function MatGetSubMatrix(arg1::Mat,arg2::IS,arg3::IS,arg4::MatReuse,arg5::Ptr{Mat})
    ccall((:MatGetSubMatrix,petsc),PetscErrorCode,(Mat,IS,IS,MatReuse,Ptr{Mat}),arg1,arg2,arg3,arg4,arg5)
end

function MatGetLocalSubMatrix(arg1::Mat,arg2::IS,arg3::IS,arg4::Ptr{Mat})
    ccall((:MatGetLocalSubMatrix,petsc),PetscErrorCode,(Mat,IS,IS,Ptr{Mat}),arg1,arg2,arg3,arg4)
end

function MatRestoreLocalSubMatrix(arg1::Mat,arg2::IS,arg3::IS,arg4::Ptr{Mat})
    ccall((:MatRestoreLocalSubMatrix,petsc),PetscErrorCode,(Mat,IS,IS,Ptr{Mat}),arg1,arg2,arg3,arg4)
end

function MatGetSeqNonzeroStructure(arg1::Mat,arg2::Ptr{Mat})
    ccall((:MatGetSeqNonzeroStructure,petsc),PetscErrorCode,(Mat,Ptr{Mat}),arg1,arg2)
end

function MatDestroySeqNonzeroStructure(arg1::Ptr{Mat})
    ccall((:MatDestroySeqNonzeroStructure,petsc),PetscErrorCode,(Ptr{Mat},),arg1)
end

function MatCreateMPIAIJSumSeqAIJ(arg1::MPI_Comm,arg2::Mat,arg3::PetscInt,arg4::PetscInt,arg5::MatReuse,arg6::Ptr{Mat})
    ccall((:MatCreateMPIAIJSumSeqAIJ,petsc),PetscErrorCode,(MPI_Comm,Mat,PetscInt,PetscInt,MatReuse,Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatCreateMPIAIJSumSeqAIJSymbolic(arg1::MPI_Comm,arg2::Mat,arg3::PetscInt,arg4::PetscInt,arg5::Ptr{Mat})
    ccall((:MatCreateMPIAIJSumSeqAIJSymbolic,petsc),PetscErrorCode,(MPI_Comm,Mat,PetscInt,PetscInt,Ptr{Mat}),arg1,arg2,arg3,arg4,arg5)
end

function MatCreateMPIAIJSumSeqAIJNumeric(arg1::Mat,arg2::Mat)
    ccall((:MatCreateMPIAIJSumSeqAIJNumeric,petsc),PetscErrorCode,(Mat,Mat),arg1,arg2)
end

function MatMPIAIJGetLocalMat(arg1::Mat,arg2::MatReuse,arg3::Ptr{Mat})
    ccall((:MatMPIAIJGetLocalMat,petsc),PetscErrorCode,(Mat,MatReuse,Ptr{Mat}),arg1,arg2,arg3)
end

function MatMPIAIJGetLocalMatCondensed(arg1::Mat,arg2::MatReuse,arg3::Ptr{IS},arg4::Ptr{IS},arg5::Ptr{Mat})
    ccall((:MatMPIAIJGetLocalMatCondensed,petsc),PetscErrorCode,(Mat,MatReuse,Ptr{IS},Ptr{IS},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5)
end

function MatGetBrowsOfAcols(arg1::Mat,arg2::Mat,arg3::MatReuse,arg4::Ptr{IS},arg5::Ptr{IS},arg6::Ptr{Mat})
    ccall((:MatGetBrowsOfAcols,petsc),PetscErrorCode,(Mat,Mat,MatReuse,Ptr{IS},Ptr{IS},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatGetGhosts(arg1::Mat,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{PetscInt}})
    ccall((:MatGetGhosts,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3)
end

function MatIncreaseOverlap(arg1::Mat,arg2::PetscInt,arg3::Ptr{IS},arg4::PetscInt)
    ccall((:MatIncreaseOverlap,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{IS},PetscInt),arg1,arg2,arg3,arg4)
end

function MatMatMult(arg1::Mat,arg2::Mat,arg3::MatReuse,PetscReal::Cint,arg4::Ptr{Mat})
    ccall((:MatMatMult,petsc),PetscErrorCode,(Mat,Mat,MatReuse,Cint,Ptr{Mat}),arg1,arg2,arg3,PetscReal,arg4)
end

function MatMatMultSymbolic(arg1::Mat,arg2::Mat,PetscReal::Cint,arg3::Ptr{Mat})
    ccall((:MatMatMultSymbolic,petsc),PetscErrorCode,(Mat,Mat,Cint,Ptr{Mat}),arg1,arg2,PetscReal,arg3)
end

function MatMatMultNumeric(arg1::Mat,arg2::Mat,arg3::Mat)
    ccall((:MatMatMultNumeric,petsc),PetscErrorCode,(Mat,Mat,Mat),arg1,arg2,arg3)
end

function MatMatMatMult(arg1::Mat,arg2::Mat,arg3::Mat,arg4::MatReuse,PetscReal::Cint,arg5::Ptr{Mat})
    ccall((:MatMatMatMult,petsc),PetscErrorCode,(Mat,Mat,Mat,MatReuse,Cint,Ptr{Mat}),arg1,arg2,arg3,arg4,PetscReal,arg5)
end

function MatMatMatMultSymbolic(arg1::Mat,arg2::Mat,arg3::Mat,PetscReal::Cint,arg4::Ptr{Mat})
    ccall((:MatMatMatMultSymbolic,petsc),PetscErrorCode,(Mat,Mat,Mat,Cint,Ptr{Mat}),arg1,arg2,arg3,PetscReal,arg4)
end

function MatMatMatMultNumeric(arg1::Mat,arg2::Mat,arg3::Mat,arg4::Mat)
    ccall((:MatMatMatMultNumeric,petsc),PetscErrorCode,(Mat,Mat,Mat,Mat),arg1,arg2,arg3,arg4)
end

function MatPtAP(arg1::Mat,arg2::Mat,arg3::MatReuse,PetscReal::Cint,arg4::Ptr{Mat})
    ccall((:MatPtAP,petsc),PetscErrorCode,(Mat,Mat,MatReuse,Cint,Ptr{Mat}),arg1,arg2,arg3,PetscReal,arg4)
end

function MatPtAPSymbolic(arg1::Mat,arg2::Mat,PetscReal::Cint,arg3::Ptr{Mat})
    ccall((:MatPtAPSymbolic,petsc),PetscErrorCode,(Mat,Mat,Cint,Ptr{Mat}),arg1,arg2,PetscReal,arg3)
end

function MatPtAPNumeric(arg1::Mat,arg2::Mat,arg3::Mat)
    ccall((:MatPtAPNumeric,petsc),PetscErrorCode,(Mat,Mat,Mat),arg1,arg2,arg3)
end

function MatRARt(arg1::Mat,arg2::Mat,arg3::MatReuse,PetscReal::Cint,arg4::Ptr{Mat})
    ccall((:MatRARt,petsc),PetscErrorCode,(Mat,Mat,MatReuse,Cint,Ptr{Mat}),arg1,arg2,arg3,PetscReal,arg4)
end

function MatRARtSymbolic(arg1::Mat,arg2::Mat,PetscReal::Cint,arg3::Ptr{Mat})
    ccall((:MatRARtSymbolic,petsc),PetscErrorCode,(Mat,Mat,Cint,Ptr{Mat}),arg1,arg2,PetscReal,arg3)
end

function MatRARtNumeric(arg1::Mat,arg2::Mat,arg3::Mat)
    ccall((:MatRARtNumeric,petsc),PetscErrorCode,(Mat,Mat,Mat),arg1,arg2,arg3)
end

function MatTransposeMatMult(arg1::Mat,arg2::Mat,arg3::MatReuse,PetscReal::Cint,arg4::Ptr{Mat})
    ccall((:MatTransposeMatMult,petsc),PetscErrorCode,(Mat,Mat,MatReuse,Cint,Ptr{Mat}),arg1,arg2,arg3,PetscReal,arg4)
end

function MatTransposetMatMultSymbolic(arg1::Mat,arg2::Mat,PetscReal::Cint,arg3::Ptr{Mat})
    ccall((:MatTransposetMatMultSymbolic,petsc),PetscErrorCode,(Mat,Mat,Cint,Ptr{Mat}),arg1,arg2,PetscReal,arg3)
end

function MatTransposetMatMultNumeric(arg1::Mat,arg2::Mat,arg3::Mat)
    ccall((:MatTransposetMatMultNumeric,petsc),PetscErrorCode,(Mat,Mat,Mat),arg1,arg2,arg3)
end

function MatMatTransposeMult(arg1::Mat,arg2::Mat,arg3::MatReuse,PetscReal::Cint,arg4::Ptr{Mat})
    ccall((:MatMatTransposeMult,petsc),PetscErrorCode,(Mat,Mat,MatReuse,Cint,Ptr{Mat}),arg1,arg2,arg3,PetscReal,arg4)
end

function MatMatTransposeMultSymbolic(arg1::Mat,arg2::Mat,PetscReal::Cint,arg3::Ptr{Mat})
    ccall((:MatMatTransposeMultSymbolic,petsc),PetscErrorCode,(Mat,Mat,Cint,Ptr{Mat}),arg1,arg2,PetscReal,arg3)
end

function MatMatTransposeMultNumeric(arg1::Mat,arg2::Mat,arg3::Mat)
    ccall((:MatMatTransposeMultNumeric,petsc),PetscErrorCode,(Mat,Mat,Mat),arg1,arg2,arg3)
end

function MatAXPY(arg1::Mat,arg2::PetscScalar,arg3::Mat,arg4::MatStructure)
    ccall((:MatAXPY,petsc),PetscErrorCode,(Mat,PetscScalar,Mat,MatStructure),arg1,arg2,arg3,arg4)
end

function MatAYPX(arg1::Mat,arg2::PetscScalar,arg3::Mat,arg4::MatStructure)
    ccall((:MatAYPX,petsc),PetscErrorCode,(Mat,PetscScalar,Mat,MatStructure),arg1,arg2,arg3,arg4)
end

function MatScale(arg1::Mat,arg2::PetscScalar)
    ccall((:MatScale,petsc),PetscErrorCode,(Mat,PetscScalar),arg1,arg2)
end

function MatShift(arg1::Mat,arg2::PetscScalar)
    ccall((:MatShift,petsc),PetscErrorCode,(Mat,PetscScalar),arg1,arg2)
end

function MatSetLocalToGlobalMapping(arg1::Mat,arg2::ISLocalToGlobalMapping,arg3::ISLocalToGlobalMapping)
    ccall((:MatSetLocalToGlobalMapping,petsc),PetscErrorCode,(Mat,ISLocalToGlobalMapping,ISLocalToGlobalMapping),arg1,arg2,arg3)
end

function MatGetLocalToGlobalMapping(arg1::Mat,arg2::Ptr{ISLocalToGlobalMapping},arg3::Ptr{ISLocalToGlobalMapping})
    ccall((:MatGetLocalToGlobalMapping,petsc),PetscErrorCode,(Mat,Ptr{ISLocalToGlobalMapping},Ptr{ISLocalToGlobalMapping}),arg1,arg2,arg3)
end

function MatGetLayouts(arg1::Mat,arg2::Ptr{PetscLayout},arg3::Ptr{PetscLayout})
    ccall((:MatGetLayouts,petsc),PetscErrorCode,(Mat,Ptr{PetscLayout},Ptr{PetscLayout}),arg1,arg2,arg3)
end

function MatZeroRowsLocal(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::PetscScalar,arg5::Vec,arg6::Vec)
    ccall((:MatZeroRowsLocal,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscScalar,Vec,Vec),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatZeroRowsLocalIS(arg1::Mat,arg2::IS,arg3::PetscScalar,arg4::Vec,arg5::Vec)
    ccall((:MatZeroRowsLocalIS,petsc),PetscErrorCode,(Mat,IS,PetscScalar,Vec,Vec),arg1,arg2,arg3,arg4,arg5)
end

function MatZeroRowsColumnsLocal(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::PetscScalar,arg5::Vec,arg6::Vec)
    ccall((:MatZeroRowsColumnsLocal,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscScalar,Vec,Vec),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatZeroRowsColumnsLocalIS(arg1::Mat,arg2::IS,arg3::PetscScalar,arg4::Vec,arg5::Vec)
    ccall((:MatZeroRowsColumnsLocalIS,petsc),PetscErrorCode,(Mat,IS,PetscScalar,Vec,Vec),arg1,arg2,arg3,arg4,arg5)
end

function MatSetValuesLocal(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::PetscInt,arg5::Ptr{PetscInt},arg6::Ptr{PetscScalar},arg7::InsertMode)
    ccall((:MatSetValuesLocal,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{PetscScalar},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatSetValuesBlockedLocal(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::PetscInt,arg5::Ptr{PetscInt},arg6::Ptr{PetscScalar},arg7::InsertMode)
    ccall((:MatSetValuesBlockedLocal,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{PetscScalar},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatStashSetInitialSize(arg1::Mat,arg2::PetscInt,arg3::PetscInt)
    ccall((:MatStashSetInitialSize,petsc),PetscErrorCode,(Mat,PetscInt,PetscInt),arg1,arg2,arg3)
end

function MatStashGetInfo(arg1::Mat,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{PetscInt})
    ccall((:MatStashGetInfo,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end

function MatInterpolate(arg1::Mat,arg2::Vec,arg3::Vec)
    ccall((:MatInterpolate,petsc),PetscErrorCode,(Mat,Vec,Vec),arg1,arg2,arg3)
end

function MatInterpolateAdd(arg1::Mat,arg2::Vec,arg3::Vec,arg4::Vec)
    ccall((:MatInterpolateAdd,petsc),PetscErrorCode,(Mat,Vec,Vec,Vec),arg1,arg2,arg3,arg4)
end

function MatRestrict(arg1::Mat,arg2::Vec,arg3::Vec)
    ccall((:MatRestrict,petsc),PetscErrorCode,(Mat,Vec,Vec),arg1,arg2,arg3)
end

function MatCreateVecs(arg1::Mat,arg2::Ptr{Vec},arg3::Ptr{Vec})
    ccall((:MatCreateVecs,petsc),PetscErrorCode,(Mat,Ptr{Vec},Ptr{Vec}),arg1,arg2,arg3)
end

function MatGetMultiProcBlock(arg1::Mat,arg2::MPI_Comm,arg3::MatReuse,arg4::Ptr{Mat})
    ccall((:MatGetMultiProcBlock,petsc),PetscErrorCode,(Mat,MPI_Comm,MatReuse,Ptr{Mat}),arg1,arg2,arg3,arg4)
end

function MatFindZeroDiagonals(arg1::Mat,arg2::Ptr{IS})
    ccall((:MatFindZeroDiagonals,petsc),PetscErrorCode,(Mat,Ptr{IS}),arg1,arg2)
end

function MatFindOffBlockDiagonalEntries(arg1::Mat,arg2::Ptr{IS})
    ccall((:MatFindOffBlockDiagonalEntries,petsc),PetscErrorCode,(Mat,Ptr{IS}),arg1,arg2)
end

function MatCreateMPIMatConcatenateSeqMat(arg1::MPI_Comm,arg2::Mat,arg3::PetscInt,arg4::MatReuse,arg5::Ptr{Mat})
    ccall((:MatCreateMPIMatConcatenateSeqMat,petsc),PetscErrorCode,(MPI_Comm,Mat,PetscInt,MatReuse,Ptr{Mat}),arg1,arg2,arg3,arg4,arg5)
end

function MatInodeAdjustForInodes(arg1::Mat,arg2::Ptr{IS},arg3::Ptr{IS})
    ccall((:MatInodeAdjustForInodes,petsc),PetscErrorCode,(Mat,Ptr{IS},Ptr{IS}),arg1,arg2,arg3)
end

function MatInodeGetInodeSizes(arg1::Mat,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{PetscInt}},arg4::Ptr{PetscInt})
    ccall((:MatInodeGetInodeSizes,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function MatSeqAIJSetColumnIndices(arg1::Mat,arg2::Ptr{PetscInt})
    ccall((:MatSeqAIJSetColumnIndices,petsc),PetscErrorCode,(Mat,Ptr{PetscInt}),arg1,arg2)
end

function MatSeqBAIJSetColumnIndices(arg1::Mat,arg2::Ptr{PetscInt})
    ccall((:MatSeqBAIJSetColumnIndices,petsc),PetscErrorCode,(Mat,Ptr{PetscInt}),arg1,arg2)
end

function MatCreateSeqAIJWithArrays(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscInt},arg5::Ptr{PetscInt},arg6::Ptr{PetscScalar},arg7::Ptr{Mat})
    ccall((:MatCreateSeqAIJWithArrays,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscScalar},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateSeqBAIJWithArrays(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::Ptr{PetscInt},arg6::Ptr{PetscInt},arg7::Ptr{PetscScalar},arg8::Ptr{Mat})
    ccall((:MatCreateSeqBAIJWithArrays,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscScalar},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatCreateSeqSBAIJWithArrays(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::Ptr{PetscInt},arg6::Ptr{PetscInt},arg7::Ptr{PetscScalar},arg8::Ptr{Mat})
    ccall((:MatCreateSeqSBAIJWithArrays,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscScalar},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatCreateSeqAIJFromTriple(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscInt},arg5::Ptr{PetscInt},arg6::Ptr{PetscScalar},arg7::Ptr{Mat},arg8::PetscInt,arg9::PetscBool)
    ccall((:MatCreateSeqAIJFromTriple,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscScalar},Ptr{Mat},PetscInt,PetscBool),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function MatSeqBAIJSetPreallocation(arg1::Mat,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscInt})
    ccall((:MatSeqBAIJSetPreallocation,petsc),PetscErrorCode,(Mat,PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function MatSeqSBAIJSetPreallocation(arg1::Mat,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscInt})
    ccall((:MatSeqSBAIJSetPreallocation,petsc),PetscErrorCode,(Mat,PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function MatSeqAIJSetPreallocation(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt})
    ccall((:MatSeqAIJSetPreallocation,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatMPIBAIJSetPreallocation(arg1::Mat,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscInt},arg5::PetscInt,arg6::Ptr{PetscInt})
    ccall((:MatMPIBAIJSetPreallocation,petsc),PetscErrorCode,(Mat,PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatMPISBAIJSetPreallocation(arg1::Mat,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscInt},arg5::PetscInt,arg6::Ptr{PetscInt})
    ccall((:MatMPISBAIJSetPreallocation,petsc),PetscErrorCode,(Mat,PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatMPIAIJSetPreallocation(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::PetscInt,arg5::Ptr{PetscInt})
    ccall((:MatMPIAIJSetPreallocation,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end

function MatSeqAIJSetPreallocationCSR(arg1::Mat,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{PetscScalar})
    ccall((:MatSeqAIJSetPreallocationCSR,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscScalar}),arg1,arg2,arg3,arg4)
end

function MatSeqBAIJSetPreallocationCSR(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{PetscScalar})
    ccall((:MatSeqBAIJSetPreallocationCSR,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscScalar}),arg1,arg2,arg3,arg4,arg5)
end

function MatMPIAIJSetPreallocationCSR(arg1::Mat,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{PetscScalar})
    ccall((:MatMPIAIJSetPreallocationCSR,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscScalar}),arg1,arg2,arg3,arg4)
end

function MatMPIBAIJSetPreallocationCSR(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{PetscScalar})
    ccall((:MatMPIBAIJSetPreallocationCSR,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscScalar}),arg1,arg2,arg3,arg4,arg5)
end

function MatMPIAdjSetPreallocation(arg1::Mat,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:MatMPIAdjSetPreallocation,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function MatMPIDenseSetPreallocation(arg1::Mat,arg2::Ptr{PetscScalar})
    ccall((:MatMPIDenseSetPreallocation,petsc),PetscErrorCode,(Mat,Ptr{PetscScalar}),arg1,arg2)
end

function MatSeqDenseSetPreallocation(arg1::Mat,arg2::Ptr{PetscScalar})
    ccall((:MatSeqDenseSetPreallocation,petsc),PetscErrorCode,(Mat,Ptr{PetscScalar}),arg1,arg2)
end

function MatMPIAIJGetSeqAIJ(arg1::Mat,arg2::Ptr{Mat},arg3::Ptr{Mat},arg4::Ptr{Ptr{PetscInt}})
    ccall((:MatMPIAIJGetSeqAIJ,petsc),PetscErrorCode,(Mat,Ptr{Mat},Ptr{Mat},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4)
end

function MatMPIBAIJGetSeqBAIJ(arg1::Mat,arg2::Ptr{Mat},arg3::Ptr{Mat},arg4::Ptr{Ptr{PetscInt}})
    ccall((:MatMPIBAIJGetSeqBAIJ,petsc),PetscErrorCode,(Mat,Ptr{Mat},Ptr{Mat},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4)
end

function MatMPIAdjCreateNonemptySubcommMat(arg1::Mat,arg2::Ptr{Mat})
    ccall((:MatMPIAdjCreateNonemptySubcommMat,petsc),PetscErrorCode,(Mat,Ptr{Mat}),arg1,arg2)
end

function MatISSetPreallocation(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::PetscInt,arg5::Ptr{PetscInt})
    ccall((:MatISSetPreallocation,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end

function MatSeqDenseSetLDA(arg1::Mat,arg2::PetscInt)
    ccall((:MatSeqDenseSetLDA,petsc),PetscErrorCode,(Mat,PetscInt),arg1,arg2)
end

function MatDenseGetLocalMatrix(arg1::Mat,arg2::Ptr{Mat})
    ccall((:MatDenseGetLocalMatrix,petsc),PetscErrorCode,(Mat,Ptr{Mat}),arg1,arg2)
end

function MatStoreValues(arg1::Mat)
    ccall((:MatStoreValues,petsc),PetscErrorCode,(Mat,),arg1)
end

function MatRetrieveValues(arg1::Mat)
    ccall((:MatRetrieveValues,petsc),PetscErrorCode,(Mat,),arg1)
end

function MatDAADSetCtx(arg1::Mat,arg2::Ptr{Void})
    ccall((:MatDAADSetCtx,petsc),PetscErrorCode,(Mat,Ptr{Void}),arg1,arg2)
end

function MatFindNonzeroRows(arg1::Mat,arg2::Ptr{IS})
    ccall((:MatFindNonzeroRows,petsc),PetscErrorCode,(Mat,Ptr{IS}),arg1,arg2)
end

function MatGetOrdering(arg1::Mat,arg2::MatOrderingType,arg3::Ptr{IS},arg4::Ptr{IS})
    ccall((:MatGetOrdering,petsc),PetscErrorCode,(Mat,MatOrderingType,Ptr{IS},Ptr{IS}),arg1,arg2,arg3,arg4)
end

function MatGetOrderingList(arg1::Ptr{PetscFunctionList})
    ccall((:MatGetOrderingList,petsc),PetscErrorCode,(Ptr{PetscFunctionList},),arg1)
end

function MatOrderingRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:MatOrderingRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function MatReorderForNonzeroDiagonal(arg1::Mat,PetscReal::Cint,arg2::IS,arg3::IS)
    ccall((:MatReorderForNonzeroDiagonal,petsc),PetscErrorCode,(Mat,Cint,IS,IS),arg1,PetscReal,arg2,arg3)
end

function MatCreateLaplacian(arg1::Mat,PetscReal::Cint,arg2::PetscBool,arg3::Ptr{Mat})
    ccall((:MatCreateLaplacian,petsc),PetscErrorCode,(Mat,Cint,PetscBool,Ptr{Mat}),arg1,PetscReal,arg2,arg3)
end

function MatFactorInfoInitialize(arg1::Ptr{MatFactorInfo})
    ccall((:MatFactorInfoInitialize,petsc),PetscErrorCode,(Ptr{MatFactorInfo},),arg1)
end

function MatCholeskyFactor(arg1::Mat,arg2::IS,arg3::Ptr{MatFactorInfo})
    ccall((:MatCholeskyFactor,petsc),PetscErrorCode,(Mat,IS,Ptr{MatFactorInfo}),arg1,arg2,arg3)
end

function MatCholeskyFactorSymbolic(arg1::Mat,arg2::Mat,arg3::IS,arg4::Ptr{MatFactorInfo})
    ccall((:MatCholeskyFactorSymbolic,petsc),PetscErrorCode,(Mat,Mat,IS,Ptr{MatFactorInfo}),arg1,arg2,arg3,arg4)
end

function MatCholeskyFactorNumeric(arg1::Mat,arg2::Mat,arg3::Ptr{MatFactorInfo})
    ccall((:MatCholeskyFactorNumeric,petsc),PetscErrorCode,(Mat,Mat,Ptr{MatFactorInfo}),arg1,arg2,arg3)
end

function MatLUFactor(arg1::Mat,arg2::IS,arg3::IS,arg4::Ptr{MatFactorInfo})
    ccall((:MatLUFactor,petsc),PetscErrorCode,(Mat,IS,IS,Ptr{MatFactorInfo}),arg1,arg2,arg3,arg4)
end

function MatILUFactor(arg1::Mat,arg2::IS,arg3::IS,arg4::Ptr{MatFactorInfo})
    ccall((:MatILUFactor,petsc),PetscErrorCode,(Mat,IS,IS,Ptr{MatFactorInfo}),arg1,arg2,arg3,arg4)
end

function MatLUFactorSymbolic(arg1::Mat,arg2::Mat,arg3::IS,arg4::IS,arg5::Ptr{MatFactorInfo})
    ccall((:MatLUFactorSymbolic,petsc),PetscErrorCode,(Mat,Mat,IS,IS,Ptr{MatFactorInfo}),arg1,arg2,arg3,arg4,arg5)
end

function MatILUFactorSymbolic(arg1::Mat,arg2::Mat,arg3::IS,arg4::IS,arg5::Ptr{MatFactorInfo})
    ccall((:MatILUFactorSymbolic,petsc),PetscErrorCode,(Mat,Mat,IS,IS,Ptr{MatFactorInfo}),arg1,arg2,arg3,arg4,arg5)
end

function MatICCFactorSymbolic(arg1::Mat,arg2::Mat,arg3::IS,arg4::Ptr{MatFactorInfo})
    ccall((:MatICCFactorSymbolic,petsc),PetscErrorCode,(Mat,Mat,IS,Ptr{MatFactorInfo}),arg1,arg2,arg3,arg4)
end

function MatICCFactor(arg1::Mat,arg2::IS,arg3::Ptr{MatFactorInfo})
    ccall((:MatICCFactor,petsc),PetscErrorCode,(Mat,IS,Ptr{MatFactorInfo}),arg1,arg2,arg3)
end

function MatLUFactorNumeric(arg1::Mat,arg2::Mat,arg3::Ptr{MatFactorInfo})
    ccall((:MatLUFactorNumeric,petsc),PetscErrorCode,(Mat,Mat,Ptr{MatFactorInfo}),arg1,arg2,arg3)
end

function MatGetInertia(arg1::Mat,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:MatGetInertia,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function MatSolve(arg1::Mat,arg2::Vec,arg3::Vec)
    ccall((:MatSolve,petsc),PetscErrorCode,(Mat,Vec,Vec),arg1,arg2,arg3)
end

function MatForwardSolve(arg1::Mat,arg2::Vec,arg3::Vec)
    ccall((:MatForwardSolve,petsc),PetscErrorCode,(Mat,Vec,Vec),arg1,arg2,arg3)
end

function MatBackwardSolve(arg1::Mat,arg2::Vec,arg3::Vec)
    ccall((:MatBackwardSolve,petsc),PetscErrorCode,(Mat,Vec,Vec),arg1,arg2,arg3)
end

function MatSolveAdd(arg1::Mat,arg2::Vec,arg3::Vec,arg4::Vec)
    ccall((:MatSolveAdd,petsc),PetscErrorCode,(Mat,Vec,Vec,Vec),arg1,arg2,arg3,arg4)
end

function MatSolveTranspose(arg1::Mat,arg2::Vec,arg3::Vec)
    ccall((:MatSolveTranspose,petsc),PetscErrorCode,(Mat,Vec,Vec),arg1,arg2,arg3)
end

function MatSolveTransposeAdd(arg1::Mat,arg2::Vec,arg3::Vec,arg4::Vec)
    ccall((:MatSolveTransposeAdd,petsc),PetscErrorCode,(Mat,Vec,Vec,Vec),arg1,arg2,arg3,arg4)
end

function MatSolves(arg1::Mat,arg2::Vecs,arg3::Vecs)
    ccall((:MatSolves,petsc),PetscErrorCode,(Mat,Vecs,Vecs),arg1,arg2,arg3)
end

function MatSetUnfactored(arg1::Mat)
    ccall((:MatSetUnfactored,petsc),PetscErrorCode,(Mat,),arg1)
end

function MatSOR(arg1::Mat,arg2::Vec,PetscReal::Cint,arg3::MatSORType,arg4::Cint,arg5::PetscInt,arg6::PetscInt,arg7::Vec)
    ccall((:MatSOR,petsc),PetscErrorCode,(Mat,Vec,Cint,MatSORType,Cint,PetscInt,PetscInt,Vec),arg1,arg2,PetscReal,arg3,arg4,arg5,arg6,arg7)
end

function MatColoringCreate(arg1::Mat,arg2::Ptr{MatColoring})
    ccall((:MatColoringCreate,petsc),PetscErrorCode,(Mat,Ptr{MatColoring}),arg1,arg2)
end

function MatColoringGetDegrees(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt})
    ccall((:MatColoringGetDegrees,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatColoringDestroy(arg1::Ptr{MatColoring})
    ccall((:MatColoringDestroy,petsc),PetscErrorCode,(Ptr{MatColoring},),arg1)
end

function MatColoringView(arg1::MatColoring,arg2::PetscViewer)
    ccall((:MatColoringView,petsc),PetscErrorCode,(MatColoring,PetscViewer),arg1,arg2)
end

function MatColoringSetType(arg1::MatColoring,arg2::MatColoringType)
    ccall((:MatColoringSetType,petsc),PetscErrorCode,(MatColoring,MatColoringType),arg1,arg2)
end

function MatColoringSetFromOptions(arg1::MatColoring)
    ccall((:MatColoringSetFromOptions,petsc),PetscErrorCode,(MatColoring,),arg1)
end

function MatColoringSetDistance(arg1::MatColoring,arg2::PetscInt)
    ccall((:MatColoringSetDistance,petsc),PetscErrorCode,(MatColoring,PetscInt),arg1,arg2)
end

function MatColoringGetDistance(arg1::MatColoring,arg2::Ptr{PetscInt})
    ccall((:MatColoringGetDistance,petsc),PetscErrorCode,(MatColoring,Ptr{PetscInt}),arg1,arg2)
end

function MatColoringSetMaxColors(arg1::MatColoring,arg2::PetscInt)
    ccall((:MatColoringSetMaxColors,petsc),PetscErrorCode,(MatColoring,PetscInt),arg1,arg2)
end

function MatColoringGetMaxColors(arg1::MatColoring,arg2::Ptr{PetscInt})
    ccall((:MatColoringGetMaxColors,petsc),PetscErrorCode,(MatColoring,Ptr{PetscInt}),arg1,arg2)
end

function MatColoringApply(arg1::MatColoring,arg2::Ptr{ISColoring})
    ccall((:MatColoringApply,petsc),PetscErrorCode,(MatColoring,Ptr{ISColoring}),arg1,arg2)
end

function MatColoringRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:MatColoringRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function MatColoringPatch(arg1::Mat,arg2::PetscInt,arg3::PetscInt,ISColoringValue::Ptr{Cint},arg4::Ptr{ISColoring})
    ccall((:MatColoringPatch,petsc),PetscErrorCode,(Mat,PetscInt,PetscInt,Ptr{Cint},Ptr{ISColoring}),arg1,arg2,arg3,ISColoringValue,arg4)
end

function MatColoringSetWeightType(arg1::MatColoring,arg2::MatColoringWeightType)
    ccall((:MatColoringSetWeightType,petsc),PetscErrorCode,(MatColoring,MatColoringWeightType),arg1,arg2)
end

function MatColoringSetWeights(arg1::MatColoring,arg2::Ptr{Cint},arg3::Ptr{PetscInt})
    ccall((:MatColoringSetWeights,petsc),PetscErrorCode,(MatColoring,Ptr{Cint},Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatColoringCreateWeights(arg1::MatColoring,arg2::Ptr{Ptr{Cint}},lperm::Ptr{Ptr{PetscInt}})
    ccall((:MatColoringCreateWeights,petsc),PetscErrorCode,(MatColoring,Ptr{Ptr{Cint}},Ptr{Ptr{PetscInt}}),arg1,arg2,lperm)
end

function MatFDColoringCreate(arg1::Mat,arg2::ISColoring,arg3::Ptr{MatFDColoring})
    ccall((:MatFDColoringCreate,petsc),PetscErrorCode,(Mat,ISColoring,Ptr{MatFDColoring}),arg1,arg2,arg3)
end

function MatFDColoringDestroy(arg1::Ptr{MatFDColoring})
    ccall((:MatFDColoringDestroy,petsc),PetscErrorCode,(Ptr{MatFDColoring},),arg1)
end

function MatFDColoringView(arg1::MatFDColoring,arg2::PetscViewer)
    ccall((:MatFDColoringView,petsc),PetscErrorCode,(MatFDColoring,PetscViewer),arg1,arg2)
end

function MatFDColoringSetFunction(arg1::MatFDColoring,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:MatFDColoringSetFunction,petsc),PetscErrorCode,(MatFDColoring,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function MatFDColoringGetFunction(arg1::MatFDColoring,arg2::Ptr{Ptr{Void}},arg3::Ptr{Ptr{Void}})
    ccall((:MatFDColoringGetFunction,petsc),PetscErrorCode,(MatFDColoring,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end

function MatFDColoringSetParameters(arg1::MatFDColoring,PetscReal::Cint,arg2::Cint)
    ccall((:MatFDColoringSetParameters,petsc),PetscErrorCode,(MatFDColoring,Cint,Cint),arg1,PetscReal,arg2)
end

function MatFDColoringSetFromOptions(arg1::MatFDColoring)
    ccall((:MatFDColoringSetFromOptions,petsc),PetscErrorCode,(MatFDColoring,),arg1)
end

function MatFDColoringApply(arg1::Mat,arg2::MatFDColoring,arg3::Vec,arg4::Ptr{Void})
    ccall((:MatFDColoringApply,petsc),PetscErrorCode,(Mat,MatFDColoring,Vec,Ptr{Void}),arg1,arg2,arg3,arg4)
end

function MatFDColoringSetF(arg1::MatFDColoring,arg2::Vec)
    ccall((:MatFDColoringSetF,petsc),PetscErrorCode,(MatFDColoring,Vec),arg1,arg2)
end

function MatFDColoringGetPerturbedColumns(arg1::MatFDColoring,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{PetscInt}})
    ccall((:MatFDColoringGetPerturbedColumns,petsc),PetscErrorCode,(MatFDColoring,Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3)
end

function MatFDColoringSetUp(arg1::Mat,arg2::ISColoring,arg3::MatFDColoring)
    ccall((:MatFDColoringSetUp,petsc),PetscErrorCode,(Mat,ISColoring,MatFDColoring),arg1,arg2,arg3)
end

function MatFDColoringSetBlockSize(arg1::MatFDColoring,arg2::PetscInt,arg3::PetscInt)
    ccall((:MatFDColoringSetBlockSize,petsc),PetscErrorCode,(MatFDColoring,PetscInt,PetscInt),arg1,arg2,arg3)
end

function MatTransposeColoringCreate(arg1::Mat,arg2::ISColoring,arg3::Ptr{MatTransposeColoring})
    ccall((:MatTransposeColoringCreate,petsc),PetscErrorCode,(Mat,ISColoring,Ptr{MatTransposeColoring}),arg1,arg2,arg3)
end

function MatTransColoringApplySpToDen(arg1::MatTransposeColoring,arg2::Mat,arg3::Mat)
    ccall((:MatTransColoringApplySpToDen,petsc),PetscErrorCode,(MatTransposeColoring,Mat,Mat),arg1,arg2,arg3)
end

function MatTransColoringApplyDenToSp(arg1::MatTransposeColoring,arg2::Mat,arg3::Mat)
    ccall((:MatTransColoringApplyDenToSp,petsc),PetscErrorCode,(MatTransposeColoring,Mat,Mat),arg1,arg2,arg3)
end

function MatTransposeColoringDestroy(arg1::Ptr{MatTransposeColoring})
    ccall((:MatTransposeColoringDestroy,petsc),PetscErrorCode,(Ptr{MatTransposeColoring},),arg1)
end

function MatPartitioningCreate(arg1::MPI_Comm,arg2::Ptr{MatPartitioning})
    ccall((:MatPartitioningCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{MatPartitioning}),arg1,arg2)
end

function MatPartitioningSetType(arg1::MatPartitioning,arg2::MatPartitioningType)
    ccall((:MatPartitioningSetType,petsc),PetscErrorCode,(MatPartitioning,MatPartitioningType),arg1,arg2)
end

function MatPartitioningSetNParts(arg1::MatPartitioning,arg2::PetscInt)
    ccall((:MatPartitioningSetNParts,petsc),PetscErrorCode,(MatPartitioning,PetscInt),arg1,arg2)
end

function MatPartitioningSetAdjacency(arg1::MatPartitioning,arg2::Mat)
    ccall((:MatPartitioningSetAdjacency,petsc),PetscErrorCode,(MatPartitioning,Mat),arg1,arg2)
end

function MatPartitioningSetVertexWeights(arg1::MatPartitioning,arg2::Ptr{PetscInt})
    ccall((:MatPartitioningSetVertexWeights,petsc),PetscErrorCode,(MatPartitioning,Ptr{PetscInt}),arg1,arg2)
end

function MatPartitioningSetPartitionWeights(arg1::MatPartitioning,PetscReal::Ptr{Cint})
    ccall((:MatPartitioningSetPartitionWeights,petsc),PetscErrorCode,(MatPartitioning,Ptr{Cint}),arg1,PetscReal)
end

function MatPartitioningApply(arg1::MatPartitioning,arg2::Ptr{IS})
    ccall((:MatPartitioningApply,petsc),PetscErrorCode,(MatPartitioning,Ptr{IS}),arg1,arg2)
end

function MatPartitioningDestroy(arg1::Ptr{MatPartitioning})
    ccall((:MatPartitioningDestroy,petsc),PetscErrorCode,(Ptr{MatPartitioning},),arg1)
end

function MatPartitioningRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:MatPartitioningRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function MatPartitioningView(arg1::MatPartitioning,arg2::PetscViewer)
    ccall((:MatPartitioningView,petsc),PetscErrorCode,(MatPartitioning,PetscViewer),arg1,arg2)
end

function MatPartitioningSetFromOptions(arg1::MatPartitioning)
    ccall((:MatPartitioningSetFromOptions,petsc),PetscErrorCode,(MatPartitioning,),arg1)
end

function MatPartitioningGetType(arg1::MatPartitioning,arg2::Ptr{MatPartitioningType})
    ccall((:MatPartitioningGetType,petsc),PetscErrorCode,(MatPartitioning,Ptr{MatPartitioningType}),arg1,arg2)
end

function MatPartitioningParmetisSetCoarseSequential(arg1::MatPartitioning)
    ccall((:MatPartitioningParmetisSetCoarseSequential,petsc),PetscErrorCode,(MatPartitioning,),arg1)
end

function MatPartitioningParmetisGetEdgeCut(arg1::MatPartitioning,arg2::Ptr{PetscInt})
    ccall((:MatPartitioningParmetisGetEdgeCut,petsc),PetscErrorCode,(MatPartitioning,Ptr{PetscInt}),arg1,arg2)
end

function MatPartitioningChacoSetGlobal(arg1::MatPartitioning,arg2::MPChacoGlobalType)
    ccall((:MatPartitioningChacoSetGlobal,petsc),PetscErrorCode,(MatPartitioning,MPChacoGlobalType),arg1,arg2)
end

function MatPartitioningChacoGetGlobal(arg1::MatPartitioning,arg2::Ptr{MPChacoGlobalType})
    ccall((:MatPartitioningChacoGetGlobal,petsc),PetscErrorCode,(MatPartitioning,Ptr{MPChacoGlobalType}),arg1,arg2)
end

function MatPartitioningChacoSetLocal(arg1::MatPartitioning,arg2::MPChacoLocalType)
    ccall((:MatPartitioningChacoSetLocal,petsc),PetscErrorCode,(MatPartitioning,MPChacoLocalType),arg1,arg2)
end

function MatPartitioningChacoGetLocal(arg1::MatPartitioning,arg2::Ptr{MPChacoLocalType})
    ccall((:MatPartitioningChacoGetLocal,petsc),PetscErrorCode,(MatPartitioning,Ptr{MPChacoLocalType}),arg1,arg2)
end

function MatPartitioningChacoSetCoarseLevel(arg1::MatPartitioning,PetscReal::Cint)
    ccall((:MatPartitioningChacoSetCoarseLevel,petsc),PetscErrorCode,(MatPartitioning,Cint),arg1,PetscReal)
end

function MatPartitioningChacoSetEigenSolver(arg1::MatPartitioning,arg2::MPChacoEigenType)
    ccall((:MatPartitioningChacoSetEigenSolver,petsc),PetscErrorCode,(MatPartitioning,MPChacoEigenType),arg1,arg2)
end

function MatPartitioningChacoGetEigenSolver(arg1::MatPartitioning,arg2::Ptr{MPChacoEigenType})
    ccall((:MatPartitioningChacoGetEigenSolver,petsc),PetscErrorCode,(MatPartitioning,Ptr{MPChacoEigenType}),arg1,arg2)
end

function MatPartitioningChacoSetEigenTol(arg1::MatPartitioning,PetscReal::Cint)
    ccall((:MatPartitioningChacoSetEigenTol,petsc),PetscErrorCode,(MatPartitioning,Cint),arg1,PetscReal)
end

function MatPartitioningChacoGetEigenTol(arg1::MatPartitioning,arg2::Ptr{Cint})
    ccall((:MatPartitioningChacoGetEigenTol,petsc),PetscErrorCode,(MatPartitioning,Ptr{Cint}),arg1,arg2)
end

function MatPartitioningChacoSetEigenNumber(arg1::MatPartitioning,arg2::PetscInt)
    ccall((:MatPartitioningChacoSetEigenNumber,petsc),PetscErrorCode,(MatPartitioning,PetscInt),arg1,arg2)
end

function MatPartitioningChacoGetEigenNumber(arg1::MatPartitioning,arg2::Ptr{PetscInt})
    ccall((:MatPartitioningChacoGetEigenNumber,petsc),PetscErrorCode,(MatPartitioning,Ptr{PetscInt}),arg1,arg2)
end

function MatPartitioningPartySetGlobal(arg1::MatPartitioning,arg2::Ptr{Uint8})
    ccall((:MatPartitioningPartySetGlobal,petsc),PetscErrorCode,(MatPartitioning,Ptr{Uint8}),arg1,arg2)
end

function MatPartitioningPartySetLocal(arg1::MatPartitioning,arg2::Ptr{Uint8})
    ccall((:MatPartitioningPartySetLocal,petsc),PetscErrorCode,(MatPartitioning,Ptr{Uint8}),arg1,arg2)
end

function MatPartitioningPartySetCoarseLevel(arg1::MatPartitioning,PetscReal::Cint)
    ccall((:MatPartitioningPartySetCoarseLevel,petsc),PetscErrorCode,(MatPartitioning,Cint),arg1,PetscReal)
end

function MatPartitioningPartySetBipart(arg1::MatPartitioning,arg2::PetscBool)
    ccall((:MatPartitioningPartySetBipart,petsc),PetscErrorCode,(MatPartitioning,PetscBool),arg1,arg2)
end

function MatPartitioningPartySetMatchOptimization(arg1::MatPartitioning,arg2::PetscBool)
    ccall((:MatPartitioningPartySetMatchOptimization,petsc),PetscErrorCode,(MatPartitioning,PetscBool),arg1,arg2)
end

function MatPartitioningPTScotchSetImbalance(arg1::MatPartitioning,PetscReal::Cint)
    ccall((:MatPartitioningPTScotchSetImbalance,petsc),PetscErrorCode,(MatPartitioning,Cint),arg1,PetscReal)
end

function MatPartitioningPTScotchGetImbalance(arg1::MatPartitioning,arg2::Ptr{Cint})
    ccall((:MatPartitioningPTScotchGetImbalance,petsc),PetscErrorCode,(MatPartitioning,Ptr{Cint}),arg1,arg2)
end

function MatPartitioningPTScotchSetStrategy(arg1::MatPartitioning,arg2::MPPTScotchStrategyType)
    ccall((:MatPartitioningPTScotchSetStrategy,petsc),PetscErrorCode,(MatPartitioning,MPPTScotchStrategyType),arg1,arg2)
end

function MatPartitioningPTScotchGetStrategy(arg1::MatPartitioning,arg2::Ptr{MPPTScotchStrategyType})
    ccall((:MatPartitioningPTScotchGetStrategy,petsc),PetscErrorCode,(MatPartitioning,Ptr{MPPTScotchStrategyType}),arg1,arg2)
end

function MatCoarsenCreate(arg1::MPI_Comm,arg2::Ptr{MatCoarsen})
    ccall((:MatCoarsenCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{MatCoarsen}),arg1,arg2)
end

function MatCoarsenSetType(arg1::MatCoarsen,arg2::MatCoarsenType)
    ccall((:MatCoarsenSetType,petsc),PetscErrorCode,(MatCoarsen,MatCoarsenType),arg1,arg2)
end

function MatCoarsenSetAdjacency(arg1::MatCoarsen,arg2::Mat)
    ccall((:MatCoarsenSetAdjacency,petsc),PetscErrorCode,(MatCoarsen,Mat),arg1,arg2)
end

function MatCoarsenSetGreedyOrdering(arg1::MatCoarsen,arg2::IS)
    ccall((:MatCoarsenSetGreedyOrdering,petsc),PetscErrorCode,(MatCoarsen,IS),arg1,arg2)
end

function MatCoarsenSetStrictAggs(arg1::MatCoarsen,arg2::PetscBool)
    ccall((:MatCoarsenSetStrictAggs,petsc),PetscErrorCode,(MatCoarsen,PetscBool),arg1,arg2)
end

function MatCoarsenGetData(arg1::MatCoarsen,arg2::Ptr{Ptr{PetscCoarsenData}})
    ccall((:MatCoarsenGetData,petsc),PetscErrorCode,(MatCoarsen,Ptr{Ptr{PetscCoarsenData}}),arg1,arg2)
end

function MatCoarsenApply(arg1::MatCoarsen)
    ccall((:MatCoarsenApply,petsc),PetscErrorCode,(MatCoarsen,),arg1)
end

function MatCoarsenDestroy(arg1::Ptr{MatCoarsen})
    ccall((:MatCoarsenDestroy,petsc),PetscErrorCode,(Ptr{MatCoarsen},),arg1)
end

function MatCoarsenRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:MatCoarsenRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function MatCoarsenView(arg1::MatCoarsen,arg2::PetscViewer)
    ccall((:MatCoarsenView,petsc),PetscErrorCode,(MatCoarsen,PetscViewer),arg1,arg2)
end

function MatCoarsenSetFromOptions(arg1::MatCoarsen)
    ccall((:MatCoarsenSetFromOptions,petsc),PetscErrorCode,(MatCoarsen,),arg1)
end

function MatCoarsenGetType(arg1::MatCoarsen,arg2::Ptr{MatCoarsenType})
    ccall((:MatCoarsenGetType,petsc),PetscErrorCode,(MatCoarsen,Ptr{MatCoarsenType}),arg1,arg2)
end

function MatMeshToCellGraph(arg1::Mat,arg2::PetscInt,arg3::Ptr{Mat})
    ccall((:MatMeshToCellGraph,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{Mat}),arg1,arg2,arg3)
end

function MatHasOperation(arg1::Mat,arg2::MatOperation,arg3::Ptr{PetscBool})
    ccall((:MatHasOperation,petsc),PetscErrorCode,(Mat,MatOperation,Ptr{PetscBool}),arg1,arg2,arg3)
end

function MatShellSetOperation(arg1::Mat,arg2::MatOperation,arg3::Ptr{Void})
    ccall((:MatShellSetOperation,petsc),PetscErrorCode,(Mat,MatOperation,Ptr{Void}),arg1,arg2,arg3)
end

function MatShellGetOperation(arg1::Mat,arg2::MatOperation,arg3::Ptr{Ptr{Void}})
    ccall((:MatShellGetOperation,petsc),PetscErrorCode,(Mat,MatOperation,Ptr{Ptr{Void}}),arg1,arg2,arg3)
end

function MatShellSetContext(arg1::Mat,arg2::Ptr{Void})
    ccall((:MatShellSetContext,petsc),PetscErrorCode,(Mat,Ptr{Void}),arg1,arg2)
end

function MatMPIBAIJSetHashTableFactor(arg1::Mat,PetscReal::Cint)
    ccall((:MatMPIBAIJSetHashTableFactor,petsc),PetscErrorCode,(Mat,Cint),arg1,PetscReal)
end

function MatISGetLocalMat(arg1::Mat,arg2::Ptr{Mat})
    ccall((:MatISGetLocalMat,petsc),PetscErrorCode,(Mat,Ptr{Mat}),arg1,arg2)
end

function MatISSetLocalMat(arg1::Mat,arg2::Mat)
    ccall((:MatISSetLocalMat,petsc),PetscErrorCode,(Mat,Mat),arg1,arg2)
end

function MatISGetMPIXAIJ(arg1::Mat,arg2::MatReuse,arg3::Ptr{Mat})
    ccall((:MatISGetMPIXAIJ,petsc),PetscErrorCode,(Mat,MatReuse,Ptr{Mat}),arg1,arg2,arg3)
end

function MatNullSpaceCreate(arg1::MPI_Comm,arg2::PetscBool,arg3::PetscInt,arg4::Ptr{Vec},arg5::Ptr{MatNullSpace})
    ccall((:MatNullSpaceCreate,petsc),PetscErrorCode,(MPI_Comm,PetscBool,PetscInt,Ptr{Vec},Ptr{MatNullSpace}),arg1,arg2,arg3,arg4,arg5)
end

function MatNullSpaceSetFunction(arg1::MatNullSpace,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:MatNullSpaceSetFunction,petsc),PetscErrorCode,(MatNullSpace,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function MatNullSpaceDestroy(arg1::Ptr{MatNullSpace})
    ccall((:MatNullSpaceDestroy,petsc),PetscErrorCode,(Ptr{MatNullSpace},),arg1)
end

function MatNullSpaceRemove(arg1::MatNullSpace,arg2::Vec)
    ccall((:MatNullSpaceRemove,petsc),PetscErrorCode,(MatNullSpace,Vec),arg1,arg2)
end

function MatGetNullSpace(arg1::Mat,arg2::Ptr{MatNullSpace})
    ccall((:MatGetNullSpace,petsc),PetscErrorCode,(Mat,Ptr{MatNullSpace}),arg1,arg2)
end

function MatGetTransposeNullSpace(arg1::Mat,arg2::Ptr{MatNullSpace})
    ccall((:MatGetTransposeNullSpace,petsc),PetscErrorCode,(Mat,Ptr{MatNullSpace}),arg1,arg2)
end

function MatSetTransposeNullSpace(arg1::Mat,arg2::MatNullSpace)
    ccall((:MatSetTransposeNullSpace,petsc),PetscErrorCode,(Mat,MatNullSpace),arg1,arg2)
end

function MatSetNullSpace(arg1::Mat,arg2::MatNullSpace)
    ccall((:MatSetNullSpace,petsc),PetscErrorCode,(Mat,MatNullSpace),arg1,arg2)
end

function MatSetNearNullSpace(arg1::Mat,arg2::MatNullSpace)
    ccall((:MatSetNearNullSpace,petsc),PetscErrorCode,(Mat,MatNullSpace),arg1,arg2)
end

function MatGetNearNullSpace(arg1::Mat,arg2::Ptr{MatNullSpace})
    ccall((:MatGetNearNullSpace,petsc),PetscErrorCode,(Mat,Ptr{MatNullSpace}),arg1,arg2)
end

function MatNullSpaceTest(arg1::MatNullSpace,arg2::Mat,arg3::Ptr{PetscBool})
    ccall((:MatNullSpaceTest,petsc),PetscErrorCode,(MatNullSpace,Mat,Ptr{PetscBool}),arg1,arg2,arg3)
end

function MatNullSpaceView(arg1::MatNullSpace,arg2::PetscViewer)
    ccall((:MatNullSpaceView,petsc),PetscErrorCode,(MatNullSpace,PetscViewer),arg1,arg2)
end

function MatNullSpaceGetVecs(arg1::MatNullSpace,arg2::Ptr{PetscBool},arg3::Ptr{PetscInt},arg4::Ptr{Ptr{Vec}})
    ccall((:MatNullSpaceGetVecs,petsc),PetscErrorCode,(MatNullSpace,Ptr{PetscBool},Ptr{PetscInt},Ptr{Ptr{Vec}}),arg1,arg2,arg3,arg4)
end

function MatNullSpaceCreateRigidBody(arg1::Vec,arg2::Ptr{MatNullSpace})
    ccall((:MatNullSpaceCreateRigidBody,petsc),PetscErrorCode,(Vec,Ptr{MatNullSpace}),arg1,arg2)
end

function MatReorderingSeqSBAIJ(arg1::Mat,arg2::IS)
    ccall((:MatReorderingSeqSBAIJ,petsc),PetscErrorCode,(Mat,IS),arg1,arg2)
end

function MatMPISBAIJSetHashTableFactor(arg1::Mat,PetscReal::Cint)
    ccall((:MatMPISBAIJSetHashTableFactor,petsc),PetscErrorCode,(Mat,Cint),arg1,PetscReal)
end

function MatSeqSBAIJSetColumnIndices(arg1::Mat,arg2::Ptr{PetscInt})
    ccall((:MatSeqSBAIJSetColumnIndices,petsc),PetscErrorCode,(Mat,Ptr{PetscInt}),arg1,arg2)
end

function MatSeqBAIJInvertBlockDiagonal(arg1::Mat)
    ccall((:MatSeqBAIJInvertBlockDiagonal,petsc),PetscErrorCode,(Mat,),arg1)
end

function MatCreateMAIJ(arg1::Mat,arg2::PetscInt,arg3::Ptr{Mat})
    ccall((:MatCreateMAIJ,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{Mat}),arg1,arg2,arg3)
end

function MatMAIJRedimension(arg1::Mat,arg2::PetscInt,arg3::Ptr{Mat})
    ccall((:MatMAIJRedimension,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{Mat}),arg1,arg2,arg3)
end

function MatMAIJGetAIJ(arg1::Mat,arg2::Ptr{Mat})
    ccall((:MatMAIJGetAIJ,petsc),PetscErrorCode,(Mat,Ptr{Mat}),arg1,arg2)
end

function MatComputeExplicitOperator(arg1::Mat,arg2::Ptr{Mat})
    ccall((:MatComputeExplicitOperator,petsc),PetscErrorCode,(Mat,Ptr{Mat}),arg1,arg2)
end

function MatDiagonalScaleLocal(arg1::Mat,arg2::Vec)
    ccall((:MatDiagonalScaleLocal,petsc),PetscErrorCode,(Mat,Vec),arg1,arg2)
end

function MatCreateMFFD(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::Ptr{Mat})
    ccall((:MatCreateMFFD,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatMFFDSetBase(arg1::Mat,arg2::Vec,arg3::Vec)
    ccall((:MatMFFDSetBase,petsc),PetscErrorCode,(Mat,Vec,Vec),arg1,arg2,arg3)
end

function MatMFFDSetFunction(arg1::Mat,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:MatMFFDSetFunction,petsc),PetscErrorCode,(Mat,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function MatMFFDSetFunctioni(arg1::Mat,arg2::Ptr{Void})
    ccall((:MatMFFDSetFunctioni,petsc),PetscErrorCode,(Mat,Ptr{Void}),arg1,arg2)
end

function MatMFFDSetFunctioniBase(arg1::Mat,arg2::Ptr{Void})
    ccall((:MatMFFDSetFunctioniBase,petsc),PetscErrorCode,(Mat,Ptr{Void}),arg1,arg2)
end

function MatMFFDAddNullSpace(arg1::Mat,arg2::MatNullSpace)
    ccall((:MatMFFDAddNullSpace,petsc),PetscErrorCode,(Mat,MatNullSpace),arg1,arg2)
end

function MatMFFDSetHHistory(arg1::Mat,arg2::Ptr{PetscScalar},arg3::PetscInt)
    ccall((:MatMFFDSetHHistory,petsc),PetscErrorCode,(Mat,Ptr{PetscScalar},PetscInt),arg1,arg2,arg3)
end

function MatMFFDResetHHistory(arg1::Mat)
    ccall((:MatMFFDResetHHistory,petsc),PetscErrorCode,(Mat,),arg1)
end

function MatMFFDSetFunctionError(arg1::Mat,PetscReal::Cint)
    ccall((:MatMFFDSetFunctionError,petsc),PetscErrorCode,(Mat,Cint),arg1,PetscReal)
end

function MatMFFDSetPeriod(arg1::Mat,arg2::PetscInt)
    ccall((:MatMFFDSetPeriod,petsc),PetscErrorCode,(Mat,PetscInt),arg1,arg2)
end

function MatMFFDGetH(arg1::Mat,arg2::Ptr{PetscScalar})
    ccall((:MatMFFDGetH,petsc),PetscErrorCode,(Mat,Ptr{PetscScalar}),arg1,arg2)
end

function MatMFFDSetOptionsPrefix(arg1::Mat,arg2::Ptr{Uint8})
    ccall((:MatMFFDSetOptionsPrefix,petsc),PetscErrorCode,(Mat,Ptr{Uint8}),arg1,arg2)
end

function MatMFFDCheckPositivity(arg1::Ptr{Void},arg2::Vec,arg3::Vec,arg4::Ptr{PetscScalar})
    ccall((:MatMFFDCheckPositivity,petsc),PetscErrorCode,(Ptr{Void},Vec,Vec,Ptr{PetscScalar}),arg1,arg2,arg3,arg4)
end

function MatMFFDSetCheckh(arg1::Mat,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:MatMFFDSetCheckh,petsc),PetscErrorCode,(Mat,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function MatMFFDSetType(arg1::Mat,arg2::MatMFFDType)
    ccall((:MatMFFDSetType,petsc),PetscErrorCode,(Mat,MatMFFDType),arg1,arg2)
end

function MatMFFDRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:MatMFFDRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function MatMFFDDSSetUmin(arg1::Mat,PetscReal::Cint)
    ccall((:MatMFFDDSSetUmin,petsc),PetscErrorCode,(Mat,Cint),arg1,PetscReal)
end

function MatMFFDWPSetComputeNormU(arg1::Mat,arg2::PetscBool)
    ccall((:MatMFFDWPSetComputeNormU,petsc),PetscErrorCode,(Mat,PetscBool),arg1,arg2)
end

function PetscViewerMathematicaPutMatrix(arg1::PetscViewer,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{Cint})
    ccall((:PetscViewerMathematicaPutMatrix,petsc),PetscErrorCode,(PetscViewer,PetscInt,PetscInt,Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function PetscViewerMathematicaPutCSRMatrix(arg1::PetscViewer,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscInt},arg5::Ptr{PetscInt},arg6::Ptr{Cint})
    ccall((:PetscViewerMathematicaPutCSRMatrix,petsc),PetscErrorCode,(PetscViewer,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatCreateNest(arg1::MPI_Comm,arg2::PetscInt,arg3::Ptr{IS},arg4::PetscInt,arg5::Ptr{IS},arg6::Ptr{Mat},arg7::Ptr{Mat})
    ccall((:MatCreateNest,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{IS},PetscInt,Ptr{IS},Ptr{Mat},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatNestGetSize(arg1::Mat,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt})
    ccall((:MatNestGetSize,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatNestGetISs(arg1::Mat,arg2::Ptr{IS},arg3::Ptr{IS})
    ccall((:MatNestGetISs,petsc),PetscErrorCode,(Mat,Ptr{IS},Ptr{IS}),arg1,arg2,arg3)
end

function MatNestGetLocalISs(arg1::Mat,arg2::Ptr{IS},arg3::Ptr{IS})
    ccall((:MatNestGetLocalISs,petsc),PetscErrorCode,(Mat,Ptr{IS},Ptr{IS}),arg1,arg2,arg3)
end

function MatNestGetSubMats(arg1::Mat,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{Ptr{Ptr{Mat}}})
    ccall((:MatNestGetSubMats,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{Ptr{Mat}}}),arg1,arg2,arg3,arg4)
end

function MatNestGetSubMat(arg1::Mat,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{Mat})
    ccall((:MatNestGetSubMat,petsc),PetscErrorCode,(Mat,PetscInt,PetscInt,Ptr{Mat}),arg1,arg2,arg3,arg4)
end

function MatNestSetVecType(arg1::Mat,arg2::VecType)
    ccall((:MatNestSetVecType,petsc),PetscErrorCode,(Mat,VecType),arg1,arg2)
end

function MatNestSetSubMats(arg1::Mat,arg2::PetscInt,arg3::Ptr{IS},arg4::PetscInt,arg5::Ptr{IS},arg6::Ptr{Mat})
    ccall((:MatNestSetSubMats,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{IS},PetscInt,Ptr{IS},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatNestSetSubMat(arg1::Mat,arg2::PetscInt,arg3::PetscInt,arg4::Mat)
    ccall((:MatNestSetSubMat,petsc),PetscErrorCode,(Mat,PetscInt,PetscInt,Mat),arg1,arg2,arg3,arg4)
end

function MatChop(arg1::Mat,PetscReal::Cint)
    ccall((:MatChop,petsc),PetscErrorCode,(Mat,Cint),arg1,PetscReal)
end

function MatComputeBandwidth(arg1::Mat,PetscReal::Cint,arg2::Ptr{PetscInt})
    ccall((:MatComputeBandwidth,petsc),PetscErrorCode,(Mat,Cint,Ptr{PetscInt}),arg1,PetscReal,arg2)
end

function MatSubdomainsCreateCoalesce(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{Ptr{IS}})
    ccall((:MatSubdomainsCreateCoalesce,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{Ptr{IS}}),arg1,arg2,arg3,arg4)
end

function DMInitializePackage()
    ccall((:DMInitializePackage,petsc),PetscErrorCode,())
end

function DMCreate(arg1::MPI_Comm,arg2::Ptr{DM})
    ccall((:DMCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{DM}),arg1,arg2)
end

function DMClone(arg1::DM,arg2::Ptr{DM})
    ccall((:DMClone,petsc),PetscErrorCode,(DM,Ptr{DM}),arg1,arg2)
end

function DMSetType(arg1::DM,arg2::DMType)
    ccall((:DMSetType,petsc),PetscErrorCode,(DM,DMType),arg1,arg2)
end

function DMGetType(arg1::DM,arg2::Ptr{DMType})
    ccall((:DMGetType,petsc),PetscErrorCode,(DM,Ptr{DMType}),arg1,arg2)
end

function DMRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:DMRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function DMRegisterDestroy()
    ccall((:DMRegisterDestroy,petsc),PetscErrorCode,())
end

function DMView(arg1::DM,arg2::PetscViewer)
    ccall((:DMView,petsc),PetscErrorCode,(DM,PetscViewer),arg1,arg2)
end

function DMLoad(arg1::DM,arg2::PetscViewer)
    ccall((:DMLoad,petsc),PetscErrorCode,(DM,PetscViewer),arg1,arg2)
end

function DMDestroy(arg1::Ptr{DM})
    ccall((:DMDestroy,petsc),PetscErrorCode,(Ptr{DM},),arg1)
end

function DMCreateGlobalVector(arg1::DM,arg2::Ptr{Vec})
    ccall((:DMCreateGlobalVector,petsc),PetscErrorCode,(DM,Ptr{Vec}),arg1,arg2)
end

function DMCreateLocalVector(arg1::DM,arg2::Ptr{Vec})
    ccall((:DMCreateLocalVector,petsc),PetscErrorCode,(DM,Ptr{Vec}),arg1,arg2)
end

function DMGetLocalVector(arg1::DM,arg2::Ptr{Vec})
    ccall((:DMGetLocalVector,petsc),PetscErrorCode,(DM,Ptr{Vec}),arg1,arg2)
end

function DMRestoreLocalVector(arg1::DM,arg2::Ptr{Vec})
    ccall((:DMRestoreLocalVector,petsc),PetscErrorCode,(DM,Ptr{Vec}),arg1,arg2)
end

function DMGetGlobalVector(arg1::DM,arg2::Ptr{Vec})
    ccall((:DMGetGlobalVector,petsc),PetscErrorCode,(DM,Ptr{Vec}),arg1,arg2)
end

function DMRestoreGlobalVector(arg1::DM,arg2::Ptr{Vec})
    ccall((:DMRestoreGlobalVector,petsc),PetscErrorCode,(DM,Ptr{Vec}),arg1,arg2)
end

function DMClearGlobalVectors(arg1::DM)
    ccall((:DMClearGlobalVectors,petsc),PetscErrorCode,(DM,),arg1)
end

function DMGetNamedGlobalVector(arg1::DM,arg2::Ptr{Uint8},arg3::Ptr{Vec})
    ccall((:DMGetNamedGlobalVector,petsc),PetscErrorCode,(DM,Ptr{Uint8},Ptr{Vec}),arg1,arg2,arg3)
end

function DMRestoreNamedGlobalVector(arg1::DM,arg2::Ptr{Uint8},arg3::Ptr{Vec})
    ccall((:DMRestoreNamedGlobalVector,petsc),PetscErrorCode,(DM,Ptr{Uint8},Ptr{Vec}),arg1,arg2,arg3)
end

function DMGetNamedLocalVector(arg1::DM,arg2::Ptr{Uint8},arg3::Ptr{Vec})
    ccall((:DMGetNamedLocalVector,petsc),PetscErrorCode,(DM,Ptr{Uint8},Ptr{Vec}),arg1,arg2,arg3)
end

function DMRestoreNamedLocalVector(arg1::DM,arg2::Ptr{Uint8},arg3::Ptr{Vec})
    ccall((:DMRestoreNamedLocalVector,petsc),PetscErrorCode,(DM,Ptr{Uint8},Ptr{Vec}),arg1,arg2,arg3)
end

function DMGetLocalToGlobalMapping(arg1::DM,arg2::Ptr{ISLocalToGlobalMapping})
    ccall((:DMGetLocalToGlobalMapping,petsc),PetscErrorCode,(DM,Ptr{ISLocalToGlobalMapping}),arg1,arg2)
end

function DMCreateFieldIS(arg1::DM,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{Ptr{Uint8}}},arg4::Ptr{Ptr{IS}})
    ccall((:DMCreateFieldIS,petsc),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{Ptr{Ptr{Uint8}}},Ptr{Ptr{IS}}),arg1,arg2,arg3,arg4)
end

function DMGetBlockSize(arg1::DM,arg2::Ptr{PetscInt})
    ccall((:DMGetBlockSize,petsc),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end

function DMCreateColoring(arg1::DM,arg2::ISColoringType,arg3::Ptr{ISColoring})
    ccall((:DMCreateColoring,petsc),PetscErrorCode,(DM,ISColoringType,Ptr{ISColoring}),arg1,arg2,arg3)
end

function DMCreateMatrix(arg1::DM,arg2::Ptr{Mat})
    ccall((:DMCreateMatrix,petsc),PetscErrorCode,(DM,Ptr{Mat}),arg1,arg2)
end

function DMSetMatrixPreallocateOnly(arg1::DM,arg2::PetscBool)
    ccall((:DMSetMatrixPreallocateOnly,petsc),PetscErrorCode,(DM,PetscBool),arg1,arg2)
end

function DMCreateInterpolation(arg1::DM,arg2::DM,arg3::Ptr{Mat},arg4::Ptr{Vec})
    ccall((:DMCreateInterpolation,petsc),PetscErrorCode,(DM,DM,Ptr{Mat},Ptr{Vec}),arg1,arg2,arg3,arg4)
end

function DMCreateInjection(arg1::DM,arg2::DM,arg3::Ptr{Mat})
    ccall((:DMCreateInjection,petsc),PetscErrorCode,(DM,DM,Ptr{Mat}),arg1,arg2,arg3)
end

function DMGetWorkArray(arg1::DM,arg2::PetscInt,arg3::PetscDataType,arg4::Ptr{Void})
    ccall((:DMGetWorkArray,petsc),PetscErrorCode,(DM,PetscInt,PetscDataType,Ptr{Void}),arg1,arg2,arg3,arg4)
end

function DMRestoreWorkArray(arg1::DM,arg2::PetscInt,arg3::PetscDataType,arg4::Ptr{Void})
    ccall((:DMRestoreWorkArray,petsc),PetscErrorCode,(DM,PetscInt,PetscDataType,Ptr{Void}),arg1,arg2,arg3,arg4)
end

function DMRefine(arg1::DM,arg2::MPI_Comm,arg3::Ptr{DM})
    ccall((:DMRefine,petsc),PetscErrorCode,(DM,MPI_Comm,Ptr{DM}),arg1,arg2,arg3)
end

function DMCoarsen(arg1::DM,arg2::MPI_Comm,arg3::Ptr{DM})
    ccall((:DMCoarsen,petsc),PetscErrorCode,(DM,MPI_Comm,Ptr{DM}),arg1,arg2,arg3)
end

function DMRefineHierarchy(arg1::DM,arg2::PetscInt,arg3::Ptr{DM})
    ccall((:DMRefineHierarchy,petsc),PetscErrorCode,(DM,PetscInt,Ptr{DM}),arg1,arg2,arg3)
end

function DMCoarsenHierarchy(arg1::DM,arg2::PetscInt,arg3::Ptr{DM})
    ccall((:DMCoarsenHierarchy,petsc),PetscErrorCode,(DM,PetscInt,Ptr{DM}),arg1,arg2,arg3)
end

function DMCoarsenHookAdd(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:DMCoarsenHookAdd,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function DMRefineHookAdd(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:DMRefineHookAdd,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function DMRestrict(arg1::DM,arg2::Mat,arg3::Vec,arg4::Mat,arg5::DM)
    ccall((:DMRestrict,petsc),PetscErrorCode,(DM,Mat,Vec,Mat,DM),arg1,arg2,arg3,arg4,arg5)
end

function DMInterpolate(arg1::DM,arg2::Mat,arg3::DM)
    ccall((:DMInterpolate,petsc),PetscErrorCode,(DM,Mat,DM),arg1,arg2,arg3)
end

function DMSetFromOptions(arg1::DM)
    ccall((:DMSetFromOptions,petsc),PetscErrorCode,(DM,),arg1)
end

function DMCreateInterpolationScale(arg1::DM,arg2::DM,arg3::Mat,arg4::Ptr{Vec})
    ccall((:DMCreateInterpolationScale,petsc),PetscErrorCode,(DM,DM,Mat,Ptr{Vec}),arg1,arg2,arg3,arg4)
end

function DMCreateAggregates(arg1::DM,arg2::DM,arg3::Ptr{Mat})
    ccall((:DMCreateAggregates,petsc),PetscErrorCode,(DM,DM,Ptr{Mat}),arg1,arg2,arg3)
end

function DMGlobalToLocalHookAdd(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:DMGlobalToLocalHookAdd,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function DMLocalToGlobalHookAdd(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:DMLocalToGlobalHookAdd,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function DMGlobalToLocalBegin(arg1::DM,arg2::Vec,arg3::InsertMode,arg4::Vec)
    ccall((:DMGlobalToLocalBegin,petsc),PetscErrorCode,(DM,Vec,InsertMode,Vec),arg1,arg2,arg3,arg4)
end

function DMGlobalToLocalEnd(arg1::DM,arg2::Vec,arg3::InsertMode,arg4::Vec)
    ccall((:DMGlobalToLocalEnd,petsc),PetscErrorCode,(DM,Vec,InsertMode,Vec),arg1,arg2,arg3,arg4)
end

function DMLocalToGlobalBegin(arg1::DM,arg2::Vec,arg3::InsertMode,arg4::Vec)
    ccall((:DMLocalToGlobalBegin,petsc),PetscErrorCode,(DM,Vec,InsertMode,Vec),arg1,arg2,arg3,arg4)
end

function DMLocalToGlobalEnd(arg1::DM,arg2::Vec,arg3::InsertMode,arg4::Vec)
    ccall((:DMLocalToGlobalEnd,petsc),PetscErrorCode,(DM,Vec,InsertMode,Vec),arg1,arg2,arg3,arg4)
end

function DMLocalToLocalBegin(arg1::DM,arg2::Vec,arg3::InsertMode,arg4::Vec)
    ccall((:DMLocalToLocalBegin,petsc),PetscErrorCode,(DM,Vec,InsertMode,Vec),arg1,arg2,arg3,arg4)
end

function DMLocalToLocalEnd(arg1::DM,arg2::Vec,arg3::InsertMode,arg4::Vec)
    ccall((:DMLocalToLocalEnd,petsc),PetscErrorCode,(DM,Vec,InsertMode,Vec),arg1,arg2,arg3,arg4)
end

function DMConvert(arg1::DM,arg2::DMType,arg3::Ptr{DM})
    ccall((:DMConvert,petsc),PetscErrorCode,(DM,DMType,Ptr{DM}),arg1,arg2,arg3)
end

function DMGetDimension(arg1::DM,arg2::Ptr{PetscInt})
    ccall((:DMGetDimension,petsc),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end

function DMSetDimension(arg1::DM,arg2::PetscInt)
    ccall((:DMSetDimension,petsc),PetscErrorCode,(DM,PetscInt),arg1,arg2)
end

function DMGetDimPoints(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:DMGetDimPoints,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function DMGetCoordinateDM(arg1::DM,arg2::Ptr{DM})
    ccall((:DMGetCoordinateDM,petsc),PetscErrorCode,(DM,Ptr{DM}),arg1,arg2)
end

function DMSetCoordinateDM(arg1::DM,arg2::DM)
    ccall((:DMSetCoordinateDM,petsc),PetscErrorCode,(DM,DM),arg1,arg2)
end

function DMGetCoordinateDim(arg1::DM,arg2::Ptr{PetscInt})
    ccall((:DMGetCoordinateDim,petsc),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end

function DMSetCoordinateDim(arg1::DM,arg2::PetscInt)
    ccall((:DMSetCoordinateDim,petsc),PetscErrorCode,(DM,PetscInt),arg1,arg2)
end

function DMGetCoordinateSection(arg1::DM,arg2::Ptr{PetscSection})
    ccall((:DMGetCoordinateSection,petsc),PetscErrorCode,(DM,Ptr{PetscSection}),arg1,arg2)
end

function DMSetCoordinateSection(arg1::DM,arg2::PetscInt,arg3::PetscSection)
    ccall((:DMSetCoordinateSection,petsc),PetscErrorCode,(DM,PetscInt,PetscSection),arg1,arg2,arg3)
end

function DMGetCoordinates(arg1::DM,arg2::Ptr{Vec})
    ccall((:DMGetCoordinates,petsc),PetscErrorCode,(DM,Ptr{Vec}),arg1,arg2)
end

function DMSetCoordinates(arg1::DM,arg2::Vec)
    ccall((:DMSetCoordinates,petsc),PetscErrorCode,(DM,Vec),arg1,arg2)
end

function DMGetCoordinatesLocal(arg1::DM,arg2::Ptr{Vec})
    ccall((:DMGetCoordinatesLocal,petsc),PetscErrorCode,(DM,Ptr{Vec}),arg1,arg2)
end

function DMSetCoordinatesLocal(arg1::DM,arg2::Vec)
    ccall((:DMSetCoordinatesLocal,petsc),PetscErrorCode,(DM,Vec),arg1,arg2)
end

function DMLocatePoints(arg1::DM,arg2::Vec,arg3::Ptr{IS})
    ccall((:DMLocatePoints,petsc),PetscErrorCode,(DM,Vec,Ptr{IS}),arg1,arg2,arg3)
end

function DMGetPeriodicity(arg1::DM,arg2::Ptr{Ptr{Cint}},arg3::Ptr{Ptr{Cint}},arg4::Ptr{Ptr{DMBoundaryType}})
    ccall((:DMGetPeriodicity,petsc),PetscErrorCode,(DM,Ptr{Ptr{Cint}},Ptr{Ptr{Cint}},Ptr{Ptr{DMBoundaryType}}),arg1,arg2,arg3,arg4)
end

function DMSetPeriodicity(arg1::DM,PetscReal::Ptr{Cint},arg2::Ptr{Cint},arg3::Ptr{DMBoundaryType})
    ccall((:DMSetPeriodicity,petsc),PetscErrorCode,(DM,Ptr{Cint},Ptr{Cint},Ptr{DMBoundaryType}),arg1,PetscReal,arg2,arg3)
end

function DMSubDomainHookAdd(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:DMSubDomainHookAdd,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function DMSubDomainRestrict(arg1::DM,arg2::VecScatter,arg3::VecScatter,arg4::DM)
    ccall((:DMSubDomainRestrict,petsc),PetscErrorCode,(DM,VecScatter,VecScatter,DM),arg1,arg2,arg3,arg4)
end

function DMSetOptionsPrefix(arg1::DM,arg2::Ptr{Uint8})
    ccall((:DMSetOptionsPrefix,petsc),PetscErrorCode,(DM,Ptr{Uint8}),arg1,arg2)
end

function DMSetVecType(arg1::DM,arg2::VecType)
    ccall((:DMSetVecType,petsc),PetscErrorCode,(DM,VecType),arg1,arg2)
end

function DMGetVecType(arg1::DM,arg2::Ptr{VecType})
    ccall((:DMGetVecType,petsc),PetscErrorCode,(DM,Ptr{VecType}),arg1,arg2)
end

function DMSetMatType(arg1::DM,arg2::MatType)
    ccall((:DMSetMatType,petsc),PetscErrorCode,(DM,MatType),arg1,arg2)
end

function DMGetMatType(arg1::DM,arg2::Ptr{MatType})
    ccall((:DMGetMatType,petsc),PetscErrorCode,(DM,Ptr{MatType}),arg1,arg2)
end

function DMSetApplicationContext(arg1::DM,arg2::Ptr{Void})
    ccall((:DMSetApplicationContext,petsc),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end

function DMSetApplicationContextDestroy(arg1::DM,arg2::Ptr{Void})
    ccall((:DMSetApplicationContextDestroy,petsc),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end

function DMGetApplicationContext(arg1::DM,arg2::Ptr{Void})
    ccall((:DMGetApplicationContext,petsc),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end

function DMSetVariableBounds(arg1::DM,arg2::Ptr{Void})
    ccall((:DMSetVariableBounds,petsc),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end

function DMHasVariableBounds(arg1::DM,arg2::Ptr{PetscBool})
    ccall((:DMHasVariableBounds,petsc),PetscErrorCode,(DM,Ptr{PetscBool}),arg1,arg2)
end

function DMHasColoring(arg1::DM,arg2::Ptr{PetscBool})
    ccall((:DMHasColoring,petsc),PetscErrorCode,(DM,Ptr{PetscBool}),arg1,arg2)
end

function DMComputeVariableBounds(arg1::DM,arg2::Vec,arg3::Vec)
    ccall((:DMComputeVariableBounds,petsc),PetscErrorCode,(DM,Vec,Vec),arg1,arg2,arg3)
end

function DMCreateSubDM(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{IS},arg5::Ptr{DM})
    ccall((:DMCreateSubDM,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{IS},Ptr{DM}),arg1,arg2,arg3,arg4,arg5)
end

function DMCreateFieldDecomposition(arg1::DM,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{Ptr{Uint8}}},arg4::Ptr{Ptr{IS}},arg5::Ptr{Ptr{DM}})
    ccall((:DMCreateFieldDecomposition,petsc),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{Ptr{Ptr{Uint8}}},Ptr{Ptr{IS}},Ptr{Ptr{DM}}),arg1,arg2,arg3,arg4,arg5)
end

function DMCreateDomainDecomposition(arg1::DM,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{Ptr{Uint8}}},arg4::Ptr{Ptr{IS}},arg5::Ptr{Ptr{IS}},arg6::Ptr{Ptr{DM}})
    ccall((:DMCreateDomainDecomposition,petsc),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{Ptr{Ptr{Uint8}}},Ptr{Ptr{IS}},Ptr{Ptr{IS}},Ptr{Ptr{DM}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function DMCreateDomainDecompositionScatters(arg1::DM,arg2::PetscInt,arg3::Ptr{DM},arg4::Ptr{Ptr{VecScatter}},arg5::Ptr{Ptr{VecScatter}},arg6::Ptr{Ptr{VecScatter}})
    ccall((:DMCreateDomainDecompositionScatters,petsc),PetscErrorCode,(DM,PetscInt,Ptr{DM},Ptr{Ptr{VecScatter}},Ptr{Ptr{VecScatter}},Ptr{Ptr{VecScatter}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function DMGetRefineLevel(arg1::DM,arg2::Ptr{PetscInt})
    ccall((:DMGetRefineLevel,petsc),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end

function DMGetCoarsenLevel(arg1::DM,arg2::Ptr{PetscInt})
    ccall((:DMGetCoarsenLevel,petsc),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end

function DMFinalizePackage()
    ccall((:DMFinalizePackage,petsc),PetscErrorCode,())
end

function VecGetDM(arg1::Vec,arg2::Ptr{DM})
    ccall((:VecGetDM,petsc),PetscErrorCode,(Vec,Ptr{DM}),arg1,arg2)
end

function VecSetDM(arg1::Vec,arg2::DM)
    ccall((:VecSetDM,petsc),PetscErrorCode,(Vec,DM),arg1,arg2)
end

function MatGetDM(arg1::Mat,arg2::Ptr{DM})
    ccall((:MatGetDM,petsc),PetscErrorCode,(Mat,Ptr{DM}),arg1,arg2)
end

function MatSetDM(arg1::Mat,arg2::DM)
    ccall((:MatSetDM,petsc),PetscErrorCode,(Mat,DM),arg1,arg2)
end

function DMPrintCellVector(arg1::PetscInt,arg2::Ptr{Uint8},arg3::PetscInt,arg4::Ptr{PetscScalar})
    ccall((:DMPrintCellVector,petsc),PetscErrorCode,(PetscInt,Ptr{Uint8},PetscInt,Ptr{PetscScalar}),arg1,arg2,arg3,arg4)
end

function DMPrintCellMatrix(arg1::PetscInt,arg2::Ptr{Uint8},arg3::PetscInt,arg4::PetscInt,arg5::Ptr{PetscScalar})
    ccall((:DMPrintCellMatrix,petsc),PetscErrorCode,(PetscInt,Ptr{Uint8},PetscInt,PetscInt,Ptr{PetscScalar}),arg1,arg2,arg3,arg4,arg5)
end

function DMPrintLocalVec(arg1::DM,arg2::Ptr{Uint8},PetscReal::Cint,arg3::Vec)
    ccall((:DMPrintLocalVec,petsc),PetscErrorCode,(DM,Ptr{Uint8},Cint,Vec),arg1,arg2,PetscReal,arg3)
end

function DMGetDefaultSection(arg1::DM,arg2::Ptr{PetscSection})
    ccall((:DMGetDefaultSection,petsc),PetscErrorCode,(DM,Ptr{PetscSection}),arg1,arg2)
end

function DMSetDefaultSection(arg1::DM,arg2::PetscSection)
    ccall((:DMSetDefaultSection,petsc),PetscErrorCode,(DM,PetscSection),arg1,arg2)
end

function DMGetDefaultConstraints(arg1::DM,arg2::Ptr{PetscSection},arg3::Ptr{Mat})
    ccall((:DMGetDefaultConstraints,petsc),PetscErrorCode,(DM,Ptr{PetscSection},Ptr{Mat}),arg1,arg2,arg3)
end

function DMSetDefaultConstraints(arg1::DM,arg2::PetscSection,arg3::Mat)
    ccall((:DMSetDefaultConstraints,petsc),PetscErrorCode,(DM,PetscSection,Mat),arg1,arg2,arg3)
end

function DMGetDefaultGlobalSection(arg1::DM,arg2::Ptr{PetscSection})
    ccall((:DMGetDefaultGlobalSection,petsc),PetscErrorCode,(DM,Ptr{PetscSection}),arg1,arg2)
end

function DMSetDefaultGlobalSection(arg1::DM,arg2::PetscSection)
    ccall((:DMSetDefaultGlobalSection,petsc),PetscErrorCode,(DM,PetscSection),arg1,arg2)
end

function DMGetDefaultSF(arg1::DM,arg2::Ptr{PetscSF})
    ccall((:DMGetDefaultSF,petsc),PetscErrorCode,(DM,Ptr{PetscSF}),arg1,arg2)
end

function DMSetDefaultSF(arg1::DM,arg2::PetscSF)
    ccall((:DMSetDefaultSF,petsc),PetscErrorCode,(DM,PetscSF),arg1,arg2)
end

function DMCreateDefaultSF(arg1::DM,arg2::PetscSection,arg3::PetscSection)
    ccall((:DMCreateDefaultSF,petsc),PetscErrorCode,(DM,PetscSection,PetscSection),arg1,arg2,arg3)
end

function DMGetPointSF(arg1::DM,arg2::Ptr{PetscSF})
    ccall((:DMGetPointSF,petsc),PetscErrorCode,(DM,Ptr{PetscSF}),arg1,arg2)
end

function DMSetPointSF(arg1::DM,arg2::PetscSF)
    ccall((:DMSetPointSF,petsc),PetscErrorCode,(DM,PetscSF),arg1,arg2)
end

function DMGetOutputDM(arg1::DM,arg2::Ptr{DM})
    ccall((:DMGetOutputDM,petsc),PetscErrorCode,(DM,Ptr{DM}),arg1,arg2)
end

function DMGetOutputSequenceNumber(arg1::DM,arg2::Ptr{PetscInt},arg3::Ptr{Cint})
    ccall((:DMGetOutputSequenceNumber,petsc),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{Cint}),arg1,arg2,arg3)
end

function DMSetOutputSequenceNumber(arg1::DM,arg2::PetscInt,PetscReal::Cint)
    ccall((:DMSetOutputSequenceNumber,petsc),PetscErrorCode,(DM,PetscInt,Cint),arg1,arg2,PetscReal)
end

function DMOutputSequenceLoad(arg1::DM,arg2::PetscViewer,arg3::Ptr{Uint8},arg4::PetscInt,arg5::Ptr{Cint})
    ccall((:DMOutputSequenceLoad,petsc),PetscErrorCode,(DM,PetscViewer,Ptr{Uint8},PetscInt,Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
end

function DMGetDS(arg1::DM,arg2::Ptr{PetscDS})
    ccall((:DMGetDS,petsc),PetscErrorCode,(DM,Ptr{PetscDS}),arg1,arg2)
end

function DMSetDS(arg1::DM,arg2::PetscDS)
    ccall((:DMSetDS,petsc),PetscErrorCode,(DM,PetscDS),arg1,arg2)
end

function DMGetNumFields(arg1::DM,arg2::Ptr{PetscInt})
    ccall((:DMGetNumFields,petsc),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end

function DMSetNumFields(arg1::DM,arg2::PetscInt)
    ccall((:DMSetNumFields,petsc),PetscErrorCode,(DM,PetscInt),arg1,arg2)
end

function DMGetField(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscObject})
    ccall((:DMGetField,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscObject}),arg1,arg2,arg3)
end

function DMSetField(arg1::DM,arg2::PetscInt,arg3::PetscObject)
    ccall((:DMSetField,petsc),PetscErrorCode,(DM,PetscInt,PetscObject),arg1,arg2,arg3)
end

function DMInterpolationCreate(arg1::MPI_Comm,arg2::Ptr{DMInterpolationInfo})
    ccall((:DMInterpolationCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{DMInterpolationInfo}),arg1,arg2)
end

function DMInterpolationSetDim(arg1::DMInterpolationInfo,arg2::PetscInt)
    ccall((:DMInterpolationSetDim,petsc),PetscErrorCode,(DMInterpolationInfo,PetscInt),arg1,arg2)
end

function DMInterpolationGetDim(arg1::DMInterpolationInfo,arg2::Ptr{PetscInt})
    ccall((:DMInterpolationGetDim,petsc),PetscErrorCode,(DMInterpolationInfo,Ptr{PetscInt}),arg1,arg2)
end

function DMInterpolationSetDof(arg1::DMInterpolationInfo,arg2::PetscInt)
    ccall((:DMInterpolationSetDof,petsc),PetscErrorCode,(DMInterpolationInfo,PetscInt),arg1,arg2)
end

function DMInterpolationGetDof(arg1::DMInterpolationInfo,arg2::Ptr{PetscInt})
    ccall((:DMInterpolationGetDof,petsc),PetscErrorCode,(DMInterpolationInfo,Ptr{PetscInt}),arg1,arg2)
end

function DMInterpolationAddPoints(arg1::DMInterpolationInfo,arg2::PetscInt,PetscReal::Ptr{Cint})
    ccall((:DMInterpolationAddPoints,petsc),PetscErrorCode,(DMInterpolationInfo,PetscInt,Ptr{Cint}),arg1,arg2,PetscReal)
end

function DMInterpolationSetUp(arg1::DMInterpolationInfo,arg2::DM,arg3::PetscBool)
    ccall((:DMInterpolationSetUp,petsc),PetscErrorCode,(DMInterpolationInfo,DM,PetscBool),arg1,arg2,arg3)
end

function DMInterpolationGetCoordinates(arg1::DMInterpolationInfo,arg2::Ptr{Vec})
    ccall((:DMInterpolationGetCoordinates,petsc),PetscErrorCode,(DMInterpolationInfo,Ptr{Vec}),arg1,arg2)
end

function DMInterpolationGetVector(arg1::DMInterpolationInfo,arg2::Ptr{Vec})
    ccall((:DMInterpolationGetVector,petsc),PetscErrorCode,(DMInterpolationInfo,Ptr{Vec}),arg1,arg2)
end

function DMInterpolationRestoreVector(arg1::DMInterpolationInfo,arg2::Ptr{Vec})
    ccall((:DMInterpolationRestoreVector,petsc),PetscErrorCode,(DMInterpolationInfo,Ptr{Vec}),arg1,arg2)
end

function DMInterpolationEvaluate(arg1::DMInterpolationInfo,arg2::DM,arg3::Vec,arg4::Vec)
    ccall((:DMInterpolationEvaluate,petsc),PetscErrorCode,(DMInterpolationInfo,DM,Vec,Vec),arg1,arg2,arg3,arg4)
end

function DMInterpolationDestroy(arg1::Ptr{DMInterpolationInfo})
    ccall((:DMInterpolationDestroy,petsc),PetscErrorCode,(Ptr{DMInterpolationInfo},),arg1)
end

function PFCreate(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PF})
    ccall((:PFCreate,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,Ptr{PF}),arg1,arg2,arg3,arg4)
end

function PFSetType(arg1::PF,arg2::PFType,arg3::Ptr{Void})
    ccall((:PFSetType,petsc),PetscErrorCode,(PF,PFType,Ptr{Void}),arg1,arg2,arg3)
end

function PFSet(arg1::PF,arg2::Ptr{Void},arg3::Ptr{Void},arg4::Ptr{Void},arg5::Ptr{Void},arg6::Ptr{Void})
    ccall((:PFSet,petsc),PetscErrorCode,(PF,Ptr{Void},Ptr{Void},Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PFApply(arg1::PF,arg2::PetscInt,arg3::Ptr{PetscScalar},arg4::Ptr{PetscScalar})
    ccall((:PFApply,petsc),PetscErrorCode,(PF,PetscInt,Ptr{PetscScalar},Ptr{PetscScalar}),arg1,arg2,arg3,arg4)
end

function PFApplyVec(arg1::PF,arg2::Vec,arg3::Vec)
    ccall((:PFApplyVec,petsc),PetscErrorCode,(PF,Vec,Vec),arg1,arg2,arg3)
end

function PFInitializePackage()
    ccall((:PFInitializePackage,petsc),PetscErrorCode,())
end

function PFRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:PFRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function PFDestroy(arg1::Ptr{PF})
    ccall((:PFDestroy,petsc),PetscErrorCode,(Ptr{PF},),arg1)
end

function PFSetFromOptions(arg1::PF)
    ccall((:PFSetFromOptions,petsc),PetscErrorCode,(PF,),arg1)
end

function PFGetType(arg1::PF,arg2::Ptr{PFType})
    ccall((:PFGetType,petsc),PetscErrorCode,(PF,Ptr{PFType}),arg1,arg2)
end

function PFView(arg1::PF,arg2::PetscViewer)
    ccall((:PFView,petsc),PetscErrorCode,(PF,PetscViewer),arg1,arg2)
end

function AOInitializePackage()
    ccall((:AOInitializePackage,petsc),PetscErrorCode,())
end

function AOCreate(arg1::MPI_Comm,arg2::Ptr{Cint})
    ccall((:AOCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{Cint}),arg1,arg2)
end

function AOSetIS()
    ccall((:AOSetIS,petsc),PetscErrorCode,())
end

function AOSetFromOptions()
    ccall((:AOSetFromOptions,petsc),PetscErrorCode,())
end

function AOCreateBasic(arg1::MPI_Comm,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{Cint})
    ccall((:AOCreateBasic,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
end

function AOCreateBasicIS(arg1::IS,arg2::IS,arg3::Ptr{Cint})
    ccall((:AOCreateBasicIS,petsc),PetscErrorCode,(IS,IS,Ptr{Cint}),arg1,arg2,arg3)
end

function AOCreateMemoryScalable(arg1::MPI_Comm,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{Cint})
    ccall((:AOCreateMemoryScalable,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
end

function AOCreateMemoryScalableIS(arg1::IS,arg2::IS,arg3::Ptr{Cint})
    ccall((:AOCreateMemoryScalableIS,petsc),PetscErrorCode,(IS,IS,Ptr{Cint}),arg1,arg2,arg3)
end

function AOCreateMapping(arg1::MPI_Comm,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{Cint})
    ccall((:AOCreateMapping,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
end

function AOCreateMappingIS(arg1::IS,arg2::IS,arg3::Ptr{Cint})
    ccall((:AOCreateMappingIS,petsc),PetscErrorCode,(IS,IS,Ptr{Cint}),arg1,arg2,arg3)
end

function AOView()
    ccall((:AOView,petsc),PetscErrorCode,())
end

function AOSetType()
    ccall((:AOSetType,petsc),PetscErrorCode,())
end

function AOGetType()
    ccall((:AOGetType,petsc),PetscErrorCode,())
end

function AORegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:AORegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function AOPetscToApplication()
    ccall((:AOPetscToApplication,petsc),PetscErrorCode,())
end

function AOApplicationToPetsc()
    ccall((:AOApplicationToPetsc,petsc),PetscErrorCode,())
end

function AOPetscToApplicationIS()
    ccall((:AOPetscToApplicationIS,petsc),PetscErrorCode,())
end

function AOApplicationToPetscIS()
    ccall((:AOApplicationToPetscIS,petsc),PetscErrorCode,())
end

function AOPetscToApplicationPermuteInt()
    ccall((:AOPetscToApplicationPermuteInt,petsc),PetscErrorCode,())
end

function AOApplicationToPetscPermuteInt()
    ccall((:AOApplicationToPetscPermuteInt,petsc),PetscErrorCode,())
end

function AOPetscToApplicationPermuteReal()
    ccall((:AOPetscToApplicationPermuteReal,petsc),PetscErrorCode,())
end

function AOApplicationToPetscPermuteReal()
    ccall((:AOApplicationToPetscPermuteReal,petsc),PetscErrorCode,())
end

function AOMappingHasApplicationIndex()
    ccall((:AOMappingHasApplicationIndex,petsc),PetscErrorCode,())
end

function AOMappingHasPetscIndex()
    ccall((:AOMappingHasPetscIndex,petsc),PetscErrorCode,())
end

function PetscQuadratureCreate(arg1::MPI_Comm,arg2::Ptr{PetscQuadrature})
    ccall((:PetscQuadratureCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscQuadrature}),arg1,arg2)
end

function PetscQuadratureDuplicate(arg1::PetscQuadrature,arg2::Ptr{PetscQuadrature})
    ccall((:PetscQuadratureDuplicate,petsc),PetscErrorCode,(PetscQuadrature,Ptr{PetscQuadrature}),arg1,arg2)
end

function PetscQuadratureGetOrder(arg1::PetscQuadrature,arg2::Ptr{PetscInt})
    ccall((:PetscQuadratureGetOrder,petsc),PetscErrorCode,(PetscQuadrature,Ptr{PetscInt}),arg1,arg2)
end

function PetscQuadratureSetOrder(arg1::PetscQuadrature,arg2::PetscInt)
    ccall((:PetscQuadratureSetOrder,petsc),PetscErrorCode,(PetscQuadrature,PetscInt),arg1,arg2)
end

function PetscQuadratureGetData(arg1::PetscQuadrature,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{Ptr{Cint}},arg5::Ptr{Ptr{Cint}})
    ccall((:PetscQuadratureGetData,petsc),PetscErrorCode,(PetscQuadrature,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}}),arg1,arg2,arg3,arg4,arg5)
end

function PetscQuadratureSetData(arg1::PetscQuadrature,arg2::PetscInt,arg3::PetscInt,PetscReal::Ptr{Cint},arg4::Ptr{Cint})
    ccall((:PetscQuadratureSetData,petsc),PetscErrorCode,(PetscQuadrature,PetscInt,PetscInt,Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,PetscReal,arg4)
end

function PetscQuadratureView(arg1::PetscQuadrature,arg2::PetscViewer)
    ccall((:PetscQuadratureView,petsc),PetscErrorCode,(PetscQuadrature,PetscViewer),arg1,arg2)
end

function PetscQuadratureDestroy(arg1::Ptr{PetscQuadrature})
    ccall((:PetscQuadratureDestroy,petsc),PetscErrorCode,(Ptr{PetscQuadrature},),arg1)
end

function PetscQuadratureExpandComposite(arg1::PetscQuadrature,arg2::PetscInt,PetscReal::Ptr{Cint},arg3::Ptr{Cint},arg4::Ptr{PetscQuadrature})
    ccall((:PetscQuadratureExpandComposite,petsc),PetscErrorCode,(PetscQuadrature,PetscInt,Ptr{Cint},Ptr{Cint},Ptr{PetscQuadrature}),arg1,arg2,PetscReal,arg3,arg4)
end

function PetscDTLegendreEval(arg1::PetscInt,arg2::Ptr{Cint},arg3::PetscInt,arg4::Ptr{PetscInt},arg5::Ptr{Cint},arg6::Ptr{Cint},arg7::Ptr{Cint})
    ccall((:PetscDTLegendreEval,petsc),PetscErrorCode,(PetscInt,Ptr{Cint},PetscInt,Ptr{PetscInt},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscDTGaussQuadrature(arg1::PetscInt,PetscReal::Cint,arg2::Cint,arg3::Ptr{Cint},arg4::Ptr{Cint})
    ccall((:PetscDTGaussQuadrature,petsc),PetscErrorCode,(PetscInt,Cint,Cint,Ptr{Cint},Ptr{Cint}),arg1,PetscReal,arg2,arg3,arg4)
end

function PetscDTReconstructPoly(arg1::PetscInt,arg2::PetscInt,arg3::Ptr{Cint},arg4::PetscInt,arg5::Ptr{Cint},arg6::Ptr{Cint})
    ccall((:PetscDTReconstructPoly,petsc),PetscErrorCode,(PetscInt,PetscInt,Ptr{Cint},PetscInt,Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscDTGaussTensorQuadrature(arg1::PetscInt,arg2::PetscInt,PetscReal::Cint,arg3::Cint,arg4::Ptr{PetscQuadrature})
    ccall((:PetscDTGaussTensorQuadrature,petsc),PetscErrorCode,(PetscInt,PetscInt,Cint,Cint,Ptr{PetscQuadrature}),arg1,arg2,PetscReal,arg3,arg4)
end

function PetscDTGaussJacobiQuadrature(arg1::PetscInt,arg2::PetscInt,PetscReal::Cint,arg3::Cint,arg4::Ptr{PetscQuadrature})
    ccall((:PetscDTGaussJacobiQuadrature,petsc),PetscErrorCode,(PetscInt,PetscInt,Cint,Cint,Ptr{PetscQuadrature}),arg1,arg2,PetscReal,arg3,arg4)
end

function PetscFEInitializePackage()
    ccall((:PetscFEInitializePackage,petsc),PetscErrorCode,())
end

function PetscSpaceCreate(arg1::MPI_Comm,arg2::Ptr{PetscSpace})
    ccall((:PetscSpaceCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscSpace}),arg1,arg2)
end

function PetscSpaceDestroy(arg1::Ptr{PetscSpace})
    ccall((:PetscSpaceDestroy,petsc),PetscErrorCode,(Ptr{PetscSpace},),arg1)
end

function PetscSpaceSetType(arg1::PetscSpace,arg2::PetscSpaceType)
    ccall((:PetscSpaceSetType,petsc),PetscErrorCode,(PetscSpace,PetscSpaceType),arg1,arg2)
end

function PetscSpaceGetType(arg1::PetscSpace,arg2::Ptr{PetscSpaceType})
    ccall((:PetscSpaceGetType,petsc),PetscErrorCode,(PetscSpace,Ptr{PetscSpaceType}),arg1,arg2)
end

function PetscSpaceSetUp(arg1::PetscSpace)
    ccall((:PetscSpaceSetUp,petsc),PetscErrorCode,(PetscSpace,),arg1)
end

function PetscSpaceSetFromOptions(arg1::PetscSpace)
    ccall((:PetscSpaceSetFromOptions,petsc),PetscErrorCode,(PetscSpace,),arg1)
end

function PetscSpaceRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:PetscSpaceRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function PetscSpaceRegisterDestroy()
    ccall((:PetscSpaceRegisterDestroy,petsc),PetscErrorCode,())
end

function PetscSpaceGetDimension(arg1::PetscSpace,arg2::Ptr{PetscInt})
    ccall((:PetscSpaceGetDimension,petsc),PetscErrorCode,(PetscSpace,Ptr{PetscInt}),arg1,arg2)
end

function PetscSpaceSetOrder(arg1::PetscSpace,arg2::PetscInt)
    ccall((:PetscSpaceSetOrder,petsc),PetscErrorCode,(PetscSpace,PetscInt),arg1,arg2)
end

function PetscSpaceGetOrder(arg1::PetscSpace,arg2::Ptr{PetscInt})
    ccall((:PetscSpaceGetOrder,petsc),PetscErrorCode,(PetscSpace,Ptr{PetscInt}),arg1,arg2)
end

function PetscSpaceEvaluate(arg1::PetscSpace,arg2::PetscInt,PetscReal::Ptr{Cint},arg3::Ptr{Cint},arg4::Ptr{Cint},arg5::Ptr{Cint})
    ccall((:PetscSpaceEvaluate,petsc),PetscErrorCode,(PetscSpace,PetscInt,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,PetscReal,arg3,arg4,arg5)
end

function PetscSpacePolynomialSetNumVariables(arg1::PetscSpace,arg2::PetscInt)
    ccall((:PetscSpacePolynomialSetNumVariables,petsc),PetscErrorCode,(PetscSpace,PetscInt),arg1,arg2)
end

function PetscSpacePolynomialGetNumVariables(arg1::PetscSpace,arg2::Ptr{PetscInt})
    ccall((:PetscSpacePolynomialGetNumVariables,petsc),PetscErrorCode,(PetscSpace,Ptr{PetscInt}),arg1,arg2)
end

function PetscSpacePolynomialSetSymmetric(arg1::PetscSpace,arg2::PetscBool)
    ccall((:PetscSpacePolynomialSetSymmetric,petsc),PetscErrorCode,(PetscSpace,PetscBool),arg1,arg2)
end

function PetscSpacePolynomialGetSymmetric(arg1::PetscSpace,arg2::Ptr{PetscBool})
    ccall((:PetscSpacePolynomialGetSymmetric,petsc),PetscErrorCode,(PetscSpace,Ptr{PetscBool}),arg1,arg2)
end

function PetscSpacePolynomialSetTensor(arg1::PetscSpace,arg2::PetscBool)
    ccall((:PetscSpacePolynomialSetTensor,petsc),PetscErrorCode,(PetscSpace,PetscBool),arg1,arg2)
end

function PetscSpacePolynomialGetTensor(arg1::PetscSpace,arg2::Ptr{PetscBool})
    ccall((:PetscSpacePolynomialGetTensor,petsc),PetscErrorCode,(PetscSpace,Ptr{PetscBool}),arg1,arg2)
end

function PetscSpaceDGSetQuadrature(arg1::PetscSpace,arg2::PetscQuadrature)
    ccall((:PetscSpaceDGSetQuadrature,petsc),PetscErrorCode,(PetscSpace,PetscQuadrature),arg1,arg2)
end

function PetscSpaceDGGetQuadrature(arg1::PetscSpace,arg2::Ptr{PetscQuadrature})
    ccall((:PetscSpaceDGGetQuadrature,petsc),PetscErrorCode,(PetscSpace,Ptr{PetscQuadrature}),arg1,arg2)
end

function PetscDualSpaceCreate(arg1::MPI_Comm,arg2::Ptr{PetscDualSpace})
    ccall((:PetscDualSpaceCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscDualSpace}),arg1,arg2)
end

function PetscDualSpaceDestroy(arg1::Ptr{PetscDualSpace})
    ccall((:PetscDualSpaceDestroy,petsc),PetscErrorCode,(Ptr{PetscDualSpace},),arg1)
end

function PetscDualSpaceDuplicate(arg1::PetscDualSpace,arg2::Ptr{PetscDualSpace})
    ccall((:PetscDualSpaceDuplicate,petsc),PetscErrorCode,(PetscDualSpace,Ptr{PetscDualSpace}),arg1,arg2)
end

function PetscDualSpaceSetType(arg1::PetscDualSpace,arg2::PetscDualSpaceType)
    ccall((:PetscDualSpaceSetType,petsc),PetscErrorCode,(PetscDualSpace,PetscDualSpaceType),arg1,arg2)
end

function PetscDualSpaceGetType(arg1::PetscDualSpace,arg2::Ptr{PetscDualSpaceType})
    ccall((:PetscDualSpaceGetType,petsc),PetscErrorCode,(PetscDualSpace,Ptr{PetscDualSpaceType}),arg1,arg2)
end

function PetscDualSpaceSetUp(arg1::PetscDualSpace)
    ccall((:PetscDualSpaceSetUp,petsc),PetscErrorCode,(PetscDualSpace,),arg1)
end

function PetscDualSpaceSetFromOptions(arg1::PetscDualSpace)
    ccall((:PetscDualSpaceSetFromOptions,petsc),PetscErrorCode,(PetscDualSpace,),arg1)
end

function PetscDualSpaceRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:PetscDualSpaceRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function PetscDualSpaceRegisterDestroy()
    ccall((:PetscDualSpaceRegisterDestroy,petsc),PetscErrorCode,())
end

function PetscDualSpaceGetDimension(arg1::PetscDualSpace,arg2::Ptr{PetscInt})
    ccall((:PetscDualSpaceGetDimension,petsc),PetscErrorCode,(PetscDualSpace,Ptr{PetscInt}),arg1,arg2)
end

function PetscDualSpaceSetOrder(arg1::PetscDualSpace,arg2::PetscInt)
    ccall((:PetscDualSpaceSetOrder,petsc),PetscErrorCode,(PetscDualSpace,PetscInt),arg1,arg2)
end

function PetscDualSpaceGetOrder(arg1::PetscDualSpace,arg2::Ptr{PetscInt})
    ccall((:PetscDualSpaceGetOrder,petsc),PetscErrorCode,(PetscDualSpace,Ptr{PetscInt}),arg1,arg2)
end

function PetscDualSpaceSetDM(arg1::PetscDualSpace,arg2::DM)
    ccall((:PetscDualSpaceSetDM,petsc),PetscErrorCode,(PetscDualSpace,DM),arg1,arg2)
end

function PetscDualSpaceGetDM(arg1::PetscDualSpace,arg2::Ptr{DM})
    ccall((:PetscDualSpaceGetDM,petsc),PetscErrorCode,(PetscDualSpace,Ptr{DM}),arg1,arg2)
end

function PetscDualSpaceGetFunctional(arg1::PetscDualSpace,arg2::PetscInt,arg3::Ptr{PetscQuadrature})
    ccall((:PetscDualSpaceGetFunctional,petsc),PetscErrorCode,(PetscDualSpace,PetscInt,Ptr{PetscQuadrature}),arg1,arg2,arg3)
end

function PetscDualSpaceCreateReferenceCell(arg1::PetscDualSpace,arg2::PetscInt,arg3::PetscBool,arg4::Ptr{DM})
    ccall((:PetscDualSpaceCreateReferenceCell,petsc),PetscErrorCode,(PetscDualSpace,PetscInt,PetscBool,Ptr{DM}),arg1,arg2,arg3,arg4)
end

function PetscDualSpaceApply(arg1::PetscDualSpace,arg2::PetscInt,arg3::Ptr{PetscFECellGeom},arg4::PetscInt,arg5::Ptr{Void},arg6::Ptr{Void},arg7::Ptr{PetscScalar})
    ccall((:PetscDualSpaceApply,petsc),PetscErrorCode,(PetscDualSpace,PetscInt,Ptr{PetscFECellGeom},PetscInt,Ptr{Void},Ptr{Void},Ptr{PetscScalar}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscDualSpaceLagrangeGetContinuity(arg1::PetscDualSpace,arg2::Ptr{PetscBool})
    ccall((:PetscDualSpaceLagrangeGetContinuity,petsc),PetscErrorCode,(PetscDualSpace,Ptr{PetscBool}),arg1,arg2)
end

function PetscDualSpaceLagrangeSetContinuity(arg1::PetscDualSpace,arg2::PetscBool)
    ccall((:PetscDualSpaceLagrangeSetContinuity,petsc),PetscErrorCode,(PetscDualSpace,PetscBool),arg1,arg2)
end

function PetscDualSpaceGetHeightSubspace(arg1::PetscDualSpace,arg2::PetscInt,arg3::Ptr{PetscDualSpace})
    ccall((:PetscDualSpaceGetHeightSubspace,petsc),PetscErrorCode,(PetscDualSpace,PetscInt,Ptr{PetscDualSpace}),arg1,arg2,arg3)
end

function PetscDualSpaceSimpleSetDimension(arg1::PetscDualSpace,arg2::PetscInt)
    ccall((:PetscDualSpaceSimpleSetDimension,petsc),PetscErrorCode,(PetscDualSpace,PetscInt),arg1,arg2)
end

function PetscDualSpaceSimpleSetFunctional(arg1::PetscDualSpace,arg2::PetscInt,arg3::PetscQuadrature)
    ccall((:PetscDualSpaceSimpleSetFunctional,petsc),PetscErrorCode,(PetscDualSpace,PetscInt,PetscQuadrature),arg1,arg2,arg3)
end

function PetscFECreate(arg1::MPI_Comm,arg2::Ptr{PetscFE})
    ccall((:PetscFECreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscFE}),arg1,arg2)
end

function PetscFEDestroy(arg1::Ptr{PetscFE})
    ccall((:PetscFEDestroy,petsc),PetscErrorCode,(Ptr{PetscFE},),arg1)
end

function PetscFESetType(arg1::PetscFE,arg2::PetscFEType)
    ccall((:PetscFESetType,petsc),PetscErrorCode,(PetscFE,PetscFEType),arg1,arg2)
end

function PetscFEGetType(arg1::PetscFE,arg2::Ptr{PetscFEType})
    ccall((:PetscFEGetType,petsc),PetscErrorCode,(PetscFE,Ptr{PetscFEType}),arg1,arg2)
end

function PetscFESetUp(arg1::PetscFE)
    ccall((:PetscFESetUp,petsc),PetscErrorCode,(PetscFE,),arg1)
end

function PetscFESetFromOptions(arg1::PetscFE)
    ccall((:PetscFESetFromOptions,petsc),PetscErrorCode,(PetscFE,),arg1)
end

function PetscFERegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:PetscFERegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function PetscFERegisterDestroy()
    ccall((:PetscFERegisterDestroy,petsc),PetscErrorCode,())
end

function PetscFECreateDefault(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::PetscBool,arg5::Ptr{Uint8},arg6::PetscInt,arg7::Ptr{PetscFE})
    ccall((:PetscFECreateDefault,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,PetscBool,Ptr{Uint8},PetscInt,Ptr{PetscFE}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscFEGetDimension(arg1::PetscFE,arg2::Ptr{PetscInt})
    ccall((:PetscFEGetDimension,petsc),PetscErrorCode,(PetscFE,Ptr{PetscInt}),arg1,arg2)
end

function PetscFEGetSpatialDimension(arg1::PetscFE,arg2::Ptr{PetscInt})
    ccall((:PetscFEGetSpatialDimension,petsc),PetscErrorCode,(PetscFE,Ptr{PetscInt}),arg1,arg2)
end

function PetscFESetNumComponents(arg1::PetscFE,arg2::PetscInt)
    ccall((:PetscFESetNumComponents,petsc),PetscErrorCode,(PetscFE,PetscInt),arg1,arg2)
end

function PetscFEGetNumComponents(arg1::PetscFE,arg2::Ptr{PetscInt})
    ccall((:PetscFEGetNumComponents,petsc),PetscErrorCode,(PetscFE,Ptr{PetscInt}),arg1,arg2)
end

function PetscFEGetTileSizes(arg1::PetscFE,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{PetscInt})
    ccall((:PetscFEGetTileSizes,petsc),PetscErrorCode,(PetscFE,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end

function PetscFESetTileSizes(arg1::PetscFE,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt)
    ccall((:PetscFESetTileSizes,petsc),PetscErrorCode,(PetscFE,PetscInt,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4,arg5)
end

function PetscFESetBasisSpace(arg1::PetscFE,arg2::PetscSpace)
    ccall((:PetscFESetBasisSpace,petsc),PetscErrorCode,(PetscFE,PetscSpace),arg1,arg2)
end

function PetscFEGetBasisSpace(arg1::PetscFE,arg2::Ptr{PetscSpace})
    ccall((:PetscFEGetBasisSpace,petsc),PetscErrorCode,(PetscFE,Ptr{PetscSpace}),arg1,arg2)
end

function PetscFESetDualSpace(arg1::PetscFE,arg2::PetscDualSpace)
    ccall((:PetscFESetDualSpace,petsc),PetscErrorCode,(PetscFE,PetscDualSpace),arg1,arg2)
end

function PetscFEGetDualSpace(arg1::PetscFE,arg2::Ptr{PetscDualSpace})
    ccall((:PetscFEGetDualSpace,petsc),PetscErrorCode,(PetscFE,Ptr{PetscDualSpace}),arg1,arg2)
end

function PetscFESetQuadrature(arg1::PetscFE,arg2::PetscQuadrature)
    ccall((:PetscFESetQuadrature,petsc),PetscErrorCode,(PetscFE,PetscQuadrature),arg1,arg2)
end

function PetscFEGetQuadrature(arg1::PetscFE,arg2::Ptr{PetscQuadrature})
    ccall((:PetscFEGetQuadrature,petsc),PetscErrorCode,(PetscFE,Ptr{PetscQuadrature}),arg1,arg2)
end

function PetscFEGetNumDof(arg1::PetscFE,arg2::Ptr{Ptr{PetscInt}})
    ccall((:PetscFEGetNumDof,petsc),PetscErrorCode,(PetscFE,Ptr{Ptr{PetscInt}}),arg1,arg2)
end

function PetscFEGetDefaultTabulation(arg1::PetscFE,arg2::Ptr{Ptr{Cint}},arg3::Ptr{Ptr{Cint}},arg4::Ptr{Ptr{Cint}})
    ccall((:PetscFEGetDefaultTabulation,petsc),PetscErrorCode,(PetscFE,Ptr{Ptr{Cint}},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}}),arg1,arg2,arg3,arg4)
end

function PetscFEGetFaceTabulation(arg1::PetscFE,arg2::Ptr{Ptr{Cint}})
    ccall((:PetscFEGetFaceTabulation,petsc),PetscErrorCode,(PetscFE,Ptr{Ptr{Cint}}),arg1,arg2)
end

function PetscFEGetTabulation(arg1::PetscFE,arg2::PetscInt,PetscReal::Ptr{Cint},arg3::Ptr{Ptr{Cint}},arg4::Ptr{Ptr{Cint}},arg5::Ptr{Ptr{Cint}})
    ccall((:PetscFEGetTabulation,petsc),PetscErrorCode,(PetscFE,PetscInt,Ptr{Cint},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}}),arg1,arg2,PetscReal,arg3,arg4,arg5)
end

function PetscFERestoreTabulation(arg1::PetscFE,arg2::PetscInt,PetscReal::Ptr{Cint},arg3::Ptr{Ptr{Cint}},arg4::Ptr{Ptr{Cint}},arg5::Ptr{Ptr{Cint}})
    ccall((:PetscFERestoreTabulation,petsc),PetscErrorCode,(PetscFE,PetscInt,Ptr{Cint},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}}),arg1,arg2,PetscReal,arg3,arg4,arg5)
end

function PetscFERefine(arg1::PetscFE,arg2::Ptr{PetscFE})
    ccall((:PetscFERefine,petsc),PetscErrorCode,(PetscFE,Ptr{PetscFE}),arg1,arg2)
end

function PetscFEIntegrate(arg1::PetscFE,arg2::PetscDS,arg3::PetscInt,arg4::PetscInt,arg5::Ptr{PetscFECellGeom},arg6::Ptr{PetscScalar},arg7::PetscDS,arg8::Ptr{PetscScalar},PetscReal::Ptr{Cint})
    ccall((:PetscFEIntegrate,petsc),PetscErrorCode,(PetscFE,PetscDS,PetscInt,PetscInt,Ptr{PetscFECellGeom},Ptr{PetscScalar},PetscDS,Ptr{PetscScalar},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,PetscReal)
end

function PetscFEIntegrateResidual(arg1::PetscFE,arg2::PetscDS,arg3::PetscInt,arg4::PetscInt,arg5::Ptr{PetscFECellGeom},arg6::Ptr{PetscScalar},arg7::Ptr{PetscScalar},arg8::PetscDS,arg9::Ptr{PetscScalar},arg10::Ptr{PetscScalar})
    ccall((:PetscFEIntegrateResidual,petsc),PetscErrorCode,(PetscFE,PetscDS,PetscInt,PetscInt,Ptr{PetscFECellGeom},Ptr{PetscScalar},Ptr{PetscScalar},PetscDS,Ptr{PetscScalar},Ptr{PetscScalar}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function PetscFEIntegrateBdResidual(arg1::PetscFE,arg2::PetscDS,arg3::PetscInt,arg4::PetscInt,arg5::Ptr{PetscFECellGeom},arg6::Ptr{PetscScalar},arg7::Ptr{PetscScalar},arg8::PetscDS,arg9::Ptr{PetscScalar},arg10::Ptr{PetscScalar})
    ccall((:PetscFEIntegrateBdResidual,petsc),PetscErrorCode,(PetscFE,PetscDS,PetscInt,PetscInt,Ptr{PetscFECellGeom},Ptr{PetscScalar},Ptr{PetscScalar},PetscDS,Ptr{PetscScalar},Ptr{PetscScalar}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function PetscFEIntegrateJacobian(arg1::PetscFE,arg2::PetscDS,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::Ptr{PetscFECellGeom},arg7::Ptr{PetscScalar},arg8::Ptr{PetscScalar},arg9::PetscDS,arg10::Ptr{PetscScalar},arg11::Ptr{PetscScalar})
    ccall((:PetscFEIntegrateJacobian,petsc),PetscErrorCode,(PetscFE,PetscDS,PetscInt,PetscInt,PetscInt,Ptr{PetscFECellGeom},Ptr{PetscScalar},Ptr{PetscScalar},PetscDS,Ptr{PetscScalar},Ptr{PetscScalar}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
end

function PetscFEIntegrateBdJacobian(arg1::PetscFE,arg2::PetscDS,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::Ptr{PetscFECellGeom},arg7::Ptr{PetscScalar},arg8::Ptr{PetscScalar},arg9::PetscDS,arg10::Ptr{PetscScalar},arg11::Ptr{PetscScalar})
    ccall((:PetscFEIntegrateBdJacobian,petsc),PetscErrorCode,(PetscFE,PetscDS,PetscInt,PetscInt,PetscInt,Ptr{PetscFECellGeom},Ptr{PetscScalar},Ptr{PetscScalar},PetscDS,Ptr{PetscScalar},Ptr{PetscScalar}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
end

function PetscFECompositeGetMapping(arg1::PetscFE,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{Cint}},arg4::Ptr{Ptr{Cint}},arg5::Ptr{Ptr{Cint}})
    ccall((:PetscFECompositeGetMapping,petsc),PetscErrorCode,(PetscFE,Ptr{PetscInt},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}}),arg1,arg2,arg3,arg4,arg5)
end

function PetscFEOpenCLSetRealType(arg1::PetscFE,arg2::PetscDataType)
    ccall((:PetscFEOpenCLSetRealType,petsc),PetscErrorCode,(PetscFE,PetscDataType),arg1,arg2)
end

function PetscFEOpenCLGetRealType(arg1::PetscFE,arg2::Ptr{PetscDataType})
    ccall((:PetscFEOpenCLGetRealType,petsc),PetscErrorCode,(PetscFE,Ptr{PetscDataType}),arg1,arg2)
end

function DMDASetInterpolationType(arg1::DM,arg2::DMDAInterpolationType)
    ccall((:DMDASetInterpolationType,petsc),PetscErrorCode,(DM,DMDAInterpolationType),arg1,arg2)
end

function DMDAGetInterpolationType(arg1::DM,arg2::Ptr{DMDAInterpolationType})
    ccall((:DMDAGetInterpolationType,petsc),PetscErrorCode,(DM,Ptr{DMDAInterpolationType}),arg1,arg2)
end

function DMDASetElementType(arg1::DM,arg2::DMDAElementType)
    ccall((:DMDASetElementType,petsc),PetscErrorCode,(DM,DMDAElementType),arg1,arg2)
end

function DMDAGetElementType(arg1::DM,arg2::Ptr{DMDAElementType})
    ccall((:DMDAGetElementType,petsc),PetscErrorCode,(DM,Ptr{DMDAElementType}),arg1,arg2)
end

function DMDAGetElements(arg1::DM,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{Ptr{PetscInt}})
    ccall((:DMDAGetElements,petsc),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4)
end

function DMDARestoreElements(arg1::DM,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{Ptr{PetscInt}})
    ccall((:DMDARestoreElements,petsc),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4)
end

function DMDACreate(arg1::MPI_Comm,arg2::Ptr{DM})
    ccall((:DMDACreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{DM}),arg1,arg2)
end

function DMDASetSizes(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt)
    ccall((:DMDASetSizes,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end

function DMDACreate1d(arg1::MPI_Comm,arg2::DMBoundaryType,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::Ptr{PetscInt},arg7::Ptr{DM})
    ccall((:DMDACreate1d,petsc),PetscErrorCode,(MPI_Comm,DMBoundaryType,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{DM}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function DMDACreate2d(arg1::MPI_Comm,arg2::DMBoundaryType,arg3::DMBoundaryType,arg4::DMDAStencilType,arg5::PetscInt,arg6::PetscInt,arg7::PetscInt,arg8::PetscInt,arg9::PetscInt,arg10::PetscInt,arg11::Ptr{PetscInt},arg12::Ptr{PetscInt},arg13::Ptr{DM})
    ccall((:DMDACreate2d,petsc),PetscErrorCode,(MPI_Comm,DMBoundaryType,DMBoundaryType,DMDAStencilType,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{DM}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13)
end

function DMDACreate3d(arg1::MPI_Comm,arg2::DMBoundaryType,arg3::DMBoundaryType,arg4::DMBoundaryType,arg5::DMDAStencilType,arg6::PetscInt,arg7::PetscInt,arg8::PetscInt,arg9::PetscInt,arg10::PetscInt,arg11::PetscInt,arg12::PetscInt,arg13::PetscInt,arg14::Ptr{PetscInt},arg15::Ptr{PetscInt},arg16::Ptr{PetscInt},arg17::Ptr{DM})
    ccall((:DMDACreate3d,petsc),PetscErrorCode,(MPI_Comm,DMBoundaryType,DMBoundaryType,DMBoundaryType,DMDAStencilType,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{DM}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17)
end

function DMDAGlobalToNaturalBegin(arg1::DM,arg2::Vec,arg3::InsertMode,arg4::Vec)
    ccall((:DMDAGlobalToNaturalBegin,petsc),PetscErrorCode,(DM,Vec,InsertMode,Vec),arg1,arg2,arg3,arg4)
end

function DMDAGlobalToNaturalEnd(arg1::DM,arg2::Vec,arg3::InsertMode,arg4::Vec)
    ccall((:DMDAGlobalToNaturalEnd,petsc),PetscErrorCode,(DM,Vec,InsertMode,Vec),arg1,arg2,arg3,arg4)
end

function DMDANaturalToGlobalBegin(arg1::DM,arg2::Vec,arg3::InsertMode,arg4::Vec)
    ccall((:DMDANaturalToGlobalBegin,petsc),PetscErrorCode,(DM,Vec,InsertMode,Vec),arg1,arg2,arg3,arg4)
end

function DMDANaturalToGlobalEnd(arg1::DM,arg2::Vec,arg3::InsertMode,arg4::Vec)
    ccall((:DMDANaturalToGlobalEnd,petsc),PetscErrorCode,(DM,Vec,InsertMode,Vec),arg1,arg2,arg3,arg4)
end

function DMDAGetCorners(arg1::DM,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{PetscInt},arg6::Ptr{PetscInt},arg7::Ptr{PetscInt})
    ccall((:DMDAGetCorners,petsc),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function DMDAGetGhostCorners(arg1::DM,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{PetscInt},arg6::Ptr{PetscInt},arg7::Ptr{PetscInt})
    ccall((:DMDAGetGhostCorners,petsc),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function DMDAGetInfo(arg1::DM,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{PetscInt},arg6::Ptr{PetscInt},arg7::Ptr{PetscInt},arg8::Ptr{PetscInt},arg9::Ptr{PetscInt},arg10::Ptr{PetscInt},arg11::Ptr{DMBoundaryType},arg12::Ptr{DMBoundaryType},arg13::Ptr{DMBoundaryType},arg14::Ptr{DMDAStencilType})
    ccall((:DMDAGetInfo,petsc),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{DMBoundaryType},Ptr{DMBoundaryType},Ptr{DMBoundaryType},Ptr{DMDAStencilType}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14)
end

function DMDAGetProcessorSubset(arg1::DM,arg2::DMDADirection,arg3::PetscInt,arg4::Ptr{MPI_Comm})
    ccall((:DMDAGetProcessorSubset,petsc),PetscErrorCode,(DM,DMDADirection,PetscInt,Ptr{MPI_Comm}),arg1,arg2,arg3,arg4)
end

function DMDAGetProcessorSubsets(arg1::DM,arg2::DMDADirection,arg3::Ptr{MPI_Comm})
    ccall((:DMDAGetProcessorSubsets,petsc),PetscErrorCode,(DM,DMDADirection,Ptr{MPI_Comm}),arg1,arg2,arg3)
end

function DMDAGetRay(arg1::DM,arg2::DMDADirection,arg3::PetscInt,arg4::Ptr{Vec},arg5::Ptr{VecScatter})
    ccall((:DMDAGetRay,petsc),PetscErrorCode,(DM,DMDADirection,PetscInt,Ptr{Vec},Ptr{VecScatter}),arg1,arg2,arg3,arg4,arg5)
end

function DMDAGlobalToNaturalAllCreate(arg1::DM,arg2::Ptr{VecScatter})
    ccall((:DMDAGlobalToNaturalAllCreate,petsc),PetscErrorCode,(DM,Ptr{VecScatter}),arg1,arg2)
end

function DMDANaturalAllToGlobalCreate(arg1::DM,arg2::Ptr{VecScatter})
    ccall((:DMDANaturalAllToGlobalCreate,petsc),PetscErrorCode,(DM,Ptr{VecScatter}),arg1,arg2)
end

function DMDAGetScatter(arg1::DM,arg2::Ptr{VecScatter},arg3::Ptr{VecScatter})
    ccall((:DMDAGetScatter,petsc),PetscErrorCode,(DM,Ptr{VecScatter},Ptr{VecScatter}),arg1,arg2,arg3)
end

function DMDAGetNeighbors(arg1::DM,arg2::Ptr{Ptr{PetscMPIInt}})
    ccall((:DMDAGetNeighbors,petsc),PetscErrorCode,(DM,Ptr{Ptr{PetscMPIInt}}),arg1,arg2)
end

function DMDASetAOType(arg1::DM,arg2::AOType)
    ccall((:DMDASetAOType,petsc),PetscErrorCode,(DM,AOType),arg1,arg2)
end

function DMDAGetAO(arg1::DM,arg2::Ptr{Cint})
    ccall((:DMDAGetAO,petsc),PetscErrorCode,(DM,Ptr{Cint}),arg1,arg2)
end

function DMDASetUniformCoordinates(arg1::DM,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::Cint,arg5::Cint,arg6::Cint)
    ccall((:DMDASetUniformCoordinates,petsc),PetscErrorCode,(DM,Cint,Cint,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4,arg5,arg6)
end

function DMDAGetCoordinateArray(arg1::DM,arg2::Ptr{Void})
    ccall((:DMDAGetCoordinateArray,petsc),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end

function DMDARestoreCoordinateArray(arg1::DM,arg2::Ptr{Void})
    ccall((:DMDARestoreCoordinateArray,petsc),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end

function DMDAGetBoundingBox(arg1::DM,PetscReal::Ptr{Cint},arg2::Ptr{Cint})
    ccall((:DMDAGetBoundingBox,petsc),PetscErrorCode,(DM,Ptr{Cint},Ptr{Cint}),arg1,PetscReal,arg2)
end

function DMDAGetLocalBoundingBox(arg1::DM,PetscReal::Ptr{Cint},arg2::Ptr{Cint})
    ccall((:DMDAGetLocalBoundingBox,petsc),PetscErrorCode,(DM,Ptr{Cint},Ptr{Cint}),arg1,PetscReal,arg2)
end

function DMDAGetLogicalCoordinate(arg1::DM,arg2::PetscScalar,arg3::PetscScalar,arg4::PetscScalar,arg5::Ptr{PetscInt},arg6::Ptr{PetscInt},arg7::Ptr{PetscInt},arg8::Ptr{PetscScalar},arg9::Ptr{PetscScalar},arg10::Ptr{PetscScalar})
    ccall((:DMDAGetLogicalCoordinate,petsc),PetscErrorCode,(DM,PetscScalar,PetscScalar,PetscScalar,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscScalar},Ptr{PetscScalar},Ptr{PetscScalar}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function DMDAMapCoordsToPeriodicDomain(arg1::DM,arg2::Ptr{PetscScalar},arg3::Ptr{PetscScalar})
    ccall((:DMDAMapCoordsToPeriodicDomain,petsc),PetscErrorCode,(DM,Ptr{PetscScalar},Ptr{PetscScalar}),arg1,arg2,arg3)
end

function DMDAGetReducedDMDA(arg1::DM,arg2::PetscInt,arg3::Ptr{DM})
    ccall((:DMDAGetReducedDMDA,petsc),PetscErrorCode,(DM,PetscInt,Ptr{DM}),arg1,arg2,arg3)
end

function DMDASetFieldName(arg1::DM,arg2::PetscInt,arg3::Ptr{Uint8})
    ccall((:DMDASetFieldName,petsc),PetscErrorCode,(DM,PetscInt,Ptr{Uint8}),arg1,arg2,arg3)
end

function DMDAGetFieldName(arg1::DM,arg2::PetscInt,arg3::Ptr{Ptr{Uint8}})
    ccall((:DMDAGetFieldName,petsc),PetscErrorCode,(DM,PetscInt,Ptr{Ptr{Uint8}}),arg1,arg2,arg3)
end

function DMDASetFieldNames(arg1::DM,arg2::Ptr{Ptr{Uint8}})
    ccall((:DMDASetFieldNames,petsc),PetscErrorCode,(DM,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function DMDAGetFieldNames(arg1::DM,arg2::Ptr{Ptr{Ptr{Uint8}}})
    ccall((:DMDAGetFieldNames,petsc),PetscErrorCode,(DM,Ptr{Ptr{Ptr{Uint8}}}),arg1,arg2)
end

function DMDASetCoordinateName(arg1::DM,arg2::PetscInt,arg3::Ptr{Uint8})
    ccall((:DMDASetCoordinateName,petsc),PetscErrorCode,(DM,PetscInt,Ptr{Uint8}),arg1,arg2,arg3)
end

function DMDAGetCoordinateName(arg1::DM,arg2::PetscInt,arg3::Ptr{Ptr{Uint8}})
    ccall((:DMDAGetCoordinateName,petsc),PetscErrorCode,(DM,PetscInt,Ptr{Ptr{Uint8}}),arg1,arg2,arg3)
end

function DMDASetBoundaryType(arg1::DM,arg2::DMBoundaryType,arg3::DMBoundaryType,arg4::DMBoundaryType)
    ccall((:DMDASetBoundaryType,petsc),PetscErrorCode,(DM,DMBoundaryType,DMBoundaryType,DMBoundaryType),arg1,arg2,arg3,arg4)
end

function DMDASetDof(arg1::DM,arg2::PetscInt)
    ccall((:DMDASetDof,petsc),PetscErrorCode,(DM,PetscInt),arg1,arg2)
end

function DMDASetOverlap(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt)
    ccall((:DMDASetOverlap,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end

function DMDAGetOverlap(arg1::DM,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:DMDAGetOverlap,petsc),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function DMDASetNumLocalSubDomains(arg1::DM,arg2::PetscInt)
    ccall((:DMDASetNumLocalSubDomains,petsc),PetscErrorCode,(DM,PetscInt),arg1,arg2)
end

function DMDAGetNumLocalSubDomains(arg1::DM,arg2::Ptr{PetscInt})
    ccall((:DMDAGetNumLocalSubDomains,petsc),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end

function DMDAGetOffset(arg1::DM,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{PetscInt},arg6::Ptr{PetscInt},arg7::Ptr{PetscInt})
    ccall((:DMDAGetOffset,petsc),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function DMDASetOffset(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::PetscInt,arg7::PetscInt)
    ccall((:DMDASetOffset,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function DMDAGetNonOverlappingRegion(arg1::DM,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{PetscInt},arg6::Ptr{PetscInt},arg7::Ptr{PetscInt})
    ccall((:DMDAGetNonOverlappingRegion,petsc),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function DMDASetNonOverlappingRegion(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::PetscInt,arg7::PetscInt)
    ccall((:DMDASetNonOverlappingRegion,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function DMDASetStencilWidth(arg1::DM,arg2::PetscInt)
    ccall((:DMDASetStencilWidth,petsc),PetscErrorCode,(DM,PetscInt),arg1,arg2)
end

function DMDASetOwnershipRanges(arg1::DM,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:DMDASetOwnershipRanges,petsc),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function DMDAGetOwnershipRanges(arg1::DM,arg2::Ptr{Ptr{PetscInt}},arg3::Ptr{Ptr{PetscInt}},arg4::Ptr{Ptr{PetscInt}})
    ccall((:DMDAGetOwnershipRanges,petsc),PetscErrorCode,(DM,Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4)
end

function DMDASetNumProcs(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt)
    ccall((:DMDASetNumProcs,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end

function DMDASetStencilType(arg1::DM,arg2::DMDAStencilType)
    ccall((:DMDASetStencilType,petsc),PetscErrorCode,(DM,DMDAStencilType),arg1,arg2)
end

function DMDAVecGetArray(arg1::DM,arg2::Vec,arg3::Ptr{Void})
    ccall((:DMDAVecGetArray,petsc),PetscErrorCode,(DM,Vec,Ptr{Void}),arg1,arg2,arg3)
end

function DMDAVecRestoreArray(arg1::DM,arg2::Vec,arg3::Ptr{Void})
    ccall((:DMDAVecRestoreArray,petsc),PetscErrorCode,(DM,Vec,Ptr{Void}),arg1,arg2,arg3)
end

function DMDAVecGetArrayDOF(arg1::DM,arg2::Vec,arg3::Ptr{Void})
    ccall((:DMDAVecGetArrayDOF,petsc),PetscErrorCode,(DM,Vec,Ptr{Void}),arg1,arg2,arg3)
end

function DMDAVecRestoreArrayDOF(arg1::DM,arg2::Vec,arg3::Ptr{Void})
    ccall((:DMDAVecRestoreArrayDOF,petsc),PetscErrorCode,(DM,Vec,Ptr{Void}),arg1,arg2,arg3)
end

function DMDAVecGetArrayRead(arg1::DM,arg2::Vec,arg3::Ptr{Void})
    ccall((:DMDAVecGetArrayRead,petsc),PetscErrorCode,(DM,Vec,Ptr{Void}),arg1,arg2,arg3)
end

function DMDAVecRestoreArrayRead(arg1::DM,arg2::Vec,arg3::Ptr{Void})
    ccall((:DMDAVecRestoreArrayRead,petsc),PetscErrorCode,(DM,Vec,Ptr{Void}),arg1,arg2,arg3)
end

function DMDAVecGetArrayDOFRead(arg1::DM,arg2::Vec,arg3::Ptr{Void})
    ccall((:DMDAVecGetArrayDOFRead,petsc),PetscErrorCode,(DM,Vec,Ptr{Void}),arg1,arg2,arg3)
end

function DMDAVecRestoreArrayDOFRead(arg1::DM,arg2::Vec,arg3::Ptr{Void})
    ccall((:DMDAVecRestoreArrayDOFRead,petsc),PetscErrorCode,(DM,Vec,Ptr{Void}),arg1,arg2,arg3)
end

function DMDASplitComm2d(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::Ptr{MPI_Comm})
    ccall((:DMDASplitComm2d,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,Ptr{MPI_Comm}),arg1,arg2,arg3,arg4,arg5)
end

function DMDACreatePatchIS(arg1::DM,arg2::Ptr{MatStencil},arg3::Ptr{MatStencil},arg4::Ptr{IS})
    ccall((:DMDACreatePatchIS,petsc),PetscErrorCode,(DM,Ptr{MatStencil},Ptr{MatStencil},Ptr{IS}),arg1,arg2,arg3,arg4)
end

function DMDAGetLocalInfo(arg1::DM,arg2::Ptr{DMDALocalInfo})
    ccall((:DMDAGetLocalInfo,petsc),PetscErrorCode,(DM,Ptr{DMDALocalInfo}),arg1,arg2)
end

function MatRegisterDAAD()
    ccall((:MatRegisterDAAD,petsc),PetscErrorCode,())
end

function MatCreateDAAD(arg1::DM,arg2::Ptr{Mat})
    ccall((:MatCreateDAAD,petsc),PetscErrorCode,(DM,Ptr{Mat}),arg1,arg2)
end

function MatCreateSeqUSFFT(arg1::Vec,arg2::DM,arg3::Ptr{Mat})
    ccall((:MatCreateSeqUSFFT,petsc),PetscErrorCode,(Vec,DM,Ptr{Mat}),arg1,arg2,arg3)
end

function DMDASetGetMatrix(arg1::DM,arg2::Ptr{Void})
    ccall((:DMDASetGetMatrix,petsc),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end

function DMDASetBlockFills(arg1::DM,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt})
    ccall((:DMDASetBlockFills,petsc),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function DMDASetRefinementFactor(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt)
    ccall((:DMDASetRefinementFactor,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end

function DMDAGetRefinementFactor(arg1::DM,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:DMDAGetRefinementFactor,petsc),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function DMDAGetArray(arg1::DM,arg2::PetscBool,arg3::Ptr{Void})
    ccall((:DMDAGetArray,petsc),PetscErrorCode,(DM,PetscBool,Ptr{Void}),arg1,arg2,arg3)
end

function DMDARestoreArray(arg1::DM,arg2::PetscBool,arg3::Ptr{Void})
    ccall((:DMDARestoreArray,petsc),PetscErrorCode,(DM,PetscBool,Ptr{Void}),arg1,arg2,arg3)
end

function DMDACreatePF(arg1::DM,arg2::Ptr{PF})
    ccall((:DMDACreatePF,petsc),PetscErrorCode,(DM,Ptr{PF}),arg1,arg2)
end

function DMDAGetNumCells(arg1::DM,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{PetscInt})
    ccall((:DMDAGetNumCells,petsc),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end

function DMDAGetCellPoint(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::Ptr{PetscInt})
    ccall((:DMDAGetCellPoint,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end

function DMDAGetNumVertices(arg1::DM,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{PetscInt})
    ccall((:DMDAGetNumVertices,petsc),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end

function DMDAGetNumFaces(arg1::DM,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{PetscInt},arg6::Ptr{PetscInt},arg7::Ptr{PetscInt})
    ccall((:DMDAGetNumFaces,petsc),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function DMDAGetHeightStratum(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:DMDAGetHeightStratum,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function DMDAGetDepthStratum(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:DMDAGetDepthStratum,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function DMDACreateSection(arg1::DM,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{PetscSection})
    ccall((:DMDACreateSection,petsc),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscSection}),arg1,arg2,arg3,arg4,arg5)
end

function DMDAComputeCellGeometryFEM(arg1::DM,arg2::PetscInt,arg3::PetscQuadrature,PetscReal::Ptr{Cint},arg4::Ptr{Cint},arg5::Ptr{Cint},arg6::Ptr{Cint})
    ccall((:DMDAComputeCellGeometryFEM,petsc),PetscErrorCode,(DM,PetscInt,PetscQuadrature,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,PetscReal,arg4,arg5,arg6)
end

function DMDAGetTransitiveClosure(arg1::DM,arg2::PetscInt,arg3::PetscBool,arg4::Ptr{PetscInt},arg5::Ptr{Ptr{PetscInt}})
    ccall((:DMDAGetTransitiveClosure,petsc),PetscErrorCode,(DM,PetscInt,PetscBool,Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end

function DMDARestoreTransitiveClosure(arg1::DM,arg2::PetscInt,arg3::PetscBool,arg4::Ptr{PetscInt},arg5::Ptr{Ptr{PetscInt}})
    ccall((:DMDARestoreTransitiveClosure,petsc),PetscErrorCode,(DM,PetscInt,PetscBool,Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end

function DMDAVecGetClosure(arg1::DM,arg2::PetscSection,arg3::Vec,arg4::PetscInt,arg5::Ptr{PetscInt},arg6::Ptr{Ptr{PetscScalar}})
    ccall((:DMDAVecGetClosure,petsc),PetscErrorCode,(DM,PetscSection,Vec,PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscScalar}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function DMDAVecRestoreClosure(arg1::DM,arg2::PetscSection,arg3::Vec,arg4::PetscInt,arg5::Ptr{PetscInt},arg6::Ptr{Ptr{PetscScalar}})
    ccall((:DMDAVecRestoreClosure,petsc),PetscErrorCode,(DM,PetscSection,Vec,PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscScalar}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function DMDAVecSetClosure(arg1::DM,arg2::PetscSection,arg3::Vec,arg4::PetscInt,arg5::Ptr{PetscScalar},arg6::InsertMode)
    ccall((:DMDAVecSetClosure,petsc),PetscErrorCode,(DM,PetscSection,Vec,PetscInt,Ptr{PetscScalar},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6)
end

function DMDAGetClosure(arg1::DM,arg2::PetscSection,arg3::PetscInt,arg4::Ptr{PetscInt},arg5::Ptr{Ptr{PetscInt}})
    ccall((:DMDAGetClosure,petsc),PetscErrorCode,(DM,PetscSection,PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end

function DMDARestoreClosure(arg1::DM,arg2::PetscSection,arg3::PetscInt,arg4::Ptr{PetscInt},arg5::Ptr{Ptr{PetscInt}})
    ccall((:DMDARestoreClosure,petsc),PetscErrorCode,(DM,PetscSection,PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end

function DMDAGetClosureScalar(arg1::DM,arg2::PetscSection,arg3::PetscInt,arg4::Ptr{PetscScalar},arg5::Ptr{PetscInt},arg6::Ptr{Ptr{PetscScalar}})
    ccall((:DMDAGetClosureScalar,petsc),PetscErrorCode,(DM,PetscSection,PetscInt,Ptr{PetscScalar},Ptr{PetscInt},Ptr{Ptr{PetscScalar}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function DMDARestoreClosureScalar(arg1::DM,arg2::PetscSection,arg3::PetscInt,arg4::Ptr{PetscScalar},arg5::Ptr{PetscInt},arg6::Ptr{Ptr{PetscScalar}})
    ccall((:DMDARestoreClosureScalar,petsc),PetscErrorCode,(DM,PetscSection,PetscInt,Ptr{PetscScalar},Ptr{PetscInt},Ptr{Ptr{PetscScalar}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function DMDASetClosureScalar(arg1::DM,arg2::PetscSection,arg3::PetscInt,arg4::Ptr{PetscScalar},arg5::Ptr{PetscScalar},arg6::InsertMode)
    ccall((:DMDASetClosureScalar,petsc),PetscErrorCode,(DM,PetscSection,PetscInt,Ptr{PetscScalar},Ptr{PetscScalar},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6)
end

function DMDAConvertToCell(arg1::DM,arg2::MatStencil,arg3::Ptr{PetscInt})
    ccall((:DMDAConvertToCell,petsc),PetscErrorCode,(DM,MatStencil,Ptr{PetscInt}),arg1,arg2,arg3)
end

function DMDASetVertexCoordinates(arg1::DM,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::Cint,arg5::Cint,arg6::Cint)
    ccall((:DMDASetVertexCoordinates,petsc),PetscErrorCode,(DM,Cint,Cint,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4,arg5,arg6)
end

function DMDASetPreallocationCenterDimension(arg1::DM,arg2::PetscInt)
    ccall((:DMDASetPreallocationCenterDimension,petsc),PetscErrorCode,(DM,PetscInt),arg1,arg2)
end

function DMDAGetPreallocationCenterDimension(arg1::DM,arg2::Ptr{PetscInt})
    ccall((:DMDAGetPreallocationCenterDimension,petsc),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end

function DMDAProjectFunction(arg1::DM,arg2::Ptr{Ptr{Void}},arg3::Ptr{Ptr{Void}},arg4::InsertMode,arg5::Vec)
    ccall((:DMDAProjectFunction,petsc),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}},InsertMode,Vec),arg1,arg2,arg3,arg4,arg5)
end

function DMDAProjectFunctionLocal(arg1::DM,arg2::Ptr{Ptr{Void}},arg3::Ptr{Ptr{Void}},arg4::InsertMode,arg5::Vec)
    ccall((:DMDAProjectFunctionLocal,petsc),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}},InsertMode,Vec),arg1,arg2,arg3,arg4,arg5)
end

function DMDAComputeL2Diff(arg1::DM,arg2::Ptr{Ptr{Void}},arg3::Ptr{Ptr{Void}},arg4::Vec,arg5::Ptr{Cint})
    ccall((:DMDAComputeL2Diff,petsc),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}},Vec,Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
end

function DMDAComputeL2GradientDiff(arg1::DM,arg2::Ptr{Ptr{Void}},arg3::Ptr{Ptr{Void}},arg4::Vec,PetscReal::Ptr{Cint},arg5::Ptr{Cint})
    ccall((:DMDAComputeL2GradientDiff,petsc),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}},Vec,Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,PetscReal,arg5)
end

function DMCompositeCreate(arg1::MPI_Comm,arg2::Ptr{DM})
    ccall((:DMCompositeCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{DM}),arg1,arg2)
end

function DMCompositeAddDM(arg1::DM,arg2::DM)
    ccall((:DMCompositeAddDM,petsc),PetscErrorCode,(DM,DM),arg1,arg2)
end

function DMCompositeSetCoupling(arg1::DM,arg2::Ptr{Void})
    ccall((:DMCompositeSetCoupling,petsc),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end

function DMCompositeAddVecScatter(arg1::DM,arg2::VecScatter)
    ccall((:DMCompositeAddVecScatter,petsc),PetscErrorCode,(DM,VecScatter),arg1,arg2)
end

function DMCompositeScatterArray(arg1::DM,arg2::Vec,arg3::Ptr{Vec})
    ccall((:DMCompositeScatterArray,petsc),PetscErrorCode,(DM,Vec,Ptr{Vec}),arg1,arg2,arg3)
end

function DMCompositeGatherArray(arg1::DM,arg2::Vec,arg3::InsertMode,arg4::Ptr{Vec})
    ccall((:DMCompositeGatherArray,petsc),PetscErrorCode,(DM,Vec,InsertMode,Ptr{Vec}),arg1,arg2,arg3,arg4)
end

function DMCompositeGetNumberDM(arg1::DM,arg2::Ptr{PetscInt})
    ccall((:DMCompositeGetNumberDM,petsc),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end

function DMCompositeGetAccessArray(arg1::DM,arg2::Vec,arg3::PetscInt,arg4::Ptr{PetscInt},arg5::Ptr{Vec})
    ccall((:DMCompositeGetAccessArray,petsc),PetscErrorCode,(DM,Vec,PetscInt,Ptr{PetscInt},Ptr{Vec}),arg1,arg2,arg3,arg4,arg5)
end

function DMCompositeRestoreAccessArray(arg1::DM,arg2::Vec,arg3::PetscInt,arg4::Ptr{PetscInt},arg5::Ptr{Vec})
    ccall((:DMCompositeRestoreAccessArray,petsc),PetscErrorCode,(DM,Vec,PetscInt,Ptr{PetscInt},Ptr{Vec}),arg1,arg2,arg3,arg4,arg5)
end

function DMCompositeGetEntriesArray(arg1::DM,arg2::Ptr{DM})
    ccall((:DMCompositeGetEntriesArray,petsc),PetscErrorCode,(DM,Ptr{DM}),arg1,arg2)
end

function DMCompositeGetGlobalISs(arg1::DM,arg2::Ptr{Ptr{IS}})
    ccall((:DMCompositeGetGlobalISs,petsc),PetscErrorCode,(DM,Ptr{Ptr{IS}}),arg1,arg2)
end

function DMCompositeGetLocalISs(arg1::DM,arg2::Ptr{Ptr{IS}})
    ccall((:DMCompositeGetLocalISs,petsc),PetscErrorCode,(DM,Ptr{Ptr{IS}}),arg1,arg2)
end

function DMCompositeGetISLocalToGlobalMappings(arg1::DM,arg2::Ptr{Ptr{ISLocalToGlobalMapping}})
    ccall((:DMCompositeGetISLocalToGlobalMappings,petsc),PetscErrorCode,(DM,Ptr{Ptr{ISLocalToGlobalMapping}}),arg1,arg2)
end

function DMPatchCreate(arg1::MPI_Comm,arg2::Ptr{DM})
    ccall((:DMPatchCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{DM}),arg1,arg2)
end

function DMPatchZoom(arg1::DM,arg2::Vec,arg3::MatStencil,arg4::MatStencil,arg5::MPI_Comm,arg6::Ptr{DM},arg7::Ptr{PetscSF},arg8::Ptr{PetscSF})
    ccall((:DMPatchZoom,petsc),PetscErrorCode,(DM,Vec,MatStencil,MatStencil,MPI_Comm,Ptr{DM},Ptr{PetscSF},Ptr{PetscSF}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function DMPatchSolve(arg1::DM)
    ccall((:DMPatchSolve,petsc),PetscErrorCode,(DM,),arg1)
end

function DMPatchGetPatchSize(arg1::DM,arg2::Ptr{MatStencil})
    ccall((:DMPatchGetPatchSize,petsc),PetscErrorCode,(DM,Ptr{MatStencil}),arg1,arg2)
end

function DMPatchSetPatchSize(arg1::DM,arg2::MatStencil)
    ccall((:DMPatchSetPatchSize,petsc),PetscErrorCode,(DM,MatStencil),arg1,arg2)
end

function DMPatchGetCommSize(arg1::DM,arg2::Ptr{MatStencil})
    ccall((:DMPatchGetCommSize,petsc),PetscErrorCode,(DM,Ptr{MatStencil}),arg1,arg2)
end

function DMPatchSetCommSize(arg1::DM,arg2::MatStencil)
    ccall((:DMPatchSetCommSize,petsc),PetscErrorCode,(DM,MatStencil),arg1,arg2)
end

function DMPatchGetCoarse(arg1::DM,arg2::Ptr{DM})
    ccall((:DMPatchGetCoarse,petsc),PetscErrorCode,(DM,Ptr{DM}),arg1,arg2)
end

function DMPatchCreateGrid(arg1::MPI_Comm,arg2::PetscInt,arg3::MatStencil,arg4::MatStencil,arg5::MatStencil,arg6::Ptr{DM})
    ccall((:DMPatchCreateGrid,petsc),PetscErrorCode,(MPI_Comm,PetscInt,MatStencil,MatStencil,MatStencil,Ptr{DM}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscLimiterCreate(arg1::MPI_Comm,arg2::Ptr{PetscLimiter})
    ccall((:PetscLimiterCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscLimiter}),arg1,arg2)
end

function PetscLimiterDestroy(arg1::Ptr{PetscLimiter})
    ccall((:PetscLimiterDestroy,petsc),PetscErrorCode,(Ptr{PetscLimiter},),arg1)
end

function PetscLimiterSetType(arg1::PetscLimiter,arg2::PetscLimiterType)
    ccall((:PetscLimiterSetType,petsc),PetscErrorCode,(PetscLimiter,PetscLimiterType),arg1,arg2)
end

function PetscLimiterGetType(arg1::PetscLimiter,arg2::Ptr{PetscLimiterType})
    ccall((:PetscLimiterGetType,petsc),PetscErrorCode,(PetscLimiter,Ptr{PetscLimiterType}),arg1,arg2)
end

function PetscLimiterSetUp(arg1::PetscLimiter)
    ccall((:PetscLimiterSetUp,petsc),PetscErrorCode,(PetscLimiter,),arg1)
end

function PetscLimiterSetFromOptions(arg1::PetscLimiter)
    ccall((:PetscLimiterSetFromOptions,petsc),PetscErrorCode,(PetscLimiter,),arg1)
end

function PetscLimiterRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:PetscLimiterRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function PetscLimiterRegisterDestroy()
    ccall((:PetscLimiterRegisterDestroy,petsc),PetscErrorCode,())
end

function PetscLimiterLimit(arg1::PetscLimiter,PetscReal::Cint,arg2::Ptr{Cint})
    ccall((:PetscLimiterLimit,petsc),PetscErrorCode,(PetscLimiter,Cint,Ptr{Cint}),arg1,PetscReal,arg2)
end

function PetscFVInitializePackage()
    ccall((:PetscFVInitializePackage,petsc),PetscErrorCode,())
end

function PetscFVCreate(arg1::MPI_Comm,arg2::Ptr{PetscFV})
    ccall((:PetscFVCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscFV}),arg1,arg2)
end

function PetscFVDestroy(arg1::Ptr{PetscFV})
    ccall((:PetscFVDestroy,petsc),PetscErrorCode,(Ptr{PetscFV},),arg1)
end

function PetscFVSetType(arg1::PetscFV,arg2::PetscFVType)
    ccall((:PetscFVSetType,petsc),PetscErrorCode,(PetscFV,PetscFVType),arg1,arg2)
end

function PetscFVGetType(arg1::PetscFV,arg2::Ptr{PetscFVType})
    ccall((:PetscFVGetType,petsc),PetscErrorCode,(PetscFV,Ptr{PetscFVType}),arg1,arg2)
end

function PetscFVSetUp(arg1::PetscFV)
    ccall((:PetscFVSetUp,petsc),PetscErrorCode,(PetscFV,),arg1)
end

function PetscFVSetFromOptions(arg1::PetscFV)
    ccall((:PetscFVSetFromOptions,petsc),PetscErrorCode,(PetscFV,),arg1)
end

function PetscFVRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:PetscFVRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function PetscFVRegisterDestroy()
    ccall((:PetscFVRegisterDestroy,petsc),PetscErrorCode,())
end

function PetscFVSetLimiter(arg1::PetscFV,arg2::PetscLimiter)
    ccall((:PetscFVSetLimiter,petsc),PetscErrorCode,(PetscFV,PetscLimiter),arg1,arg2)
end

function PetscFVGetLimiter(arg1::PetscFV,arg2::Ptr{PetscLimiter})
    ccall((:PetscFVGetLimiter,petsc),PetscErrorCode,(PetscFV,Ptr{PetscLimiter}),arg1,arg2)
end

function PetscFVSetNumComponents(arg1::PetscFV,arg2::PetscInt)
    ccall((:PetscFVSetNumComponents,petsc),PetscErrorCode,(PetscFV,PetscInt),arg1,arg2)
end

function PetscFVGetNumComponents(arg1::PetscFV,arg2::Ptr{PetscInt})
    ccall((:PetscFVGetNumComponents,petsc),PetscErrorCode,(PetscFV,Ptr{PetscInt}),arg1,arg2)
end

function PetscFVSetSpatialDimension(arg1::PetscFV,arg2::PetscInt)
    ccall((:PetscFVSetSpatialDimension,petsc),PetscErrorCode,(PetscFV,PetscInt),arg1,arg2)
end

function PetscFVGetSpatialDimension(arg1::PetscFV,arg2::Ptr{PetscInt})
    ccall((:PetscFVGetSpatialDimension,petsc),PetscErrorCode,(PetscFV,Ptr{PetscInt}),arg1,arg2)
end

function PetscFVSetComputeGradients(arg1::PetscFV,arg2::PetscBool)
    ccall((:PetscFVSetComputeGradients,petsc),PetscErrorCode,(PetscFV,PetscBool),arg1,arg2)
end

function PetscFVGetComputeGradients(arg1::PetscFV,arg2::Ptr{PetscBool})
    ccall((:PetscFVGetComputeGradients,petsc),PetscErrorCode,(PetscFV,Ptr{PetscBool}),arg1,arg2)
end

function PetscFVSetQuadrature(arg1::PetscFV,arg2::PetscQuadrature)
    ccall((:PetscFVSetQuadrature,petsc),PetscErrorCode,(PetscFV,PetscQuadrature),arg1,arg2)
end

function PetscFVGetQuadrature(arg1::PetscFV,arg2::Ptr{PetscQuadrature})
    ccall((:PetscFVGetQuadrature,petsc),PetscErrorCode,(PetscFV,Ptr{PetscQuadrature}),arg1,arg2)
end

function PetscFVSetDualSpace(arg1::PetscFV,arg2::PetscDualSpace)
    ccall((:PetscFVSetDualSpace,petsc),PetscErrorCode,(PetscFV,PetscDualSpace),arg1,arg2)
end

function PetscFVGetDualSpace(arg1::PetscFV,arg2::Ptr{PetscDualSpace})
    ccall((:PetscFVGetDualSpace,petsc),PetscErrorCode,(PetscFV,Ptr{PetscDualSpace}),arg1,arg2)
end

function PetscFVRefine(arg1::PetscFV,arg2::Ptr{PetscFV})
    ccall((:PetscFVRefine,petsc),PetscErrorCode,(PetscFV,Ptr{PetscFV}),arg1,arg2)
end

function PetscFVGetDefaultTabulation(arg1::PetscFV,arg2::Ptr{Ptr{Cint}},arg3::Ptr{Ptr{Cint}},arg4::Ptr{Ptr{Cint}})
    ccall((:PetscFVGetDefaultTabulation,petsc),PetscErrorCode,(PetscFV,Ptr{Ptr{Cint}},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}}),arg1,arg2,arg3,arg4)
end

function PetscFVGetTabulation(arg1::PetscFV,arg2::PetscInt,PetscReal::Ptr{Cint},arg3::Ptr{Ptr{Cint}},arg4::Ptr{Ptr{Cint}},arg5::Ptr{Ptr{Cint}})
    ccall((:PetscFVGetTabulation,petsc),PetscErrorCode,(PetscFV,PetscInt,Ptr{Cint},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}}),arg1,arg2,PetscReal,arg3,arg4,arg5)
end

function PetscFVRestoreTabulation(arg1::PetscFV,arg2::PetscInt,PetscReal::Ptr{Cint},arg3::Ptr{Ptr{Cint}},arg4::Ptr{Ptr{Cint}},arg5::Ptr{Ptr{Cint}})
    ccall((:PetscFVRestoreTabulation,petsc),PetscErrorCode,(PetscFV,PetscInt,Ptr{Cint},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}}),arg1,arg2,PetscReal,arg3,arg4,arg5)
end

function PetscFVComputeGradient(arg1::PetscFV,arg2::PetscInt,arg3::Ptr{PetscScalar},arg4::Ptr{PetscScalar})
    ccall((:PetscFVComputeGradient,petsc),PetscErrorCode,(PetscFV,PetscInt,Ptr{PetscScalar},Ptr{PetscScalar}),arg1,arg2,arg3,arg4)
end

function PetscFVIntegrateRHSFunction(arg1::PetscFV,arg2::PetscDS,arg3::PetscInt,arg4::PetscInt,arg5::Ptr{PetscFVFaceGeom},arg6::Ptr{Cint},arg7::Ptr{PetscScalar},arg8::Ptr{PetscScalar},arg9::Ptr{PetscScalar},arg10::Ptr{PetscScalar})
    ccall((:PetscFVIntegrateRHSFunction,petsc),PetscErrorCode,(PetscFV,PetscDS,PetscInt,PetscInt,Ptr{PetscFVFaceGeom},Ptr{Cint},Ptr{PetscScalar},Ptr{PetscScalar},Ptr{PetscScalar},Ptr{PetscScalar}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function PetscFVLeastSquaresSetMaxFaces(arg1::PetscFV,arg2::PetscInt)
    ccall((:PetscFVLeastSquaresSetMaxFaces,petsc),PetscErrorCode,(PetscFV,PetscInt),arg1,arg2)
end

function PetscPartitionerCreate(arg1::MPI_Comm,arg2::Ptr{PetscPartitioner})
    ccall((:PetscPartitionerCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscPartitioner}),arg1,arg2)
end

function PetscPartitionerDestroy(arg1::Ptr{PetscPartitioner})
    ccall((:PetscPartitionerDestroy,petsc),PetscErrorCode,(Ptr{PetscPartitioner},),arg1)
end

function PetscPartitionerSetType(arg1::PetscPartitioner,arg2::PetscPartitionerType)
    ccall((:PetscPartitionerSetType,petsc),PetscErrorCode,(PetscPartitioner,PetscPartitionerType),arg1,arg2)
end

function PetscPartitionerGetType(arg1::PetscPartitioner,arg2::Ptr{PetscPartitionerType})
    ccall((:PetscPartitionerGetType,petsc),PetscErrorCode,(PetscPartitioner,Ptr{PetscPartitionerType}),arg1,arg2)
end

function PetscPartitionerSetUp(arg1::PetscPartitioner)
    ccall((:PetscPartitionerSetUp,petsc),PetscErrorCode,(PetscPartitioner,),arg1)
end

function PetscPartitionerSetFromOptions(arg1::PetscPartitioner)
    ccall((:PetscPartitionerSetFromOptions,petsc),PetscErrorCode,(PetscPartitioner,),arg1)
end

function PetscPartitionerRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:PetscPartitionerRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function PetscPartitionerRegisterDestroy()
    ccall((:PetscPartitionerRegisterDestroy,petsc),PetscErrorCode,())
end

function PetscPartitionerPartition(arg1::PetscPartitioner,arg2::DM,arg3::PetscSection,arg4::Ptr{IS})
    ccall((:PetscPartitionerPartition,petsc),PetscErrorCode,(PetscPartitioner,DM,PetscSection,Ptr{IS}),arg1,arg2,arg3,arg4)
end

function PetscPartitionerShellSetPartition(arg1::PetscPartitioner,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:PetscPartitionerShellSetPartition,petsc),PetscErrorCode,(PetscPartitioner,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function DMPlexCreate(arg1::MPI_Comm,arg2::Ptr{DM})
    ccall((:DMPlexCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{DM}),arg1,arg2)
end

function DMPlexCreateCohesiveSubmesh(arg1::DM,arg2::PetscBool,arg3::Ptr{Uint8},arg4::PetscInt,arg5::Ptr{DM})
    ccall((:DMPlexCreateCohesiveSubmesh,petsc),PetscErrorCode,(DM,PetscBool,Ptr{Uint8},PetscInt,Ptr{DM}),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexCreateFromCellList(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::PetscBool,arg7::Ptr{Cint},arg8::PetscInt,arg9::Ptr{Cdouble},arg10::Ptr{DM})
    ccall((:DMPlexCreateFromCellList,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscBool,Ptr{Cint},PetscInt,Ptr{Cdouble},Ptr{DM}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function DMPlexCreateFromDAG(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{PetscInt},arg6::Ptr{PetscInt},arg7::Ptr{PetscScalar})
    ccall((:DMPlexCreateFromDAG,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscScalar}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function DMPlexCreateReferenceCell(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscBool,arg4::Ptr{DM})
    ccall((:DMPlexCreateReferenceCell,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscBool,Ptr{DM}),arg1,arg2,arg3,arg4)
end

function DMPlexGetChart(arg1::DM,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt})
    ccall((:DMPlexGetChart,petsc),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function DMPlexSetChart(arg1::DM,arg2::PetscInt,arg3::PetscInt)
    ccall((:DMPlexSetChart,petsc),PetscErrorCode,(DM,PetscInt,PetscInt),arg1,arg2,arg3)
end

function DMPlexGetConeSize(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt})
    ccall((:DMPlexGetConeSize,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function DMPlexSetConeSize(arg1::DM,arg2::PetscInt,arg3::PetscInt)
    ccall((:DMPlexSetConeSize,petsc),PetscErrorCode,(DM,PetscInt,PetscInt),arg1,arg2,arg3)
end

function DMPlexAddConeSize(arg1::DM,arg2::PetscInt,arg3::PetscInt)
    ccall((:DMPlexAddConeSize,petsc),PetscErrorCode,(DM,PetscInt,PetscInt),arg1,arg2,arg3)
end

function DMPlexGetCone(arg1::DM,arg2::PetscInt,arg3::Ptr{Ptr{PetscInt}})
    ccall((:DMPlexGetCone,petsc),PetscErrorCode,(DM,PetscInt,Ptr{Ptr{PetscInt}}),arg1,arg2,arg3)
end

function DMPlexSetCone(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt})
    ccall((:DMPlexSetCone,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function DMPlexInsertCone(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt)
    ccall((:DMPlexInsertCone,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end

function DMPlexInsertConeOrientation(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt)
    ccall((:DMPlexInsertConeOrientation,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end

function DMPlexGetConeOrientation(arg1::DM,arg2::PetscInt,arg3::Ptr{Ptr{PetscInt}})
    ccall((:DMPlexGetConeOrientation,petsc),PetscErrorCode,(DM,PetscInt,Ptr{Ptr{PetscInt}}),arg1,arg2,arg3)
end

function DMPlexSetConeOrientation(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt})
    ccall((:DMPlexSetConeOrientation,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function DMPlexGetSupportSize(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt})
    ccall((:DMPlexGetSupportSize,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function DMPlexSetSupportSize(arg1::DM,arg2::PetscInt,arg3::PetscInt)
    ccall((:DMPlexSetSupportSize,petsc),PetscErrorCode,(DM,PetscInt,PetscInt),arg1,arg2,arg3)
end

function DMPlexGetSupport(arg1::DM,arg2::PetscInt,arg3::Ptr{Ptr{PetscInt}})
    ccall((:DMPlexGetSupport,petsc),PetscErrorCode,(DM,PetscInt,Ptr{Ptr{PetscInt}}),arg1,arg2,arg3)
end

function DMPlexSetSupport(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt})
    ccall((:DMPlexSetSupport,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function DMPlexInsertSupport(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt)
    ccall((:DMPlexInsertSupport,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end

function DMPlexGetConeSection(arg1::DM,arg2::Ptr{PetscSection})
    ccall((:DMPlexGetConeSection,petsc),PetscErrorCode,(DM,Ptr{PetscSection}),arg1,arg2)
end

function DMPlexGetSupportSection(arg1::DM,arg2::Ptr{PetscSection})
    ccall((:DMPlexGetSupportSection,petsc),PetscErrorCode,(DM,Ptr{PetscSection}),arg1,arg2)
end

function DMPlexGetCones(arg1::DM,arg2::Ptr{Ptr{PetscInt}})
    ccall((:DMPlexGetCones,petsc),PetscErrorCode,(DM,Ptr{Ptr{PetscInt}}),arg1,arg2)
end

function DMPlexGetConeOrientations(arg1::DM,arg2::Ptr{Ptr{PetscInt}})
    ccall((:DMPlexGetConeOrientations,petsc),PetscErrorCode,(DM,Ptr{Ptr{PetscInt}}),arg1,arg2)
end

function DMPlexGetMaxSizes(arg1::DM,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt})
    ccall((:DMPlexGetMaxSizes,petsc),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function DMPlexSymmetrize(arg1::DM)
    ccall((:DMPlexSymmetrize,petsc),PetscErrorCode,(DM,),arg1)
end

function DMPlexStratify(arg1::DM)
    ccall((:DMPlexStratify,petsc),PetscErrorCode,(DM,),arg1)
end

function DMPlexEqual(arg1::DM,arg2::DM,arg3::Ptr{PetscBool})
    ccall((:DMPlexEqual,petsc),PetscErrorCode,(DM,DM,Ptr{PetscBool}),arg1,arg2,arg3)
end

function DMPlexReverseCell(arg1::DM,arg2::PetscInt)
    ccall((:DMPlexReverseCell,petsc),PetscErrorCode,(DM,PetscInt),arg1,arg2)
end

function DMPlexOrient(arg1::DM)
    ccall((:DMPlexOrient,petsc),PetscErrorCode,(DM,),arg1)
end

function DMPlexInterpolate(arg1::DM,arg2::Ptr{DM})
    ccall((:DMPlexInterpolate,petsc),PetscErrorCode,(DM,Ptr{DM}),arg1,arg2)
end

function DMPlexUninterpolate(arg1::DM,arg2::Ptr{DM})
    ccall((:DMPlexUninterpolate,petsc),PetscErrorCode,(DM,Ptr{DM}),arg1,arg2)
end

function DMPlexLoad(arg1::PetscViewer,arg2::DM)
    ccall((:DMPlexLoad,petsc),PetscErrorCode,(PetscViewer,DM),arg1,arg2)
end

function DMPlexPreallocateOperator(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{PetscInt},arg6::Ptr{PetscInt},arg7::Mat,arg8::PetscBool)
    ccall((:DMPlexPreallocateOperator,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Mat,PetscBool),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function DMPlexGetPointLocal(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:DMPlexGetPointLocal,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function DMPlexPointLocalRead(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscScalar},arg4::Ptr{Void})
    ccall((:DMPlexPointLocalRead,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscScalar},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function DMPlexPointLocalRef(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscScalar},arg4::Ptr{Void})
    ccall((:DMPlexPointLocalRef,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscScalar},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function DMPlexGetPointLocalField(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscInt},arg5::Ptr{PetscInt})
    ccall((:DMPlexGetPointLocalField,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexPointLocalFieldRef(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscScalar},arg5::Ptr{Void})
    ccall((:DMPlexPointLocalFieldRef,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,Ptr{PetscScalar},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexPointLocalFieldRead(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscScalar},arg5::Ptr{Void})
    ccall((:DMPlexPointLocalFieldRead,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,Ptr{PetscScalar},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexGetPointGlobal(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:DMPlexGetPointGlobal,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function DMPlexPointGlobalRead(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscScalar},arg4::Ptr{Void})
    ccall((:DMPlexPointGlobalRead,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscScalar},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function DMPlexPointGlobalRef(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscScalar},arg4::Ptr{Void})
    ccall((:DMPlexPointGlobalRef,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscScalar},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function DMPlexGetPointGlobalField(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscInt},arg5::Ptr{PetscInt})
    ccall((:DMPlexGetPointGlobalField,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexPointGlobalFieldRef(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscScalar},arg5::Ptr{Void})
    ccall((:DMPlexPointGlobalFieldRef,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,Ptr{PetscScalar},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexPointGlobalFieldRead(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscScalar},arg5::Ptr{Void})
    ccall((:DMPlexPointGlobalFieldRead,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,Ptr{PetscScalar},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end

function DMLabelCreate(arg1::Ptr{Uint8},arg2::Ptr{DMLabel})
    ccall((:DMLabelCreate,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{DMLabel}),arg1,arg2)
end

function DMLabelView(arg1::DMLabel,arg2::PetscViewer)
    ccall((:DMLabelView,petsc),PetscErrorCode,(DMLabel,PetscViewer),arg1,arg2)
end

function DMLabelDestroy(arg1::Ptr{DMLabel})
    ccall((:DMLabelDestroy,petsc),PetscErrorCode,(Ptr{DMLabel},),arg1)
end

function DMLabelDuplicate(arg1::DMLabel,arg2::Ptr{DMLabel})
    ccall((:DMLabelDuplicate,petsc),PetscErrorCode,(DMLabel,Ptr{DMLabel}),arg1,arg2)
end

function DMLabelGetName(arg1::DMLabel,arg2::Ptr{Ptr{Uint8}})
    ccall((:DMLabelGetName,petsc),PetscErrorCode,(DMLabel,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function DMLabelGetValue(arg1::DMLabel,arg2::PetscInt,arg3::Ptr{PetscInt})
    ccall((:DMLabelGetValue,petsc),PetscErrorCode,(DMLabel,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function DMLabelSetValue(arg1::DMLabel,arg2::PetscInt,arg3::PetscInt)
    ccall((:DMLabelSetValue,petsc),PetscErrorCode,(DMLabel,PetscInt,PetscInt),arg1,arg2,arg3)
end

function DMLabelClearValue(arg1::DMLabel,arg2::PetscInt,arg3::PetscInt)
    ccall((:DMLabelClearValue,petsc),PetscErrorCode,(DMLabel,PetscInt,PetscInt),arg1,arg2,arg3)
end

function DMLabelInsertIS(arg1::DMLabel,arg2::IS,arg3::PetscInt)
    ccall((:DMLabelInsertIS,petsc),PetscErrorCode,(DMLabel,IS,PetscInt),arg1,arg2,arg3)
end

function DMLabelGetNumValues(arg1::DMLabel,arg2::Ptr{PetscInt})
    ccall((:DMLabelGetNumValues,petsc),PetscErrorCode,(DMLabel,Ptr{PetscInt}),arg1,arg2)
end

function DMLabelGetStratumBounds(arg1::DMLabel,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:DMLabelGetStratumBounds,petsc),PetscErrorCode,(DMLabel,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function DMLabelGetValueIS(arg1::DMLabel,arg2::Ptr{IS})
    ccall((:DMLabelGetValueIS,petsc),PetscErrorCode,(DMLabel,Ptr{IS}),arg1,arg2)
end

function DMLabelStratumHasPoint(arg1::DMLabel,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscBool})
    ccall((:DMLabelStratumHasPoint,petsc),PetscErrorCode,(DMLabel,PetscInt,PetscInt,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function DMLabelGetStratumSize(arg1::DMLabel,arg2::PetscInt,arg3::Ptr{PetscInt})
    ccall((:DMLabelGetStratumSize,petsc),PetscErrorCode,(DMLabel,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function DMLabelGetStratumIS(arg1::DMLabel,arg2::PetscInt,arg3::Ptr{IS})
    ccall((:DMLabelGetStratumIS,petsc),PetscErrorCode,(DMLabel,PetscInt,Ptr{IS}),arg1,arg2,arg3)
end

function DMLabelClearStratum(arg1::DMLabel,arg2::PetscInt)
    ccall((:DMLabelClearStratum,petsc),PetscErrorCode,(DMLabel,PetscInt),arg1,arg2)
end

function DMLabelCreateIndex(arg1::DMLabel,arg2::PetscInt,arg3::PetscInt)
    ccall((:DMLabelCreateIndex,petsc),PetscErrorCode,(DMLabel,PetscInt,PetscInt),arg1,arg2,arg3)
end

function DMLabelDestroyIndex(arg1::DMLabel)
    ccall((:DMLabelDestroyIndex,petsc),PetscErrorCode,(DMLabel,),arg1)
end

function DMLabelHasValue(arg1::DMLabel,arg2::PetscInt,arg3::Ptr{PetscBool})
    ccall((:DMLabelHasValue,petsc),PetscErrorCode,(DMLabel,PetscInt,Ptr{PetscBool}),arg1,arg2,arg3)
end

function DMLabelHasPoint(arg1::DMLabel,arg2::PetscInt,arg3::Ptr{PetscBool})
    ccall((:DMLabelHasPoint,petsc),PetscErrorCode,(DMLabel,PetscInt,Ptr{PetscBool}),arg1,arg2,arg3)
end

function DMLabelFilter(arg1::DMLabel,arg2::PetscInt,arg3::PetscInt)
    ccall((:DMLabelFilter,petsc),PetscErrorCode,(DMLabel,PetscInt,PetscInt),arg1,arg2,arg3)
end

function DMLabelPermute(arg1::DMLabel,arg2::IS,arg3::Ptr{DMLabel})
    ccall((:DMLabelPermute,petsc),PetscErrorCode,(DMLabel,IS,Ptr{DMLabel}),arg1,arg2,arg3)
end

function DMLabelDistribute(arg1::DMLabel,arg2::PetscSF,arg3::Ptr{DMLabel})
    ccall((:DMLabelDistribute,petsc),PetscErrorCode,(DMLabel,PetscSF,Ptr{DMLabel}),arg1,arg2,arg3)
end

function DMPlexCreateLabel(arg1::DM,arg2::Ptr{Uint8})
    ccall((:DMPlexCreateLabel,petsc),PetscErrorCode,(DM,Ptr{Uint8}),arg1,arg2)
end

function DMPlexGetLabelValue(arg1::DM,arg2::Ptr{Uint8},arg3::PetscInt,arg4::Ptr{PetscInt})
    ccall((:DMPlexGetLabelValue,petsc),PetscErrorCode,(DM,Ptr{Uint8},PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function DMPlexSetLabelValue(arg1::DM,arg2::Ptr{Uint8},arg3::PetscInt,arg4::PetscInt)
    ccall((:DMPlexSetLabelValue,petsc),PetscErrorCode,(DM,Ptr{Uint8},PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end

function DMPlexClearLabelValue(arg1::DM,arg2::Ptr{Uint8},arg3::PetscInt,arg4::PetscInt)
    ccall((:DMPlexClearLabelValue,petsc),PetscErrorCode,(DM,Ptr{Uint8},PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end

function DMPlexGetLabelSize(arg1::DM,arg2::Ptr{Uint8},arg3::Ptr{PetscInt})
    ccall((:DMPlexGetLabelSize,petsc),PetscErrorCode,(DM,Ptr{Uint8},Ptr{PetscInt}),arg1,arg2,arg3)
end

function DMPlexGetLabelIdIS(arg1::DM,arg2::Ptr{Uint8},arg3::Ptr{IS})
    ccall((:DMPlexGetLabelIdIS,petsc),PetscErrorCode,(DM,Ptr{Uint8},Ptr{IS}),arg1,arg2,arg3)
end

function DMPlexGetStratumSize(arg1::DM,arg2::Ptr{Uint8},arg3::PetscInt,arg4::Ptr{PetscInt})
    ccall((:DMPlexGetStratumSize,petsc),PetscErrorCode,(DM,Ptr{Uint8},PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function DMPlexGetStratumIS(arg1::DM,arg2::Ptr{Uint8},arg3::PetscInt,arg4::Ptr{IS})
    ccall((:DMPlexGetStratumIS,petsc),PetscErrorCode,(DM,Ptr{Uint8},PetscInt,Ptr{IS}),arg1,arg2,arg3,arg4)
end

function DMPlexClearLabelStratum(arg1::DM,arg2::Ptr{Uint8},arg3::PetscInt)
    ccall((:DMPlexClearLabelStratum,petsc),PetscErrorCode,(DM,Ptr{Uint8},PetscInt),arg1,arg2,arg3)
end

function DMPlexGetLabelOutput(arg1::DM,arg2::Ptr{Uint8},arg3::Ptr{PetscBool})
    ccall((:DMPlexGetLabelOutput,petsc),PetscErrorCode,(DM,Ptr{Uint8},Ptr{PetscBool}),arg1,arg2,arg3)
end

function DMPlexSetLabelOutput(arg1::DM,arg2::Ptr{Uint8},arg3::PetscBool)
    ccall((:DMPlexSetLabelOutput,petsc),PetscErrorCode,(DM,Ptr{Uint8},PetscBool),arg1,arg2,arg3)
end

function PetscSectionCreateGlobalSectionLabel(arg1::PetscSection,arg2::PetscSF,arg3::PetscBool,arg4::DMLabel,arg5::PetscInt,arg6::Ptr{PetscSection})
    ccall((:PetscSectionCreateGlobalSectionLabel,petsc),PetscErrorCode,(PetscSection,PetscSF,PetscBool,DMLabel,PetscInt,Ptr{PetscSection}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function DMPlexGetNumLabels(arg1::DM,arg2::Ptr{PetscInt})
    ccall((:DMPlexGetNumLabels,petsc),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end

function DMPlexGetLabelName(arg1::DM,arg2::PetscInt,arg3::Ptr{Ptr{Uint8}})
    ccall((:DMPlexGetLabelName,petsc),PetscErrorCode,(DM,PetscInt,Ptr{Ptr{Uint8}}),arg1,arg2,arg3)
end

function DMPlexHasLabel(arg1::DM,arg2::Ptr{Uint8},arg3::Ptr{PetscBool})
    ccall((:DMPlexHasLabel,petsc),PetscErrorCode,(DM,Ptr{Uint8},Ptr{PetscBool}),arg1,arg2,arg3)
end

function DMPlexGetLabel(arg1::DM,arg2::Ptr{Uint8},arg3::Ptr{DMLabel})
    ccall((:DMPlexGetLabel,petsc),PetscErrorCode,(DM,Ptr{Uint8},Ptr{DMLabel}),arg1,arg2,arg3)
end

function DMPlexGetLabelByNum(arg1::DM,arg2::PetscInt,arg3::Ptr{DMLabel})
    ccall((:DMPlexGetLabelByNum,petsc),PetscErrorCode,(DM,PetscInt,Ptr{DMLabel}),arg1,arg2,arg3)
end

function DMPlexAddLabel(arg1::DM,arg2::DMLabel)
    ccall((:DMPlexAddLabel,petsc),PetscErrorCode,(DM,DMLabel),arg1,arg2)
end

function DMPlexRemoveLabel(arg1::DM,arg2::Ptr{Uint8},arg3::Ptr{DMLabel})
    ccall((:DMPlexRemoveLabel,petsc),PetscErrorCode,(DM,Ptr{Uint8},Ptr{DMLabel}),arg1,arg2,arg3)
end

function DMPlexGetCellNumbering(arg1::DM,arg2::Ptr{IS})
    ccall((:DMPlexGetCellNumbering,petsc),PetscErrorCode,(DM,Ptr{IS}),arg1,arg2)
end

function DMPlexGetVertexNumbering(arg1::DM,arg2::Ptr{IS})
    ccall((:DMPlexGetVertexNumbering,petsc),PetscErrorCode,(DM,Ptr{IS}),arg1,arg2)
end

function DMPlexCreatePointNumbering(arg1::DM,arg2::Ptr{IS})
    ccall((:DMPlexCreatePointNumbering,petsc),PetscErrorCode,(DM,Ptr{IS}),arg1,arg2)
end

function DMPlexGetDepth(arg1::DM,arg2::Ptr{PetscInt})
    ccall((:DMPlexGetDepth,petsc),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end

function DMPlexGetDepthLabel(arg1::DM,arg2::Ptr{DMLabel})
    ccall((:DMPlexGetDepthLabel,petsc),PetscErrorCode,(DM,Ptr{DMLabel}),arg1,arg2)
end

function DMPlexGetDepthStratum(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:DMPlexGetDepthStratum,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function DMPlexGetHeightStratum(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:DMPlexGetHeightStratum,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function DMPlexGetMeet(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{Ptr{PetscInt}})
    ccall((:DMPlexGetMeet,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexGetFullMeet(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{Ptr{PetscInt}})
    ccall((:DMPlexGetFullMeet,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexRestoreMeet(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{Ptr{PetscInt}})
    ccall((:DMPlexRestoreMeet,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexGetJoin(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{Ptr{PetscInt}})
    ccall((:DMPlexGetJoin,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexGetFullJoin(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{Ptr{PetscInt}})
    ccall((:DMPlexGetFullJoin,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexRestoreJoin(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{Ptr{PetscInt}})
    ccall((:DMPlexRestoreJoin,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexGetTransitiveClosure(arg1::DM,arg2::PetscInt,arg3::PetscBool,arg4::Ptr{PetscInt},arg5::Ptr{Ptr{PetscInt}})
    ccall((:DMPlexGetTransitiveClosure,petsc),PetscErrorCode,(DM,PetscInt,PetscBool,Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexRestoreTransitiveClosure(arg1::DM,arg2::PetscInt,arg3::PetscBool,arg4::Ptr{PetscInt},arg5::Ptr{Ptr{PetscInt}})
    ccall((:DMPlexRestoreTransitiveClosure,petsc),PetscErrorCode,(DM,PetscInt,PetscBool,Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexGenerate(arg1::DM,arg2::Ptr{Uint8},arg3::PetscBool,arg4::Ptr{DM})
    ccall((:DMPlexGenerate,petsc),PetscErrorCode,(DM,Ptr{Uint8},PetscBool,Ptr{DM}),arg1,arg2,arg3,arg4)
end

function DMPlexCopyCoordinates(arg1::DM,arg2::DM)
    ccall((:DMPlexCopyCoordinates,petsc),PetscErrorCode,(DM,DM),arg1,arg2)
end

function DMPlexCopyLabels(arg1::DM,arg2::DM)
    ccall((:DMPlexCopyLabels,petsc),PetscErrorCode,(DM,DM),arg1,arg2)
end

function DMPlexCreateDoublet(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscBool,arg4::PetscBool,arg5::PetscBool,PetscReal::Cint,arg6::Ptr{DM})
    ccall((:DMPlexCreateDoublet,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscBool,PetscBool,PetscBool,Cint,Ptr{DM}),arg1,arg2,arg3,arg4,arg5,PetscReal,arg6)
end

function DMPlexCreateSquareBoundary(arg1::DM,PetscReal::Ptr{Cint},arg2::Ptr{Cint},arg3::Ptr{PetscInt})
    ccall((:DMPlexCreateSquareBoundary,petsc),PetscErrorCode,(DM,Ptr{Cint},Ptr{Cint},Ptr{PetscInt}),arg1,PetscReal,arg2,arg3)
end

function DMPlexCreateCubeBoundary(arg1::DM,PetscReal::Ptr{Cint},arg2::Ptr{Cint},arg3::Ptr{PetscInt})
    ccall((:DMPlexCreateCubeBoundary,petsc),PetscErrorCode,(DM,Ptr{Cint},Ptr{Cint},Ptr{PetscInt}),arg1,PetscReal,arg2,arg3)
end

function DMPlexCreateSquareMesh(arg1::DM,PetscReal::Ptr{Cint},arg2::Ptr{Cint},arg3::Ptr{PetscInt},arg4::DMBoundaryType,arg5::DMBoundaryType)
    ccall((:DMPlexCreateSquareMesh,petsc),PetscErrorCode,(DM,Ptr{Cint},Ptr{Cint},Ptr{PetscInt},DMBoundaryType,DMBoundaryType),arg1,PetscReal,arg2,arg3,arg4,arg5)
end

function DMPlexCreateBoxMesh(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscBool,arg4::Ptr{DM})
    ccall((:DMPlexCreateBoxMesh,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscBool,Ptr{DM}),arg1,arg2,arg3,arg4)
end

function DMPlexCreateHexBoxMesh(arg1::MPI_Comm,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::DMBoundaryType,arg5::DMBoundaryType,arg6::DMBoundaryType,arg7::Ptr{DM})
    ccall((:DMPlexCreateHexBoxMesh,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{PetscInt},DMBoundaryType,DMBoundaryType,DMBoundaryType,Ptr{DM}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function DMPlexCreateConeSection(arg1::DM,arg2::Ptr{PetscSection})
    ccall((:DMPlexCreateConeSection,petsc),PetscErrorCode,(DM,Ptr{PetscSection}),arg1,arg2)
end

function DMPlexInvertCell(arg1::PetscInt,arg2::PetscInt,arg3::Ptr{Cint})
    ccall((:DMPlexInvertCell,petsc),PetscErrorCode,(PetscInt,PetscInt,Ptr{Cint}),arg1,arg2,arg3)
end

function DMPlexLocalizeCoordinate(arg1::DM,arg2::Ptr{PetscScalar},arg3::Ptr{PetscScalar})
    ccall((:DMPlexLocalizeCoordinate,petsc),PetscErrorCode,(DM,Ptr{PetscScalar},Ptr{PetscScalar}),arg1,arg2,arg3)
end

function DMPlexLocalizeCoordinates(arg1::DM)
    ccall((:DMPlexLocalizeCoordinates,petsc),PetscErrorCode,(DM,),arg1)
end

function DMPlexCheckSymmetry(arg1::DM)
    ccall((:DMPlexCheckSymmetry,petsc),PetscErrorCode,(DM,),arg1)
end

function DMPlexCheckSkeleton(arg1::DM,arg2::PetscBool,arg3::PetscInt)
    ccall((:DMPlexCheckSkeleton,petsc),PetscErrorCode,(DM,PetscBool,PetscInt),arg1,arg2,arg3)
end

function DMPlexCheckFaces(arg1::DM,arg2::PetscBool,arg3::PetscInt)
    ccall((:DMPlexCheckFaces,petsc),PetscErrorCode,(DM,PetscBool,PetscInt),arg1,arg2,arg3)
end

function DMPlexTriangleSetOptions(arg1::DM,arg2::Ptr{Uint8})
    ccall((:DMPlexTriangleSetOptions,petsc),PetscErrorCode,(DM,Ptr{Uint8}),arg1,arg2)
end

function DMPlexTetgenSetOptions(arg1::DM,arg2::Ptr{Uint8})
    ccall((:DMPlexTetgenSetOptions,petsc),PetscErrorCode,(DM,Ptr{Uint8}),arg1,arg2)
end

function DMPlexCreateNeighborCSR(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{Ptr{PetscInt}},arg5::Ptr{Ptr{PetscInt}})
    ccall((:DMPlexCreateNeighborCSR,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexGetPartitioner(arg1::DM,arg2::Ptr{PetscPartitioner})
    ccall((:DMPlexGetPartitioner,petsc),PetscErrorCode,(DM,Ptr{PetscPartitioner}),arg1,arg2)
end

function DMPlexSetPartitioner(arg1::DM,arg2::PetscPartitioner)
    ccall((:DMPlexSetPartitioner,petsc),PetscErrorCode,(DM,PetscPartitioner),arg1,arg2)
end

function DMPlexCreatePartition(arg1::DM,arg2::Ptr{Uint8},arg3::PetscInt,arg4::PetscBool,arg5::Ptr{PetscSection},arg6::Ptr{IS},arg7::Ptr{PetscSection},arg8::Ptr{IS})
    ccall((:DMPlexCreatePartition,petsc),PetscErrorCode,(DM,Ptr{Uint8},PetscInt,PetscBool,Ptr{PetscSection},Ptr{IS},Ptr{PetscSection},Ptr{IS}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function DMPlexCreatePartitionerGraph(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{Ptr{PetscInt}},arg5::Ptr{Ptr{PetscInt}})
    ccall((:DMPlexCreatePartitionerGraph,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexCreatePartitionClosure(arg1::DM,arg2::PetscSection,arg3::IS,arg4::Ptr{PetscSection},arg5::Ptr{IS})
    ccall((:DMPlexCreatePartitionClosure,petsc),PetscErrorCode,(DM,PetscSection,IS,Ptr{PetscSection},Ptr{IS}),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexPartitionLabelInvert(arg1::DM,arg2::DMLabel,arg3::PetscSF,arg4::DMLabel)
    ccall((:DMPlexPartitionLabelInvert,petsc),PetscErrorCode,(DM,DMLabel,PetscSF,DMLabel),arg1,arg2,arg3,arg4)
end

function DMPlexPartitionLabelClosure(arg1::DM,arg2::DMLabel)
    ccall((:DMPlexPartitionLabelClosure,petsc),PetscErrorCode,(DM,DMLabel),arg1,arg2)
end

function DMPlexPartitionLabelAdjacency(arg1::DM,arg2::DMLabel)
    ccall((:DMPlexPartitionLabelAdjacency,petsc),PetscErrorCode,(DM,DMLabel),arg1,arg2)
end

function DMPlexPartitionLabelCreateSF(arg1::DM,arg2::DMLabel,arg3::Ptr{PetscSF})
    ccall((:DMPlexPartitionLabelCreateSF,petsc),PetscErrorCode,(DM,DMLabel,Ptr{PetscSF}),arg1,arg2,arg3)
end

function DMPlexDistribute(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscSF},arg4::Ptr{DM})
    ccall((:DMPlexDistribute,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscSF},Ptr{DM}),arg1,arg2,arg3,arg4)
end

function DMPlexDistributeOverlap(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscSF},arg4::Ptr{DM})
    ccall((:DMPlexDistributeOverlap,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscSF},Ptr{DM}),arg1,arg2,arg3,arg4)
end

function DMPlexDistributeField(arg1::DM,arg2::PetscSF,arg3::PetscSection,arg4::Vec,arg5::PetscSection,arg6::Vec)
    ccall((:DMPlexDistributeField,petsc),PetscErrorCode,(DM,PetscSF,PetscSection,Vec,PetscSection,Vec),arg1,arg2,arg3,arg4,arg5,arg6)
end

function DMPlexDistributeFieldIS(arg1::DM,arg2::PetscSF,arg3::PetscSection,arg4::IS,arg5::PetscSection,arg6::Ptr{IS})
    ccall((:DMPlexDistributeFieldIS,petsc),PetscErrorCode,(DM,PetscSF,PetscSection,IS,PetscSection,Ptr{IS}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function DMPlexDistributeData(arg1::DM,arg2::PetscSF,arg3::PetscSection,arg4::MPI_Datatype,arg5::Ptr{Void},arg6::PetscSection,arg7::Ptr{Ptr{Void}})
    ccall((:DMPlexDistributeData,petsc),PetscErrorCode,(DM,PetscSF,PetscSection,MPI_Datatype,Ptr{Void},PetscSection,Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function DMPlexMigrate(arg1::DM,arg2::PetscSF,arg3::DM)
    ccall((:DMPlexMigrate,petsc),PetscErrorCode,(DM,PetscSF,DM),arg1,arg2,arg3)
end

function DMPlexSetAdjacencyUseCone(arg1::DM,arg2::PetscBool)
    ccall((:DMPlexSetAdjacencyUseCone,petsc),PetscErrorCode,(DM,PetscBool),arg1,arg2)
end

function DMPlexGetAdjacencyUseCone(arg1::DM,arg2::Ptr{PetscBool})
    ccall((:DMPlexGetAdjacencyUseCone,petsc),PetscErrorCode,(DM,Ptr{PetscBool}),arg1,arg2)
end

function DMPlexSetAdjacencyUseClosure(arg1::DM,arg2::PetscBool)
    ccall((:DMPlexSetAdjacencyUseClosure,petsc),PetscErrorCode,(DM,PetscBool),arg1,arg2)
end

function DMPlexGetAdjacencyUseClosure(arg1::DM,arg2::Ptr{PetscBool})
    ccall((:DMPlexGetAdjacencyUseClosure,petsc),PetscErrorCode,(DM,Ptr{PetscBool}),arg1,arg2)
end

function DMPlexSetAdjacencyUseAnchors(arg1::DM,arg2::PetscBool)
    ccall((:DMPlexSetAdjacencyUseAnchors,petsc),PetscErrorCode,(DM,PetscBool),arg1,arg2)
end

function DMPlexGetAdjacencyUseAnchors(arg1::DM,arg2::Ptr{PetscBool})
    ccall((:DMPlexGetAdjacencyUseAnchors,petsc),PetscErrorCode,(DM,Ptr{PetscBool}),arg1,arg2)
end

function DMPlexGetAdjacency(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{Ptr{PetscInt}})
    ccall((:DMPlexGetAdjacency,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4)
end

function DMPlexGetOrdering(arg1::DM,arg2::MatOrderingType,arg3::Ptr{IS})
    ccall((:DMPlexGetOrdering,petsc),PetscErrorCode,(DM,MatOrderingType,Ptr{IS}),arg1,arg2,arg3)
end

function DMPlexPermute(arg1::DM,arg2::IS,arg3::Ptr{DM})
    ccall((:DMPlexPermute,petsc),PetscErrorCode,(DM,IS,Ptr{DM}),arg1,arg2,arg3)
end

function DMPlexCreateProcessSF(arg1::DM,arg2::PetscSF,arg3::Ptr{IS},arg4::Ptr{PetscSF})
    ccall((:DMPlexCreateProcessSF,petsc),PetscErrorCode,(DM,PetscSF,Ptr{IS},Ptr{PetscSF}),arg1,arg2,arg3,arg4)
end

function DMPlexCreateTwoSidedProcessSF(arg1::DM,arg2::PetscSF,arg3::PetscSection,arg4::IS,arg5::PetscSection,arg6::IS,arg7::Ptr{IS},arg8::Ptr{PetscSF})
    ccall((:DMPlexCreateTwoSidedProcessSF,petsc),PetscErrorCode,(DM,PetscSF,PetscSection,IS,PetscSection,IS,Ptr{IS},Ptr{PetscSF}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function DMPlexDistributeOwnership(arg1::DM,arg2::PetscSection,arg3::Ptr{IS},arg4::PetscSection,arg5::Ptr{IS})
    ccall((:DMPlexDistributeOwnership,petsc),PetscErrorCode,(DM,PetscSection,Ptr{IS},PetscSection,Ptr{IS}),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexCreateOverlap(arg1::DM,arg2::PetscInt,arg3::PetscSection,arg4::IS,arg5::PetscSection,arg6::IS,arg7::Ptr{DMLabel})
    ccall((:DMPlexCreateOverlap,petsc),PetscErrorCode,(DM,PetscInt,PetscSection,IS,PetscSection,IS,Ptr{DMLabel}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function DMPlexCreateOverlapMigrationSF(arg1::DM,arg2::PetscSF,arg3::Ptr{PetscSF})
    ccall((:DMPlexCreateOverlapMigrationSF,petsc),PetscErrorCode,(DM,PetscSF,Ptr{PetscSF}),arg1,arg2,arg3)
end

function DMPlexStratifyMigrationSF(arg1::DM,arg2::PetscSF,arg3::Ptr{PetscSF})
    ccall((:DMPlexStratifyMigrationSF,petsc),PetscErrorCode,(DM,PetscSF,Ptr{PetscSF}),arg1,arg2,arg3)
end

function DMPlexCreateSubmesh(arg1::DM,arg2::DMLabel,arg3::PetscInt,arg4::Ptr{DM})
    ccall((:DMPlexCreateSubmesh,petsc),PetscErrorCode,(DM,DMLabel,PetscInt,Ptr{DM}),arg1,arg2,arg3,arg4)
end

function DMPlexCreateHybridMesh(arg1::DM,arg2::DMLabel,arg3::Ptr{DMLabel},arg4::Ptr{DM})
    ccall((:DMPlexCreateHybridMesh,petsc),PetscErrorCode,(DM,DMLabel,Ptr{DMLabel},Ptr{DM}),arg1,arg2,arg3,arg4)
end

function DMPlexGetSubpointMap(arg1::DM,arg2::Ptr{DMLabel})
    ccall((:DMPlexGetSubpointMap,petsc),PetscErrorCode,(DM,Ptr{DMLabel}),arg1,arg2)
end

function DMPlexSetSubpointMap(arg1::DM,arg2::DMLabel)
    ccall((:DMPlexSetSubpointMap,petsc),PetscErrorCode,(DM,DMLabel),arg1,arg2)
end

function DMPlexCreateSubpointIS(arg1::DM,arg2::Ptr{IS})
    ccall((:DMPlexCreateSubpointIS,petsc),PetscErrorCode,(DM,Ptr{IS}),arg1,arg2)
end

function DMPlexMarkBoundaryFaces(arg1::DM,arg2::DMLabel)
    ccall((:DMPlexMarkBoundaryFaces,petsc),PetscErrorCode,(DM,DMLabel),arg1,arg2)
end

function DMPlexLabelComplete(arg1::DM,arg2::DMLabel)
    ccall((:DMPlexLabelComplete,petsc),PetscErrorCode,(DM,DMLabel),arg1,arg2)
end

function DMPlexLabelCohesiveComplete(arg1::DM,arg2::DMLabel,arg3::DMLabel,arg4::PetscBool,arg5::DM)
    ccall((:DMPlexLabelCohesiveComplete,petsc),PetscErrorCode,(DM,DMLabel,DMLabel,PetscBool,DM),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexLabelAddCells(arg1::DM,arg2::DMLabel)
    ccall((:DMPlexLabelAddCells,petsc),PetscErrorCode,(DM,DMLabel),arg1,arg2)
end

function DMPlexGetRefinementLimit(arg1::DM,arg2::Ptr{Cint})
    ccall((:DMPlexGetRefinementLimit,petsc),PetscErrorCode,(DM,Ptr{Cint}),arg1,arg2)
end

function DMPlexSetRefinementLimit(arg1::DM,PetscReal::Cint)
    ccall((:DMPlexSetRefinementLimit,petsc),PetscErrorCode,(DM,Cint),arg1,PetscReal)
end

function DMPlexGetRefinementUniform(arg1::DM,arg2::Ptr{PetscBool})
    ccall((:DMPlexGetRefinementUniform,petsc),PetscErrorCode,(DM,Ptr{PetscBool}),arg1,arg2)
end

function DMPlexSetRefinementUniform(arg1::DM,arg2::PetscBool)
    ccall((:DMPlexSetRefinementUniform,petsc),PetscErrorCode,(DM,PetscBool),arg1,arg2)
end

function DMPlexGetCoarseDM(arg1::DM,arg2::Ptr{DM})
    ccall((:DMPlexGetCoarseDM,petsc),PetscErrorCode,(DM,Ptr{DM}),arg1,arg2)
end

function DMPlexSetCoarseDM(arg1::DM,arg2::DM)
    ccall((:DMPlexSetCoarseDM,petsc),PetscErrorCode,(DM,DM),arg1,arg2)
end

function DMPlexCreateCoarsePointIS(arg1::DM,arg2::Ptr{IS})
    ccall((:DMPlexCreateCoarsePointIS,petsc),PetscErrorCode,(DM,Ptr{IS}),arg1,arg2)
end

function DMPlexGetNumFaceVertices(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscInt})
    ccall((:DMPlexGetNumFaceVertices,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function DMPlexGetOrientedFace(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscInt},arg5::PetscInt,arg6::Ptr{PetscInt},arg7::Ptr{PetscInt},arg8::Ptr{PetscInt},arg9::Ptr{PetscBool})
    ccall((:DMPlexGetOrientedFace,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function DMPlexGetMinRadius(arg1::DM,arg2::Ptr{Cint})
    ccall((:DMPlexGetMinRadius,petsc),PetscErrorCode,(DM,Ptr{Cint}),arg1,arg2)
end

function DMPlexSetMinRadius(arg1::DM,PetscReal::Cint)
    ccall((:DMPlexSetMinRadius,petsc),PetscErrorCode,(DM,Cint),arg1,PetscReal)
end

function DMPlexComputeCellGeometryFVM(arg1::DM,arg2::PetscInt,arg3::Ptr{Cint},PetscReal::Ptr{Cint},arg4::Ptr{Cint})
    ccall((:DMPlexComputeCellGeometryFVM,petsc),PetscErrorCode,(DM,PetscInt,Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,PetscReal,arg4)
end

function DMPlexComputeGeometryFVM(arg1::DM,arg2::Ptr{Vec},arg3::Ptr{Vec})
    ccall((:DMPlexComputeGeometryFVM,petsc),PetscErrorCode,(DM,Ptr{Vec},Ptr{Vec}),arg1,arg2,arg3)
end

function DMPlexComputeGradientFVM(arg1::DM,arg2::PetscFV,arg3::Vec,arg4::Vec,arg5::Ptr{DM})
    ccall((:DMPlexComputeGradientFVM,petsc),PetscErrorCode,(DM,PetscFV,Vec,Vec,Ptr{DM}),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexInsertBoundaryValues(arg1::DM,arg2::Vec,PetscReal::Cint,arg3::Vec,arg4::Vec,arg5::Vec)
    ccall((:DMPlexInsertBoundaryValues,petsc),PetscErrorCode,(DM,Vec,Cint,Vec,Vec,Vec),arg1,arg2,PetscReal,arg3,arg4,arg5)
end

function DMPlexCreateSection(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{PetscInt},arg5::Ptr{PetscInt},arg6::PetscInt,arg7::Ptr{PetscInt},arg8::Ptr{IS},arg9::Ptr{IS},arg10::IS,arg11::Ptr{PetscSection})
    ccall((:DMPlexCreateSection,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{IS},Ptr{IS},IS,Ptr{PetscSection}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
end

function DMPlexComputeCellGeometryAffineFEM(arg1::DM,arg2::PetscInt,arg3::Ptr{Cint},arg4::Ptr{Cint},arg5::Ptr{Cint},arg6::Ptr{Cint})
    ccall((:DMPlexComputeCellGeometryAffineFEM,petsc),PetscErrorCode,(DM,PetscInt,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function DMPlexComputeCellGeometryFEM(arg1::DM,arg2::PetscInt,arg3::PetscFE,arg4::Ptr{Cint},arg5::Ptr{Cint},arg6::Ptr{Cint},arg7::Ptr{Cint})
    ccall((:DMPlexComputeCellGeometryFEM,petsc),PetscErrorCode,(DM,PetscInt,PetscFE,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function DMPlexComputeGeometryFEM(arg1::DM,arg2::Ptr{Vec})
    ccall((:DMPlexComputeGeometryFEM,petsc),PetscErrorCode,(DM,Ptr{Vec}),arg1,arg2)
end

function DMPlexVecGetClosure(arg1::DM,arg2::PetscSection,arg3::Vec,arg4::PetscInt,arg5::Ptr{PetscInt},arg6::Ptr{Ptr{PetscScalar}})
    ccall((:DMPlexVecGetClosure,petsc),PetscErrorCode,(DM,PetscSection,Vec,PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscScalar}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function DMPlexVecRestoreClosure(arg1::DM,arg2::PetscSection,arg3::Vec,arg4::PetscInt,arg5::Ptr{PetscInt},arg6::Ptr{Ptr{PetscScalar}})
    ccall((:DMPlexVecRestoreClosure,petsc),PetscErrorCode,(DM,PetscSection,Vec,PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscScalar}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function DMPlexVecSetClosure(arg1::DM,arg2::PetscSection,arg3::Vec,arg4::PetscInt,arg5::Ptr{PetscScalar},arg6::InsertMode)
    ccall((:DMPlexVecSetClosure,petsc),PetscErrorCode,(DM,PetscSection,Vec,PetscInt,Ptr{PetscScalar},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6)
end

function DMPlexMatSetClosure(arg1::DM,arg2::PetscSection,arg3::PetscSection,arg4::Mat,arg5::PetscInt,arg6::Ptr{PetscScalar},arg7::InsertMode)
    ccall((:DMPlexMatSetClosure,petsc),PetscErrorCode,(DM,PetscSection,PetscSection,Mat,PetscInt,Ptr{PetscScalar},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function DMPlexMatSetClosureRefined(arg1::DM,arg2::PetscSection,arg3::PetscSection,arg4::DM,arg5::PetscSection,arg6::PetscSection,arg7::Mat,arg8::PetscInt,arg9::Ptr{PetscScalar},arg10::InsertMode)
    ccall((:DMPlexMatSetClosureRefined,petsc),PetscErrorCode,(DM,PetscSection,PetscSection,DM,PetscSection,PetscSection,Mat,PetscInt,Ptr{PetscScalar},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function DMPlexMatGetClosureIndicesRefined(arg1::DM,arg2::PetscSection,arg3::PetscSection,arg4::DM,arg5::PetscSection,arg6::PetscSection,arg7::PetscInt,arg8::Ptr{PetscInt},arg9::Ptr{PetscInt})
    ccall((:DMPlexMatGetClosureIndicesRefined,petsc),PetscErrorCode,(DM,PetscSection,PetscSection,DM,PetscSection,PetscSection,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function DMPlexCreateClosureIndex(arg1::DM,arg2::PetscSection)
    ccall((:DMPlexCreateClosureIndex,petsc),PetscErrorCode,(DM,PetscSection),arg1,arg2)
end

function DMPlexCreateFromFile(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::PetscBool,arg4::Ptr{DM})
    ccall((:DMPlexCreateFromFile,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},PetscBool,Ptr{DM}),arg1,arg2,arg3,arg4)
end

function DMPlexCreateExodus(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscBool,arg4::Ptr{DM})
    ccall((:DMPlexCreateExodus,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscBool,Ptr{DM}),arg1,arg2,arg3,arg4)
end

function DMPlexCreateExodusFromFile(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::PetscBool,arg4::Ptr{DM})
    ccall((:DMPlexCreateExodusFromFile,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},PetscBool,Ptr{DM}),arg1,arg2,arg3,arg4)
end

function DMPlexCreateCGNS(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscBool,arg4::Ptr{DM})
    ccall((:DMPlexCreateCGNS,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscBool,Ptr{DM}),arg1,arg2,arg3,arg4)
end

function DMPlexCreateCGNSFromFile(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::PetscBool,arg4::Ptr{DM})
    ccall((:DMPlexCreateCGNSFromFile,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},PetscBool,Ptr{DM}),arg1,arg2,arg3,arg4)
end

function DMPlexCreateGmsh(arg1::MPI_Comm,arg2::PetscViewer,arg3::PetscBool,arg4::Ptr{DM})
    ccall((:DMPlexCreateGmsh,petsc),PetscErrorCode,(MPI_Comm,PetscViewer,PetscBool,Ptr{DM}),arg1,arg2,arg3,arg4)
end

function DMPlexCreateGmshFromFile(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::PetscBool,arg4::Ptr{DM})
    ccall((:DMPlexCreateGmshFromFile,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},PetscBool,Ptr{DM}),arg1,arg2,arg3,arg4)
end

function DMPlexCreateFluent(arg1::MPI_Comm,arg2::PetscViewer,arg3::PetscBool,arg4::Ptr{DM})
    ccall((:DMPlexCreateFluent,petsc),PetscErrorCode,(MPI_Comm,PetscViewer,PetscBool,Ptr{DM}),arg1,arg2,arg3,arg4)
end

function DMPlexCreateFluentFromFile(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::PetscBool,arg4::Ptr{DM})
    ccall((:DMPlexCreateFluentFromFile,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},PetscBool,Ptr{DM}),arg1,arg2,arg3,arg4)
end

function DMPlexConstructGhostCells(arg1::DM,arg2::Ptr{Uint8},arg3::Ptr{PetscInt},arg4::Ptr{DM})
    ccall((:DMPlexConstructGhostCells,petsc),PetscErrorCode,(DM,Ptr{Uint8},Ptr{PetscInt},Ptr{DM}),arg1,arg2,arg3,arg4)
end

function DMPlexConstructCohesiveCells(arg1::DM,arg2::DMLabel,arg3::Ptr{DM})
    ccall((:DMPlexConstructCohesiveCells,petsc),PetscErrorCode,(DM,DMLabel,Ptr{DM}),arg1,arg2,arg3)
end

function DMPlexGetHybridBounds(arg1::DM,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},arg5::Ptr{PetscInt})
    ccall((:DMPlexGetHybridBounds,petsc),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexSetHybridBounds(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt)
    ccall((:DMPlexSetHybridBounds,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexGetVTKCellHeight(arg1::DM,arg2::Ptr{PetscInt})
    ccall((:DMPlexGetVTKCellHeight,petsc),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end

function DMPlexSetVTKCellHeight(arg1::DM,arg2::PetscInt)
    ccall((:DMPlexSetVTKCellHeight,petsc),PetscErrorCode,(DM,PetscInt),arg1,arg2)
end

function DMPlexVTKWriteAll(arg1::PetscObject,arg2::PetscViewer)
    ccall((:DMPlexVTKWriteAll,petsc),PetscErrorCode,(PetscObject,PetscViewer),arg1,arg2)
end

function DMPlexGetScale(arg1::DM,arg2::PetscUnit,arg3::Ptr{Cint})
    ccall((:DMPlexGetScale,petsc),PetscErrorCode,(DM,PetscUnit,Ptr{Cint}),arg1,arg2,arg3)
end

function DMPlexSetScale(arg1::DM,arg2::PetscUnit,PetscReal::Cint)
    ccall((:DMPlexSetScale,petsc),PetscErrorCode,(DM,PetscUnit,Cint),arg1,arg2,PetscReal)
end

function DMPlexAddBoundary(arg1::DM,arg2::PetscBool,arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::PetscInt,arg6::PetscInt,arg7::Ptr{PetscInt},arg8::Ptr{Void},arg9::PetscInt,arg10::Ptr{PetscInt},arg11::Ptr{Void})
    ccall((:DMPlexAddBoundary,petsc),PetscErrorCode,(DM,PetscBool,Ptr{Uint8},Ptr{Uint8},PetscInt,PetscInt,Ptr{PetscInt},Ptr{Void},PetscInt,Ptr{PetscInt},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
end

function DMPlexGetNumBoundary(arg1::DM,arg2::Ptr{PetscInt})
    ccall((:DMPlexGetNumBoundary,petsc),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end

function DMPlexGetBoundary(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscBool},arg4::Ptr{Ptr{Uint8}},arg5::Ptr{Ptr{Uint8}},arg6::Ptr{PetscInt},arg7::Ptr{PetscInt},arg8::Ptr{Ptr{PetscInt}},arg9::Ptr{Ptr{Void}},arg10::Ptr{PetscInt},arg11::Ptr{Ptr{PetscInt}},arg12::Ptr{Ptr{Void}})
    ccall((:DMPlexGetBoundary,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscBool},Ptr{Ptr{Uint8}},Ptr{Ptr{Uint8}},Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{Void}},Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12)
end

function DMPlexIsBoundaryPoint(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscBool})
    ccall((:DMPlexIsBoundaryPoint,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscBool}),arg1,arg2,arg3)
end

function DMPlexCopyBoundary(arg1::DM,arg2::DM)
    ccall((:DMPlexCopyBoundary,petsc),PetscErrorCode,(DM,DM),arg1,arg2)
end

function DMPlexInsertBoundaryValuesFEM(arg1::DM,arg2::Vec)
    ccall((:DMPlexInsertBoundaryValuesFEM,petsc),PetscErrorCode,(DM,Vec),arg1,arg2)
end

function DMPlexSetMaxProjectionHeight(arg1::DM,arg2::PetscInt)
    ccall((:DMPlexSetMaxProjectionHeight,petsc),PetscErrorCode,(DM,PetscInt),arg1,arg2)
end

function DMPlexGetMaxProjectionHeight(arg1::DM,arg2::Ptr{PetscInt})
    ccall((:DMPlexGetMaxProjectionHeight,petsc),PetscErrorCode,(DM,Ptr{PetscInt}),arg1,arg2)
end

function DMPlexProjectFunction(arg1::DM,arg2::Ptr{Ptr{Void}},arg3::Ptr{Ptr{Void}},arg4::InsertMode,arg5::Vec)
    ccall((:DMPlexProjectFunction,petsc),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}},InsertMode,Vec),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexProjectFunctionLocal(arg1::DM,arg2::Ptr{Ptr{Void}},arg3::Ptr{Ptr{Void}},arg4::InsertMode,arg5::Vec)
    ccall((:DMPlexProjectFunctionLocal,petsc),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}},InsertMode,Vec),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexProjectFieldLocal(arg1::DM,arg2::Vec,arg3::Ptr{Ptr{Void}},arg4::InsertMode,arg5::Vec)
    ccall((:DMPlexProjectFieldLocal,petsc),PetscErrorCode,(DM,Vec,Ptr{Ptr{Void}},InsertMode,Vec),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexComputeL2Diff(arg1::DM,arg2::Ptr{Ptr{Void}},arg3::Ptr{Ptr{Void}},arg4::Vec,arg5::Ptr{Cint})
    ccall((:DMPlexComputeL2Diff,petsc),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}},Vec,Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexComputeL2GradientDiff(arg1::DM,arg2::Ptr{Ptr{Void}},arg3::Ptr{Ptr{Void}},arg4::Vec,PetscReal::Ptr{Cint},arg5::Ptr{Cint})
    ccall((:DMPlexComputeL2GradientDiff,petsc),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}},Vec,Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,PetscReal,arg5)
end

function DMPlexComputeL2FieldDiff(arg1::DM,arg2::Ptr{Ptr{Void}},arg3::Ptr{Ptr{Void}},arg4::Vec,PetscReal::Ptr{Cint})
    ccall((:DMPlexComputeL2FieldDiff,petsc),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}},Vec,Ptr{Cint}),arg1,arg2,arg3,arg4,PetscReal)
end

function DMPlexComputeIntegralFEM(arg1::DM,arg2::Vec,arg3::Ptr{Cint},arg4::Ptr{Void})
    ccall((:DMPlexComputeIntegralFEM,petsc),PetscErrorCode,(DM,Vec,Ptr{Cint},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function DMPlexComputeInterpolatorFEM(arg1::DM,arg2::DM,arg3::Mat,arg4::Ptr{Void})
    ccall((:DMPlexComputeInterpolatorFEM,petsc),PetscErrorCode,(DM,DM,Mat,Ptr{Void}),arg1,arg2,arg3,arg4)
end

function DMPlexComputeInjectorFEM(arg1::DM,arg2::DM,arg3::Ptr{VecScatter},arg4::Ptr{Void})
    ccall((:DMPlexComputeInjectorFEM,petsc),PetscErrorCode,(DM,DM,Ptr{VecScatter},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function DMPlexCreateRigidBody(arg1::DM,arg2::Ptr{MatNullSpace})
    ccall((:DMPlexCreateRigidBody,petsc),PetscErrorCode,(DM,Ptr{MatNullSpace}),arg1,arg2)
end

function DMPlexSNESComputeResidualFEM(arg1::DM,arg2::Vec,arg3::Vec,arg4::Ptr{Void})
    ccall((:DMPlexSNESComputeResidualFEM,petsc),PetscErrorCode,(DM,Vec,Vec,Ptr{Void}),arg1,arg2,arg3,arg4)
end

function DMPlexSNESComputeJacobianFEM(arg1::DM,arg2::Vec,arg3::Mat,arg4::Mat,arg5::Ptr{Void})
    ccall((:DMPlexSNESComputeJacobianFEM,petsc),PetscErrorCode,(DM,Vec,Mat,Mat,Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexTSComputeRHSFunctionFVM(arg1::DM,PetscReal::Cint,arg2::Vec,arg3::Vec,arg4::Ptr{Void})
    ccall((:DMPlexTSComputeRHSFunctionFVM,petsc),PetscErrorCode,(DM,Cint,Vec,Vec,Ptr{Void}),arg1,PetscReal,arg2,arg3,arg4)
end

function DMPlexTSComputeIFunctionFEM(arg1::DM,PetscReal::Cint,arg2::Vec,arg3::Vec,arg4::Vec,arg5::Ptr{Void})
    ccall((:DMPlexTSComputeIFunctionFEM,petsc),PetscErrorCode,(DM,Cint,Vec,Vec,Vec,Ptr{Void}),arg1,PetscReal,arg2,arg3,arg4,arg5)
end

function DMPlexComputeRHSFunctionFVM(arg1::DM,PetscReal::Cint,arg2::Vec,arg3::Vec,arg4::Ptr{Void})
    ccall((:DMPlexComputeRHSFunctionFVM,petsc),PetscErrorCode,(DM,Cint,Vec,Vec,Ptr{Void}),arg1,PetscReal,arg2,arg3,arg4)
end

function DMPlexGetAnchors(arg1::DM,arg2::Ptr{PetscSection},arg3::Ptr{IS})
    ccall((:DMPlexGetAnchors,petsc),PetscErrorCode,(DM,Ptr{PetscSection},Ptr{IS}),arg1,arg2,arg3)
end

function DMPlexSetAnchors(arg1::DM,arg2::PetscSection,arg3::IS)
    ccall((:DMPlexSetAnchors,petsc),PetscErrorCode,(DM,PetscSection,IS),arg1,arg2,arg3)
end

function DMPlexSetReferenceTree(arg1::DM,arg2::DM)
    ccall((:DMPlexSetReferenceTree,petsc),PetscErrorCode,(DM,DM),arg1,arg2)
end

function DMPlexGetReferenceTree(arg1::DM,arg2::Ptr{DM})
    ccall((:DMPlexGetReferenceTree,petsc),PetscErrorCode,(DM,Ptr{DM}),arg1,arg2)
end

function DMPlexReferenceTreeGetChildSymmetry(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::PetscInt,arg7::Ptr{PetscInt},arg8::Ptr{PetscInt})
    ccall((:DMPlexReferenceTreeGetChildSymmetry,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function DMPlexCreateDefaultReferenceTree(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscBool,arg4::Ptr{DM})
    ccall((:DMPlexCreateDefaultReferenceTree,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscBool,Ptr{DM}),arg1,arg2,arg3,arg4)
end

function DMPlexSetTree(arg1::DM,arg2::PetscSection,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:DMPlexSetTree,petsc),PetscErrorCode,(DM,PetscSection,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function DMPlexGetTree(arg1::DM,arg2::Ptr{PetscSection},arg3::Ptr{Ptr{PetscInt}},arg4::Ptr{Ptr{PetscInt}},arg5::Ptr{PetscSection},arg6::Ptr{Ptr{PetscInt}})
    ccall((:DMPlexGetTree,petsc),PetscErrorCode,(DM,Ptr{PetscSection},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{PetscSection},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function DMPlexGetTreeParent(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:DMPlexGetTreeParent,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function DMPlexGetTreeChildren(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{Ptr{PetscInt}})
    ccall((:DMPlexGetTreeChildren,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3,arg4)
end

function DMPlexTreeRefineCell(arg1::DM,arg2::PetscInt,arg3::Ptr{DM})
    ccall((:DMPlexTreeRefineCell,petsc),PetscErrorCode,(DM,PetscInt,Ptr{DM}),arg1,arg2,arg3)
end

function DMRedundantCreate(arg1::MPI_Comm,arg2::PetscMPIInt,arg3::PetscInt,arg4::Ptr{DM})
    ccall((:DMRedundantCreate,petsc),PetscErrorCode,(MPI_Comm,PetscMPIInt,PetscInt,Ptr{DM}),arg1,arg2,arg3,arg4)
end

function DMRedundantSetSize(arg1::DM,arg2::PetscMPIInt,arg3::PetscInt)
    ccall((:DMRedundantSetSize,petsc),PetscErrorCode,(DM,PetscMPIInt,PetscInt),arg1,arg2,arg3)
end

function DMRedundantGetSize(arg1::DM,arg2::Ptr{PetscMPIInt},arg3::Ptr{PetscInt})
    ccall((:DMRedundantGetSize,petsc),PetscErrorCode,(DM,Ptr{PetscMPIInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function DMShellCreate(arg1::MPI_Comm,arg2::Ptr{DM})
    ccall((:DMShellCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{DM}),arg1,arg2)
end

function DMShellSetMatrix(arg1::DM,arg2::Mat)
    ccall((:DMShellSetMatrix,petsc),PetscErrorCode,(DM,Mat),arg1,arg2)
end

function DMShellSetGlobalVector(arg1::DM,arg2::Vec)
    ccall((:DMShellSetGlobalVector,petsc),PetscErrorCode,(DM,Vec),arg1,arg2)
end

function DMShellSetLocalVector(arg1::DM,arg2::Vec)
    ccall((:DMShellSetLocalVector,petsc),PetscErrorCode,(DM,Vec),arg1,arg2)
end

function DMShellSetCreateGlobalVector(arg1::DM,arg2::Ptr{Void})
    ccall((:DMShellSetCreateGlobalVector,petsc),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end

function DMShellSetCreateLocalVector(arg1::DM,arg2::Ptr{Void})
    ccall((:DMShellSetCreateLocalVector,petsc),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end

function DMShellSetGlobalToLocal(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:DMShellSetGlobalToLocal,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function DMShellSetGlobalToLocalVecScatter(arg1::DM,arg2::VecScatter)
    ccall((:DMShellSetGlobalToLocalVecScatter,petsc),PetscErrorCode,(DM,VecScatter),arg1,arg2)
end

function DMShellSetLocalToGlobal(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:DMShellSetLocalToGlobal,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function DMShellSetLocalToGlobalVecScatter(arg1::DM,arg2::VecScatter)
    ccall((:DMShellSetLocalToGlobalVecScatter,petsc),PetscErrorCode,(DM,VecScatter),arg1,arg2)
end

function DMShellSetLocalToLocal(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:DMShellSetLocalToLocal,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function DMShellSetLocalToLocalVecScatter(arg1::DM,arg2::VecScatter)
    ccall((:DMShellSetLocalToLocalVecScatter,petsc),PetscErrorCode,(DM,VecScatter),arg1,arg2)
end

function DMShellSetCreateMatrix(arg1::DM,arg2::Ptr{Void})
    ccall((:DMShellSetCreateMatrix,petsc),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end

function DMShellSetCoarsen(arg1::DM,arg2::Ptr{Void})
    ccall((:DMShellSetCoarsen,petsc),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end

function DMShellSetRefine(arg1::DM,arg2::Ptr{Void})
    ccall((:DMShellSetRefine,petsc),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end

function DMShellSetCreateInterpolation(arg1::DM,arg2::Ptr{Void})
    ccall((:DMShellSetCreateInterpolation,petsc),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end

function DMShellSetCreateInjection(arg1::DM,arg2::Ptr{Void})
    ccall((:DMShellSetCreateInjection,petsc),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end

function DMShellSetCreateFieldDecomposition(arg1::DM,arg2::Ptr{Void})
    ccall((:DMShellSetCreateFieldDecomposition,petsc),PetscErrorCode,(DM,Ptr{Void}),arg1,arg2)
end

function DMGlobalToLocalBeginDefaultShell(arg1::DM,arg2::Vec,arg3::InsertMode,arg4::Vec)
    ccall((:DMGlobalToLocalBeginDefaultShell,petsc),PetscErrorCode,(DM,Vec,InsertMode,Vec),arg1,arg2,arg3,arg4)
end

function DMGlobalToLocalEndDefaultShell(arg1::DM,arg2::Vec,arg3::InsertMode,arg4::Vec)
    ccall((:DMGlobalToLocalEndDefaultShell,petsc),PetscErrorCode,(DM,Vec,InsertMode,Vec),arg1,arg2,arg3,arg4)
end

function DMLocalToGlobalBeginDefaultShell(arg1::DM,arg2::Vec,arg3::InsertMode,arg4::Vec)
    ccall((:DMLocalToGlobalBeginDefaultShell,petsc),PetscErrorCode,(DM,Vec,InsertMode,Vec),arg1,arg2,arg3,arg4)
end

function DMLocalToGlobalEndDefaultShell(arg1::DM,arg2::Vec,arg3::InsertMode,arg4::Vec)
    ccall((:DMLocalToGlobalEndDefaultShell,petsc),PetscErrorCode,(DM,Vec,InsertMode,Vec),arg1,arg2,arg3,arg4)
end

function DMLocalToLocalBeginDefaultShell(arg1::DM,arg2::Vec,arg3::InsertMode,arg4::Vec)
    ccall((:DMLocalToLocalBeginDefaultShell,petsc),PetscErrorCode,(DM,Vec,InsertMode,Vec),arg1,arg2,arg3,arg4)
end

function DMLocalToLocalEndDefaultShell(arg1::DM,arg2::Vec,arg3::InsertMode,arg4::Vec)
    ccall((:DMLocalToLocalEndDefaultShell,petsc),PetscErrorCode,(DM,Vec,InsertMode,Vec),arg1,arg2,arg3,arg4)
end

function DMSlicedCreate(arg1::MPI_Comm,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::Ptr{PetscInt},arg6::Ptr{PetscInt},arg7::Ptr{PetscInt},arg8::Ptr{DM})
    ccall((:DMSlicedCreate,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{DM}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function DMSlicedSetPreallocation(arg1::DM,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::PetscInt,arg5::Ptr{PetscInt})
    ccall((:DMSlicedSetPreallocation,petsc),PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end

function DMSlicedSetBlockFills(arg1::DM,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt})
    ccall((:DMSlicedSetBlockFills,petsc),PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function DMSlicedSetGhosts(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::Ptr{PetscInt})
    ccall((:DMSlicedSetGhosts,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end

function PetscDSInitializePackage()
    ccall((:PetscDSInitializePackage,petsc),PetscErrorCode,())
end

function PetscDSCreate(arg1::MPI_Comm,arg2::Ptr{PetscDS})
    ccall((:PetscDSCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscDS}),arg1,arg2)
end

function PetscDSDestroy(arg1::Ptr{PetscDS})
    ccall((:PetscDSDestroy,petsc),PetscErrorCode,(Ptr{PetscDS},),arg1)
end

function PetscDSSetType(arg1::PetscDS,arg2::PetscDSType)
    ccall((:PetscDSSetType,petsc),PetscErrorCode,(PetscDS,PetscDSType),arg1,arg2)
end

function PetscDSGetType(arg1::PetscDS,arg2::Ptr{PetscDSType})
    ccall((:PetscDSGetType,petsc),PetscErrorCode,(PetscDS,Ptr{PetscDSType}),arg1,arg2)
end

function PetscDSSetUp(arg1::PetscDS)
    ccall((:PetscDSSetUp,petsc),PetscErrorCode,(PetscDS,),arg1)
end

function PetscDSSetFromOptions(arg1::PetscDS)
    ccall((:PetscDSSetFromOptions,petsc),PetscErrorCode,(PetscDS,),arg1)
end

function PetscDSRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:PetscDSRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function PetscDSRegisterDestroy()
    ccall((:PetscDSRegisterDestroy,petsc),PetscErrorCode,())
end

function PetscDSGetSpatialDimension(arg1::PetscDS,arg2::Ptr{PetscInt})
    ccall((:PetscDSGetSpatialDimension,petsc),PetscErrorCode,(PetscDS,Ptr{PetscInt}),arg1,arg2)
end

function PetscDSGetNumFields(arg1::PetscDS,arg2::Ptr{PetscInt})
    ccall((:PetscDSGetNumFields,petsc),PetscErrorCode,(PetscDS,Ptr{PetscInt}),arg1,arg2)
end

function PetscDSGetTotalDimension(arg1::PetscDS,arg2::Ptr{PetscInt})
    ccall((:PetscDSGetTotalDimension,petsc),PetscErrorCode,(PetscDS,Ptr{PetscInt}),arg1,arg2)
end

function PetscDSGetTotalBdDimension(arg1::PetscDS,arg2::Ptr{PetscInt})
    ccall((:PetscDSGetTotalBdDimension,petsc),PetscErrorCode,(PetscDS,Ptr{PetscInt}),arg1,arg2)
end

function PetscDSGetTotalComponents(arg1::PetscDS,arg2::Ptr{PetscInt})
    ccall((:PetscDSGetTotalComponents,petsc),PetscErrorCode,(PetscDS,Ptr{PetscInt}),arg1,arg2)
end

function PetscDSGetFieldOffset(arg1::PetscDS,arg2::PetscInt,arg3::Ptr{PetscInt})
    ccall((:PetscDSGetFieldOffset,petsc),PetscErrorCode,(PetscDS,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscDSGetBdFieldOffset(arg1::PetscDS,arg2::PetscInt,arg3::Ptr{PetscInt})
    ccall((:PetscDSGetBdFieldOffset,petsc),PetscErrorCode,(PetscDS,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscDSGetComponentOffset(arg1::PetscDS,arg2::PetscInt,arg3::Ptr{PetscInt})
    ccall((:PetscDSGetComponentOffset,petsc),PetscErrorCode,(PetscDS,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscDSGetComponentOffsets(arg1::PetscDS,arg2::Ptr{Ptr{PetscInt}})
    ccall((:PetscDSGetComponentOffsets,petsc),PetscErrorCode,(PetscDS,Ptr{Ptr{PetscInt}}),arg1,arg2)
end

function PetscDSGetComponentBdOffsets(arg1::PetscDS,arg2::Ptr{Ptr{PetscInt}})
    ccall((:PetscDSGetComponentBdOffsets,petsc),PetscErrorCode,(PetscDS,Ptr{Ptr{PetscInt}}),arg1,arg2)
end

function PetscDSGetComponentDerivativeOffsets(arg1::PetscDS,arg2::Ptr{Ptr{PetscInt}})
    ccall((:PetscDSGetComponentDerivativeOffsets,petsc),PetscErrorCode,(PetscDS,Ptr{Ptr{PetscInt}}),arg1,arg2)
end

function PetscDSGetComponentBdDerivativeOffsets(arg1::PetscDS,arg2::Ptr{Ptr{PetscInt}})
    ccall((:PetscDSGetComponentBdDerivativeOffsets,petsc),PetscErrorCode,(PetscDS,Ptr{Ptr{PetscInt}}),arg1,arg2)
end

function PetscDSGetDiscretization(arg1::PetscDS,arg2::PetscInt,arg3::Ptr{PetscObject})
    ccall((:PetscDSGetDiscretization,petsc),PetscErrorCode,(PetscDS,PetscInt,Ptr{PetscObject}),arg1,arg2,arg3)
end

function PetscDSSetDiscretization(arg1::PetscDS,arg2::PetscInt,arg3::PetscObject)
    ccall((:PetscDSSetDiscretization,petsc),PetscErrorCode,(PetscDS,PetscInt,PetscObject),arg1,arg2,arg3)
end

function PetscDSAddDiscretization(arg1::PetscDS,arg2::PetscObject)
    ccall((:PetscDSAddDiscretization,petsc),PetscErrorCode,(PetscDS,PetscObject),arg1,arg2)
end

function PetscDSGetBdDiscretization(arg1::PetscDS,arg2::PetscInt,arg3::Ptr{PetscObject})
    ccall((:PetscDSGetBdDiscretization,petsc),PetscErrorCode,(PetscDS,PetscInt,Ptr{PetscObject}),arg1,arg2,arg3)
end

function PetscDSSetBdDiscretization(arg1::PetscDS,arg2::PetscInt,arg3::PetscObject)
    ccall((:PetscDSSetBdDiscretization,petsc),PetscErrorCode,(PetscDS,PetscInt,PetscObject),arg1,arg2,arg3)
end

function PetscDSAddBdDiscretization(arg1::PetscDS,arg2::PetscObject)
    ccall((:PetscDSAddBdDiscretization,petsc),PetscErrorCode,(PetscDS,PetscObject),arg1,arg2)
end

function PetscDSGetImplicit(arg1::PetscDS,arg2::PetscInt,arg3::Ptr{PetscBool})
    ccall((:PetscDSGetImplicit,petsc),PetscErrorCode,(PetscDS,PetscInt,Ptr{PetscBool}),arg1,arg2,arg3)
end

function PetscDSSetImplicit(arg1::PetscDS,arg2::PetscInt,arg3::PetscBool)
    ccall((:PetscDSSetImplicit,petsc),PetscErrorCode,(PetscDS,PetscInt,PetscBool),arg1,arg2,arg3)
end

function PetscDSGetAdjacency(arg1::PetscDS,arg2::PetscInt,arg3::Ptr{PetscBool},arg4::Ptr{PetscBool})
    ccall((:PetscDSGetAdjacency,petsc),PetscErrorCode,(PetscDS,PetscInt,Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function PetscDSSetAdjacency(arg1::PetscDS,arg2::PetscInt,arg3::PetscBool,arg4::PetscBool)
    ccall((:PetscDSSetAdjacency,petsc),PetscErrorCode,(PetscDS,PetscInt,PetscBool,PetscBool),arg1,arg2,arg3,arg4)
end

function PetscDSGetObjective(arg1::PetscDS,arg2::PetscInt,arg3::Ptr{Ptr{Void}})
    ccall((:PetscDSGetObjective,petsc),PetscErrorCode,(PetscDS,PetscInt,Ptr{Ptr{Void}}),arg1,arg2,arg3)
end

function PetscDSSetObjective(arg1::PetscDS,arg2::PetscInt,arg3::Ptr{Void})
    ccall((:PetscDSSetObjective,petsc),PetscErrorCode,(PetscDS,PetscInt,Ptr{Void}),arg1,arg2,arg3)
end

function PetscDSGetResidual(arg1::PetscDS,arg2::PetscInt,arg3::Ptr{Ptr{Void}},arg4::Ptr{Ptr{Void}})
    ccall((:PetscDSGetResidual,petsc),PetscErrorCode,(PetscDS,PetscInt,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4)
end

function PetscDSSetResidual(arg1::PetscDS,arg2::PetscInt,arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:PetscDSSetResidual,petsc),PetscErrorCode,(PetscDS,PetscInt,Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function PetscDSGetJacobian(arg1::PetscDS,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{Ptr{Void}},arg5::Ptr{Ptr{Void}},arg6::Ptr{Ptr{Void}},arg7::Ptr{Ptr{Void}})
    ccall((:PetscDSGetJacobian,petsc),PetscErrorCode,(PetscDS,PetscInt,PetscInt,Ptr{Ptr{Void}},Ptr{Ptr{Void}},Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscDSSetJacobian(arg1::PetscDS,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{Void},arg5::Ptr{Void},arg6::Ptr{Void},arg7::Ptr{Void})
    ccall((:PetscDSSetJacobian,petsc),PetscErrorCode,(PetscDS,PetscInt,PetscInt,Ptr{Void},Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscDSGetRiemannSolver(arg1::PetscDS,arg2::PetscInt,arg3::Ptr{Ptr{Void}})
    ccall((:PetscDSGetRiemannSolver,petsc),PetscErrorCode,(PetscDS,PetscInt,Ptr{Ptr{Void}}),arg1,arg2,arg3)
end

function PetscDSSetRiemannSolver(arg1::PetscDS,arg2::PetscInt,arg3::Ptr{Void})
    ccall((:PetscDSSetRiemannSolver,petsc),PetscErrorCode,(PetscDS,PetscInt,Ptr{Void}),arg1,arg2,arg3)
end

function PetscDSGetContext(arg1::PetscDS,arg2::PetscInt,arg3::Ptr{Ptr{Void}})
    ccall((:PetscDSGetContext,petsc),PetscErrorCode,(PetscDS,PetscInt,Ptr{Ptr{Void}}),arg1,arg2,arg3)
end

function PetscDSSetContext(arg1::PetscDS,arg2::PetscInt,arg3::Ptr{Void})
    ccall((:PetscDSSetContext,petsc),PetscErrorCode,(PetscDS,PetscInt,Ptr{Void}),arg1,arg2,arg3)
end

function PetscDSGetBdResidual(arg1::PetscDS,arg2::PetscInt,arg3::Ptr{Ptr{Void}},arg4::Ptr{Ptr{Void}})
    ccall((:PetscDSGetBdResidual,petsc),PetscErrorCode,(PetscDS,PetscInt,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4)
end

function PetscDSSetBdResidual(arg1::PetscDS,arg2::PetscInt,arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:PetscDSSetBdResidual,petsc),PetscErrorCode,(PetscDS,PetscInt,Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function PetscDSGetBdJacobian(arg1::PetscDS,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{Ptr{Void}},arg5::Ptr{Ptr{Void}},arg6::Ptr{Ptr{Void}},arg7::Ptr{Ptr{Void}})
    ccall((:PetscDSGetBdJacobian,petsc),PetscErrorCode,(PetscDS,PetscInt,PetscInt,Ptr{Ptr{Void}},Ptr{Ptr{Void}},Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscDSSetBdJacobian(arg1::PetscDS,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{Void},arg5::Ptr{Void},arg6::Ptr{Void},arg7::Ptr{Void})
    ccall((:PetscDSSetBdJacobian,petsc),PetscErrorCode,(PetscDS,PetscInt,PetscInt,Ptr{Void},Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscDSGetTabulation(arg1::PetscDS,arg2::Ptr{Ptr{Ptr{Cint}}},arg3::Ptr{Ptr{Ptr{Cint}}})
    ccall((:PetscDSGetTabulation,petsc),PetscErrorCode,(PetscDS,Ptr{Ptr{Ptr{Cint}}},Ptr{Ptr{Ptr{Cint}}}),arg1,arg2,arg3)
end

function PetscDSGetBdTabulation(arg1::PetscDS,arg2::Ptr{Ptr{Ptr{Cint}}},arg3::Ptr{Ptr{Ptr{Cint}}})
    ccall((:PetscDSGetBdTabulation,petsc),PetscErrorCode,(PetscDS,Ptr{Ptr{Ptr{Cint}}},Ptr{Ptr{Ptr{Cint}}}),arg1,arg2,arg3)
end

function PetscDSGetEvaluationArrays(arg1::PetscDS,arg2::Ptr{Ptr{PetscScalar}},arg3::Ptr{Ptr{PetscScalar}},arg4::Ptr{Ptr{PetscScalar}})
    ccall((:PetscDSGetEvaluationArrays,petsc),PetscErrorCode,(PetscDS,Ptr{Ptr{PetscScalar}},Ptr{Ptr{PetscScalar}},Ptr{Ptr{PetscScalar}}),arg1,arg2,arg3,arg4)
end

function PetscDSGetWeakFormArrays(arg1::PetscDS,arg2::Ptr{Ptr{PetscScalar}},arg3::Ptr{Ptr{PetscScalar}},arg4::Ptr{Ptr{PetscScalar}},arg5::Ptr{Ptr{PetscScalar}},arg6::Ptr{Ptr{PetscScalar}},arg7::Ptr{Ptr{PetscScalar}})
    ccall((:PetscDSGetWeakFormArrays,petsc),PetscErrorCode,(PetscDS,Ptr{Ptr{PetscScalar}},Ptr{Ptr{PetscScalar}},Ptr{Ptr{PetscScalar}},Ptr{Ptr{PetscScalar}},Ptr{Ptr{PetscScalar}},Ptr{Ptr{PetscScalar}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscDSGetRefCoordArrays(arg1::PetscDS,arg2::Ptr{Ptr{Cint}},arg3::Ptr{Ptr{PetscScalar}})
    ccall((:PetscDSGetRefCoordArrays,petsc),PetscErrorCode,(PetscDS,Ptr{Ptr{Cint}},Ptr{Ptr{PetscScalar}}),arg1,arg2,arg3)
end

function CharacteristicInitializePackage()
    ccall((:CharacteristicInitializePackage,petsc),PetscErrorCode,())
end

function CharacteristicCreate(arg1::MPI_Comm,arg2::Ptr{Characteristic})
    ccall((:CharacteristicCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{Characteristic}),arg1,arg2)
end

function CharacteristicSetType(arg1::Characteristic,arg2::CharacteristicType)
    ccall((:CharacteristicSetType,petsc),PetscErrorCode,(Characteristic,CharacteristicType),arg1,arg2)
end

function CharacteristicSetUp(arg1::Characteristic)
    ccall((:CharacteristicSetUp,petsc),PetscErrorCode,(Characteristic,),arg1)
end

function CharacteristicSetVelocityInterpolation(arg1::Characteristic,arg2::DM,arg3::Vec,arg4::Vec,arg5::PetscInt,arg6::Ptr{PetscInt},arg7::Ptr{Void},arg8::Ptr{Void})
    ccall((:CharacteristicSetVelocityInterpolation,petsc),PetscErrorCode,(Characteristic,DM,Vec,Vec,PetscInt,Ptr{PetscInt},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function CharacteristicSetVelocityInterpolationLocal(arg1::Characteristic,arg2::DM,arg3::Vec,arg4::Vec,arg5::PetscInt,arg6::Ptr{PetscInt},arg7::Ptr{Void},arg8::Ptr{Void})
    ccall((:CharacteristicSetVelocityInterpolationLocal,petsc),PetscErrorCode,(Characteristic,DM,Vec,Vec,PetscInt,Ptr{PetscInt},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function CharacteristicSetFieldInterpolation(arg1::Characteristic,arg2::DM,arg3::Vec,arg4::PetscInt,arg5::Ptr{PetscInt},arg6::Ptr{Void},arg7::Ptr{Void})
    ccall((:CharacteristicSetFieldInterpolation,petsc),PetscErrorCode,(Characteristic,DM,Vec,PetscInt,Ptr{PetscInt},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function CharacteristicSetFieldInterpolationLocal(arg1::Characteristic,arg2::DM,arg3::Vec,arg4::PetscInt,arg5::Ptr{PetscInt},arg6::Ptr{Void},arg7::Ptr{Void})
    ccall((:CharacteristicSetFieldInterpolationLocal,petsc),PetscErrorCode,(Characteristic,DM,Vec,PetscInt,Ptr{PetscInt},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function CharacteristicSolve(arg1::Characteristic,PetscReal::Cint,arg2::Vec)
    ccall((:CharacteristicSolve,petsc),PetscErrorCode,(Characteristic,Cint,Vec),arg1,PetscReal,arg2)
end

function CharacteristicDestroy(arg1::Ptr{Characteristic})
    ccall((:CharacteristicDestroy,petsc),PetscErrorCode,(Ptr{Characteristic},),arg1)
end

function CharacteristicRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:CharacteristicRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function PCExoticSetType(arg1::PC,arg2::PCExoticType)
    ccall((:PCExoticSetType,petsc),PetscErrorCode,(PC,PCExoticType),arg1,arg2)
end

function PCInitializePackage()
    ccall((:PCInitializePackage,petsc),PetscErrorCode,())
end

function PCCreate(arg1::MPI_Comm,arg2::Ptr{PC})
    ccall((:PCCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{PC}),arg1,arg2)
end

function PCSetType(arg1::PC,arg2::PCType)
    ccall((:PCSetType,petsc),PetscErrorCode,(PC,PCType),arg1,arg2)
end

function PCGetType(arg1::PC,arg2::Ptr{PCType})
    ccall((:PCGetType,petsc),PetscErrorCode,(PC,Ptr{PCType}),arg1,arg2)
end

function PCSetUp(arg1::PC)
    ccall((:PCSetUp,petsc),PetscErrorCode,(PC,),arg1)
end

function PCGetSetUpFailedReason(arg1::PC,arg2::Ptr{PetscInt})
    ccall((:PCGetSetUpFailedReason,petsc),PetscErrorCode,(PC,Ptr{PetscInt}),arg1,arg2)
end

function PCSetUpOnBlocks(arg1::PC)
    ccall((:PCSetUpOnBlocks,petsc),PetscErrorCode,(PC,),arg1)
end

function PCApply(arg1::PC,arg2::Vec,arg3::Vec)
    ccall((:PCApply,petsc),PetscErrorCode,(PC,Vec,Vec),arg1,arg2,arg3)
end

function PCApplySymmetricLeft(arg1::PC,arg2::Vec,arg3::Vec)
    ccall((:PCApplySymmetricLeft,petsc),PetscErrorCode,(PC,Vec,Vec),arg1,arg2,arg3)
end

function PCApplySymmetricRight(arg1::PC,arg2::Vec,arg3::Vec)
    ccall((:PCApplySymmetricRight,petsc),PetscErrorCode,(PC,Vec,Vec),arg1,arg2,arg3)
end

function PCApplyBAorAB(arg1::PC,arg2::PCSide,arg3::Vec,arg4::Vec,arg5::Vec)
    ccall((:PCApplyBAorAB,petsc),PetscErrorCode,(PC,PCSide,Vec,Vec,Vec),arg1,arg2,arg3,arg4,arg5)
end

function PCApplyTranspose(arg1::PC,arg2::Vec,arg3::Vec)
    ccall((:PCApplyTranspose,petsc),PetscErrorCode,(PC,Vec,Vec),arg1,arg2,arg3)
end

function PCApplyTransposeExists(arg1::PC,arg2::Ptr{PetscBool})
    ccall((:PCApplyTransposeExists,petsc),PetscErrorCode,(PC,Ptr{PetscBool}),arg1,arg2)
end

function PCApplyBAorABTranspose(arg1::PC,arg2::PCSide,arg3::Vec,arg4::Vec,arg5::Vec)
    ccall((:PCApplyBAorABTranspose,petsc),PetscErrorCode,(PC,PCSide,Vec,Vec,Vec),arg1,arg2,arg3,arg4,arg5)
end

function PCSetReusePreconditioner(arg1::PC,arg2::PetscBool)
    ccall((:PCSetReusePreconditioner,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCGetReusePreconditioner(arg1::PC,arg2::Ptr{PetscBool})
    ccall((:PCGetReusePreconditioner,petsc),PetscErrorCode,(PC,Ptr{PetscBool}),arg1,arg2)
end

function PCSetErrorIfFailure(arg1::PC,arg2::PetscBool)
    ccall((:PCSetErrorIfFailure,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCApplyRichardson(arg1::PC,arg2::Vec,arg3::Vec,arg4::Vec,PetscReal::Cint,arg5::Cint,arg6::Cint,arg7::PetscInt,arg8::PetscBool,arg9::Ptr{PetscInt},arg10::Ptr{PCRichardsonConvergedReason})
    ccall((:PCApplyRichardson,petsc),PetscErrorCode,(PC,Vec,Vec,Vec,Cint,Cint,Cint,PetscInt,PetscBool,Ptr{PetscInt},Ptr{PCRichardsonConvergedReason}),arg1,arg2,arg3,arg4,PetscReal,arg5,arg6,arg7,arg8,arg9,arg10)
end

function PCApplyRichardsonExists(arg1::PC,arg2::Ptr{PetscBool})
    ccall((:PCApplyRichardsonExists,petsc),PetscErrorCode,(PC,Ptr{PetscBool}),arg1,arg2)
end

function PCSetInitialGuessNonzero(arg1::PC,arg2::PetscBool)
    ccall((:PCSetInitialGuessNonzero,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCGetInitialGuessNonzero(arg1::PC,arg2::Ptr{PetscBool})
    ccall((:PCGetInitialGuessNonzero,petsc),PetscErrorCode,(PC,Ptr{PetscBool}),arg1,arg2)
end

function PCSetUseAmat(arg1::PC,arg2::PetscBool)
    ccall((:PCSetUseAmat,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCGetUseAmat(arg1::PC,arg2::Ptr{PetscBool})
    ccall((:PCGetUseAmat,petsc),PetscErrorCode,(PC,Ptr{PetscBool}),arg1,arg2)
end

function PCRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:PCRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function PCReset(arg1::PC)
    ccall((:PCReset,petsc),PetscErrorCode,(PC,),arg1)
end

function PCDestroy(arg1::Ptr{PC})
    ccall((:PCDestroy,petsc),PetscErrorCode,(Ptr{PC},),arg1)
end

function PCSetFromOptions(arg1::PC)
    ccall((:PCSetFromOptions,petsc),PetscErrorCode,(PC,),arg1)
end

function PCFactorGetMatrix(arg1::PC,arg2::Ptr{Mat})
    ccall((:PCFactorGetMatrix,petsc),PetscErrorCode,(PC,Ptr{Mat}),arg1,arg2)
end

function PCSetModifySubMatrices(arg1::PC,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:PCSetModifySubMatrices,petsc),PetscErrorCode,(PC,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function PCModifySubMatrices(arg1::PC,arg2::PetscInt,arg3::Ptr{IS},arg4::Ptr{IS},arg5::Ptr{Mat},arg6::Ptr{Void})
    ccall((:PCModifySubMatrices,petsc),PetscErrorCode,(PC,PetscInt,Ptr{IS},Ptr{IS},Ptr{Mat},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PCSetOperators(arg1::PC,arg2::Mat,arg3::Mat)
    ccall((:PCSetOperators,petsc),PetscErrorCode,(PC,Mat,Mat),arg1,arg2,arg3)
end

function PCGetOperators(arg1::PC,arg2::Ptr{Mat},arg3::Ptr{Mat})
    ccall((:PCGetOperators,petsc),PetscErrorCode,(PC,Ptr{Mat},Ptr{Mat}),arg1,arg2,arg3)
end

function PCGetOperatorsSet(arg1::PC,arg2::Ptr{PetscBool},arg3::Ptr{PetscBool})
    ccall((:PCGetOperatorsSet,petsc),PetscErrorCode,(PC,Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3)
end

function PCView(arg1::PC,arg2::PetscViewer)
    ccall((:PCView,petsc),PetscErrorCode,(PC,PetscViewer),arg1,arg2)
end

function PCLoad(arg1::PC,arg2::PetscViewer)
    ccall((:PCLoad,petsc),PetscErrorCode,(PC,PetscViewer),arg1,arg2)
end

function PCAppendOptionsPrefix(arg1::PC,arg2::Ptr{Uint8})
    ccall((:PCAppendOptionsPrefix,petsc),PetscErrorCode,(PC,Ptr{Uint8}),arg1,arg2)
end

function PCGetOptionsPrefix(arg1::PC,arg2::Ptr{Ptr{Uint8}})
    ccall((:PCGetOptionsPrefix,petsc),PetscErrorCode,(PC,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function PCComputeExplicitOperator(arg1::PC,arg2::Ptr{Mat})
    ccall((:PCComputeExplicitOperator,petsc),PetscErrorCode,(PC,Ptr{Mat}),arg1,arg2)
end

function PCGetDiagonalScale(arg1::PC,arg2::Ptr{PetscBool})
    ccall((:PCGetDiagonalScale,petsc),PetscErrorCode,(PC,Ptr{PetscBool}),arg1,arg2)
end

function PCDiagonalScaleLeft(arg1::PC,arg2::Vec,arg3::Vec)
    ccall((:PCDiagonalScaleLeft,petsc),PetscErrorCode,(PC,Vec,Vec),arg1,arg2,arg3)
end

function PCDiagonalScaleRight(arg1::PC,arg2::Vec,arg3::Vec)
    ccall((:PCDiagonalScaleRight,petsc),PetscErrorCode,(PC,Vec,Vec),arg1,arg2,arg3)
end

function PCSetDiagonalScale(arg1::PC,arg2::Vec)
    ccall((:PCSetDiagonalScale,petsc),PetscErrorCode,(PC,Vec),arg1,arg2)
end

function PCJacobiSetType(arg1::PC,arg2::PCJacobiType)
    ccall((:PCJacobiSetType,petsc),PetscErrorCode,(PC,PCJacobiType),arg1,arg2)
end

function PCJacobiGetType(arg1::PC,arg2::Ptr{PCJacobiType})
    ccall((:PCJacobiGetType,petsc),PetscErrorCode,(PC,Ptr{PCJacobiType}),arg1,arg2)
end

function PCJacobiSetUseAbs(arg1::PC,arg2::PetscBool)
    ccall((:PCJacobiSetUseAbs,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCJacobiGetUseAbs(arg1::PC,arg2::Ptr{PetscBool})
    ccall((:PCJacobiGetUseAbs,petsc),PetscErrorCode,(PC,Ptr{PetscBool}),arg1,arg2)
end

function PCSORSetSymmetric(arg1::PC,arg2::MatSORType)
    ccall((:PCSORSetSymmetric,petsc),PetscErrorCode,(PC,MatSORType),arg1,arg2)
end

function PCSORGetSymmetric(arg1::PC,arg2::Ptr{MatSORType})
    ccall((:PCSORGetSymmetric,petsc),PetscErrorCode,(PC,Ptr{MatSORType}),arg1,arg2)
end

function PCSORSetOmega(arg1::PC,PetscReal::Cint)
    ccall((:PCSORSetOmega,petsc),PetscErrorCode,(PC,Cint),arg1,PetscReal)
end

function PCSORGetOmega(arg1::PC,arg2::Ptr{Cint})
    ccall((:PCSORGetOmega,petsc),PetscErrorCode,(PC,Ptr{Cint}),arg1,arg2)
end

function PCSORSetIterations(arg1::PC,arg2::PetscInt,arg3::PetscInt)
    ccall((:PCSORSetIterations,petsc),PetscErrorCode,(PC,PetscInt,PetscInt),arg1,arg2,arg3)
end

function PCSORGetIterations(arg1::PC,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt})
    ccall((:PCSORGetIterations,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function PCEisenstatSetOmega(arg1::PC,PetscReal::Cint)
    ccall((:PCEisenstatSetOmega,petsc),PetscErrorCode,(PC,Cint),arg1,PetscReal)
end

function PCEisenstatGetOmega(arg1::PC,arg2::Ptr{Cint})
    ccall((:PCEisenstatGetOmega,petsc),PetscErrorCode,(PC,Ptr{Cint}),arg1,arg2)
end

function PCEisenstatSetNoDiagonalScaling(arg1::PC,arg2::PetscBool)
    ccall((:PCEisenstatSetNoDiagonalScaling,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCEisenstatGetNoDiagonalScaling(arg1::PC,arg2::Ptr{PetscBool})
    ccall((:PCEisenstatGetNoDiagonalScaling,petsc),PetscErrorCode,(PC,Ptr{PetscBool}),arg1,arg2)
end

function PCBJacobiSetTotalBlocks(arg1::PC,arg2::PetscInt,arg3::Ptr{PetscInt})
    ccall((:PCBJacobiSetTotalBlocks,petsc),PetscErrorCode,(PC,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function PCBJacobiGetTotalBlocks(arg1::PC,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{PetscInt}})
    ccall((:PCBJacobiGetTotalBlocks,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3)
end

function PCBJacobiSetLocalBlocks(arg1::PC,arg2::PetscInt,arg3::Ptr{PetscInt})
    ccall((:PCBJacobiSetLocalBlocks,petsc),PetscErrorCode,(PC,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function PCBJacobiGetLocalBlocks(arg1::PC,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{PetscInt}})
    ccall((:PCBJacobiGetLocalBlocks,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3)
end

function PCShellSetApply(arg1::PC,arg2::Ptr{Void})
    ccall((:PCShellSetApply,petsc),PetscErrorCode,(PC,Ptr{Void}),arg1,arg2)
end

function PCShellSetApplyBA(arg1::PC,arg2::Ptr{Void})
    ccall((:PCShellSetApplyBA,petsc),PetscErrorCode,(PC,Ptr{Void}),arg1,arg2)
end

function PCShellSetApplyTranspose(arg1::PC,arg2::Ptr{Void})
    ccall((:PCShellSetApplyTranspose,petsc),PetscErrorCode,(PC,Ptr{Void}),arg1,arg2)
end

function PCShellSetSetUp(arg1::PC,arg2::Ptr{Void})
    ccall((:PCShellSetSetUp,petsc),PetscErrorCode,(PC,Ptr{Void}),arg1,arg2)
end

function PCShellSetApplyRichardson(arg1::PC,arg2::Ptr{Void})
    ccall((:PCShellSetApplyRichardson,petsc),PetscErrorCode,(PC,Ptr{Void}),arg1,arg2)
end

function PCShellSetView(arg1::PC,arg2::Ptr{Void})
    ccall((:PCShellSetView,petsc),PetscErrorCode,(PC,Ptr{Void}),arg1,arg2)
end

function PCShellSetDestroy(arg1::PC,arg2::Ptr{Void})
    ccall((:PCShellSetDestroy,petsc),PetscErrorCode,(PC,Ptr{Void}),arg1,arg2)
end

function PCShellSetContext(arg1::PC,arg2::Ptr{Void})
    ccall((:PCShellSetContext,petsc),PetscErrorCode,(PC,Ptr{Void}),arg1,arg2)
end

function PCShellGetContext(arg1::PC,arg2::Ptr{Ptr{Void}})
    ccall((:PCShellGetContext,petsc),PetscErrorCode,(PC,Ptr{Ptr{Void}}),arg1,arg2)
end

function PCShellSetName(arg1::PC,arg2::Ptr{Uint8})
    ccall((:PCShellSetName,petsc),PetscErrorCode,(PC,Ptr{Uint8}),arg1,arg2)
end

function PCShellGetName(arg1::PC,arg2::Ptr{Ptr{Uint8}})
    ccall((:PCShellGetName,petsc),PetscErrorCode,(PC,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function PCFactorSetZeroPivot(arg1::PC,PetscReal::Cint)
    ccall((:PCFactorSetZeroPivot,petsc),PetscErrorCode,(PC,Cint),arg1,PetscReal)
end

function PCFactorSetShiftType(arg1::PC,arg2::MatFactorShiftType)
    ccall((:PCFactorSetShiftType,petsc),PetscErrorCode,(PC,MatFactorShiftType),arg1,arg2)
end

function PCFactorSetShiftAmount(arg1::PC,PetscReal::Cint)
    ccall((:PCFactorSetShiftAmount,petsc),PetscErrorCode,(PC,Cint),arg1,PetscReal)
end

function PCFactorSetMatSolverPackage(arg1::PC,arg2::Ptr{Uint8})
    ccall((:PCFactorSetMatSolverPackage,petsc),PetscErrorCode,(PC,Ptr{Uint8}),arg1,arg2)
end

function PCFactorGetMatSolverPackage(arg1::PC,arg2::Ptr{Ptr{Uint8}})
    ccall((:PCFactorGetMatSolverPackage,petsc),PetscErrorCode,(PC,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function PCFactorSetUpMatSolverPackage(arg1::PC)
    ccall((:PCFactorSetUpMatSolverPackage,petsc),PetscErrorCode,(PC,),arg1)
end

function PCFactorSetFill(arg1::PC,PetscReal::Cint)
    ccall((:PCFactorSetFill,petsc),PetscErrorCode,(PC,Cint),arg1,PetscReal)
end

function PCFactorSetColumnPivot(arg1::PC,PetscReal::Cint)
    ccall((:PCFactorSetColumnPivot,petsc),PetscErrorCode,(PC,Cint),arg1,PetscReal)
end

function PCFactorReorderForNonzeroDiagonal(arg1::PC,PetscReal::Cint)
    ccall((:PCFactorReorderForNonzeroDiagonal,petsc),PetscErrorCode,(PC,Cint),arg1,PetscReal)
end

function PCFactorSetMatOrderingType(arg1::PC,arg2::MatOrderingType)
    ccall((:PCFactorSetMatOrderingType,petsc),PetscErrorCode,(PC,MatOrderingType),arg1,arg2)
end

function PCFactorSetReuseOrdering(arg1::PC,arg2::PetscBool)
    ccall((:PCFactorSetReuseOrdering,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCFactorSetReuseFill(arg1::PC,arg2::PetscBool)
    ccall((:PCFactorSetReuseFill,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCFactorSetUseInPlace(arg1::PC,arg2::PetscBool)
    ccall((:PCFactorSetUseInPlace,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCFactorGetUseInPlace(arg1::PC,arg2::Ptr{PetscBool})
    ccall((:PCFactorGetUseInPlace,petsc),PetscErrorCode,(PC,Ptr{PetscBool}),arg1,arg2)
end

function PCFactorSetAllowDiagonalFill(arg1::PC,arg2::PetscBool)
    ccall((:PCFactorSetAllowDiagonalFill,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCFactorGetAllowDiagonalFill(arg1::PC,arg2::Ptr{PetscBool})
    ccall((:PCFactorGetAllowDiagonalFill,petsc),PetscErrorCode,(PC,Ptr{PetscBool}),arg1,arg2)
end

function PCFactorSetPivotInBlocks(arg1::PC,arg2::PetscBool)
    ccall((:PCFactorSetPivotInBlocks,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCFactorSetLevels(arg1::PC,arg2::PetscInt)
    ccall((:PCFactorSetLevels,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCFactorGetLevels(arg1::PC,arg2::Ptr{PetscInt})
    ccall((:PCFactorGetLevels,petsc),PetscErrorCode,(PC,Ptr{PetscInt}),arg1,arg2)
end

function PCFactorSetDropTolerance(arg1::PC,PetscReal::Cint,arg2::Cint,arg3::PetscInt)
    ccall((:PCFactorSetDropTolerance,petsc),PetscErrorCode,(PC,Cint,Cint,PetscInt),arg1,PetscReal,arg2,arg3)
end

function PCASMSetLocalSubdomains(arg1::PC,arg2::PetscInt,arg3::Ptr{IS},arg4::Ptr{IS})
    ccall((:PCASMSetLocalSubdomains,petsc),PetscErrorCode,(PC,PetscInt,Ptr{IS},Ptr{IS}),arg1,arg2,arg3,arg4)
end

function PCASMSetTotalSubdomains(arg1::PC,arg2::PetscInt,arg3::Ptr{IS},arg4::Ptr{IS})
    ccall((:PCASMSetTotalSubdomains,petsc),PetscErrorCode,(PC,PetscInt,Ptr{IS},Ptr{IS}),arg1,arg2,arg3,arg4)
end

function PCASMSetOverlap(arg1::PC,arg2::PetscInt)
    ccall((:PCASMSetOverlap,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCASMSetDMSubdomains(arg1::PC,arg2::PetscBool)
    ccall((:PCASMSetDMSubdomains,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCASMGetDMSubdomains(arg1::PC,arg2::Ptr{PetscBool})
    ccall((:PCASMGetDMSubdomains,petsc),PetscErrorCode,(PC,Ptr{PetscBool}),arg1,arg2)
end

function PCASMSetSortIndices(arg1::PC,arg2::PetscBool)
    ccall((:PCASMSetSortIndices,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCASMSetType(arg1::PC,arg2::PCASMType)
    ccall((:PCASMSetType,petsc),PetscErrorCode,(PC,PCASMType),arg1,arg2)
end

function PCASMGetType(arg1::PC,arg2::Ptr{PCASMType})
    ccall((:PCASMGetType,petsc),PetscErrorCode,(PC,Ptr{PCASMType}),arg1,arg2)
end

function PCASMSetLocalType(arg1::PC,arg2::PCCompositeType)
    ccall((:PCASMSetLocalType,petsc),PetscErrorCode,(PC,PCCompositeType),arg1,arg2)
end

function PCASMGetLocalType(arg1::PC,arg2::Ptr{PCCompositeType})
    ccall((:PCASMGetLocalType,petsc),PetscErrorCode,(PC,Ptr{PCCompositeType}),arg1,arg2)
end

function PCASMCreateSubdomains(arg1::Mat,arg2::PetscInt,arg3::Ptr{Ptr{IS}})
    ccall((:PCASMCreateSubdomains,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{Ptr{IS}}),arg1,arg2,arg3)
end

function PCASMDestroySubdomains(arg1::PetscInt,arg2::Ptr{IS},arg3::Ptr{IS})
    ccall((:PCASMDestroySubdomains,petsc),PetscErrorCode,(PetscInt,Ptr{IS},Ptr{IS}),arg1,arg2,arg3)
end

function PCASMCreateSubdomains2D(arg1::PetscInt,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::PetscInt,arg7::Ptr{PetscInt},arg8::Ptr{Ptr{IS}},arg9::Ptr{Ptr{IS}})
    ccall((:PCASMCreateSubdomains2D,petsc),PetscErrorCode,(PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Ptr{IS}},Ptr{Ptr{IS}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function PCASMGetLocalSubdomains(arg1::PC,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{IS}},arg4::Ptr{Ptr{IS}})
    ccall((:PCASMGetLocalSubdomains,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{Ptr{IS}},Ptr{Ptr{IS}}),arg1,arg2,arg3,arg4)
end

function PCASMGetLocalSubmatrices(arg1::PC,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{Mat}})
    ccall((:PCASMGetLocalSubmatrices,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{Ptr{Mat}}),arg1,arg2,arg3)
end

function PCGASMSetTotalSubdomains(arg1::PC,arg2::PetscInt)
    ccall((:PCGASMSetTotalSubdomains,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCGASMSetSubdomains(arg1::PC,arg2::PetscInt,arg3::Ptr{IS},arg4::Ptr{IS})
    ccall((:PCGASMSetSubdomains,petsc),PetscErrorCode,(PC,PetscInt,Ptr{IS},Ptr{IS}),arg1,arg2,arg3,arg4)
end

function PCGASMSetOverlap(arg1::PC,arg2::PetscInt)
    ccall((:PCGASMSetOverlap,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCGASMSetUseDMSubdomains(arg1::PC,arg2::PetscBool)
    ccall((:PCGASMSetUseDMSubdomains,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCGASMGetUseDMSubdomains(arg1::PC,arg2::Ptr{PetscBool})
    ccall((:PCGASMGetUseDMSubdomains,petsc),PetscErrorCode,(PC,Ptr{PetscBool}),arg1,arg2)
end

function PCGASMSetSortIndices(arg1::PC,arg2::PetscBool)
    ccall((:PCGASMSetSortIndices,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCGASMSetType(arg1::PC,arg2::PCGASMType)
    ccall((:PCGASMSetType,petsc),PetscErrorCode,(PC,PCGASMType),arg1,arg2)
end

function PCGASMCreateSubdomains(arg1::Mat,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{Ptr{IS}})
    ccall((:PCGASMCreateSubdomains,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{Ptr{IS}}),arg1,arg2,arg3,arg4)
end

function PCGASMDestroySubdomains(arg1::PetscInt,arg2::Ptr{Ptr{IS}},arg3::Ptr{Ptr{IS}})
    ccall((:PCGASMDestroySubdomains,petsc),PetscErrorCode,(PetscInt,Ptr{Ptr{IS}},Ptr{Ptr{IS}}),arg1,arg2,arg3)
end

function PCGASMCreateSubdomains2D(arg1::PC,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt,arg5::PetscInt,arg6::PetscInt,arg7::PetscInt,arg8::Ptr{PetscInt},arg9::Ptr{Ptr{IS}},arg10::Ptr{Ptr{IS}})
    ccall((:PCGASMCreateSubdomains2D,petsc),PetscErrorCode,(PC,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Ptr{IS}},Ptr{Ptr{IS}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function PCGASMGetSubdomains(arg1::PC,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{IS}},arg4::Ptr{Ptr{IS}})
    ccall((:PCGASMGetSubdomains,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{Ptr{IS}},Ptr{Ptr{IS}}),arg1,arg2,arg3,arg4)
end

function PCGASMGetSubmatrices(arg1::PC,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{Mat}})
    ccall((:PCGASMGetSubmatrices,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{Ptr{Mat}}),arg1,arg2,arg3)
end

function PCCompositeSetType(arg1::PC,arg2::PCCompositeType)
    ccall((:PCCompositeSetType,petsc),PetscErrorCode,(PC,PCCompositeType),arg1,arg2)
end

function PCCompositeGetType(arg1::PC,arg2::Ptr{PCCompositeType})
    ccall((:PCCompositeGetType,petsc),PetscErrorCode,(PC,Ptr{PCCompositeType}),arg1,arg2)
end

function PCCompositeAddPC(arg1::PC,arg2::PCType)
    ccall((:PCCompositeAddPC,petsc),PetscErrorCode,(PC,PCType),arg1,arg2)
end

function PCCompositeGetPC(arg1::PC,arg2::PetscInt,arg3::Ptr{PC})
    ccall((:PCCompositeGetPC,petsc),PetscErrorCode,(PC,PetscInt,Ptr{PC}),arg1,arg2,arg3)
end

function PCCompositeSpecialSetAlpha(arg1::PC,arg2::PetscScalar)
    ccall((:PCCompositeSpecialSetAlpha,petsc),PetscErrorCode,(PC,PetscScalar),arg1,arg2)
end

function PCRedundantSetNumber(arg1::PC,arg2::PetscInt)
    ccall((:PCRedundantSetNumber,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCRedundantSetScatter(arg1::PC,arg2::VecScatter,arg3::VecScatter)
    ccall((:PCRedundantSetScatter,petsc),PetscErrorCode,(PC,VecScatter,VecScatter),arg1,arg2,arg3)
end

function PCRedundantGetOperators(arg1::PC,arg2::Ptr{Mat},arg3::Ptr{Mat})
    ccall((:PCRedundantGetOperators,petsc),PetscErrorCode,(PC,Ptr{Mat},Ptr{Mat}),arg1,arg2,arg3)
end

function PCSPAISetEpsilon(arg1::PC,arg2::Cdouble)
    ccall((:PCSPAISetEpsilon,petsc),PetscErrorCode,(PC,Cdouble),arg1,arg2)
end

function PCSPAISetNBSteps(arg1::PC,arg2::PetscInt)
    ccall((:PCSPAISetNBSteps,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCSPAISetMax(arg1::PC,arg2::PetscInt)
    ccall((:PCSPAISetMax,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCSPAISetMaxNew(arg1::PC,arg2::PetscInt)
    ccall((:PCSPAISetMaxNew,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCSPAISetBlockSize(arg1::PC,arg2::PetscInt)
    ccall((:PCSPAISetBlockSize,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCSPAISetCacheSize(arg1::PC,arg2::PetscInt)
    ccall((:PCSPAISetCacheSize,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCSPAISetVerbose(arg1::PC,arg2::PetscInt)
    ccall((:PCSPAISetVerbose,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCSPAISetSp(arg1::PC,arg2::PetscInt)
    ccall((:PCSPAISetSp,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCHYPRESetType(arg1::PC,arg2::Ptr{Uint8})
    ccall((:PCHYPRESetType,petsc),PetscErrorCode,(PC,Ptr{Uint8}),arg1,arg2)
end

function PCHYPREGetType(arg1::PC,arg2::Ptr{Ptr{Uint8}})
    ccall((:PCHYPREGetType,petsc),PetscErrorCode,(PC,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function PCHYPRESetDiscreteGradient(arg1::PC,arg2::Mat)
    ccall((:PCHYPRESetDiscreteGradient,petsc),PetscErrorCode,(PC,Mat),arg1,arg2)
end

function PCHYPRESetDiscreteCurl(arg1::PC,arg2::Mat)
    ccall((:PCHYPRESetDiscreteCurl,petsc),PetscErrorCode,(PC,Mat),arg1,arg2)
end

function PCHYPRESetEdgeConstantVectors(arg1::PC,arg2::Vec,arg3::Vec,arg4::Vec)
    ccall((:PCHYPRESetEdgeConstantVectors,petsc),PetscErrorCode,(PC,Vec,Vec,Vec),arg1,arg2,arg3,arg4)
end

function PCHYPRESetAlphaPoissonMatrix(arg1::PC,arg2::Mat)
    ccall((:PCHYPRESetAlphaPoissonMatrix,petsc),PetscErrorCode,(PC,Mat),arg1,arg2)
end

function PCHYPRESetBetaPoissonMatrix(arg1::PC,arg2::Mat)
    ccall((:PCHYPRESetBetaPoissonMatrix,petsc),PetscErrorCode,(PC,Mat),arg1,arg2)
end

function PCBJacobiGetLocalBlocks(arg1::PC,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{PetscInt}})
    ccall((:PCBJacobiGetLocalBlocks,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3)
end

function PCBJacobiGetTotalBlocks(arg1::PC,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{PetscInt}})
    ccall((:PCBJacobiGetTotalBlocks,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{Ptr{PetscInt}}),arg1,arg2,arg3)
end

function PCFieldSplitSetFields(arg1::PC,arg2::Ptr{Uint8},arg3::PetscInt,arg4::Ptr{PetscInt},arg5::Ptr{PetscInt})
    ccall((:PCFieldSplitSetFields,petsc),PetscErrorCode,(PC,Ptr{Uint8},PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end

function PCFieldSplitSetType(arg1::PC,arg2::PCCompositeType)
    ccall((:PCFieldSplitSetType,petsc),PetscErrorCode,(PC,PCCompositeType),arg1,arg2)
end

function PCFieldSplitGetType(arg1::PC,arg2::Ptr{PCCompositeType})
    ccall((:PCFieldSplitGetType,petsc),PetscErrorCode,(PC,Ptr{PCCompositeType}),arg1,arg2)
end

function PCFieldSplitSetBlockSize(arg1::PC,arg2::PetscInt)
    ccall((:PCFieldSplitSetBlockSize,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCFieldSplitSetIS(arg1::PC,arg2::Ptr{Uint8},arg3::IS)
    ccall((:PCFieldSplitSetIS,petsc),PetscErrorCode,(PC,Ptr{Uint8},IS),arg1,arg2,arg3)
end

function PCFieldSplitGetIS(arg1::PC,arg2::Ptr{Uint8},arg3::Ptr{IS})
    ccall((:PCFieldSplitGetIS,petsc),PetscErrorCode,(PC,Ptr{Uint8},Ptr{IS}),arg1,arg2,arg3)
end

function PCFieldSplitSetDMSplits(arg1::PC,arg2::PetscBool)
    ccall((:PCFieldSplitSetDMSplits,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCFieldSplitGetDMSplits(arg1::PC,arg2::Ptr{PetscBool})
    ccall((:PCFieldSplitGetDMSplits,petsc),PetscErrorCode,(PC,Ptr{PetscBool}),arg1,arg2)
end

function PCFieldSplitSetDiagUseAmat(arg1::PC,arg2::PetscBool)
    ccall((:PCFieldSplitSetDiagUseAmat,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCFieldSplitGetDiagUseAmat(arg1::PC,arg2::Ptr{PetscBool})
    ccall((:PCFieldSplitGetDiagUseAmat,petsc),PetscErrorCode,(PC,Ptr{PetscBool}),arg1,arg2)
end

function PCFieldSplitSetOffDiagUseAmat(arg1::PC,arg2::PetscBool)
    ccall((:PCFieldSplitSetOffDiagUseAmat,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCFieldSplitGetOffDiagUseAmat(arg1::PC,arg2::Ptr{PetscBool})
    ccall((:PCFieldSplitGetOffDiagUseAmat,petsc),PetscErrorCode,(PC,Ptr{PetscBool}),arg1,arg2)
end

function PETSC_DEPRECATED()
    ccall((:PETSC_DEPRECATED,petsc),Cint,())
end

function PCFieldSplitSchurPrecondition(arg1::PC,arg2::PCFieldSplitSchurPreType,arg3::Mat)
    ccall((:PCFieldSplitSchurPrecondition,petsc),PetscErrorCode,(PC,PCFieldSplitSchurPreType,Mat),arg1,arg2,arg3)
end

function PCFieldSplitSetSchurPre(arg1::PC,arg2::PCFieldSplitSchurPreType,arg3::Mat)
    ccall((:PCFieldSplitSetSchurPre,petsc),PetscErrorCode,(PC,PCFieldSplitSchurPreType,Mat),arg1,arg2,arg3)
end

function PCFieldSplitGetSchurPre(arg1::PC,arg2::Ptr{PCFieldSplitSchurPreType},arg3::Ptr{Mat})
    ccall((:PCFieldSplitGetSchurPre,petsc),PetscErrorCode,(PC,Ptr{PCFieldSplitSchurPreType},Ptr{Mat}),arg1,arg2,arg3)
end

function PCFieldSplitSetSchurFactType(arg1::PC,arg2::PCFieldSplitSchurFactType)
    ccall((:PCFieldSplitSetSchurFactType,petsc),PetscErrorCode,(PC,PCFieldSplitSchurFactType),arg1,arg2)
end

function PCFieldSplitGetSchurBlocks(arg1::PC,arg2::Ptr{Mat},arg3::Ptr{Mat},arg4::Ptr{Mat},arg5::Ptr{Mat})
    ccall((:PCFieldSplitGetSchurBlocks,petsc),PetscErrorCode,(PC,Ptr{Mat},Ptr{Mat},Ptr{Mat},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5)
end

function PCFieldSplitSchurGetS(arg1::PC,S::Ptr{Mat})
    ccall((:PCFieldSplitSchurGetS,petsc),PetscErrorCode,(PC,Ptr{Mat}),arg1,S)
end

function PCFieldSplitSchurRestoreS(arg1::PC,S::Ptr{Mat})
    ccall((:PCFieldSplitSchurRestoreS,petsc),PetscErrorCode,(PC,Ptr{Mat}),arg1,S)
end

function PCGalerkinSetRestriction(arg1::PC,arg2::Mat)
    ccall((:PCGalerkinSetRestriction,petsc),PetscErrorCode,(PC,Mat),arg1,arg2)
end

function PCGalerkinSetInterpolation(arg1::PC,arg2::Mat)
    ccall((:PCGalerkinSetInterpolation,petsc),PetscErrorCode,(PC,Mat),arg1,arg2)
end

function PCSetCoordinates(arg1::PC,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{Cint})
    ccall((:PCSetCoordinates,petsc),PetscErrorCode,(PC,PetscInt,PetscInt,Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function PCPythonSetType(arg1::PC,arg2::Ptr{Uint8})
    ccall((:PCPythonSetType,petsc),PetscErrorCode,(PC,Ptr{Uint8}),arg1,arg2)
end

function PCSetDM(arg1::PC,arg2::DM)
    ccall((:PCSetDM,petsc),PetscErrorCode,(PC,DM),arg1,arg2)
end

function PCGetDM(arg1::PC,arg2::Ptr{DM})
    ccall((:PCGetDM,petsc),PetscErrorCode,(PC,Ptr{DM}),arg1,arg2)
end

function PCSetApplicationContext(arg1::PC,arg2::Ptr{Void})
    ccall((:PCSetApplicationContext,petsc),PetscErrorCode,(PC,Ptr{Void}),arg1,arg2)
end

function PCGetApplicationContext(arg1::PC,arg2::Ptr{Void})
    ccall((:PCGetApplicationContext,petsc),PetscErrorCode,(PC,Ptr{Void}),arg1,arg2)
end

function PCBiCGStabCUSPSetTolerance(arg1::PC,PetscReal::Cint)
    ccall((:PCBiCGStabCUSPSetTolerance,petsc),PetscErrorCode,(PC,Cint),arg1,PetscReal)
end

function PCBiCGStabCUSPSetIterations(arg1::PC,arg2::PetscInt)
    ccall((:PCBiCGStabCUSPSetIterations,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCBiCGStabCUSPSetUseVerboseMonitor(arg1::PC,arg2::PetscBool)
    ccall((:PCBiCGStabCUSPSetUseVerboseMonitor,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCAINVCUSPSetDropTolerance(arg1::PC,PetscReal::Cint)
    ccall((:PCAINVCUSPSetDropTolerance,petsc),PetscErrorCode,(PC,Cint),arg1,PetscReal)
end

function PCAINVCUSPUseScaling(arg1::PC,arg2::PetscBool)
    ccall((:PCAINVCUSPUseScaling,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCAINVCUSPSetNonzeros(arg1::PC,arg2::PetscInt)
    ccall((:PCAINVCUSPSetNonzeros,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCAINVCUSPSetLinParameter(arg1::PC,arg2::PetscInt)
    ccall((:PCAINVCUSPSetLinParameter,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCPARMSSetGlobal(arg1::PC,arg2::PCPARMSGlobalType)
    ccall((:PCPARMSSetGlobal,petsc),PetscErrorCode,(PC,PCPARMSGlobalType),arg1,arg2)
end

function PCPARMSSetLocal(arg1::PC,arg2::PCPARMSLocalType)
    ccall((:PCPARMSSetLocal,petsc),PetscErrorCode,(PC,PCPARMSLocalType),arg1,arg2)
end

function PCPARMSSetSolveTolerances(arg1::PC,PetscReal::Cint,arg2::PetscInt)
    ccall((:PCPARMSSetSolveTolerances,petsc),PetscErrorCode,(PC,Cint,PetscInt),arg1,PetscReal,arg2)
end

function PCPARMSSetSolveRestart(arg1::PC,arg2::PetscInt)
    ccall((:PCPARMSSetSolveRestart,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCPARMSSetNonsymPerm(arg1::PC,arg2::PetscBool)
    ccall((:PCPARMSSetNonsymPerm,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCPARMSSetFill(arg1::PC,arg2::PetscInt,arg3::PetscInt,arg4::PetscInt)
    ccall((:PCPARMSSetFill,petsc),PetscErrorCode,(PC,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end

function PCGAMGSetType(arg1::PC,arg2::PCGAMGType)
    ccall((:PCGAMGSetType,petsc),PetscErrorCode,(PC,PCGAMGType),arg1,arg2)
end

function PCGAMGGetType(arg1::PC,arg2::Ptr{PCGAMGType})
    ccall((:PCGAMGGetType,petsc),PetscErrorCode,(PC,Ptr{PCGAMGType}),arg1,arg2)
end

function PCGAMGSetProcEqLim(arg1::PC,arg2::PetscInt)
    ccall((:PCGAMGSetProcEqLim,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCGAMGSetRepartitioning(arg1::PC,arg2::PetscBool)
    ccall((:PCGAMGSetRepartitioning,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCGAMGSetUseASMAggs(arg1::PC,arg2::PetscBool)
    ccall((:PCGAMGSetUseASMAggs,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCGAMGSetSolverType(arg1::PC,arg2::Ptr{Uint8},arg3::PetscInt)
    ccall((:PCGAMGSetSolverType,petsc),PetscErrorCode,(PC,Ptr{Uint8},PetscInt),arg1,arg2,arg3)
end

function PCGAMGSetThreshold(arg1::PC,PetscReal::Cint)
    ccall((:PCGAMGSetThreshold,petsc),PetscErrorCode,(PC,Cint),arg1,PetscReal)
end

function PCGAMGSetCoarseEqLim(arg1::PC,arg2::PetscInt)
    ccall((:PCGAMGSetCoarseEqLim,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCGAMGSetNlevels(arg1::PC,arg2::PetscInt)
    ccall((:PCGAMGSetNlevels,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCGAMGSetNSmooths(arg1::PC,arg2::PetscInt)
    ccall((:PCGAMGSetNSmooths,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCGAMGSetSymGraph(arg1::PC,arg2::PetscBool)
    ccall((:PCGAMGSetSymGraph,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCGAMGSetSquareGraph(arg1::PC,arg2::PetscInt)
    ccall((:PCGAMGSetSquareGraph,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCGAMGSetReuseInterpolation(arg1::PC,arg2::PetscBool)
    ccall((:PCGAMGSetReuseInterpolation,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCGAMGFinalizePackage()
    ccall((:PCGAMGFinalizePackage,petsc),PetscErrorCode,())
end

function PCGAMGInitializePackage()
    ccall((:PCGAMGInitializePackage,petsc),PetscErrorCode,())
end

function PCGAMGRegister(arg1::PCGAMGType,arg2::Ptr{Void})
    ccall((:PCGAMGRegister,petsc),PetscErrorCode,(PCGAMGType,Ptr{Void}),arg1,arg2)
end

function PCGAMGClassicalSetType(arg1::PC,arg2::PCGAMGClassicalType)
    ccall((:PCGAMGClassicalSetType,petsc),PetscErrorCode,(PC,PCGAMGClassicalType),arg1,arg2)
end

function PCGAMGClassicalGetType(arg1::PC,arg2::Ptr{PCGAMGClassicalType})
    ccall((:PCGAMGClassicalGetType,petsc),PetscErrorCode,(PC,Ptr{PCGAMGClassicalType}),arg1,arg2)
end

function PCBDDCSetChangeOfBasisMat(arg1::PC,arg2::Mat)
    ccall((:PCBDDCSetChangeOfBasisMat,petsc),PetscErrorCode,(PC,Mat),arg1,arg2)
end

function PCBDDCSetPrimalVerticesLocalIS(arg1::PC,arg2::IS)
    ccall((:PCBDDCSetPrimalVerticesLocalIS,petsc),PetscErrorCode,(PC,IS),arg1,arg2)
end

function PCBDDCSetCoarseningRatio(arg1::PC,arg2::PetscInt)
    ccall((:PCBDDCSetCoarseningRatio,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCBDDCSetLevels(arg1::PC,arg2::PetscInt)
    ccall((:PCBDDCSetLevels,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCBDDCSetNullSpace(arg1::PC,arg2::MatNullSpace)
    ccall((:PCBDDCSetNullSpace,petsc),PetscErrorCode,(PC,MatNullSpace),arg1,arg2)
end

function PCBDDCSetDirichletBoundaries(arg1::PC,arg2::IS)
    ccall((:PCBDDCSetDirichletBoundaries,petsc),PetscErrorCode,(PC,IS),arg1,arg2)
end

function PCBDDCSetDirichletBoundariesLocal(arg1::PC,arg2::IS)
    ccall((:PCBDDCSetDirichletBoundariesLocal,petsc),PetscErrorCode,(PC,IS),arg1,arg2)
end

function PCBDDCGetDirichletBoundaries(arg1::PC,arg2::Ptr{IS})
    ccall((:PCBDDCGetDirichletBoundaries,petsc),PetscErrorCode,(PC,Ptr{IS}),arg1,arg2)
end

function PCBDDCGetDirichletBoundariesLocal(arg1::PC,arg2::Ptr{IS})
    ccall((:PCBDDCGetDirichletBoundariesLocal,petsc),PetscErrorCode,(PC,Ptr{IS}),arg1,arg2)
end

function PCBDDCSetNeumannBoundaries(arg1::PC,arg2::IS)
    ccall((:PCBDDCSetNeumannBoundaries,petsc),PetscErrorCode,(PC,IS),arg1,arg2)
end

function PCBDDCSetNeumannBoundariesLocal(arg1::PC,arg2::IS)
    ccall((:PCBDDCSetNeumannBoundariesLocal,petsc),PetscErrorCode,(PC,IS),arg1,arg2)
end

function PCBDDCGetNeumannBoundaries(arg1::PC,arg2::Ptr{IS})
    ccall((:PCBDDCGetNeumannBoundaries,petsc),PetscErrorCode,(PC,Ptr{IS}),arg1,arg2)
end

function PCBDDCGetNeumannBoundariesLocal(arg1::PC,arg2::Ptr{IS})
    ccall((:PCBDDCGetNeumannBoundariesLocal,petsc),PetscErrorCode,(PC,Ptr{IS}),arg1,arg2)
end

function PCBDDCSetDofsSplitting(arg1::PC,arg2::PetscInt,arg3::Ptr{IS})
    ccall((:PCBDDCSetDofsSplitting,petsc),PetscErrorCode,(PC,PetscInt,Ptr{IS}),arg1,arg2,arg3)
end

function PCBDDCSetDofsSplittingLocal(arg1::PC,arg2::PetscInt,arg3::Ptr{IS})
    ccall((:PCBDDCSetDofsSplittingLocal,petsc),PetscErrorCode,(PC,PetscInt,Ptr{IS}),arg1,arg2,arg3)
end

function PCBDDCSetLocalAdjacencyGraph(arg1::PC,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscInt},PetscCopyMode::Cint)
    ccall((:PCBDDCSetLocalAdjacencyGraph,petsc),PetscErrorCode,(PC,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Cint),arg1,arg2,arg3,arg4,PetscCopyMode)
end

function PCBDDCCreateFETIDPOperators(arg1::PC,arg2::Ptr{Mat},arg3::Ptr{PC})
    ccall((:PCBDDCCreateFETIDPOperators,petsc),PetscErrorCode,(PC,Ptr{Mat},Ptr{PC}),arg1,arg2,arg3)
end

function PCBDDCMatFETIDPGetRHS(arg1::Mat,arg2::Vec,arg3::Vec)
    ccall((:PCBDDCMatFETIDPGetRHS,petsc),PetscErrorCode,(Mat,Vec,Vec),arg1,arg2,arg3)
end

function PCBDDCMatFETIDPGetSolution(arg1::Mat,arg2::Vec,arg3::Vec)
    ccall((:PCBDDCMatFETIDPGetSolution,petsc),PetscErrorCode,(Mat,Vec,Vec),arg1,arg2,arg3)
end

function PCISSetUseStiffnessScaling(arg1::PC,arg2::PetscBool)
    ccall((:PCISSetUseStiffnessScaling,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCISSetSubdomainScalingFactor(arg1::PC,arg2::PetscScalar)
    ccall((:PCISSetSubdomainScalingFactor,petsc),PetscErrorCode,(PC,PetscScalar),arg1,arg2)
end

function PCISSetSubdomainDiagonalScaling(arg1::PC,arg2::Vec)
    ccall((:PCISSetSubdomainDiagonalScaling,petsc),PetscErrorCode,(PC,Vec),arg1,arg2)
end

function PCMGSetType(arg1::PC,arg2::PCMGType)
    ccall((:PCMGSetType,petsc),PetscErrorCode,(PC,PCMGType),arg1,arg2)
end

function PCMGGetType(arg1::PC,arg2::Ptr{PCMGType})
    ccall((:PCMGGetType,petsc),PetscErrorCode,(PC,Ptr{PCMGType}),arg1,arg2)
end

function PCMGSetLevels(arg1::PC,arg2::PetscInt,arg3::Ptr{MPI_Comm})
    ccall((:PCMGSetLevels,petsc),PetscErrorCode,(PC,PetscInt,Ptr{MPI_Comm}),arg1,arg2,arg3)
end

function PCMGGetLevels(arg1::PC,arg2::Ptr{PetscInt})
    ccall((:PCMGGetLevels,petsc),PetscErrorCode,(PC,Ptr{PetscInt}),arg1,arg2)
end

function PCMGSetNumberSmoothUp(arg1::PC,arg2::PetscInt)
    ccall((:PCMGSetNumberSmoothUp,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCMGSetNumberSmoothDown(arg1::PC,arg2::PetscInt)
    ccall((:PCMGSetNumberSmoothDown,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCMGSetCycleType(arg1::PC,arg2::PCMGCycleType)
    ccall((:PCMGSetCycleType,petsc),PetscErrorCode,(PC,PCMGCycleType),arg1,arg2)
end

function PCMGSetCycleTypeOnLevel(arg1::PC,arg2::PetscInt,arg3::PCMGCycleType)
    ccall((:PCMGSetCycleTypeOnLevel,petsc),PetscErrorCode,(PC,PetscInt,PCMGCycleType),arg1,arg2,arg3)
end

function PCMGSetCyclesOnLevel(arg1::PC,arg2::PetscInt,arg3::PetscInt)
    ccall((:PCMGSetCyclesOnLevel,petsc),PetscErrorCode,(PC,PetscInt,PetscInt),arg1,arg2,arg3)
end

function PCMGMultiplicativeSetCycles(arg1::PC,arg2::PetscInt)
    ccall((:PCMGMultiplicativeSetCycles,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCMGSetGalerkin(arg1::PC,arg2::PetscBool)
    ccall((:PCMGSetGalerkin,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCMGGetGalerkin(arg1::PC,arg2::Ptr{PetscBool})
    ccall((:PCMGGetGalerkin,petsc),PetscErrorCode,(PC,Ptr{PetscBool}),arg1,arg2)
end

function PCMGSetRhs(arg1::PC,arg2::PetscInt,arg3::Vec)
    ccall((:PCMGSetRhs,petsc),PetscErrorCode,(PC,PetscInt,Vec),arg1,arg2,arg3)
end

function PCMGSetX(arg1::PC,arg2::PetscInt,arg3::Vec)
    ccall((:PCMGSetX,petsc),PetscErrorCode,(PC,PetscInt,Vec),arg1,arg2,arg3)
end

function PCMGSetR(arg1::PC,arg2::PetscInt,arg3::Vec)
    ccall((:PCMGSetR,petsc),PetscErrorCode,(PC,PetscInt,Vec),arg1,arg2,arg3)
end

function PCMGSetRestriction(arg1::PC,arg2::PetscInt,arg3::Mat)
    ccall((:PCMGSetRestriction,petsc),PetscErrorCode,(PC,PetscInt,Mat),arg1,arg2,arg3)
end

function PCMGGetRestriction(arg1::PC,arg2::PetscInt,arg3::Ptr{Mat})
    ccall((:PCMGGetRestriction,petsc),PetscErrorCode,(PC,PetscInt,Ptr{Mat}),arg1,arg2,arg3)
end

function PCMGSetInterpolation(arg1::PC,arg2::PetscInt,arg3::Mat)
    ccall((:PCMGSetInterpolation,petsc),PetscErrorCode,(PC,PetscInt,Mat),arg1,arg2,arg3)
end

function PCMGGetInterpolation(arg1::PC,arg2::PetscInt,arg3::Ptr{Mat})
    ccall((:PCMGGetInterpolation,petsc),PetscErrorCode,(PC,PetscInt,Ptr{Mat}),arg1,arg2,arg3)
end

function PCMGSetRScale(arg1::PC,arg2::PetscInt,arg3::Vec)
    ccall((:PCMGSetRScale,petsc),PetscErrorCode,(PC,PetscInt,Vec),arg1,arg2,arg3)
end

function PCMGGetRScale(arg1::PC,arg2::PetscInt,arg3::Ptr{Vec})
    ccall((:PCMGGetRScale,petsc),PetscErrorCode,(PC,PetscInt,Ptr{Vec}),arg1,arg2,arg3)
end

function PCMGSetResidual(arg1::PC,arg2::PetscInt,arg3::Ptr{Void},arg4::Mat)
    ccall((:PCMGSetResidual,petsc),PetscErrorCode,(PC,PetscInt,Ptr{Void},Mat),arg1,arg2,arg3,arg4)
end

function PCMGResidualDefault(arg1::Mat,arg2::Vec,arg3::Vec,arg4::Vec)
    ccall((:PCMGResidualDefault,petsc),PetscErrorCode,(Mat,Vec,Vec,Vec),arg1,arg2,arg3,arg4)
end

function KSPInitializePackage()
    ccall((:KSPInitializePackage,petsc),PetscErrorCode,())
end

function KSPCreate(arg1::MPI_Comm,arg2::Ptr{KSP})
    ccall((:KSPCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{KSP}),arg1,arg2)
end

function KSPSetType(arg1::KSP,arg2::KSPType)
    ccall((:KSPSetType,petsc),PetscErrorCode,(KSP,KSPType),arg1,arg2)
end

function KSPGetType(arg1::KSP,arg2::Ptr{KSPType})
    ccall((:KSPGetType,petsc),PetscErrorCode,(KSP,Ptr{KSPType}),arg1,arg2)
end

function KSPSetUp(arg1::KSP)
    ccall((:KSPSetUp,petsc),PetscErrorCode,(KSP,),arg1)
end

function KSPSetUpOnBlocks(arg1::KSP)
    ccall((:KSPSetUpOnBlocks,petsc),PetscErrorCode,(KSP,),arg1)
end

function KSPSolve(arg1::KSP,arg2::Vec,arg3::Vec)
    ccall((:KSPSolve,petsc),PetscErrorCode,(KSP,Vec,Vec),arg1,arg2,arg3)
end

function KSPSolveTranspose(arg1::KSP,arg2::Vec,arg3::Vec)
    ccall((:KSPSolveTranspose,petsc),PetscErrorCode,(KSP,Vec,Vec),arg1,arg2,arg3)
end

function KSPReset(arg1::KSP)
    ccall((:KSPReset,petsc),PetscErrorCode,(KSP,),arg1)
end

function KSPDestroy(arg1::Ptr{KSP})
    ccall((:KSPDestroy,petsc),PetscErrorCode,(Ptr{KSP},),arg1)
end

function KSPSetReusePreconditioner(arg1::KSP,arg2::PetscBool)
    ccall((:KSPSetReusePreconditioner,petsc),PetscErrorCode,(KSP,PetscBool),arg1,arg2)
end

function KSPSetSkipPCSetFromOptions(arg1::KSP,arg2::PetscBool)
    ccall((:KSPSetSkipPCSetFromOptions,petsc),PetscErrorCode,(KSP,PetscBool),arg1,arg2)
end

function KSPRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:KSPRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function KSPSetPCSide(arg1::KSP,arg2::PCSide)
    ccall((:KSPSetPCSide,petsc),PetscErrorCode,(KSP,PCSide),arg1,arg2)
end

function KSPGetPCSide(arg1::KSP,arg2::Ptr{PCSide})
    ccall((:KSPGetPCSide,petsc),PetscErrorCode,(KSP,Ptr{PCSide}),arg1,arg2)
end

function KSPSetTolerances(arg1::KSP,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::PetscInt)
    ccall((:KSPSetTolerances,petsc),PetscErrorCode,(KSP,Cint,Cint,Cint,PetscInt),arg1,PetscReal,arg2,arg3,arg4)
end

function KSPGetTolerances(arg1::KSP,arg2::Ptr{Cint},arg3::Ptr{Cint},arg4::Ptr{Cint},arg5::Ptr{PetscInt})
    ccall((:KSPGetTolerances,petsc),PetscErrorCode,(KSP,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end

function KSPSetInitialGuessNonzero(arg1::KSP,arg2::PetscBool)
    ccall((:KSPSetInitialGuessNonzero,petsc),PetscErrorCode,(KSP,PetscBool),arg1,arg2)
end

function KSPGetInitialGuessNonzero(arg1::KSP,arg2::Ptr{PetscBool})
    ccall((:KSPGetInitialGuessNonzero,petsc),PetscErrorCode,(KSP,Ptr{PetscBool}),arg1,arg2)
end

function KSPSetInitialGuessKnoll(arg1::KSP,arg2::PetscBool)
    ccall((:KSPSetInitialGuessKnoll,petsc),PetscErrorCode,(KSP,PetscBool),arg1,arg2)
end

function KSPGetInitialGuessKnoll(arg1::KSP,arg2::Ptr{PetscBool})
    ccall((:KSPGetInitialGuessKnoll,petsc),PetscErrorCode,(KSP,Ptr{PetscBool}),arg1,arg2)
end

function KSPSetErrorIfNotConverged(arg1::KSP,arg2::PetscBool)
    ccall((:KSPSetErrorIfNotConverged,petsc),PetscErrorCode,(KSP,PetscBool),arg1,arg2)
end

function KSPGetErrorIfNotConverged(arg1::KSP,arg2::Ptr{PetscBool})
    ccall((:KSPGetErrorIfNotConverged,petsc),PetscErrorCode,(KSP,Ptr{PetscBool}),arg1,arg2)
end

function KSPSetComputeEigenvalues(arg1::KSP,arg2::PetscBool)
    ccall((:KSPSetComputeEigenvalues,petsc),PetscErrorCode,(KSP,PetscBool),arg1,arg2)
end

function KSPGetComputeEigenvalues(arg1::KSP,arg2::Ptr{PetscBool})
    ccall((:KSPGetComputeEigenvalues,petsc),PetscErrorCode,(KSP,Ptr{PetscBool}),arg1,arg2)
end

function KSPSetComputeSingularValues(arg1::KSP,arg2::PetscBool)
    ccall((:KSPSetComputeSingularValues,petsc),PetscErrorCode,(KSP,PetscBool),arg1,arg2)
end

function KSPGetComputeSingularValues(arg1::KSP,arg2::Ptr{PetscBool})
    ccall((:KSPGetComputeSingularValues,petsc),PetscErrorCode,(KSP,Ptr{PetscBool}),arg1,arg2)
end

function KSPGetRhs(arg1::KSP,arg2::Ptr{Vec})
    ccall((:KSPGetRhs,petsc),PetscErrorCode,(KSP,Ptr{Vec}),arg1,arg2)
end

function KSPGetSolution(arg1::KSP,arg2::Ptr{Vec})
    ccall((:KSPGetSolution,petsc),PetscErrorCode,(KSP,Ptr{Vec}),arg1,arg2)
end

function KSPGetResidualNorm(arg1::KSP,arg2::Ptr{Cint})
    ccall((:KSPGetResidualNorm,petsc),PetscErrorCode,(KSP,Ptr{Cint}),arg1,arg2)
end

function KSPGetIterationNumber(arg1::KSP,arg2::Ptr{PetscInt})
    ccall((:KSPGetIterationNumber,petsc),PetscErrorCode,(KSP,Ptr{PetscInt}),arg1,arg2)
end

function KSPGetTotalIterations(arg1::KSP,arg2::Ptr{PetscInt})
    ccall((:KSPGetTotalIterations,petsc),PetscErrorCode,(KSP,Ptr{PetscInt}),arg1,arg2)
end

function KSPCreateVecs(arg1::KSP,arg2::PetscInt,arg3::Ptr{Ptr{Vec}},arg4::PetscInt,arg5::Ptr{Ptr{Vec}})
    ccall((:KSPCreateVecs,petsc),PetscErrorCode,(KSP,PetscInt,Ptr{Ptr{Vec}},PetscInt,Ptr{Ptr{Vec}}),arg1,arg2,arg3,arg4,arg5)
end

function KSPSetPostSolve(arg1::KSP,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:KSPSetPostSolve,petsc),PetscErrorCode,(KSP,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function KSPSetPC(arg1::KSP,arg2::PC)
    ccall((:KSPSetPC,petsc),PetscErrorCode,(KSP,PC),arg1,arg2)
end

function KSPGetPC(arg1::KSP,arg2::Ptr{PC})
    ccall((:KSPGetPC,petsc),PetscErrorCode,(KSP,Ptr{PC}),arg1,arg2)
end

function KSPMonitor(arg1::KSP,arg2::PetscInt,PetscReal::Cint)
    ccall((:KSPMonitor,petsc),PetscErrorCode,(KSP,PetscInt,Cint),arg1,arg2,PetscReal)
end

function KSPMonitorSet(arg1::KSP,arg2::Ptr{Void},arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:KSPMonitorSet,petsc),PetscErrorCode,(KSP,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function KSPMonitorCancel(arg1::KSP)
    ccall((:KSPMonitorCancel,petsc),PetscErrorCode,(KSP,),arg1)
end

function KSPGetMonitorContext(arg1::KSP,arg2::Ptr{Ptr{Void}})
    ccall((:KSPGetMonitorContext,petsc),PetscErrorCode,(KSP,Ptr{Ptr{Void}}),arg1,arg2)
end

function KSPGetResidualHistory(arg1::KSP,arg2::Ptr{Ptr{Cint}},arg3::Ptr{PetscInt})
    ccall((:KSPGetResidualHistory,petsc),PetscErrorCode,(KSP,Ptr{Ptr{Cint}},Ptr{PetscInt}),arg1,arg2,arg3)
end

function KSPSetResidualHistory(arg1::KSP,PetscReal::Ptr{Cint},arg2::PetscInt,arg3::PetscBool)
    ccall((:KSPSetResidualHistory,petsc),PetscErrorCode,(KSP,Ptr{Cint},PetscInt,PetscBool),arg1,PetscReal,arg2,arg3)
end

function KSPBuildSolutionDefault(arg1::KSP,arg2::Vec,arg3::Ptr{Vec})
    ccall((:KSPBuildSolutionDefault,petsc),PetscErrorCode,(KSP,Vec,Ptr{Vec}),arg1,arg2,arg3)
end

function KSPBuildResidualDefault(arg1::KSP,arg2::Vec,arg3::Vec,arg4::Ptr{Vec})
    ccall((:KSPBuildResidualDefault,petsc),PetscErrorCode,(KSP,Vec,Vec,Ptr{Vec}),arg1,arg2,arg3,arg4)
end

function KSPDestroyDefault(arg1::KSP)
    ccall((:KSPDestroyDefault,petsc),PetscErrorCode,(KSP,),arg1)
end

function KSPSetWorkVecs(arg1::KSP,arg2::PetscInt)
    ccall((:KSPSetWorkVecs,petsc),PetscErrorCode,(KSP,PetscInt),arg1,arg2)
end

function PCKSPGetKSP(arg1::PC,arg2::Ptr{KSP})
    ccall((:PCKSPGetKSP,petsc),PetscErrorCode,(PC,Ptr{KSP}),arg1,arg2)
end

function PCBJacobiGetSubKSP(arg1::PC,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{Ptr{KSP}})
    ccall((:PCBJacobiGetSubKSP,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{KSP}}),arg1,arg2,arg3,arg4)
end

function PCASMGetSubKSP(arg1::PC,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{Ptr{KSP}})
    ccall((:PCASMGetSubKSP,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{KSP}}),arg1,arg2,arg3,arg4)
end

function PCGASMGetSubKSP(arg1::PC,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{Ptr{KSP}})
    ccall((:PCGASMGetSubKSP,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{KSP}}),arg1,arg2,arg3,arg4)
end

function PCFieldSplitGetSubKSP(arg1::PC,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{KSP}})
    ccall((:PCFieldSplitGetSubKSP,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{Ptr{KSP}}),arg1,arg2,arg3)
end

function PCMGGetSmoother(arg1::PC,arg2::PetscInt,arg3::Ptr{KSP})
    ccall((:PCMGGetSmoother,petsc),PetscErrorCode,(PC,PetscInt,Ptr{KSP}),arg1,arg2,arg3)
end

function PCMGGetSmootherDown(arg1::PC,arg2::PetscInt,arg3::Ptr{KSP})
    ccall((:PCMGGetSmootherDown,petsc),PetscErrorCode,(PC,PetscInt,Ptr{KSP}),arg1,arg2,arg3)
end

function PCMGGetSmootherUp(arg1::PC,arg2::PetscInt,arg3::Ptr{KSP})
    ccall((:PCMGGetSmootherUp,petsc),PetscErrorCode,(PC,PetscInt,Ptr{KSP}),arg1,arg2,arg3)
end

function PCMGGetCoarseSolve(arg1::PC,arg2::Ptr{KSP})
    ccall((:PCMGGetCoarseSolve,petsc),PetscErrorCode,(PC,Ptr{KSP}),arg1,arg2)
end

function PCGalerkinGetKSP(arg1::PC,arg2::Ptr{KSP})
    ccall((:PCGalerkinGetKSP,petsc),PetscErrorCode,(PC,Ptr{KSP}),arg1,arg2)
end

function KSPBuildSolution(arg1::KSP,arg2::Vec,arg3::Ptr{Vec})
    ccall((:KSPBuildSolution,petsc),PetscErrorCode,(KSP,Vec,Ptr{Vec}),arg1,arg2,arg3)
end

function KSPBuildResidual(arg1::KSP,arg2::Vec,arg3::Vec,arg4::Ptr{Vec})
    ccall((:KSPBuildResidual,petsc),PetscErrorCode,(KSP,Vec,Vec,Ptr{Vec}),arg1,arg2,arg3,arg4)
end

function KSPRichardsonSetScale(arg1::KSP,PetscReal::Cint)
    ccall((:KSPRichardsonSetScale,petsc),PetscErrorCode,(KSP,Cint),arg1,PetscReal)
end

function KSPRichardsonSetSelfScale(arg1::KSP,arg2::PetscBool)
    ccall((:KSPRichardsonSetSelfScale,petsc),PetscErrorCode,(KSP,PetscBool),arg1,arg2)
end

function KSPChebyshevSetEigenvalues(arg1::KSP,PetscReal::Cint,arg2::Cint)
    ccall((:KSPChebyshevSetEigenvalues,petsc),PetscErrorCode,(KSP,Cint,Cint),arg1,PetscReal,arg2)
end

function KSPChebyshevEstEigSet(arg1::KSP,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::Cint)
    ccall((:KSPChebyshevEstEigSet,petsc),PetscErrorCode,(KSP,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4)
end

function KSPChebyshevEstEigSetRandom(arg1::KSP,arg2::PetscRandom)
    ccall((:KSPChebyshevEstEigSetRandom,petsc),PetscErrorCode,(KSP,PetscRandom),arg1,arg2)
end

function KSPChebyshevEstEigGetKSP(arg1::KSP,arg2::Ptr{KSP})
    ccall((:KSPChebyshevEstEigGetKSP,petsc),PetscErrorCode,(KSP,Ptr{KSP}),arg1,arg2)
end

function KSPComputeExtremeSingularValues(arg1::KSP,arg2::Ptr{Cint},arg3::Ptr{Cint})
    ccall((:KSPComputeExtremeSingularValues,petsc),PetscErrorCode,(KSP,Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3)
end

function KSPComputeEigenvalues(arg1::KSP,arg2::PetscInt,PetscReal::Ptr{Cint},arg3::Ptr{Cint},arg4::Ptr{PetscInt})
    ccall((:KSPComputeEigenvalues,petsc),PetscErrorCode,(KSP,PetscInt,Ptr{Cint},Ptr{Cint},Ptr{PetscInt}),arg1,arg2,PetscReal,arg3,arg4)
end

function KSPComputeEigenvaluesExplicitly(arg1::KSP,arg2::PetscInt,PetscReal::Ptr{Cint},arg3::Ptr{Cint})
    ccall((:KSPComputeEigenvaluesExplicitly,petsc),PetscErrorCode,(KSP,PetscInt,Ptr{Cint},Ptr{Cint}),arg1,arg2,PetscReal,arg3)
end

function KSPFCGSetMmax(arg1::KSP,arg2::PetscInt)
    ccall((:KSPFCGSetMmax,petsc),PetscErrorCode,(KSP,PetscInt),arg1,arg2)
end

function KSPFCGGetMmax(arg1::KSP,arg2::Ptr{PetscInt})
    ccall((:KSPFCGGetMmax,petsc),PetscErrorCode,(KSP,Ptr{PetscInt}),arg1,arg2)
end

function KSPFCGSetNprealloc(arg1::KSP,arg2::PetscInt)
    ccall((:KSPFCGSetNprealloc,petsc),PetscErrorCode,(KSP,PetscInt),arg1,arg2)
end

function KSPFCGGetNprealloc(arg1::KSP,arg2::Ptr{PetscInt})
    ccall((:KSPFCGGetNprealloc,petsc),PetscErrorCode,(KSP,Ptr{PetscInt}),arg1,arg2)
end

function KSPFCGSetTruncationType(arg1::KSP,arg2::KSPFCGTruncationType)
    ccall((:KSPFCGSetTruncationType,petsc),PetscErrorCode,(KSP,KSPFCGTruncationType),arg1,arg2)
end

function KSPFCGGetTruncationType(arg1::KSP,arg2::Ptr{KSPFCGTruncationType})
    ccall((:KSPFCGGetTruncationType,petsc),PetscErrorCode,(KSP,Ptr{KSPFCGTruncationType}),arg1,arg2)
end

function KSPGMRESSetRestart(arg1::KSP,arg2::PetscInt)
    ccall((:KSPGMRESSetRestart,petsc),PetscErrorCode,(KSP,PetscInt),arg1,arg2)
end

function KSPGMRESGetRestart(arg1::KSP,arg2::Ptr{PetscInt})
    ccall((:KSPGMRESGetRestart,petsc),PetscErrorCode,(KSP,Ptr{PetscInt}),arg1,arg2)
end

function KSPGMRESSetHapTol(arg1::KSP,PetscReal::Cint)
    ccall((:KSPGMRESSetHapTol,petsc),PetscErrorCode,(KSP,Cint),arg1,PetscReal)
end

function KSPGMRESSetPreAllocateVectors(arg1::KSP)
    ccall((:KSPGMRESSetPreAllocateVectors,petsc),PetscErrorCode,(KSP,),arg1)
end

function KSPGMRESSetOrthogonalization(arg1::KSP,arg2::Ptr{Void})
    ccall((:KSPGMRESSetOrthogonalization,petsc),PetscErrorCode,(KSP,Ptr{Void}),arg1,arg2)
end

function KSPGMRESGetOrthogonalization(arg1::KSP,arg2::Ptr{Ptr{Void}})
    ccall((:KSPGMRESGetOrthogonalization,petsc),PetscErrorCode,(KSP,Ptr{Ptr{Void}}),arg1,arg2)
end

function KSPGMRESModifiedGramSchmidtOrthogonalization(arg1::KSP,arg2::PetscInt)
    ccall((:KSPGMRESModifiedGramSchmidtOrthogonalization,petsc),PetscErrorCode,(KSP,PetscInt),arg1,arg2)
end

function KSPGMRESClassicalGramSchmidtOrthogonalization(arg1::KSP,arg2::PetscInt)
    ccall((:KSPGMRESClassicalGramSchmidtOrthogonalization,petsc),PetscErrorCode,(KSP,PetscInt),arg1,arg2)
end

function KSPLGMRESSetAugDim(arg1::KSP,arg2::PetscInt)
    ccall((:KSPLGMRESSetAugDim,petsc),PetscErrorCode,(KSP,PetscInt),arg1,arg2)
end

function KSPLGMRESSetConstant(arg1::KSP)
    ccall((:KSPLGMRESSetConstant,petsc),PetscErrorCode,(KSP,),arg1)
end

function KSPGCRSetRestart(arg1::KSP,arg2::PetscInt)
    ccall((:KSPGCRSetRestart,petsc),PetscErrorCode,(KSP,PetscInt),arg1,arg2)
end

function KSPGCRGetRestart(arg1::KSP,arg2::Ptr{PetscInt})
    ccall((:KSPGCRGetRestart,petsc),PetscErrorCode,(KSP,Ptr{PetscInt}),arg1,arg2)
end

function KSPGCRSetModifyPC(arg1::KSP,arg2::Ptr{Void},arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:KSPGCRSetModifyPC,petsc),PetscErrorCode,(KSP,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function KSPGMRESSetCGSRefinementType(arg1::KSP,arg2::KSPGMRESCGSRefinementType)
    ccall((:KSPGMRESSetCGSRefinementType,petsc),PetscErrorCode,(KSP,KSPGMRESCGSRefinementType),arg1,arg2)
end

function KSPGMRESGetCGSRefinementType(arg1::KSP,arg2::Ptr{KSPGMRESCGSRefinementType})
    ccall((:KSPGMRESGetCGSRefinementType,petsc),PetscErrorCode,(KSP,Ptr{KSPGMRESCGSRefinementType}),arg1,arg2)
end

function KSPFGMRESModifyPCNoChange(arg1::KSP,arg2::PetscInt,arg3::PetscInt,PetscReal::Cint,arg4::Ptr{Void})
    ccall((:KSPFGMRESModifyPCNoChange,petsc),PetscErrorCode,(KSP,PetscInt,PetscInt,Cint,Ptr{Void}),arg1,arg2,arg3,PetscReal,arg4)
end

function KSPFGMRESModifyPCKSP(arg1::KSP,arg2::PetscInt,arg3::PetscInt,PetscReal::Cint,arg4::Ptr{Void})
    ccall((:KSPFGMRESModifyPCKSP,petsc),PetscErrorCode,(KSP,PetscInt,PetscInt,Cint,Ptr{Void}),arg1,arg2,arg3,PetscReal,arg4)
end

function KSPFGMRESSetModifyPC(arg1::KSP,arg2::Ptr{Void},arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:KSPFGMRESSetModifyPC,petsc),PetscErrorCode,(KSP,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function KSPQCGSetTrustRegionRadius(arg1::KSP,PetscReal::Cint)
    ccall((:KSPQCGSetTrustRegionRadius,petsc),PetscErrorCode,(KSP,Cint),arg1,PetscReal)
end

function KSPQCGGetQuadratic(arg1::KSP,arg2::Ptr{Cint})
    ccall((:KSPQCGGetQuadratic,petsc),PetscErrorCode,(KSP,Ptr{Cint}),arg1,arg2)
end

function KSPQCGGetTrialStepNorm(arg1::KSP,arg2::Ptr{Cint})
    ccall((:KSPQCGGetTrialStepNorm,petsc),PetscErrorCode,(KSP,Ptr{Cint}),arg1,arg2)
end

function KSPBCGSLSetXRes(arg1::KSP,PetscReal::Cint)
    ccall((:KSPBCGSLSetXRes,petsc),PetscErrorCode,(KSP,Cint),arg1,PetscReal)
end

function KSPBCGSLSetPol(arg1::KSP,arg2::PetscBool)
    ccall((:KSPBCGSLSetPol,petsc),PetscErrorCode,(KSP,PetscBool),arg1,arg2)
end

function KSPBCGSLSetEll(arg1::KSP,arg2::PetscInt)
    ccall((:KSPBCGSLSetEll,petsc),PetscErrorCode,(KSP,PetscInt),arg1,arg2)
end

function KSPBCGSLSetUsePseudoinverse(arg1::KSP,arg2::PetscBool)
    ccall((:KSPBCGSLSetUsePseudoinverse,petsc),PetscErrorCode,(KSP,PetscBool),arg1,arg2)
end

function KSPSetFromOptions(arg1::KSP)
    ccall((:KSPSetFromOptions,petsc),PetscErrorCode,(KSP,),arg1)
end

function KSPAddOptionsChecker(arg1::Ptr{Void})
    ccall((:KSPAddOptionsChecker,petsc),PetscErrorCode,(Ptr{Void},),arg1)
end

function KSPMonitorSingularValue(arg1::KSP,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:KSPMonitorSingularValue,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorDefault(arg1::KSP,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:KSPMonitorDefault,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPLSQRMonitorDefault(arg1::KSP,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:KSPLSQRMonitorDefault,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorRange(arg1::KSP,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:KSPMonitorRange,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorDynamicTolerance(ksp::KSP,its::PetscInt,fnorm::Cint,dummy::Ptr{Void})
    ccall((:KSPMonitorDynamicTolerance,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),ksp,its,fnorm,dummy)
end

function KSPMonitorDynamicToleranceDestroy(dummy::Ptr{Ptr{Void}})
    ccall((:KSPMonitorDynamicToleranceDestroy,petsc),PetscErrorCode,(Ptr{Ptr{Void}},),dummy)
end

function KSPMonitorTrueResidualNorm(arg1::KSP,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:KSPMonitorTrueResidualNorm,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorTrueResidualMaxNorm(arg1::KSP,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:KSPMonitorTrueResidualMaxNorm,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorDefaultShort(arg1::KSP,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:KSPMonitorDefaultShort,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorSolution(arg1::KSP,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:KSPMonitorSolution,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorSAWs(arg1::KSP,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:KSPMonitorSAWs,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorSAWsCreate(arg1::KSP,arg2::Ptr{Ptr{Void}})
    ccall((:KSPMonitorSAWsCreate,petsc),PetscErrorCode,(KSP,Ptr{Ptr{Void}}),arg1,arg2)
end

function KSPMonitorSAWsDestroy(arg1::Ptr{Ptr{Void}})
    ccall((:KSPMonitorSAWsDestroy,petsc),PetscErrorCode,(Ptr{Ptr{Void}},),arg1)
end

function KSPGMRESMonitorKrylov(arg1::KSP,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:KSPGMRESMonitorKrylov,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPUnwindPreconditioner(arg1::KSP,arg2::Vec,arg3::Vec)
    ccall((:KSPUnwindPreconditioner,petsc),PetscErrorCode,(KSP,Vec,Vec),arg1,arg2,arg3)
end

function KSPInitialResidual(arg1::KSP,arg2::Vec,arg3::Vec,arg4::Vec,arg5::Vec,arg6::Vec)
    ccall((:KSPInitialResidual,petsc),PetscErrorCode,(KSP,Vec,Vec,Vec,Vec,Vec),arg1,arg2,arg3,arg4,arg5,arg6)
end

function KSPSetOperators(arg1::KSP,arg2::Mat,arg3::Mat)
    ccall((:KSPSetOperators,petsc),PetscErrorCode,(KSP,Mat,Mat),arg1,arg2,arg3)
end

function KSPGetOperators(arg1::KSP,arg2::Ptr{Mat},arg3::Ptr{Mat})
    ccall((:KSPGetOperators,petsc),PetscErrorCode,(KSP,Ptr{Mat},Ptr{Mat}),arg1,arg2,arg3)
end

function KSPGetOperatorsSet(arg1::KSP,arg2::Ptr{PetscBool},arg3::Ptr{PetscBool})
    ccall((:KSPGetOperatorsSet,petsc),PetscErrorCode,(KSP,Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3)
end

function KSPSetOptionsPrefix(arg1::KSP,arg2::Ptr{Uint8})
    ccall((:KSPSetOptionsPrefix,petsc),PetscErrorCode,(KSP,Ptr{Uint8}),arg1,arg2)
end

function KSPAppendOptionsPrefix(arg1::KSP,arg2::Ptr{Uint8})
    ccall((:KSPAppendOptionsPrefix,petsc),PetscErrorCode,(KSP,Ptr{Uint8}),arg1,arg2)
end

function KSPGetOptionsPrefix(arg1::KSP,arg2::Ptr{Ptr{Uint8}})
    ccall((:KSPGetOptionsPrefix,petsc),PetscErrorCode,(KSP,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function KSPSetTabLevel(arg1::KSP,arg2::PetscInt)
    ccall((:KSPSetTabLevel,petsc),PetscErrorCode,(KSP,PetscInt),arg1,arg2)
end

function KSPGetTabLevel(arg1::KSP,arg2::Ptr{PetscInt})
    ccall((:KSPGetTabLevel,petsc),PetscErrorCode,(KSP,Ptr{PetscInt}),arg1,arg2)
end

function KSPSetDiagonalScale(arg1::KSP,arg2::PetscBool)
    ccall((:KSPSetDiagonalScale,petsc),PetscErrorCode,(KSP,PetscBool),arg1,arg2)
end

function KSPGetDiagonalScale(arg1::KSP,arg2::Ptr{PetscBool})
    ccall((:KSPGetDiagonalScale,petsc),PetscErrorCode,(KSP,Ptr{PetscBool}),arg1,arg2)
end

function KSPSetDiagonalScaleFix(arg1::KSP,arg2::PetscBool)
    ccall((:KSPSetDiagonalScaleFix,petsc),PetscErrorCode,(KSP,PetscBool),arg1,arg2)
end

function KSPGetDiagonalScaleFix(arg1::KSP,arg2::Ptr{PetscBool})
    ccall((:KSPGetDiagonalScaleFix,petsc),PetscErrorCode,(KSP,Ptr{PetscBool}),arg1,arg2)
end

function KSPView(arg1::KSP,arg2::PetscViewer)
    ccall((:KSPView,petsc),PetscErrorCode,(KSP,PetscViewer),arg1,arg2)
end

function KSPLoad(arg1::KSP,arg2::PetscViewer)
    ccall((:KSPLoad,petsc),PetscErrorCode,(KSP,PetscViewer),arg1,arg2)
end

function KSPReasonViewFromOptions(arg1::KSP)
    ccall((:KSPReasonViewFromOptions,petsc),PetscErrorCode,(KSP,),arg1)
end

function KSPLSQRSetStandardErrorVec(arg1::KSP,arg2::Vec)
    ccall((:KSPLSQRSetStandardErrorVec,petsc),PetscErrorCode,(KSP,Vec),arg1,arg2)
end

function KSPLSQRGetStandardErrorVec(arg1::KSP,arg2::Ptr{Vec})
    ccall((:KSPLSQRGetStandardErrorVec,petsc),PetscErrorCode,(KSP,Ptr{Vec}),arg1,arg2)
end

function PCRedundantGetKSP(arg1::PC,arg2::Ptr{KSP})
    ccall((:PCRedundantGetKSP,petsc),PetscErrorCode,(PC,Ptr{KSP}),arg1,arg2)
end

function PCRedistributeGetKSP(arg1::PC,arg2::Ptr{KSP})
    ccall((:PCRedistributeGetKSP,petsc),PetscErrorCode,(PC,Ptr{KSP}),arg1,arg2)
end

function KSPSetNormType(arg1::KSP,arg2::KSPNormType)
    ccall((:KSPSetNormType,petsc),PetscErrorCode,(KSP,KSPNormType),arg1,arg2)
end

function KSPGetNormType(arg1::KSP,arg2::Ptr{KSPNormType})
    ccall((:KSPGetNormType,petsc),PetscErrorCode,(KSP,Ptr{KSPNormType}),arg1,arg2)
end

function KSPSetSupportedNorm(ksp::KSP,arg1::KSPNormType,arg2::PCSide,arg3::PetscInt)
    ccall((:KSPSetSupportedNorm,petsc),PetscErrorCode,(KSP,KSPNormType,PCSide,PetscInt),ksp,arg1,arg2,arg3)
end

function KSPSetCheckNormIteration(arg1::KSP,arg2::PetscInt)
    ccall((:KSPSetCheckNormIteration,petsc),PetscErrorCode,(KSP,PetscInt),arg1,arg2)
end

function KSPSetLagNorm(arg1::KSP,arg2::PetscBool)
    ccall((:KSPSetLagNorm,petsc),PetscErrorCode,(KSP,PetscBool),arg1,arg2)
end

function KSPSetConvergenceTest(arg1::KSP,arg2::Ptr{Void},arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:KSPSetConvergenceTest,petsc),PetscErrorCode,(KSP,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function KSPGetConvergenceContext(arg1::KSP,arg2::Ptr{Ptr{Void}})
    ccall((:KSPGetConvergenceContext,petsc),PetscErrorCode,(KSP,Ptr{Ptr{Void}}),arg1,arg2)
end

function KSPConvergedDefault(arg1::KSP,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{KSPConvergedReason},arg4::Ptr{Void})
    ccall((:KSPConvergedDefault,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{KSPConvergedReason},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end

function KSPConvergedLSQR(arg1::KSP,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{KSPConvergedReason},arg4::Ptr{Void})
    ccall((:KSPConvergedLSQR,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{KSPConvergedReason},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end

function KSPConvergedDefaultDestroy(arg1::Ptr{Void})
    ccall((:KSPConvergedDefaultDestroy,petsc),PetscErrorCode,(Ptr{Void},),arg1)
end

function KSPConvergedDefaultCreate(arg1::Ptr{Ptr{Void}})
    ccall((:KSPConvergedDefaultCreate,petsc),PetscErrorCode,(Ptr{Ptr{Void}},),arg1)
end

function KSPConvergedDefaultSetUIRNorm(arg1::KSP)
    ccall((:KSPConvergedDefaultSetUIRNorm,petsc),PetscErrorCode,(KSP,),arg1)
end

function KSPConvergedDefaultSetUMIRNorm(arg1::KSP)
    ccall((:KSPConvergedDefaultSetUMIRNorm,petsc),PetscErrorCode,(KSP,),arg1)
end

function KSPConvergedSkip(arg1::KSP,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{KSPConvergedReason},arg4::Ptr{Void})
    ccall((:KSPConvergedSkip,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{KSPConvergedReason},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end

function KSPGetConvergedReason(arg1::KSP,arg2::Ptr{KSPConvergedReason})
    ccall((:KSPGetConvergedReason,petsc),PetscErrorCode,(KSP,Ptr{KSPConvergedReason}),arg1,arg2)
end

function KSPCGSetType(arg1::KSP,arg2::KSPCGType)
    ccall((:KSPCGSetType,petsc),PetscErrorCode,(KSP,KSPCGType),arg1,arg2)
end

function KSPCGUseSingleReduction(arg1::KSP,arg2::PetscBool)
    ccall((:KSPCGUseSingleReduction,petsc),PetscErrorCode,(KSP,PetscBool),arg1,arg2)
end

function KSPNASHSetRadius(arg1::KSP,PetscReal::Cint)
    ccall((:KSPNASHSetRadius,petsc),PetscErrorCode,(KSP,Cint),arg1,PetscReal)
end

function KSPNASHGetNormD(arg1::KSP,arg2::Ptr{Cint})
    ccall((:KSPNASHGetNormD,petsc),PetscErrorCode,(KSP,Ptr{Cint}),arg1,arg2)
end

function KSPNASHGetObjFcn(arg1::KSP,arg2::Ptr{Cint})
    ccall((:KSPNASHGetObjFcn,petsc),PetscErrorCode,(KSP,Ptr{Cint}),arg1,arg2)
end

function KSPSTCGSetRadius(arg1::KSP,PetscReal::Cint)
    ccall((:KSPSTCGSetRadius,petsc),PetscErrorCode,(KSP,Cint),arg1,PetscReal)
end

function KSPSTCGGetNormD(arg1::KSP,arg2::Ptr{Cint})
    ccall((:KSPSTCGGetNormD,petsc),PetscErrorCode,(KSP,Ptr{Cint}),arg1,arg2)
end

function KSPSTCGGetObjFcn(arg1::KSP,arg2::Ptr{Cint})
    ccall((:KSPSTCGGetObjFcn,petsc),PetscErrorCode,(KSP,Ptr{Cint}),arg1,arg2)
end

function KSPGLTRSetRadius(arg1::KSP,PetscReal::Cint)
    ccall((:KSPGLTRSetRadius,petsc),PetscErrorCode,(KSP,Cint),arg1,PetscReal)
end

function KSPGLTRGetNormD(arg1::KSP,arg2::Ptr{Cint})
    ccall((:KSPGLTRGetNormD,petsc),PetscErrorCode,(KSP,Ptr{Cint}),arg1,arg2)
end

function KSPGLTRGetObjFcn(arg1::KSP,arg2::Ptr{Cint})
    ccall((:KSPGLTRGetObjFcn,petsc),PetscErrorCode,(KSP,Ptr{Cint}),arg1,arg2)
end

function KSPGLTRGetMinEig(arg1::KSP,arg2::Ptr{Cint})
    ccall((:KSPGLTRGetMinEig,petsc),PetscErrorCode,(KSP,Ptr{Cint}),arg1,arg2)
end

function KSPGLTRGetLambda(arg1::KSP,arg2::Ptr{Cint})
    ccall((:KSPGLTRGetLambda,petsc),PetscErrorCode,(KSP,Ptr{Cint}),arg1,arg2)
end

function KSPPythonSetType(arg1::KSP,arg2::Ptr{Uint8})
    ccall((:KSPPythonSetType,petsc),PetscErrorCode,(KSP,Ptr{Uint8}),arg1,arg2)
end

function PCPreSolve(arg1::PC,arg2::KSP)
    ccall((:PCPreSolve,petsc),PetscErrorCode,(PC,KSP),arg1,arg2)
end

function PCPostSolve(arg1::PC,arg2::KSP)
    ccall((:PCPostSolve,petsc),PetscErrorCode,(PC,KSP),arg1,arg2)
end

function KSPMonitorLGResidualNormCreate(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Cint,arg4::Cint,arg5::Cint,arg6::Cint,arg7::Ptr{Ptr{PetscObject}})
    ccall((:KSPMonitorLGResidualNormCreate,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Cint,Cint,Cint,Cint,Ptr{Ptr{PetscObject}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function KSPMonitorLGResidualNorm(arg1::KSP,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{PetscObject})
    ccall((:KSPMonitorLGResidualNorm,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{PetscObject}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorLGResidualNormDestroy(arg1::Ptr{Ptr{PetscObject}})
    ccall((:KSPMonitorLGResidualNormDestroy,petsc),PetscErrorCode,(Ptr{Ptr{PetscObject}},),arg1)
end

function KSPMonitorLGTrueResidualNormCreate(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Cint,arg4::Cint,arg5::Cint,arg6::Cint,arg7::Ptr{Ptr{PetscObject}})
    ccall((:KSPMonitorLGTrueResidualNormCreate,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Cint,Cint,Cint,Cint,Ptr{Ptr{PetscObject}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function KSPMonitorLGTrueResidualNorm(arg1::KSP,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{PetscObject})
    ccall((:KSPMonitorLGTrueResidualNorm,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{PetscObject}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorLGTrueResidualNormDestroy(arg1::Ptr{Ptr{PetscObject}})
    ccall((:KSPMonitorLGTrueResidualNormDestroy,petsc),PetscErrorCode,(Ptr{Ptr{PetscObject}},),arg1)
end

function KSPMonitorLGRange(arg1::KSP,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:KSPMonitorLGRange,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function PCShellSetPreSolve(arg1::PC,arg2::Ptr{Void})
    ccall((:PCShellSetPreSolve,petsc),PetscErrorCode,(PC,Ptr{Void}),arg1,arg2)
end

function PCShellSetPostSolve(arg1::PC,arg2::Ptr{Void})
    ccall((:PCShellSetPostSolve,petsc),PetscErrorCode,(PC,Ptr{Void}),arg1,arg2)
end

function KSPFischerGuessCreate(arg1::KSP,arg2::PetscInt,arg3::PetscInt,arg4::Ptr{KSPFischerGuess})
    ccall((:KSPFischerGuessCreate,petsc),PetscErrorCode,(KSP,PetscInt,PetscInt,Ptr{KSPFischerGuess}),arg1,arg2,arg3,arg4)
end

function KSPFischerGuessDestroy(arg1::Ptr{KSPFischerGuess})
    ccall((:KSPFischerGuessDestroy,petsc),PetscErrorCode,(Ptr{KSPFischerGuess},),arg1)
end

function KSPFischerGuessReset(arg1::KSPFischerGuess)
    ccall((:KSPFischerGuessReset,petsc),PetscErrorCode,(KSPFischerGuess,),arg1)
end

function KSPFischerGuessUpdate(arg1::KSPFischerGuess,arg2::Vec)
    ccall((:KSPFischerGuessUpdate,petsc),PetscErrorCode,(KSPFischerGuess,Vec),arg1,arg2)
end

function KSPFischerGuessFormGuess(arg1::KSPFischerGuess,arg2::Vec,arg3::Vec)
    ccall((:KSPFischerGuessFormGuess,petsc),PetscErrorCode,(KSPFischerGuess,Vec,Vec),arg1,arg2,arg3)
end

function KSPFischerGuessSetFromOptions(arg1::KSPFischerGuess)
    ccall((:KSPFischerGuessSetFromOptions,petsc),PetscErrorCode,(KSPFischerGuess,),arg1)
end

function KSPSetUseFischerGuess(arg1::KSP,arg2::PetscInt,arg3::PetscInt)
    ccall((:KSPSetUseFischerGuess,petsc),PetscErrorCode,(KSP,PetscInt,PetscInt),arg1,arg2,arg3)
end

function KSPSetFischerGuess(arg1::KSP,arg2::KSPFischerGuess)
    ccall((:KSPSetFischerGuess,petsc),PetscErrorCode,(KSP,KSPFischerGuess),arg1,arg2)
end

function KSPGetFischerGuess(arg1::KSP,arg2::Ptr{KSPFischerGuess})
    ccall((:KSPGetFischerGuess,petsc),PetscErrorCode,(KSP,Ptr{KSPFischerGuess}),arg1,arg2)
end

function MatCreateSchurComplement(arg1::Mat,arg2::Mat,arg3::Mat,arg4::Mat,arg5::Mat,arg6::Ptr{Mat})
    ccall((:MatCreateSchurComplement,petsc),PetscErrorCode,(Mat,Mat,Mat,Mat,Mat,Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatSchurComplementGetKSP(arg1::Mat,arg2::Ptr{KSP})
    ccall((:MatSchurComplementGetKSP,petsc),PetscErrorCode,(Mat,Ptr{KSP}),arg1,arg2)
end

function MatSchurComplementSetKSP(arg1::Mat,arg2::KSP)
    ccall((:MatSchurComplementSetKSP,petsc),PetscErrorCode,(Mat,KSP),arg1,arg2)
end

function MatSchurComplementSetSubMatrices(arg1::Mat,arg2::Mat,arg3::Mat,arg4::Mat,arg5::Mat,arg6::Mat)
    ccall((:MatSchurComplementSetSubMatrices,petsc),PetscErrorCode,(Mat,Mat,Mat,Mat,Mat,Mat),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatSchurComplementUpdateSubMatrices(arg1::Mat,arg2::Mat,arg3::Mat,arg4::Mat,arg5::Mat,arg6::Mat)
    ccall((:MatSchurComplementUpdateSubMatrices,petsc),PetscErrorCode,(Mat,Mat,Mat,Mat,Mat,Mat),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatSchurComplementGetSubMatrices(arg1::Mat,arg2::Ptr{Mat},arg3::Ptr{Mat},arg4::Ptr{Mat},arg5::Ptr{Mat},arg6::Ptr{Mat})
    ccall((:MatSchurComplementGetSubMatrices,petsc),PetscErrorCode,(Mat,Ptr{Mat},Ptr{Mat},Ptr{Mat},Ptr{Mat},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatSchurComplementSetAinvType(arg1::Mat,arg2::MatSchurComplementAinvType)
    ccall((:MatSchurComplementSetAinvType,petsc),PetscErrorCode,(Mat,MatSchurComplementAinvType),arg1,arg2)
end

function MatSchurComplementGetAinvType(arg1::Mat,arg2::Ptr{MatSchurComplementAinvType})
    ccall((:MatSchurComplementGetAinvType,petsc),PetscErrorCode,(Mat,Ptr{MatSchurComplementAinvType}),arg1,arg2)
end

function MatSchurComplementGetPmat(arg1::Mat,arg2::MatReuse,arg3::Ptr{Mat})
    ccall((:MatSchurComplementGetPmat,petsc),PetscErrorCode,(Mat,MatReuse,Ptr{Mat}),arg1,arg2,arg3)
end

function MatSchurComplementComputeExplicitOperator(arg1::Mat,arg2::Ptr{Mat})
    ccall((:MatSchurComplementComputeExplicitOperator,petsc),PetscErrorCode,(Mat,Ptr{Mat}),arg1,arg2)
end

function MatGetSchurComplement(arg1::Mat,arg2::IS,arg3::IS,arg4::IS,arg5::IS,arg6::MatReuse,arg7::Ptr{Mat},arg8::MatSchurComplementAinvType,arg9::MatReuse,arg10::Ptr{Mat})
    ccall((:MatGetSchurComplement,petsc),PetscErrorCode,(Mat,IS,IS,IS,IS,MatReuse,Ptr{Mat},MatSchurComplementAinvType,MatReuse,Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function MatCreateSchurComplementPmat(arg1::Mat,arg2::Mat,arg3::Mat,arg4::Mat,arg5::MatSchurComplementAinvType,arg6::MatReuse,arg7::Ptr{Mat})
    ccall((:MatCreateSchurComplementPmat,petsc),PetscErrorCode,(Mat,Mat,Mat,Mat,MatSchurComplementAinvType,MatReuse,Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function KSPSetDM(arg1::KSP,arg2::DM)
    ccall((:KSPSetDM,petsc),PetscErrorCode,(KSP,DM),arg1,arg2)
end

function KSPSetDMActive(arg1::KSP,arg2::PetscBool)
    ccall((:KSPSetDMActive,petsc),PetscErrorCode,(KSP,PetscBool),arg1,arg2)
end

function KSPGetDM(arg1::KSP,arg2::Ptr{DM})
    ccall((:KSPGetDM,petsc),PetscErrorCode,(KSP,Ptr{DM}),arg1,arg2)
end

function KSPSetApplicationContext(arg1::KSP,arg2::Ptr{Void})
    ccall((:KSPSetApplicationContext,petsc),PetscErrorCode,(KSP,Ptr{Void}),arg1,arg2)
end

function KSPGetApplicationContext(arg1::KSP,arg2::Ptr{Void})
    ccall((:KSPGetApplicationContext,petsc),PetscErrorCode,(KSP,Ptr{Void}),arg1,arg2)
end

function KSPSetComputeRHS(arg1::KSP,func::Ptr{Void},arg2::Ptr{Void})
    ccall((:KSPSetComputeRHS,petsc),PetscErrorCode,(KSP,Ptr{Void},Ptr{Void}),arg1,func,arg2)
end

function KSPSetComputeOperators(arg1::KSP,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:KSPSetComputeOperators,petsc),PetscErrorCode,(KSP,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function KSPSetComputeInitialGuess(arg1::KSP,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:KSPSetComputeInitialGuess,petsc),PetscErrorCode,(KSP,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function DMKSPSetComputeOperators(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:DMKSPSetComputeOperators,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function DMKSPGetComputeOperators(arg1::DM,arg2::Ptr{Ptr{Void}},arg3::Ptr{Void})
    ccall((:DMKSPGetComputeOperators,petsc),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Void}),arg1,arg2,arg3)
end

function DMKSPSetComputeRHS(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:DMKSPSetComputeRHS,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function DMKSPGetComputeRHS(arg1::DM,arg2::Ptr{Ptr{Void}},arg3::Ptr{Void})
    ccall((:DMKSPGetComputeRHS,petsc),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Void}),arg1,arg2,arg3)
end

function DMKSPSetComputeInitialGuess(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:DMKSPSetComputeInitialGuess,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function DMKSPGetComputeInitialGuess(arg1::DM,arg2::Ptr{Ptr{Void}},arg3::Ptr{Void})
    ccall((:DMKSPGetComputeInitialGuess,petsc),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Void}),arg1,arg2,arg3)
end

function DMGlobalToLocalSolve(arg1::DM,arg2::Vec,arg3::Vec)
    ccall((:DMGlobalToLocalSolve,petsc),PetscErrorCode,(DM,Vec,Vec),arg1,arg2,arg3)
end

function DMPlexProjectField(arg1::DM,arg2::Vec,arg3::Ptr{Ptr{Void}},arg4::InsertMode,arg5::Vec)
    ccall((:DMPlexProjectField,petsc),PetscErrorCode,(DM,Vec,Ptr{Ptr{Void}},InsertMode,Vec),arg1,arg2,arg3,arg4,arg5)
end

function SNESInitializePackage()
    ccall((:SNESInitializePackage,petsc),PetscErrorCode,())
end

function SNESCreate(arg1::MPI_Comm,arg2::Ptr{SNES})
    ccall((:SNESCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{SNES}),arg1,arg2)
end

function SNESReset(arg1::SNES)
    ccall((:SNESReset,petsc),PetscErrorCode,(SNES,),arg1)
end

function SNESDestroy(arg1::Ptr{SNES})
    ccall((:SNESDestroy,petsc),PetscErrorCode,(Ptr{SNES},),arg1)
end

function SNESSetType(arg1::SNES,arg2::SNESType)
    ccall((:SNESSetType,petsc),PetscErrorCode,(SNES,SNESType),arg1,arg2)
end

function SNESMonitor(arg1::SNES,arg2::PetscInt,PetscReal::Cint)
    ccall((:SNESMonitor,petsc),PetscErrorCode,(SNES,PetscInt,Cint),arg1,arg2,PetscReal)
end

function SNESMonitorSet(arg1::SNES,arg2::Ptr{Void},arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:SNESMonitorSet,petsc),PetscErrorCode,(SNES,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function SNESMonitorCancel(arg1::SNES)
    ccall((:SNESMonitorCancel,petsc),PetscErrorCode,(SNES,),arg1)
end

function SNESMonitorSAWs(arg1::SNES,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:SNESMonitorSAWs,petsc),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function SNESMonitorSAWsCreate(arg1::SNES,arg2::Ptr{Ptr{Void}})
    ccall((:SNESMonitorSAWsCreate,petsc),PetscErrorCode,(SNES,Ptr{Ptr{Void}}),arg1,arg2)
end

function SNESMonitorSAWsDestroy(arg1::Ptr{Ptr{Void}})
    ccall((:SNESMonitorSAWsDestroy,petsc),PetscErrorCode,(Ptr{Ptr{Void}},),arg1)
end

function SNESSetConvergenceHistory(arg1::SNES,PetscReal::Ptr{Cint},arg2::Ptr{PetscInt},arg3::PetscInt,arg4::PetscBool)
    ccall((:SNESSetConvergenceHistory,petsc),PetscErrorCode,(SNES,Ptr{Cint},Ptr{PetscInt},PetscInt,PetscBool),arg1,PetscReal,arg2,arg3,arg4)
end

function SNESGetConvergenceHistory(arg1::SNES,arg2::Ptr{Ptr{Cint}},arg3::Ptr{Ptr{PetscInt}},arg4::Ptr{PetscInt})
    ccall((:SNESGetConvergenceHistory,petsc),PetscErrorCode,(SNES,Ptr{Ptr{Cint}},Ptr{Ptr{PetscInt}},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function SNESSetUp(arg1::SNES)
    ccall((:SNESSetUp,petsc),PetscErrorCode,(SNES,),arg1)
end

function SNESSolve(arg1::SNES,arg2::Vec,arg3::Vec)
    ccall((:SNESSolve,petsc),PetscErrorCode,(SNES,Vec,Vec),arg1,arg2,arg3)
end

function SNESSetErrorIfNotConverged(arg1::SNES,arg2::PetscBool)
    ccall((:SNESSetErrorIfNotConverged,petsc),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end

function SNESGetErrorIfNotConverged(arg1::SNES,arg2::Ptr{PetscBool})
    ccall((:SNESGetErrorIfNotConverged,petsc),PetscErrorCode,(SNES,Ptr{PetscBool}),arg1,arg2)
end

function SNESSetWorkVecs(arg1::SNES,arg2::PetscInt)
    ccall((:SNESSetWorkVecs,petsc),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end

function SNESAddOptionsChecker(arg1::Ptr{Void})
    ccall((:SNESAddOptionsChecker,petsc),PetscErrorCode,(Ptr{Void},),arg1)
end

function SNESSetUpdate(arg1::SNES,arg2::Ptr{Void})
    ccall((:SNESSetUpdate,petsc),PetscErrorCode,(SNES,Ptr{Void}),arg1,arg2)
end

function SNESRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:SNESRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function SNESGetKSP(arg1::SNES,arg2::Ptr{KSP})
    ccall((:SNESGetKSP,petsc),PetscErrorCode,(SNES,Ptr{KSP}),arg1,arg2)
end

function SNESSetKSP(arg1::SNES,arg2::KSP)
    ccall((:SNESSetKSP,petsc),PetscErrorCode,(SNES,KSP),arg1,arg2)
end

function SNESSetSolution(arg1::SNES,arg2::Vec)
    ccall((:SNESSetSolution,petsc),PetscErrorCode,(SNES,Vec),arg1,arg2)
end

function SNESGetSolution(arg1::SNES,arg2::Ptr{Vec})
    ccall((:SNESGetSolution,petsc),PetscErrorCode,(SNES,Ptr{Vec}),arg1,arg2)
end

function SNESGetSolutionUpdate(arg1::SNES,arg2::Ptr{Vec})
    ccall((:SNESGetSolutionUpdate,petsc),PetscErrorCode,(SNES,Ptr{Vec}),arg1,arg2)
end

function SNESGetRhs(arg1::SNES,arg2::Ptr{Vec})
    ccall((:SNESGetRhs,petsc),PetscErrorCode,(SNES,Ptr{Vec}),arg1,arg2)
end

function SNESView(arg1::SNES,arg2::PetscViewer)
    ccall((:SNESView,petsc),PetscErrorCode,(SNES,PetscViewer),arg1,arg2)
end

function SNESLoad(arg1::SNES,arg2::PetscViewer)
    ccall((:SNESLoad,petsc),PetscErrorCode,(SNES,PetscViewer),arg1,arg2)
end

function SNESReasonViewFromOptions(arg1::SNES)
    ccall((:SNESReasonViewFromOptions,petsc),PetscErrorCode,(SNES,),arg1)
end

function SNESSetOptionsPrefix(arg1::SNES,arg2::Ptr{Uint8})
    ccall((:SNESSetOptionsPrefix,petsc),PetscErrorCode,(SNES,Ptr{Uint8}),arg1,arg2)
end

function SNESAppendOptionsPrefix(arg1::SNES,arg2::Ptr{Uint8})
    ccall((:SNESAppendOptionsPrefix,petsc),PetscErrorCode,(SNES,Ptr{Uint8}),arg1,arg2)
end

function SNESGetOptionsPrefix(arg1::SNES,arg2::Ptr{Ptr{Uint8}})
    ccall((:SNESGetOptionsPrefix,petsc),PetscErrorCode,(SNES,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function SNESSetFromOptions(arg1::SNES)
    ccall((:SNESSetFromOptions,petsc),PetscErrorCode,(SNES,),arg1)
end

function MatCreateSNESMF(arg1::SNES,arg2::Ptr{Mat})
    ccall((:MatCreateSNESMF,petsc),PetscErrorCode,(SNES,Ptr{Mat}),arg1,arg2)
end

function MatMFFDComputeJacobian(arg1::SNES,arg2::Vec,arg3::Mat,arg4::Mat,arg5::Ptr{Void})
    ccall((:MatMFFDComputeJacobian,petsc),PetscErrorCode,(SNES,Vec,Mat,Mat,Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end

function MatDAADSetSNES(arg1::Mat,arg2::SNES)
    ccall((:MatDAADSetSNES,petsc),PetscErrorCode,(Mat,SNES),arg1,arg2)
end

function SNESGetType(arg1::SNES,arg2::Ptr{SNESType})
    ccall((:SNESGetType,petsc),PetscErrorCode,(SNES,Ptr{SNESType}),arg1,arg2)
end

function SNESMonitorDefault(arg1::SNES,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:SNESMonitorDefault,petsc),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function SNESMonitorRange(arg1::SNES,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:SNESMonitorRange,petsc),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function SNESMonitorRatio(arg1::SNES,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:SNESMonitorRatio,petsc),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function SNESMonitorSetRatio(arg1::SNES,arg2::PetscViewer)
    ccall((:SNESMonitorSetRatio,petsc),PetscErrorCode,(SNES,PetscViewer),arg1,arg2)
end

function SNESMonitorSolution(arg1::SNES,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:SNESMonitorSolution,petsc),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function SNESMonitorResidual(arg1::SNES,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:SNESMonitorResidual,petsc),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function SNESMonitorSolutionUpdate(arg1::SNES,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:SNESMonitorSolutionUpdate,petsc),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function SNESMonitorDefaultShort(arg1::SNES,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:SNESMonitorDefaultShort,petsc),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function SNESMonitorDefaultField(arg1::SNES,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:SNESMonitorDefaultField,petsc),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function SNESMonitorJacUpdateSpectrum(arg1::SNES,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:SNESMonitorJacUpdateSpectrum,petsc),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function SNESMonitorFields(arg1::SNES,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:SNESMonitorFields,petsc),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorSNES(arg1::KSP,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:KSPMonitorSNES,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorSNESLGResidualNormCreate(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Cint,arg4::Cint,arg5::Cint,arg6::Cint,arg7::Ptr{Ptr{PetscObject}})
    ccall((:KSPMonitorSNESLGResidualNormCreate,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Cint,Cint,Cint,Cint,Ptr{Ptr{PetscObject}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function KSPMonitorSNESLGResidualNorm(arg1::KSP,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{PetscObject})
    ccall((:KSPMonitorSNESLGResidualNorm,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{PetscObject}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorSNESLGResidualNormDestroy(arg1::Ptr{Ptr{PetscObject}})
    ccall((:KSPMonitorSNESLGResidualNormDestroy,petsc),PetscErrorCode,(Ptr{Ptr{PetscObject}},),arg1)
end

function SNESSetTolerances(arg1::SNES,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::PetscInt,arg5::PetscInt)
    ccall((:SNESSetTolerances,petsc),PetscErrorCode,(SNES,Cint,Cint,Cint,PetscInt,PetscInt),arg1,PetscReal,arg2,arg3,arg4,arg5)
end

function SNESGetTolerances(arg1::SNES,arg2::Ptr{Cint},arg3::Ptr{Cint},arg4::Ptr{Cint},arg5::Ptr{PetscInt},arg6::Ptr{PetscInt})
    ccall((:SNESGetTolerances,petsc),PetscErrorCode,(SNES,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function SNESSetTrustRegionTolerance(arg1::SNES,PetscReal::Cint)
    ccall((:SNESSetTrustRegionTolerance,petsc),PetscErrorCode,(SNES,Cint),arg1,PetscReal)
end

function SNESGetIterationNumber(arg1::SNES,arg2::Ptr{PetscInt})
    ccall((:SNESGetIterationNumber,petsc),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end

function SNESSetIterationNumber(arg1::SNES,arg2::PetscInt)
    ccall((:SNESSetIterationNumber,petsc),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end

function SNESGetNonlinearStepFailures(arg1::SNES,arg2::Ptr{PetscInt})
    ccall((:SNESGetNonlinearStepFailures,petsc),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end

function SNESSetMaxNonlinearStepFailures(arg1::SNES,arg2::PetscInt)
    ccall((:SNESSetMaxNonlinearStepFailures,petsc),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end

function SNESGetMaxNonlinearStepFailures(arg1::SNES,arg2::Ptr{PetscInt})
    ccall((:SNESGetMaxNonlinearStepFailures,petsc),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end

function SNESGetNumberFunctionEvals(arg1::SNES,arg2::Ptr{PetscInt})
    ccall((:SNESGetNumberFunctionEvals,petsc),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end

function SNESSetLagPreconditioner(arg1::SNES,arg2::PetscInt)
    ccall((:SNESSetLagPreconditioner,petsc),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end

function SNESGetLagPreconditioner(arg1::SNES,arg2::Ptr{PetscInt})
    ccall((:SNESGetLagPreconditioner,petsc),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end

function SNESSetLagJacobian(arg1::SNES,arg2::PetscInt)
    ccall((:SNESSetLagJacobian,petsc),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end

function SNESGetLagJacobian(arg1::SNES,arg2::Ptr{PetscInt})
    ccall((:SNESGetLagJacobian,petsc),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end

function SNESSetLagPreconditionerPersists(arg1::SNES,arg2::PetscBool)
    ccall((:SNESSetLagPreconditionerPersists,petsc),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end

function SNESSetLagJacobianPersists(arg1::SNES,arg2::PetscBool)
    ccall((:SNESSetLagJacobianPersists,petsc),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end

function SNESSetGridSequence(arg1::SNES,arg2::PetscInt)
    ccall((:SNESSetGridSequence,petsc),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end

function SNESGetGridSequence(arg1::SNES,arg2::Ptr{PetscInt})
    ccall((:SNESGetGridSequence,petsc),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end

function SNESGetLinearSolveIterations(arg1::SNES,arg2::Ptr{PetscInt})
    ccall((:SNESGetLinearSolveIterations,petsc),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end

function SNESGetLinearSolveFailures(arg1::SNES,arg2::Ptr{PetscInt})
    ccall((:SNESGetLinearSolveFailures,petsc),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end

function SNESSetMaxLinearSolveFailures(arg1::SNES,arg2::PetscInt)
    ccall((:SNESSetMaxLinearSolveFailures,petsc),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end

function SNESGetMaxLinearSolveFailures(arg1::SNES,arg2::Ptr{PetscInt})
    ccall((:SNESGetMaxLinearSolveFailures,petsc),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end

function SNESSetCountersReset(arg1::SNES,arg2::PetscBool)
    ccall((:SNESSetCountersReset,petsc),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end

function SNESKSPSetUseEW(arg1::SNES,arg2::PetscBool)
    ccall((:SNESKSPSetUseEW,petsc),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end

function SNESKSPGetUseEW(arg1::SNES,arg2::Ptr{PetscBool})
    ccall((:SNESKSPGetUseEW,petsc),PetscErrorCode,(SNES,Ptr{PetscBool}),arg1,arg2)
end

function SNESKSPSetParametersEW(arg1::SNES,arg2::PetscInt,PetscReal::Cint,arg3::Cint,arg4::Cint,arg5::Cint,arg6::Cint,arg7::Cint)
    ccall((:SNESKSPSetParametersEW,petsc),PetscErrorCode,(SNES,PetscInt,Cint,Cint,Cint,Cint,Cint,Cint),arg1,arg2,PetscReal,arg3,arg4,arg5,arg6,arg7)
end

function SNESKSPGetParametersEW(arg1::SNES,arg2::Ptr{PetscInt},arg3::Ptr{Cint},arg4::Ptr{Cint},arg5::Ptr{Cint},arg6::Ptr{Cint},arg7::Ptr{Cint},arg8::Ptr{Cint})
    ccall((:SNESKSPGetParametersEW,petsc),PetscErrorCode,(SNES,Ptr{PetscInt},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function SNESMonitorLGCreate(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Cint,arg4::Cint,arg5::Cint,arg6::Cint,arg7::Ptr{Ptr{PetscObject}})
    ccall((:SNESMonitorLGCreate,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Cint,Cint,Cint,Cint,Ptr{Ptr{PetscObject}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function SNESMonitorLGResidualNorm(arg1::SNES,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{PetscObject})
    ccall((:SNESMonitorLGResidualNorm,petsc),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{PetscObject}),arg1,arg2,PetscReal,arg3)
end

function SNESMonitorLGDestroy(arg1::Ptr{Ptr{PetscObject}})
    ccall((:SNESMonitorLGDestroy,petsc),PetscErrorCode,(Ptr{Ptr{PetscObject}},),arg1)
end

function SNESMonitorLGRange(arg1::SNES,arg2::PetscInt,PetscReal::Cint,arg3::Ptr{Void})
    ccall((:SNESMonitorLGRange,petsc),PetscErrorCode,(SNES,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function SNESSetApplicationContext(arg1::SNES,arg2::Ptr{Void})
    ccall((:SNESSetApplicationContext,petsc),PetscErrorCode,(SNES,Ptr{Void}),arg1,arg2)
end

function SNESGetApplicationContext(arg1::SNES,arg2::Ptr{Void})
    ccall((:SNESGetApplicationContext,petsc),PetscErrorCode,(SNES,Ptr{Void}),arg1,arg2)
end

function SNESSetComputeApplicationContext(arg1::SNES,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:SNESSetComputeApplicationContext,petsc),PetscErrorCode,(SNES,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function SNESPythonSetType(arg1::SNES,arg2::Ptr{Uint8})
    ccall((:SNESPythonSetType,petsc),PetscErrorCode,(SNES,Ptr{Uint8}),arg1,arg2)
end

function SNESSetFunctionDomainError(arg1::SNES)
    ccall((:SNESSetFunctionDomainError,petsc),PetscErrorCode,(SNES,),arg1)
end

function SNESGetFunctionDomainError(arg1::SNES,arg2::Ptr{PetscBool})
    ccall((:SNESGetFunctionDomainError,petsc),PetscErrorCode,(SNES,Ptr{PetscBool}),arg1,arg2)
end

function SNESSetConvergenceTest(arg1::SNES,arg2::Ptr{Void},arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:SNESSetConvergenceTest,petsc),PetscErrorCode,(SNES,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function SNESConvergedDefault(arg1::SNES,arg2::PetscInt,PetscReal::Cint,arg3::Cint,arg4::Cint,arg5::Ptr{SNESConvergedReason},arg6::Ptr{Void})
    ccall((:SNESConvergedDefault,petsc),PetscErrorCode,(SNES,PetscInt,Cint,Cint,Cint,Ptr{SNESConvergedReason},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4,arg5,arg6)
end

function SNESConvergedSkip(arg1::SNES,arg2::PetscInt,PetscReal::Cint,arg3::Cint,arg4::Cint,arg5::Ptr{SNESConvergedReason},arg6::Ptr{Void})
    ccall((:SNESConvergedSkip,petsc),PetscErrorCode,(SNES,PetscInt,Cint,Cint,Cint,Ptr{SNESConvergedReason},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4,arg5,arg6)
end

function SNESGetConvergedReason(arg1::SNES,arg2::Ptr{SNESConvergedReason})
    ccall((:SNESGetConvergedReason,petsc),PetscErrorCode,(SNES,Ptr{SNESConvergedReason}),arg1,arg2)
end

function SNESGetFunction(arg1::SNES,arg2::Ptr{Vec},arg3::Ptr{Ptr{Void}},arg4::Ptr{Ptr{Void}})
    ccall((:SNESGetFunction,petsc),PetscErrorCode,(SNES,Ptr{Vec},Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4)
end

function SNESComputeFunction(arg1::SNES,arg2::Vec,arg3::Vec)
    ccall((:SNESComputeFunction,petsc),PetscErrorCode,(SNES,Vec,Vec),arg1,arg2,arg3)
end

function SNESSetJacobian(arg1::SNES,arg2::Mat,arg3::Mat,arg4::Ptr{Void},arg5::Ptr{Void})
    ccall((:SNESSetJacobian,petsc),PetscErrorCode,(SNES,Mat,Mat,Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end

function SNESGetJacobian(arg1::SNES,arg2::Ptr{Mat},arg3::Ptr{Mat},arg4::Ptr{Ptr{Void}},arg5::Ptr{Ptr{Void}})
    ccall((:SNESGetJacobian,petsc),PetscErrorCode,(SNES,Ptr{Mat},Ptr{Mat},Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4,arg5)
end

function SNESObjectiveComputeFunctionDefaultFD(arg1::SNES,arg2::Vec,arg3::Vec,arg4::Ptr{Void})
    ccall((:SNESObjectiveComputeFunctionDefaultFD,petsc),PetscErrorCode,(SNES,Vec,Vec,Ptr{Void}),arg1,arg2,arg3,arg4)
end

function SNESComputeJacobianDefault(arg1::SNES,arg2::Vec,arg3::Mat,arg4::Mat,arg5::Ptr{Void})
    ccall((:SNESComputeJacobianDefault,petsc),PetscErrorCode,(SNES,Vec,Mat,Mat,Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end

function SNESComputeJacobianDefaultColor(arg1::SNES,arg2::Vec,arg3::Mat,arg4::Mat,arg5::Ptr{Void})
    ccall((:SNESComputeJacobianDefaultColor,petsc),PetscErrorCode,(SNES,Vec,Mat,Mat,Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end

function SNESSetComputeInitialGuess(arg1::SNES,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:SNESSetComputeInitialGuess,petsc),PetscErrorCode,(SNES,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function SNESSetPicard(arg1::SNES,arg2::Vec,arg3::Ptr{Void},arg4::Mat,arg5::Mat,arg6::Ptr{Void},arg7::Ptr{Void})
    ccall((:SNESSetPicard,petsc),PetscErrorCode,(SNES,Vec,Ptr{Void},Mat,Mat,Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function SNESGetPicard(arg1::SNES,arg2::Ptr{Vec},arg3::Ptr{Ptr{Void}},arg4::Ptr{Mat},arg5::Ptr{Mat},arg6::Ptr{Ptr{Void}},arg7::Ptr{Ptr{Void}})
    ccall((:SNESGetPicard,petsc),PetscErrorCode,(SNES,Ptr{Vec},Ptr{Ptr{Void}},Ptr{Mat},Ptr{Mat},Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function SNESSetInitialFunction(arg1::SNES,arg2::Vec)
    ccall((:SNESSetInitialFunction,petsc),PetscErrorCode,(SNES,Vec),arg1,arg2)
end

function SNESSetObjective(arg1::SNES,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:SNESSetObjective,petsc),PetscErrorCode,(SNES,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function SNESGetObjective(arg1::SNES,arg2::Ptr{Ptr{Void}},arg3::Ptr{Ptr{Void}})
    ccall((:SNESGetObjective,petsc),PetscErrorCode,(SNES,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end

function SNESComputeObjective(arg1::SNES,arg2::Vec,arg3::Ptr{Cint})
    ccall((:SNESComputeObjective,petsc),PetscErrorCode,(SNES,Vec,Ptr{Cint}),arg1,arg2,arg3)
end

function SNESSetNormSchedule(arg1::SNES,arg2::SNESNormSchedule)
    ccall((:SNESSetNormSchedule,petsc),PetscErrorCode,(SNES,SNESNormSchedule),arg1,arg2)
end

function SNESGetNormSchedule(arg1::SNES,arg2::Ptr{SNESNormSchedule})
    ccall((:SNESGetNormSchedule,petsc),PetscErrorCode,(SNES,Ptr{SNESNormSchedule}),arg1,arg2)
end

function SNESSetFunctionType(arg1::SNES,arg2::SNESFunctionType)
    ccall((:SNESSetFunctionType,petsc),PetscErrorCode,(SNES,SNESFunctionType),arg1,arg2)
end

function SNESGetFunctionType(arg1::SNES,arg2::Ptr{SNESFunctionType})
    ccall((:SNESGetFunctionType,petsc),PetscErrorCode,(SNES,Ptr{SNESFunctionType}),arg1,arg2)
end

function SNESSetNGS(arg1::SNES,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:SNESSetNGS,petsc),PetscErrorCode,(SNES,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function SNESGetNGS(arg1::SNES,arg2::Ptr{Ptr{Void}},arg3::Ptr{Ptr{Void}})
    ccall((:SNESGetNGS,petsc),PetscErrorCode,(SNES,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end

function SNESSetUseNGS(arg1::SNES,arg2::PetscBool)
    ccall((:SNESSetUseNGS,petsc),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end

function SNESGetUseNGS(arg1::SNES,arg2::Ptr{PetscBool})
    ccall((:SNESGetUseNGS,petsc),PetscErrorCode,(SNES,Ptr{PetscBool}),arg1,arg2)
end

function SNESComputeNGS(arg1::SNES,arg2::Vec,arg3::Vec)
    ccall((:SNESComputeNGS,petsc),PetscErrorCode,(SNES,Vec,Vec),arg1,arg2,arg3)
end

function SNESNGSSetSweeps(arg1::SNES,arg2::PetscInt)
    ccall((:SNESNGSSetSweeps,petsc),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end

function SNESNGSGetSweeps(arg1::SNES,arg2::Ptr{PetscInt})
    ccall((:SNESNGSGetSweeps,petsc),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end

function SNESNGSSetTolerances(arg1::SNES,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::PetscInt)
    ccall((:SNESNGSSetTolerances,petsc),PetscErrorCode,(SNES,Cint,Cint,Cint,PetscInt),arg1,PetscReal,arg2,arg3,arg4)
end

function SNESNGSGetTolerances(arg1::SNES,arg2::Ptr{Cint},arg3::Ptr{Cint},arg4::Ptr{Cint},arg5::Ptr{PetscInt})
    ccall((:SNESNGSGetTolerances,petsc),PetscErrorCode,(SNES,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end

function SNESUpdateCheckJacobian(arg1::SNES,arg2::PetscInt)
    ccall((:SNESUpdateCheckJacobian,petsc),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end

function SNESShellGetContext(arg1::SNES,arg2::Ptr{Ptr{Void}})
    ccall((:SNESShellGetContext,petsc),PetscErrorCode,(SNES,Ptr{Ptr{Void}}),arg1,arg2)
end

function SNESShellSetContext(arg1::SNES,arg2::Ptr{Void})
    ccall((:SNESShellSetContext,petsc),PetscErrorCode,(SNES,Ptr{Void}),arg1,arg2)
end

function SNESShellSetSolve(arg1::SNES,arg2::Ptr{Void})
    ccall((:SNESShellSetSolve,petsc),PetscErrorCode,(SNES,Ptr{Void}),arg1,arg2)
end

function SNESLineSearchCreate(arg1::MPI_Comm,arg2::Ptr{SNESLineSearch})
    ccall((:SNESLineSearchCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{SNESLineSearch}),arg1,arg2)
end

function SNESLineSearchReset(arg1::SNESLineSearch)
    ccall((:SNESLineSearchReset,petsc),PetscErrorCode,(SNESLineSearch,),arg1)
end

function SNESLineSearchView(arg1::SNESLineSearch,arg2::PetscViewer)
    ccall((:SNESLineSearchView,petsc),PetscErrorCode,(SNESLineSearch,PetscViewer),arg1,arg2)
end

function SNESLineSearchDestroy(arg1::Ptr{SNESLineSearch})
    ccall((:SNESLineSearchDestroy,petsc),PetscErrorCode,(Ptr{SNESLineSearch},),arg1)
end

function SNESLineSearchSetType(arg1::SNESLineSearch,arg2::SNESLineSearchType)
    ccall((:SNESLineSearchSetType,petsc),PetscErrorCode,(SNESLineSearch,SNESLineSearchType),arg1,arg2)
end

function SNESLineSearchSetFromOptions(arg1::SNESLineSearch)
    ccall((:SNESLineSearchSetFromOptions,petsc),PetscErrorCode,(SNESLineSearch,),arg1)
end

function SNESLineSearchSetFunction(arg1::SNESLineSearch,arg2::Ptr{Void})
    ccall((:SNESLineSearchSetFunction,petsc),PetscErrorCode,(SNESLineSearch,Ptr{Void}),arg1,arg2)
end

function SNESLineSearchSetUp(arg1::SNESLineSearch)
    ccall((:SNESLineSearchSetUp,petsc),PetscErrorCode,(SNESLineSearch,),arg1)
end

function SNESLineSearchApply(arg1::SNESLineSearch,arg2::Vec,arg3::Vec,arg4::Ptr{Cint},arg5::Vec)
    ccall((:SNESLineSearchApply,petsc),PetscErrorCode,(SNESLineSearch,Vec,Vec,Ptr{Cint},Vec),arg1,arg2,arg3,arg4,arg5)
end

function SNESLineSearchPreCheck(arg1::SNESLineSearch,arg2::Vec,arg3::Vec,arg4::Ptr{PetscBool})
    ccall((:SNESLineSearchPreCheck,petsc),PetscErrorCode,(SNESLineSearch,Vec,Vec,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function SNESLineSearchPostCheck(arg1::SNESLineSearch,arg2::Vec,arg3::Vec,arg4::Vec,arg5::Ptr{PetscBool},arg6::Ptr{PetscBool})
    ccall((:SNESLineSearchPostCheck,petsc),PetscErrorCode,(SNESLineSearch,Vec,Vec,Vec,Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function SNESLineSearchSetWorkVecs(arg1::SNESLineSearch,arg2::PetscInt)
    ccall((:SNESLineSearchSetWorkVecs,petsc),PetscErrorCode,(SNESLineSearch,PetscInt),arg1,arg2)
end

function SNESLineSearchSetPreCheck(arg1::SNESLineSearch,arg2::Ptr{Void},ctx::Ptr{Void})
    ccall((:SNESLineSearchSetPreCheck,petsc),PetscErrorCode,(SNESLineSearch,Ptr{Void},Ptr{Void}),arg1,arg2,ctx)
end

function SNESLineSearchSetPostCheck(arg1::SNESLineSearch,arg2::Ptr{Void},ctx::Ptr{Void})
    ccall((:SNESLineSearchSetPostCheck,petsc),PetscErrorCode,(SNESLineSearch,Ptr{Void},Ptr{Void}),arg1,arg2,ctx)
end

function SNESLineSearchGetPreCheck(arg1::SNESLineSearch,arg2::Ptr{Ptr{Void}},ctx::Ptr{Ptr{Void}})
    ccall((:SNESLineSearchGetPreCheck,petsc),PetscErrorCode,(SNESLineSearch,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,ctx)
end

function SNESLineSearchGetPostCheck(arg1::SNESLineSearch,arg2::Ptr{Ptr{Void}},ctx::Ptr{Ptr{Void}})
    ccall((:SNESLineSearchGetPostCheck,petsc),PetscErrorCode,(SNESLineSearch,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,ctx)
end

function SNESLineSearchSetVIFunctions(arg1::SNESLineSearch,arg2::SNESLineSearchVIProjectFunc,arg3::SNESLineSearchVINormFunc)
    ccall((:SNESLineSearchSetVIFunctions,petsc),PetscErrorCode,(SNESLineSearch,SNESLineSearchVIProjectFunc,SNESLineSearchVINormFunc),arg1,arg2,arg3)
end

function SNESLineSearchGetVIFunctions(arg1::SNESLineSearch,arg2::Ptr{SNESLineSearchVIProjectFunc},arg3::Ptr{SNESLineSearchVINormFunc})
    ccall((:SNESLineSearchGetVIFunctions,petsc),PetscErrorCode,(SNESLineSearch,Ptr{SNESLineSearchVIProjectFunc},Ptr{SNESLineSearchVINormFunc}),arg1,arg2,arg3)
end

function SNESLineSearchSetSNES(arg1::SNESLineSearch,arg2::SNES)
    ccall((:SNESLineSearchSetSNES,petsc),PetscErrorCode,(SNESLineSearch,SNES),arg1,arg2)
end

function SNESLineSearchGetSNES(arg1::SNESLineSearch,arg2::Ptr{SNES})
    ccall((:SNESLineSearchGetSNES,petsc),PetscErrorCode,(SNESLineSearch,Ptr{SNES}),arg1,arg2)
end

function SNESLineSearchGetTolerances(arg1::SNESLineSearch,arg2::Ptr{Cint},arg3::Ptr{Cint},arg4::Ptr{Cint},arg5::Ptr{Cint},arg6::Ptr{Cint},arg7::Ptr{PetscInt})
    ccall((:SNESLineSearchGetTolerances,petsc),PetscErrorCode,(SNESLineSearch,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function SNESLineSearchSetTolerances(arg1::SNESLineSearch,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::Cint,arg5::Cint,arg6::PetscInt)
    ccall((:SNESLineSearchSetTolerances,petsc),PetscErrorCode,(SNESLineSearch,Cint,Cint,Cint,Cint,Cint,PetscInt),arg1,PetscReal,arg2,arg3,arg4,arg5,arg6)
end

function SNESLineSearchPreCheckPicard(arg1::SNESLineSearch,arg2::Vec,arg3::Vec,arg4::Ptr{PetscBool},arg5::Ptr{Void})
    ccall((:SNESLineSearchPreCheckPicard,petsc),PetscErrorCode,(SNESLineSearch,Vec,Vec,Ptr{PetscBool},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end

function SNESLineSearchGetLambda(arg1::SNESLineSearch,arg2::Ptr{Cint})
    ccall((:SNESLineSearchGetLambda,petsc),PetscErrorCode,(SNESLineSearch,Ptr{Cint}),arg1,arg2)
end

function SNESLineSearchSetLambda(arg1::SNESLineSearch,PetscReal::Cint)
    ccall((:SNESLineSearchSetLambda,petsc),PetscErrorCode,(SNESLineSearch,Cint),arg1,PetscReal)
end

function SNESLineSearchGetDamping(arg1::SNESLineSearch,arg2::Ptr{Cint})
    ccall((:SNESLineSearchGetDamping,petsc),PetscErrorCode,(SNESLineSearch,Ptr{Cint}),arg1,arg2)
end

function SNESLineSearchSetDamping(arg1::SNESLineSearch,PetscReal::Cint)
    ccall((:SNESLineSearchSetDamping,petsc),PetscErrorCode,(SNESLineSearch,Cint),arg1,PetscReal)
end

function SNESLineSearchGetOrder(arg1::SNESLineSearch,order::Ptr{PetscInt})
    ccall((:SNESLineSearchGetOrder,petsc),PetscErrorCode,(SNESLineSearch,Ptr{PetscInt}),arg1,order)
end

function SNESLineSearchSetOrder(arg1::SNESLineSearch,order::PetscInt)
    ccall((:SNESLineSearchSetOrder,petsc),PetscErrorCode,(SNESLineSearch,PetscInt),arg1,order)
end

function SNESLineSearchGetReason(arg1::SNESLineSearch,arg2::Ptr{SNESLineSearchReason})
    ccall((:SNESLineSearchGetReason,petsc),PetscErrorCode,(SNESLineSearch,Ptr{SNESLineSearchReason}),arg1,arg2)
end

function SNESLineSearchSetReason(arg1::SNESLineSearch,arg2::SNESLineSearchReason)
    ccall((:SNESLineSearchSetReason,petsc),PetscErrorCode,(SNESLineSearch,SNESLineSearchReason),arg1,arg2)
end

function SNESLineSearchGetVecs(arg1::SNESLineSearch,arg2::Ptr{Vec},arg3::Ptr{Vec},arg4::Ptr{Vec},arg5::Ptr{Vec},arg6::Ptr{Vec})
    ccall((:SNESLineSearchGetVecs,petsc),PetscErrorCode,(SNESLineSearch,Ptr{Vec},Ptr{Vec},Ptr{Vec},Ptr{Vec},Ptr{Vec}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function SNESLineSearchSetVecs(arg1::SNESLineSearch,arg2::Vec,arg3::Vec,arg4::Vec,arg5::Vec,arg6::Vec)
    ccall((:SNESLineSearchSetVecs,petsc),PetscErrorCode,(SNESLineSearch,Vec,Vec,Vec,Vec,Vec),arg1,arg2,arg3,arg4,arg5,arg6)
end

function SNESLineSearchGetNorms(arg1::SNESLineSearch,arg2::Ptr{Cint},arg3::Ptr{Cint},arg4::Ptr{Cint})
    ccall((:SNESLineSearchGetNorms,petsc),PetscErrorCode,(SNESLineSearch,Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function SNESLineSearchSetNorms(arg1::SNESLineSearch,PetscReal::Cint,arg2::Cint,arg3::Cint)
    ccall((:SNESLineSearchSetNorms,petsc),PetscErrorCode,(SNESLineSearch,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3)
end

function SNESLineSearchComputeNorms(arg1::SNESLineSearch)
    ccall((:SNESLineSearchComputeNorms,petsc),PetscErrorCode,(SNESLineSearch,),arg1)
end

function SNESLineSearchSetComputeNorms(arg1::SNESLineSearch,arg2::PetscBool)
    ccall((:SNESLineSearchSetComputeNorms,petsc),PetscErrorCode,(SNESLineSearch,PetscBool),arg1,arg2)
end

function SNESLineSearchSetMonitor(arg1::SNESLineSearch,arg2::PetscBool)
    ccall((:SNESLineSearchSetMonitor,petsc),PetscErrorCode,(SNESLineSearch,PetscBool),arg1,arg2)
end

function SNESLineSearchGetMonitor(arg1::SNESLineSearch,arg2::Ptr{PetscViewer})
    ccall((:SNESLineSearchGetMonitor,petsc),PetscErrorCode,(SNESLineSearch,Ptr{PetscViewer}),arg1,arg2)
end

function SNESLineSearchAppendOptionsPrefix(arg1::SNESLineSearch,prefix::Ptr{Uint8})
    ccall((:SNESLineSearchAppendOptionsPrefix,petsc),PetscErrorCode,(SNESLineSearch,Ptr{Uint8}),arg1,prefix)
end

function SNESLineSearchGetOptionsPrefix(arg1::SNESLineSearch,prefix::Ptr{Ptr{Uint8}})
    ccall((:SNESLineSearchGetOptionsPrefix,petsc),PetscErrorCode,(SNESLineSearch,Ptr{Ptr{Uint8}}),arg1,prefix)
end

function SNESLineSearchShellSetUserFunc(arg1::SNESLineSearch,arg2::SNESLineSearchUserFunc,arg3::Ptr{Void})
    ccall((:SNESLineSearchShellSetUserFunc,petsc),PetscErrorCode,(SNESLineSearch,SNESLineSearchUserFunc,Ptr{Void}),arg1,arg2,arg3)
end

function SNESLineSearchShellGetUserFunc(arg1::SNESLineSearch,arg2::Ptr{SNESLineSearchUserFunc},arg3::Ptr{Ptr{Void}})
    ccall((:SNESLineSearchShellGetUserFunc,petsc),PetscErrorCode,(SNESLineSearch,Ptr{SNESLineSearchUserFunc},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end

function SNESLineSearchBTSetAlpha(arg1::SNESLineSearch,PetscReal::Cint)
    ccall((:SNESLineSearchBTSetAlpha,petsc),PetscErrorCode,(SNESLineSearch,Cint),arg1,PetscReal)
end

function SNESLineSearchBTGetAlpha(arg1::SNESLineSearch,arg2::Ptr{Cint})
    ccall((:SNESLineSearchBTGetAlpha,petsc),PetscErrorCode,(SNESLineSearch,Ptr{Cint}),arg1,arg2)
end

function SNESLineSearchRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:SNESLineSearchRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function SNESVISetVariableBounds(arg1::SNES,arg2::Vec,arg3::Vec)
    ccall((:SNESVISetVariableBounds,petsc),PetscErrorCode,(SNES,Vec,Vec),arg1,arg2,arg3)
end

function SNESVISetComputeVariableBounds(arg1::SNES,arg2::Ptr{Void})
    ccall((:SNESVISetComputeVariableBounds,petsc),PetscErrorCode,(SNES,Ptr{Void}),arg1,arg2)
end

function SNESVIGetInactiveSet(arg1::SNES,arg2::Ptr{IS})
    ccall((:SNESVIGetInactiveSet,petsc),PetscErrorCode,(SNES,Ptr{IS}),arg1,arg2)
end

function SNESVIGetActiveSetIS(arg1::SNES,arg2::Vec,arg3::Vec,arg4::Ptr{IS})
    ccall((:SNESVIGetActiveSetIS,petsc),PetscErrorCode,(SNES,Vec,Vec,Ptr{IS}),arg1,arg2,arg3,arg4)
end

function SNESVIComputeInactiveSetFnorm(arg1::SNES,arg2::Vec,arg3::Vec,arg4::Ptr{Cint})
    ccall((:SNESVIComputeInactiveSetFnorm,petsc),PetscErrorCode,(SNES,Vec,Vec,Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function SNESVISetRedundancyCheck(arg1::SNES,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:SNESVISetRedundancyCheck,petsc),PetscErrorCode,(SNES,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function SNESTestLocalMin(arg1::SNES)
    ccall((:SNESTestLocalMin,petsc),PetscErrorCode,(SNES,),arg1)
end

function SNESComputeJacobian(arg1::SNES,arg2::Vec,arg3::Mat,arg4::Mat)
    ccall((:SNESComputeJacobian,petsc),PetscErrorCode,(SNES,Vec,Mat,Mat),arg1,arg2,arg3,arg4)
end

function SNESSetDM(arg1::SNES,arg2::DM)
    ccall((:SNESSetDM,petsc),PetscErrorCode,(SNES,DM),arg1,arg2)
end

function SNESGetDM(arg1::SNES,arg2::Ptr{DM})
    ccall((:SNESGetDM,petsc),PetscErrorCode,(SNES,Ptr{DM}),arg1,arg2)
end

function SNESSetNPC(arg1::SNES,arg2::SNES)
    ccall((:SNESSetNPC,petsc),PetscErrorCode,(SNES,SNES),arg1,arg2)
end

function SNESGetNPC(arg1::SNES,arg2::Ptr{SNES})
    ccall((:SNESGetNPC,petsc),PetscErrorCode,(SNES,Ptr{SNES}),arg1,arg2)
end

function SNESHasNPC(arg1::SNES,arg2::Ptr{PetscBool})
    ccall((:SNESHasNPC,petsc),PetscErrorCode,(SNES,Ptr{PetscBool}),arg1,arg2)
end

function SNESApplyNPC(arg1::SNES,arg2::Vec,arg3::Vec,arg4::Vec)
    ccall((:SNESApplyNPC,petsc),PetscErrorCode,(SNES,Vec,Vec,Vec),arg1,arg2,arg3,arg4)
end

function SNESGetNPCFunction(arg1::SNES,arg2::Vec,arg3::Ptr{Cint})
    ccall((:SNESGetNPCFunction,petsc),PetscErrorCode,(SNES,Vec,Ptr{Cint}),arg1,arg2,arg3)
end

function SNESComputeFunctionDefaultNPC(arg1::SNES,arg2::Vec,arg3::Vec)
    ccall((:SNESComputeFunctionDefaultNPC,petsc),PetscErrorCode,(SNES,Vec,Vec),arg1,arg2,arg3)
end

function SNESSetNPCSide(arg1::SNES,arg2::PCSide)
    ccall((:SNESSetNPCSide,petsc),PetscErrorCode,(SNES,PCSide),arg1,arg2)
end

function SNESGetNPCSide(arg1::SNES,arg2::Ptr{PCSide})
    ccall((:SNESGetNPCSide,petsc),PetscErrorCode,(SNES,Ptr{PCSide}),arg1,arg2)
end

function SNESSetLineSearch(arg1::SNES,arg2::SNESLineSearch)
    ccall((:SNESSetLineSearch,petsc),PetscErrorCode,(SNES,SNESLineSearch),arg1,arg2)
end

function SNESGetLineSearch(arg1::SNES,arg2::Ptr{SNESLineSearch})
    ccall((:SNESGetLineSearch,petsc),PetscErrorCode,(SNES,Ptr{SNESLineSearch}),arg1,arg2)
end

function SNESRestrictHookAdd(arg1::SNES,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:SNESRestrictHookAdd,petsc),PetscErrorCode,(SNES,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function SNESRestrictHooksRun(arg1::SNES,arg2::SNES)
    ccall((:SNESRestrictHooksRun,petsc),PetscErrorCode,(SNES,SNES),arg1,arg2)
end

function DMSNESSetFunction(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:DMSNESSetFunction,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function DMSNESGetFunction(arg1::DM,arg2::Ptr{Ptr{Void}},arg3::Ptr{Ptr{Void}})
    ccall((:DMSNESGetFunction,petsc),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end

function DMSNESSetNGS(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:DMSNESSetNGS,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function DMSNESGetNGS(arg1::DM,arg2::Ptr{Ptr{Void}},arg3::Ptr{Ptr{Void}})
    ccall((:DMSNESGetNGS,petsc),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end

function DMSNESSetJacobian(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:DMSNESSetJacobian,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function DMSNESGetJacobian(arg1::DM,arg2::Ptr{Ptr{Void}},arg3::Ptr{Ptr{Void}})
    ccall((:DMSNESGetJacobian,petsc),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end

function DMSNESSetPicard(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:DMSNESSetPicard,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function DMSNESGetPicard(arg1::DM,arg2::Ptr{Ptr{Void}},arg3::Ptr{Ptr{Void}},arg4::Ptr{Ptr{Void}})
    ccall((:DMSNESGetPicard,petsc),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4)
end

function DMSNESSetObjective(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:DMSNESSetObjective,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function DMSNESGetObjective(arg1::DM,arg2::Ptr{Ptr{Void}},arg3::Ptr{Ptr{Void}})
    ccall((:DMSNESGetObjective,petsc),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end

function DMDASNESSetFunctionLocal(arg1::DM,arg2::InsertMode,arg3::DMDASNESFunction,arg4::Ptr{Void})
    ccall((:DMDASNESSetFunctionLocal,petsc),PetscErrorCode,(DM,InsertMode,DMDASNESFunction,Ptr{Void}),arg1,arg2,arg3,arg4)
end

function DMDASNESSetJacobianLocal(arg1::DM,arg2::DMDASNESJacobian,arg3::Ptr{Void})
    ccall((:DMDASNESSetJacobianLocal,petsc),PetscErrorCode,(DM,DMDASNESJacobian,Ptr{Void}),arg1,arg2,arg3)
end

function DMDASNESSetObjectiveLocal(arg1::DM,arg2::DMDASNESObjective,arg3::Ptr{Void})
    ccall((:DMDASNESSetObjectiveLocal,petsc),PetscErrorCode,(DM,DMDASNESObjective,Ptr{Void}),arg1,arg2,arg3)
end

function DMDASNESSetPicardLocal(arg1::DM,arg2::InsertMode,arg3::Ptr{Void},arg4::Ptr{Void},arg5::Ptr{Void})
    ccall((:DMDASNESSetPicardLocal,petsc),PetscErrorCode,(DM,InsertMode,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end

function DMPlexSNESGetGeometryFEM(arg1::DM,arg2::Ptr{Vec})
    ccall((:DMPlexSNESGetGeometryFEM,petsc),PetscErrorCode,(DM,Ptr{Vec}),arg1,arg2)
end

function DMPlexSNESGetGeometryFVM(arg1::DM,arg2::Ptr{Vec},arg3::Ptr{Vec},arg4::Ptr{Cint})
    ccall((:DMPlexSNESGetGeometryFVM,petsc),PetscErrorCode,(DM,Ptr{Vec},Ptr{Vec},Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function DMPlexSNESGetGradientDM(arg1::DM,arg2::PetscFV,arg3::Ptr{DM})
    ccall((:DMPlexSNESGetGradientDM,petsc),PetscErrorCode,(DM,PetscFV,Ptr{DM}),arg1,arg2,arg3)
end

function DMPlexGetCellFields(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::Vec,arg5::Vec,arg6::Vec,arg7::Ptr{Ptr{PetscScalar}},arg8::Ptr{Ptr{PetscScalar}},arg9::Ptr{Ptr{PetscScalar}})
    ccall((:DMPlexGetCellFields,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,Vec,Vec,Vec,Ptr{Ptr{PetscScalar}},Ptr{Ptr{PetscScalar}},Ptr{Ptr{PetscScalar}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function DMPlexRestoreCellFields(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::Vec,arg5::Vec,arg6::Vec,arg7::Ptr{Ptr{PetscScalar}},arg8::Ptr{Ptr{PetscScalar}},arg9::Ptr{Ptr{PetscScalar}})
    ccall((:DMPlexRestoreCellFields,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,Vec,Vec,Vec,Ptr{Ptr{PetscScalar}},Ptr{Ptr{PetscScalar}},Ptr{Ptr{PetscScalar}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function DMPlexGetFaceFields(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::Vec,arg5::Vec,arg6::Vec,arg7::Vec,arg8::Vec,arg9::Ptr{Ptr{PetscScalar}},arg10::Ptr{Ptr{PetscScalar}})
    ccall((:DMPlexGetFaceFields,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,Vec,Vec,Vec,Vec,Vec,Ptr{Ptr{PetscScalar}},Ptr{Ptr{PetscScalar}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function DMPlexRestoreFaceFields(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::Vec,arg5::Vec,arg6::Vec,arg7::Vec,arg8::Vec,arg9::Ptr{Ptr{PetscScalar}},arg10::Ptr{Ptr{PetscScalar}})
    ccall((:DMPlexRestoreFaceFields,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,Vec,Vec,Vec,Vec,Vec,Ptr{Ptr{PetscScalar}},Ptr{Ptr{PetscScalar}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function DMPlexGetFaceGeometry(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::Vec,arg5::Vec,arg6::Ptr{Ptr{PetscFVFaceGeom}},arg7::Ptr{Ptr{Cint}})
    ccall((:DMPlexGetFaceGeometry,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,Vec,Vec,Ptr{Ptr{PetscFVFaceGeom}},Ptr{Ptr{Cint}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function DMPlexRestoreFaceGeometry(arg1::DM,arg2::PetscInt,arg3::PetscInt,arg4::Vec,arg5::Vec,arg6::Ptr{Ptr{PetscFVFaceGeom}},arg7::Ptr{Ptr{Cint}})
    ccall((:DMPlexRestoreFaceGeometry,petsc),PetscErrorCode,(DM,PetscInt,PetscInt,Vec,Vec,Ptr{Ptr{PetscFVFaceGeom}},Ptr{Ptr{Cint}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function DMSNESSetFunctionLocal(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:DMSNESSetFunctionLocal,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function DMSNESSetJacobianLocal(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:DMSNESSetJacobianLocal,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function SNESMultiblockSetFields(arg1::SNES,arg2::Ptr{Uint8},arg3::PetscInt,arg4::Ptr{PetscInt})
    ccall((:SNESMultiblockSetFields,petsc),PetscErrorCode,(SNES,Ptr{Uint8},PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function SNESMultiblockSetIS(arg1::SNES,arg2::Ptr{Uint8},arg3::IS)
    ccall((:SNESMultiblockSetIS,petsc),PetscErrorCode,(SNES,Ptr{Uint8},IS),arg1,arg2,arg3)
end

function SNESMultiblockSetBlockSize(arg1::SNES,arg2::PetscInt)
    ccall((:SNESMultiblockSetBlockSize,petsc),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end

function SNESMultiblockSetType(arg1::SNES,arg2::PCCompositeType)
    ccall((:SNESMultiblockSetType,petsc),PetscErrorCode,(SNES,PCCompositeType),arg1,arg2)
end

function SNESMSRegister(arg1::SNESMSType,arg2::PetscInt,arg3::PetscInt,PetscReal::Cint,arg4::Ptr{Cint},arg5::Ptr{Cint},arg6::Ptr{Cint})
    ccall((:SNESMSRegister,petsc),PetscErrorCode,(SNESMSType,PetscInt,PetscInt,Cint,Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,PetscReal,arg4,arg5,arg6)
end

function SNESMSSetType(arg1::SNES,arg2::SNESMSType)
    ccall((:SNESMSSetType,petsc),PetscErrorCode,(SNES,SNESMSType),arg1,arg2)
end

function SNESMSFinalizePackage()
    ccall((:SNESMSFinalizePackage,petsc),PetscErrorCode,())
end

function SNESMSInitializePackage()
    ccall((:SNESMSInitializePackage,petsc),PetscErrorCode,())
end

function SNESMSRegisterDestroy()
    ccall((:SNESMSRegisterDestroy,petsc),PetscErrorCode,())
end

function SNESNGMRESSetRestartType(arg1::SNES,arg2::SNESNGMRESRestartType)
    ccall((:SNESNGMRESSetRestartType,petsc),PetscErrorCode,(SNES,SNESNGMRESRestartType),arg1,arg2)
end

function SNESNGMRESSetSelectType(arg1::SNES,arg2::SNESNGMRESSelectType)
    ccall((:SNESNGMRESSetSelectType,petsc),PetscErrorCode,(SNES,SNESNGMRESSelectType),arg1,arg2)
end

function SNESNCGSetType(arg1::SNES,arg2::SNESNCGType)
    ccall((:SNESNCGSetType,petsc),PetscErrorCode,(SNES,SNESNCGType),arg1,arg2)
end

function SNESQNSetType(arg1::SNES,arg2::SNESQNType)
    ccall((:SNESQNSetType,petsc),PetscErrorCode,(SNES,SNESQNType),arg1,arg2)
end

function SNESQNSetScaleType(arg1::SNES,arg2::SNESQNScaleType)
    ccall((:SNESQNSetScaleType,petsc),PetscErrorCode,(SNES,SNESQNScaleType),arg1,arg2)
end

function SNESQNSetRestartType(arg1::SNES,arg2::SNESQNRestartType)
    ccall((:SNESQNSetRestartType,petsc),PetscErrorCode,(SNES,SNESQNRestartType),arg1,arg2)
end

function SNESNASMGetType(arg1::SNES,arg2::Ptr{PCASMType})
    ccall((:SNESNASMGetType,petsc),PetscErrorCode,(SNES,Ptr{PCASMType}),arg1,arg2)
end

function SNESNASMSetType(arg1::SNES,arg2::PCASMType)
    ccall((:SNESNASMSetType,petsc),PetscErrorCode,(SNES,PCASMType),arg1,arg2)
end

function SNESNASMGetSubdomains(arg1::SNES,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{SNES}},arg4::Ptr{Ptr{VecScatter}},arg5::Ptr{Ptr{VecScatter}},arg6::Ptr{Ptr{VecScatter}})
    ccall((:SNESNASMGetSubdomains,petsc),PetscErrorCode,(SNES,Ptr{PetscInt},Ptr{Ptr{SNES}},Ptr{Ptr{VecScatter}},Ptr{Ptr{VecScatter}},Ptr{Ptr{VecScatter}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function SNESNASMSetSubdomains(arg1::SNES,arg2::PetscInt,arg3::Ptr{SNES},arg4::Ptr{VecScatter},arg5::Ptr{VecScatter},arg6::Ptr{VecScatter})
    ccall((:SNESNASMSetSubdomains,petsc),PetscErrorCode,(SNES,PetscInt,Ptr{SNES},Ptr{VecScatter},Ptr{VecScatter},Ptr{VecScatter}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function SNESNASMSetDamping(arg1::SNES,PetscReal::Cint)
    ccall((:SNESNASMSetDamping,petsc),PetscErrorCode,(SNES,Cint),arg1,PetscReal)
end

function SNESNASMGetDamping(arg1::SNES,arg2::Ptr{Cint})
    ccall((:SNESNASMGetDamping,petsc),PetscErrorCode,(SNES,Ptr{Cint}),arg1,arg2)
end

function SNESNASMGetSubdomainVecs(arg1::SNES,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{Vec}},arg4::Ptr{Ptr{Vec}},arg5::Ptr{Ptr{Vec}},arg6::Ptr{Ptr{Vec}})
    ccall((:SNESNASMGetSubdomainVecs,petsc),PetscErrorCode,(SNES,Ptr{PetscInt},Ptr{Ptr{Vec}},Ptr{Ptr{Vec}},Ptr{Ptr{Vec}},Ptr{Ptr{Vec}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function SNESNASMSetComputeFinalJacobian(arg1::SNES,arg2::PetscBool)
    ccall((:SNESNASMSetComputeFinalJacobian,petsc),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end

function SNESCompositeSetType(arg1::SNES,arg2::SNESCompositeType)
    ccall((:SNESCompositeSetType,petsc),PetscErrorCode,(SNES,SNESCompositeType),arg1,arg2)
end

function SNESCompositeAddSNES(arg1::SNES,arg2::SNESType)
    ccall((:SNESCompositeAddSNES,petsc),PetscErrorCode,(SNES,SNESType),arg1,arg2)
end

function SNESCompositeGetSNES(arg1::SNES,arg2::PetscInt,arg3::Ptr{SNES})
    ccall((:SNESCompositeGetSNES,petsc),PetscErrorCode,(SNES,PetscInt,Ptr{SNES}),arg1,arg2,arg3)
end

function SNESCompositeGetNumber(arg1::SNES,arg2::Ptr{PetscInt})
    ccall((:SNESCompositeGetNumber,petsc),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end

function SNESCompositeSetDamping(arg1::SNES,arg2::PetscInt,PetscReal::Cint)
    ccall((:SNESCompositeSetDamping,petsc),PetscErrorCode,(SNES,PetscInt,Cint),arg1,arg2,PetscReal)
end

function SNESFASSetType(arg1::SNES,arg2::SNESFASType)
    ccall((:SNESFASSetType,petsc),PetscErrorCode,(SNES,SNESFASType),arg1,arg2)
end

function SNESFASGetType(arg1::SNES,arg2::Ptr{SNESFASType})
    ccall((:SNESFASGetType,petsc),PetscErrorCode,(SNES,Ptr{SNESFASType}),arg1,arg2)
end

function SNESFASSetLevels(arg1::SNES,arg2::PetscInt,arg3::Ptr{MPI_Comm})
    ccall((:SNESFASSetLevels,petsc),PetscErrorCode,(SNES,PetscInt,Ptr{MPI_Comm}),arg1,arg2,arg3)
end

function SNESFASGetLevels(arg1::SNES,arg2::Ptr{PetscInt})
    ccall((:SNESFASGetLevels,petsc),PetscErrorCode,(SNES,Ptr{PetscInt}),arg1,arg2)
end

function SNESFASGetCycleSNES(arg1::SNES,arg2::PetscInt,arg3::Ptr{SNES})
    ccall((:SNESFASGetCycleSNES,petsc),PetscErrorCode,(SNES,PetscInt,Ptr{SNES}),arg1,arg2,arg3)
end

function SNESFASSetNumberSmoothUp(arg1::SNES,arg2::PetscInt)
    ccall((:SNESFASSetNumberSmoothUp,petsc),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end

function SNESFASSetNumberSmoothDown(arg1::SNES,arg2::PetscInt)
    ccall((:SNESFASSetNumberSmoothDown,petsc),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end

function SNESFASSetCycles(arg1::SNES,arg2::PetscInt)
    ccall((:SNESFASSetCycles,petsc),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end

function SNESFASSetMonitor(arg1::SNES,arg2::PetscBool)
    ccall((:SNESFASSetMonitor,petsc),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end

function SNESFASSetLog(arg1::SNES,arg2::PetscBool)
    ccall((:SNESFASSetLog,petsc),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end

function SNESFASSetGalerkin(arg1::SNES,arg2::PetscBool)
    ccall((:SNESFASSetGalerkin,petsc),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end

function SNESFASGetGalerkin(arg1::SNES,arg2::Ptr{PetscBool})
    ccall((:SNESFASGetGalerkin,petsc),PetscErrorCode,(SNES,Ptr{PetscBool}),arg1,arg2)
end

function SNESFASCycleGetSmoother(arg1::SNES,arg2::Ptr{SNES})
    ccall((:SNESFASCycleGetSmoother,petsc),PetscErrorCode,(SNES,Ptr{SNES}),arg1,arg2)
end

function SNESFASCycleGetSmootherUp(arg1::SNES,arg2::Ptr{SNES})
    ccall((:SNESFASCycleGetSmootherUp,petsc),PetscErrorCode,(SNES,Ptr{SNES}),arg1,arg2)
end

function SNESFASCycleGetSmootherDown(arg1::SNES,arg2::Ptr{SNES})
    ccall((:SNESFASCycleGetSmootherDown,petsc),PetscErrorCode,(SNES,Ptr{SNES}),arg1,arg2)
end

function SNESFASCycleGetCorrection(arg1::SNES,arg2::Ptr{SNES})
    ccall((:SNESFASCycleGetCorrection,petsc),PetscErrorCode,(SNES,Ptr{SNES}),arg1,arg2)
end

function SNESFASCycleGetInterpolation(arg1::SNES,arg2::Ptr{Mat})
    ccall((:SNESFASCycleGetInterpolation,petsc),PetscErrorCode,(SNES,Ptr{Mat}),arg1,arg2)
end

function SNESFASCycleGetRestriction(arg1::SNES,arg2::Ptr{Mat})
    ccall((:SNESFASCycleGetRestriction,petsc),PetscErrorCode,(SNES,Ptr{Mat}),arg1,arg2)
end

function SNESFASCycleGetInjection(arg1::SNES,arg2::Ptr{Mat})
    ccall((:SNESFASCycleGetInjection,petsc),PetscErrorCode,(SNES,Ptr{Mat}),arg1,arg2)
end

function SNESFASCycleGetRScale(arg1::SNES,arg2::Ptr{Vec})
    ccall((:SNESFASCycleGetRScale,petsc),PetscErrorCode,(SNES,Ptr{Vec}),arg1,arg2)
end

function SNESFASCycleSetCycles(arg1::SNES,arg2::PetscInt)
    ccall((:SNESFASCycleSetCycles,petsc),PetscErrorCode,(SNES,PetscInt),arg1,arg2)
end

function SNESFASCycleIsFine(arg1::SNES,arg2::Ptr{PetscBool})
    ccall((:SNESFASCycleIsFine,petsc),PetscErrorCode,(SNES,Ptr{PetscBool}),arg1,arg2)
end

function SNESFASSetInterpolation(arg1::SNES,arg2::PetscInt,arg3::Mat)
    ccall((:SNESFASSetInterpolation,petsc),PetscErrorCode,(SNES,PetscInt,Mat),arg1,arg2,arg3)
end

function SNESFASGetInterpolation(arg1::SNES,arg2::PetscInt,arg3::Ptr{Mat})
    ccall((:SNESFASGetInterpolation,petsc),PetscErrorCode,(SNES,PetscInt,Ptr{Mat}),arg1,arg2,arg3)
end

function SNESFASSetRestriction(arg1::SNES,arg2::PetscInt,arg3::Mat)
    ccall((:SNESFASSetRestriction,petsc),PetscErrorCode,(SNES,PetscInt,Mat),arg1,arg2,arg3)
end

function SNESFASGetRestriction(arg1::SNES,arg2::PetscInt,arg3::Ptr{Mat})
    ccall((:SNESFASGetRestriction,petsc),PetscErrorCode,(SNES,PetscInt,Ptr{Mat}),arg1,arg2,arg3)
end

function SNESFASSetInjection(arg1::SNES,arg2::PetscInt,arg3::Mat)
    ccall((:SNESFASSetInjection,petsc),PetscErrorCode,(SNES,PetscInt,Mat),arg1,arg2,arg3)
end

function SNESFASGetInjection(arg1::SNES,arg2::PetscInt,arg3::Ptr{Mat})
    ccall((:SNESFASGetInjection,petsc),PetscErrorCode,(SNES,PetscInt,Ptr{Mat}),arg1,arg2,arg3)
end

function SNESFASSetRScale(arg1::SNES,arg2::PetscInt,arg3::Vec)
    ccall((:SNESFASSetRScale,petsc),PetscErrorCode,(SNES,PetscInt,Vec),arg1,arg2,arg3)
end

function SNESFASGetRScale(arg1::SNES,arg2::PetscInt,arg3::Ptr{Vec})
    ccall((:SNESFASGetRScale,petsc),PetscErrorCode,(SNES,PetscInt,Ptr{Vec}),arg1,arg2,arg3)
end

function SNESFASSetContinuation(arg1::SNES,arg2::PetscBool)
    ccall((:SNESFASSetContinuation,petsc),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end

function SNESFASGetSmoother(arg1::SNES,arg2::PetscInt,arg3::Ptr{SNES})
    ccall((:SNESFASGetSmoother,petsc),PetscErrorCode,(SNES,PetscInt,Ptr{SNES}),arg1,arg2,arg3)
end

function SNESFASGetSmootherUp(arg1::SNES,arg2::PetscInt,arg3::Ptr{SNES})
    ccall((:SNESFASGetSmootherUp,petsc),PetscErrorCode,(SNES,PetscInt,Ptr{SNES}),arg1,arg2,arg3)
end

function SNESFASGetSmootherDown(arg1::SNES,arg2::PetscInt,arg3::Ptr{SNES})
    ccall((:SNESFASGetSmootherDown,petsc),PetscErrorCode,(SNES,PetscInt,Ptr{SNES}),arg1,arg2,arg3)
end

function SNESFASGetCoarseSolve(arg1::SNES,arg2::Ptr{SNES})
    ccall((:SNESFASGetCoarseSolve,petsc),PetscErrorCode,(SNES,Ptr{SNES}),arg1,arg2)
end

function SNESFASFullSetDownSweep(arg1::SNES,arg2::PetscBool)
    ccall((:SNESFASFullSetDownSweep,petsc),PetscErrorCode,(SNES,PetscBool),arg1,arg2)
end

function SNESFASCreateCoarseVec(arg1::SNES,arg2::Ptr{Vec})
    ccall((:SNESFASCreateCoarseVec,petsc),PetscErrorCode,(SNES,Ptr{Vec}),arg1,arg2)
end

function SNESFASRestrict(arg1::SNES,arg2::Vec,arg3::Vec)
    ccall((:SNESFASRestrict,petsc),PetscErrorCode,(SNES,Vec,Vec),arg1,arg2,arg3)
end

function DMSNESCheckFromOptions(arg1::SNES,arg2::Vec,arg3::Ptr{Ptr{Void}},arg4::Ptr{Ptr{Void}})
    ccall((:DMSNESCheckFromOptions,petsc),PetscErrorCode,(SNES,Vec,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4)
end

function TSInitializePackage()
    ccall((:TSInitializePackage,petsc),PetscErrorCode,())
end

function TSCreate(arg1::MPI_Comm,arg2::Ptr{TS})
    ccall((:TSCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{TS}),arg1,arg2)
end

function TSClone(arg1::TS,arg2::Ptr{TS})
    ccall((:TSClone,petsc),PetscErrorCode,(TS,Ptr{TS}),arg1,arg2)
end

function TSDestroy(arg1::Ptr{TS})
    ccall((:TSDestroy,petsc),PetscErrorCode,(Ptr{TS},),arg1)
end

function TSSetProblemType(arg1::TS,arg2::TSProblemType)
    ccall((:TSSetProblemType,petsc),PetscErrorCode,(TS,TSProblemType),arg1,arg2)
end

function TSGetProblemType(arg1::TS,arg2::Ptr{TSProblemType})
    ccall((:TSGetProblemType,petsc),PetscErrorCode,(TS,Ptr{TSProblemType}),arg1,arg2)
end

function TSMonitor(arg1::TS,arg2::PetscInt,PetscReal::Cint,arg3::Vec)
    ccall((:TSMonitor,petsc),PetscErrorCode,(TS,PetscInt,Cint,Vec),arg1,arg2,PetscReal,arg3)
end

function TSMonitorSet(arg1::TS,arg2::Ptr{Void},arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:TSMonitorSet,petsc),PetscErrorCode,(TS,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function TSMonitorCancel(arg1::TS)
    ccall((:TSMonitorCancel,petsc),PetscErrorCode,(TS,),arg1)
end

function TSSetOptionsPrefix(arg1::TS,arg2::Ptr{Uint8})
    ccall((:TSSetOptionsPrefix,petsc),PetscErrorCode,(TS,Ptr{Uint8}),arg1,arg2)
end

function TSAppendOptionsPrefix(arg1::TS,arg2::Ptr{Uint8})
    ccall((:TSAppendOptionsPrefix,petsc),PetscErrorCode,(TS,Ptr{Uint8}),arg1,arg2)
end

function TSGetOptionsPrefix(arg1::TS,arg2::Ptr{Ptr{Uint8}})
    ccall((:TSGetOptionsPrefix,petsc),PetscErrorCode,(TS,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function TSSetFromOptions(arg1::TS)
    ccall((:TSSetFromOptions,petsc),PetscErrorCode,(TS,),arg1)
end

function TSSetUp(arg1::TS)
    ccall((:TSSetUp,petsc),PetscErrorCode,(TS,),arg1)
end

function TSReset(arg1::TS)
    ccall((:TSReset,petsc),PetscErrorCode,(TS,),arg1)
end

function TSSetSolution(arg1::TS,arg2::Vec)
    ccall((:TSSetSolution,petsc),PetscErrorCode,(TS,Vec),arg1,arg2)
end

function TSGetSolution(arg1::TS,arg2::Ptr{Vec})
    ccall((:TSGetSolution,petsc),PetscErrorCode,(TS,Ptr{Vec}),arg1,arg2)
end

function TSSetSaveTrajectory(arg1::TS)
    ccall((:TSSetSaveTrajectory,petsc),PetscErrorCode,(TS,),arg1)
end

function TSTrajectoryCreate(arg1::MPI_Comm,arg2::Ptr{TSTrajectory})
    ccall((:TSTrajectoryCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{TSTrajectory}),arg1,arg2)
end

function TSTrajectoryDestroy(arg1::Ptr{TSTrajectory})
    ccall((:TSTrajectoryDestroy,petsc),PetscErrorCode,(Ptr{TSTrajectory},),arg1)
end

function TSTrajectorySetType(arg1::TSTrajectory,arg2::TSTrajectoryType)
    ccall((:TSTrajectorySetType,petsc),PetscErrorCode,(TSTrajectory,TSTrajectoryType),arg1,arg2)
end

function TSTrajectorySet(arg1::TSTrajectory,arg2::TS,arg3::PetscInt,PetscReal::Cint,arg4::Vec)
    ccall((:TSTrajectorySet,petsc),PetscErrorCode,(TSTrajectory,TS,PetscInt,Cint,Vec),arg1,arg2,arg3,PetscReal,arg4)
end

function TSTrajectoryGet(arg1::TSTrajectory,arg2::TS,arg3::PetscInt,arg4::Ptr{Cint})
    ccall((:TSTrajectoryGet,petsc),PetscErrorCode,(TSTrajectory,TS,PetscInt,Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function TSTrajectorySetFromOptions(arg1::TSTrajectory)
    ccall((:TSTrajectorySetFromOptions,petsc),PetscErrorCode,(TSTrajectory,),arg1)
end

function TSTrajectoryRegisterAll()
    ccall((:TSTrajectoryRegisterAll,petsc),PetscErrorCode,())
end

function TSSetCostGradients(arg1::TS,arg2::PetscInt,arg3::Ptr{Vec},arg4::Ptr{Vec})
    ccall((:TSSetCostGradients,petsc),PetscErrorCode,(TS,PetscInt,Ptr{Vec},Ptr{Vec}),arg1,arg2,arg3,arg4)
end

function TSGetCostGradients(arg1::TS,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{Vec}},arg4::Ptr{Ptr{Vec}})
    ccall((:TSGetCostGradients,petsc),PetscErrorCode,(TS,Ptr{PetscInt},Ptr{Ptr{Vec}},Ptr{Ptr{Vec}}),arg1,arg2,arg3,arg4)
end

function TSSetCostIntegrand(arg1::TS,arg2::PetscInt,arg3::Ptr{Void},arg4::Ptr{Void},arg5::Ptr{Void},arg6::Ptr{Void})
    ccall((:TSSetCostIntegrand,petsc),PetscErrorCode,(TS,PetscInt,Ptr{Void},Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function TSGetCostIntegral(arg1::TS,arg2::Ptr{Vec})
    ccall((:TSGetCostIntegral,petsc),PetscErrorCode,(TS,Ptr{Vec}),arg1,arg2)
end

function TSAdjointSetRHSJacobian(arg1::TS,arg2::Mat,arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:TSAdjointSetRHSJacobian,petsc),PetscErrorCode,(TS,Mat,Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function TSAdjointSolve(arg1::TS)
    ccall((:TSAdjointSolve,petsc),PetscErrorCode,(TS,),arg1)
end

function TSAdjointSetSteps(arg1::TS,arg2::PetscInt)
    ccall((:TSAdjointSetSteps,petsc),PetscErrorCode,(TS,PetscInt),arg1,arg2)
end

function TSAdjointComputeRHSJacobian(arg1::TS,PetscReal::Cint,arg2::Vec,arg3::Mat)
    ccall((:TSAdjointComputeRHSJacobian,petsc),PetscErrorCode,(TS,Cint,Vec,Mat),arg1,PetscReal,arg2,arg3)
end

function TSAdjointStep(arg1::TS)
    ccall((:TSAdjointStep,petsc),PetscErrorCode,(TS,),arg1)
end

function TSAdjointSetUp(arg1::TS)
    ccall((:TSAdjointSetUp,petsc),PetscErrorCode,(TS,),arg1)
end

function TSAdjointComputeDRDPFunction(arg1::TS,PetscReal::Cint,arg2::Vec,arg3::Ptr{Vec})
    ccall((:TSAdjointComputeDRDPFunction,petsc),PetscErrorCode,(TS,Cint,Vec,Ptr{Vec}),arg1,PetscReal,arg2,arg3)
end

function TSAdjointComputeDRDYFunction(arg1::TS,PetscReal::Cint,arg2::Vec,arg3::Ptr{Vec})
    ccall((:TSAdjointComputeDRDYFunction,petsc),PetscErrorCode,(TS,Cint,Vec,Ptr{Vec}),arg1,PetscReal,arg2,arg3)
end

function TSAdjointComputeCostIntegrand(arg1::TS,PetscReal::Cint,arg2::Vec,arg3::Vec)
    ccall((:TSAdjointComputeCostIntegrand,petsc),PetscErrorCode,(TS,Cint,Vec,Vec),arg1,PetscReal,arg2,arg3)
end

function TSSetDuration(arg1::TS,arg2::PetscInt,PetscReal::Cint)
    ccall((:TSSetDuration,petsc),PetscErrorCode,(TS,PetscInt,Cint),arg1,arg2,PetscReal)
end

function TSGetDuration(arg1::TS,arg2::Ptr{PetscInt},arg3::Ptr{Cint})
    ccall((:TSGetDuration,petsc),PetscErrorCode,(TS,Ptr{PetscInt},Ptr{Cint}),arg1,arg2,arg3)
end

function TSSetExactFinalTime(arg1::TS,arg2::TSExactFinalTimeOption)
    ccall((:TSSetExactFinalTime,petsc),PetscErrorCode,(TS,TSExactFinalTimeOption),arg1,arg2)
end

function TSMonitorDefault(arg1::TS,arg2::PetscInt,PetscReal::Cint,arg3::Vec,arg4::Ptr{Void})
    ccall((:TSMonitorDefault,petsc),PetscErrorCode,(TS,PetscInt,Cint,Vec,Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end

function TSMonitorDrawCtxCreate(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Cint,arg5::Cint,arg6::Cint,arg7::Cint,arg8::PetscInt,arg9::Ptr{TSMonitorDrawCtx})
    ccall((:TSMonitorDrawCtxCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Cint,Cint,Cint,Cint,PetscInt,Ptr{TSMonitorDrawCtx}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function TSMonitorDrawCtxDestroy(arg1::Ptr{TSMonitorDrawCtx})
    ccall((:TSMonitorDrawCtxDestroy,petsc),PetscErrorCode,(Ptr{TSMonitorDrawCtx},),arg1)
end

function TSMonitorDrawSolution(arg1::TS,arg2::PetscInt,PetscReal::Cint,arg3::Vec,arg4::Ptr{Void})
    ccall((:TSMonitorDrawSolution,petsc),PetscErrorCode,(TS,PetscInt,Cint,Vec,Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end

function TSMonitorDrawSolutionPhase(arg1::TS,arg2::PetscInt,PetscReal::Cint,arg3::Vec,arg4::Ptr{Void})
    ccall((:TSMonitorDrawSolutionPhase,petsc),PetscErrorCode,(TS,PetscInt,Cint,Vec,Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end

function TSMonitorDrawError(arg1::TS,arg2::PetscInt,PetscReal::Cint,arg3::Vec,arg4::Ptr{Void})
    ccall((:TSMonitorDrawError,petsc),PetscErrorCode,(TS,PetscInt,Cint,Vec,Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end

function TSMonitorSolutionBinary(arg1::TS,arg2::PetscInt,PetscReal::Cint,arg3::Vec,arg4::Ptr{Void})
    ccall((:TSMonitorSolutionBinary,petsc),PetscErrorCode,(TS,PetscInt,Cint,Vec,Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end

function TSMonitorSolutionVTK(arg1::TS,arg2::PetscInt,PetscReal::Cint,arg3::Vec,arg4::Ptr{Void})
    ccall((:TSMonitorSolutionVTK,petsc),PetscErrorCode,(TS,PetscInt,Cint,Vec,Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end

function TSMonitorSolutionVTKDestroy(arg1::Ptr{Void})
    ccall((:TSMonitorSolutionVTKDestroy,petsc),PetscErrorCode,(Ptr{Void},),arg1)
end

function TSStep(arg1::TS)
    ccall((:TSStep,petsc),PetscErrorCode,(TS,),arg1)
end

function TSEvaluateStep(arg1::TS,arg2::PetscInt,arg3::Vec,arg4::Ptr{PetscBool})
    ccall((:TSEvaluateStep,petsc),PetscErrorCode,(TS,PetscInt,Vec,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function TSSolve(arg1::TS,arg2::Vec)
    ccall((:TSSolve,petsc),PetscErrorCode,(TS,Vec),arg1,arg2)
end

function TSGetEquationType(arg1::TS,arg2::Ptr{TSEquationType})
    ccall((:TSGetEquationType,petsc),PetscErrorCode,(TS,Ptr{TSEquationType}),arg1,arg2)
end

function TSSetEquationType(arg1::TS,arg2::TSEquationType)
    ccall((:TSSetEquationType,petsc),PetscErrorCode,(TS,TSEquationType),arg1,arg2)
end

function TSGetConvergedReason(arg1::TS,arg2::Ptr{TSConvergedReason})
    ccall((:TSGetConvergedReason,petsc),PetscErrorCode,(TS,Ptr{TSConvergedReason}),arg1,arg2)
end

function TSSetConvergedReason(arg1::TS,arg2::TSConvergedReason)
    ccall((:TSSetConvergedReason,petsc),PetscErrorCode,(TS,TSConvergedReason),arg1,arg2)
end

function TSGetSolveTime(arg1::TS,arg2::Ptr{Cint})
    ccall((:TSGetSolveTime,petsc),PetscErrorCode,(TS,Ptr{Cint}),arg1,arg2)
end

function TSGetSNESIterations(arg1::TS,arg2::Ptr{PetscInt})
    ccall((:TSGetSNESIterations,petsc),PetscErrorCode,(TS,Ptr{PetscInt}),arg1,arg2)
end

function TSGetKSPIterations(arg1::TS,arg2::Ptr{PetscInt})
    ccall((:TSGetKSPIterations,petsc),PetscErrorCode,(TS,Ptr{PetscInt}),arg1,arg2)
end

function TSGetStepRejections(arg1::TS,arg2::Ptr{PetscInt})
    ccall((:TSGetStepRejections,petsc),PetscErrorCode,(TS,Ptr{PetscInt}),arg1,arg2)
end

function TSSetMaxStepRejections(arg1::TS,arg2::PetscInt)
    ccall((:TSSetMaxStepRejections,petsc),PetscErrorCode,(TS,PetscInt),arg1,arg2)
end

function TSGetSNESFailures(arg1::TS,arg2::Ptr{PetscInt})
    ccall((:TSGetSNESFailures,petsc),PetscErrorCode,(TS,Ptr{PetscInt}),arg1,arg2)
end

function TSSetMaxSNESFailures(arg1::TS,arg2::PetscInt)
    ccall((:TSSetMaxSNESFailures,petsc),PetscErrorCode,(TS,PetscInt),arg1,arg2)
end

function TSSetErrorIfStepFails(arg1::TS,arg2::PetscBool)
    ccall((:TSSetErrorIfStepFails,petsc),PetscErrorCode,(TS,PetscBool),arg1,arg2)
end

function TSRollBack(arg1::TS)
    ccall((:TSRollBack,petsc),PetscErrorCode,(TS,),arg1)
end

function TSGetTotalSteps(arg1::TS,arg2::Ptr{PetscInt})
    ccall((:TSGetTotalSteps,petsc),PetscErrorCode,(TS,Ptr{PetscInt}),arg1,arg2)
end

function TSGetStages(arg1::TS,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{Vec}})
    ccall((:TSGetStages,petsc),PetscErrorCode,(TS,Ptr{PetscInt},Ptr{Ptr{Vec}}),arg1,arg2,arg3)
end

function TSSetInitialTimeStep(arg1::TS,PetscReal::Cint,arg2::Cint)
    ccall((:TSSetInitialTimeStep,petsc),PetscErrorCode,(TS,Cint,Cint),arg1,PetscReal,arg2)
end

function TSGetTimeStep(arg1::TS,arg2::Ptr{Cint})
    ccall((:TSGetTimeStep,petsc),PetscErrorCode,(TS,Ptr{Cint}),arg1,arg2)
end

function TSGetTime(arg1::TS,arg2::Ptr{Cint})
    ccall((:TSGetTime,petsc),PetscErrorCode,(TS,Ptr{Cint}),arg1,arg2)
end

function TSSetTime(arg1::TS,PetscReal::Cint)
    ccall((:TSSetTime,petsc),PetscErrorCode,(TS,Cint),arg1,PetscReal)
end

function TSGetTimeStepNumber(arg1::TS,arg2::Ptr{PetscInt})
    ccall((:TSGetTimeStepNumber,petsc),PetscErrorCode,(TS,Ptr{PetscInt}),arg1,arg2)
end

function TSSetTimeStep(arg1::TS,PetscReal::Cint)
    ccall((:TSSetTimeStep,petsc),PetscErrorCode,(TS,Cint),arg1,PetscReal)
end

function TSGetPrevTime(arg1::TS,arg2::Ptr{Cint})
    ccall((:TSGetPrevTime,petsc),PetscErrorCode,(TS,Ptr{Cint}),arg1,arg2)
end

function TSSetRHSFunction(arg1::TS,arg2::Vec,arg3::TSRHSFunction,arg4::Ptr{Void})
    ccall((:TSSetRHSFunction,petsc),PetscErrorCode,(TS,Vec,TSRHSFunction,Ptr{Void}),arg1,arg2,arg3,arg4)
end

function TSGetRHSFunction(arg1::TS,arg2::Ptr{Vec},arg3::Ptr{TSRHSFunction},arg4::Ptr{Ptr{Void}})
    ccall((:TSGetRHSFunction,petsc),PetscErrorCode,(TS,Ptr{Vec},Ptr{TSRHSFunction},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4)
end

function TSSetRHSJacobian(arg1::TS,arg2::Mat,arg3::Mat,arg4::TSRHSJacobian,arg5::Ptr{Void})
    ccall((:TSSetRHSJacobian,petsc),PetscErrorCode,(TS,Mat,Mat,TSRHSJacobian,Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end

function TSGetRHSJacobian(arg1::TS,arg2::Ptr{Mat},arg3::Ptr{Mat},arg4::Ptr{TSRHSJacobian},arg5::Ptr{Ptr{Void}})
    ccall((:TSGetRHSJacobian,petsc),PetscErrorCode,(TS,Ptr{Mat},Ptr{Mat},Ptr{TSRHSJacobian},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4,arg5)
end

function TSRHSJacobianSetReuse(arg1::TS,arg2::PetscBool)
    ccall((:TSRHSJacobianSetReuse,petsc),PetscErrorCode,(TS,PetscBool),arg1,arg2)
end

function TSSetSolutionFunction(arg1::TS,arg2::TSSolutionFunction,arg3::Ptr{Void})
    ccall((:TSSetSolutionFunction,petsc),PetscErrorCode,(TS,TSSolutionFunction,Ptr{Void}),arg1,arg2,arg3)
end

function TSSetForcingFunction(arg1::TS,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:TSSetForcingFunction,petsc),PetscErrorCode,(TS,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function TSSetIFunction(arg1::TS,arg2::Vec,arg3::TSIFunction,arg4::Ptr{Void})
    ccall((:TSSetIFunction,petsc),PetscErrorCode,(TS,Vec,TSIFunction,Ptr{Void}),arg1,arg2,arg3,arg4)
end

function TSGetIFunction(arg1::TS,arg2::Ptr{Vec},arg3::Ptr{TSIFunction},arg4::Ptr{Ptr{Void}})
    ccall((:TSGetIFunction,petsc),PetscErrorCode,(TS,Ptr{Vec},Ptr{TSIFunction},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4)
end

function TSSetIJacobian(arg1::TS,arg2::Mat,arg3::Mat,arg4::TSIJacobian,arg5::Ptr{Void})
    ccall((:TSSetIJacobian,petsc),PetscErrorCode,(TS,Mat,Mat,TSIJacobian,Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end

function TSGetIJacobian(arg1::TS,arg2::Ptr{Mat},arg3::Ptr{Mat},arg4::Ptr{TSIJacobian},arg5::Ptr{Ptr{Void}})
    ccall((:TSGetIJacobian,petsc),PetscErrorCode,(TS,Ptr{Mat},Ptr{Mat},Ptr{TSIJacobian},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4,arg5)
end

function TSComputeRHSFunctionLinear(arg1::TS,PetscReal::Cint,arg2::Vec,arg3::Vec,arg4::Ptr{Void})
    ccall((:TSComputeRHSFunctionLinear,petsc),PetscErrorCode,(TS,Cint,Vec,Vec,Ptr{Void}),arg1,PetscReal,arg2,arg3,arg4)
end

function TSComputeRHSJacobianConstant(arg1::TS,PetscReal::Cint,arg2::Vec,arg3::Mat,arg4::Mat,arg5::Ptr{Void})
    ccall((:TSComputeRHSJacobianConstant,petsc),PetscErrorCode,(TS,Cint,Vec,Mat,Mat,Ptr{Void}),arg1,PetscReal,arg2,arg3,arg4,arg5)
end

function TSComputeIFunctionLinear(arg1::TS,PetscReal::Cint,arg2::Vec,arg3::Vec,arg4::Vec,arg5::Ptr{Void})
    ccall((:TSComputeIFunctionLinear,petsc),PetscErrorCode,(TS,Cint,Vec,Vec,Vec,Ptr{Void}),arg1,PetscReal,arg2,arg3,arg4,arg5)
end

function TSComputeIJacobianConstant(arg1::TS,PetscReal::Cint,arg2::Vec,arg3::Vec,arg4::Cint,arg5::Mat,arg6::Mat,arg7::Ptr{Void})
    ccall((:TSComputeIJacobianConstant,petsc),PetscErrorCode,(TS,Cint,Vec,Vec,Cint,Mat,Mat,Ptr{Void}),arg1,PetscReal,arg2,arg3,arg4,arg5,arg6,arg7)
end

function TSComputeSolutionFunction(arg1::TS,PetscReal::Cint,arg2::Vec)
    ccall((:TSComputeSolutionFunction,petsc),PetscErrorCode,(TS,Cint,Vec),arg1,PetscReal,arg2)
end

function TSComputeForcingFunction(arg1::TS,PetscReal::Cint,arg2::Vec)
    ccall((:TSComputeForcingFunction,petsc),PetscErrorCode,(TS,Cint,Vec),arg1,PetscReal,arg2)
end

function TSComputeIJacobianDefaultColor(arg1::TS,PetscReal::Cint,arg2::Vec,arg3::Vec,arg4::Cint,arg5::Mat,arg6::Mat,arg7::Ptr{Void})
    ccall((:TSComputeIJacobianDefaultColor,petsc),PetscErrorCode,(TS,Cint,Vec,Vec,Cint,Mat,Mat,Ptr{Void}),arg1,PetscReal,arg2,arg3,arg4,arg5,arg6,arg7)
end

function TSSetPreStep(arg1::TS,arg2::Ptr{Void})
    ccall((:TSSetPreStep,petsc),PetscErrorCode,(TS,Ptr{Void}),arg1,arg2)
end

function TSSetPreStage(arg1::TS,arg2::Ptr{Void})
    ccall((:TSSetPreStage,petsc),PetscErrorCode,(TS,Ptr{Void}),arg1,arg2)
end

function TSSetPostStage(arg1::TS,arg2::Ptr{Void})
    ccall((:TSSetPostStage,petsc),PetscErrorCode,(TS,Ptr{Void}),arg1,arg2)
end

function TSSetPostStep(arg1::TS,arg2::Ptr{Void})
    ccall((:TSSetPostStep,petsc),PetscErrorCode,(TS,Ptr{Void}),arg1,arg2)
end

function TSPreStep(arg1::TS)
    ccall((:TSPreStep,petsc),PetscErrorCode,(TS,),arg1)
end

function TSPreStage(arg1::TS,PetscReal::Cint)
    ccall((:TSPreStage,petsc),PetscErrorCode,(TS,Cint),arg1,PetscReal)
end

function TSPostStage(arg1::TS,PetscReal::Cint,arg2::PetscInt,arg3::Ptr{Vec})
    ccall((:TSPostStage,petsc),PetscErrorCode,(TS,Cint,PetscInt,Ptr{Vec}),arg1,PetscReal,arg2,arg3)
end

function TSPostStep(arg1::TS)
    ccall((:TSPostStep,petsc),PetscErrorCode,(TS,),arg1)
end

function TSSetRetainStages(arg1::TS,arg2::PetscBool)
    ccall((:TSSetRetainStages,petsc),PetscErrorCode,(TS,PetscBool),arg1,arg2)
end

function TSInterpolate(arg1::TS,PetscReal::Cint,arg2::Vec)
    ccall((:TSInterpolate,petsc),PetscErrorCode,(TS,Cint,Vec),arg1,PetscReal,arg2)
end

function TSSetTolerances(arg1::TS,PetscReal::Cint,arg2::Vec,arg3::Cint,arg4::Vec)
    ccall((:TSSetTolerances,petsc),PetscErrorCode,(TS,Cint,Vec,Cint,Vec),arg1,PetscReal,arg2,arg3,arg4)
end

function TSGetTolerances(arg1::TS,arg2::Ptr{Cint},arg3::Ptr{Vec},arg4::Ptr{Cint},arg5::Ptr{Vec})
    ccall((:TSGetTolerances,petsc),PetscErrorCode,(TS,Ptr{Cint},Ptr{Vec},Ptr{Cint},Ptr{Vec}),arg1,arg2,arg3,arg4,arg5)
end

function TSErrorWeightedNormInfinity(arg1::TS,arg2::Vec,arg3::Vec,arg4::Ptr{Cint})
    ccall((:TSErrorWeightedNormInfinity,petsc),PetscErrorCode,(TS,Vec,Vec,Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function TSErrorWeightedNorm2(arg1::TS,arg2::Vec,arg3::Vec,arg4::Ptr{Cint})
    ccall((:TSErrorWeightedNorm2,petsc),PetscErrorCode,(TS,Vec,Vec,Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function TSErrorWeightedNorm(arg1::TS,arg2::Vec,arg3::Vec,arg4::NormType,arg5::Ptr{Cint})
    ccall((:TSErrorWeightedNorm,petsc),PetscErrorCode,(TS,Vec,Vec,NormType,Ptr{Cint}),arg1,arg2,arg3,arg4,arg5)
end

function TSSetCFLTimeLocal(arg1::TS,PetscReal::Cint)
    ccall((:TSSetCFLTimeLocal,petsc),PetscErrorCode,(TS,Cint),arg1,PetscReal)
end

function TSGetCFLTime(arg1::TS,arg2::Ptr{Cint})
    ccall((:TSGetCFLTime,petsc),PetscErrorCode,(TS,Ptr{Cint}),arg1,arg2)
end

function TSPseudoSetTimeStep(arg1::TS,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:TSPseudoSetTimeStep,petsc),PetscErrorCode,(TS,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function TSPseudoTimeStepDefault(arg1::TS,arg2::Ptr{Cint},arg3::Ptr{Void})
    ccall((:TSPseudoTimeStepDefault,petsc),PetscErrorCode,(TS,Ptr{Cint},Ptr{Void}),arg1,arg2,arg3)
end

function TSPseudoComputeTimeStep(arg1::TS,arg2::Ptr{Cint})
    ccall((:TSPseudoComputeTimeStep,petsc),PetscErrorCode,(TS,Ptr{Cint}),arg1,arg2)
end

function TSPseudoSetMaxTimeStep(arg1::TS,PetscReal::Cint)
    ccall((:TSPseudoSetMaxTimeStep,petsc),PetscErrorCode,(TS,Cint),arg1,PetscReal)
end

function TSPseudoSetVerifyTimeStep(arg1::TS,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:TSPseudoSetVerifyTimeStep,petsc),PetscErrorCode,(TS,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function TSPseudoVerifyTimeStepDefault(arg1::TS,arg2::Vec,arg3::Ptr{Void},arg4::Ptr{Cint},arg5::Ptr{PetscBool})
    ccall((:TSPseudoVerifyTimeStepDefault,petsc),PetscErrorCode,(TS,Vec,Ptr{Void},Ptr{Cint},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function TSPseudoVerifyTimeStep(arg1::TS,arg2::Vec,arg3::Ptr{Cint},arg4::Ptr{PetscBool})
    ccall((:TSPseudoVerifyTimeStep,petsc),PetscErrorCode,(TS,Vec,Ptr{Cint},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function TSPseudoSetTimeStepIncrement(arg1::TS,PetscReal::Cint)
    ccall((:TSPseudoSetTimeStepIncrement,petsc),PetscErrorCode,(TS,Cint),arg1,PetscReal)
end

function TSPseudoIncrementDtFromInitialDt(arg1::TS)
    ccall((:TSPseudoIncrementDtFromInitialDt,petsc),PetscErrorCode,(TS,),arg1)
end

function TSPythonSetType(arg1::TS,arg2::Ptr{Uint8})
    ccall((:TSPythonSetType,petsc),PetscErrorCode,(TS,Ptr{Uint8}),arg1,arg2)
end

function TSComputeRHSFunction(arg1::TS,PetscReal::Cint,arg2::Vec,arg3::Vec)
    ccall((:TSComputeRHSFunction,petsc),PetscErrorCode,(TS,Cint,Vec,Vec),arg1,PetscReal,arg2,arg3)
end

function TSComputeRHSJacobian(arg1::TS,PetscReal::Cint,arg2::Vec,arg3::Mat,arg4::Mat)
    ccall((:TSComputeRHSJacobian,petsc),PetscErrorCode,(TS,Cint,Vec,Mat,Mat),arg1,PetscReal,arg2,arg3,arg4)
end

function TSComputeIFunction(arg1::TS,PetscReal::Cint,arg2::Vec,arg3::Vec,arg4::Vec,arg5::PetscBool)
    ccall((:TSComputeIFunction,petsc),PetscErrorCode,(TS,Cint,Vec,Vec,Vec,PetscBool),arg1,PetscReal,arg2,arg3,arg4,arg5)
end

function TSComputeIJacobian(arg1::TS,PetscReal::Cint,arg2::Vec,arg3::Vec,arg4::Cint,arg5::Mat,arg6::Mat,arg7::PetscBool)
    ccall((:TSComputeIJacobian,petsc),PetscErrorCode,(TS,Cint,Vec,Vec,Cint,Mat,Mat,PetscBool),arg1,PetscReal,arg2,arg3,arg4,arg5,arg6,arg7)
end

function TSComputeLinearStability(arg1::TS,PetscReal::Cint,arg2::Cint,arg3::Ptr{Cint},arg4::Ptr{Cint})
    ccall((:TSComputeLinearStability,petsc),PetscErrorCode,(TS,Cint,Cint,Ptr{Cint},Ptr{Cint}),arg1,PetscReal,arg2,arg3,arg4)
end

function TSVISetVariableBounds(arg1::TS,arg2::Vec,arg3::Vec)
    ccall((:TSVISetVariableBounds,petsc),PetscErrorCode,(TS,Vec,Vec),arg1,arg2,arg3)
end

function DMTSSetRHSFunction(arg1::DM,arg2::TSRHSFunction,arg3::Ptr{Void})
    ccall((:DMTSSetRHSFunction,petsc),PetscErrorCode,(DM,TSRHSFunction,Ptr{Void}),arg1,arg2,arg3)
end

function DMTSGetRHSFunction(arg1::DM,arg2::Ptr{TSRHSFunction},arg3::Ptr{Ptr{Void}})
    ccall((:DMTSGetRHSFunction,petsc),PetscErrorCode,(DM,Ptr{TSRHSFunction},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end

function DMTSSetRHSJacobian(arg1::DM,arg2::TSRHSJacobian,arg3::Ptr{Void})
    ccall((:DMTSSetRHSJacobian,petsc),PetscErrorCode,(DM,TSRHSJacobian,Ptr{Void}),arg1,arg2,arg3)
end

function DMTSGetRHSJacobian(arg1::DM,arg2::Ptr{TSRHSJacobian},arg3::Ptr{Ptr{Void}})
    ccall((:DMTSGetRHSJacobian,petsc),PetscErrorCode,(DM,Ptr{TSRHSJacobian},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end

function DMTSSetIFunction(arg1::DM,arg2::TSIFunction,arg3::Ptr{Void})
    ccall((:DMTSSetIFunction,petsc),PetscErrorCode,(DM,TSIFunction,Ptr{Void}),arg1,arg2,arg3)
end

function DMTSGetIFunction(arg1::DM,arg2::Ptr{TSIFunction},arg3::Ptr{Ptr{Void}})
    ccall((:DMTSGetIFunction,petsc),PetscErrorCode,(DM,Ptr{TSIFunction},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end

function DMTSSetIJacobian(arg1::DM,arg2::TSIJacobian,arg3::Ptr{Void})
    ccall((:DMTSSetIJacobian,petsc),PetscErrorCode,(DM,TSIJacobian,Ptr{Void}),arg1,arg2,arg3)
end

function DMTSGetIJacobian(arg1::DM,arg2::Ptr{TSIJacobian},arg3::Ptr{Ptr{Void}})
    ccall((:DMTSGetIJacobian,petsc),PetscErrorCode,(DM,Ptr{TSIJacobian},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end

function DMTSSetSolutionFunction(arg1::DM,arg2::TSSolutionFunction,arg3::Ptr{Void})
    ccall((:DMTSSetSolutionFunction,petsc),PetscErrorCode,(DM,TSSolutionFunction,Ptr{Void}),arg1,arg2,arg3)
end

function DMTSGetSolutionFunction(arg1::DM,arg2::Ptr{TSSolutionFunction},arg3::Ptr{Ptr{Void}})
    ccall((:DMTSGetSolutionFunction,petsc),PetscErrorCode,(DM,Ptr{TSSolutionFunction},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end

function DMTSSetForcingFunction(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:DMTSSetForcingFunction,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function DMTSGetForcingFunction(arg1::DM,arg2::Ptr{Ptr{Void}},arg3::Ptr{Ptr{Void}})
    ccall((:DMTSGetForcingFunction,petsc),PetscErrorCode,(DM,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3)
end

function DMTSGetMinRadius(arg1::DM,arg2::Ptr{Cint})
    ccall((:DMTSGetMinRadius,petsc),PetscErrorCode,(DM,Ptr{Cint}),arg1,arg2)
end

function DMTSSetMinRadius(arg1::DM,PetscReal::Cint)
    ccall((:DMTSSetMinRadius,petsc),PetscErrorCode,(DM,Cint),arg1,PetscReal)
end

function DMTSCheckFromOptions(arg1::TS,arg2::Vec,arg3::Ptr{Ptr{Void}},arg4::Ptr{Ptr{Void}})
    ccall((:DMTSCheckFromOptions,petsc),PetscErrorCode,(TS,Vec,Ptr{Ptr{Void}},Ptr{Ptr{Void}}),arg1,arg2,arg3,arg4)
end

function DMTSSetIFunctionLocal(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:DMTSSetIFunctionLocal,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function DMTSSetIJacobianLocal(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:DMTSSetIJacobianLocal,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function DMTSSetRHSFunctionLocal(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:DMTSSetRHSFunctionLocal,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function DMTSSetIFunctionSerialize(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:DMTSSetIFunctionSerialize,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function DMTSSetIJacobianSerialize(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:DMTSSetIJacobianSerialize,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function DMDATSSetRHSFunctionLocal(arg1::DM,arg2::InsertMode,arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:DMDATSSetRHSFunctionLocal,petsc),PetscErrorCode,(DM,InsertMode,Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function DMDATSSetRHSJacobianLocal(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:DMDATSSetRHSJacobianLocal,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function DMDATSSetIFunctionLocal(arg1::DM,arg2::InsertMode,arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:DMDATSSetIFunctionLocal,petsc),PetscErrorCode,(DM,InsertMode,Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function DMDATSSetIJacobianLocal(arg1::DM,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:DMDATSSetIJacobianLocal,petsc),PetscErrorCode,(DM,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function DMPlexTSGetGeometryFVM(arg1::DM,arg2::Ptr{Vec},arg3::Ptr{Vec},arg4::Ptr{Cint})
    ccall((:DMPlexTSGetGeometryFVM,petsc),PetscErrorCode,(DM,Ptr{Vec},Ptr{Vec},Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function TSMonitorDMDARayDestroy(arg1::Ptr{Ptr{Void}})
    ccall((:TSMonitorDMDARayDestroy,petsc),PetscErrorCode,(Ptr{Ptr{Void}},),arg1)
end

function TSMonitorDMDARay(arg1::TS,arg2::PetscInt,PetscReal::Cint,arg3::Vec,arg4::Ptr{Void})
    ccall((:TSMonitorDMDARay,petsc),PetscErrorCode,(TS,PetscInt,Cint,Vec,Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end

function TSMonitorLGDMDARay(arg1::TS,arg2::PetscInt,PetscReal::Cint,arg3::Vec,arg4::Ptr{Void})
    ccall((:TSMonitorLGDMDARay,petsc),PetscErrorCode,(TS,PetscInt,Cint,Vec,Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end

function TSGetType(arg1::TS,arg2::Ptr{TSType})
    ccall((:TSGetType,petsc),PetscErrorCode,(TS,Ptr{TSType}),arg1,arg2)
end

function TSSetType(arg1::TS,arg2::TSType)
    ccall((:TSSetType,petsc),PetscErrorCode,(TS,TSType),arg1,arg2)
end

function TSRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:TSRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function TSGetSNES(arg1::TS,arg2::Ptr{SNES})
    ccall((:TSGetSNES,petsc),PetscErrorCode,(TS,Ptr{SNES}),arg1,arg2)
end

function TSSetSNES(arg1::TS,arg2::SNES)
    ccall((:TSSetSNES,petsc),PetscErrorCode,(TS,SNES),arg1,arg2)
end

function TSGetKSP(arg1::TS,arg2::Ptr{KSP})
    ccall((:TSGetKSP,petsc),PetscErrorCode,(TS,Ptr{KSP}),arg1,arg2)
end

function TSView(arg1::TS,arg2::PetscViewer)
    ccall((:TSView,petsc),PetscErrorCode,(TS,PetscViewer),arg1,arg2)
end

function TSLoad(arg1::TS,arg2::PetscViewer)
    ccall((:TSLoad,petsc),PetscErrorCode,(TS,PetscViewer),arg1,arg2)
end

function TSGetApplicationContext(arg1::TS,arg2::Ptr{Void})
    ccall((:TSGetApplicationContext,petsc),PetscErrorCode,(TS,Ptr{Void}),arg1,arg2)
end

function TSMonitorLGCtxCreate(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Cint,arg5::Cint,arg6::Cint,arg7::Cint,arg8::PetscInt,arg9::Ptr{TSMonitorLGCtx})
    ccall((:TSMonitorLGCtxCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Cint,Cint,Cint,Cint,PetscInt,Ptr{TSMonitorLGCtx}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function TSMonitorLGCtxDestroy(arg1::Ptr{TSMonitorLGCtx})
    ccall((:TSMonitorLGCtxDestroy,petsc),PetscErrorCode,(Ptr{TSMonitorLGCtx},),arg1)
end

function TSMonitorLGTimeStep(arg1::TS,arg2::PetscInt,PetscReal::Cint,arg3::Vec,arg4::Ptr{Void})
    ccall((:TSMonitorLGTimeStep,petsc),PetscErrorCode,(TS,PetscInt,Cint,Vec,Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end

function TSMonitorLGSolution(arg1::TS,arg2::PetscInt,PetscReal::Cint,arg3::Vec,arg4::Ptr{Void})
    ccall((:TSMonitorLGSolution,petsc),PetscErrorCode,(TS,PetscInt,Cint,Vec,Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end

function TSMonitorLGSetVariableNames(arg1::TS,arg2::Ptr{Ptr{Uint8}})
    ccall((:TSMonitorLGSetVariableNames,petsc),PetscErrorCode,(TS,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function TSMonitorLGGetVariableNames(arg1::TS,arg2::Ptr{Ptr{Ptr{Uint8}}})
    ccall((:TSMonitorLGGetVariableNames,petsc),PetscErrorCode,(TS,Ptr{Ptr{Ptr{Uint8}}}),arg1,arg2)
end

function TSMonitorLGCtxSetVariableNames(arg1::TSMonitorLGCtx,arg2::Ptr{Ptr{Uint8}})
    ccall((:TSMonitorLGCtxSetVariableNames,petsc),PetscErrorCode,(TSMonitorLGCtx,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function TSMonitorLGSetDisplayVariables(arg1::TS,arg2::Ptr{Ptr{Uint8}})
    ccall((:TSMonitorLGSetDisplayVariables,petsc),PetscErrorCode,(TS,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function TSMonitorLGCtxSetDisplayVariables(arg1::TSMonitorLGCtx,arg2::Ptr{Ptr{Uint8}})
    ccall((:TSMonitorLGCtxSetDisplayVariables,petsc),PetscErrorCode,(TSMonitorLGCtx,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function TSMonitorLGSetTransform(arg1::TS,arg2::Ptr{Void},arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:TSMonitorLGSetTransform,petsc),PetscErrorCode,(TS,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function TSMonitorLGCtxSetTransform(arg1::TSMonitorLGCtx,arg2::Ptr{Void},arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:TSMonitorLGCtxSetTransform,petsc),PetscErrorCode,(TSMonitorLGCtx,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function TSMonitorLGError(arg1::TS,arg2::PetscInt,PetscReal::Cint,arg3::Vec,arg4::Ptr{Void})
    ccall((:TSMonitorLGError,petsc),PetscErrorCode,(TS,PetscInt,Cint,Vec,Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end

function TSMonitorLGSNESIterations(arg1::TS,arg2::PetscInt,PetscReal::Cint,arg3::Vec,arg4::Ptr{Void})
    ccall((:TSMonitorLGSNESIterations,petsc),PetscErrorCode,(TS,PetscInt,Cint,Vec,Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end

function TSMonitorLGKSPIterations(arg1::TS,arg2::PetscInt,PetscReal::Cint,arg3::Vec,arg4::Ptr{Void})
    ccall((:TSMonitorLGKSPIterations,petsc),PetscErrorCode,(TS,PetscInt,Cint,Vec,Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end

function TSMonitorEnvelopeCtxCreate(arg1::TS,arg2::Ptr{TSMonitorEnvelopeCtx})
    ccall((:TSMonitorEnvelopeCtxCreate,petsc),PetscErrorCode,(TS,Ptr{TSMonitorEnvelopeCtx}),arg1,arg2)
end

function TSMonitorEnvelope(arg1::TS,arg2::PetscInt,PetscReal::Cint,arg3::Vec,arg4::Ptr{Void})
    ccall((:TSMonitorEnvelope,petsc),PetscErrorCode,(TS,PetscInt,Cint,Vec,Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end

function TSMonitorEnvelopeGetBounds(arg1::TS,arg2::Ptr{Vec},arg3::Ptr{Vec})
    ccall((:TSMonitorEnvelopeGetBounds,petsc),PetscErrorCode,(TS,Ptr{Vec},Ptr{Vec}),arg1,arg2,arg3)
end

function TSMonitorEnvelopeCtxDestroy(arg1::Ptr{TSMonitorEnvelopeCtx})
    ccall((:TSMonitorEnvelopeCtxDestroy,petsc),PetscErrorCode,(Ptr{TSMonitorEnvelopeCtx},),arg1)
end

function TSMonitorSPEigCtxCreate(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Cint,arg5::Cint,arg6::Cint,arg7::Cint,arg8::PetscInt,arg9::Ptr{TSMonitorSPEigCtx})
    ccall((:TSMonitorSPEigCtxCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Cint,Cint,Cint,Cint,PetscInt,Ptr{TSMonitorSPEigCtx}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function TSMonitorSPEigCtxDestroy(arg1::Ptr{TSMonitorSPEigCtx})
    ccall((:TSMonitorSPEigCtxDestroy,petsc),PetscErrorCode,(Ptr{TSMonitorSPEigCtx},),arg1)
end

function TSMonitorSPEig(arg1::TS,arg2::PetscInt,PetscReal::Cint,arg3::Vec,arg4::Ptr{Void})
    ccall((:TSMonitorSPEig,petsc),PetscErrorCode,(TS,PetscInt,Cint,Vec,Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end

function TSSetEventMonitor(arg1::TS,arg2::PetscInt,arg3::Ptr{PetscInt},arg4::Ptr{PetscBool},arg5::Ptr{Void},arg6::Ptr{Void},arg7::Ptr{Void})
    ccall((:TSSetEventMonitor,petsc),PetscErrorCode,(TS,PetscInt,Ptr{PetscInt},Ptr{PetscBool},Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function TSSetEventTolerances(arg1::TS,PetscReal::Cint,arg2::Ptr{Cint})
    ccall((:TSSetEventTolerances,petsc),PetscErrorCode,(TS,Cint,Ptr{Cint}),arg1,PetscReal,arg2)
end

function TSSSPSetType(arg1::TS,arg2::TSSSPType)
    ccall((:TSSSPSetType,petsc),PetscErrorCode,(TS,TSSSPType),arg1,arg2)
end

function TSSSPGetType(arg1::TS,arg2::Ptr{TSSSPType})
    ccall((:TSSSPGetType,petsc),PetscErrorCode,(TS,Ptr{TSSSPType}),arg1,arg2)
end

function TSSSPSetNumStages(arg1::TS,arg2::PetscInt)
    ccall((:TSSSPSetNumStages,petsc),PetscErrorCode,(TS,PetscInt),arg1,arg2)
end

function TSSSPGetNumStages(arg1::TS,arg2::Ptr{PetscInt})
    ccall((:TSSSPGetNumStages,petsc),PetscErrorCode,(TS,Ptr{PetscInt}),arg1,arg2)
end

function TSSSPFinalizePackage()
    ccall((:TSSSPFinalizePackage,petsc),PetscErrorCode,())
end

function TSSSPInitializePackage()
    ccall((:TSSSPInitializePackage,petsc),PetscErrorCode,())
end

function TSGetAdapt(arg1::TS,arg2::Ptr{TSAdapt})
    ccall((:TSGetAdapt,petsc),PetscErrorCode,(TS,Ptr{TSAdapt}),arg1,arg2)
end

function TSAdaptRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:TSAdaptRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function TSAdaptInitializePackage()
    ccall((:TSAdaptInitializePackage,petsc),PetscErrorCode,())
end

function TSAdaptFinalizePackage()
    ccall((:TSAdaptFinalizePackage,petsc),PetscErrorCode,())
end

function TSAdaptCreate(arg1::MPI_Comm,arg2::Ptr{TSAdapt})
    ccall((:TSAdaptCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{TSAdapt}),arg1,arg2)
end

function TSAdaptSetType(arg1::TSAdapt,arg2::TSAdaptType)
    ccall((:TSAdaptSetType,petsc),PetscErrorCode,(TSAdapt,TSAdaptType),arg1,arg2)
end

function TSAdaptSetOptionsPrefix(arg1::TSAdapt,arg2::Ptr{Uint8})
    ccall((:TSAdaptSetOptionsPrefix,petsc),PetscErrorCode,(TSAdapt,Ptr{Uint8}),arg1,arg2)
end

function TSAdaptCandidatesClear(arg1::TSAdapt)
    ccall((:TSAdaptCandidatesClear,petsc),PetscErrorCode,(TSAdapt,),arg1)
end

function TSAdaptCandidateAdd(arg1::TSAdapt,arg2::Ptr{Uint8},arg3::PetscInt,arg4::PetscInt,PetscReal::Cint,arg5::Cint,arg6::PetscBool)
    ccall((:TSAdaptCandidateAdd,petsc),PetscErrorCode,(TSAdapt,Ptr{Uint8},PetscInt,PetscInt,Cint,Cint,PetscBool),arg1,arg2,arg3,arg4,PetscReal,arg5,arg6)
end

function TSAdaptCandidatesGet(arg1::TSAdapt,arg2::Ptr{PetscInt},arg3::Ptr{Ptr{PetscInt}},arg4::Ptr{Ptr{PetscInt}},arg5::Ptr{Ptr{Cint}},arg6::Ptr{Ptr{Cint}})
    ccall((:TSAdaptCandidatesGet,petsc),PetscErrorCode,(TSAdapt,Ptr{PetscInt},Ptr{Ptr{PetscInt}},Ptr{Ptr{PetscInt}},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function TSAdaptChoose(arg1::TSAdapt,arg2::TS,PetscReal::Cint,arg3::Ptr{PetscInt},arg4::Ptr{Cint},arg5::Ptr{PetscBool})
    ccall((:TSAdaptChoose,petsc),PetscErrorCode,(TSAdapt,TS,Cint,Ptr{PetscInt},Ptr{Cint},Ptr{PetscBool}),arg1,arg2,PetscReal,arg3,arg4,arg5)
end

function TSAdaptCheckStage(arg1::TSAdapt,arg2::TS,arg3::Ptr{PetscBool})
    ccall((:TSAdaptCheckStage,petsc),PetscErrorCode,(TSAdapt,TS,Ptr{PetscBool}),arg1,arg2,arg3)
end

function TSAdaptView(arg1::TSAdapt,arg2::PetscViewer)
    ccall((:TSAdaptView,petsc),PetscErrorCode,(TSAdapt,PetscViewer),arg1,arg2)
end

function TSAdaptLoad(arg1::TSAdapt,arg2::PetscViewer)
    ccall((:TSAdaptLoad,petsc),PetscErrorCode,(TSAdapt,PetscViewer),arg1,arg2)
end

function TSAdaptSetFromOptions(arg1::Ptr{PetscOptions},arg2::TSAdapt)
    ccall((:TSAdaptSetFromOptions,petsc),PetscErrorCode,(Ptr{PetscOptions},TSAdapt),arg1,arg2)
end

function TSAdaptReset(arg1::TSAdapt)
    ccall((:TSAdaptReset,petsc),PetscErrorCode,(TSAdapt,),arg1)
end

function TSAdaptDestroy(arg1::Ptr{TSAdapt})
    ccall((:TSAdaptDestroy,petsc),PetscErrorCode,(Ptr{TSAdapt},),arg1)
end

function TSAdaptSetMonitor(arg1::TSAdapt,arg2::PetscBool)
    ccall((:TSAdaptSetMonitor,petsc),PetscErrorCode,(TSAdapt,PetscBool),arg1,arg2)
end

function TSAdaptSetStepLimits(arg1::TSAdapt,PetscReal::Cint,arg2::Cint)
    ccall((:TSAdaptSetStepLimits,petsc),PetscErrorCode,(TSAdapt,Cint,Cint),arg1,PetscReal,arg2)
end

function TSAdaptSetCheckStage(arg1::TSAdapt,arg2::Ptr{Void})
    ccall((:TSAdaptSetCheckStage,petsc),PetscErrorCode,(TSAdapt,Ptr{Void}),arg1,arg2)
end

function TSGLAdaptRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:TSGLAdaptRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function TSGLAdaptInitializePackage()
    ccall((:TSGLAdaptInitializePackage,petsc),PetscErrorCode,())
end

function TSGLAdaptFinalizePackage()
    ccall((:TSGLAdaptFinalizePackage,petsc),PetscErrorCode,())
end

function TSGLAdaptCreate(arg1::MPI_Comm,arg2::Ptr{TSGLAdapt})
    ccall((:TSGLAdaptCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{TSGLAdapt}),arg1,arg2)
end

function TSGLAdaptSetType(arg1::TSGLAdapt,arg2::TSGLAdaptType)
    ccall((:TSGLAdaptSetType,petsc),PetscErrorCode,(TSGLAdapt,TSGLAdaptType),arg1,arg2)
end

function TSGLAdaptSetOptionsPrefix(arg1::TSGLAdapt,arg2::Ptr{Uint8})
    ccall((:TSGLAdaptSetOptionsPrefix,petsc),PetscErrorCode,(TSGLAdapt,Ptr{Uint8}),arg1,arg2)
end

function TSGLAdaptChoose(arg1::TSGLAdapt,arg2::PetscInt,arg3::Ptr{PetscInt},PetscReal::Ptr{Cint},arg4::Ptr{Cint},arg5::PetscInt,arg6::Cint,arg7::Cint,arg8::Ptr{PetscInt},arg9::Ptr{Cint},arg10::Ptr{PetscBool})
    ccall((:TSGLAdaptChoose,petsc),PetscErrorCode,(TSGLAdapt,PetscInt,Ptr{PetscInt},Ptr{Cint},Ptr{Cint},PetscInt,Cint,Cint,Ptr{PetscInt},Ptr{Cint},Ptr{PetscBool}),arg1,arg2,arg3,PetscReal,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function TSGLAdaptView(arg1::TSGLAdapt,arg2::PetscViewer)
    ccall((:TSGLAdaptView,petsc),PetscErrorCode,(TSGLAdapt,PetscViewer),arg1,arg2)
end

function TSGLAdaptSetFromOptions(arg1::Ptr{PetscOptions},arg2::TSGLAdapt)
    ccall((:TSGLAdaptSetFromOptions,petsc),PetscErrorCode,(Ptr{PetscOptions},TSGLAdapt),arg1,arg2)
end

function TSGLAdaptDestroy(arg1::Ptr{TSGLAdapt})
    ccall((:TSGLAdaptDestroy,petsc),PetscErrorCode,(Ptr{TSGLAdapt},),arg1)
end

function TSGLAcceptRegister(arg1::Ptr{Uint8},arg2::TSGLAcceptFunction)
    ccall((:TSGLAcceptRegister,petsc),PetscErrorCode,(Ptr{Uint8},TSGLAcceptFunction),arg1,arg2)
end

function TSGLRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:TSGLRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function TSGLInitializePackage()
    ccall((:TSGLInitializePackage,petsc),PetscErrorCode,())
end

function TSGLFinalizePackage()
    ccall((:TSGLFinalizePackage,petsc),PetscErrorCode,())
end

function TSGLSetType(arg1::TS,arg2::TSGLType)
    ccall((:TSGLSetType,petsc),PetscErrorCode,(TS,TSGLType),arg1,arg2)
end

function TSGLGetAdapt(arg1::TS,arg2::Ptr{TSGLAdapt})
    ccall((:TSGLGetAdapt,petsc),PetscErrorCode,(TS,Ptr{TSGLAdapt}),arg1,arg2)
end

function TSGLSetAcceptType(arg1::TS,arg2::TSGLAcceptType)
    ccall((:TSGLSetAcceptType,petsc),PetscErrorCode,(TS,TSGLAcceptType),arg1,arg2)
end

function TSEIMEXSetMaxRows(ts::TS,arg1::PetscInt)
    ccall((:TSEIMEXSetMaxRows,petsc),PetscErrorCode,(TS,PetscInt),ts,arg1)
end

function TSEIMEXSetRowCol(ts::TS,arg1::PetscInt,arg2::PetscInt)
    ccall((:TSEIMEXSetRowCol,petsc),PetscErrorCode,(TS,PetscInt,PetscInt),ts,arg1,arg2)
end

function TSEIMEXSetOrdAdapt(arg1::TS,arg2::PetscBool)
    ccall((:TSEIMEXSetOrdAdapt,petsc),PetscErrorCode,(TS,PetscBool),arg1,arg2)
end

function TSRKGetType(ts::TS,arg1::Ptr{TSRKType})
    ccall((:TSRKGetType,petsc),PetscErrorCode,(TS,Ptr{TSRKType}),ts,arg1)
end

function TSRKSetType(ts::TS,arg1::TSRKType)
    ccall((:TSRKSetType,petsc),PetscErrorCode,(TS,TSRKType),ts,arg1)
end

function TSRKSetFullyImplicit(arg1::TS,arg2::PetscBool)
    ccall((:TSRKSetFullyImplicit,petsc),PetscErrorCode,(TS,PetscBool),arg1,arg2)
end

function TSRKRegister(arg1::TSRKType,arg2::PetscInt,arg3::PetscInt,PetscReal::Ptr{Cint},arg4::Ptr{Cint},arg5::Ptr{Cint},arg6::Ptr{Cint},arg7::PetscInt,arg8::Ptr{Cint})
    ccall((:TSRKRegister,petsc),PetscErrorCode,(TSRKType,PetscInt,PetscInt,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},PetscInt,Ptr{Cint}),arg1,arg2,arg3,PetscReal,arg4,arg5,arg6,arg7,arg8)
end

function TSRKFinalizePackage()
    ccall((:TSRKFinalizePackage,petsc),PetscErrorCode,())
end

function TSRKInitializePackage()
    ccall((:TSRKInitializePackage,petsc),PetscErrorCode,())
end

function TSRKRegisterDestroy()
    ccall((:TSRKRegisterDestroy,petsc),PetscErrorCode,())
end

function TSARKIMEXGetType(ts::TS,arg1::Ptr{TSARKIMEXType})
    ccall((:TSARKIMEXGetType,petsc),PetscErrorCode,(TS,Ptr{TSARKIMEXType}),ts,arg1)
end

function TSARKIMEXSetType(ts::TS,arg1::TSARKIMEXType)
    ccall((:TSARKIMEXSetType,petsc),PetscErrorCode,(TS,TSARKIMEXType),ts,arg1)
end

function TSARKIMEXSetFullyImplicit(arg1::TS,arg2::PetscBool)
    ccall((:TSARKIMEXSetFullyImplicit,petsc),PetscErrorCode,(TS,PetscBool),arg1,arg2)
end

function TSARKIMEXRegister(arg1::TSARKIMEXType,arg2::PetscInt,arg3::PetscInt,PetscReal::Ptr{Cint},arg4::Ptr{Cint},arg5::Ptr{Cint},arg6::Ptr{Cint},arg7::Ptr{Cint},arg8::Ptr{Cint},arg9::Ptr{Cint},arg10::Ptr{Cint},arg11::PetscInt,arg12::Ptr{Cint},arg13::Ptr{Cint})
    ccall((:TSARKIMEXRegister,petsc),PetscErrorCode,(TSARKIMEXType,PetscInt,PetscInt,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},PetscInt,Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,PetscReal,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13)
end

function TSARKIMEXFinalizePackage()
    ccall((:TSARKIMEXFinalizePackage,petsc),PetscErrorCode,())
end

function TSARKIMEXInitializePackage()
    ccall((:TSARKIMEXInitializePackage,petsc),PetscErrorCode,())
end

function TSARKIMEXRegisterDestroy()
    ccall((:TSARKIMEXRegisterDestroy,petsc),PetscErrorCode,())
end

function TSRosWGetType(ts::TS,arg1::Ptr{TSRosWType})
    ccall((:TSRosWGetType,petsc),PetscErrorCode,(TS,Ptr{TSRosWType}),ts,arg1)
end

function TSRosWSetType(ts::TS,arg1::TSRosWType)
    ccall((:TSRosWSetType,petsc),PetscErrorCode,(TS,TSRosWType),ts,arg1)
end

function TSRosWSetRecomputeJacobian(arg1::TS,arg2::PetscBool)
    ccall((:TSRosWSetRecomputeJacobian,petsc),PetscErrorCode,(TS,PetscBool),arg1,arg2)
end

function TSRosWRegister(arg1::TSRosWType,arg2::PetscInt,arg3::PetscInt,PetscReal::Ptr{Cint},arg4::Ptr{Cint},arg5::Ptr{Cint},arg6::Ptr{Cint},arg7::PetscInt,arg8::Ptr{Cint})
    ccall((:TSRosWRegister,petsc),PetscErrorCode,(TSRosWType,PetscInt,PetscInt,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},PetscInt,Ptr{Cint}),arg1,arg2,arg3,PetscReal,arg4,arg5,arg6,arg7,arg8)
end

function TSRosWRegisterRos4(arg1::TSRosWType,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::Cint,arg5::Cint)
    ccall((:TSRosWRegisterRos4,petsc),PetscErrorCode,(TSRosWType,Cint,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4,arg5)
end

function TSRosWFinalizePackage()
    ccall((:TSRosWFinalizePackage,petsc),PetscErrorCode,())
end

function TSRosWInitializePackage()
    ccall((:TSRosWInitializePackage,petsc),PetscErrorCode,())
end

function TSRosWRegisterDestroy()
    ccall((:TSRosWRegisterDestroy,petsc),PetscErrorCode,())
end

function TSThetaSetTheta(arg1::TS,PetscReal::Cint)
    ccall((:TSThetaSetTheta,petsc),PetscErrorCode,(TS,Cint),arg1,PetscReal)
end

function TSThetaGetTheta(arg1::TS,arg2::Ptr{Cint})
    ccall((:TSThetaGetTheta,petsc),PetscErrorCode,(TS,Ptr{Cint}),arg1,arg2)
end

function TSThetaGetEndpoint(arg1::TS,arg2::Ptr{PetscBool})
    ccall((:TSThetaGetEndpoint,petsc),PetscErrorCode,(TS,Ptr{PetscBool}),arg1,arg2)
end

function TSThetaSetEndpoint(arg1::TS,arg2::PetscBool)
    ccall((:TSThetaSetEndpoint,petsc),PetscErrorCode,(TS,PetscBool),arg1,arg2)
end

function TSAlphaSetAdapt(arg1::TS,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:TSAlphaSetAdapt,petsc),PetscErrorCode,(TS,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function TSAlphaAdaptDefault(arg1::TS,PetscReal::Cint,arg2::Vec,arg3::Vec,arg4::Ptr{Cint},arg5::Ptr{PetscBool},arg6::Ptr{Void})
    ccall((:TSAlphaAdaptDefault,petsc),PetscErrorCode,(TS,Cint,Vec,Vec,Ptr{Cint},Ptr{PetscBool},Ptr{Void}),arg1,PetscReal,arg2,arg3,arg4,arg5,arg6)
end

function TSAlphaSetRadius(arg1::TS,PetscReal::Cint)
    ccall((:TSAlphaSetRadius,petsc),PetscErrorCode,(TS,Cint),arg1,PetscReal)
end

function TSAlphaSetParams(arg1::TS,PetscReal::Cint,arg2::Cint,arg3::Cint)
    ccall((:TSAlphaSetParams,petsc),PetscErrorCode,(TS,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3)
end

function TSAlphaGetParams(arg1::TS,arg2::Ptr{Cint},arg3::Ptr{Cint},arg4::Ptr{Cint})
    ccall((:TSAlphaGetParams,petsc),PetscErrorCode,(TS,Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function TSSetDM(arg1::TS,arg2::DM)
    ccall((:TSSetDM,petsc),PetscErrorCode,(TS,DM),arg1,arg2)
end

function TSGetDM(arg1::TS,arg2::Ptr{DM})
    ccall((:TSGetDM,petsc),PetscErrorCode,(TS,Ptr{DM}),arg1,arg2)
end

function SNESTSFormFunction(arg1::SNES,arg2::Vec,arg3::Vec,arg4::Ptr{Void})
    ccall((:SNESTSFormFunction,petsc),PetscErrorCode,(SNES,Vec,Vec,Ptr{Void}),arg1,arg2,arg3,arg4)
end

function SNESTSFormJacobian(arg1::SNES,arg2::Vec,arg3::Mat,arg4::Mat,arg5::Ptr{Void})
    ccall((:SNESTSFormJacobian,petsc),PetscErrorCode,(SNES,Vec,Mat,Mat,Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end

function VecFischer(arg1::Vec,arg2::Vec,arg3::Vec,arg4::Vec,arg5::Vec)
    ccall((:VecFischer,petsc),PetscErrorCode,(Vec,Vec,Vec,Vec,Vec),arg1,arg2,arg3,arg4,arg5)
end

function VecSFischer(arg1::Vec,arg2::Vec,arg3::Vec,arg4::Vec,PetscReal::Cint,arg5::Vec)
    ccall((:VecSFischer,petsc),PetscErrorCode,(Vec,Vec,Vec,Vec,Cint,Vec),arg1,arg2,arg3,arg4,PetscReal,arg5)
end

function MatDFischer(arg1::Mat,arg2::Vec,arg3::Vec,arg4::Vec,arg5::Vec,arg6::Vec,arg7::Vec,arg8::Vec,arg9::Vec)
    ccall((:MatDFischer,petsc),PetscErrorCode,(Mat,Vec,Vec,Vec,Vec,Vec,Vec,Vec,Vec),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function MatDSFischer(arg1::Mat,arg2::Vec,arg3::Vec,arg4::Vec,arg5::Vec,PetscReal::Cint,arg6::Vec,arg7::Vec,arg8::Vec,arg9::Vec,arg10::Vec)
    ccall((:MatDSFischer,petsc),PetscErrorCode,(Mat,Vec,Vec,Vec,Vec,Cint,Vec,Vec,Vec,Vec,Vec),arg1,arg2,arg3,arg4,arg5,PetscReal,arg6,arg7,arg8,arg9,arg10)
end

function TaoInitializePackage()
    ccall((:TaoInitializePackage,petsc),PetscErrorCode,())
end

function TaoFinalizePackage()
    ccall((:TaoFinalizePackage,petsc),PetscErrorCode,())
end

function TaoCreate(arg1::MPI_Comm,arg2::Ptr{Tao})
    ccall((:TaoCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{Tao}),arg1,arg2)
end

function TaoSetFromOptions(arg1::Tao)
    ccall((:TaoSetFromOptions,petsc),PetscErrorCode,(Tao,),arg1)
end

function TaoSetUp(arg1::Tao)
    ccall((:TaoSetUp,petsc),PetscErrorCode,(Tao,),arg1)
end

function TaoSetType(arg1::Tao,arg2::Ptr{Uint8})
    ccall((:TaoSetType,petsc),PetscErrorCode,(Tao,Ptr{Uint8}),arg1,arg2)
end

function TaoGetType(arg1::Tao,arg2::Ptr{Ptr{Uint8}})
    ccall((:TaoGetType,petsc),PetscErrorCode,(Tao,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function TaoSetApplicationContext(arg1::Tao,arg2::Ptr{Void})
    ccall((:TaoSetApplicationContext,petsc),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end

function TaoGetApplicationContext(arg1::Tao,arg2::Ptr{Void})
    ccall((:TaoGetApplicationContext,petsc),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end

function TaoDestroy(arg1::Ptr{Tao})
    ccall((:TaoDestroy,petsc),PetscErrorCode,(Ptr{Tao},),arg1)
end

function TaoSetOptionsPrefix(arg1::Tao,arg2::Ptr{Uint8})
    ccall((:TaoSetOptionsPrefix,petsc),PetscErrorCode,(Tao,Ptr{Uint8}),arg1,arg2)
end

function TaoView(arg1::Tao,arg2::PetscViewer)
    ccall((:TaoView,petsc),PetscErrorCode,(Tao,PetscViewer),arg1,arg2)
end

function TaoRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:TaoRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function TaoRegisterDestroy()
    ccall((:TaoRegisterDestroy,petsc),PetscErrorCode,())
end

function TaoGetConvergedReason(arg1::Tao,arg2::Ptr{TaoConvergedReason})
    ccall((:TaoGetConvergedReason,petsc),PetscErrorCode,(Tao,Ptr{TaoConvergedReason}),arg1,arg2)
end

function TaoGetSolutionStatus(arg1::Tao,arg2::Ptr{PetscInt},arg3::Ptr{Cint},arg4::Ptr{Cint},arg5::Ptr{Cint},arg6::Ptr{Cint},arg7::Ptr{TaoConvergedReason})
    ccall((:TaoGetSolutionStatus,petsc),PetscErrorCode,(Tao,Ptr{PetscInt},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{TaoConvergedReason}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function TaoSetConvergedReason(arg1::Tao,arg2::TaoConvergedReason)
    ccall((:TaoSetConvergedReason,petsc),PetscErrorCode,(Tao,TaoConvergedReason),arg1,arg2)
end

function TaoSetInitialVector(arg1::Tao,arg2::Vec)
    ccall((:TaoSetInitialVector,petsc),PetscErrorCode,(Tao,Vec),arg1,arg2)
end

function TaoGetSolutionVector(arg1::Tao,arg2::Ptr{Vec})
    ccall((:TaoGetSolutionVector,petsc),PetscErrorCode,(Tao,Ptr{Vec}),arg1,arg2)
end

function TaoGetGradientVector(arg1::Tao,arg2::Ptr{Vec})
    ccall((:TaoGetGradientVector,petsc),PetscErrorCode,(Tao,Ptr{Vec}),arg1,arg2)
end

function TaoSetObjectiveRoutine(arg1::Tao,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:TaoSetObjectiveRoutine,petsc),PetscErrorCode,(Tao,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function TaoSetGradientRoutine(arg1::Tao,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:TaoSetGradientRoutine,petsc),PetscErrorCode,(Tao,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function TaoSetObjectiveAndGradientRoutine(arg1::Tao,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:TaoSetObjectiveAndGradientRoutine,petsc),PetscErrorCode,(Tao,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function TaoSetHessianRoutine(arg1::Tao,arg2::Mat,arg3::Mat,arg4::Ptr{Void},arg5::Ptr{Void})
    ccall((:TaoSetHessianRoutine,petsc),PetscErrorCode,(Tao,Mat,Mat,Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end

function TaoSetSeparableObjectiveRoutine(arg1::Tao,arg2::Vec,arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:TaoSetSeparableObjectiveRoutine,petsc),PetscErrorCode,(Tao,Vec,Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function TaoSetConstraintsRoutine(arg1::Tao,arg2::Vec,arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:TaoSetConstraintsRoutine,petsc),PetscErrorCode,(Tao,Vec,Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function TaoSetInequalityConstraintsRoutine(arg1::Tao,arg2::Vec,arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:TaoSetInequalityConstraintsRoutine,petsc),PetscErrorCode,(Tao,Vec,Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function TaoSetEqualityConstraintsRoutine(arg1::Tao,arg2::Vec,arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:TaoSetEqualityConstraintsRoutine,petsc),PetscErrorCode,(Tao,Vec,Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function TaoSetJacobianRoutine(arg1::Tao,arg2::Mat,arg3::Mat,arg4::Ptr{Void},arg5::Ptr{Void})
    ccall((:TaoSetJacobianRoutine,petsc),PetscErrorCode,(Tao,Mat,Mat,Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end

function TaoSetJacobianStateRoutine(arg1::Tao,arg2::Mat,arg3::Mat,arg4::Mat,arg5::Ptr{Void},arg6::Ptr{Void})
    ccall((:TaoSetJacobianStateRoutine,petsc),PetscErrorCode,(Tao,Mat,Mat,Mat,Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function TaoSetJacobianDesignRoutine(arg1::Tao,arg2::Mat,arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:TaoSetJacobianDesignRoutine,petsc),PetscErrorCode,(Tao,Mat,Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function TaoSetJacobianInequalityRoutine(arg1::Tao,arg2::Mat,arg3::Mat,arg4::Ptr{Void},arg5::Ptr{Void})
    ccall((:TaoSetJacobianInequalityRoutine,petsc),PetscErrorCode,(Tao,Mat,Mat,Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end

function TaoSetJacobianEqualityRoutine(arg1::Tao,arg2::Mat,arg3::Mat,arg4::Ptr{Void},arg5::Ptr{Void})
    ccall((:TaoSetJacobianEqualityRoutine,petsc),PetscErrorCode,(Tao,Mat,Mat,Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end

function TaoSetStateDesignIS(arg1::Tao,arg2::IS,arg3::IS)
    ccall((:TaoSetStateDesignIS,petsc),PetscErrorCode,(Tao,IS,IS),arg1,arg2,arg3)
end

function TaoComputeObjective(arg1::Tao,arg2::Vec,arg3::Ptr{Cint})
    ccall((:TaoComputeObjective,petsc),PetscErrorCode,(Tao,Vec,Ptr{Cint}),arg1,arg2,arg3)
end

function TaoComputeSeparableObjective(arg1::Tao,arg2::Vec,arg3::Vec)
    ccall((:TaoComputeSeparableObjective,petsc),PetscErrorCode,(Tao,Vec,Vec),arg1,arg2,arg3)
end

function TaoComputeGradient(arg1::Tao,arg2::Vec,arg3::Vec)
    ccall((:TaoComputeGradient,petsc),PetscErrorCode,(Tao,Vec,Vec),arg1,arg2,arg3)
end

function TaoComputeObjectiveAndGradient(arg1::Tao,arg2::Vec,arg3::Ptr{Cint},arg4::Vec)
    ccall((:TaoComputeObjectiveAndGradient,petsc),PetscErrorCode,(Tao,Vec,Ptr{Cint},Vec),arg1,arg2,arg3,arg4)
end

function TaoComputeConstraints(arg1::Tao,arg2::Vec,arg3::Vec)
    ccall((:TaoComputeConstraints,petsc),PetscErrorCode,(Tao,Vec,Vec),arg1,arg2,arg3)
end

function TaoComputeInequalityConstraints(arg1::Tao,arg2::Vec,arg3::Vec)
    ccall((:TaoComputeInequalityConstraints,petsc),PetscErrorCode,(Tao,Vec,Vec),arg1,arg2,arg3)
end

function TaoComputeEqualityConstraints(arg1::Tao,arg2::Vec,arg3::Vec)
    ccall((:TaoComputeEqualityConstraints,petsc),PetscErrorCode,(Tao,Vec,Vec),arg1,arg2,arg3)
end

function TaoDefaultComputeGradient(arg1::Tao,arg2::Vec,arg3::Vec,arg4::Ptr{Void})
    ccall((:TaoDefaultComputeGradient,petsc),PetscErrorCode,(Tao,Vec,Vec,Ptr{Void}),arg1,arg2,arg3,arg4)
end

function TaoIsObjectiveDefined(arg1::Tao,arg2::Ptr{PetscBool})
    ccall((:TaoIsObjectiveDefined,petsc),PetscErrorCode,(Tao,Ptr{PetscBool}),arg1,arg2)
end

function TaoIsGradientDefined(arg1::Tao,arg2::Ptr{PetscBool})
    ccall((:TaoIsGradientDefined,petsc),PetscErrorCode,(Tao,Ptr{PetscBool}),arg1,arg2)
end

function TaoIsObjectiveAndGradientDefined(arg1::Tao,arg2::Ptr{PetscBool})
    ccall((:TaoIsObjectiveAndGradientDefined,petsc),PetscErrorCode,(Tao,Ptr{PetscBool}),arg1,arg2)
end

function TaoComputeHessian(arg1::Tao,arg2::Vec,arg3::Mat,arg4::Mat)
    ccall((:TaoComputeHessian,petsc),PetscErrorCode,(Tao,Vec,Mat,Mat),arg1,arg2,arg3,arg4)
end

function TaoComputeJacobian(arg1::Tao,arg2::Vec,arg3::Mat,arg4::Mat)
    ccall((:TaoComputeJacobian,petsc),PetscErrorCode,(Tao,Vec,Mat,Mat),arg1,arg2,arg3,arg4)
end

function TaoComputeJacobianState(arg1::Tao,arg2::Vec,arg3::Mat,arg4::Mat,arg5::Mat)
    ccall((:TaoComputeJacobianState,petsc),PetscErrorCode,(Tao,Vec,Mat,Mat,Mat),arg1,arg2,arg3,arg4,arg5)
end

function TaoComputeJacobianEquality(arg1::Tao,arg2::Vec,arg3::Mat,arg4::Mat)
    ccall((:TaoComputeJacobianEquality,petsc),PetscErrorCode,(Tao,Vec,Mat,Mat),arg1,arg2,arg3,arg4)
end

function TaoComputeJacobianInequality(arg1::Tao,arg2::Vec,arg3::Mat,arg4::Mat)
    ccall((:TaoComputeJacobianInequality,petsc),PetscErrorCode,(Tao,Vec,Mat,Mat),arg1,arg2,arg3,arg4)
end

function TaoComputeJacobianDesign(arg1::Tao,arg2::Vec,arg3::Mat)
    ccall((:TaoComputeJacobianDesign,petsc),PetscErrorCode,(Tao,Vec,Mat),arg1,arg2,arg3)
end

function TaoDefaultComputeHessian(arg1::Tao,arg2::Vec,arg3::Mat,arg4::Mat,arg5::Ptr{Void})
    ccall((:TaoDefaultComputeHessian,petsc),PetscErrorCode,(Tao,Vec,Mat,Mat,Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end

function TaoDefaultComputeHessianColor(arg1::Tao,arg2::Vec,arg3::Mat,arg4::Mat,arg5::Ptr{Void})
    ccall((:TaoDefaultComputeHessianColor,petsc),PetscErrorCode,(Tao,Vec,Mat,Mat,Ptr{Void}),arg1,arg2,arg3,arg4,arg5)
end

function TaoComputeDualVariables(arg1::Tao,arg2::Vec,arg3::Vec)
    ccall((:TaoComputeDualVariables,petsc),PetscErrorCode,(Tao,Vec,Vec),arg1,arg2,arg3)
end

function TaoComputeDualVariables(arg1::Tao,arg2::Vec,arg3::Vec)
    ccall((:TaoComputeDualVariables,petsc),PetscErrorCode,(Tao,Vec,Vec),arg1,arg2,arg3)
end

function TaoSetVariableBounds(arg1::Tao,arg2::Vec,arg3::Vec)
    ccall((:TaoSetVariableBounds,petsc),PetscErrorCode,(Tao,Vec,Vec),arg1,arg2,arg3)
end

function TaoGetVariableBounds(arg1::Tao,arg2::Ptr{Vec},arg3::Ptr{Vec})
    ccall((:TaoGetVariableBounds,petsc),PetscErrorCode,(Tao,Ptr{Vec},Ptr{Vec}),arg1,arg2,arg3)
end

function TaoGetDualVariables(arg1::Tao,arg2::Ptr{Vec},arg3::Ptr{Vec})
    ccall((:TaoGetDualVariables,petsc),PetscErrorCode,(Tao,Ptr{Vec},Ptr{Vec}),arg1,arg2,arg3)
end

function TaoSetInequalityBounds(arg1::Tao,arg2::Vec,arg3::Vec)
    ccall((:TaoSetInequalityBounds,petsc),PetscErrorCode,(Tao,Vec,Vec),arg1,arg2,arg3)
end

function TaoGetInequalityBounds(arg1::Tao,arg2::Ptr{Vec},arg3::Ptr{Vec})
    ccall((:TaoGetInequalityBounds,petsc),PetscErrorCode,(Tao,Ptr{Vec},Ptr{Vec}),arg1,arg2,arg3)
end

function TaoSetVariableBoundsRoutine(arg1::Tao,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:TaoSetVariableBoundsRoutine,petsc),PetscErrorCode,(Tao,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function TaoComputeVariableBounds(arg1::Tao)
    ccall((:TaoComputeVariableBounds,petsc),PetscErrorCode,(Tao,),arg1)
end

function TaoGetTolerances(arg1::Tao,arg2::Ptr{Cint},arg3::Ptr{Cint},arg4::Ptr{Cint},arg5::Ptr{Cint},arg6::Ptr{Cint})
    ccall((:TaoGetTolerances,petsc),PetscErrorCode,(Tao,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function TaoSetTolerances(arg1::Tao,PetscReal::Cint,arg2::Cint,arg3::Cint,arg4::Cint,arg5::Cint)
    ccall((:TaoSetTolerances,petsc),PetscErrorCode,(Tao,Cint,Cint,Cint,Cint,Cint),arg1,PetscReal,arg2,arg3,arg4,arg5)
end

function TaoGetConstraintTolerances(arg1::Tao,arg2::Ptr{Cint},arg3::Ptr{Cint})
    ccall((:TaoGetConstraintTolerances,petsc),PetscErrorCode,(Tao,Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3)
end

function TaoSetConstraintTolerances(arg1::Tao,PetscReal::Cint,arg2::Cint)
    ccall((:TaoSetConstraintTolerances,petsc),PetscErrorCode,(Tao,Cint,Cint),arg1,PetscReal,arg2)
end

function TaoSetFunctionLowerBound(arg1::Tao,PetscReal::Cint)
    ccall((:TaoSetFunctionLowerBound,petsc),PetscErrorCode,(Tao,Cint),arg1,PetscReal)
end

function TaoSetInitialTrustRegionRadius(arg1::Tao,PetscReal::Cint)
    ccall((:TaoSetInitialTrustRegionRadius,petsc),PetscErrorCode,(Tao,Cint),arg1,PetscReal)
end

function TaoSetMaximumIterations(arg1::Tao,arg2::PetscInt)
    ccall((:TaoSetMaximumIterations,petsc),PetscErrorCode,(Tao,PetscInt),arg1,arg2)
end

function TaoSetMaximumFunctionEvaluations(arg1::Tao,arg2::PetscInt)
    ccall((:TaoSetMaximumFunctionEvaluations,petsc),PetscErrorCode,(Tao,PetscInt),arg1,arg2)
end

function TaoGetFunctionLowerBound(arg1::Tao,arg2::Ptr{Cint})
    ccall((:TaoGetFunctionLowerBound,petsc),PetscErrorCode,(Tao,Ptr{Cint}),arg1,arg2)
end

function TaoGetInitialTrustRegionRadius(arg1::Tao,arg2::Ptr{Cint})
    ccall((:TaoGetInitialTrustRegionRadius,petsc),PetscErrorCode,(Tao,Ptr{Cint}),arg1,arg2)
end

function TaoGetCurrentTrustRegionRadius(arg1::Tao,arg2::Ptr{Cint})
    ccall((:TaoGetCurrentTrustRegionRadius,petsc),PetscErrorCode,(Tao,Ptr{Cint}),arg1,arg2)
end

function TaoGetMaximumIterations(arg1::Tao,arg2::Ptr{PetscInt})
    ccall((:TaoGetMaximumIterations,petsc),PetscErrorCode,(Tao,Ptr{PetscInt}),arg1,arg2)
end

function TaoGetCurrentFunctionEvaluations(arg1::Tao,arg2::Ptr{PetscInt})
    ccall((:TaoGetCurrentFunctionEvaluations,petsc),PetscErrorCode,(Tao,Ptr{PetscInt}),arg1,arg2)
end

function TaoGetMaximumFunctionEvaluations(arg1::Tao,arg2::Ptr{PetscInt})
    ccall((:TaoGetMaximumFunctionEvaluations,petsc),PetscErrorCode,(Tao,Ptr{PetscInt}),arg1,arg2)
end

function TaoGetIterationNumber(arg1::Tao,arg2::Ptr{PetscInt})
    ccall((:TaoGetIterationNumber,petsc),PetscErrorCode,(Tao,Ptr{PetscInt}),arg1,arg2)
end

function TaoSetIterationNumber(arg1::Tao,arg2::PetscInt)
    ccall((:TaoSetIterationNumber,petsc),PetscErrorCode,(Tao,PetscInt),arg1,arg2)
end

function TaoGetTotalIterationNumber(arg1::Tao,arg2::Ptr{PetscInt})
    ccall((:TaoGetTotalIterationNumber,petsc),PetscErrorCode,(Tao,Ptr{PetscInt}),arg1,arg2)
end

function TaoSetTotalIterationNumber(arg1::Tao,arg2::PetscInt)
    ccall((:TaoSetTotalIterationNumber,petsc),PetscErrorCode,(Tao,PetscInt),arg1,arg2)
end

function TaoSetOptionsPrefix(arg1::Tao,p::Ptr{Uint8})
    ccall((:TaoSetOptionsPrefix,petsc),PetscErrorCode,(Tao,Ptr{Uint8}),arg1,p)
end

function TaoAppendOptionsPrefix(arg1::Tao,p::Ptr{Uint8})
    ccall((:TaoAppendOptionsPrefix,petsc),PetscErrorCode,(Tao,Ptr{Uint8}),arg1,p)
end

function TaoGetOptionsPrefix(arg1::Tao,p::Ptr{Ptr{Uint8}})
    ccall((:TaoGetOptionsPrefix,petsc),PetscErrorCode,(Tao,Ptr{Ptr{Uint8}}),arg1,p)
end

function TaoResetStatistics(arg1::Tao)
    ccall((:TaoResetStatistics,petsc),PetscErrorCode,(Tao,),arg1)
end

function TaoGetKSP(arg1::Tao,arg2::Ptr{KSP})
    ccall((:TaoGetKSP,petsc),PetscErrorCode,(Tao,Ptr{KSP}),arg1,arg2)
end

function TaoGetLinearSolveIterations(arg1::Tao,arg2::Ptr{PetscInt})
    ccall((:TaoGetLinearSolveIterations,petsc),PetscErrorCode,(Tao,Ptr{PetscInt}),arg1,arg2)
end

function TaoLineSearchCreate(arg1::MPI_Comm,arg2::Ptr{TaoLineSearch})
    ccall((:TaoLineSearchCreate,petsc),PetscErrorCode,(MPI_Comm,Ptr{TaoLineSearch}),arg1,arg2)
end

function TaoLineSearchSetFromOptions(arg1::TaoLineSearch)
    ccall((:TaoLineSearchSetFromOptions,petsc),PetscErrorCode,(TaoLineSearch,),arg1)
end

function TaoLineSearchSetUp(arg1::TaoLineSearch)
    ccall((:TaoLineSearchSetUp,petsc),PetscErrorCode,(TaoLineSearch,),arg1)
end

function TaoLineSearchDestroy(arg1::Ptr{TaoLineSearch})
    ccall((:TaoLineSearchDestroy,petsc),PetscErrorCode,(Ptr{TaoLineSearch},),arg1)
end

function TaoLineSearchView(arg1::TaoLineSearch,arg2::PetscViewer)
    ccall((:TaoLineSearchView,petsc),PetscErrorCode,(TaoLineSearch,PetscViewer),arg1,arg2)
end

function TaoLineSearchReset(arg1::TaoLineSearch)
    ccall((:TaoLineSearchReset,petsc),PetscErrorCode,(TaoLineSearch,),arg1)
end

function TaoLineSearchAppendOptionsPrefix(arg1::TaoLineSearch,prefix::Ptr{Uint8})
    ccall((:TaoLineSearchAppendOptionsPrefix,petsc),PetscErrorCode,(TaoLineSearch,Ptr{Uint8}),arg1,prefix)
end

function TaoLineSearchGetOptionsPrefix(arg1::TaoLineSearch,prefix::Ptr{Ptr{Uint8}})
    ccall((:TaoLineSearchGetOptionsPrefix,petsc),PetscErrorCode,(TaoLineSearch,Ptr{Ptr{Uint8}}),arg1,prefix)
end

function TaoLineSearchApply(arg1::TaoLineSearch,arg2::Vec,arg3::Ptr{Cint},arg4::Vec,arg5::Vec,arg6::Ptr{Cint},arg7::Ptr{TaoLineSearchConvergedReason})
    ccall((:TaoLineSearchApply,petsc),PetscErrorCode,(TaoLineSearch,Vec,Ptr{Cint},Vec,Vec,Ptr{Cint},Ptr{TaoLineSearchConvergedReason}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function TaoLineSearchGetStepLength(arg1::TaoLineSearch,arg2::Ptr{Cint})
    ccall((:TaoLineSearchGetStepLength,petsc),PetscErrorCode,(TaoLineSearch,Ptr{Cint}),arg1,arg2)
end

function TaoLineSearchGetStartingVector(arg1::TaoLineSearch,arg2::Ptr{Vec})
    ccall((:TaoLineSearchGetStartingVector,petsc),PetscErrorCode,(TaoLineSearch,Ptr{Vec}),arg1,arg2)
end

function TaoLineSearchGetStepDirection(arg1::TaoLineSearch,arg2::Ptr{Vec})
    ccall((:TaoLineSearchGetStepDirection,petsc),PetscErrorCode,(TaoLineSearch,Ptr{Vec}),arg1,arg2)
end

function TaoLineSearchSetInitialStepLength(arg1::TaoLineSearch,PetscReal::Cint)
    ccall((:TaoLineSearchSetInitialStepLength,petsc),PetscErrorCode,(TaoLineSearch,Cint),arg1,PetscReal)
end

function TaoLineSearchGetSolution(arg1::TaoLineSearch,arg2::Vec,arg3::Ptr{Cint},arg4::Vec,arg5::Ptr{Cint},arg6::Ptr{TaoLineSearchConvergedReason})
    ccall((:TaoLineSearchGetSolution,petsc),PetscErrorCode,(TaoLineSearch,Vec,Ptr{Cint},Vec,Ptr{Cint},Ptr{TaoLineSearchConvergedReason}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function TaoLineSearchGetFullStepObjective(arg1::TaoLineSearch,arg2::Ptr{Cint})
    ccall((:TaoLineSearchGetFullStepObjective,petsc),PetscErrorCode,(TaoLineSearch,Ptr{Cint}),arg1,arg2)
end

function TaoLineSearchGetNumberFunctionEvaluations(arg1::TaoLineSearch,arg2::Ptr{PetscInt},arg3::Ptr{PetscInt},arg4::Ptr{PetscInt})
    ccall((:TaoLineSearchGetNumberFunctionEvaluations,petsc),PetscErrorCode,(TaoLineSearch,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function TaoLineSearchGetType(arg1::TaoLineSearch,arg2::Ptr{Ptr{Uint8}})
    ccall((:TaoLineSearchGetType,petsc),PetscErrorCode,(TaoLineSearch,Ptr{Ptr{Uint8}}),arg1,arg2)
end

function TaoLineSearchSetType(arg1::TaoLineSearch,arg2::Ptr{Uint8})
    ccall((:TaoLineSearchSetType,petsc),PetscErrorCode,(TaoLineSearch,Ptr{Uint8}),arg1,arg2)
end

function TaoLineSearchIsUsingTaoRoutines(arg1::TaoLineSearch,arg2::Ptr{PetscBool})
    ccall((:TaoLineSearchIsUsingTaoRoutines,petsc),PetscErrorCode,(TaoLineSearch,Ptr{PetscBool}),arg1,arg2)
end

function TaoLineSearchSetObjectiveAndGTSRoutine(arg1::TaoLineSearch,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:TaoLineSearchSetObjectiveAndGTSRoutine,petsc),PetscErrorCode,(TaoLineSearch,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function TaoLineSearchSetObjectiveRoutine(arg1::TaoLineSearch,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:TaoLineSearchSetObjectiveRoutine,petsc),PetscErrorCode,(TaoLineSearch,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function TaoLineSearchSetGradientRoutine(arg1::TaoLineSearch,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:TaoLineSearchSetGradientRoutine,petsc),PetscErrorCode,(TaoLineSearch,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function TaoLineSearchSetObjectiveAndGradientRoutine(arg1::TaoLineSearch,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:TaoLineSearchSetObjectiveAndGradientRoutine,petsc),PetscErrorCode,(TaoLineSearch,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function TaoLineSearchComputeObjective(arg1::TaoLineSearch,arg2::Vec,arg3::Ptr{Cint})
    ccall((:TaoLineSearchComputeObjective,petsc),PetscErrorCode,(TaoLineSearch,Vec,Ptr{Cint}),arg1,arg2,arg3)
end

function TaoLineSearchComputeGradient(arg1::TaoLineSearch,arg2::Vec,arg3::Vec)
    ccall((:TaoLineSearchComputeGradient,petsc),PetscErrorCode,(TaoLineSearch,Vec,Vec),arg1,arg2,arg3)
end

function TaoLineSearchComputeObjectiveAndGradient(arg1::TaoLineSearch,arg2::Vec,arg3::Ptr{Cint},arg4::Vec)
    ccall((:TaoLineSearchComputeObjectiveAndGradient,petsc),PetscErrorCode,(TaoLineSearch,Vec,Ptr{Cint},Vec),arg1,arg2,arg3,arg4)
end

function TaoLineSearchComputeObjectiveAndGTS(arg1::TaoLineSearch,arg2::Vec,arg3::Ptr{Cint},arg4::Ptr{Cint})
    ccall((:TaoLineSearchComputeObjectiveAndGTS,petsc),PetscErrorCode,(TaoLineSearch,Vec,Ptr{Cint},Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function TaoLineSearchSetVariableBounds(arg1::TaoLineSearch,arg2::Vec,arg3::Vec)
    ccall((:TaoLineSearchSetVariableBounds,petsc),PetscErrorCode,(TaoLineSearch,Vec,Vec),arg1,arg2,arg3)
end

function TaoLineSearchInitializePackage()
    ccall((:TaoLineSearchInitializePackage,petsc),PetscErrorCode,())
end

function TaoLineSearchFinalizePackage()
    ccall((:TaoLineSearchFinalizePackage,petsc),PetscErrorCode,())
end

function TaoLineSearchRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:TaoLineSearchRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function TaoLineSearchUseTaoRoutines(arg1::TaoLineSearch,arg2::Tao)
    ccall((:TaoLineSearchUseTaoRoutines,petsc),PetscErrorCode,(TaoLineSearch,Tao),arg1,arg2)
end

function TaoGetLineSearch(arg1::Tao,arg2::Ptr{TaoLineSearch})
    ccall((:TaoGetLineSearch,petsc),PetscErrorCode,(Tao,Ptr{TaoLineSearch}),arg1,arg2)
end

function TaoSetConvergenceHistory(arg1::Tao,arg2::Ptr{Cint},arg3::Ptr{Cint},arg4::Ptr{Cint},arg5::Ptr{PetscInt},arg6::PetscInt,arg7::PetscBool)
    ccall((:TaoSetConvergenceHistory,petsc),PetscErrorCode,(Tao,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{PetscInt},PetscInt,PetscBool),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function TaoGetConvergenceHistory(arg1::Tao,arg2::Ptr{Ptr{Cint}},arg3::Ptr{Ptr{Cint}},arg4::Ptr{Ptr{Cint}},arg5::Ptr{Ptr{PetscInt}},arg6::Ptr{PetscInt})
    ccall((:TaoGetConvergenceHistory,petsc),PetscErrorCode,(Tao,Ptr{Ptr{Cint}},Ptr{Ptr{Cint}},Ptr{Ptr{Cint}},Ptr{Ptr{PetscInt}},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function TaoSetMonitor(arg1::Tao,arg2::Ptr{Void},arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:TaoSetMonitor,petsc),PetscErrorCode,(Tao,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function TaoCancelMonitors(arg1::Tao)
    ccall((:TaoCancelMonitors,petsc),PetscErrorCode,(Tao,),arg1)
end

function TaoDefaultMonitor(arg1::Tao,arg2::Ptr{Void})
    ccall((:TaoDefaultMonitor,petsc),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end

function TaoDefaultSMonitor(arg1::Tao,arg2::Ptr{Void})
    ccall((:TaoDefaultSMonitor,petsc),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end

function TaoDefaultCMonitor(arg1::Tao,arg2::Ptr{Void})
    ccall((:TaoDefaultCMonitor,petsc),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end

function TaoSolutionMonitor(arg1::Tao,arg2::Ptr{Void})
    ccall((:TaoSolutionMonitor,petsc),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end

function TaoSeparableObjectiveMonitor(arg1::Tao,arg2::Ptr{Void})
    ccall((:TaoSeparableObjectiveMonitor,petsc),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end

function TaoGradientMonitor(arg1::Tao,arg2::Ptr{Void})
    ccall((:TaoGradientMonitor,petsc),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end

function TaoStepDirectionMonitor(arg1::Tao,arg2::Ptr{Void})
    ccall((:TaoStepDirectionMonitor,petsc),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end

function TaoDrawSolutionMonitor(arg1::Tao,arg2::Ptr{Void})
    ccall((:TaoDrawSolutionMonitor,petsc),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end

function TaoDrawStepMonitor(arg1::Tao,arg2::Ptr{Void})
    ccall((:TaoDrawStepMonitor,petsc),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end

function TaoDrawGradientMonitor(arg1::Tao,arg2::Ptr{Void})
    ccall((:TaoDrawGradientMonitor,petsc),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end

function TaoAddLineSearchCounts(arg1::Tao)
    ccall((:TaoAddLineSearchCounts,petsc),PetscErrorCode,(Tao,),arg1)
end

function TaoDefaultConvergenceTest(arg1::Tao,arg2::Ptr{Void})
    ccall((:TaoDefaultConvergenceTest,petsc),PetscErrorCode,(Tao,Ptr{Void}),arg1,arg2)
end

function TaoSetConvergenceTest(arg1::Tao,arg2::Ptr{Void},arg3::Ptr{Void})
    ccall((:TaoSetConvergenceTest,petsc),PetscErrorCode,(Tao,Ptr{Void},Ptr{Void}),arg1,arg2,arg3)
end

function TaoSQPCONSetStateDesignIS(arg1::Tao,arg2::IS,arg3::IS)
    ccall((:TaoSQPCONSetStateDesignIS,petsc),PetscErrorCode,(Tao,IS,IS),arg1,arg2,arg3)
end

function TaoLCLSetStateDesignIS(arg1::Tao,arg2::IS,arg3::IS)
    ccall((:TaoLCLSetStateDesignIS,petsc),PetscErrorCode,(Tao,IS,IS),arg1,arg2,arg3)
end

function TaoMonitor(arg1::Tao,arg2::PetscInt,PetscReal::Cint,arg3::Cint,arg4::Cint,arg5::Cint,arg6::Ptr{TaoConvergedReason})
    ccall((:TaoMonitor,petsc),PetscErrorCode,(Tao,PetscInt,Cint,Cint,Cint,Cint,Ptr{TaoConvergedReason}),arg1,arg2,PetscReal,arg3,arg4,arg5,arg6)
end

end # end module
