# Julia wrapper for header: /home/jared/build/petsc-3.6.0/include/petscksp.h
# Automatically generated using Clang.jl wrap_c, version 0.0.0


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

function PetscMallocDump(arg1::Ptr{FILE})
    ccall((:PetscMallocDump,petsc),PetscErrorCode,(Ptr{FILE},),arg1)
end

function PetscMallocDumpLog(arg1::Ptr{FILE})
    ccall((:PetscMallocDumpLog,petsc),PetscErrorCode,(Ptr{FILE},),arg1)
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

function PetscBitMemcpy(arg1::Ptr{Void},arg2::Int32,arg3::Ptr{Void},arg4::Int32,arg5::Int32,arg6::PetscDataType)
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

function PetscStrendswithwhich(arg1::Ptr{Uint8},arg2::Ptr{Ptr{Uint8}},arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
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

function PetscStrNArrayallocpy(arg1::Int32,arg2::Ptr{Ptr{Uint8}},arg3::Ptr{Ptr{Ptr{Uint8}}})
    ccall((:PetscStrNArrayallocpy,petsc),PetscErrorCode,(PetscInt,Ptr{Ptr{Uint8}},Ptr{Ptr{Ptr{Uint8}}}),arg1,arg2,arg3)
end

function PetscStrNArrayDestroy(arg1::Int32,arg2::Ptr{Ptr{Ptr{Uint8}}})
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

function PetscEListFind(arg1::Int32,arg2::Ptr{Ptr{Uint8}},arg3::Ptr{Uint8},arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Ptr{PetscBool})
    ccall((:PetscEListFind,petsc),PetscErrorCode,(PetscInt,Ptr{Ptr{Uint8}},Ptr{Uint8},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function PetscEnumFind(arg1::Ptr{Ptr{Uint8}},arg2::Ptr{Uint8},arg3::Ptr{PetscEnum},arg4::Ptr{PetscBool})
    ccall((:PetscEnumFind,petsc),PetscErrorCode,(Ptr{Ptr{Uint8}},Ptr{Uint8},Ptr{PetscEnum},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function PetscMaxSum(arg1::MPI_Comm,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscMaxSum,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function MPIULong_Send(arg1::Ptr{Void},arg2::Int32,arg3::MPI_Datatype,arg4::PetscMPIInt,arg5::PetscMPIInt,arg6::MPI_Comm)
    ccall((:MPIULong_Send,petsc),PetscErrorCode,(Ptr{Void},PetscInt,MPI_Datatype,PetscMPIInt,PetscMPIInt,MPI_Comm),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MPIULong_Recv(arg1::Ptr{Void},arg2::Int32,arg3::MPI_Datatype,arg4::PetscMPIInt,arg5::PetscMPIInt,arg6::MPI_Comm)
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

function PetscCheckPointerSetIntensity(arg1::Int32)
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

function PetscStackPrint(arg1::Ptr{PetscStack},arg2::Ptr{FILE})
    ccall((:PetscStackPrint,petsc),PetscErrorCode,(Ptr{PetscStack},Ptr{FILE}),arg1,arg2)
end

function PetscStackView(arg1::Ptr{FILE})
    ccall((:PetscStackView,petsc),PetscErrorCode,(Ptr{FILE},),arg1)
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

function PetscObjectSetTabLevel(arg1::PetscObject,arg2::Int32)
    ccall((:PetscObjectSetTabLevel,petsc),PetscErrorCode,(PetscObject,PetscInt),arg1,arg2)
end

function PetscObjectGetTabLevel(arg1::PetscObject,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscObjectGetTabLevel,petsc),PetscErrorCode,(PetscObject,Ptr{PetscInt}),arg1,arg2)
end

function PetscObjectIncrementTabLevel(arg1::PetscObject,arg2::PetscObject,arg3::Int32)
    ccall((:PetscObjectIncrementTabLevel,petsc),PetscErrorCode,(PetscObject,PetscObject,PetscInt),arg1,arg2,arg3)
end

function PetscObjectReference(arg1::PetscObject)
    ccall((:PetscObjectReference,petsc),PetscErrorCode,(PetscObject,),arg1)
end

function PetscObjectGetReference(arg1::PetscObject,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
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

function PetscObjectsListGetGlobalNumbering(arg1::MPI_Comm,arg2::Int32,arg3::Ptr{PetscObject},arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscObjectsListGetGlobalNumbering,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{PetscObject},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end

function PetscOptionsHasName(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{PetscBool})
    ccall((:PetscOptionsHasName,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{PetscBool}),arg1,arg2,arg3)
end

function PetscOptionsGetInt(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{PetscBool})
    ccall((:PetscOptionsGetInt,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function PetscOptionsGetBool(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{PetscBool},arg4::Ptr{PetscBool})
    ccall((:PetscOptionsGetBool,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function PetscOptionsGetReal(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{Cint},arg4::Ptr{PetscBool})
    ccall((:PetscOptionsGetReal,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{Cint},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function PetscOptionsGetScalar(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{Float64},arg4::Ptr{PetscBool})
    ccall((:PetscOptionsGetScalar,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{Float64},Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function PetscOptionsGetIntArray(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Ptr{PetscBool})
    ccall((:PetscOptionsGetIntArray,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function PetscOptionsGetRealArray(arg1::Ptr{Uint8},arg2::Ptr{Uint8},Float64::Ptr{Cint},arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{PetscBool})
    ccall((:PetscOptionsGetRealArray,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{Cint},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,PetscReal,arg3,arg4)
end

function PetscOptionsGetScalarArray(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{Float64},arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Ptr{PetscBool})
    ccall((:PetscOptionsGetScalarArray,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{Float64},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function PetscOptionsGetBoolArray(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{PetscBool},arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Ptr{PetscBool})
    ccall((:PetscOptionsGetBoolArray,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{PetscBool},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function PetscOptionsGetString(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{Uint8},size_t::Cint,arg4::Ptr{PetscBool})
    ccall((:PetscOptionsGetString,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Cint,Ptr{PetscBool}),arg1,arg2,arg3,size_t,arg4)
end

function PetscOptionsGetStringArray(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{Ptr{Uint8}},arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Ptr{PetscBool})
    ccall((:PetscOptionsGetStringArray,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{Ptr{Uint8}},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function PetscOptionsGetEList(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{Ptr{Uint8}},arg4::Int32,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Ptr{PetscBool})
    ccall((:PetscOptionsGetEList,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{Ptr{Uint8}},PetscInt,Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscOptionsGetEnum(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{Ptr{Uint8}},arg4::Ptr{PetscEnum},arg5::Ptr{PetscBool})
    ccall((:PetscOptionsGetEnum,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Ptr{Ptr{Uint8}},Ptr{PetscEnum},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function PetscOptionsGetEnumArray(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Ptr{Ptr{Uint8}},arg4::Ptr{PetscEnum},arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Ptr{PetscBool})
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

function PetscOptionsAllUsed(arg1::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
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

function PetscOptionsStringToInt(arg1::Ptr{Uint8},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
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

function PetscOptionsInt_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Int32,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Ptr{PetscBool})
    ccall((:PetscOptionsInt_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},PetscInt,Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscOptionsReal_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},Float64::Cint,arg5::Ptr{Cint},arg6::Ptr{PetscBool})
    ccall((:PetscOptionsReal_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Cint,Ptr{Cint},Ptr{PetscBool}),arg1,arg2,arg3,arg4,PetscReal,arg5,arg6)
end

function PetscOptionsScalar_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Float64,arg6::Ptr{Float64},arg7::Ptr{PetscBool})
    ccall((:PetscOptionsScalar_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},PetscScalar,Ptr{Float64},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
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

function PetscOptionsEList_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Ptr{Ptr{Uint8}},arg6::Int32,arg7::Ptr{Uint8},arg8::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg9::Ptr{PetscBool})
    ccall((:PetscOptionsEList_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Ptr{Uint8}},PetscInt,Ptr{Uint8},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function PetscOptionsRealArray_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},Float64::Ptr{Cint},arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Ptr{PetscBool})
    ccall((:PetscOptionsRealArray_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Cint},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,PetscReal,arg5,arg6)
end

function PetscOptionsScalarArray_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Ptr{Float64},arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Ptr{PetscBool})
    ccall((:PetscOptionsScalarArray_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Float64},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscOptionsIntArray_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Ptr{PetscBool})
    ccall((:PetscOptionsIntArray_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscOptionsStringArray_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Ptr{Ptr{Uint8}},arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Ptr{PetscBool})
    ccall((:PetscOptionsStringArray_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Ptr{Uint8}},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscOptionsBoolArray_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Ptr{PetscBool},arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Ptr{PetscBool})
    ccall((:PetscOptionsBoolArray_Private,petsc),PetscErrorCode,(Ptr{PetscOptions},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{PetscBool},Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscOptionsEnumArray_Private(arg1::Ptr{PetscOptions},arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Ptr{Ptr{Uint8}},arg6::Ptr{PetscEnum},arg7::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg8::Ptr{PetscBool})
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

function PetscObjectsDump(arg1::Ptr{FILE},arg2::PetscBool)
    ccall((:PetscObjectsDump,petsc),PetscErrorCode,(Ptr{FILE},PetscBool),arg1,arg2)
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

function PetscFunctionListPrintTypes(arg1::MPI_Comm,arg2::Ptr{FILE},arg3::Ptr{Uint8},arg4::Ptr{Uint8},arg5::Ptr{Uint8},arg6::Ptr{Uint8},arg7::PetscFunctionList,arg8::Ptr{Uint8})
    ccall((:PetscFunctionListPrintTypes,petsc),PetscErrorCode,(MPI_Comm,Ptr{FILE},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},PetscFunctionList,Ptr{Uint8}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
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

function PetscSplitOwnership(arg1::MPI_Comm,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSplitOwnership,petsc),PetscErrorCode,(MPI_Comm,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscSplitOwnershipBlock(arg1::MPI_Comm,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
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

function PetscMPIDump(arg1::Ptr{FILE})
    ccall((:PetscMPIDump,petsc),PetscErrorCode,(Ptr{FILE},),arg1)
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

function PetscFOpen(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Ptr{FILE}})
    ccall((:PetscFOpen,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Ptr{Ptr{FILE}}),arg1,arg2,arg3,arg4)
end

function PetscFClose(arg1::MPI_Comm,arg2::Ptr{FILE})
    ccall((:PetscFClose,petsc),PetscErrorCode,(MPI_Comm,Ptr{FILE}),arg1,arg2)
end

function PetscVSNPrintf(arg1::Ptr{Uint8},size_t::Cint,arg2::Ptr{Uint8},arg3::Ptr{Cint},va_list::Cint)
    ccall((:PetscVSNPrintf,petsc),PetscErrorCode,(Ptr{Uint8},Cint,Ptr{Uint8},Ptr{Cint},Cint),arg1,size_t,arg2,arg3,va_list)
end

function PetscVFPrintfDefault(arg1::Ptr{FILE},arg2::Ptr{Uint8},va_list::Cint)
    ccall((:PetscVFPrintfDefault,petsc),PetscErrorCode,(Ptr{FILE},Ptr{Uint8},Cint),arg1,arg2,va_list)
end

function PetscSynchronizedFlush(arg1::MPI_Comm,arg2::Ptr{FILE})
    ccall((:PetscSynchronizedFlush,petsc),PetscErrorCode,(MPI_Comm,Ptr{FILE}),arg1,arg2)
end

function PetscSynchronizedFGets(arg1::MPI_Comm,arg2::Ptr{FILE},size_t::Cint,arg3::Ptr{Uint8})
    ccall((:PetscSynchronizedFGets,petsc),PetscErrorCode,(MPI_Comm,Ptr{FILE},Cint,Ptr{Uint8}),arg1,arg2,size_t,arg3)
end

function PetscStartMatlab(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Ptr{FILE}})
    ccall((:PetscStartMatlab,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Ptr{Ptr{FILE}}),arg1,arg2,arg3,arg4)
end

function PetscStartJava(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{Uint8},arg4::Ptr{Ptr{FILE}})
    ccall((:PetscStartJava,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{Uint8},Ptr{Ptr{FILE}}),arg1,arg2,arg3,arg4)
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

function PetscIntView(arg1::Int32,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::PetscViewer)
    ccall((:PetscIntView,petsc),PetscErrorCode,(PetscInt,Ptr{PetscInt},PetscViewer),arg1,arg2,arg3)
end

function PetscRealView(arg1::Int32,Float64::Ptr{Cint},arg2::PetscViewer)
    ccall((:PetscRealView,petsc),PetscErrorCode,(PetscInt,Ptr{Cint},PetscViewer),arg1,PetscReal,arg2)
end

function PetscScalarView(arg1::Int32,arg2::Ptr{Float64},arg3::PetscViewer)
    ccall((:PetscScalarView,petsc),PetscErrorCode,(PetscInt,Ptr{Float64},PetscViewer),arg1,arg2,arg3)
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

function PetscSortInt(arg1::Int32,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSortInt,petsc),PetscErrorCode,(PetscInt,Ptr{PetscInt}),arg1,arg2)
end

function PetscSortRemoveDupsInt(arg1::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSortRemoveDupsInt,petsc),PetscErrorCode,(Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2)
end

function PetscFindInt(arg1::Int32,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscFindInt,petsc),PetscErrorCode,(PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function PetscSortIntWithPermutation(arg1::Int32,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSortIntWithPermutation,petsc),PetscErrorCode,(PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscSortStrWithPermutation(arg1::Int32,arg2::Ptr{Ptr{Uint8}},arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSortStrWithPermutation,petsc),PetscErrorCode,(PetscInt,Ptr{Ptr{Uint8}},Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscSortIntWithArray(arg1::Int32,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSortIntWithArray,petsc),PetscErrorCode,(PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscSortIntWithArrayPair(arg1::Int32,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSortIntWithArrayPair,petsc),PetscErrorCode,(PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function PetscSortMPIInt(arg1::Int32,arg2::Ptr{PetscMPIInt})
    ccall((:PetscSortMPIInt,petsc),PetscErrorCode,(PetscInt,Ptr{PetscMPIInt}),arg1,arg2)
end

function PetscSortRemoveDupsMPIInt(arg1::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg2::Ptr{PetscMPIInt})
    ccall((:PetscSortRemoveDupsMPIInt,petsc),PetscErrorCode,(Ptr{PetscInt},Ptr{PetscMPIInt}),arg1,arg2)
end

function PetscSortMPIIntWithArray(arg1::PetscMPIInt,arg2::Ptr{PetscMPIInt},arg3::Ptr{PetscMPIInt})
    ccall((:PetscSortMPIIntWithArray,petsc),PetscErrorCode,(PetscMPIInt,Ptr{PetscMPIInt},Ptr{PetscMPIInt}),arg1,arg2,arg3)
end

function PetscSortIntWithScalarArray(arg1::Int32,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Float64})
    ccall((:PetscSortIntWithScalarArray,petsc),PetscErrorCode,(PetscInt,Ptr{PetscInt},Ptr{Float64}),arg1,arg2,arg3)
end

function PetscSortIntWithDataArray(arg1::Int32,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Void},size_t::Cint,arg4::Ptr{Void})
    ccall((:PetscSortIntWithDataArray,petsc),PetscErrorCode,(PetscInt,Ptr{PetscInt},Ptr{Void},Cint,Ptr{Void}),arg1,arg2,arg3,size_t,arg4)
end

function PetscSortReal(arg1::Int32,Float64::Ptr{Cint})
    ccall((:PetscSortReal,petsc),PetscErrorCode,(PetscInt,Ptr{Cint}),arg1,PetscReal)
end

function PetscSortRealWithPermutation(arg1::Int32,Float64::Ptr{Cint},arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSortRealWithPermutation,petsc),PetscErrorCode,(PetscInt,Ptr{Cint},Ptr{PetscInt}),arg1,PetscReal,arg2)
end

function PetscSortRemoveDupsReal(arg1::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),Float64::Ptr{Cint})
    ccall((:PetscSortRemoveDupsReal,petsc),PetscErrorCode,(Ptr{PetscInt},Ptr{Cint}),arg1,PetscReal)
end

function PetscSortSplit(arg1::Int32,arg2::Int32,arg3::Ptr{Float64},arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSortSplit,petsc),PetscErrorCode,(PetscInt,PetscInt,Ptr{Float64},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function PetscSortSplitReal(arg1::Int32,arg2::Int32,Float64::Ptr{Cint},arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSortSplitReal,petsc),PetscErrorCode,(PetscInt,PetscInt,Ptr{Cint},Ptr{PetscInt}),arg1,arg2,PetscReal,arg3)
end

function PetscProcessTree(arg1::Int32,arg2::Ptr{PetscBool},arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Ptr{Ptr{Int32}},arg6::Ptr{Ptr{Int32}},arg7::Ptr{Ptr{Int32}},arg8::Ptr{Ptr{Int32}})
    ccall((:PetscProcessTree,petsc),PetscErrorCode,(PetscInt,Ptr{PetscBool},Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function PetscMergeIntArrayPair(arg1::Int32,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Int32,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg8::Ptr{Ptr{Int32}},arg9::Ptr{Ptr{Int32}})
    ccall((:PetscMergeIntArrayPair,petsc),PetscErrorCode,(PetscInt,Ptr{PetscInt},Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function PetscMergeIntArray(arg1::Int32,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Ptr{Ptr{Int32}})
    ccall((:PetscMergeIntArray,petsc),PetscErrorCode,(PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{Int32}}),arg1,arg2,arg3,arg4,arg5,arg6)
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

function PetscRandomGetValue(arg1::PetscRandom,arg2::Ptr{Float64})
    ccall((:PetscRandomGetValue,petsc),PetscErrorCode,(PetscRandom,Ptr{Float64}),arg1,arg2)
end

function PetscRandomGetValueReal(arg1::PetscRandom,arg2::Ptr{Cint})
    ccall((:PetscRandomGetValueReal,petsc),PetscErrorCode,(PetscRandom,Ptr{Cint}),arg1,arg2)
end

function PetscRandomGetInterval(arg1::PetscRandom,arg2::Ptr{Float64},arg3::Ptr{Float64})
    ccall((:PetscRandomGetInterval,petsc),PetscErrorCode,(PetscRandom,Ptr{Float64},Ptr{Float64}),arg1,arg2,arg3)
end

function PetscRandomSetInterval(arg1::PetscRandom,arg2::Float64,arg3::Float64)
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

function PetscBinaryRead(arg1::Cint,arg2::Ptr{Void},arg3::Int32,arg4::PetscDataType)
    ccall((:PetscBinaryRead,petsc),PetscErrorCode,(Cint,Ptr{Void},PetscInt,PetscDataType),arg1,arg2,arg3,arg4)
end

function PetscBinarySynchronizedRead(arg1::MPI_Comm,arg2::Cint,arg3::Ptr{Void},arg4::Int32,arg5::PetscDataType)
    ccall((:PetscBinarySynchronizedRead,petsc),PetscErrorCode,(MPI_Comm,Cint,Ptr{Void},PetscInt,PetscDataType),arg1,arg2,arg3,arg4,arg5)
end

function PetscBinarySynchronizedWrite(arg1::MPI_Comm,arg2::Cint,arg3::Ptr{Void},arg4::Int32,arg5::PetscDataType,arg6::PetscBool)
    ccall((:PetscBinarySynchronizedWrite,petsc),PetscErrorCode,(MPI_Comm,Cint,Ptr{Void},PetscInt,PetscDataType,PetscBool),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscBinaryWrite(arg1::Cint,arg2::Ptr{Void},arg3::Int32,arg4::PetscDataType,arg5::PetscBool)
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

function PetscBinarySeek(arg1::Cint,off_t::Cint,arg2::PetscBinarySeekType,arg3::Ptr{Cint})
    ccall((:PetscBinarySeek,petsc),PetscErrorCode,(Cint,Cint,PetscBinarySeekType,Ptr{Cint}),arg1,off_t,arg2,arg3)
end

function PetscBinarySynchronizedSeek(arg1::MPI_Comm,arg2::Cint,off_t::Cint,arg3::PetscBinarySeekType,arg4::Ptr{Cint})
    ccall((:PetscBinarySynchronizedSeek,petsc),PetscErrorCode,(MPI_Comm,Cint,Cint,PetscBinarySeekType,Ptr{Cint}),arg1,arg2,off_t,arg3,arg4)
end

function PetscByteSwap(arg1::Ptr{Void},arg2::PetscDataType,arg3::Int32)
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

function PetscPostIrecvInt(arg1::MPI_Comm,arg2::PetscMPIInt,arg3::PetscMPIInt,arg4::Ptr{PetscMPIInt},arg5::Ptr{PetscMPIInt},arg6::Ptr{Ptr{Ptr{Int32}}},arg7::Ptr{Ptr{MPI_Request}})
    ccall((:PetscPostIrecvInt,petsc),PetscErrorCode,(MPI_Comm,PetscMPIInt,PetscMPIInt,Ptr{PetscMPIInt},Ptr{PetscMPIInt},Ptr{Ptr{Ptr{Int32}}},Ptr{Ptr{MPI_Request}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscPostIrecvScalar(arg1::MPI_Comm,arg2::PetscMPIInt,arg3::PetscMPIInt,arg4::Ptr{PetscMPIInt},arg5::Ptr{PetscMPIInt},arg6::Ptr{Ptr{Ptr{Float64}}},arg7::Ptr{Ptr{MPI_Request}})
    ccall((:PetscPostIrecvScalar,petsc),PetscErrorCode,(MPI_Comm,PetscMPIInt,PetscMPIInt,Ptr{PetscMPIInt},Ptr{PetscMPIInt},Ptr{Ptr{Ptr{Float64}}},Ptr{Ptr{MPI_Request}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function PetscCommBuildTwoSided(arg1::MPI_Comm,arg2::PetscMPIInt,arg3::MPI_Datatype,arg4::Int32,arg5::Ptr{PetscMPIInt},arg6::Ptr{Void},arg7::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg8::Ptr{Ptr{PetscMPIInt}},arg9::Ptr{Void})
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

function PetscSubcommSetNumber(arg1::PetscSubcomm,arg2::Int32)
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

function ISCreateGeneral(arg1::MPI_Comm,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),PetscCopyMode::Cint,arg4::Ptr{IS})
    ccall((:ISCreateGeneral,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{PetscInt},Cint,Ptr{IS}),arg1,arg2,arg3,PetscCopyMode,arg4)
end

function ISGeneralSetIndices(arg1::IS,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),PetscCopyMode::Cint)
    ccall((:ISGeneralSetIndices,petsc),PetscErrorCode,(IS,PetscInt,Ptr{PetscInt},Cint),arg1,arg2,arg3,PetscCopyMode)
end

function ISCreateBlock(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),PetscCopyMode::Cint,arg5::Ptr{IS})
    ccall((:ISCreateBlock,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,Ptr{PetscInt},Cint,Ptr{IS}),arg1,arg2,arg3,arg4,PetscCopyMode,arg5)
end

function ISBlockSetIndices(arg1::IS,arg2::Int32,arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),PetscCopyMode::Cint)
    ccall((:ISBlockSetIndices,petsc),PetscErrorCode,(IS,PetscInt,PetscInt,Ptr{PetscInt},Cint),arg1,arg2,arg3,arg4,PetscCopyMode)
end

function ISCreateStride(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Ptr{IS})
    ccall((:ISCreateStride,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,Ptr{IS}),arg1,arg2,arg3,arg4,arg5)
end

function ISStrideSetStride(arg1::IS,arg2::Int32,arg3::Int32,arg4::Int32)
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

function ISContiguousLocal(arg1::IS,arg2::Int32,arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Ptr{PetscBool})
    ccall((:ISContiguousLocal,petsc),PetscErrorCode,(IS,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5)
end

function ISGetIndices(arg1::IS,arg2::Ptr{Ptr{Int32}})
    ccall((:ISGetIndices,petsc),PetscErrorCode,(IS,Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISRestoreIndices(arg1::IS,arg2::Ptr{Ptr{Int32}})
    ccall((:ISRestoreIndices,petsc),PetscErrorCode,(IS,Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISGetTotalIndices(arg1::IS,arg2::Ptr{Ptr{Int32}})
    ccall((:ISGetTotalIndices,petsc),PetscErrorCode,(IS,Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISRestoreTotalIndices(arg1::IS,arg2::Ptr{Ptr{Int32}})
    ccall((:ISRestoreTotalIndices,petsc),PetscErrorCode,(IS,Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISGetNonlocalIndices(arg1::IS,arg2::Ptr{Ptr{Int32}})
    ccall((:ISGetNonlocalIndices,petsc),PetscErrorCode,(IS,Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISRestoreNonlocalIndices(arg1::IS,arg2::Ptr{Ptr{Int32}})
    ccall((:ISRestoreNonlocalIndices,petsc),PetscErrorCode,(IS,Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISGetNonlocalIS(arg1::IS,is::Ptr{IS})
    ccall((:ISGetNonlocalIS,petsc),PetscErrorCode,(IS,Ptr{IS}),arg1,is)
end

function ISRestoreNonlocalIS(arg1::IS,is::Ptr{IS})
    ccall((:ISRestoreNonlocalIS,petsc),PetscErrorCode,(IS,Ptr{IS}),arg1,is)
end

function ISGetSize(arg1::IS,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISGetSize,petsc),PetscErrorCode,(IS,Ptr{PetscInt}),arg1,arg2)
end

function ISGetLocalSize(arg1::IS,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISGetLocalSize,petsc),PetscErrorCode,(IS,Ptr{PetscInt}),arg1,arg2)
end

function ISInvertPermutation(arg1::IS,arg2::Int32,arg3::Ptr{IS})
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

function ISGetMinMax(arg1::IS,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISGetMinMax,petsc),PetscErrorCode,(IS,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function ISBlockGetIndices(arg1::IS,arg2::Ptr{Ptr{Int32}})
    ccall((:ISBlockGetIndices,petsc),PetscErrorCode,(IS,Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISBlockRestoreIndices(arg1::IS,arg2::Ptr{Ptr{Int32}})
    ccall((:ISBlockRestoreIndices,petsc),PetscErrorCode,(IS,Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISBlockGetLocalSize(arg1::IS,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISBlockGetLocalSize,petsc),PetscErrorCode,(IS,Ptr{PetscInt}),arg1,arg2)
end

function ISBlockGetSize(arg1::IS,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISBlockGetSize,petsc),PetscErrorCode,(IS,Ptr{PetscInt}),arg1,arg2)
end

function ISGetBlockSize(arg1::IS,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISGetBlockSize,petsc),PetscErrorCode,(IS,Ptr{PetscInt}),arg1,arg2)
end

function ISSetBlockSize(arg1::IS,arg2::Int32)
    ccall((:ISSetBlockSize,petsc),PetscErrorCode,(IS,PetscInt),arg1,arg2)
end

function ISStrideGetInfo(arg1::IS,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
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

function ISComplement(arg1::IS,arg2::Int32,arg3::Int32,arg4::Ptr{IS})
    ccall((:ISComplement,petsc),PetscErrorCode,(IS,PetscInt,PetscInt,Ptr{IS}),arg1,arg2,arg3,arg4)
end

function ISConcatenate(arg1::MPI_Comm,arg2::Int32,arg3::Ptr{IS},arg4::Ptr{IS})
    ccall((:ISConcatenate,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{IS},Ptr{IS}),arg1,arg2,arg3,arg4)
end

function ISListToPair(arg1::MPI_Comm,arg2::Int32,arg3::Ptr{IS},arg4::Ptr{IS},arg5::Ptr{IS})
    ccall((:ISListToPair,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{IS},Ptr{IS},Ptr{IS}),arg1,arg2,arg3,arg4,arg5)
end

function ISPairToList(arg1::IS,arg2::IS,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{Ptr{IS}})
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

function ISLocalToGlobalMappingCreate(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),PetscCopyMode::Cint,arg5::Ptr{ISLocalToGlobalMapping})
    ccall((:ISLocalToGlobalMappingCreate,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,Ptr{PetscInt},Cint,Ptr{ISLocalToGlobalMapping}),arg1,arg2,arg3,arg4,PetscCopyMode,arg5)
end

function ISLocalToGlobalMappingCreateIS(arg1::IS,arg2::Ptr{ISLocalToGlobalMapping})
    ccall((:ISLocalToGlobalMappingCreateIS,petsc),PetscErrorCode,(IS,Ptr{ISLocalToGlobalMapping}),arg1,arg2)
end

function ISLocalToGlobalMappingCreateSF(arg1::PetscSF,arg2::Int32,arg3::Ptr{ISLocalToGlobalMapping})
    ccall((:ISLocalToGlobalMappingCreateSF,petsc),PetscErrorCode,(PetscSF,PetscInt,Ptr{ISLocalToGlobalMapping}),arg1,arg2,arg3)
end

function ISLocalToGlobalMappingView(arg1::ISLocalToGlobalMapping,arg2::PetscViewer)
    ccall((:ISLocalToGlobalMappingView,petsc),PetscErrorCode,(ISLocalToGlobalMapping,PetscViewer),arg1,arg2)
end

function ISLocalToGlobalMappingDestroy(arg1::Ptr{ISLocalToGlobalMapping})
    ccall((:ISLocalToGlobalMappingDestroy,petsc),PetscErrorCode,(Ptr{ISLocalToGlobalMapping},),arg1)
end

function ISLocalToGlobalMappingApply(arg1::ISLocalToGlobalMapping,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISLocalToGlobalMappingApply,petsc),PetscErrorCode,(ISLocalToGlobalMapping,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function ISLocalToGlobalMappingApplyBlock(arg1::ISLocalToGlobalMapping,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISLocalToGlobalMappingApplyBlock,petsc),PetscErrorCode,(ISLocalToGlobalMapping,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function ISLocalToGlobalMappingApplyIS(arg1::ISLocalToGlobalMapping,arg2::IS,arg3::Ptr{IS})
    ccall((:ISLocalToGlobalMappingApplyIS,petsc),PetscErrorCode,(ISLocalToGlobalMapping,IS,Ptr{IS}),arg1,arg2,arg3)
end

function ISGlobalToLocalMappingApply(arg1::ISLocalToGlobalMapping,arg2::ISGlobalToLocalMappingType,arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISGlobalToLocalMappingApply,petsc),PetscErrorCode,(ISLocalToGlobalMapping,ISGlobalToLocalMappingType,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function ISGlobalToLocalMappingApplyBlock(arg1::ISLocalToGlobalMapping,arg2::ISGlobalToLocalMappingType,arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISGlobalToLocalMappingApplyBlock,petsc),PetscErrorCode,(ISLocalToGlobalMapping,ISGlobalToLocalMappingType,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function ISGlobalToLocalMappingApplyIS(arg1::ISLocalToGlobalMapping,arg2::ISGlobalToLocalMappingType,arg3::IS,arg4::Ptr{IS})
    ccall((:ISGlobalToLocalMappingApplyIS,petsc),PetscErrorCode,(ISLocalToGlobalMapping,ISGlobalToLocalMappingType,IS,Ptr{IS}),arg1,arg2,arg3,arg4)
end

function ISLocalToGlobalMappingGetSize(arg1::ISLocalToGlobalMapping,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISLocalToGlobalMappingGetSize,petsc),PetscErrorCode,(ISLocalToGlobalMapping,Ptr{PetscInt}),arg1,arg2)
end

function ISLocalToGlobalMappingGetInfo(arg1::ISLocalToGlobalMapping,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Ptr{Int32}},arg4::Ptr{Ptr{Int32}},arg5::Ptr{Ptr{Ptr{Int32}}})
    ccall((:ISLocalToGlobalMappingGetInfo,petsc),PetscErrorCode,(ISLocalToGlobalMapping,Ptr{PetscInt},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}},Ptr{Ptr{Ptr{Int32}}}),arg1,arg2,arg3,arg4,arg5)
end

function ISLocalToGlobalMappingRestoreInfo(arg1::ISLocalToGlobalMapping,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Ptr{Int32}},arg4::Ptr{Ptr{Int32}},arg5::Ptr{Ptr{Ptr{Int32}}})
    ccall((:ISLocalToGlobalMappingRestoreInfo,petsc),PetscErrorCode,(ISLocalToGlobalMapping,Ptr{PetscInt},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}},Ptr{Ptr{Ptr{Int32}}}),arg1,arg2,arg3,arg4,arg5)
end

function ISLocalToGlobalMappingGetBlockInfo(arg1::ISLocalToGlobalMapping,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Ptr{Int32}},arg4::Ptr{Ptr{Int32}},arg5::Ptr{Ptr{Ptr{Int32}}})
    ccall((:ISLocalToGlobalMappingGetBlockInfo,petsc),PetscErrorCode,(ISLocalToGlobalMapping,Ptr{PetscInt},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}},Ptr{Ptr{Ptr{Int32}}}),arg1,arg2,arg3,arg4,arg5)
end

function ISLocalToGlobalMappingRestoreBlockInfo(arg1::ISLocalToGlobalMapping,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Ptr{Int32}},arg4::Ptr{Ptr{Int32}},arg5::Ptr{Ptr{Ptr{Int32}}})
    ccall((:ISLocalToGlobalMappingRestoreBlockInfo,petsc),PetscErrorCode,(ISLocalToGlobalMapping,Ptr{PetscInt},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}},Ptr{Ptr{Ptr{Int32}}}),arg1,arg2,arg3,arg4,arg5)
end

function ISLocalToGlobalMappingGetIndices(arg1::ISLocalToGlobalMapping,arg2::Ptr{Ptr{Int32}})
    ccall((:ISLocalToGlobalMappingGetIndices,petsc),PetscErrorCode,(ISLocalToGlobalMapping,Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISLocalToGlobalMappingRestoreIndices(arg1::ISLocalToGlobalMapping,arg2::Ptr{Ptr{Int32}})
    ccall((:ISLocalToGlobalMappingRestoreIndices,petsc),PetscErrorCode,(ISLocalToGlobalMapping,Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISLocalToGlobalMappingGetBlockIndices(arg1::ISLocalToGlobalMapping,arg2::Ptr{Ptr{Int32}})
    ccall((:ISLocalToGlobalMappingGetBlockIndices,petsc),PetscErrorCode,(ISLocalToGlobalMapping,Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISLocalToGlobalMappingRestoreBlockIndices(arg1::ISLocalToGlobalMapping,arg2::Ptr{Ptr{Int32}})
    ccall((:ISLocalToGlobalMappingRestoreBlockIndices,petsc),PetscErrorCode,(ISLocalToGlobalMapping,Ptr{Ptr{Int32}}),arg1,arg2)
end

function ISLocalToGlobalMappingConcatenate(arg1::MPI_Comm,arg2::Int32,arg3::Ptr{ISLocalToGlobalMapping},arg4::Ptr{ISLocalToGlobalMapping})
    ccall((:ISLocalToGlobalMappingConcatenate,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{ISLocalToGlobalMapping},Ptr{ISLocalToGlobalMapping}),arg1,arg2,arg3,arg4)
end

function ISG2LMapApply(arg1::ISLocalToGlobalMapping,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISG2LMapApply,petsc),PetscErrorCode,(ISLocalToGlobalMapping,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function ISLocalToGlobalMappingGetBlockSize(arg1::ISLocalToGlobalMapping,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISLocalToGlobalMappingGetBlockSize,petsc),PetscErrorCode,(ISLocalToGlobalMapping,Ptr{PetscInt}),arg1,arg2)
end

function ISAllGatherColors(arg1::MPI_Comm,arg2::Int32,arg3::Ptr{Cint},arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Ptr{Ptr{Cint}})
    ccall((:ISAllGatherColors,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{Cint},Ptr{PetscInt},Ptr{Ptr{Cint}}),arg1,arg2,arg3,arg4,arg5)
end

function ISColoringCreate(arg1::MPI_Comm,arg2::Int32,arg3::Int32,ISColoringValue::Ptr{Cint},PetscCopyMode::Cint,arg4::Ptr{ISColoring})
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

function ISColoringGetIS(arg1::ISColoring,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Ptr{IS}})
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

function ISPartitioningCount(arg1::IS,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:ISPartitioningCount,petsc),PetscErrorCode,(IS,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function ISCompressIndicesGeneral(arg1::Int32,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Ptr{IS},arg6::Ptr{IS})
    ccall((:ISCompressIndicesGeneral,petsc),PetscErrorCode,(PetscInt,PetscInt,PetscInt,PetscInt,Ptr{IS},Ptr{IS}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function ISCompressIndicesSorted(arg1::Int32,arg2::Int32,arg3::Int32,arg4::Ptr{IS},arg5::Ptr{IS})
    ccall((:ISCompressIndicesSorted,petsc),PetscErrorCode,(PetscInt,PetscInt,PetscInt,Ptr{IS},Ptr{IS}),arg1,arg2,arg3,arg4,arg5)
end

function ISExpandIndicesGeneral(arg1::Int32,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Ptr{IS},arg6::Ptr{IS})
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

function PetscLayoutSetLocalSize(arg1::PetscLayout,arg2::Int32)
    ccall((:PetscLayoutSetLocalSize,petsc),PetscErrorCode,(PetscLayout,PetscInt),arg1,arg2)
end

function PetscLayoutGetLocalSize(arg1::PetscLayout,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscLayoutGetLocalSize,petsc),PetscErrorCode,(PetscLayout,Ptr{PetscInt}),arg1,arg2)
end

function PetscLayoutSetSize(arg1::PetscLayout,arg2::Int32)
    ccall((:PetscLayoutSetSize,petsc),PetscErrorCode,(PetscLayout,PetscInt),arg1,arg2)
end

function PetscLayoutGetSize(arg1::PetscLayout,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscLayoutGetSize,petsc),PetscErrorCode,(PetscLayout,Ptr{PetscInt}),arg1,arg2)
end

function PetscLayoutSetBlockSize(arg1::PetscLayout,arg2::Int32)
    ccall((:PetscLayoutSetBlockSize,petsc),PetscErrorCode,(PetscLayout,PetscInt),arg1,arg2)
end

function PetscLayoutGetBlockSize(arg1::PetscLayout,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscLayoutGetBlockSize,petsc),PetscErrorCode,(PetscLayout,Ptr{PetscInt}),arg1,arg2)
end

function PetscLayoutGetRange(arg1::PetscLayout,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscLayoutGetRange,petsc),PetscErrorCode,(PetscLayout,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscLayoutGetRanges(arg1::PetscLayout,arg2::Ptr{Ptr{Int32}})
    ccall((:PetscLayoutGetRanges,petsc),PetscErrorCode,(PetscLayout,Ptr{Ptr{Int32}}),arg1,arg2)
end

function PetscLayoutSetISLocalToGlobalMapping(arg1::PetscLayout,arg2::ISLocalToGlobalMapping)
    ccall((:PetscLayoutSetISLocalToGlobalMapping,petsc),PetscErrorCode,(PetscLayout,ISLocalToGlobalMapping),arg1,arg2)
end

function PetscSFSetGraphLayout(arg1::PetscSF,arg2::PetscLayout,arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),PetscCopyMode::Cint,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
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

function PetscSectionGetNumFields(arg1::PetscSection,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetNumFields,petsc),PetscErrorCode,(PetscSection,Ptr{PetscInt}),arg1,arg2)
end

function PetscSectionSetNumFields(arg1::PetscSection,arg2::Int32)
    ccall((:PetscSectionSetNumFields,petsc),PetscErrorCode,(PetscSection,PetscInt),arg1,arg2)
end

function PetscSectionGetFieldName(arg1::PetscSection,arg2::Int32,arg3::Ptr{Ptr{Uint8}})
    ccall((:PetscSectionGetFieldName,petsc),PetscErrorCode,(PetscSection,PetscInt,Ptr{Ptr{Uint8}}),arg1,arg2,arg3)
end

function PetscSectionSetFieldName(arg1::PetscSection,arg2::Int32,arg3::Ptr{Uint8})
    ccall((:PetscSectionSetFieldName,petsc),PetscErrorCode,(PetscSection,PetscInt,Ptr{Uint8}),arg1,arg2,arg3)
end

function PetscSectionGetFieldComponents(arg1::PetscSection,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetFieldComponents,petsc),PetscErrorCode,(PetscSection,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscSectionSetFieldComponents(arg1::PetscSection,arg2::Int32,arg3::Int32)
    ccall((:PetscSectionSetFieldComponents,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt),arg1,arg2,arg3)
end

function PetscSectionGetChart(arg1::PetscSection,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetChart,petsc),PetscErrorCode,(PetscSection,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscSectionSetChart(arg1::PetscSection,arg2::Int32,arg3::Int32)
    ccall((:PetscSectionSetChart,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt),arg1,arg2,arg3)
end

function PetscSectionGetPermutation(arg1::PetscSection,arg2::Ptr{IS})
    ccall((:PetscSectionGetPermutation,petsc),PetscErrorCode,(PetscSection,Ptr{IS}),arg1,arg2)
end

function PetscSectionSetPermutation(arg1::PetscSection,arg2::IS)
    ccall((:PetscSectionSetPermutation,petsc),PetscErrorCode,(PetscSection,IS),arg1,arg2)
end

function PetscSectionGetDof(arg1::PetscSection,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetDof,petsc),PetscErrorCode,(PetscSection,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscSectionSetDof(arg1::PetscSection,arg2::Int32,arg3::Int32)
    ccall((:PetscSectionSetDof,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt),arg1,arg2,arg3)
end

function PetscSectionAddDof(arg1::PetscSection,arg2::Int32,arg3::Int32)
    ccall((:PetscSectionAddDof,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt),arg1,arg2,arg3)
end

function PetscSectionGetFieldDof(arg1::PetscSection,arg2::Int32,arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetFieldDof,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function PetscSectionSetFieldDof(arg1::PetscSection,arg2::Int32,arg3::Int32,arg4::Int32)
    ccall((:PetscSectionSetFieldDof,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end

function PetscSectionAddFieldDof(arg1::PetscSection,arg2::Int32,arg3::Int32,arg4::Int32)
    ccall((:PetscSectionAddFieldDof,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end

function PetscSectionHasConstraints(arg1::PetscSection,arg2::Ptr{PetscBool})
    ccall((:PetscSectionHasConstraints,petsc),PetscErrorCode,(PetscSection,Ptr{PetscBool}),arg1,arg2)
end

function PetscSectionGetConstraintDof(arg1::PetscSection,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetConstraintDof,petsc),PetscErrorCode,(PetscSection,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscSectionSetConstraintDof(arg1::PetscSection,arg2::Int32,arg3::Int32)
    ccall((:PetscSectionSetConstraintDof,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt),arg1,arg2,arg3)
end

function PetscSectionAddConstraintDof(arg1::PetscSection,arg2::Int32,arg3::Int32)
    ccall((:PetscSectionAddConstraintDof,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt),arg1,arg2,arg3)
end

function PetscSectionGetFieldConstraintDof(arg1::PetscSection,arg2::Int32,arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetFieldConstraintDof,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function PetscSectionSetFieldConstraintDof(arg1::PetscSection,arg2::Int32,arg3::Int32,arg4::Int32)
    ccall((:PetscSectionSetFieldConstraintDof,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end

function PetscSectionAddFieldConstraintDof(arg1::PetscSection,arg2::Int32,arg3::Int32,arg4::Int32)
    ccall((:PetscSectionAddFieldConstraintDof,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end

function PetscSectionGetConstraintIndices(arg1::PetscSection,arg2::Int32,arg3::Ptr{Ptr{Int32}})
    ccall((:PetscSectionGetConstraintIndices,petsc),PetscErrorCode,(PetscSection,PetscInt,Ptr{Ptr{Int32}}),arg1,arg2,arg3)
end

function PetscSectionSetConstraintIndices(arg1::PetscSection,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionSetConstraintIndices,petsc),PetscErrorCode,(PetscSection,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscSectionGetFieldConstraintIndices(arg1::PetscSection,arg2::Int32,arg3::Int32,arg4::Ptr{Ptr{Int32}})
    ccall((:PetscSectionGetFieldConstraintIndices,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt,Ptr{Ptr{Int32}}),arg1,arg2,arg3,arg4)
end

function PetscSectionSetFieldConstraintIndices(arg1::PetscSection,arg2::Int32,arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionSetFieldConstraintIndices,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function PetscSectionSetUpBC(arg1::PetscSection)
    ccall((:PetscSectionSetUpBC,petsc),PetscErrorCode,(PetscSection,),arg1)
end

function PetscSectionSetUp(arg1::PetscSection)
    ccall((:PetscSectionSetUp,petsc),PetscErrorCode,(PetscSection,),arg1)
end

function PetscSectionGetMaxDof(arg1::PetscSection,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetMaxDof,petsc),PetscErrorCode,(PetscSection,Ptr{PetscInt}),arg1,arg2)
end

function PetscSectionGetStorageSize(arg1::PetscSection,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetStorageSize,petsc),PetscErrorCode,(PetscSection,Ptr{PetscInt}),arg1,arg2)
end

function PetscSectionGetConstrainedStorageSize(arg1::PetscSection,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetConstrainedStorageSize,petsc),PetscErrorCode,(PetscSection,Ptr{PetscInt}),arg1,arg2)
end

function PetscSectionGetOffset(arg1::PetscSection,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetOffset,petsc),PetscErrorCode,(PetscSection,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function PetscSectionSetOffset(arg1::PetscSection,arg2::Int32,arg3::Int32)
    ccall((:PetscSectionSetOffset,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt),arg1,arg2,arg3)
end

function PetscSectionGetFieldOffset(arg1::PetscSection,arg2::Int32,arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscSectionGetFieldOffset,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function PetscSectionSetFieldOffset(arg1::PetscSection,arg2::Int32,arg3::Int32,arg4::Int32)
    ccall((:PetscSectionSetFieldOffset,petsc),PetscErrorCode,(PetscSection,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end

function PetscSectionGetOffsetRange(arg1::PetscSection,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
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

function PetscSectionCreateGlobalSectionCensored(arg1::PetscSection,arg2::PetscSF,arg3::PetscBool,arg4::Int32,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Ptr{PetscSection})
    ccall((:PetscSectionCreateGlobalSectionCensored,petsc),PetscErrorCode,(PetscSection,PetscSF,PetscBool,PetscInt,Ptr{PetscInt},Ptr{PetscSection}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function PetscSectionCreateSubsection(arg1::PetscSection,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{PetscSection})
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

function PetscSectionGetField(arg1::PetscSection,arg2::Int32,arg3::Ptr{PetscSection})
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

function PetscSFCreateRemoteOffsets(arg1::PetscSF,arg2::PetscSection,arg3::PetscSection,arg4::Ptr{Ptr{Int32}})
    ccall((:PetscSFCreateRemoteOffsets,petsc),PetscErrorCode,(PetscSF,PetscSection,PetscSection,Ptr{Ptr{Int32}}),arg1,arg2,arg3,arg4)
end

function PetscSFDistributeSection(arg1::PetscSF,arg2::PetscSection,arg3::Ptr{Ptr{Int32}},arg4::PetscSection)
    ccall((:PetscSFDistributeSection,petsc),PetscErrorCode,(PetscSF,PetscSection,Ptr{Ptr{Int32}},PetscSection),arg1,arg2,arg3,arg4)
end

function PetscSFCreateSectionSF(arg1::PetscSF,arg2::PetscSection,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::PetscSection,arg5::Ptr{PetscSF})
    ccall((:PetscSFCreateSectionSF,petsc),PetscErrorCode,(PetscSF,PetscSection,Ptr{PetscInt},PetscSection,Ptr{PetscSF}),arg1,arg2,arg3,arg4,arg5)
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

function PetscViewerASCIIOpenWithFILE(arg1::MPI_Comm,arg2::Ptr{FILE},arg3::Ptr{PetscViewer})
    ccall((:PetscViewerASCIIOpenWithFILE,petsc),PetscErrorCode,(MPI_Comm,Ptr{FILE},Ptr{PetscViewer}),arg1,arg2,arg3)
end

function PetscViewerASCIIOpen(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::Ptr{PetscViewer})
    ccall((:PetscViewerASCIIOpen,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},Ptr{PetscViewer}),arg1,arg2,arg3)
end

function PetscViewerASCIISetFILE(arg1::PetscViewer,arg2::Ptr{FILE})
    ccall((:PetscViewerASCIISetFILE,petsc),PetscErrorCode,(PetscViewer,Ptr{FILE}),arg1,arg2)
end

function PetscViewerBinaryOpen(arg1::MPI_Comm,arg2::Ptr{Uint8},arg3::PetscFileMode,arg4::Ptr{PetscViewer})
    ccall((:PetscViewerBinaryOpen,petsc),PetscErrorCode,(MPI_Comm,Ptr{Uint8},PetscFileMode,Ptr{PetscViewer}),arg1,arg2,arg3,arg4)
end

function PetscViewerBinaryGetFlowControl(arg1::PetscViewer,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscViewerBinaryGetFlowControl,petsc),PetscErrorCode,(PetscViewer,Ptr{PetscInt}),arg1,arg2)
end

function PetscViewerBinarySetFlowControl(arg1::PetscViewer,arg2::Int32)
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

function PetscViewerASCIIGetPointer(arg1::PetscViewer,arg2::Ptr{Ptr{FILE}})
    ccall((:PetscViewerASCIIGetPointer,petsc),PetscErrorCode,(PetscViewer,Ptr{Ptr{FILE}}),arg1,arg2)
end

function PetscViewerFileGetMode(arg1::PetscViewer,arg2::Ptr{PetscFileMode})
    ccall((:PetscViewerFileGetMode,petsc),PetscErrorCode,(PetscViewer,Ptr{PetscFileMode}),arg1,arg2)
end

function PetscViewerFileSetMode(arg1::PetscViewer,arg2::PetscFileMode)
    ccall((:PetscViewerFileSetMode,petsc),PetscErrorCode,(PetscViewer,PetscFileMode),arg1,arg2)
end

function PetscViewerRead(arg1::PetscViewer,arg2::Ptr{Void},arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::PetscDataType)
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

function PetscViewerASCIISetTab(arg1::PetscViewer,arg2::Int32)
    ccall((:PetscViewerASCIISetTab,petsc),PetscErrorCode,(PetscViewer,PetscInt),arg1,arg2)
end

function PetscViewerASCIIGetTab(arg1::PetscViewer,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PetscViewerASCIIGetTab,petsc),PetscErrorCode,(PetscViewer,Ptr{PetscInt}),arg1,arg2)
end

function PetscViewerASCIIAddTab(arg1::PetscViewer,arg2::Int32)
    ccall((:PetscViewerASCIIAddTab,petsc),PetscErrorCode,(PetscViewer,PetscInt),arg1,arg2)
end

function PetscViewerASCIISubtractTab(arg1::PetscViewer,arg2::Int32)
    ccall((:PetscViewerASCIISubtractTab,petsc),PetscErrorCode,(PetscViewer,PetscInt),arg1,arg2)
end

function PetscViewerASCIIRead(arg1::PetscViewer,arg2::Ptr{Void},arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::PetscDataType)
    ccall((:PetscViewerASCIIRead,petsc),PetscErrorCode,(PetscViewer,Ptr{Void},PetscInt,Ptr{PetscInt},PetscDataType),arg1,arg2,arg3,arg4,arg5)
end

function PetscViewerBinaryGetDescriptor(arg1::PetscViewer,arg2::Ptr{Cint})
    ccall((:PetscViewerBinaryGetDescriptor,petsc),PetscErrorCode,(PetscViewer,Ptr{Cint}),arg1,arg2)
end

function PetscViewerBinaryGetInfoPointer(arg1::PetscViewer,arg2::Ptr{Ptr{FILE}})
    ccall((:PetscViewerBinaryGetInfoPointer,petsc),PetscErrorCode,(PetscViewer,Ptr{Ptr{FILE}}),arg1,arg2)
end

function PetscViewerBinaryRead(arg1::PetscViewer,arg2::Ptr{Void},arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::PetscDataType)
    ccall((:PetscViewerBinaryRead,petsc),PetscErrorCode,(PetscViewer,Ptr{Void},PetscInt,Ptr{PetscInt},PetscDataType),arg1,arg2,arg3,arg4,arg5)
end

function PetscViewerBinaryWrite(arg1::PetscViewer,arg2::Ptr{Void},arg3::Int32,arg4::PetscDataType,arg5::PetscBool)
    ccall((:PetscViewerBinaryWrite,petsc),PetscErrorCode,(PetscViewer,Ptr{Void},PetscInt,PetscDataType,PetscBool),arg1,arg2,arg3,arg4,arg5)
end

function PetscViewerStringSetString(arg1::PetscViewer,arg2::Ptr{Uint8},arg3::Int32)
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

function PetscViewerDrawSetPause(arg1::PetscViewer,Float64::Cint)
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

function PetscViewerDrawSetBounds(arg1::PetscViewer,arg2::Int32,arg3::Ptr{Cint})
    ccall((:PetscViewerDrawSetBounds,petsc),PetscErrorCode,(PetscViewer,PetscInt,Ptr{Cint}),arg1,arg2,arg3)
end

function PetscViewerDrawGetBounds(arg1::PetscViewer,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Ptr{Cint}})
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

function PetscViewerVUGetPointer(arg1::PetscViewer,arg2::Ptr{Ptr{FILE}})
    ccall((:PetscViewerVUGetPointer,petsc),PetscErrorCode,(PetscViewer,Ptr{Ptr{FILE}}),arg1,arg2)
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

function PetscViewerMatlabGetArray(arg1::PetscViewer,arg2::Cint,arg3::Cint,arg4::Ptr{Float64},arg5::Ptr{Uint8})
    ccall((:PetscViewerMatlabGetArray,petsc),PetscErrorCode,(PetscViewer,Cint,Cint,Ptr{Float64},Ptr{Uint8}),arg1,arg2,arg3,arg4,arg5)
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

function PetscViewersGetViewer(arg1::PetscViewers,arg2::Int32,arg3::Ptr{PetscViewer})
    ccall((:PetscViewersGetViewer,petsc),PetscErrorCode,(PetscViewers,PetscInt,Ptr{PetscViewer}),arg1,arg2,arg3)
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

function VecCreateSeq(arg1::MPI_Comm,arg2::Int32,arg3::Ptr{Vec})
    ccall((:VecCreateSeq,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{Vec}),arg1,arg2,arg3)
end

function VecCreateMPI(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Ptr{Vec})
    ccall((:VecCreateMPI,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,Ptr{Vec}),arg1,arg2,arg3,arg4)
end

function VecCreateSeqWithArray(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Ptr{Float64},arg5::Ptr{Vec})
    ccall((:VecCreateSeqWithArray,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,Ptr{Float64},Ptr{Vec}),arg1,arg2,arg3,arg4,arg5)
end

function VecCreateMPIWithArray(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Ptr{Float64},arg6::Ptr{Vec})
    ccall((:VecCreateMPIWithArray,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,Ptr{Float64},Ptr{Vec}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecCreateShared(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Ptr{Vec})
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

function VecSetSizes(arg1::Vec,arg2::Int32,arg3::Int32)
    ccall((:VecSetSizes,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt),arg1,arg2,arg3)
end

function VecDotNorm2(arg1::Vec,arg2::Vec,arg3::Ptr{Float64},arg4::Ptr{Cint})
    ccall((:VecDotNorm2,petsc),PetscErrorCode,(Vec,Vec,Ptr{Float64},Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function VecDot(arg1::Vec,arg2::Vec,arg3::Ptr{Float64})
    ccall((:VecDot,petsc),PetscErrorCode,(Vec,Vec,Ptr{Float64}),arg1,arg2,arg3)
end

function VecDotRealPart(arg1::Vec,arg2::Vec,arg3::Ptr{Cint})
    ccall((:VecDotRealPart,petsc),PetscErrorCode,(Vec,Vec,Ptr{Cint}),arg1,arg2,arg3)
end

function VecTDot(arg1::Vec,arg2::Vec,arg3::Ptr{Float64})
    ccall((:VecTDot,petsc),PetscErrorCode,(Vec,Vec,Ptr{Float64}),arg1,arg2,arg3)
end

function VecMDot(arg1::Vec,arg2::Int32,arg3::Ptr{Vec},arg4::Ptr{Float64})
    ccall((:VecMDot,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{Vec},Ptr{Float64}),arg1,arg2,arg3,arg4)
end

function VecMTDot(arg1::Vec,arg2::Int32,arg3::Ptr{Vec},arg4::Ptr{Float64})
    ccall((:VecMTDot,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{Vec},Ptr{Float64}),arg1,arg2,arg3,arg4)
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

function VecSum(arg1::Vec,arg2::Ptr{Float64})
    ccall((:VecSum,petsc),PetscErrorCode,(Vec,Ptr{Float64}),arg1,arg2)
end

function VecMax(arg1::Vec,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Cint})
    ccall((:VecMax,petsc),PetscErrorCode,(Vec,Ptr{PetscInt},Ptr{Cint}),arg1,arg2,arg3)
end

function VecMin(arg1::Vec,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Cint})
    ccall((:VecMin,petsc),PetscErrorCode,(Vec,Ptr{PetscInt},Ptr{Cint}),arg1,arg2,arg3)
end

function VecScale(arg1::Vec,arg2::Float64)
    ccall((:VecScale,petsc),PetscErrorCode,(Vec,PetscScalar),arg1,arg2)
end

function VecCopy(arg1::Vec,arg2::Vec)
    ccall((:VecCopy,petsc),PetscErrorCode,(Vec,Vec),arg1,arg2)
end

function VecSetRandom(arg1::Vec,arg2::PetscRandom)
    ccall((:VecSetRandom,petsc),PetscErrorCode,(Vec,PetscRandom),arg1,arg2)
end

function VecSet(arg1::Vec,arg2::Float64)
    ccall((:VecSet,petsc),PetscErrorCode,(Vec,PetscScalar),arg1,arg2)
end

function VecSetInf(arg1::Vec)
    ccall((:VecSetInf,petsc),PetscErrorCode,(Vec,),arg1)
end

function VecSwap(arg1::Vec,arg2::Vec)
    ccall((:VecSwap,petsc),PetscErrorCode,(Vec,Vec),arg1,arg2)
end

function VecAXPY(arg1::Vec,arg2::Float64,arg3::Vec)
    ccall((:VecAXPY,petsc),PetscErrorCode,(Vec,PetscScalar,Vec),arg1,arg2,arg3)
end

function VecAXPBY(arg1::Vec,arg2::Float64,arg3::Float64,arg4::Vec)
    ccall((:VecAXPBY,petsc),PetscErrorCode,(Vec,PetscScalar,PetscScalar,Vec),arg1,arg2,arg3,arg4)
end

function VecMAXPY(arg1::Vec,arg2::Int32,arg3::Ptr{Float64},arg4::Ptr{Vec})
    ccall((:VecMAXPY,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{Float64},Ptr{Vec}),arg1,arg2,arg3,arg4)
end

function VecAYPX(arg1::Vec,arg2::Float64,arg3::Vec)
    ccall((:VecAYPX,petsc),PetscErrorCode,(Vec,PetscScalar,Vec),arg1,arg2,arg3)
end

function VecWAXPY(arg1::Vec,arg2::Float64,arg3::Vec,arg4::Vec)
    ccall((:VecWAXPY,petsc),PetscErrorCode,(Vec,PetscScalar,Vec,Vec),arg1,arg2,arg3,arg4)
end

function VecAXPBYPCZ(arg1::Vec,arg2::Float64,arg3::Float64,arg4::Float64,arg5::Vec,arg6::Vec)
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

function VecShift(arg1::Vec,arg2::Float64)
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

function VecDuplicateVecs(arg1::Vec,arg2::Int32,arg3::Ptr{Ptr{Vec}})
    ccall((:VecDuplicateVecs,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{Ptr{Vec}}),arg1,arg2,arg3)
end

function VecDestroyVecs(arg1::Int32,arg2::Ptr{Ptr{Vec}})
    ccall((:VecDestroyVecs,petsc),PetscErrorCode,(PetscInt,Ptr{Ptr{Vec}}),arg1,arg2)
end

function VecStrideNormAll(arg1::Vec,arg2::NormType,Float64::Ptr{Cint})
    ccall((:VecStrideNormAll,petsc),PetscErrorCode,(Vec,NormType,Ptr{Cint}),arg1,arg2,PetscReal)
end

function VecStrideMaxAll(arg1::Vec,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),Float64::Ptr{Cint})
    ccall((:VecStrideMaxAll,petsc),PetscErrorCode,(Vec,Ptr{PetscInt},Ptr{Cint}),arg1,arg2,PetscReal)
end

function VecStrideMinAll(arg1::Vec,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),Float64::Ptr{Cint})
    ccall((:VecStrideMinAll,petsc),PetscErrorCode,(Vec,Ptr{PetscInt},Ptr{Cint}),arg1,arg2,PetscReal)
end

function VecStrideScaleAll(arg1::Vec,arg2::Ptr{Float64})
    ccall((:VecStrideScaleAll,petsc),PetscErrorCode,(Vec,Ptr{Float64}),arg1,arg2)
end

function VecUniqueEntries(arg1::Vec,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Ptr{Float64}})
    ccall((:VecUniqueEntries,petsc),PetscErrorCode,(Vec,Ptr{PetscInt},Ptr{Ptr{Float64}}),arg1,arg2,arg3)
end

function VecStrideNorm(arg1::Vec,arg2::Int32,arg3::NormType,arg4::Ptr{Cint})
    ccall((:VecStrideNorm,petsc),PetscErrorCode,(Vec,PetscInt,NormType,Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function VecStrideMax(arg1::Vec,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{Cint})
    ccall((:VecStrideMax,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{PetscInt},Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function VecStrideMin(arg1::Vec,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{Cint})
    ccall((:VecStrideMin,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{PetscInt},Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function VecStrideScale(arg1::Vec,arg2::Int32,arg3::Float64)
    ccall((:VecStrideScale,petsc),PetscErrorCode,(Vec,PetscInt,PetscScalar),arg1,arg2,arg3)
end

function VecStrideSet(arg1::Vec,arg2::Int32,arg3::Float64)
    ccall((:VecStrideSet,petsc),PetscErrorCode,(Vec,PetscInt,PetscScalar),arg1,arg2,arg3)
end

function VecStrideGather(arg1::Vec,arg2::Int32,arg3::Vec,arg4::InsertMode)
    ccall((:VecStrideGather,petsc),PetscErrorCode,(Vec,PetscInt,Vec,InsertMode),arg1,arg2,arg3,arg4)
end

function VecStrideScatter(arg1::Vec,arg2::Int32,arg3::Vec,arg4::InsertMode)
    ccall((:VecStrideScatter,petsc),PetscErrorCode,(Vec,PetscInt,Vec,InsertMode),arg1,arg2,arg3,arg4)
end

function VecStrideGatherAll(arg1::Vec,arg2::Ptr{Vec},arg3::InsertMode)
    ccall((:VecStrideGatherAll,petsc),PetscErrorCode,(Vec,Ptr{Vec},InsertMode),arg1,arg2,arg3)
end

function VecStrideScatterAll(arg1::Ptr{Vec},arg2::Vec,arg3::InsertMode)
    ccall((:VecStrideScatterAll,petsc),PetscErrorCode,(Ptr{Vec},Vec,InsertMode),arg1,arg2,arg3)
end

function VecStrideSubSetScatter(arg1::Vec,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Vec,arg6::InsertMode)
    ccall((:VecStrideSubSetScatter,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Vec,InsertMode),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecStrideSubSetGather(arg1::Vec,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Vec,arg6::InsertMode)
    ccall((:VecStrideSubSetGather,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Vec,InsertMode),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecSetValues(arg1::Vec,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{Float64},arg5::InsertMode)
    ccall((:VecSetValues,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{PetscInt},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5)
end

function VecGetValues(arg1::Vec,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{Float64})
    ccall((:VecGetValues,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{PetscInt},Ptr{Float64}),arg1,arg2,arg3,arg4)
end

function VecAssemblyBegin(arg1::Vec)
    ccall((:VecAssemblyBegin,petsc),PetscErrorCode,(Vec,),arg1)
end

function VecAssemblyEnd(arg1::Vec)
    ccall((:VecAssemblyEnd,petsc),PetscErrorCode,(Vec,),arg1)
end

function VecStashSetInitialSize(arg1::Vec,arg2::Int32,arg3::Int32)
    ccall((:VecStashSetInitialSize,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt),arg1,arg2,arg3)
end

function VecStashView(arg1::Vec,arg2::PetscViewer)
    ccall((:VecStashView,petsc),PetscErrorCode,(Vec,PetscViewer),arg1,arg2)
end

function VecStashViewFromOptions(arg1::Vec,arg2::PetscObject,arg3::Ptr{Uint8})
    ccall((:VecStashViewFromOptions,petsc),PetscErrorCode,(Vec,PetscObject,Ptr{Uint8}),arg1,arg2,arg3)
end

function VecStashGetInfo(arg1::Vec,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:VecStashGetInfo,petsc),PetscErrorCode,(Vec,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end

function VecGetBlockSize(arg1::Vec,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:VecGetBlockSize,petsc),PetscErrorCode,(Vec,Ptr{PetscInt}),arg1,arg2)
end

function VecSetValuesBlocked(arg1::Vec,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{Float64},arg5::InsertMode)
    ccall((:VecSetValuesBlocked,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{PetscInt},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5)
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

function VecScatterCreateLocal(arg1::VecScatter,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Int32,arg7::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg8::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg9::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg10::Int32)
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

function VecGetArray4d(arg1::Vec,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Int32,arg7::Int32,arg8::Int32,arg9::Int32,arg10::Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}})
    ccall((:VecGetArray4d,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function VecRestoreArray4d(arg1::Vec,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Int32,arg7::Int32,arg8::Int32,arg9::Int32,arg10::Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}})
    ccall((:VecRestoreArray4d,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function VecGetArray3d(arg1::Vec,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Int32,arg7::Int32,arg8::Ptr{Ptr{Ptr{Ptr{Float64}}}})
    ccall((:VecGetArray3d,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{Float64}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function VecRestoreArray3d(arg1::Vec,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Int32,arg7::Int32,arg8::Ptr{Ptr{Ptr{Ptr{Float64}}}})
    ccall((:VecRestoreArray3d,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{Float64}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function VecGetArray2d(arg1::Vec,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Ptr{Ptr{Ptr{Float64}}})
    ccall((:VecGetArray2d,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Float64}}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecRestoreArray2d(arg1::Vec,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Ptr{Ptr{Ptr{Float64}}})
    ccall((:VecRestoreArray2d,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Float64}}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecGetArray1d(arg1::Vec,arg2::Int32,arg3::Int32,arg4::Ptr{Ptr{Float64}})
    ccall((:VecGetArray1d,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4)
end

function VecRestoreArray1d(arg1::Vec,arg2::Int32,arg3::Int32,arg4::Ptr{Ptr{Float64}})
    ccall((:VecRestoreArray1d,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4)
end

function VecGetArray4dRead(arg1::Vec,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Int32,arg7::Int32,arg8::Int32,arg9::Int32,arg10::Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}})
    ccall((:VecGetArray4dRead,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function VecRestoreArray4dRead(arg1::Vec,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Int32,arg7::Int32,arg8::Int32,arg9::Int32,arg10::Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}})
    ccall((:VecRestoreArray4dRead,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{Ptr{Float64}}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function VecGetArray3dRead(arg1::Vec,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Int32,arg7::Int32,arg8::Ptr{Ptr{Ptr{Ptr{Float64}}}})
    ccall((:VecGetArray3dRead,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{Float64}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function VecRestoreArray3dRead(arg1::Vec,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Int32,arg7::Int32,arg8::Ptr{Ptr{Ptr{Ptr{Float64}}}})
    ccall((:VecRestoreArray3dRead,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Ptr{Float64}}}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function VecGetArray2dRead(arg1::Vec,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Ptr{Ptr{Ptr{Float64}}})
    ccall((:VecGetArray2dRead,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Float64}}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecRestoreArray2dRead(arg1::Vec,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Ptr{Ptr{Ptr{Float64}}})
    ccall((:VecRestoreArray2dRead,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Ptr{Ptr{Float64}}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecGetArray1dRead(arg1::Vec,arg2::Int32,arg3::Int32,arg4::Ptr{Ptr{Float64}})
    ccall((:VecGetArray1dRead,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4)
end

function VecRestoreArray1dRead(arg1::Vec,arg2::Int32,arg3::Int32,arg4::Ptr{Ptr{Float64}})
    ccall((:VecRestoreArray1dRead,petsc),PetscErrorCode,(Vec,PetscInt,PetscInt,Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4)
end

function VecPlaceArray(arg1::Vec,arg2::Ptr{Float64})
    ccall((:VecPlaceArray,petsc),PetscErrorCode,(Vec,Ptr{Float64}),arg1,arg2)
end

function VecResetArray(arg1::Vec)
    ccall((:VecResetArray,petsc),PetscErrorCode,(Vec,),arg1)
end

function VecReplaceArray(arg1::Vec,arg2::Ptr{Float64})
    ccall((:VecReplaceArray,petsc),PetscErrorCode,(Vec,Ptr{Float64}),arg1,arg2)
end

function VecGetArrays(arg1::Ptr{Vec},arg2::Int32,arg3::Ptr{Ptr{Ptr{Float64}}})
    ccall((:VecGetArrays,petsc),PetscErrorCode,(Ptr{Vec},PetscInt,Ptr{Ptr{Ptr{Float64}}}),arg1,arg2,arg3)
end

function VecRestoreArrays(arg1::Ptr{Vec},arg2::Int32,arg3::Ptr{Ptr{Ptr{Float64}}})
    ccall((:VecRestoreArrays,petsc),PetscErrorCode,(Ptr{Vec},PetscInt,Ptr{Ptr{Ptr{Float64}}}),arg1,arg2,arg3)
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

function VecGetSize(arg1::Vec,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:VecGetSize,petsc),PetscErrorCode,(Vec,Ptr{PetscInt}),arg1,arg2)
end

function VecGetLocalSize(arg1::Vec,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:VecGetLocalSize,petsc),PetscErrorCode,(Vec,Ptr{PetscInt}),arg1,arg2)
end

function VecGetOwnershipRange(arg1::Vec,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:VecGetOwnershipRange,petsc),PetscErrorCode,(Vec,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function VecGetOwnershipRanges(arg1::Vec,arg2::Ptr{Ptr{Int32}})
    ccall((:VecGetOwnershipRanges,petsc),PetscErrorCode,(Vec,Ptr{Ptr{Int32}}),arg1,arg2)
end

function VecSetLocalToGlobalMapping(arg1::Vec,arg2::ISLocalToGlobalMapping)
    ccall((:VecSetLocalToGlobalMapping,petsc),PetscErrorCode,(Vec,ISLocalToGlobalMapping),arg1,arg2)
end

function VecSetValuesLocal(arg1::Vec,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{Float64},arg5::InsertMode)
    ccall((:VecSetValuesLocal,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{PetscInt},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5)
end

function VecGetLocalToGlobalMapping(arg1::Vec,arg2::Ptr{ISLocalToGlobalMapping})
    ccall((:VecGetLocalToGlobalMapping,petsc),PetscErrorCode,(Vec,Ptr{ISLocalToGlobalMapping}),arg1,arg2)
end

function VecDotBegin(arg1::Vec,arg2::Vec,arg3::Ptr{Float64})
    ccall((:VecDotBegin,petsc),PetscErrorCode,(Vec,Vec,Ptr{Float64}),arg1,arg2,arg3)
end

function VecDotEnd(arg1::Vec,arg2::Vec,arg3::Ptr{Float64})
    ccall((:VecDotEnd,petsc),PetscErrorCode,(Vec,Vec,Ptr{Float64}),arg1,arg2,arg3)
end

function VecTDotBegin(arg1::Vec,arg2::Vec,arg3::Ptr{Float64})
    ccall((:VecTDotBegin,petsc),PetscErrorCode,(Vec,Vec,Ptr{Float64}),arg1,arg2,arg3)
end

function VecTDotEnd(arg1::Vec,arg2::Vec,arg3::Ptr{Float64})
    ccall((:VecTDotEnd,petsc),PetscErrorCode,(Vec,Vec,Ptr{Float64}),arg1,arg2,arg3)
end

function VecNormBegin(arg1::Vec,arg2::NormType,arg3::Ptr{Cint})
    ccall((:VecNormBegin,petsc),PetscErrorCode,(Vec,NormType,Ptr{Cint}),arg1,arg2,arg3)
end

function VecNormEnd(arg1::Vec,arg2::NormType,arg3::Ptr{Cint})
    ccall((:VecNormEnd,petsc),PetscErrorCode,(Vec,NormType,Ptr{Cint}),arg1,arg2,arg3)
end

function VecMDotBegin(arg1::Vec,arg2::Int32,arg3::Ptr{Vec},arg4::Ptr{Float64})
    ccall((:VecMDotBegin,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{Vec},Ptr{Float64}),arg1,arg2,arg3,arg4)
end

function VecMDotEnd(arg1::Vec,arg2::Int32,arg3::Ptr{Vec},arg4::Ptr{Float64})
    ccall((:VecMDotEnd,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{Vec},Ptr{Float64}),arg1,arg2,arg3,arg4)
end

function VecMTDotBegin(arg1::Vec,arg2::Int32,arg3::Ptr{Vec},arg4::Ptr{Float64})
    ccall((:VecMTDotBegin,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{Vec},Ptr{Float64}),arg1,arg2,arg3,arg4)
end

function VecMTDotEnd(arg1::Vec,arg2::Int32,arg3::Ptr{Vec},arg4::Ptr{Float64})
    ccall((:VecMTDotEnd,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{Vec},Ptr{Float64}),arg1,arg2,arg3,arg4)
end

function PetscCommSplitReductionBegin(arg1::MPI_Comm)
    ccall((:PetscCommSplitReductionBegin,petsc),PetscErrorCode,(MPI_Comm,),arg1)
end

function VecSetOption(arg1::Vec,arg2::VecOption,arg3::PetscBool)
    ccall((:VecSetOption,petsc),PetscErrorCode,(Vec,VecOption,PetscBool),arg1,arg2,arg3)
end

function VecGetArray(arg1::Vec,arg2::Ptr{Ptr{Float64}})
    ccall((:VecGetArray,petsc),PetscErrorCode,(Vec,Ptr{Ptr{Float64}}),arg1,arg2)
end

function VecGetArrayRead(arg1::Vec,arg2::Ptr{Ptr{Float64}})
    ccall((:VecGetArrayRead,petsc),PetscErrorCode,(Vec,Ptr{Ptr{Float64}}),arg1,arg2)
end

function VecRestoreArray(arg1::Vec,arg2::Ptr{Ptr{Float64}})
    ccall((:VecRestoreArray,petsc),PetscErrorCode,(Vec,Ptr{Ptr{Float64}}),arg1,arg2)
end

function VecRestoreArrayRead(arg1::Vec,arg2::Ptr{Ptr{Float64}})
    ccall((:VecRestoreArrayRead,petsc),PetscErrorCode,(Vec,Ptr{Ptr{Float64}}),arg1,arg2)
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

function VecContourScale(arg1::Vec,Float64::Cint,arg2::Cint)
    ccall((:VecContourScale,petsc),PetscErrorCode,(Vec,Cint,Cint),arg1,PetscReal,arg2)
end

function VecSetOperation(arg1::Vec,arg2::VecOperation,arg3::Ptr{Void})
    ccall((:VecSetOperation,petsc),PetscErrorCode,(Vec,VecOperation,Ptr{Void}),arg1,arg2,arg3)
end

function VecMPISetGhost(arg1::Vec,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:VecMPISetGhost,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function VecCreateGhost(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Ptr{Vec})
    ccall((:VecCreateGhost,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Vec}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function VecCreateGhostWithArray(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Ptr{Float64},arg7::Ptr{Vec})
    ccall((:VecCreateGhostWithArray,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Float64},Ptr{Vec}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function VecCreateGhostBlock(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Ptr{Vec})
    ccall((:VecCreateGhostBlock,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Vec}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function VecCreateGhostBlockWithArray(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Ptr{Float64},arg8::Ptr{Vec})
    ccall((:VecCreateGhostBlockWithArray,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Float64},Ptr{Vec}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
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

function VecPow(arg1::Vec,arg2::Float64)
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

function VecISAXPY(arg1::Vec,arg2::IS,arg3::Float64,arg4::Vec)
    ccall((:VecISAXPY,petsc),PetscErrorCode,(Vec,IS,PetscScalar,Vec),arg1,arg2,arg3,arg4)
end

function VecISSet(arg1::Vec,arg2::IS,arg3::Float64)
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

function VecsCreateSeq(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Ptr{Vecs})
    ccall((:VecsCreateSeq,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,Ptr{Vecs}),arg1,arg2,arg3,arg4)
end

function VecsCreateSeqWithArray(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Ptr{Float64},arg5::Ptr{Vecs})
    ccall((:VecsCreateSeqWithArray,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,Ptr{Float64},Ptr{Vecs}),arg1,arg2,arg3,arg4,arg5)
end

function VecsDuplicate(arg1::Vecs,arg2::Ptr{Vecs})
    ccall((:VecsDuplicate,petsc),PetscErrorCode,(Vecs,Ptr{Vecs}),arg1,arg2)
end

function VecNestGetSubVecs(arg1::Vec,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Ptr{Vec}})
    ccall((:VecNestGetSubVecs,petsc),PetscErrorCode,(Vec,Ptr{PetscInt},Ptr{Ptr{Vec}}),arg1,arg2,arg3)
end

function VecNestGetSubVec(arg1::Vec,arg2::Int32,arg3::Ptr{Vec})
    ccall((:VecNestGetSubVec,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{Vec}),arg1,arg2,arg3)
end

function VecNestSetSubVecs(arg1::Vec,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{Vec})
    ccall((:VecNestSetSubVecs,petsc),PetscErrorCode,(Vec,PetscInt,Ptr{PetscInt},Ptr{Vec}),arg1,arg2,arg3,arg4)
end

function VecNestSetSubVec(arg1::Vec,arg2::Int32,arg3::Vec)
    ccall((:VecNestSetSubVec,petsc),PetscErrorCode,(Vec,PetscInt,Vec),arg1,arg2,arg3)
end

function VecCreateNest(arg1::MPI_Comm,arg2::Int32,arg3::Ptr{IS},arg4::Ptr{Vec},arg5::Ptr{Vec})
    ccall((:VecCreateNest,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{IS},Ptr{Vec},Ptr{Vec}),arg1,arg2,arg3,arg4,arg5)
end

function VecNestGetSize(arg1::Vec,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:VecNestGetSize,petsc),PetscErrorCode,(Vec,Ptr{PetscInt}),arg1,arg2)
end

function PetscOptionsGetVec(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Vec,arg4::Ptr{PetscBool})
    ccall((:PetscOptionsGetVec,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Vec,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function VecChop(arg1::Vec,Float64::Cint)
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

function VecGetValuesSection(arg1::Vec,arg2::PetscSection,arg3::Int32,arg4::Ptr{Ptr{Float64}})
    ccall((:VecGetValuesSection,petsc),PetscErrorCode,(Vec,PetscSection,PetscInt,Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4)
end

function VecSetValuesSection(arg1::Vec,arg2::PetscSection,arg3::Int32,arg4::Ptr{Float64},arg5::InsertMode)
    ccall((:VecSetValuesSection,petsc),PetscErrorCode,(Vec,PetscSection,PetscInt,Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5)
end

function PetscSectionVecNorm(arg1::PetscSection,arg2::PetscSection,arg3::Vec,arg4::NormType,Float64::Ptr{Cint})
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

function MatSetSizes(arg1::Mat,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32)
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

function MatCreateSeqDense(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Ptr{Float64},arg5::Ptr{Mat})
    ccall((:MatCreateSeqDense,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,Ptr{Float64},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5)
end

function MatCreateDense(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Ptr{Float64},arg7::Ptr{Mat})
    ccall((:MatCreateDense,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Float64},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateSeqAIJ(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Ptr{Mat})
    ccall((:MatCreateSeqAIJ,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatCreateAIJ(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Int32,arg7::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg8::Int32,arg9::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg10::Ptr{Mat})
    ccall((:MatCreateAIJ,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function MatCreateMPIAIJWithArrays(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg8::Ptr{Float64},arg9::Ptr{Mat})
    ccall((:MatCreateMPIAIJWithArrays,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function MatCreateMPIAIJWithSplitArrays(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg8::Ptr{Float64},arg9::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg10::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg11::Ptr{Float64},arg12::Ptr{Mat})
    ccall((:MatCreateMPIAIJWithSplitArrays,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64},Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12)
end

function MatCreateSeqBAIJ(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Ptr{Mat})
    ccall((:MatCreateSeqBAIJ,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateBAIJ(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Int32,arg7::Int32,arg8::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg9::Int32,arg10::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg11::Ptr{Mat})
    ccall((:MatCreateBAIJ,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
end

function MatCreateMPIBAIJWithArrays(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Int32,arg7::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg8::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg9::Ptr{Float64},arg10::Ptr{Mat})
    ccall((:MatCreateMPIBAIJWithArrays,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function MatCreateMPIAdj(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Ptr{Mat})
    ccall((:MatCreateMPIAdj,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateSeqSBAIJ(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Ptr{Mat})
    ccall((:MatCreateSeqSBAIJ,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateSBAIJ(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Int32,arg7::Int32,arg8::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg9::Int32,arg10::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg11::Ptr{Mat})
    ccall((:MatCreateSBAIJ,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
end

function MatCreateMPISBAIJWithArrays(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Int32,arg7::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg8::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg9::Ptr{Float64},arg10::Ptr{Mat})
    ccall((:MatCreateMPISBAIJWithArrays,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function MatSeqSBAIJSetPreallocationCSR(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Ptr{Float64})
    ccall((:MatSeqSBAIJSetPreallocationCSR,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5)
end

function MatMPISBAIJSetPreallocationCSR(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Ptr{Float64})
    ccall((:MatMPISBAIJSetPreallocationCSR,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5)
end

function MatXAIJSetPreallocation(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatXAIJSetPreallocation,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatCreateShell(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Ptr{Void},arg7::Ptr{Mat})
    ccall((:MatCreateShell,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Void},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateNormal(arg1::Mat,arg2::Ptr{Mat})
    ccall((:MatCreateNormal,petsc),PetscErrorCode,(Mat,Ptr{Mat}),arg1,arg2)
end

function MatCreateLRC(arg1::Mat,arg2::Mat,arg3::Mat,arg4::Ptr{Mat})
    ccall((:MatCreateLRC,petsc),PetscErrorCode,(Mat,Mat,Mat,Ptr{Mat}),arg1,arg2,arg3,arg4)
end

function MatCreateIS(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Int32,arg7::ISLocalToGlobalMapping,arg8::Ptr{Mat})
    ccall((:MatCreateIS,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,ISLocalToGlobalMapping,Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatCreateSeqAIJCRL(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Ptr{Mat})
    ccall((:MatCreateSeqAIJCRL,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatCreateMPIAIJCRL(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Int32,arg7::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg8::Ptr{Mat})
    ccall((:MatCreateMPIAIJCRL,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatCreateSeqBSTRM(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Ptr{Mat})
    ccall((:MatCreateSeqBSTRM,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateMPIBSTRM(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Int32,arg7::Int32,arg8::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg9::Int32,arg10::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg11::Ptr{Mat})
    ccall((:MatCreateMPIBSTRM,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
end

function MatCreateSeqSBSTRM(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Ptr{Mat})
    ccall((:MatCreateSeqSBSTRM,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateMPISBSTRM(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Int32,arg7::Int32,arg8::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg9::Int32,arg10::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg11::Ptr{Mat})
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

function MatCreateBlockMat(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Ptr{Mat})
    ccall((:MatCreateBlockMat,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCompositeAddMat(arg1::Mat,arg2::Mat)
    ccall((:MatCompositeAddMat,petsc),PetscErrorCode,(Mat,Mat),arg1,arg2)
end

function MatCompositeMerge(arg1::Mat)
    ccall((:MatCompositeMerge,petsc),PetscErrorCode,(Mat,),arg1)
end

function MatCreateComposite(arg1::MPI_Comm,arg2::Int32,arg3::Ptr{Mat},arg4::Ptr{Mat})
    ccall((:MatCreateComposite,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{Mat},Ptr{Mat}),arg1,arg2,arg3,arg4)
end

function MatCompositeSetType(arg1::Mat,arg2::MatCompositeType)
    ccall((:MatCompositeSetType,petsc),PetscErrorCode,(Mat,MatCompositeType),arg1,arg2)
end

function MatCreateFFT(arg1::MPI_Comm,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::MatType,arg5::Ptr{Mat})
    ccall((:MatCreateFFT,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{PetscInt},MatType,Ptr{Mat}),arg1,arg2,arg3,arg4,arg5)
end

function MatCreateSeqCUFFT(arg1::MPI_Comm,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{Mat})
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

function MatGetTrace(arg1::Mat,arg2::Ptr{Float64})
    ccall((:MatGetTrace,petsc),PetscErrorCode,(Mat,Ptr{Float64}),arg1,arg2)
end

function MatInvertBlockDiagonal(arg1::Mat,arg2::Ptr{Ptr{Float64}})
    ccall((:MatInvertBlockDiagonal,petsc),PetscErrorCode,(Mat,Ptr{Ptr{Float64}}),arg1,arg2)
end

function MatSetValues(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Int32,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Ptr{Float64},arg7::InsertMode)
    ccall((:MatSetValues,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatSetValuesBlocked(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Int32,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Ptr{Float64},arg7::InsertMode)
    ccall((:MatSetValuesBlocked,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatSetValuesRow(arg1::Mat,arg2::Int32,arg3::Ptr{Float64})
    ccall((:MatSetValuesRow,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{Float64}),arg1,arg2,arg3)
end

function MatSetValuesRowLocal(arg1::Mat,arg2::Int32,arg3::Ptr{Float64})
    ccall((:MatSetValuesRowLocal,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{Float64}),arg1,arg2,arg3)
end

function MatSetValuesBatch(arg1::Mat,arg2::Int32,arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Ptr{Float64})
    ccall((:MatSetValuesBatch,petsc),PetscErrorCode,(Mat,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5)
end

function MatSetRandom(arg1::Mat,arg2::PetscRandom)
    ccall((:MatSetRandom,petsc),PetscErrorCode,(Mat,PetscRandom),arg1,arg2)
end

function MatSetValuesStencil(arg1::Mat,arg2::Int32,arg3::Ptr{MatStencil},arg4::Int32,arg5::Ptr{MatStencil},arg6::Ptr{Float64},arg7::InsertMode)
    ccall((:MatSetValuesStencil,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{MatStencil},PetscInt,Ptr{MatStencil},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatSetValuesBlockedStencil(arg1::Mat,arg2::Int32,arg3::Ptr{MatStencil},arg4::Int32,arg5::Ptr{MatStencil},arg6::Ptr{Float64},arg7::InsertMode)
    ccall((:MatSetValuesBlockedStencil,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{MatStencil},PetscInt,Ptr{MatStencil},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatSetStencil(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Int32)
    ccall((:MatSetStencil,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{PetscInt},PetscInt),arg1,arg2,arg3,arg4,arg5)
end

function MatSetColoring(arg1::Mat,arg2::ISColoring)
    ccall((:MatSetColoring,petsc),PetscErrorCode,(Mat,ISColoring),arg1,arg2)
end

function MatSetValuesAdifor(arg1::Mat,arg2::Int32,arg3::Ptr{Void})
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

function MatGetValues(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Int32,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Ptr{Float64})
    ccall((:MatGetValues,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatGetRow(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{Ptr{Int32}},arg5::Ptr{Ptr{Float64}})
    ccall((:MatGetRow,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{Ptr{Int32}},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4,arg5)
end

function MatRestoreRow(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{Ptr{Int32}},arg5::Ptr{Ptr{Float64}})
    ccall((:MatRestoreRow,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{Ptr{Int32}},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4,arg5)
end

function MatGetRowUpperTriangular(arg1::Mat)
    ccall((:MatGetRowUpperTriangular,petsc),PetscErrorCode,(Mat,),arg1)
end

function MatRestoreRowUpperTriangular(arg1::Mat)
    ccall((:MatRestoreRowUpperTriangular,petsc),PetscErrorCode,(Mat,),arg1)
end

function MatGetColumn(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{Ptr{Int32}},arg5::Ptr{Ptr{Float64}})
    ccall((:MatGetColumn,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{Ptr{Int32}},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4,arg5)
end

function MatRestoreColumn(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{Ptr{Int32}},arg5::Ptr{Ptr{Float64}})
    ccall((:MatRestoreColumn,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{Ptr{Int32}},Ptr{Ptr{Float64}}),arg1,arg2,arg3,arg4,arg5)
end

function MatGetColumnVector(arg1::Mat,arg2::Vec,arg3::Int32)
    ccall((:MatGetColumnVector,petsc),PetscErrorCode,(Mat,Vec,PetscInt),arg1,arg2,arg3)
end

function MatSeqAIJGetArray(arg1::Mat,arg2::Ptr{Ptr{Float64}})
    ccall((:MatSeqAIJGetArray,petsc),PetscErrorCode,(Mat,Ptr{Ptr{Float64}}),arg1,arg2)
end

function MatSeqAIJRestoreArray(arg1::Mat,arg2::Ptr{Ptr{Float64}})
    ccall((:MatSeqAIJRestoreArray,petsc),PetscErrorCode,(Mat,Ptr{Ptr{Float64}}),arg1,arg2)
end

function MatSeqAIJGetMaxRowNonzeros(arg1::Mat,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatSeqAIJGetMaxRowNonzeros,petsc),PetscErrorCode,(Mat,Ptr{PetscInt}),arg1,arg2)
end

function MatSeqAIJSetValuesLocalFast(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Int32,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Ptr{Float64},arg7::InsertMode)
    ccall((:MatSeqAIJSetValuesLocalFast,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatDenseGetArray(arg1::Mat,arg2::Ptr{Ptr{Float64}})
    ccall((:MatDenseGetArray,petsc),PetscErrorCode,(Mat,Ptr{Ptr{Float64}}),arg1,arg2)
end

function MatDenseRestoreArray(arg1::Mat,arg2::Ptr{Ptr{Float64}})
    ccall((:MatDenseRestoreArray,petsc),PetscErrorCode,(Mat,Ptr{Ptr{Float64}}),arg1,arg2)
end

function MatGetBlockSize(arg1::Mat,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatGetBlockSize,petsc),PetscErrorCode,(Mat,Ptr{PetscInt}),arg1,arg2)
end

function MatSetBlockSize(arg1::Mat,arg2::Int32)
    ccall((:MatSetBlockSize,petsc),PetscErrorCode,(Mat,PetscInt),arg1,arg2)
end

function MatGetBlockSizes(arg1::Mat,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatGetBlockSizes,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatSetBlockSizes(arg1::Mat,arg2::Int32,arg3::Int32)
    ccall((:MatSetBlockSizes,petsc),PetscErrorCode,(Mat,PetscInt,PetscInt),arg1,arg2,arg3)
end

function MatSetBlockSizesFromMats(arg1::Mat,arg2::Mat,arg3::Mat)
    ccall((:MatSetBlockSizesFromMats,petsc),PetscErrorCode,(Mat,Mat,Mat),arg1,arg2,arg3)
end

function MatSetNThreads(arg1::Mat,arg2::Int32)
    ccall((:MatSetNThreads,petsc),PetscErrorCode,(Mat,PetscInt),arg1,arg2)
end

function MatGetNThreads(arg1::Mat,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
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

function MatIsTranspose(arg1::Mat,arg2::Mat,Float64::Cint,arg3::Ptr{PetscBool})
    ccall((:MatIsTranspose,petsc),PetscErrorCode,(Mat,Mat,Cint,Ptr{PetscBool}),arg1,arg2,PetscReal,arg3)
end

function MatIsHermitianTranspose(arg1::Mat,arg2::Mat,Float64::Cint,arg3::Ptr{PetscBool})
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

function MatIsSymmetric(arg1::Mat,Float64::Cint,arg2::Ptr{PetscBool})
    ccall((:MatIsSymmetric,petsc),PetscErrorCode,(Mat,Cint,Ptr{PetscBool}),arg1,PetscReal,arg2)
end

function MatIsStructurallySymmetric(arg1::Mat,arg2::Ptr{PetscBool})
    ccall((:MatIsStructurallySymmetric,petsc),PetscErrorCode,(Mat,Ptr{PetscBool}),arg1,arg2)
end

function MatIsHermitian(arg1::Mat,Float64::Cint,arg2::Ptr{PetscBool})
    ccall((:MatIsHermitian,petsc),PetscErrorCode,(Mat,Cint,Ptr{PetscBool}),arg1,PetscReal,arg2)
end

function MatIsSymmetricKnown(arg1::Mat,arg2::Ptr{PetscBool},arg3::Ptr{PetscBool})
    ccall((:MatIsSymmetricKnown,petsc),PetscErrorCode,(Mat,Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3)
end

function MatIsHermitianKnown(arg1::Mat,arg2::Ptr{PetscBool},arg3::Ptr{PetscBool})
    ccall((:MatIsHermitianKnown,petsc),PetscErrorCode,(Mat,Ptr{PetscBool},Ptr{PetscBool}),arg1,arg2,arg3)
end

function MatMissingDiagonal(arg1::Mat,arg2::Ptr{PetscBool},arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatMissingDiagonal,petsc),PetscErrorCode,(Mat,Ptr{PetscBool},Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatLoad(arg1::Mat,arg2::PetscViewer)
    ccall((:MatLoad,petsc),PetscErrorCode,(Mat,PetscViewer),arg1,arg2)
end

function MatGetRowIJ(arg1::Mat,arg2::Int32,arg3::PetscBool,arg4::PetscBool,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Ptr{Ptr{Int32}},arg7::Ptr{Ptr{Int32}},arg8::Ptr{PetscBool})
    ccall((:MatGetRowIJ,petsc),PetscErrorCode,(Mat,PetscInt,PetscBool,PetscBool,Ptr{PetscInt},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatRestoreRowIJ(arg1::Mat,arg2::Int32,arg3::PetscBool,arg4::PetscBool,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Ptr{Ptr{Int32}},arg7::Ptr{Ptr{Int32}},arg8::Ptr{PetscBool})
    ccall((:MatRestoreRowIJ,petsc),PetscErrorCode,(Mat,PetscInt,PetscBool,PetscBool,Ptr{PetscInt},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatGetColumnIJ(arg1::Mat,arg2::Int32,arg3::PetscBool,arg4::PetscBool,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Ptr{Ptr{Int32}},arg7::Ptr{Ptr{Int32}},arg8::Ptr{PetscBool})
    ccall((:MatGetColumnIJ,petsc),PetscErrorCode,(Mat,PetscInt,PetscBool,PetscBool,Ptr{PetscInt},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatRestoreColumnIJ(arg1::Mat,arg2::Int32,arg3::PetscBool,arg4::PetscBool,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Ptr{Ptr{Int32}},arg7::Ptr{Ptr{Int32}},arg8::Ptr{PetscBool})
    ccall((:MatRestoreColumnIJ,petsc),PetscErrorCode,(Mat,PetscInt,PetscBool,PetscBool,Ptr{PetscInt},Ptr{Ptr{Int32}},Ptr{Ptr{Int32}},Ptr{PetscBool}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatGetInfo(arg1::Mat,arg2::MatInfoType,arg3::Ptr{MatInfo})
    ccall((:MatGetInfo,petsc),PetscErrorCode,(Mat,MatInfoType,Ptr{MatInfo}),arg1,arg2,arg3)
end

function MatGetDiagonal(arg1::Mat,arg2::Vec)
    ccall((:MatGetDiagonal,petsc),PetscErrorCode,(Mat,Vec),arg1,arg2)
end

function MatGetRowMax(arg1::Mat,arg2::Vec,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatGetRowMax,petsc),PetscErrorCode,(Mat,Vec,Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatGetRowMin(arg1::Mat,arg2::Vec,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatGetRowMin,petsc),PetscErrorCode,(Mat,Vec,Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatGetRowMaxAbs(arg1::Mat,arg2::Vec,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatGetRowMaxAbs,petsc),PetscErrorCode,(Mat,Vec,Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatGetRowMinAbs(arg1::Mat,arg2::Vec,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
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

function MatMultEqual(arg1::Mat,arg2::Mat,arg3::Int32,arg4::Ptr{PetscBool})
    ccall((:MatMultEqual,petsc),PetscErrorCode,(Mat,Mat,PetscInt,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function MatMultAddEqual(arg1::Mat,arg2::Mat,arg3::Int32,arg4::Ptr{PetscBool})
    ccall((:MatMultAddEqual,petsc),PetscErrorCode,(Mat,Mat,PetscInt,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function MatMultTransposeEqual(arg1::Mat,arg2::Mat,arg3::Int32,arg4::Ptr{PetscBool})
    ccall((:MatMultTransposeEqual,petsc),PetscErrorCode,(Mat,Mat,PetscInt,Ptr{PetscBool}),arg1,arg2,arg3,arg4)
end

function MatMultTransposeAddEqual(arg1::Mat,arg2::Mat,arg3::Int32,arg4::Ptr{PetscBool})
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

function MatZeroRows(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Float64,arg5::Vec,arg6::Vec)
    ccall((:MatZeroRows,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscScalar,Vec,Vec),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatZeroRowsIS(arg1::Mat,arg2::IS,arg3::Float64,arg4::Vec,arg5::Vec)
    ccall((:MatZeroRowsIS,petsc),PetscErrorCode,(Mat,IS,PetscScalar,Vec,Vec),arg1,arg2,arg3,arg4,arg5)
end

function MatZeroRowsStencil(arg1::Mat,arg2::Int32,arg3::Ptr{MatStencil},arg4::Float64,arg5::Vec,arg6::Vec)
    ccall((:MatZeroRowsStencil,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{MatStencil},PetscScalar,Vec,Vec),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatZeroRowsColumnsStencil(arg1::Mat,arg2::Int32,arg3::Ptr{MatStencil},arg4::Float64,arg5::Vec,arg6::Vec)
    ccall((:MatZeroRowsColumnsStencil,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{MatStencil},PetscScalar,Vec,Vec),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatZeroRowsColumns(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Float64,arg5::Vec,arg6::Vec)
    ccall((:MatZeroRowsColumns,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscScalar,Vec,Vec),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatZeroRowsColumnsIS(arg1::Mat,arg2::IS,arg3::Float64,arg4::Vec,arg5::Vec)
    ccall((:MatZeroRowsColumnsIS,petsc),PetscErrorCode,(Mat,IS,PetscScalar,Vec,Vec),arg1,arg2,arg3,arg4,arg5)
end

function MatGetSize(arg1::Mat,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatGetSize,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatGetLocalSize(arg1::Mat,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatGetLocalSize,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatGetOwnershipRange(arg1::Mat,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatGetOwnershipRange,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatGetOwnershipRanges(arg1::Mat,arg2::Ptr{Ptr{Int32}})
    ccall((:MatGetOwnershipRanges,petsc),PetscErrorCode,(Mat,Ptr{Ptr{Int32}}),arg1,arg2)
end

function MatGetOwnershipRangeColumn(arg1::Mat,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatGetOwnershipRangeColumn,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatGetOwnershipRangesColumn(arg1::Mat,arg2::Ptr{Ptr{Int32}})
    ccall((:MatGetOwnershipRangesColumn,petsc),PetscErrorCode,(Mat,Ptr{Ptr{Int32}}),arg1,arg2)
end

function MatGetOwnershipIS(arg1::Mat,arg2::Ptr{IS},arg3::Ptr{IS})
    ccall((:MatGetOwnershipIS,petsc),PetscErrorCode,(Mat,Ptr{IS},Ptr{IS}),arg1,arg2,arg3)
end

function MatGetSubMatrices(arg1::Mat,arg2::Int32,arg3::Ptr{IS},arg4::Ptr{IS},arg5::MatReuse,arg6::Ptr{Ptr{Mat}})
    ccall((:MatGetSubMatrices,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{IS},Ptr{IS},MatReuse,Ptr{Ptr{Mat}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatGetSubMatricesMPI(arg1::Mat,arg2::Int32,arg3::Ptr{IS},arg4::Ptr{IS},arg5::MatReuse,arg6::Ptr{Ptr{Mat}})
    ccall((:MatGetSubMatricesMPI,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{IS},Ptr{IS},MatReuse,Ptr{Ptr{Mat}}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatDestroyMatrices(arg1::Int32,arg2::Ptr{Ptr{Mat}})
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

function MatCreateMPIAIJSumSeqAIJ(arg1::MPI_Comm,arg2::Mat,arg3::Int32,arg4::Int32,arg5::MatReuse,arg6::Ptr{Mat})
    ccall((:MatCreateMPIAIJSumSeqAIJ,petsc),PetscErrorCode,(MPI_Comm,Mat,PetscInt,PetscInt,MatReuse,Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatCreateMPIAIJSumSeqAIJSymbolic(arg1::MPI_Comm,arg2::Mat,arg3::Int32,arg4::Int32,arg5::Ptr{Mat})
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

function MatGetGhosts(arg1::Mat,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Ptr{Int32}})
    ccall((:MatGetGhosts,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{Ptr{Int32}}),arg1,arg2,arg3)
end

function MatIncreaseOverlap(arg1::Mat,arg2::Int32,arg3::Ptr{IS},arg4::Int32)
    ccall((:MatIncreaseOverlap,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{IS},PetscInt),arg1,arg2,arg3,arg4)
end

function MatMatMult(arg1::Mat,arg2::Mat,arg3::MatReuse,Float64::Cint,arg4::Ptr{Mat})
    ccall((:MatMatMult,petsc),PetscErrorCode,(Mat,Mat,MatReuse,Cint,Ptr{Mat}),arg1,arg2,arg3,PetscReal,arg4)
end

function MatMatMultSymbolic(arg1::Mat,arg2::Mat,Float64::Cint,arg3::Ptr{Mat})
    ccall((:MatMatMultSymbolic,petsc),PetscErrorCode,(Mat,Mat,Cint,Ptr{Mat}),arg1,arg2,PetscReal,arg3)
end

function MatMatMultNumeric(arg1::Mat,arg2::Mat,arg3::Mat)
    ccall((:MatMatMultNumeric,petsc),PetscErrorCode,(Mat,Mat,Mat),arg1,arg2,arg3)
end

function MatMatMatMult(arg1::Mat,arg2::Mat,arg3::Mat,arg4::MatReuse,Float64::Cint,arg5::Ptr{Mat})
    ccall((:MatMatMatMult,petsc),PetscErrorCode,(Mat,Mat,Mat,MatReuse,Cint,Ptr{Mat}),arg1,arg2,arg3,arg4,PetscReal,arg5)
end

function MatMatMatMultSymbolic(arg1::Mat,arg2::Mat,arg3::Mat,Float64::Cint,arg4::Ptr{Mat})
    ccall((:MatMatMatMultSymbolic,petsc),PetscErrorCode,(Mat,Mat,Mat,Cint,Ptr{Mat}),arg1,arg2,arg3,PetscReal,arg4)
end

function MatMatMatMultNumeric(arg1::Mat,arg2::Mat,arg3::Mat,arg4::Mat)
    ccall((:MatMatMatMultNumeric,petsc),PetscErrorCode,(Mat,Mat,Mat,Mat),arg1,arg2,arg3,arg4)
end

function MatPtAP(arg1::Mat,arg2::Mat,arg3::MatReuse,Float64::Cint,arg4::Ptr{Mat})
    ccall((:MatPtAP,petsc),PetscErrorCode,(Mat,Mat,MatReuse,Cint,Ptr{Mat}),arg1,arg2,arg3,PetscReal,arg4)
end

function MatPtAPSymbolic(arg1::Mat,arg2::Mat,Float64::Cint,arg3::Ptr{Mat})
    ccall((:MatPtAPSymbolic,petsc),PetscErrorCode,(Mat,Mat,Cint,Ptr{Mat}),arg1,arg2,PetscReal,arg3)
end

function MatPtAPNumeric(arg1::Mat,arg2::Mat,arg3::Mat)
    ccall((:MatPtAPNumeric,petsc),PetscErrorCode,(Mat,Mat,Mat),arg1,arg2,arg3)
end

function MatRARt(arg1::Mat,arg2::Mat,arg3::MatReuse,Float64::Cint,arg4::Ptr{Mat})
    ccall((:MatRARt,petsc),PetscErrorCode,(Mat,Mat,MatReuse,Cint,Ptr{Mat}),arg1,arg2,arg3,PetscReal,arg4)
end

function MatRARtSymbolic(arg1::Mat,arg2::Mat,Float64::Cint,arg3::Ptr{Mat})
    ccall((:MatRARtSymbolic,petsc),PetscErrorCode,(Mat,Mat,Cint,Ptr{Mat}),arg1,arg2,PetscReal,arg3)
end

function MatRARtNumeric(arg1::Mat,arg2::Mat,arg3::Mat)
    ccall((:MatRARtNumeric,petsc),PetscErrorCode,(Mat,Mat,Mat),arg1,arg2,arg3)
end

function MatTransposeMatMult(arg1::Mat,arg2::Mat,arg3::MatReuse,Float64::Cint,arg4::Ptr{Mat})
    ccall((:MatTransposeMatMult,petsc),PetscErrorCode,(Mat,Mat,MatReuse,Cint,Ptr{Mat}),arg1,arg2,arg3,PetscReal,arg4)
end

function MatTransposetMatMultSymbolic(arg1::Mat,arg2::Mat,Float64::Cint,arg3::Ptr{Mat})
    ccall((:MatTransposetMatMultSymbolic,petsc),PetscErrorCode,(Mat,Mat,Cint,Ptr{Mat}),arg1,arg2,PetscReal,arg3)
end

function MatTransposetMatMultNumeric(arg1::Mat,arg2::Mat,arg3::Mat)
    ccall((:MatTransposetMatMultNumeric,petsc),PetscErrorCode,(Mat,Mat,Mat),arg1,arg2,arg3)
end

function MatMatTransposeMult(arg1::Mat,arg2::Mat,arg3::MatReuse,Float64::Cint,arg4::Ptr{Mat})
    ccall((:MatMatTransposeMult,petsc),PetscErrorCode,(Mat,Mat,MatReuse,Cint,Ptr{Mat}),arg1,arg2,arg3,PetscReal,arg4)
end

function MatMatTransposeMultSymbolic(arg1::Mat,arg2::Mat,Float64::Cint,arg3::Ptr{Mat})
    ccall((:MatMatTransposeMultSymbolic,petsc),PetscErrorCode,(Mat,Mat,Cint,Ptr{Mat}),arg1,arg2,PetscReal,arg3)
end

function MatMatTransposeMultNumeric(arg1::Mat,arg2::Mat,arg3::Mat)
    ccall((:MatMatTransposeMultNumeric,petsc),PetscErrorCode,(Mat,Mat,Mat),arg1,arg2,arg3)
end

function MatAXPY(arg1::Mat,arg2::Float64,arg3::Mat,arg4::MatStructure)
    ccall((:MatAXPY,petsc),PetscErrorCode,(Mat,PetscScalar,Mat,MatStructure),arg1,arg2,arg3,arg4)
end

function MatAYPX(arg1::Mat,arg2::Float64,arg3::Mat,arg4::MatStructure)
    ccall((:MatAYPX,petsc),PetscErrorCode,(Mat,PetscScalar,Mat,MatStructure),arg1,arg2,arg3,arg4)
end

function MatScale(arg1::Mat,arg2::Float64)
    ccall((:MatScale,petsc),PetscErrorCode,(Mat,PetscScalar),arg1,arg2)
end

function MatShift(arg1::Mat,arg2::Float64)
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

function MatZeroRowsLocal(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Float64,arg5::Vec,arg6::Vec)
    ccall((:MatZeroRowsLocal,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscScalar,Vec,Vec),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatZeroRowsLocalIS(arg1::Mat,arg2::IS,arg3::Float64,arg4::Vec,arg5::Vec)
    ccall((:MatZeroRowsLocalIS,petsc),PetscErrorCode,(Mat,IS,PetscScalar,Vec,Vec),arg1,arg2,arg3,arg4,arg5)
end

function MatZeroRowsColumnsLocal(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Float64,arg5::Vec,arg6::Vec)
    ccall((:MatZeroRowsColumnsLocal,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscScalar,Vec,Vec),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatZeroRowsColumnsLocalIS(arg1::Mat,arg2::IS,arg3::Float64,arg4::Vec,arg5::Vec)
    ccall((:MatZeroRowsColumnsLocalIS,petsc),PetscErrorCode,(Mat,IS,PetscScalar,Vec,Vec),arg1,arg2,arg3,arg4,arg5)
end

function MatSetValuesLocal(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Int32,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Ptr{Float64},arg7::InsertMode)
    ccall((:MatSetValuesLocal,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatSetValuesBlockedLocal(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Int32,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Ptr{Float64},arg7::InsertMode)
    ccall((:MatSetValuesBlockedLocal,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{Float64},InsertMode),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatStashSetInitialSize(arg1::Mat,arg2::Int32,arg3::Int32)
    ccall((:MatStashSetInitialSize,petsc),PetscErrorCode,(Mat,PetscInt,PetscInt),arg1,arg2,arg3)
end

function MatStashGetInfo(arg1::Mat,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
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

function MatCreateMPIMatConcatenateSeqMat(arg1::MPI_Comm,arg2::Mat,arg3::Int32,arg4::MatReuse,arg5::Ptr{Mat})
    ccall((:MatCreateMPIMatConcatenateSeqMat,petsc),PetscErrorCode,(MPI_Comm,Mat,PetscInt,MatReuse,Ptr{Mat}),arg1,arg2,arg3,arg4,arg5)
end

function MatInodeAdjustForInodes(arg1::Mat,arg2::Ptr{IS},arg3::Ptr{IS})
    ccall((:MatInodeAdjustForInodes,petsc),PetscErrorCode,(Mat,Ptr{IS},Ptr{IS}),arg1,arg2,arg3)
end

function MatInodeGetInodeSizes(arg1::Mat,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Ptr{Int32}},arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatInodeGetInodeSizes,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{Ptr{Int32}},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function MatSeqAIJSetColumnIndices(arg1::Mat,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatSeqAIJSetColumnIndices,petsc),PetscErrorCode,(Mat,Ptr{PetscInt}),arg1,arg2)
end

function MatSeqBAIJSetColumnIndices(arg1::Mat,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatSeqBAIJSetColumnIndices,petsc),PetscErrorCode,(Mat,Ptr{PetscInt}),arg1,arg2)
end

function MatCreateSeqAIJWithArrays(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Ptr{Float64},arg7::Ptr{Mat})
    ccall((:MatCreateSeqAIJWithArrays,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatCreateSeqBAIJWithArrays(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Ptr{Float64},arg8::Ptr{Mat})
    ccall((:MatCreateSeqBAIJWithArrays,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatCreateSeqSBAIJWithArrays(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg7::Ptr{Float64},arg8::Ptr{Mat})
    ccall((:MatCreateSeqSBAIJWithArrays,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end

function MatCreateSeqAIJFromTriple(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Ptr{Float64},arg7::Ptr{Mat},arg8::Int32,arg9::PetscBool)
    ccall((:MatCreateSeqAIJFromTriple,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64},Ptr{Mat},PetscInt,PetscBool),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function MatSeqBAIJSetPreallocation(arg1::Mat,arg2::Int32,arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatSeqBAIJSetPreallocation,petsc),PetscErrorCode,(Mat,PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function MatSeqSBAIJSetPreallocation(arg1::Mat,arg2::Int32,arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatSeqSBAIJSetPreallocation,petsc),PetscErrorCode,(Mat,PetscInt,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function MatSeqAIJSetPreallocation(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatSeqAIJSetPreallocation,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatMPIBAIJSetPreallocation(arg1::Mat,arg2::Int32,arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Int32,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatMPIBAIJSetPreallocation,petsc),PetscErrorCode,(Mat,PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatMPISBAIJSetPreallocation(arg1::Mat,arg2::Int32,arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Int32,arg6::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatMPISBAIJSetPreallocation,petsc),PetscErrorCode,(Mat,PetscInt,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatMPIAIJSetPreallocation(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Int32,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatMPIAIJSetPreallocation,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end

function MatSeqAIJSetPreallocationCSR(arg1::Mat,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{Float64})
    ccall((:MatSeqAIJSetPreallocationCSR,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64}),arg1,arg2,arg3,arg4)
end

function MatSeqBAIJSetPreallocationCSR(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Ptr{Float64})
    ccall((:MatSeqBAIJSetPreallocationCSR,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5)
end

function MatMPIAIJSetPreallocationCSR(arg1::Mat,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{Float64})
    ccall((:MatMPIAIJSetPreallocationCSR,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64}),arg1,arg2,arg3,arg4)
end

function MatMPIBAIJSetPreallocationCSR(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Ptr{Float64})
    ccall((:MatMPIBAIJSetPreallocationCSR,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Float64}),arg1,arg2,arg3,arg4,arg5)
end

function MatMPIAdjSetPreallocation(arg1::Mat,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatMPIAdjSetPreallocation,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4)
end

function MatMPIDenseSetPreallocation(arg1::Mat,arg2::Ptr{Float64})
    ccall((:MatMPIDenseSetPreallocation,petsc),PetscErrorCode,(Mat,Ptr{Float64}),arg1,arg2)
end

function MatSeqDenseSetPreallocation(arg1::Mat,arg2::Ptr{Float64})
    ccall((:MatSeqDenseSetPreallocation,petsc),PetscErrorCode,(Mat,Ptr{Float64}),arg1,arg2)
end

function MatMPIAIJGetSeqAIJ(arg1::Mat,arg2::Ptr{Mat},arg3::Ptr{Mat},arg4::Ptr{Ptr{Int32}})
    ccall((:MatMPIAIJGetSeqAIJ,petsc),PetscErrorCode,(Mat,Ptr{Mat},Ptr{Mat},Ptr{Ptr{Int32}}),arg1,arg2,arg3,arg4)
end

function MatMPIBAIJGetSeqBAIJ(arg1::Mat,arg2::Ptr{Mat},arg3::Ptr{Mat},arg4::Ptr{Ptr{Int32}})
    ccall((:MatMPIBAIJGetSeqBAIJ,petsc),PetscErrorCode,(Mat,Ptr{Mat},Ptr{Mat},Ptr{Ptr{Int32}}),arg1,arg2,arg3,arg4)
end

function MatMPIAdjCreateNonemptySubcommMat(arg1::Mat,arg2::Ptr{Mat})
    ccall((:MatMPIAdjCreateNonemptySubcommMat,petsc),PetscErrorCode,(Mat,Ptr{Mat}),arg1,arg2)
end

function MatISSetPreallocation(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Int32,arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatISSetPreallocation,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end

function MatSeqDenseSetLDA(arg1::Mat,arg2::Int32)
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

function MatReorderForNonzeroDiagonal(arg1::Mat,Float64::Cint,arg2::IS,arg3::IS)
    ccall((:MatReorderForNonzeroDiagonal,petsc),PetscErrorCode,(Mat,Cint,IS,IS),arg1,PetscReal,arg2,arg3)
end

function MatCreateLaplacian(arg1::Mat,Float64::Cint,arg2::PetscBool,arg3::Ptr{Mat})
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

function MatGetInertia(arg1::Mat,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
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

function MatSOR(arg1::Mat,arg2::Vec,Float64::Cint,arg3::MatSORType,arg4::Cint,arg5::Int32,arg6::Int32,arg7::Vec)
    ccall((:MatSOR,petsc),PetscErrorCode,(Mat,Vec,Cint,MatSORType,Cint,PetscInt,PetscInt,Vec),arg1,arg2,PetscReal,arg3,arg4,arg5,arg6,arg7)
end

function MatColoringCreate(arg1::Mat,arg2::Ptr{MatColoring})
    ccall((:MatColoringCreate,petsc),PetscErrorCode,(Mat,Ptr{MatColoring}),arg1,arg2)
end

function MatColoringGetDegrees(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
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

function MatColoringSetDistance(arg1::MatColoring,arg2::Int32)
    ccall((:MatColoringSetDistance,petsc),PetscErrorCode,(MatColoring,PetscInt),arg1,arg2)
end

function MatColoringGetDistance(arg1::MatColoring,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatColoringGetDistance,petsc),PetscErrorCode,(MatColoring,Ptr{PetscInt}),arg1,arg2)
end

function MatColoringSetMaxColors(arg1::MatColoring,arg2::Int32)
    ccall((:MatColoringSetMaxColors,petsc),PetscErrorCode,(MatColoring,PetscInt),arg1,arg2)
end

function MatColoringGetMaxColors(arg1::MatColoring,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatColoringGetMaxColors,petsc),PetscErrorCode,(MatColoring,Ptr{PetscInt}),arg1,arg2)
end

function MatColoringApply(arg1::MatColoring,arg2::Ptr{ISColoring})
    ccall((:MatColoringApply,petsc),PetscErrorCode,(MatColoring,Ptr{ISColoring}),arg1,arg2)
end

function MatColoringRegister(arg1::Ptr{Uint8},arg2::Ptr{Void})
    ccall((:MatColoringRegister,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Void}),arg1,arg2)
end

function MatColoringPatch(arg1::Mat,arg2::Int32,arg3::Int32,ISColoringValue::Ptr{Cint},arg4::Ptr{ISColoring})
    ccall((:MatColoringPatch,petsc),PetscErrorCode,(Mat,PetscInt,PetscInt,Ptr{Cint},Ptr{ISColoring}),arg1,arg2,arg3,ISColoringValue,arg4)
end

function MatColoringSetWeightType(arg1::MatColoring,arg2::MatColoringWeightType)
    ccall((:MatColoringSetWeightType,petsc),PetscErrorCode,(MatColoring,MatColoringWeightType),arg1,arg2)
end

function MatColoringSetWeights(arg1::MatColoring,arg2::Ptr{Cint},arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatColoringSetWeights,petsc),PetscErrorCode,(MatColoring,Ptr{Cint},Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatColoringCreateWeights(arg1::MatColoring,arg2::Ptr{Ptr{Cint}},lperm::Ptr{Ptr{Int32}})
    ccall((:MatColoringCreateWeights,petsc),PetscErrorCode,(MatColoring,Ptr{Ptr{Cint}},Ptr{Ptr{Int32}}),arg1,arg2,lperm)
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

function MatFDColoringSetParameters(arg1::MatFDColoring,Float64::Cint,arg2::Cint)
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

function MatFDColoringGetPerturbedColumns(arg1::MatFDColoring,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Ptr{Int32}})
    ccall((:MatFDColoringGetPerturbedColumns,petsc),PetscErrorCode,(MatFDColoring,Ptr{PetscInt},Ptr{Ptr{Int32}}),arg1,arg2,arg3)
end

function MatFDColoringSetUp(arg1::Mat,arg2::ISColoring,arg3::MatFDColoring)
    ccall((:MatFDColoringSetUp,petsc),PetscErrorCode,(Mat,ISColoring,MatFDColoring),arg1,arg2,arg3)
end

function MatFDColoringSetBlockSize(arg1::MatFDColoring,arg2::Int32,arg3::Int32)
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

function MatPartitioningSetNParts(arg1::MatPartitioning,arg2::Int32)
    ccall((:MatPartitioningSetNParts,petsc),PetscErrorCode,(MatPartitioning,PetscInt),arg1,arg2)
end

function MatPartitioningSetAdjacency(arg1::MatPartitioning,arg2::Mat)
    ccall((:MatPartitioningSetAdjacency,petsc),PetscErrorCode,(MatPartitioning,Mat),arg1,arg2)
end

function MatPartitioningSetVertexWeights(arg1::MatPartitioning,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatPartitioningSetVertexWeights,petsc),PetscErrorCode,(MatPartitioning,Ptr{PetscInt}),arg1,arg2)
end

function MatPartitioningSetPartitionWeights(arg1::MatPartitioning,Float64::Ptr{Cint})
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

function MatPartitioningParmetisGetEdgeCut(arg1::MatPartitioning,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
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

function MatPartitioningChacoSetCoarseLevel(arg1::MatPartitioning,Float64::Cint)
    ccall((:MatPartitioningChacoSetCoarseLevel,petsc),PetscErrorCode,(MatPartitioning,Cint),arg1,PetscReal)
end

function MatPartitioningChacoSetEigenSolver(arg1::MatPartitioning,arg2::MPChacoEigenType)
    ccall((:MatPartitioningChacoSetEigenSolver,petsc),PetscErrorCode,(MatPartitioning,MPChacoEigenType),arg1,arg2)
end

function MatPartitioningChacoGetEigenSolver(arg1::MatPartitioning,arg2::Ptr{MPChacoEigenType})
    ccall((:MatPartitioningChacoGetEigenSolver,petsc),PetscErrorCode,(MatPartitioning,Ptr{MPChacoEigenType}),arg1,arg2)
end

function MatPartitioningChacoSetEigenTol(arg1::MatPartitioning,Float64::Cint)
    ccall((:MatPartitioningChacoSetEigenTol,petsc),PetscErrorCode,(MatPartitioning,Cint),arg1,PetscReal)
end

function MatPartitioningChacoGetEigenTol(arg1::MatPartitioning,arg2::Ptr{Cint})
    ccall((:MatPartitioningChacoGetEigenTol,petsc),PetscErrorCode,(MatPartitioning,Ptr{Cint}),arg1,arg2)
end

function MatPartitioningChacoSetEigenNumber(arg1::MatPartitioning,arg2::Int32)
    ccall((:MatPartitioningChacoSetEigenNumber,petsc),PetscErrorCode,(MatPartitioning,PetscInt),arg1,arg2)
end

function MatPartitioningChacoGetEigenNumber(arg1::MatPartitioning,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatPartitioningChacoGetEigenNumber,petsc),PetscErrorCode,(MatPartitioning,Ptr{PetscInt}),arg1,arg2)
end

function MatPartitioningPartySetGlobal(arg1::MatPartitioning,arg2::Ptr{Uint8})
    ccall((:MatPartitioningPartySetGlobal,petsc),PetscErrorCode,(MatPartitioning,Ptr{Uint8}),arg1,arg2)
end

function MatPartitioningPartySetLocal(arg1::MatPartitioning,arg2::Ptr{Uint8})
    ccall((:MatPartitioningPartySetLocal,petsc),PetscErrorCode,(MatPartitioning,Ptr{Uint8}),arg1,arg2)
end

function MatPartitioningPartySetCoarseLevel(arg1::MatPartitioning,Float64::Cint)
    ccall((:MatPartitioningPartySetCoarseLevel,petsc),PetscErrorCode,(MatPartitioning,Cint),arg1,PetscReal)
end

function MatPartitioningPartySetBipart(arg1::MatPartitioning,arg2::PetscBool)
    ccall((:MatPartitioningPartySetBipart,petsc),PetscErrorCode,(MatPartitioning,PetscBool),arg1,arg2)
end

function MatPartitioningPartySetMatchOptimization(arg1::MatPartitioning,arg2::PetscBool)
    ccall((:MatPartitioningPartySetMatchOptimization,petsc),PetscErrorCode,(MatPartitioning,PetscBool),arg1,arg2)
end

function MatPartitioningPTScotchSetImbalance(arg1::MatPartitioning,Float64::Cint)
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

function MatMeshToCellGraph(arg1::Mat,arg2::Int32,arg3::Ptr{Mat})
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

function MatMPIBAIJSetHashTableFactor(arg1::Mat,Float64::Cint)
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

function MatNullSpaceCreate(arg1::MPI_Comm,arg2::PetscBool,arg3::Int32,arg4::Ptr{Vec},arg5::Ptr{MatNullSpace})
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

function MatNullSpaceGetVecs(arg1::MatNullSpace,arg2::Ptr{PetscBool},arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{Ptr{Vec}})
    ccall((:MatNullSpaceGetVecs,petsc),PetscErrorCode,(MatNullSpace,Ptr{PetscBool},Ptr{PetscInt},Ptr{Ptr{Vec}}),arg1,arg2,arg3,arg4)
end

function MatNullSpaceCreateRigidBody(arg1::Vec,arg2::Ptr{MatNullSpace})
    ccall((:MatNullSpaceCreateRigidBody,petsc),PetscErrorCode,(Vec,Ptr{MatNullSpace}),arg1,arg2)
end

function MatReorderingSeqSBAIJ(arg1::Mat,arg2::IS)
    ccall((:MatReorderingSeqSBAIJ,petsc),PetscErrorCode,(Mat,IS),arg1,arg2)
end

function MatMPISBAIJSetHashTableFactor(arg1::Mat,Float64::Cint)
    ccall((:MatMPISBAIJSetHashTableFactor,petsc),PetscErrorCode,(Mat,Cint),arg1,PetscReal)
end

function MatSeqSBAIJSetColumnIndices(arg1::Mat,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatSeqSBAIJSetColumnIndices,petsc),PetscErrorCode,(Mat,Ptr{PetscInt}),arg1,arg2)
end

function MatSeqBAIJInvertBlockDiagonal(arg1::Mat)
    ccall((:MatSeqBAIJInvertBlockDiagonal,petsc),PetscErrorCode,(Mat,),arg1)
end

function MatCreateMAIJ(arg1::Mat,arg2::Int32,arg3::Ptr{Mat})
    ccall((:MatCreateMAIJ,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{Mat}),arg1,arg2,arg3)
end

function MatMAIJRedimension(arg1::Mat,arg2::Int32,arg3::Ptr{Mat})
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

function MatCreateMFFD(arg1::MPI_Comm,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Ptr{Mat})
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

function MatMFFDSetHHistory(arg1::Mat,arg2::Ptr{Float64},arg3::Int32)
    ccall((:MatMFFDSetHHistory,petsc),PetscErrorCode,(Mat,Ptr{Float64},PetscInt),arg1,arg2,arg3)
end

function MatMFFDResetHHistory(arg1::Mat)
    ccall((:MatMFFDResetHHistory,petsc),PetscErrorCode,(Mat,),arg1)
end

function MatMFFDSetFunctionError(arg1::Mat,Float64::Cint)
    ccall((:MatMFFDSetFunctionError,petsc),PetscErrorCode,(Mat,Cint),arg1,PetscReal)
end

function MatMFFDSetPeriod(arg1::Mat,arg2::Int32)
    ccall((:MatMFFDSetPeriod,petsc),PetscErrorCode,(Mat,PetscInt),arg1,arg2)
end

function MatMFFDGetH(arg1::Mat,arg2::Ptr{Float64})
    ccall((:MatMFFDGetH,petsc),PetscErrorCode,(Mat,Ptr{Float64}),arg1,arg2)
end

function MatMFFDSetOptionsPrefix(arg1::Mat,arg2::Ptr{Uint8})
    ccall((:MatMFFDSetOptionsPrefix,petsc),PetscErrorCode,(Mat,Ptr{Uint8}),arg1,arg2)
end

function MatMFFDCheckPositivity(arg1::Ptr{Void},arg2::Vec,arg3::Vec,arg4::Ptr{Float64})
    ccall((:MatMFFDCheckPositivity,petsc),PetscErrorCode,(Ptr{Void},Vec,Vec,Ptr{Float64}),arg1,arg2,arg3,arg4)
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

function MatMFFDDSSetUmin(arg1::Mat,Float64::Cint)
    ccall((:MatMFFDDSSetUmin,petsc),PetscErrorCode,(Mat,Cint),arg1,PetscReal)
end

function MatMFFDWPSetComputeNormU(arg1::Mat,arg2::PetscBool)
    ccall((:MatMFFDWPSetComputeNormU,petsc),PetscErrorCode,(Mat,PetscBool),arg1,arg2)
end

function PetscViewerMathematicaPutMatrix(arg1::PetscViewer,arg2::Int32,arg3::Int32,arg4::Ptr{Cint})
    ccall((:PetscViewerMathematicaPutMatrix,petsc),PetscErrorCode,(PetscViewer,PetscInt,PetscInt,Ptr{Cint}),arg1,arg2,arg3,arg4)
end

function PetscViewerMathematicaPutCSRMatrix(arg1::PetscViewer,arg2::Int32,arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg6::Ptr{Cint})
    ccall((:PetscViewerMathematicaPutCSRMatrix,petsc),PetscErrorCode,(PetscViewer,PetscInt,PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{Cint}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatCreateNest(arg1::MPI_Comm,arg2::Int32,arg3::Ptr{IS},arg4::Int32,arg5::Ptr{IS},arg6::Ptr{Mat},arg7::Ptr{Mat})
    ccall((:MatCreateNest,petsc),PetscErrorCode,(MPI_Comm,PetscInt,Ptr{IS},PetscInt,Ptr{IS},Ptr{Mat},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function MatNestGetSize(arg1::Mat,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatNestGetSize,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function MatNestGetISs(arg1::Mat,arg2::Ptr{IS},arg3::Ptr{IS})
    ccall((:MatNestGetISs,petsc),PetscErrorCode,(Mat,Ptr{IS},Ptr{IS}),arg1,arg2,arg3)
end

function MatNestGetLocalISs(arg1::Mat,arg2::Ptr{IS},arg3::Ptr{IS})
    ccall((:MatNestGetLocalISs,petsc),PetscErrorCode,(Mat,Ptr{IS},Ptr{IS}),arg1,arg2,arg3)
end

function MatNestGetSubMats(arg1::Mat,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{Ptr{Ptr{Mat}}})
    ccall((:MatNestGetSubMats,petsc),PetscErrorCode,(Mat,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{Ptr{Mat}}}),arg1,arg2,arg3,arg4)
end

function MatNestGetSubMat(arg1::Mat,arg2::Int32,arg3::Int32,arg4::Ptr{Mat})
    ccall((:MatNestGetSubMat,petsc),PetscErrorCode,(Mat,PetscInt,PetscInt,Ptr{Mat}),arg1,arg2,arg3,arg4)
end

function MatNestSetVecType(arg1::Mat,arg2::VecType)
    ccall((:MatNestSetVecType,petsc),PetscErrorCode,(Mat,VecType),arg1,arg2)
end

function MatNestSetSubMats(arg1::Mat,arg2::Int32,arg3::Ptr{IS},arg4::Int32,arg5::Ptr{IS},arg6::Ptr{Mat})
    ccall((:MatNestSetSubMats,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{IS},PetscInt,Ptr{IS},Ptr{Mat}),arg1,arg2,arg3,arg4,arg5,arg6)
end

function MatNestSetSubMat(arg1::Mat,arg2::Int32,arg3::Int32,arg4::Mat)
    ccall((:MatNestSetSubMat,petsc),PetscErrorCode,(Mat,PetscInt,PetscInt,Mat),arg1,arg2,arg3,arg4)
end

function MatChop(arg1::Mat,Float64::Cint)
    ccall((:MatChop,petsc),PetscErrorCode,(Mat,Cint),arg1,PetscReal)
end

function MatComputeBandwidth(arg1::Mat,Float64::Cint,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:MatComputeBandwidth,petsc),PetscErrorCode,(Mat,Cint,Ptr{PetscInt}),arg1,PetscReal,arg2)
end

function MatSubdomainsCreateCoalesce(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{Ptr{IS}})
    ccall((:MatSubdomainsCreateCoalesce,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{Ptr{IS}}),arg1,arg2,arg3,arg4)
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

function PCGetSetUpFailedReason(arg1::PC,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
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

function PCApplyRichardson(arg1::PC,arg2::Vec,arg3::Vec,arg4::Vec,Float64::Cint,arg5::Cint,arg6::Cint,arg7::Int32,arg8::PetscBool,arg9::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg10::Ptr{PCRichardsonConvergedReason})
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

function PCModifySubMatrices(arg1::PC,arg2::Int32,arg3::Ptr{IS},arg4::Ptr{IS},arg5::Ptr{Mat},arg6::Ptr{Void})
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

function PCSORSetOmega(arg1::PC,Float64::Cint)
    ccall((:PCSORSetOmega,petsc),PetscErrorCode,(PC,Cint),arg1,PetscReal)
end

function PCSORGetOmega(arg1::PC,arg2::Ptr{Cint})
    ccall((:PCSORGetOmega,petsc),PetscErrorCode,(PC,Ptr{Cint}),arg1,arg2)
end

function PCSORSetIterations(arg1::PC,arg2::Int32,arg3::Int32)
    ccall((:PCSORSetIterations,petsc),PetscErrorCode,(PC,PetscInt,PetscInt),arg1,arg2,arg3)
end

function PCSORGetIterations(arg1::PC,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PCSORGetIterations,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3)
end

function PCEisenstatSetOmega(arg1::PC,Float64::Cint)
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

function PCBJacobiSetTotalBlocks(arg1::PC,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PCBJacobiSetTotalBlocks,petsc),PetscErrorCode,(PC,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function PCBJacobiGetTotalBlocks(arg1::PC,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Ptr{Int32}})
    ccall((:PCBJacobiGetTotalBlocks,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{Ptr{Int32}}),arg1,arg2,arg3)
end

function PCBJacobiSetLocalBlocks(arg1::PC,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PCBJacobiSetLocalBlocks,petsc),PetscErrorCode,(PC,PetscInt,Ptr{PetscInt}),arg1,arg2,arg3)
end

function PCBJacobiGetLocalBlocks(arg1::PC,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Ptr{Int32}})
    ccall((:PCBJacobiGetLocalBlocks,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{Ptr{Int32}}),arg1,arg2,arg3)
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

function PCFactorSetZeroPivot(arg1::PC,Float64::Cint)
    ccall((:PCFactorSetZeroPivot,petsc),PetscErrorCode,(PC,Cint),arg1,PetscReal)
end

function PCFactorSetShiftType(arg1::PC,arg2::MatFactorShiftType)
    ccall((:PCFactorSetShiftType,petsc),PetscErrorCode,(PC,MatFactorShiftType),arg1,arg2)
end

function PCFactorSetShiftAmount(arg1::PC,Float64::Cint)
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

function PCFactorSetFill(arg1::PC,Float64::Cint)
    ccall((:PCFactorSetFill,petsc),PetscErrorCode,(PC,Cint),arg1,PetscReal)
end

function PCFactorSetColumnPivot(arg1::PC,Float64::Cint)
    ccall((:PCFactorSetColumnPivot,petsc),PetscErrorCode,(PC,Cint),arg1,PetscReal)
end

function PCFactorReorderForNonzeroDiagonal(arg1::PC,Float64::Cint)
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

function PCFactorSetLevels(arg1::PC,arg2::Int32)
    ccall((:PCFactorSetLevels,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCFactorGetLevels(arg1::PC,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PCFactorGetLevels,petsc),PetscErrorCode,(PC,Ptr{PetscInt}),arg1,arg2)
end

function PCFactorSetDropTolerance(arg1::PC,Float64::Cint,arg2::Cint,arg3::Int32)
    ccall((:PCFactorSetDropTolerance,petsc),PetscErrorCode,(PC,Cint,Cint,PetscInt),arg1,PetscReal,arg2,arg3)
end

function PCASMSetLocalSubdomains(arg1::PC,arg2::Int32,arg3::Ptr{IS},arg4::Ptr{IS})
    ccall((:PCASMSetLocalSubdomains,petsc),PetscErrorCode,(PC,PetscInt,Ptr{IS},Ptr{IS}),arg1,arg2,arg3,arg4)
end

function PCASMSetTotalSubdomains(arg1::PC,arg2::Int32,arg3::Ptr{IS},arg4::Ptr{IS})
    ccall((:PCASMSetTotalSubdomains,petsc),PetscErrorCode,(PC,PetscInt,Ptr{IS},Ptr{IS}),arg1,arg2,arg3,arg4)
end

function PCASMSetOverlap(arg1::PC,arg2::Int32)
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

function PCASMCreateSubdomains(arg1::Mat,arg2::Int32,arg3::Ptr{Ptr{IS}})
    ccall((:PCASMCreateSubdomains,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{Ptr{IS}}),arg1,arg2,arg3)
end

function PCASMDestroySubdomains(arg1::Int32,arg2::Ptr{IS},arg3::Ptr{IS})
    ccall((:PCASMDestroySubdomains,petsc),PetscErrorCode,(PetscInt,Ptr{IS},Ptr{IS}),arg1,arg2,arg3)
end

function PCASMCreateSubdomains2D(arg1::Int32,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Int32,arg7::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg8::Ptr{Ptr{IS}},arg9::Ptr{Ptr{IS}})
    ccall((:PCASMCreateSubdomains2D,petsc),PetscErrorCode,(PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Ptr{IS}},Ptr{Ptr{IS}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
end

function PCASMGetLocalSubdomains(arg1::PC,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Ptr{IS}},arg4::Ptr{Ptr{IS}})
    ccall((:PCASMGetLocalSubdomains,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{Ptr{IS}},Ptr{Ptr{IS}}),arg1,arg2,arg3,arg4)
end

function PCASMGetLocalSubmatrices(arg1::PC,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Ptr{Mat}})
    ccall((:PCASMGetLocalSubmatrices,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{Ptr{Mat}}),arg1,arg2,arg3)
end

function PCGASMSetTotalSubdomains(arg1::PC,arg2::Int32)
    ccall((:PCGASMSetTotalSubdomains,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCGASMSetSubdomains(arg1::PC,arg2::Int32,arg3::Ptr{IS},arg4::Ptr{IS})
    ccall((:PCGASMSetSubdomains,petsc),PetscErrorCode,(PC,PetscInt,Ptr{IS},Ptr{IS}),arg1,arg2,arg3,arg4)
end

function PCGASMSetOverlap(arg1::PC,arg2::Int32)
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

function PCGASMCreateSubdomains(arg1::Mat,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{Ptr{IS}})
    ccall((:PCGASMCreateSubdomains,petsc),PetscErrorCode,(Mat,PetscInt,Ptr{PetscInt},Ptr{Ptr{IS}}),arg1,arg2,arg3,arg4)
end

function PCGASMDestroySubdomains(arg1::Int32,arg2::Ptr{Ptr{IS}},arg3::Ptr{Ptr{IS}})
    ccall((:PCGASMDestroySubdomains,petsc),PetscErrorCode,(PetscInt,Ptr{Ptr{IS}},Ptr{Ptr{IS}}),arg1,arg2,arg3)
end

function PCGASMCreateSubdomains2D(arg1::PC,arg2::Int32,arg3::Int32,arg4::Int32,arg5::Int32,arg6::Int32,arg7::Int32,arg8::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg9::Ptr{Ptr{IS}},arg10::Ptr{Ptr{IS}})
    ccall((:PCGASMCreateSubdomains2D,petsc),PetscErrorCode,(PC,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{PetscInt},Ptr{Ptr{IS}},Ptr{Ptr{IS}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
end

function PCGASMGetSubdomains(arg1::PC,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Ptr{IS}},arg4::Ptr{Ptr{IS}})
    ccall((:PCGASMGetSubdomains,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{Ptr{IS}},Ptr{Ptr{IS}}),arg1,arg2,arg3,arg4)
end

function PCGASMGetSubmatrices(arg1::PC,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Ptr{Mat}})
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

function PCCompositeGetPC(arg1::PC,arg2::Int32,arg3::Ptr{PC})
    ccall((:PCCompositeGetPC,petsc),PetscErrorCode,(PC,PetscInt,Ptr{PC}),arg1,arg2,arg3)
end

function PCCompositeSpecialSetAlpha(arg1::PC,arg2::Float64)
    ccall((:PCCompositeSpecialSetAlpha,petsc),PetscErrorCode,(PC,PetscScalar),arg1,arg2)
end

function PCRedundantSetNumber(arg1::PC,arg2::Int32)
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

function PCSPAISetNBSteps(arg1::PC,arg2::Int32)
    ccall((:PCSPAISetNBSteps,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCSPAISetMax(arg1::PC,arg2::Int32)
    ccall((:PCSPAISetMax,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCSPAISetMaxNew(arg1::PC,arg2::Int32)
    ccall((:PCSPAISetMaxNew,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCSPAISetBlockSize(arg1::PC,arg2::Int32)
    ccall((:PCSPAISetBlockSize,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCSPAISetCacheSize(arg1::PC,arg2::Int32)
    ccall((:PCSPAISetCacheSize,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCSPAISetVerbose(arg1::PC,arg2::Int32)
    ccall((:PCSPAISetVerbose,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCSPAISetSp(arg1::PC,arg2::Int32)
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

function PCBJacobiGetLocalBlocks(arg1::PC,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Ptr{Int32}})
    ccall((:PCBJacobiGetLocalBlocks,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{Ptr{Int32}}),arg1,arg2,arg3)
end

function PCBJacobiGetTotalBlocks(arg1::PC,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Ptr{Int32}})
    ccall((:PCBJacobiGetTotalBlocks,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{Ptr{Int32}}),arg1,arg2,arg3)
end

function PCFieldSplitSetFields(arg1::PC,arg2::Ptr{Uint8},arg3::Int32,arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PCFieldSplitSetFields,petsc),PetscErrorCode,(PC,Ptr{Uint8},PetscInt,Ptr{PetscInt},Ptr{PetscInt}),arg1,arg2,arg3,arg4,arg5)
end

function PCFieldSplitSetType(arg1::PC,arg2::PCCompositeType)
    ccall((:PCFieldSplitSetType,petsc),PetscErrorCode,(PC,PCCompositeType),arg1,arg2)
end

function PCFieldSplitGetType(arg1::PC,arg2::Ptr{PCCompositeType})
    ccall((:PCFieldSplitGetType,petsc),PetscErrorCode,(PC,Ptr{PCCompositeType}),arg1,arg2)
end

function PCFieldSplitSetBlockSize(arg1::PC,arg2::Int32)
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

function PCSetCoordinates(arg1::PC,arg2::Int32,arg3::Int32,arg4::Ptr{Cint})
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

function PCBiCGStabCUSPSetTolerance(arg1::PC,Float64::Cint)
    ccall((:PCBiCGStabCUSPSetTolerance,petsc),PetscErrorCode,(PC,Cint),arg1,PetscReal)
end

function PCBiCGStabCUSPSetIterations(arg1::PC,arg2::Int32)
    ccall((:PCBiCGStabCUSPSetIterations,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCBiCGStabCUSPSetUseVerboseMonitor(arg1::PC,arg2::PetscBool)
    ccall((:PCBiCGStabCUSPSetUseVerboseMonitor,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCAINVCUSPSetDropTolerance(arg1::PC,Float64::Cint)
    ccall((:PCAINVCUSPSetDropTolerance,petsc),PetscErrorCode,(PC,Cint),arg1,PetscReal)
end

function PCAINVCUSPUseScaling(arg1::PC,arg2::PetscBool)
    ccall((:PCAINVCUSPUseScaling,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCAINVCUSPSetNonzeros(arg1::PC,arg2::Int32)
    ccall((:PCAINVCUSPSetNonzeros,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCAINVCUSPSetLinParameter(arg1::PC,arg2::Int32)
    ccall((:PCAINVCUSPSetLinParameter,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCPARMSSetGlobal(arg1::PC,arg2::PCPARMSGlobalType)
    ccall((:PCPARMSSetGlobal,petsc),PetscErrorCode,(PC,PCPARMSGlobalType),arg1,arg2)
end

function PCPARMSSetLocal(arg1::PC,arg2::PCPARMSLocalType)
    ccall((:PCPARMSSetLocal,petsc),PetscErrorCode,(PC,PCPARMSLocalType),arg1,arg2)
end

function PCPARMSSetSolveTolerances(arg1::PC,Float64::Cint,arg2::Int32)
    ccall((:PCPARMSSetSolveTolerances,petsc),PetscErrorCode,(PC,Cint,PetscInt),arg1,PetscReal,arg2)
end

function PCPARMSSetSolveRestart(arg1::PC,arg2::Int32)
    ccall((:PCPARMSSetSolveRestart,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCPARMSSetNonsymPerm(arg1::PC,arg2::PetscBool)
    ccall((:PCPARMSSetNonsymPerm,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCPARMSSetFill(arg1::PC,arg2::Int32,arg3::Int32,arg4::Int32)
    ccall((:PCPARMSSetFill,petsc),PetscErrorCode,(PC,PetscInt,PetscInt,PetscInt),arg1,arg2,arg3,arg4)
end

function PCGAMGSetType(arg1::PC,arg2::PCGAMGType)
    ccall((:PCGAMGSetType,petsc),PetscErrorCode,(PC,PCGAMGType),arg1,arg2)
end

function PCGAMGGetType(arg1::PC,arg2::Ptr{PCGAMGType})
    ccall((:PCGAMGGetType,petsc),PetscErrorCode,(PC,Ptr{PCGAMGType}),arg1,arg2)
end

function PCGAMGSetProcEqLim(arg1::PC,arg2::Int32)
    ccall((:PCGAMGSetProcEqLim,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCGAMGSetRepartitioning(arg1::PC,arg2::PetscBool)
    ccall((:PCGAMGSetRepartitioning,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCGAMGSetUseASMAggs(arg1::PC,arg2::PetscBool)
    ccall((:PCGAMGSetUseASMAggs,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCGAMGSetSolverType(arg1::PC,arg2::Ptr{Uint8},arg3::Int32)
    ccall((:PCGAMGSetSolverType,petsc),PetscErrorCode,(PC,Ptr{Uint8},PetscInt),arg1,arg2,arg3)
end

function PCGAMGSetThreshold(arg1::PC,Float64::Cint)
    ccall((:PCGAMGSetThreshold,petsc),PetscErrorCode,(PC,Cint),arg1,PetscReal)
end

function PCGAMGSetCoarseEqLim(arg1::PC,arg2::Int32)
    ccall((:PCGAMGSetCoarseEqLim,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCGAMGSetNlevels(arg1::PC,arg2::Int32)
    ccall((:PCGAMGSetNlevels,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCGAMGSetNSmooths(arg1::PC,arg2::Int32)
    ccall((:PCGAMGSetNSmooths,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCGAMGSetSymGraph(arg1::PC,arg2::PetscBool)
    ccall((:PCGAMGSetSymGraph,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCGAMGSetSquareGraph(arg1::PC,arg2::Int32)
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

function PCBDDCSetCoarseningRatio(arg1::PC,arg2::Int32)
    ccall((:PCBDDCSetCoarseningRatio,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCBDDCSetLevels(arg1::PC,arg2::Int32)
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

function PCBDDCSetDofsSplitting(arg1::PC,arg2::Int32,arg3::Ptr{IS})
    ccall((:PCBDDCSetDofsSplitting,petsc),PetscErrorCode,(PC,PetscInt,Ptr{IS}),arg1,arg2,arg3)
end

function PCBDDCSetDofsSplittingLocal(arg1::PC,arg2::Int32,arg3::Ptr{IS})
    ccall((:PCBDDCSetDofsSplittingLocal,petsc),PetscErrorCode,(PC,PetscInt,Ptr{IS}),arg1,arg2,arg3)
end

function PCBDDCSetLocalAdjacencyGraph(arg1::PC,arg2::Int32,arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),PetscCopyMode::Cint)
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

function PCISSetSubdomainScalingFactor(arg1::PC,arg2::Float64)
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

function PCMGSetLevels(arg1::PC,arg2::Int32,arg3::Ptr{MPI_Comm})
    ccall((:PCMGSetLevels,petsc),PetscErrorCode,(PC,PetscInt,Ptr{MPI_Comm}),arg1,arg2,arg3)
end

function PCMGGetLevels(arg1::PC,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:PCMGGetLevels,petsc),PetscErrorCode,(PC,Ptr{PetscInt}),arg1,arg2)
end

function PCMGSetNumberSmoothUp(arg1::PC,arg2::Int32)
    ccall((:PCMGSetNumberSmoothUp,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCMGSetNumberSmoothDown(arg1::PC,arg2::Int32)
    ccall((:PCMGSetNumberSmoothDown,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCMGSetCycleType(arg1::PC,arg2::PCMGCycleType)
    ccall((:PCMGSetCycleType,petsc),PetscErrorCode,(PC,PCMGCycleType),arg1,arg2)
end

function PCMGSetCycleTypeOnLevel(arg1::PC,arg2::Int32,arg3::PCMGCycleType)
    ccall((:PCMGSetCycleTypeOnLevel,petsc),PetscErrorCode,(PC,PetscInt,PCMGCycleType),arg1,arg2,arg3)
end

function PCMGSetCyclesOnLevel(arg1::PC,arg2::Int32,arg3::Int32)
    ccall((:PCMGSetCyclesOnLevel,petsc),PetscErrorCode,(PC,PetscInt,PetscInt),arg1,arg2,arg3)
end

function PCMGMultiplicativeSetCycles(arg1::PC,arg2::Int32)
    ccall((:PCMGMultiplicativeSetCycles,petsc),PetscErrorCode,(PC,PetscInt),arg1,arg2)
end

function PCMGSetGalerkin(arg1::PC,arg2::PetscBool)
    ccall((:PCMGSetGalerkin,petsc),PetscErrorCode,(PC,PetscBool),arg1,arg2)
end

function PCMGGetGalerkin(arg1::PC,arg2::Ptr{PetscBool})
    ccall((:PCMGGetGalerkin,petsc),PetscErrorCode,(PC,Ptr{PetscBool}),arg1,arg2)
end

function PCMGSetRhs(arg1::PC,arg2::Int32,arg3::Vec)
    ccall((:PCMGSetRhs,petsc),PetscErrorCode,(PC,PetscInt,Vec),arg1,arg2,arg3)
end

function PCMGSetX(arg1::PC,arg2::Int32,arg3::Vec)
    ccall((:PCMGSetX,petsc),PetscErrorCode,(PC,PetscInt,Vec),arg1,arg2,arg3)
end

function PCMGSetR(arg1::PC,arg2::Int32,arg3::Vec)
    ccall((:PCMGSetR,petsc),PetscErrorCode,(PC,PetscInt,Vec),arg1,arg2,arg3)
end

function PCMGSetRestriction(arg1::PC,arg2::Int32,arg3::Mat)
    ccall((:PCMGSetRestriction,petsc),PetscErrorCode,(PC,PetscInt,Mat),arg1,arg2,arg3)
end

function PCMGGetRestriction(arg1::PC,arg2::Int32,arg3::Ptr{Mat})
    ccall((:PCMGGetRestriction,petsc),PetscErrorCode,(PC,PetscInt,Ptr{Mat}),arg1,arg2,arg3)
end

function PCMGSetInterpolation(arg1::PC,arg2::Int32,arg3::Mat)
    ccall((:PCMGSetInterpolation,petsc),PetscErrorCode,(PC,PetscInt,Mat),arg1,arg2,arg3)
end

function PCMGGetInterpolation(arg1::PC,arg2::Int32,arg3::Ptr{Mat})
    ccall((:PCMGGetInterpolation,petsc),PetscErrorCode,(PC,PetscInt,Ptr{Mat}),arg1,arg2,arg3)
end

function PCMGSetRScale(arg1::PC,arg2::Int32,arg3::Vec)
    ccall((:PCMGSetRScale,petsc),PetscErrorCode,(PC,PetscInt,Vec),arg1,arg2,arg3)
end

function PCMGGetRScale(arg1::PC,arg2::Int32,arg3::Ptr{Vec})
    ccall((:PCMGGetRScale,petsc),PetscErrorCode,(PC,PetscInt,Ptr{Vec}),arg1,arg2,arg3)
end

function PCMGSetResidual(arg1::PC,arg2::Int32,arg3::Ptr{Void},arg4::Mat)
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

function KSPSetTolerances(arg1::KSP,Float64::Cint,arg2::Cint,arg3::Cint,arg4::Int32)
    ccall((:KSPSetTolerances,petsc),PetscErrorCode,(KSP,Cint,Cint,Cint,PetscInt),arg1,PetscReal,arg2,arg3,arg4)
end

function KSPGetTolerances(arg1::KSP,arg2::Ptr{Cint},arg3::Ptr{Cint},arg4::Ptr{Cint},arg5::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
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

function KSPGetIterationNumber(arg1::KSP,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:KSPGetIterationNumber,petsc),PetscErrorCode,(KSP,Ptr{PetscInt}),arg1,arg2)
end

function KSPGetTotalIterations(arg1::KSP,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:KSPGetTotalIterations,petsc),PetscErrorCode,(KSP,Ptr{PetscInt}),arg1,arg2)
end

function KSPCreateVecs(arg1::KSP,arg2::Int32,arg3::Ptr{Ptr{Vec}},arg4::Int32,arg5::Ptr{Ptr{Vec}})
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

function KSPMonitor(arg1::KSP,arg2::Int32,Float64::Cint)
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

function KSPGetResidualHistory(arg1::KSP,arg2::Ptr{Ptr{Cint}},arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:KSPGetResidualHistory,petsc),PetscErrorCode,(KSP,Ptr{Ptr{Cint}},Ptr{PetscInt}),arg1,arg2,arg3)
end

function KSPSetResidualHistory(arg1::KSP,Float64::Ptr{Cint},arg2::Int32,arg3::PetscBool)
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

function KSPSetWorkVecs(arg1::KSP,arg2::Int32)
    ccall((:KSPSetWorkVecs,petsc),PetscErrorCode,(KSP,PetscInt),arg1,arg2)
end

function PCKSPGetKSP(arg1::PC,arg2::Ptr{KSP})
    ccall((:PCKSPGetKSP,petsc),PetscErrorCode,(PC,Ptr{KSP}),arg1,arg2)
end

function PCBJacobiGetSubKSP(arg1::PC,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{Ptr{KSP}})
    ccall((:PCBJacobiGetSubKSP,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{KSP}}),arg1,arg2,arg3,arg4)
end

function PCASMGetSubKSP(arg1::PC,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{Ptr{KSP}})
    ccall((:PCASMGetSubKSP,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{KSP}}),arg1,arg2,arg3,arg4)
end

function PCGASMGetSubKSP(arg1::PC,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg4::Ptr{Ptr{KSP}})
    ccall((:PCGASMGetSubKSP,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{KSP}}),arg1,arg2,arg3,arg4)
end

function PCFieldSplitGetSubKSP(arg1::PC,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}),arg3::Ptr{Ptr{KSP}})
    ccall((:PCFieldSplitGetSubKSP,petsc),PetscErrorCode,(PC,Ptr{PetscInt},Ptr{Ptr{KSP}}),arg1,arg2,arg3)
end

function PCMGGetSmoother(arg1::PC,arg2::Int32,arg3::Ptr{KSP})
    ccall((:PCMGGetSmoother,petsc),PetscErrorCode,(PC,PetscInt,Ptr{KSP}),arg1,arg2,arg3)
end

function PCMGGetSmootherDown(arg1::PC,arg2::Int32,arg3::Ptr{KSP})
    ccall((:PCMGGetSmootherDown,petsc),PetscErrorCode,(PC,PetscInt,Ptr{KSP}),arg1,arg2,arg3)
end

function PCMGGetSmootherUp(arg1::PC,arg2::Int32,arg3::Ptr{KSP})
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

function KSPRichardsonSetScale(arg1::KSP,Float64::Cint)
    ccall((:KSPRichardsonSetScale,petsc),PetscErrorCode,(KSP,Cint),arg1,PetscReal)
end

function KSPRichardsonSetSelfScale(arg1::KSP,arg2::PetscBool)
    ccall((:KSPRichardsonSetSelfScale,petsc),PetscErrorCode,(KSP,PetscBool),arg1,arg2)
end

function KSPChebyshevSetEigenvalues(arg1::KSP,Float64::Cint,arg2::Cint)
    ccall((:KSPChebyshevSetEigenvalues,petsc),PetscErrorCode,(KSP,Cint,Cint),arg1,PetscReal,arg2)
end

function KSPChebyshevEstEigSet(arg1::KSP,Float64::Cint,arg2::Cint,arg3::Cint,arg4::Cint)
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

function KSPComputeEigenvalues(arg1::KSP,arg2::Int32,Float64::Ptr{Cint},arg3::Ptr{Cint},arg4::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:KSPComputeEigenvalues,petsc),PetscErrorCode,(KSP,PetscInt,Ptr{Cint},Ptr{Cint},Ptr{PetscInt}),arg1,arg2,PetscReal,arg3,arg4)
end

function KSPComputeEigenvaluesExplicitly(arg1::KSP,arg2::Int32,Float64::Ptr{Cint},arg3::Ptr{Cint})
    ccall((:KSPComputeEigenvaluesExplicitly,petsc),PetscErrorCode,(KSP,PetscInt,Ptr{Cint},Ptr{Cint}),arg1,arg2,PetscReal,arg3)
end

function KSPFCGSetMmax(arg1::KSP,arg2::Int32)
    ccall((:KSPFCGSetMmax,petsc),PetscErrorCode,(KSP,PetscInt),arg1,arg2)
end

function KSPFCGGetMmax(arg1::KSP,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:KSPFCGGetMmax,petsc),PetscErrorCode,(KSP,Ptr{PetscInt}),arg1,arg2)
end

function KSPFCGSetNprealloc(arg1::KSP,arg2::Int32)
    ccall((:KSPFCGSetNprealloc,petsc),PetscErrorCode,(KSP,PetscInt),arg1,arg2)
end

function KSPFCGGetNprealloc(arg1::KSP,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:KSPFCGGetNprealloc,petsc),PetscErrorCode,(KSP,Ptr{PetscInt}),arg1,arg2)
end

function KSPFCGSetTruncationType(arg1::KSP,arg2::KSPFCGTruncationType)
    ccall((:KSPFCGSetTruncationType,petsc),PetscErrorCode,(KSP,KSPFCGTruncationType),arg1,arg2)
end

function KSPFCGGetTruncationType(arg1::KSP,arg2::Ptr{KSPFCGTruncationType})
    ccall((:KSPFCGGetTruncationType,petsc),PetscErrorCode,(KSP,Ptr{KSPFCGTruncationType}),arg1,arg2)
end

function KSPGMRESSetRestart(arg1::KSP,arg2::Int32)
    ccall((:KSPGMRESSetRestart,petsc),PetscErrorCode,(KSP,PetscInt),arg1,arg2)
end

function KSPGMRESGetRestart(arg1::KSP,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
    ccall((:KSPGMRESGetRestart,petsc),PetscErrorCode,(KSP,Ptr{PetscInt}),arg1,arg2)
end

function KSPGMRESSetHapTol(arg1::KSP,Float64::Cint)
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

function KSPGMRESModifiedGramSchmidtOrthogonalization(arg1::KSP,arg2::Int32)
    ccall((:KSPGMRESModifiedGramSchmidtOrthogonalization,petsc),PetscErrorCode,(KSP,PetscInt),arg1,arg2)
end

function KSPGMRESClassicalGramSchmidtOrthogonalization(arg1::KSP,arg2::Int32)
    ccall((:KSPGMRESClassicalGramSchmidtOrthogonalization,petsc),PetscErrorCode,(KSP,PetscInt),arg1,arg2)
end

function KSPLGMRESSetAugDim(arg1::KSP,arg2::Int32)
    ccall((:KSPLGMRESSetAugDim,petsc),PetscErrorCode,(KSP,PetscInt),arg1,arg2)
end

function KSPLGMRESSetConstant(arg1::KSP)
    ccall((:KSPLGMRESSetConstant,petsc),PetscErrorCode,(KSP,),arg1)
end

function KSPGCRSetRestart(arg1::KSP,arg2::Int32)
    ccall((:KSPGCRSetRestart,petsc),PetscErrorCode,(KSP,PetscInt),arg1,arg2)
end

function KSPGCRGetRestart(arg1::KSP,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
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

function KSPFGMRESModifyPCNoChange(arg1::KSP,arg2::Int32,arg3::Int32,Float64::Cint,arg4::Ptr{Void})
    ccall((:KSPFGMRESModifyPCNoChange,petsc),PetscErrorCode,(KSP,PetscInt,PetscInt,Cint,Ptr{Void}),arg1,arg2,arg3,PetscReal,arg4)
end

function KSPFGMRESModifyPCKSP(arg1::KSP,arg2::Int32,arg3::Int32,Float64::Cint,arg4::Ptr{Void})
    ccall((:KSPFGMRESModifyPCKSP,petsc),PetscErrorCode,(KSP,PetscInt,PetscInt,Cint,Ptr{Void}),arg1,arg2,arg3,PetscReal,arg4)
end

function KSPFGMRESSetModifyPC(arg1::KSP,arg2::Ptr{Void},arg3::Ptr{Void},arg4::Ptr{Void})
    ccall((:KSPFGMRESSetModifyPC,petsc),PetscErrorCode,(KSP,Ptr{Void},Ptr{Void},Ptr{Void}),arg1,arg2,arg3,arg4)
end

function KSPQCGSetTrustRegionRadius(arg1::KSP,Float64::Cint)
    ccall((:KSPQCGSetTrustRegionRadius,petsc),PetscErrorCode,(KSP,Cint),arg1,PetscReal)
end

function KSPQCGGetQuadratic(arg1::KSP,arg2::Ptr{Cint})
    ccall((:KSPQCGGetQuadratic,petsc),PetscErrorCode,(KSP,Ptr{Cint}),arg1,arg2)
end

function KSPQCGGetTrialStepNorm(arg1::KSP,arg2::Ptr{Cint})
    ccall((:KSPQCGGetTrialStepNorm,petsc),PetscErrorCode,(KSP,Ptr{Cint}),arg1,arg2)
end

function KSPBCGSLSetXRes(arg1::KSP,Float64::Cint)
    ccall((:KSPBCGSLSetXRes,petsc),PetscErrorCode,(KSP,Cint),arg1,PetscReal)
end

function KSPBCGSLSetPol(arg1::KSP,arg2::PetscBool)
    ccall((:KSPBCGSLSetPol,petsc),PetscErrorCode,(KSP,PetscBool),arg1,arg2)
end

function KSPBCGSLSetEll(arg1::KSP,arg2::Int32)
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

function KSPMonitorSingularValue(arg1::KSP,arg2::Int32,Float64::Cint,arg3::Ptr{Void})
    ccall((:KSPMonitorSingularValue,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorDefault(arg1::KSP,arg2::Int32,Float64::Cint,arg3::Ptr{Void})
    ccall((:KSPMonitorDefault,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPLSQRMonitorDefault(arg1::KSP,arg2::Int32,Float64::Cint,arg3::Ptr{Void})
    ccall((:KSPLSQRMonitorDefault,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorRange(arg1::KSP,arg2::Int32,Float64::Cint,arg3::Ptr{Void})
    ccall((:KSPMonitorRange,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorDynamicTolerance(ksp::KSP,its::Int32,fnorm::Cint,dummy::Ptr{Void})
    ccall((:KSPMonitorDynamicTolerance,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),ksp,its,fnorm,dummy)
end

function KSPMonitorDynamicToleranceDestroy(dummy::Ptr{Ptr{Void}})
    ccall((:KSPMonitorDynamicToleranceDestroy,petsc),PetscErrorCode,(Ptr{Ptr{Void}},),dummy)
end

function KSPMonitorTrueResidualNorm(arg1::KSP,arg2::Int32,Float64::Cint,arg3::Ptr{Void})
    ccall((:KSPMonitorTrueResidualNorm,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorTrueResidualMaxNorm(arg1::KSP,arg2::Int32,Float64::Cint,arg3::Ptr{Void})
    ccall((:KSPMonitorTrueResidualMaxNorm,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorDefaultShort(arg1::KSP,arg2::Int32,Float64::Cint,arg3::Ptr{Void})
    ccall((:KSPMonitorDefaultShort,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorSolution(arg1::KSP,arg2::Int32,Float64::Cint,arg3::Ptr{Void})
    ccall((:KSPMonitorSolution,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorSAWs(arg1::KSP,arg2::Int32,Float64::Cint,arg3::Ptr{Void})
    ccall((:KSPMonitorSAWs,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorSAWsCreate(arg1::KSP,arg2::Ptr{Ptr{Void}})
    ccall((:KSPMonitorSAWsCreate,petsc),PetscErrorCode,(KSP,Ptr{Ptr{Void}}),arg1,arg2)
end

function KSPMonitorSAWsDestroy(arg1::Ptr{Ptr{Void}})
    ccall((:KSPMonitorSAWsDestroy,petsc),PetscErrorCode,(Ptr{Ptr{Void}},),arg1)
end

function KSPGMRESMonitorKrylov(arg1::KSP,arg2::Int32,Float64::Cint,arg3::Ptr{Void})
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

function KSPSetTabLevel(arg1::KSP,arg2::Int32)
    ccall((:KSPSetTabLevel,petsc),PetscErrorCode,(KSP,PetscInt),arg1,arg2)
end

function KSPGetTabLevel(arg1::KSP,arg2::Union(Ptr{Int32},StridedArray{Int32},Ptr{Void}))
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

function KSPSetSupportedNorm(ksp::KSP,arg1::KSPNormType,arg2::PCSide,arg3::Int32)
    ccall((:KSPSetSupportedNorm,petsc),PetscErrorCode,(KSP,KSPNormType,PCSide,PetscInt),ksp,arg1,arg2,arg3)
end

function KSPSetCheckNormIteration(arg1::KSP,arg2::Int32)
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

function KSPConvergedDefault(arg1::KSP,arg2::Int32,Float64::Cint,arg3::Ptr{KSPConvergedReason},arg4::Ptr{Void})
    ccall((:KSPConvergedDefault,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{KSPConvergedReason},Ptr{Void}),arg1,arg2,PetscReal,arg3,arg4)
end

function KSPConvergedLSQR(arg1::KSP,arg2::Int32,Float64::Cint,arg3::Ptr{KSPConvergedReason},arg4::Ptr{Void})
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

function KSPConvergedSkip(arg1::KSP,arg2::Int32,Float64::Cint,arg3::Ptr{KSPConvergedReason},arg4::Ptr{Void})
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

function KSPNASHSetRadius(arg1::KSP,Float64::Cint)
    ccall((:KSPNASHSetRadius,petsc),PetscErrorCode,(KSP,Cint),arg1,PetscReal)
end

function KSPNASHGetNormD(arg1::KSP,arg2::Ptr{Cint})
    ccall((:KSPNASHGetNormD,petsc),PetscErrorCode,(KSP,Ptr{Cint}),arg1,arg2)
end

function KSPNASHGetObjFcn(arg1::KSP,arg2::Ptr{Cint})
    ccall((:KSPNASHGetObjFcn,petsc),PetscErrorCode,(KSP,Ptr{Cint}),arg1,arg2)
end

function KSPSTCGSetRadius(arg1::KSP,Float64::Cint)
    ccall((:KSPSTCGSetRadius,petsc),PetscErrorCode,(KSP,Cint),arg1,PetscReal)
end

function KSPSTCGGetNormD(arg1::KSP,arg2::Ptr{Cint})
    ccall((:KSPSTCGGetNormD,petsc),PetscErrorCode,(KSP,Ptr{Cint}),arg1,arg2)
end

function KSPSTCGGetObjFcn(arg1::KSP,arg2::Ptr{Cint})
    ccall((:KSPSTCGGetObjFcn,petsc),PetscErrorCode,(KSP,Ptr{Cint}),arg1,arg2)
end

function KSPGLTRSetRadius(arg1::KSP,Float64::Cint)
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

function KSPMonitorLGResidualNorm(arg1::KSP,arg2::Int32,Float64::Cint,arg3::Ptr{PetscObject})
    ccall((:KSPMonitorLGResidualNorm,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{PetscObject}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorLGResidualNormDestroy(arg1::Ptr{Ptr{PetscObject}})
    ccall((:KSPMonitorLGResidualNormDestroy,petsc),PetscErrorCode,(Ptr{Ptr{PetscObject}},),arg1)
end

function KSPMonitorLGTrueResidualNormCreate(arg1::Ptr{Uint8},arg2::Ptr{Uint8},arg3::Cint,arg4::Cint,arg5::Cint,arg6::Cint,arg7::Ptr{Ptr{PetscObject}})
    ccall((:KSPMonitorLGTrueResidualNormCreate,petsc),PetscErrorCode,(Ptr{Uint8},Ptr{Uint8},Cint,Cint,Cint,Cint,Ptr{Ptr{PetscObject}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
end

function KSPMonitorLGTrueResidualNorm(arg1::KSP,arg2::Int32,Float64::Cint,arg3::Ptr{PetscObject})
    ccall((:KSPMonitorLGTrueResidualNorm,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{PetscObject}),arg1,arg2,PetscReal,arg3)
end

function KSPMonitorLGTrueResidualNormDestroy(arg1::Ptr{Ptr{PetscObject}})
    ccall((:KSPMonitorLGTrueResidualNormDestroy,petsc),PetscErrorCode,(Ptr{Ptr{PetscObject}},),arg1)
end

function KSPMonitorLGRange(arg1::KSP,arg2::Int32,Float64::Cint,arg3::Ptr{Void})
    ccall((:KSPMonitorLGRange,petsc),PetscErrorCode,(KSP,PetscInt,Cint,Ptr{Void}),arg1,arg2,PetscReal,arg3)
end

function PCShellSetPreSolve(arg1::PC,arg2::Ptr{Void})
    ccall((:PCShellSetPreSolve,petsc),PetscErrorCode,(PC,Ptr{Void}),arg1,arg2)
end

function PCShellSetPostSolve(arg1::PC,arg2::Ptr{Void})
    ccall((:PCShellSetPostSolve,petsc),PetscErrorCode,(PC,Ptr{Void}),arg1,arg2)
end

function KSPFischerGuessCreate(arg1::KSP,arg2::Int32,arg3::Int32,arg4::Ptr{KSPFischerGuess})
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

function KSPSetUseFischerGuess(arg1::KSP,arg2::Int32,arg3::Int32)
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
