
using Libdl
const libs = (PETSc_jll.libpetsc, )

function initialize(libhdl::Ptr{Cvoid})
  PetscInitializeNoArguments_ptr = dlsym(libhdl, :PetscInitializeNoArguments)
  @chk ccall(PetscInitializeNoArguments_ptr, PetscErrorCode, ())
end


function DataTypeFromString(libhdl::Ptr{Cvoid}, name::AbstractString)
    PetscDataTypeFromString_ptr = dlsym(libhdl, :PetscDataTypeFromString)
    dtype_ref = Ref{PetscDataType}()
    found_ref = Ref{PetscBool}()
    @chk ccall(PetscDataTypeFromString_ptr, PetscErrorCode,
             (Cstring, Ptr{PetscDataType}, Ptr{PetscBool}),
             name, dtype_ref, found_ref)
    @assert found_ref[] == PETSC_TRUE
    return dtype_ref[]
end
function PetscDataTypeGetSize(libhdl::Ptr{Cvoid}, dtype::PetscDataType)
    PetscDataTypeGetSize_ptr = dlsym(libhdl, :PetscDataTypeGetSize)
    datasize_ref = Ref{Csize_t}()
    @chk ccall(PetscDataTypeGetSize_ptr, PetscErrorCode,
             (PetscDataType, Ptr{Csize_t}), 
             dtype, datasize_ref)
    return datasize_ref[]
end

const libtypes = map(libs) do lib
    libhdl = dlopen(lib)
    initialize(libhdl)
    PETSC_REAL = DataTypeFromString(libhdl, "Real")
    PETSC_SCALAR = DataTypeFromString(libhdl, "Scalar")
    PETSC_INT_SIZE = PetscDataTypeGetSize(libhdl, PETSC_INT)

    PetscReal =
        PETSC_REAL == PETSC_DOUBLE ? Cdouble :
        PETSC_REAL == PETSC_FLOAT ? Cfloat :
        error("PETSC_REAL = $PETSC_REAL not supported.")

    PetscScalar =
        PETSC_SCALAR == PETSC_REAL ? PetscReal :
        PETSC_SCALAR == PETSC_COMPLEX ? Complex{PetscReal} :
        error("PETSC_SCALAR = $PETSC_SCALAR not supported.")
    
    PetscInt =
        PETSC_INT_SIZE == 4 ? Int32 :
        PETSC_INT_SIZE == 8 ? Int64 :
        error("PETSC_INT_SIZE = $PETSC_INT_SIZE not supported.")

    # TODO: PetscBLASInt, PetscMPIInt ?
    return (lib, PetscScalar, PetscReal, PetscInt)
end

const scalar_types = map(x -> x[2], libtypes)

macro for_libpetsc(expr)
  quote
    for (libpetsc, PetscScalar, PetscReal, PetscInt) in libtypes
      @eval esc($expr)
    end
  end
end

@for_libpetsc inttype(::Type{$PetscScalar}) = $PetscInt
