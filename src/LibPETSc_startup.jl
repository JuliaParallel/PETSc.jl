# Functions needed to find libraries
#=
function getlibs()
    libs = ()
    petsc_libs = ENV["JULIA_PETSC_LIBRARY"]

    flags = Libdl.RTLD_LAZY | Libdl.RTLD_DEEPBIND | Libdl.RTLD_GLOBAL

    for petsc_lib in Base.parse_load_path(petsc_libs)
        libs = (libs..., (petsc_lib, flags))
    end
    return libs
end
=#
const libs = @static if !haskey(ENV, "JULIA_PETSC_LIBRARY")
    using PETSc_jll
    (
        ((PETSc_jll.libpetsc_Float64_Real_Int64,), Float64, Int64),
        ((PETSc_jll.libpetsc_Float32_Real_Int64,), Float32, Int64),
        ((PETSc_jll.libpetsc_Float64_Complex_Int64,), Complex{Float64}, Int64),
        ((PETSc_jll.libpetsc_Float32_Complex_Int64,), Complex{Float32}, Int64),
        ((PETSc_jll.libpetsc_Float64_Real_Int32,), Float64, Int32),
        ((PETSc_jll.libpetsc_Float32_Real_Int32,), Float32, Int32),
        ((PETSc_jll.libpetsc_Float64_Complex_Int32,), Complex{Float64}, Int32),
        ((PETSc_jll.libpetsc_Float32_Complex_Int32,), Complex{Float32}, Int32),
    )
else
    error("JULIA_PETSC_LIBRARY not currently working")
end

const petsc_library_file =
    #get(ENV, "JULIA_PETSC_LIBRARY_PATH", "../lib/petsc_library.jl")
    get(ENV, "JULIA_PETSC_LIBRARY_PATH", "autowrapped/petsc_library.jl")    # if all is well, we should be able to use this


#=    
function DataTypeFromString(libhdl::Ptr{Cvoid}, name::AbstractString)
    PetscDataTypeFromString_ptr = dlsym(libhdl, :PetscDataTypeFromString)
    dtype_ref = Ref{PetscDataType}()
    found_ref = Ref{PetscBool}()
    @chk ccall(
        PetscDataTypeFromString_ptr,
        PetscErrorCode,
        (Cstring, Ptr{PetscDataType}, Ptr{PetscBool}),
        name,
        dtype_ref,
        found_ref,
    )
    @assert found_ref[] == PETSC_TRUE
    return dtype_ref[]
end

function PetscDataTypeGetSize(libhdl::Ptr{Cvoid}, dtype::PetscDataType)
    PetscDataTypeGetSize_ptr = dlsym(libhdl, :PetscDataTypeGetSize)
    datasize_ref = Ref{Csize_t}()
    @chk ccall(
        PetscDataTypeGetSize_ptr,
        PetscErrorCode,
        (PetscDataType, Ptr{Csize_t}),
        dtype,
        datasize_ref,
    )
    return datasize_ref[]
end
=#