# Functions needed to find libraries
using Preferences

# @load_preference is evaluated at compile time; changing a preference via
# set_library!() triggers automatic recompilation (standard Preferences.jl behaviour).
const libs = @static if @load_preference("library", nothing) == "PETSc_CUDA_jll"
    using PETSc_CUDA_jll
    PETSc_CUDA_jll.is_available() || error(
        "PETSc_CUDA_jll has no build for this machine, but `PETSc.set_library!(:cuda)` selected it. " *
        "Run `PETSc.unset_library!()` and restart Julia to use PETSc_jll.")
    # PETSc_CUDA_jll is built for PetscScalar = Float64 only
    @static if @load_preference("PetscInt", "Int64") == "Int32"
        (((PETSc_CUDA_jll.libpetsc_Float64_Real_Int32,), Float64, Int32),)
    else
        (((PETSc_CUDA_jll.libpetsc_Float64_Real_Int64,), Float64, Int64),)
    end
elseif @load_preference("library_path", nothing) === nothing
    using PETSc_jll
    # Only ONE PetscInt width is registered per process (default Int64, switchable with the
    # "PetscInt" preference / `set_petscint!`). The Int64 and Int32 libpetsc variants link
    # external packages that export identical symbols with different integer ABIs (hypre,
    # SuperLU_DIST), so loading both widths into one flat ELF namespace cross-binds them.
    @static if @load_preference("PetscInt", "Int64") == "Int32"
        (
            ((PETSc_jll.libpetsc_Float64_Real_Int32,), Float64, Int32),
            ((PETSc_jll.libpetsc_Float32_Real_Int32,), Float32, Int32),
            ((PETSc_jll.libpetsc_Float64_Complex_Int32,), Complex{Float64}, Int32),
            ((PETSc_jll.libpetsc_Float32_Complex_Int32,), Complex{Float32}, Int32),
        )
    else
        (
            ((PETSc_jll.libpetsc_Float64_Real_Int64,), Float64, Int64),
            ((PETSc_jll.libpetsc_Float32_Real_Int64,), Float32, Int64),
            ((PETSc_jll.libpetsc_Float64_Complex_Int64,), Complex{Float64}, Int64),
            ((PETSc_jll.libpetsc_Float32_Complex_Int64,), Complex{Float32}, Int64),
        )
    end
else
    let _path   = @load_preference("library_path"),
        _scalar = @load_preference("PetscScalar", "Float64"),
        _int    = @load_preference("PetscInt",    "Int64")
        _st = _scalar == "Float32"        ? Float32          :
              _scalar == "ComplexFloat64" ? Complex{Float64} :
              _scalar == "ComplexFloat32" ? Complex{Float32} : Float64
        _it = _int == "Int32" ? Int32 : Int64
        (((_path,), _st, _it),)
    end
end

const petsc_library_file =
    get(ENV, "JULIA_PETSC_LIBRARY_PATH", "autowrapped/petsc_library.jl")    # if all is well, we should be able to use this
