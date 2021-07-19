module LibPETSc

using Libdl
using MPI
MPI.Initialized() || MPI.Init()

export PetscLibType,
    petsclibs,
    PetscBool,
    PETSC_TRUE,
    PETSC_FALSE,
    UnionPetscLibType,
    getlib,
    MatAssemblyType,
    MAT_FLUSH_ASSEMBLY,
    MAT_FINAL_ASSEMBLY,
    InsertMode,
    NOT_SET_VALUES,
    INSERT_VALUES,
    ADD_VALUES,
    MAX_VALUES,
    MIN_VALUES,
    INSERT_ALL_VALUES,
    ADD_ALL_VALUES,
    INSERT_BC_VALUES,
    ADD_BC_VALUES,
    NormType,
    NORM_1,
    NORM_2,
    NORM_FROBENIUS,
    NORM_INFINITY,
    NORM_1_AND_2

include("LibPETSc_const.jl")
include("LibPETSc_startup.jl")
include("LibPETSc_lib.jl")

include(petsc_library_file)

end # module
