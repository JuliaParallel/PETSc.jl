# PetscBool is a C `bool` (one byte) since PETSc 3.24: a wider Julia type read bytes PETSc
# never wrote, so false results could come back as true (JuliaParallel/PETSc.jl#268).
using Test
using PETSc, MPI
using PETSc: LibPETSc

MPI.Initialized() || MPI.Init()

@testset "PetscBool is one byte" begin
    @test sizeof(LibPETSc.PetscBool) == 1
    @test Bool(LibPETSc.PETSC_FALSE) === false && Bool(LibPETSc.PETSC_TRUE) === true
    @test LibPETSc.PetscBool(false) == LibPETSc.PETSC_FALSE

    petsclib = PETSc.getlib()
    PETSc.initialize(petsclib)
    # false answers must read as false, many times over (fresh Refs each call)
    opts = LibPETSc.PetscOptions{typeof(petsclib)}()
    for _ in 1:200
        @test LibPETSc.PetscOptionsHasName(petsclib, opts, "", "-no_such_option") == LibPETSc.PETSC_FALSE
    end
    is = LibPETSc.ISCreateGeneral(petsclib, MPI.COMM_SELF, 3, petsclib.PetscInt[3, 1, 2],
                                  LibPETSc.PETSC_COPY_VALUES)
    @test LibPETSc.ISSorted(petsclib, is) == LibPETSc.PETSC_FALSE
    LibPETSc.ISSort(petsclib, is)
    @test LibPETSc.ISSorted(petsclib, is) == LibPETSc.PETSC_TRUE
    LibPETSc.ISDestroy(petsclib, is)
    # arrays of PetscBool have stride one
    buf = fill(LibPETSc.PETSC_TRUE, 4)
    LibPETSc.PetscOptionsSetValue(petsclib, opts, "-bools", "0,1,0")
    n, set = LibPETSc.PetscOptionsGetBoolArray(petsclib, opts, "", "-bools", buf, length(buf))
    @test n == 3 && buf[1:3] == [LibPETSc.PETSC_FALSE, LibPETSc.PETSC_TRUE, LibPETSc.PETSC_FALSE]
    PETSc.finalize(petsclib)
end

# The RegisterAllCalled reset after PetscInitialize writes the flag and nothing past it. 
# A sentinel goes into the byte after each flag, and is restored afterwards. 
# The flags are not exported on Windows, so there is nothing to check there.
@testset "RegisterAllCalled reset writes one byte" begin
    petsclib = PETSc.getlib()
    PETSc.initialize(petsclib)
    handle, _ = PETSc.ensure_library_handle(petsclib)
    lib = PETSc.library_ptr(handle)
    for sym in (:TaoRegisterAllCalled, :TaoTermRegisterAllCalled, :TSTrajectoryRegisterAllCalled)
        p = PETSc.Libdl.dlsym_e(lib, sym)
        p == C_NULL && continue
        flag = Ptr{UInt8}(p)
        saved = unsafe_load(flag + 1)
        unsafe_store!(flag + 1, 0xa5)
        PETSc._reset_stale_register_flags(petsclib)
        @test unsafe_load(flag) == 0x00
        @test unsafe_load(flag + 1) == 0xa5
        unsafe_store!(flag + 1, saved)
    end
    PETSc.finalize(petsclib)
end
