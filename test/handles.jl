# test/handles.jl
# destroy! on the handles with no high-level layer (IS, AO, PF, Tao), and the
# ownership of the index sets ISColoringGetIS hands out.

using Test
using PETSc
using MPI

MPI.Initialized() || MPI.Init()

@testset "IS, AO, PF and Tao handles" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    comm = MPI.COMM_SELF
    PetscInt = petsclib.PetscInt
    LibPETSc = PETSc.LibPETSc

    @testset "destroy!" begin
        is = LibPETSc.ISCreateStride(petsclib, comm, PetscInt(4), PetscInt(0), PetscInt(1))
        ao = LibPETSc.AOCreateBasic(petsclib, comm, PetscInt(2), PetscInt[1, 0], PetscInt[0, 1])
        pf = LibPETSc.PFCreate(petsclib, comm, PetscInt(1), PetscInt(1))
        tao = LibPETSc.TaoCreate(petsclib, comm)
        for h in (is, ao, pf, tao)
            @test PETSc.owns(h)
            @test PETSc.destroy!(h) === nothing
            @test h.ptr == C_NULL
            @test PETSc.destroy!(h) === nothing     # a second call is a no-op
        end
    end

    @testset "ISColoringGetIS: borrowed or owned by mode" begin
        da = PETSc.DMDA(petsclib, comm, (PETSc.DM_BOUNDARY_NONE,), (6,), 1, 1)
        coloring = LibPETSc.DMCreateColoring(petsclib, da, LibPETSc.IS_COLORING_GLOBAL)

        # PETSC_USE_POINTER: the coloring keeps the index sets
        nn, borrowed = LibPETSc.ISColoringGetIS(petsclib, coloring, LibPETSc.PETSC_USE_POINTER)
        @test nn == length(borrowed) > 0
        @test !any(PETSc.owns, borrowed)
        @test sum(is -> LibPETSc.ISGetSize(petsclib, is), borrowed) == 6
        PETSc.destroy!(borrowed[1])                 # does nothing
        @test borrowed[1].ptr != C_NULL
        LibPETSc.ISColoringRestoreIS(petsclib, coloring, LibPETSc.PETSC_USE_POINTER, borrowed)

        # PETSC_OWN_POINTER: the caller takes them, and holds the only reference
        nn, owned = LibPETSc.ISColoringGetIS(petsclib, coloring, LibPETSc.PETSC_OWN_POINTER)
        @test all(PETSc.owns, owned)
        @test all(is -> LibPETSc.PetscObjectGetReference(petsclib, is) == 1, owned)
        @test sum(is -> LibPETSc.ISGetSize(petsclib, is), owned) == 6
        LibPETSc.ISColoringDestroy(petsclib, coloring)
        foreach(PETSc.destroy!, owned)
        PETSc.destroy!(da)
    end

    @testset "readers hand back a borrowed PetscVec" begin
        PetscScalar = petsclib.PetscScalar
        A = PETSc.PetscMat(petsclib, PetscScalar[2 0; 0 2])
        ksp = PETSc.KSP(A)
        x = PETSc.PetscVec(petsclib, PetscScalar[0, 0])
        b = PETSc.PetscVec(petsclib, PetscScalar[2, 4])
        PETSc.solve!(x, ksp, b)

        snes = PETSc.SNES(petsclib, comm)
        LibPETSc.SNESSetSolution(petsclib, snes, x)

        ts = PETSc.TS(petsclib, comm)
        PETSc.set_solution!(ts, x)
        PETSc.set_tolerances!(ts; vatol = x)

        da = PETSc.DMDA(petsclib, comm, (PETSc.DM_BOUNDARY_NONE,), (4,), 1, 1)
        PETSc.set_uniform_coordinates!(da, (0.0,), (1.0,))

        for v in (PETSc.solution(ksp), PETSc.solution(snes), PETSc.solution(ts),
                  PETSc.tolerances(ts).vatol, PETSc.local_coordinates(da))
            @test v isa LibPETSc.PetscVec
            @test !PETSc.owns(v)
        end
        @test PETSc.solution(ksp).ptr == x.ptr
        @test PETSc.tolerances(ts).vrtol.ptr == C_NULL
        foreach(PETSc.destroy!, (da, ts, snes, ksp, b, x, A))
    end

    PETSc.finalize(petsclib)
end
