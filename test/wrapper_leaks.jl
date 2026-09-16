using Test
using PETSc, MPI
using PETSc: LibPETSc

MPI.Initialized() || MPI.Init()

# ============================================================================
#   Leak check for the generated wrappers
# ============================================================================
#
# The optimized PETSc_jll build has no malloc/object tracing, so leaks are detected by repeating a
# create/destroy or Get/Restore pair many times and checking that the process memory PETSc reports
# does not grow with the iteration count. A leaked Vec of `n` scalars per iteration would grow the
# resident set by `iters * n * sizeof(PetscScalar)`, far above the tolerance.

function rss(petsclib)
    GC.gc()
    return LibPETSc.PetscMemoryGetCurrentUsage(petsclib)
end

@testset "repeated wrapper calls do not leak" begin
    comm = LibPETSc.PETSC_COMM_SELF
    iters = 2000
    n = 4000                       # 32 KB per Float64 vector: a leak would show as 64 MB per block
    for petsclib in PETSc.petsclibs[1:2]
        PETSc.initialize(petsclib)
        PetscInt = petsclib.PetscInt
        PetscScalar = petsclib.PetscScalar
        leak_budget = 32 * 2^20    # 32 MB of allocator noise allowed; a Vec leak would be 64 MB

        pairs = Dict{String,Function}(
            "VecCreateSeq/VecDestroy" => () -> begin
                v = LibPETSc.VecCreateSeq(petsclib, comm, PetscInt(n))
                LibPETSc.VecDestroy(petsclib, v)
            end,
            "VecDuplicate/VecDestroy" => () -> begin
                v = LibPETSc.VecCreateSeq(petsclib, comm, PetscInt(n))
                w = LibPETSc.VecDuplicate(petsclib, v)
                LibPETSc.VecDestroy(petsclib, w); LibPETSc.VecDestroy(petsclib, v)
            end,
            "VecGetArray/VecRestoreArray" => () -> begin
                v = LibPETSc.VecCreateSeq(petsclib, comm, PetscInt(n))
                a = LibPETSc.VecGetArray(petsclib, v); a[1] = 1
                LibPETSc.VecRestoreArray(petsclib, v, a)
                LibPETSc.VecDestroy(petsclib, v)
            end,
            "ISCreateStride/ISGetIndices/ISRestoreIndices/ISDestroy" => () -> begin
                is = LibPETSc.ISCreateStride(petsclib, comm, PetscInt(n), PetscInt(0), PetscInt(1))
                idx = LibPETSc.ISGetIndices(petsclib, is)
                LibPETSc.ISRestoreIndices(petsclib, is, idx)
                LibPETSc.ISDestroy(petsclib, is)
            end,
            "MatCreateSeqAIJ/MatDestroy" => () -> begin
                A = LibPETSc.MatCreateSeqAIJ(petsclib, comm, PetscInt(n), PetscInt(n), PetscInt(1), C_NULL)
                LibPETSc.MatDestroy(petsclib, A)
            end,
            "KSPCreate/KSPDestroy" => () -> begin
                ksp = LibPETSc.KSPCreate(petsclib, comm)
                LibPETSc.KSPDestroy(petsclib, ksp)
            end,
            "DMDACreate1d/DMGetCoordinatesLocal/DMDestroy" => () -> begin
                dm = LibPETSc.DMDACreate1d(petsclib, comm, LibPETSc.DM_BOUNDARY_NONE, PetscInt(n), PetscInt(1), PetscInt(1), C_NULL)
                LibPETSc.DMSetUp(petsclib, dm)
                LibPETSc.DMDASetUniformCoordinates(petsclib, dm, petsclib.PetscReal.((0, 1, 0, 1, 0, 1))...)
                c = LibPETSc.DMGetCoordinatesLocal(petsclib, dm)   # borrowed, not destroyed
                LibPETSc.DMDestroy(petsclib, dm)
            end,
        )
        for (name, f) in pairs
            # first block: compilation, PETSc internal caches and allocator arenas settle;
            # second block: steady state, where a leak shows as growth proportional to `iters`
            for _ in 1:iters
                f()
            end
            before = rss(petsclib)
            for _ in 1:iters
                f()
            end
            growth = rss(petsclib) - before
            @test growth < leak_budget
            growth < leak_budget || @info "possible leak" name growth iters
        end
        PETSc.finalize(petsclib)
    end
end
