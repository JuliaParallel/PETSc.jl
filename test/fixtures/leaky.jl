# Fixture for test_audit.jl. Never executed, only parsed.
#
# Exercises what the auditor has to get right: qualified and unqualified calls,
# tuple assignment, several statements on one line, low-level destroys, and
# creations that appear only inside comments or string literals.

using PETSc

function leaky(petsclib, comm)
    # created and destroyed through the high-level name
    ksp = PETSc.KSP(petsclib, comm)
    PETSc.destroy!(ksp)

    # created and destroyed through the low-level name
    v1 = LibPETSc.VecCreateSeq(petsclib, comm, 10)
    LibPETSc.VecDestroy(petsclib, v1)

    # two statements on one line: both must register
    v2 = LibPETSc.VecCreateSeq(petsclib, comm, 10); LibPETSc.VecDestroy(petsclib, v2)

    # tuple assignment creates two objects, only one is released
    x, b = LibPETSc.MatCreateVecs(petsclib, A)
    LibPETSc.VecDestroy(petsclib, x)

    # never destroyed
    dm = PETSc.DMStag(petsclib, comm, (PETSc.DM_BOUNDARY_NONE,), (10,), 1, 1)
    mat = PETSc.MatSeqAIJ(petsclib, 10, 10, 3)

    # a comment mentioning PETSc.destroy!(dm) must not count as a release
    note = "call LibPETSc.MatDestroy(petsclib, mat) when finished"

    return note
end
