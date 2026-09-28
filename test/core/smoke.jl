# test/core/smoke.jl
# Vec, Mat, KSP and MatShell end to end on small problems, against Julia's own
# linear algebra: the quickest check that the library works at all.

using Test
using PETSc, MPI, LinearAlgebra, SparseArrays

MPI.Initialized() || MPI.Init()

@testset "Vec, Mat and KSP end to end" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    PetscInt = petsclib.PetscInt
    comm = MPI.COMM_SELF

    n = 20
    x = randn(n)
    V = LibPETSc.VecCreateSeqWithArray(petsclib, comm, PetscInt(1), PetscInt(n), x)
    @test norm(x) ≈ norm(V) rtol = 10eps()

    S = sprand(n, n, 0.1) + I
    M = PETSc.PetscMat(petsclib, comm, S)
    @test norm(S) ≈ norm(M) rtol = 10eps()

    vec_x = LibPETSc.VecCreateSeqWithArray(petsclib, comm, PetscInt(1), PetscInt(n), copy(x))
    w = M * vec_x
    @test w[:] ≈ S * x

    ksp = PETSc.KSP(M; ksp_rtol = 1e-8, pc_type = "jacobi", ksp_monitor = false)
    @test PETSc.type_name(ksp) === :gmres   # the default
    y = ksp \ w
    @test S * y[:] ≈ w[:] rtol = 1e-6

    # a matrix-free operator
    f!(y, x) = y .= 2 .* x
    shell = PETSc.MatShell(petsclib, f!, comm, 10, 10)
    z = rand(10)
    @test shell * z ≈ 2z
    @test PETSc.KSP(shell) \ z ≈ z / 2

    foreach(PETSc.destroy!, (ksp, shell, M, vec_x, V))
    PETSc.finalize(petsclib)
end
