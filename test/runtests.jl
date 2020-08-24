using Test
using PETSc, MPI, LinearAlgebra

MPI.Init()
PETSc.initialize()

x = randn(100)
v = PETSc.SeqVec(MPI.COMM_SELF, x)

@test norm(x) ≈ norm(v) rtol=eps()


PETSc.finalize()
