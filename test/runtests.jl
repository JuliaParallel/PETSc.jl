using Test
using PETSc, MPI, LinearAlgebra, SparseArrays

MPI.Init()
PETSc.initialize()


m,n = 20,20
x = randn(n)
v = PETSc.VecSeq(x)

@test norm(x) ≈ norm(v) rtol=10eps()

S = sprand(m,n,0.1) + I
M = PETSc.MatSeqAIJ(S)

@test norm(S) ≈ norm(M) rtol=10eps()

w = PETSc.VecSeq(zeros(m))

mul!(w, transpose(M), v)
@test w.array ≈ transpose(S)*x 


mul!(w, M, v)
@test w.array ≈ S*x 


ksp = PETSc.KSP(M)

pc = PETSc.PC(ksp)
PETSc.settype!(pc, "jacobi")
PETSc.settolerances!(ksp; rtol=1e-8)

u = PETSc.VecSeq(randn(n))
PETSc.solve!(u, ksp, w)

@test S*u.array ≈ w.array  rtol=1e-8



PETSc.destroy(ksp)
PETSc.destroy(M)
PETSc.destroy(u)
PETSc.destroy(v)
PETSc.destroy(w)

PETSc.finalize()
