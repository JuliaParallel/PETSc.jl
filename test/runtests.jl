using Test
using PETSc, MPI, LinearAlgebra, SparseArrays

m,n = 20,20
x = randn(n)
v = PETSc.VecSeq(x)

# PETSc.view(v)

@test norm(x) ≈ norm(v) rtol=10eps()

S = sprand(m,n,0.1) + I
M = PETSc.MatSeqAIJ(S)

@test norm(S) ≈ norm(M) rtol=10eps()

w = PETSc.VecSeq(zeros(m))

mul!(w, M, v)
@test w.array ≈ S*x 

ksp = PETSc.KSP(M; ksp_rtol=1e-8, pc_type="jacobi", ksp_monitor=true)
#PETSc.settolerances!(ksp; rtol=1e-8)

@test PETSc.gettype(ksp) == "gmres" # default

pc = PETSc.PC(ksp)
@test PETSc.gettype(pc) == "jacobi"

u = PETSc.VecSeq(randn(n))
y = ksp \ w.array
@test S*y ≈ w.array  rtol=1e-8

mul!(w, transpose(M), v)
@test w.array ≈ transpose(S)*x 

PETSc.solve!(u, transpose(ksp), w)
@test transpose(S)*u.array ≈ w.array  rtol=1e-8
