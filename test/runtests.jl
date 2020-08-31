using Test
using PETSc, MPI, LinearAlgebra, SparseArrays

m,n = 20,20
x = randn(n)
V = PETSc.VecSeq(x)


@test norm(x) ≈ norm(V) rtol=10eps()

S = sprand(m,n,0.1) + I
M = PETSc.MatSeqAIJ(S)

@test norm(S) ≈ norm(M) rtol=10eps()

w = M*x
@test w ≈ S*x 

ksp = PETSc.KSP(M; ksp_rtol=1e-8, pc_type="jacobi", ksp_monitor=true)
#PETSc.settolerances!(ksp; rtol=1e-8)

@test PETSc.gettype(ksp) == "gmres" # default

pc = PETSc.PC(ksp)
@test PETSc.gettype(pc) == "jacobi"

y = ksp \ w
@test S*y ≈ w  rtol=1e-8


w = M'*x
@test w ≈ S'*x 

y = ksp' \ w
@test S'*y ≈ w rtol=1e-8
