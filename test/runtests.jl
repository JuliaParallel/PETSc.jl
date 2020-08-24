using Test
using PETSc, MPI, LinearAlgebra, SparseArrays

MPI.Init()
PETSc.initialize()


m,n = 20,20
x = randn(n)
v = PETSc.SeqVec(MPI.COMM_SELF, x)

@test norm(x) ≈ norm(v) rtol=10eps()

S = sprand(m,n,0.1)
M = PETSc.SeqAIJMat(S)

@test norm(S) ≈ norm(M) rtol=10eps()

w = PETSc.SeqVec(MPI.COMM_SELF, zeros(m))

mul!(w, M, v)
@test w.array ≈ S*x 

mul!(w, transpose(M), v)
@test w.array ≈ transpose(S)*x 


PETSc.destroy(M)
PETSc.destroy(v)
PETSc.destroy(w)

PETSc.finalize()
