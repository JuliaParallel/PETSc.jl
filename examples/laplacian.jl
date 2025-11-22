using Test
using PETSc, MPI, LinearAlgebra, SparseArrays

MPI.Initialized() || MPI.Init()

comm = MPI.COMM_WORLD

# Constructs the central finite-difference approximation of -∇² for an lx×ly box of nx×ny points
# with Dirichlet boundary conditions. Based on the code from
# http://math.mit.edu/~stevenj/18.303/lecture-10.html
function laplace(T, nx::Integer, ny::Integer=nx, lx::Integer=1, ly::Integer=lx)
    dx = T(lx) / T(nx + 1)
    dy = T(ly) / T(ny + 1)
    Dx = [[one(T) spzeros(T, 1, nx - 1)]; spdiagm(1=>ones(T, nx - 1)) - I] / dx
    Dy = [[one(T) spzeros(T, 1, ny - 1)]; spdiagm(1=>ones(T, ny - 1)) - I] / dy
    Ax = Dx' * Dx
    Ay = Dy' * Dy
    A = kron(sparse(I, ny, ny), Ax) + kron(Ay, sparse(I, nx, nx))
    return A
end

nx = 20
for T in [Float64, Float32]
  if T==Float64
    petsclib = PETSc.petsclibs[1];
  elseif T==Float32
    petsclib = PETSc.petsclibs[2];
  end
  PETSc.initialize(petsclib)

  S = laplace(T, nx)

  m, n = size(S)
  #M = PETSc.MatSeqAIJ(S)
  M = PETSc.MatCreateSeqAIJ(petsclib, comm, S)

  ksp = PETSc.KSP(
                  M;
                  ksp_monitor_true_residual = false,
                  ksp_view = false,
                  ksp_type = "cg",
                  ksp_rtol = 1e-8,
                  pc_type = "gamg",
                  mg_levels_ksp_type = "chebyshev",
                  mg_levels_ksp_max_it = 3,
                  mg_levels_pc_type = "bjacobi",
                  mg_levels_sub_pc_type = "icc",
                  mg_coarse_ksp_type = "preonly",
                  mg_coarse_pc_type = "cholesky",
                 )

  # initial guess
  x = LibPETSc.VecCreateSeqWithArray(petsclib,comm, 1, m, zeros(T, m))

  b = LibPETSc.VecCreateSeqWithArray(petsclib,comm, 1, n, randn(T, n))
  PETSc.solve!(x, ksp, b)

  
  PETSc.destroy(ksp)
  PETSc.destroy(M)
  PETSc.destroy(b)
  PETSc.destroy(x)

  PETSc.finalize(petsclib)
  
end
