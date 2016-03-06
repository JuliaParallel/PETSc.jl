# solve Problem Set 6 Problem 1
#using ArrayViews  # non-copying subarray package
using PETSc

function driver()
  xmin = 0
  xmax = 1
  N = 50  # N+1 = # of grid points per processor
  mpi_rank = MPI.Comm_rank(MPI.COMM_WORLD) + 1
  mpi_size = MPI.Comm_size(MPI.COMM_WORLD)

  Npts_global = (mpi_size -1)*(N+1) + N + 2 - 2
  delta_x = (xmax - xmin)/(Npts_global)
  delta_t = delta_x

  tmax = 4.0

  solve(xmin, xmax, tmax, N, delta_t)
end


function solve(xmin, xmax, tmax, N, delta_t)
# solves the advection-diffusion equation using centered differences in space
# Crank-Nicolson in time
# Dirchlet boundary condition on left, Neumann on the right
# Method of Manufactured Solutions used to verify convergence rate

# Parameters:
# xmin = minimum x coordinate
# xmax = maximum x coordinate
# tmax = maximum time value
# N : N+1 =  number of x points
# r = delta_t/delta_x^2
# sigma = delta_t/delta_x

# Helper Functions:
# ICFunc = initial condition function with signature val = ICFunc(x)
# BCL = left boundary condition function with signature val = BCL(t)
# BCR = right boundary condition function
# source = source function, signature val = source(x, t)
# this function assumes no ghost point on the left, ghost point on the right

# calculate some values
mpi_rank = MPI.Comm_rank(MPI.COMM_WORLD) + 1
mpi_size = MPI.Comm_size(MPI.COMM_WORLD)

if mpi_rank != mpi_size
  mat_size = N+1
else  # last segment has ghost point
  mat_size = N+2
end

Npts_global = (mpi_size -1)*(N+1) + N + 2 - 2
delta_x = (xmax - xmin)/Npts_global

r = delta_t/(delta_x^2)
sigma = delta_t/delta_x
nStep = convert(Int, div(tmax, delta_t))

if mpi_rank == 1
println("tmax = ", tmax)
println("delta_x = ", delta_x)
println("delta_t = ", delta_t)
println("r = ", r)
println("sigma = ", sigma)
println("nStep = ", nStep)
end


# get the Petsc objects
A, u_i, u_ghost, vec_scatter, rhs, ksp, uex, ctx, ghost_offset = createPetscData(mat_size, 3)

# get the ranges to iterate over
idx_range = localpart(u_i)
idx_start = idx_range[1]
idx_end = idx_range[end]
interior_part = (idx_start+1):(idx_end-1)

# apply IC
# Not applying BCL at initial condition
for i in localpart(u_i)
  u_i[i] = ICFunc(xmin + (i-1)*delta_x)
end
# assemble so we can access the values in u_i for the Neumann BC
PETSc.AssemblyBegin(u_i)
PETSc.AssemblyEnd(u_i)

if mpi_rank == mpi_size
  u_i[idx_end] = 2*delta_x*BCR(0) + u_i[idx_end-2]
end

# scatter the initial condition u_i to u_ghost, so all procs have the values
# needed to evaluate their stencils
scatter!(vec_scatter, u_i, u_ghost)


# populate LHS matrix A
# do left BC
if mpi_rank == 1
#  println("applying left BC")
  A[1, 1] = 1
end

# do right BC
if mpi_rank == mpi_size
  A[idx_end, idx_end] = 1/(2*delta_x)
  A[idx_end, idx_end-2] = -1/(2*delta_x)  # this BC makes A not Tridiagonal
end

# to matrix interior
stencil_l = -r/2 - sigma/4
stencil_c = 1 + r
stencil_r = -r/2 + sigma/4
for i in interior_part
  A[i, i-1] = stencil_l
  A[i, i] = stencil_c
  A[i, i+1] = stencil_r
end

# do boundary nodes of each proc
if mpi_rank > 1  # all processes except leftmost
#  println("p$mpi_rank doing idx_start")
  A[idx_start, idx_start-1] = stencil_l
  A[idx_start, idx_start] = stencil_c
  A[idx_start, idx_start+1] = stencil_r
end
if mpi_rank != mpi_size # all processes except rightmost
  A[idx_end, idx_end-1] = stencil_l
  A[idx_end, idx_end] = stencil_c
  A[idx_end, idx_end+1] = stencil_r
end

# create rhs stencil
stencil_l = r/2 + sigma/4
stencil_c = 1 - r
stencil_r = r/2 - sigma/4

print("\n")
time = @elapsed for tstep=1:nStep  # loop over timesteps
# advance from timestep tstep to tstep + 1

  if mpi_rank == 1
    println("tstep = ", tstep)
  end

  for i in interior_part
    u_k = u_ghost[ i + ghost_offset ]
    u_k_1 = u_ghost[ i-1 + ghost_offset ]
    u_k_p1 = u_ghost[ i+1 + ghost_offset ]
    src_val = source(xmin + (i-1)*delta_x, (tstep - 0.5)*delta_t)
    rhs[i] = stencil_l*u_k_1 + stencil_c*u_k  + stencil_r*u_k_p1 + delta_t*src_val
  end

  if mpi_rank > 1  # all processes except leftmost
    i = idx_start
    u_k = u_ghost[i + ghost_offset]
    u_k_1 = u_ghost[i-1 + ghost_offset]
    u_k_p1 = u_ghost[i+1 + ghost_offset]
    src_val = source(xmin + (i-1)*delta_x, (tstep - 0.5)*delta_t)
    rhs[i] = stencil_l*u_k_1 + stencil_c*u_k  + stencil_r*u_k_p1 + delta_t*src_val
  end

  if mpi_rank < mpi_size  # all processes except rightmost
    i = idx_end
    u_k = u_ghost[i + ghost_offset]
    u_k_1 = u_ghost[i-1 + ghost_offset]
    u_k_p1 = u_ghost[i+1 + ghost_offset]
    src_val = source(xmin + (i-1)*delta_x, (tstep - 0.5)*delta_t)
    rhs[i] = stencil_l*u_k_1 + stencil_c*u_k  + stencil_r*u_k_p1 + delta_t*src_val
  end

  # apply BC terms to rhs
  if mpi_rank == 1
    rhs[1] = BCL(tstep*delta_t)
  end
  if mpi_rank == mpi_size
    rhs[idx_end] = BCR(tstep*delta_t)
  end


  # solve for next time step, u_i gets overwritten with new solution
  A_ldiv_B!(ksp, rhs, u_i)

  # scatter u_i into u_ghost, which is used to calculate the rhs
  # this is necessary because Petsc does not allow fetching off process values
  scatter!(vec_scatter, u_i, u_ghost)

end

println("time = ", time)

# print numerical solution
PETSc.AssemblyBegin(u_i)
PETSc.AssemblyEnd(u_i)
if mpi_rank == 1
  println("u_final = \n")
end
petscview(u_i)


# get exact solution
for i in localpart(uex)
  uex[i] = uExact( (i-1)*delta_x, delta_t*nStep)
end

# print exact solution
PETSc.AssemblyBegin(uex)
PETSc.AssemblyEnd(uex)
if mpi_rank == 1
  println("u_ex = \n")
end
petscview(uex)

# calculate error
u_err = uex - u_i
max_err = norm(u_err, Inf)
sleep(2)
if mpi_rank == 1
  println("max error = ", max_err)
end

writeU(u_i, N, nStep)
MPI.Barrier(MPI.COMM_WORLD)

return u_i, delta_t*(nStep)

end

function writeU(u_i, N, tsteps)

  mpi_rank = MPI.Comm_rank(MPI.COMM_WORLD) + 1
  larr = LocalArray(u_i)
  fname = string("u_", tsteps, "_", mpi_rank, ".dat")
  writedlm(fname, larr[1:(N+1)])
end

  


function ICFunc(x)
  return 2*cos(3*x)
end

function BCZero(x)
  return 0.0
end

function BCL(t)
  return 2*cos(t)
end

function BCR(t)
  return -6*sin(3)*cos(t)
end

function source(x, t)
  return -2*cos(3*x)*sin(t) - 6*sin(3x)*cos(t) + 18*cos(3*x)*cos(t)
end

function uExact(x, t)
  return 2*cos(3*x)*cos(t)
end

function createPetscData(mat_size, stencil_size)
# create all the Petsc objects needed

  mpi_rank = MPI.Comm_rank(MPI.COMM_WORLD) + 1
  mpi_size = MPI.Comm_size(MPI.COMM_WORLD)

  # left hand side matrix
  A = PETSc.Mat(Float64, PETSc.C.PETSC_DECIDE, PETSc.C.PETSC_DECIDE, mlocal=mat_size, nlocal=mat_size, nz=stencil_size, onz=1, comm=MPI.COMM_WORLD)
  A.assembling=true
  
  # solution at current timestep
  u_i = Vec(Float64, PETSc.C.VECMPI, comm=MPI.COMM_WORLD);
  u_i.assembling=true
  resize!(u_i, mlocal=mat_size)

  # right hand side
  rhs = Vec(Float64, PETSc.C.VECMPI, comm=MPI.COMM_WORLD);
  rhs.assembling=true
  resize!(rhs, mlocal=mat_size)

  # exact solution
  uex = Vec(Float64, PETSc.C.VECMPI, comm=MPI.COMM_WORLD);
  uex.assembling=true
  resize!(uex, mlocal=mat_size)

  # create a ghost vector that includes all points needed for each proc
  # to evaluate the rhs, duplicating the shared points.
  # Set up scattering from u_i to this vector

  # get global indices of local part of u_i
  idx_range = localpart(u_i)
  start_idx = idx_range[1] - 1
  end_idx = idx_range[end] - 1
  # figure out the indices in u_i corresponding to u_ghost
  if mpi_rank == 1  # left boundary (Dirchlet BC)
    idx = collect(start_idx:(end_idx+1))
    ghost_startidx = 0
    ghost_endidx = mat_size
    idx_ghost = collect(ghost_startidx:ghost_endidx)
  elseif mpi_rank == mpi_size  # right boundary (Neumann BC)
    idx = collect( (start_idx-1):end_idx)
    ghost_startidx = mat_size + (mpi_size-2)*(mat_size + 1)
    ghost_endidx = ghost_startidx + mat_size
    idx_ghost = collect(ghost_startidx:ghost_endidx)
  else  # interior
    idx = collect((start_idx-1):(end_idx+1))
    ghost_startidx = mat_size + 1 + (mpi_rank -2)*(mat_size + 2)
    ghost_endidx = ghost_startidx + mat_size + 1
    idx_ghost = collect(ghost_startidx:ghost_endidx)
  end
  # bump everything up to 1-based indexing
  idx += 1
  idx_ghost += 1

  ghost_offset = 2*(mpi_rank-1) # offset from u_i indices to ghost indices
  # set up the vector
  nghost_local = length(idx_ghost)
  u_ghost = Vec(Float64, PETSc.C.VECMPI, comm=MPI.COMM_WORLD);
  resize!(u_ghost, mlocal=nghost_local)
  u_ghost.assembling=true

  # create the scatter
  is_local = PETSc.IS(Float64, idx, comm=MPI.COMM_WORLD)
  is_ghost = PETSc.IS(Float64, idx_ghost, comm=MPI.COMM_WORLD)
  vec_scatter = PETSc.VecScatter(u_i, is_local, u_ghost, is_ghost)

  ctx = (is_local, is_ghost)  # avoid GC
  ksp = PETSc.KSP(A, ksp_atol=1e-12, ksp_rtol=1e-14, ksp_monitor="")

  return A, u_i, u_ghost, vec_scatter, rhs, ksp, uex, ctx, ghost_offset
end

# run
driver()
