# solve Problem Set 6 Problem 1
using ArrayViews  # non-copying subarray package
using PETSc

function driver()
  xmin = 0
  xmax = 1
  N = 110  # N+1 = # of grid points per processor
  mpi_rank = MPI.Comm_rank(MPI.COMM_WORLD) + 1
  mpi_size = MPI.Comm_size(MPI.COMM_WORLD) 

  Npts_global = (mpi_size -1)*(N+1) + N + 2 - 2
  delta_x = (xmax - xmin)/(Npts_global)
#  r = 0.5
#  sigma = 0.75
  delta_t = delta_x
#  nu = 1/6


  tmax = 1.0
  ICFunc = IC1
  BCL = BC1
  BCR = BC2
  src = SRC

#  u, tmax_ret = solve(xmin, xmax, tmax, N, delta_t, ICFunc, BCL, BCR, src)

  solve(xmin, xmax, tmax, N, delta_t, ICFunc, BCL, BCR, src)
#=
#  u = solve(xmin, xmax, tmax, N, delta_t, nu, ICFunc, BCL, BCR)
  println("tmax_ret = ", tmax_ret)
  max_err, l2_err = calcError(u, xmin, xmax, tmax_ret, N)
  vals = [xmin, xmax, tmax_ret, delta_t]
  writedlm("counts.dat", vals)
  writedlm("u.dat", u)

  f = open("convergence.dat", "a+")
  @printf(f, "%d %16.15f %16.15f\n", N, max_err, tmax_ret)
  close(f)
=#
end

function calcError(u, xmin, xmax, tmax, N)
# calculate the error
# u should not include the ghost point

  print("\n")
  delta_x = (xmax - xmin)/N
  err = zeros(size(u))
  u_ex = zeros(size(u))
  for i=1:length(u)
    x_i = xmin + (i-1)*delta_x
    u_exact_i = uExact(x_i, tmax)
    u_ex[i] = u_exact_i
    err[i] = u_exact_i - u[i]
  end

  err_max = maximum(abs(err))
  l2_err = norm(err)

  println("err = \n", err)
  println("u = \n", u)
  println("u_ex = \n", u_ex)



  println("max error = ", err_max)
  println("L2 error = ", l2_err)

  writedlm("u.dat", u)
  writedlm("uexact.dat", u_ex)
  return err_max, l2_err
end



function solve(xmin, xmax, tmax, N, delta_t,ICFunc::Function, BCL::Function, BCR::Function, source::Function)
# xmin = minimum x coordinate
# xmax = maximum x coordinate
# tmax = maximum time value
# N : N+1 =  number of x points
# r = nu*delta_t/delta_x^2
# sigma = delta_t/delta_x
# ICFunc = initial condition function with signature val = ICFunc(x)
# BCL = left boundary condition function with signature val = BCL(t)
# BCR = right boundary condition function
# source = source function, signature val = source(x, t)
# this function assumes no ghost point on the left, ghost point on the right

mpi_rank = MPI.Comm_rank(MPI.COMM_WORLD) + 1
mpi_size = MPI.Comm_size(MPI.COMM_WORLD) 

if mpi_rank != mpi_size
  mat_size = N+1
else  # last segment has ghost point
  mat_size = N+2
end

Npts_global = (mpi_size -1)*(N+1) + N + 2 - 2
delta_x = (xmax - xmin)/Npts_global

#delta_t = (r*delta_x^2)/nu  # nu*delta_t
r = delta_t/(delta_x^2)
sigma = delta_t/delta_x
nStep = convert(Int, div(tmax, delta_t))

println("tmax = ", tmax)
println("delta_x = ", delta_x)
println("delta_t = ", delta_t)
println("r = ", r)
println("sigma = ", sigma)
println("nStep = ", nStep)


# allocate storage
#A = zeros(Float64, mat_size, mat_size)  # this could be a SparseMatrixCSC
#rhs = Array(Float64, mat_size) # right hand side
#u_i = Array(Float64, mat_size)  # current timestep solution values
A, u_i, u_ghost, vec_scatter, rhs, ksp, uex, ctx, ghost_offset = createPetscData(mat_size, 3)

#createPetscData(mat_size, 3)
println("p$mpi_rank exited createPetscData")
flush(STDOUT)
#end

#sleep(3*mpi_rank)
idx_range = localpart(u_i)
#sleep(2*mpi_rank)
#println("u_i.assembling = ", u_i.assembling)
#println("p$mpi_rank finished sleeping")
#println("p$mpi_rank idx_range = ", idx_range)
idx_start = idx_range[1]
idx_end = idx_range[end]
#println("p$mpi_rank idx start, end = ", idx_start, ", ", idx_end)
interior_part = (idx_start+1):(idx_end-1)
#println("p$mpi_rank interior part = ", interior_part)

#sleep(10)
# apply IC
# Not applying BCL at initial condition
for i in localpart(u_i)
#  if mpi_rank == 3
    println("p$mpi_rank i = ", i, " x = ", (i-1)*delta_x)
#  end
  u_i[i] = val =  ICFunc(xmin + (i-1)*delta_x)
  println("  val = ", val)
end
PETSc.AssemblyBegin(u_i)
PETSc.AssemblyEnd(u_i)



#sleep(5*mpi_rank)
#println("p$mpi_rank finished applying IC")
# set ghost point value at IC
if mpi_rank == mpi_size
#  println("p$mpi_rank applying ghost point")
  u_i[idx_end] = 2*delta_x*BCR(0) + u_i[idx_end-2]
#  println("p$mpi_rank finished applying ghost point")
end
#println("p$mpi_rank finished applyc BC at IC")

PETSc.AssemblyBegin(u_i)
PETSc.AssemblyEnd(u_i)
if mpi_rank == 1
  println("u_initial = \n")
end
petscview(u_i)


#println("\nu_initial = \n", u_i)
#println("finished printing u initial")
#println("p$mpi_rank about to scatter")
scatter!(vec_scatter, u_i, u_ghost)
#sleep(1*mpi_rank)
#println("p$mpi_rank finished scatter")


#sleep(1*mpi_rank)

# apply BC to IC

# populate the matrix A
# do left BC
if mpi_rank == 1
#  println("applying left BC")
  A[1, 1] = 1
end

# do right BC
if mpi_rank == mpi_size
#  println("applying right BC")
  A[idx_end, idx_end] = 1/(2*delta_x)
  A[idx_end, idx_end-2] = -1/(2*delta_x)  # this BC makes A not Tridiagonal
end


stencil_l = -r/2 - sigma/4
stencil_c = 1 + r
stencil_r = -r/2 + sigma/4
#println("p$mpi_rank about to do interior part")
for i in interior_part  # loop over interior of matrixa
#  println("p$mpi_rank i = ", i)
  A[i, i-1] = stencil_l
  A[i, i] = stencil_c
  A[i, i+1] = stencil_r
end

#println("p$mpi_rank finished doing interior part")

if mpi_rank > 1  # all processes except leftmost
#  println("p$mpi_rank doing idx_start")
  A[idx_start, idx_start-1] = stencil_l
  A[idx_start, idx_start] = stencil_c
  A[idx_start, idx_start+1] = stencil_r
end
if mpi_rank != mpi_size # all processes except rightmost
#  println("p$mpi_rank doing idx_end = ", idx_end)
  A[idx_end, idx_end-1] = stencil_l
  A[idx_end, idx_end] = stencil_c
  A[idx_end, idx_end+1] = stencil_r
end



  println("A = \n")
  PETSc.AssemblyBegin(A, PETSc.C.MAT_FINAL_ASSEMBLY)
  PETSc.AssemblyEnd(A, PETSc.C.MAT_FINAL_ASSEMBLY)
  petscview(A)




#Af = lufact(A)
#println("Af = ", Af)
#println("typeof(Af) = ", typeof(Af))


# set up stencil for rhs
stencil_l = r/2 + sigma/4
stencil_c = 1 - r
stencil_r = r/2 - sigma/4


#println("\nstencil_l = ", stencil_l)
#println("stencil_c = ", stencil_c)
#println("stencil_r = ", stencil_r)
print("\n")

time = @elapsed for tstep=1:nStep  # loop over timesteps
# advance from timestep tstep to tstep + 1

  if mpi_rank == 1
    println("\ntstep = ", tstep)
  end
#  print("\n")
  # it shouldn't be necessary to apply the BC at subsequent time steps
  # because the BC is built into the matrix A
  # print verification here
#  uL = BCL( (tstep-1)*delta_t )
#  println("uL = ", uL, "u[1] = ", u_i[1])
#  ghost_val = 2*delta_x*BCR( (tstep-1)*delta_t ) + u_i[mat_size - 2]
#  println("ghost_val = ", ghost_val, " u[mat_size] = ", u_i[mat_size])

#  sleep(3*mpi_rank)
#  println("p$mpi_rank finished sleeping")
#  println("u_i = ", u_i)
  # calculate right hand size interior points
  ghost_i = 2  # local numbering of points aka numbering in ghost vector
  for i in interior_part
    println("i = ", i)
    
    u_k = u_ghost[ i + ghost_offset ]
    u_k_1 = u_ghost[ i-1 + ghost_offset ]
    u_k_p1 = u_ghost[ i+1 + ghost_offset ]
    src_val = source( (i-1)*delta_x, (tstep - 0.5)*delta_t) 
    rhs[i] = stencil_l*u_k_1 + stencil_c*u_k  + stencil_r*u_k_p1 + delta_t*src_val
    ghost_i += 1
  end

  if mpi_rank > 1  # all processes except leftmost
    i = idx_start
    u_k = u_ghost[i + ghost_offset]
    u_k_1 = u_ghost[i-1 + ghost_offset]
    u_k_p1 = u_ghost[i+1 + ghost_offset]
    src_val = source( (i-1)*delta_x, (tstep - 0.5)*delta_t) 
    rhs[i] = stencil_l*u_k_1 + stencil_c*u_k  + stencil_r*u_k_p1 + delta_t*src_val
  end

  if mpi_rank < mpi_size  # all processes except rightmost
    i = idx_end
    u_k = u_ghost[i + ghost_offset]
    u_k_1 = u_ghost[i-1 + ghost_offset]
    u_k_p1 = u_ghost[i+1 + ghost_offset]
    src_val = source( (i-1)*delta_x, (tstep - 0.5)*delta_t) 
    rhs[i] = stencil_l*u_k_1 + stencil_c*u_k  + stencil_r*u_k_p1 + delta_t*src_val
 
  end

  # apply BC terms to rhs
  if mpi_rank == 1
    rhs[1] = BCL(tstep*delta_t)
  end
  if mpi_rank == mpi_size
    rhs[idx_end] = BCR(tstep*delta_t)
  end

  PETSc.AssemblyBegin(rhs)
  PETSc.AssemblyEnd(rhs)
  if mpi_rank == 1
    println("rhs = \n")
  end
  petscview(rhs)





#  println("rhs = \n", rhs)
  # solve for next time step
#  A_ldiv_B!(Af, rhs)  # rhs gets overwritten with new solution values
  A_ldiv_B!(ksp, rhs, u_i)
#  u_i = A\rhs
#  println("secondard u_next = \n", u_next)

  PETSc.AssemblyBegin(u_i)
  PETSc.AssemblyEnd(u_i)
  if mpi_rank == 1
    println("u_initial = \n")
  end
  petscview(u_i)



  println("p$mpi_rank about to scatter")
  scatter!(vec_scatter, u_i, u_ghost)
  println("p$mpi_rank finished scatter")

  PETSc.AssemblyBegin(u_ghost)
  PETSc.AssemblyEnd(u_ghost)
  if mpi_rank == 1
    println("u_ghost = \n")
  end
  petscview(u_ghost)


#  println("u_next = \n", u_i)

end


println("time = ", time)

# get exact solution
for i in localpart(uex)
  uex[i] = uExact( (i-1)*delta_x, delta_t*nStep)
end

PETSc.AssemblyBegin(u_i)
PETSc.AssemblyEnd(u_i)
if mpi_rank == 1
  println("u_final = \n")
end
petscview(u_i)

PETSc.AssemblyBegin(uex)
PETSc.AssemblyEnd(uex)
if mpi_rank == 1
  println("u_ex = \n")
end
petscview(uex)




u_err = uex - u_i
max_err = norm(u_err, Inf)
sleep(2)
if mpi_rank == 1
  println("max error = ", max_err)
else
  sleep(5)
end


return u_i, delta_t*(nStep)  # plus 1 because we are at the beginning of the next timestep

end



function IC1(x)
  return 2*cos(3*x)
end

function BCZero(x)
  return 0.0
end

function BC1(t)
  return 2*cos(t)
end

function BC2(t)
  return -6*sin(3)*cos(t)
end

function SRC(x, t)
  return -2*cos(3*x)*sin(t) - 6*sin(3x)*cos(t) + 18*cos(3*x)*cos(t)
end

function uExact(x, t)
  return 2*cos(3*x)*cos(t)
end


#=
function IC1(x)
  return x*x + x + 1
end


function BC1(t)
  return t*t + t + 1
end

function BC2(t)
  return 3
end

function SRC(x, t)
  return 2*x + 2*t
end

function uExact(x, t)
  return x*x + t*t + x + t +1
end
=#

function createPetscData(mat_size, stencil_size)
  
  mpi_rank = MPI.Comm_rank(MPI.COMM_WORLD) + 1
  mpi_size = MPI.Comm_size(MPI.COMM_WORLD) 

#=
  # figure out total number of columns
  if mpi_rank == mpi_size
    ncols = (mat_size-1)*(mpi_size-1) + mat_size
  else
    ncols = mat_size*(mpi_size-1) + mat_size + 1
  end 
=#
  # left hand side matrix
  A = PETSc.Mat(Float64, PETSc.C.PETSC_DECIDE, PETSc.C.PETSC_DECIDE, mlocal=mat_size, nlocal=mat_size, nz=stencil_size, onz=1, comm=MPI.COMM_WORLD)
#  sleep(2*mpi_rank)
  println("p$mpi_rank size(A) = ", size(A))
  println("p$mpi_rank localsize(A) = ", PETSc.sizelocal(A))
  A.assembling=true

  # solution at current timestep
  u_i = Vec(Float64, PETSc.C.VECMPI, comm=MPI.COMM_WORLD); 
  u_i.assembling=true
  setsizes!(u_i, mat_size)

  # right hand side
  rhs = Vec(Float64, PETSc.C.VECMPI, comm=MPI.COMM_WORLD); 
  rhs.assembling=true
  setsizes!(rhs, mat_size)

  # exact solution
  uex = Vec(Float64, PETSc.C.VECMPI, comm=MPI.COMM_WORLD); 
  uex.assembling=true
  setsizes!(uex, mat_size)



  # create vector that includes all points needed for the stencil
  # and set up scattering from u_i to this vector

  # get global indices of local part of u_i
  idx_range = localpart(u_i)
  start_idx = idx_range[1] - 1
  end_idx = idx_range[end] - 1
  # figure out how many points are needed
  if mpi_rank == 1  # left boundary (Dirchlet BC)
    idx = collect(start_idx:(end_idx+1))
    ghost_startidx = 0
    ghost_endidx = mat_size
    idx_ghost = collect(ghost_startidx:ghost_endidx)
  elseif mpi_rank == mpi_size  # right boundary (Neumann BC)
    idx = collect( (start_idx-1):end_idx)
    ghost_startidx = mat_size + (mpi_size-2)*(mat_size + 1)  # for this point, mat_size = 1 + 
                                                # the mat size of hte other blocks
    ghost_endidx = ghost_startidx + mat_size
    idx_ghost = collect(ghost_startidx:ghost_endidx)

  else  # interior
    idx = collect((start_idx-1):(end_idx+1))
    ghost_startidx = mat_size + 1 + (mpi_rank -2)*(mat_size + 2)
    ghost_endidx = ghost_startidx + mat_size + 1
    idx_ghost = collect(ghost_startidx:ghost_endidx)
  end
  ghost_offset = 2*(mpi_rank-1) # offset from u_i indices to ghost indices
  # set up the vector  
  nghost_local = length(idx_ghost)
  u_ghost = Vec(Float64, PETSc.C.VECMPI, comm=MPI.COMM_WORLD); 
  setsizes!(u_ghost, nghost_local)
  u_ghost.assembling=true

  # create the scatter
#  sleep(2*mpi_rank)
  println("p$mpi_rank mat_size = ", mat_size)

  println("p$mpi_rank start_idx = ", start_idx)
  println("p$mpi_rank end_idx = ", end_idx)
  println("p$mpi_rank idx = ", idx)
  println("p$mpi_rank idx_ghost = ", idx_ghost)
  is_local = PETSc.IS(Float64, idx, comm=MPI.COMM_WORLD)
  is_ghost = PETSc.IS(Float64, idx_ghost, comm=MPI.COMM_WORLD)
  vec_scatter = PETSc.VecScatter(u_i, is_local, u_ghost, is_ghost)

  println("finished creating vec scatter")  
  ctx = (is_local, is_ghost)  # avoid GC
  ksp = PETSc.KSP(A)
  println("finished creating ksp")
  setoptions!(ksp; ksp_atol=1e-10, ksp_rtol=1e-8, ksp_monitor="")
  println("Finished setting ksp options")
#  sleep(5)
#  sleep(mpi_rank*5)
  println("p$mpi_rank returning")
  return A, u_i, u_ghost, vec_scatter, rhs, ksp, uex, ctx, ghost_offset
end

# run
driver()
