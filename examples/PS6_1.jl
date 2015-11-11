# solve Problem Set 6 Problem 1
using ArrayViews  # non-copying subarray package
using PETSc

function driver()
  xmin = 0
  xmax = 1
  N = 20  # N+1 = # of grid points
  delta_x = (xmax - xmin)/N
#  r = 0.5
#  sigma = 0.75
  delta_t = delta_x
#  nu = 1/6


  tmax = 1.0
  ICFunc = IC1
  BCL = BC1
  BCR = BC2
  src = SRC

  u, tmax_ret = solve(xmin, xmax, tmax, N, delta_t, ICFunc, BCL, BCR, src)
#  u = solve(xmin, xmax, tmax, N, delta_t, nu, ICFunc, BCL, BCR)
  println("tmax_ret = ", tmax_ret)
  max_err, l2_err = calcError(u, xmin, xmax, tmax_ret, N)
  vals = [xmin, xmax, tmax_ret, delta_t]
  writedlm("counts.dat", vals)
  writedlm("u.dat", u)

  f = open("convergence.dat", "a+")
  @printf(f, "%d %16.15f %16.15f\n", N, max_err, tmax_ret)
  close(f)
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
delta_x = (xmax - xmin)/N
#delta_t = (r*delta_x^2)/nu  # nu*delta_t
r = delta_t/(delta_x^2)
sigma = delta_t/delta_x
nStep = convert(Int, div(tmax, delta_t))
mat_size = N+2

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
A, u_i, rhs, ksp = createPetscData(mat_size, 3)


# apply IC
# Not applying BCL at initial condition
for i=1:(mat_size-1)
  u_i[i] = ICFunc(xmin + (i-1)*delta_x)
end

# set ghost point value at IC
u_i[mat_size] = 2*delta_x*BCR(0) + u_i[mat_size-2]


println("\nu_initial = \n", u_i)

# apply BC to IC

# construct the matrix A
# do left BC
A[1, 1] = 1

# do right BC
A[mat_size, mat_size] = 1/(2*delta_x)
A[mat_size, mat_size-2] = -1/(2*delta_x)  # this BC makes A not Tridiagonal


stencil_l = -r/2 - sigma/4
stencil_c = 1 + r
stencil_r = -r/2 + sigma/4
for i=2:(mat_size-1)  # loop over interior of matrix
  A[i, i-1] = stencil_l
  A[i, i] = stencil_c
  A[i, i+1] = stencil_r
end


#println("A = \n", A)

#Af = lufact(A)
#println("Af = ", Af)
#println("typeof(Af) = ", typeof(Af))


# set up stencil for rhs
stencil_l = r/2 + sigma/4
stencil_c = 1 - r
stencil_r = r/2 - sigma/4

println("\nstencil_l = ", stencil_l)
println("stencil_c = ", stencil_c)
println("stencil_r = ", stencil_r)
print("\n")

time = @elapsed for tstep=1:nStep  # loop over timesteps
# advance from timestep tstep to tstep + 1

  println("\ntstep = ", tstep)
#  print("\n")
  # it shouldn't be necessary to apply the BC at subsequent time steps
  # because the BC is built into the matrix A
  # print verification here
  uL = BCL( (tstep-1)*delta_t )
  println("uL = ", uL, "u[1] = ", u_i[1])
  ghost_val = 2*delta_x*BCR( (tstep-1)*delta_t ) + u_i[mat_size - 2]
  println("ghost_val = ", ghost_val, " u[mat_size] = ", u_i[mat_size])


  println("u_i = ", u_i)
  # calculate right hand size interior points
  for i=2:(mat_size-1)
    u_k = u_i[i]
    u_k_1 = u_i[i-1]
    u_k_p1 = u_i[i+1]
    src_val = source( (i-1)*delta_x, (tstep - 0.5)*delta_t) 
    rhs[i] = stencil_l*u_k_1 + stencil_c*u_k  + stencil_r*u_k_p1 + delta_t*src_val
    
  end

  # apply BC terms to rhs
  rhs[1] = BCL(tstep*delta_t)
  rhs[mat_size] = BCR(tstep*delta_t)

  println("rhs = \n", rhs)
  # solve for next time step
#  A_ldiv_B!(Af, rhs)  # rhs gets overwritten with new solution values
  A_ldiv_B!(ksp, rhs, u_i)
#  u_i = A\rhs
#  println("secondard u_next = \n", u_next)

  println("u_next = \n", u_i)



end

# apply BCs to final time
#u_i_1[1] = BCL((nStep - 1)*delta_t)
#u_i_1[N] = BCR((nStep - 1)*delta_t)


println("time = ", time)

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

  
  A = PETSc.Mat(Float64, mat_size, mat_size, nz=stencil_size)
  u_i = Vec(Float64, mat_size)
  rhs = Vec(Float64, mat_size)
  ksp = PETSc.KSP(A, A)
  setoptions!(ksp; ksp_atol=1e-10, ksp_rtol=1e-8, ksp_monitor="")

  return A, u_i, rhs, ksp
end

# run
driver()
