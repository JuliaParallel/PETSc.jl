
using PETSc
using BaseTestNext
# test timestepping (TS) functions
ST = Float64


@testset "TS" begin
  # solve the equation f =  [u1;u2] = [t + 1; 2t + 1] using forward Euler
  # This should be an exact case

  # the jacobian is always the identity for a case where u = f(t) only
  A = Mat(ST, 2, 2)
  A[1,1] = 1.0
  A[2,2] = 1.0

  # evaluate rhs u_t
  function rhs1(ts::TS, t, u::Vec, f::Vec, ctx)

    f_local = LocalArray(f)

    f_local[1] = 1.0
    f_local[2] = 2.0
    LocalArrayRestore(f_local)

    return 0
  end

  ts = TS(ST, PETSc.C.TS_LINEAR, PETSc.C.TSEULER)

  set_rhs_function(ts, NullVec[ST], rhs1)
  set_rhs_jac(ts, A, A, ComputeRHSJacobianConstant)

  set_times(ts, 0.0, 0.1, 10, 50.0)  # limited by number of timesteps

  # hold IC/final solution
  u = Vec(ST, 2)
  u[1] = 1.0
  u[2] = 1.0

  solve!(ts, u)

  # same as above, but using set_ic
  ts = TS(ST, PETSc.C.TS_LINEAR, PETSc.C.TSEULER)

  set_rhs_function(ts, NullVec[ST], rhs1)
  set_rhs_jac(ts, A, A, ComputeRHSJacobianConstant)

  set_times(ts, 0.0, 0.1, 10, 50.0)  # limited by number of timesteps

  u = Vec(ST, 2)
  u[1] = 1.0
  u[2] = 1.0
  set_ic(ts, u)
  solve!(ts)

#  petscview(u)
  @test u[1] ≈ ST(2.0)
  @test u[2] ≈ ST(3.0)

  # solve the 1D heat equation u_t = u_xx with 2nd order space, linear time
  # homogneous Dirchlet BCs

  dt = 0.05
  nsteps = 10
  dx = 0.2  # actually dx^2
  nx = 5  # number of points in domain (including boundaries
  r = dt/dx

  A_julia = zeros(ST, nx, nx)
  A_julia[1,1] = 1.0
  A_julia[nx,nx] = 1.0

  # interior discretization
  for i=2:(nx-1)
    A_julia[i, i-1] = r
    A_julia[i, i] = 1 - 2*r
    A_julia[i, i+1] = r
  end

  u_julia = ones(ST, nx)
  u_julia[1] = 0.0  # homogenous direchlet
  u_julia[nx] = 0.0

  # solve in julia
  for t = 1:(nsteps+0)
    u_next = A_julia*u_julia
    u_julia[:] = u_next
  end

  # now solve in Petsc

  function rhs_jacobian(ts::TS, t, u, A::Mat, B::Mat, ctx)
    dx = ctx[1]
    nx, ny = size(A)
    A[1,1] = 1.0
    A[nx, nx] = 1.0

    for i=2:(nx-1)
      A[i, i-1] = 1/dx
      A[i, i] = -2/dx
      A[i, i+1] = 1/dx
    end

    return 0
  end

  function rhs_func(ts::TS, t, u::Vec, F::Vec, ctx)
    dx = ctx[1]
    nx = length(u)
    u_arr = LocalArrayRead(u)
    f_arr = LocalArray(F)

    f_arr[1] = 0.0
    f_arr[nx] = 0.0

    for i=2:(nx-1)
      f_arr[i] = u_arr[i-1]/dx + -2*u_arr[i]/dx + u_arr[i+1]/dx
    end

    LocalArrayRestore(u_arr)
    LocalArrayRestore(f_arr)

    return 0
  end

  ts = TS(ST, PETSc.C.TS_LINEAR, PETSc.C.TSEULER)

  A = Mat(ST, nx, nx)
  u = Vec(ST, nx)

  # set IC
  for i=2:(nx-1)
    u[i] = 1.0
  end

  ctx1 = (dx,)
  ctx2 = (dx,)
  set_rhs_function(ts, NullVec[ST], rhs_func, ctx1)
  set_rhs_jac(ts, A, A, rhs_jacobian, ctx2)

  set_times(ts, 0.0, dt, nsteps , 50.0)  # limited by number of timesteps

  solve!(ts, u)

  for i=1:nx
    @test u[i] ≈ u_julia[i]
  end

  PETSc.PetscDestroy(ts)
  @test PETSc.isfinalized(ts)


  # now solve the same problem using implicit euler
  dt = 0.05
  nsteps = 10
  dx = 0.2  # actually dx^2
  nx = 5  # number of points in domain (including boundaries
  r = dt/dx

  A_julia = zeros(ST, nx, nx)
  A_julia[1,1] = 1.0
  A_julia[nx,nx] = 1.0

  # interior discretization
  for i=2:(nx-1)
    A_julia[i, i-1] = -r
    A_julia[i, i] = 1 + 2*r
    A_julia[i, i+1] = -r
  end

  u_julia = ones(ST, nx)
  u_julia[1] = 0.0  # homogenous direchlet
  u_julia[nx] = 0.0

  # solve in julia
  for t = 1:(nsteps+0)
    u_next = A_julia\u_julia
    u_julia[:] = u_next
  end


  function lhs_func(ts::TS, t, u::Vec, ut::Vec, F::Vec, ctx)
    dx = ctx[1]
    nx = length(u)
    u_arr = LocalArrayRead(u)
    ut_arr = LocalArrayRead(ut)
    f_arr = LocalArray(F)

    f_arr[1] = 0.0
    f_arr[nx] = 0.0

    for i=2:(nx-1)
      f_arr[i] = ut[i] - (u_arr[i-1]/dx + -2*u_arr[i]/dx + u_arr[i+1]/dx)
    end

    LocalArrayRestore(u_arr)
    LocalArrayRestore(ut_arr)
    LocalArrayRestore(f_arr)

    return 0
  end

  function lhs_jacobian(ts::TS, u::Vec, ut::Vec, t, A::Mat, B::Mat, ctx)
    dx = ctx[1]
    nx, ny = size(A)
    A[1,1] = 1.0
    A[nx, nx] = 1.0

    for i=2:(nx-1)
      A[i, i-1] = -1/dx
      A[i, i] = 2/dx
      A[i, i+1] = -1/dx
    end

    AssemblyBegin(A, PETSc.C.MAT_FINAL_ASSEMBLY)
    AssemblyEnd(A, PETSc.C.MAT_FINAL_ASSEMBLY)
    println("A.p = ", A.p)
    println("B.p = ", B.p)

    return 0
  end

  ts = TS(ST, PETSc.C.TS_LINEAR, PETSc.C.TSBEULER)

  A = Mat(ST, nx, nx, nz = 3)
  A[1,1] = 1.0
  A[nx, nx] = 1.0
  for i=2:(nx-1)
    A[i, i-1] = 1.0
    A[i, i] = 1.0
    A[i, i+1] = 1.0
  end
  AssemblyBegin(A, PETSc.C.MAT_FINAL_ASSEMBLY)
  AssemblyEnd(A, PETSc.C.MAT_FINAL_ASSEMBLY)
  u = Vec(ST, nx)

  # set IC
  for i=2:(nx-1)
    u[i] = 1.0
  end

  ctx1 = (dx,)
  ctx2 = (dx,)
  set_lhs_function(ts, NullVec[ST], lhs_func, ctx1)
  set_lhs_jac(ts, A, A, lhs_jacobian, ctx2)

  set_times(ts, 0.0, dt, nsteps , 50.0)  # limited by number of timesteps
  AssemblyBegin(u)
  AssemblyEnd(u)
  println("isfinalized(ts) = ", PETSc.isfinalized(ts))
  solve!(ts, u)

  for i=1:nx
    println("i = ", i, ", u_petsc = ", u[i], ", u_julia = ", u_julia[i])
    @test u[i] ≈ u_julia[i]
  end

  PETSc.PetscDestroy(ts)


end
