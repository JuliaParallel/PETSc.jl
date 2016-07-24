# functions for PETSc TS (time stepping) algorithms
export TS, set_ic, set_times, solve!, set_rhs_function, set_rhs_jac
export set_lhs_function, set_lhs_jac
export ComputeRHSFunctionLinear, ComputeRHSJacobianConstant

type TS{T}
  p::C.TS{T}
  data::Array{Any, 1}  # hold the various ctx tuples
  function TS(p, data=nothing; first_instance=false)
    data_arr = Array(Any, 0)
    ts = new(p, data_arr)
    push!(ts.data, data)

    if first_instance
      finalizer(ts, PetscDestroy)
    end

    return ts
  end
end


function PetscDestroy{T}(ts::TS{T})
  if !PetscFinalized(T) && !isfinalized(ts)
    ts_ref = Ref(ts.p)
    chk(C.TSDestroy(ts_ref))
    ts.p = C.TS{T}(C_NULL)
  end
end

function isfinalized(ts::TS)
  return isfinalized(ts.p)
end

function isfinalized(ts::C.TS)
  return ts.pobj == C_NULL
end

#get the internal KSP object for this TS
function KSP{T}(ts::TS{T})
  ksp_c = Ref{C.KSP{T}}()
  chk(C.TSGetKSP(ts.p, ksp_c))
  return KSP{T}(ksp_c[])
end

"""
  Most preferred constructor: take ProblemType, method from options
  database
"""
function TS{T<:Scalar}(::Type{T} ;comm=MPI.COMM_WORLD)
  ts = Ref{C.TS{T}}()
  chk(C.TSCreate(comm, ts))
  return TS{T}(ts[])
end

"""
  Preferred constructor: set problem type explicitly, get method from
  options database
"""
function TS{T<:Scalar}(::Type{T}, tsptype::C.TSProblemType; comm=MPI.COMM_WORLD)

  ts = TS(T, comm=comm)
  chk(C.TSSetProblemType(ts.p, tsptype))

  return ts
end


"""
  More explicit constructor: set problem type, method directly

`tsptype` sets the problem type and can be one of the following:
* `TS_LINEAR` - a linear set of ODEs 
* `TS_NONLINEAR` - a nonlinear set of ODEs or DEAs

`tstype` sets the method used to solve the problem. More information
about the possible methods is available at the official PETSc [docs](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSType.html).
"""
function TS{T<:Scalar}(::Type{T}, tsptype::C.TSProblemType, tstype::C.TSType; 
                       comm=MPI.COMM_WORLD)

  ts = TS(T, tsptype, comm=comm)
  chk(C.TSSetType(ts.p, tstype))

  return ts
end


function set_ic{T<:Scalar}(ts::TS{T}, u::Vec{T}) 
  chk(C.TSSetSolution(ts.p, u.p))
end

"""
  Set the times related quantities:
    t0 : initial time value
    dt0: initial time step
    nsteps: maximum number of steps
    tmax: maximum time value
"""
function set_times{T<:Scalar}(ts::TS{T}, t0, dt0,  nsteps::Integer, tmax)
  TR = real(T)  # PetscReal

  chk(C.TSSetInitialTimeStep(ts.p, TR(t0), TR(dt0)))
  chk(C.TSSetDuration(ts.p, nsteps, TR(tmax)))
end

function petscview{T}(ts::TS{T})
  viewer = C.PetscViewer{T}(C_NULL)
  chk(C.TSView(ts.p, viewer))
end

"""
  Solve the system using the initial condition provided in vec
"""
function solve!{T}(ts::TS{T}, vec::Vec{T})

  chk(C.TSSolve(ts.p, vec.p))
end

"""
  Solve the system using the intitial condition proived by set_ic
"""
function solve!{T}(ts::TS{T})

  vecp = C.Vec{T}(C_NULL)
  chk(C.TSSolve(ts.p, vecp))
end


###############################################################################
# right hand side function

"""
  Sets the function that evalutes u_t = g(u, t) for an ODE.
  The function must have the signature:

  f(TS, t, U, F, ctx)

  where TS is a TS object,
  t is the current time
  u is the current state vector
  F is the vector to be populated with u_t
  ctx is the user supplied context tuple (empty tuple if not provided)
"""
function set_rhs_function{T}(ts::TS{T}, r::Vec{T}, f::Function, ctx=())

  ctx_outer = (f, ctx)
  push!(ts.data, ctx_outer)
  ctx_ptr = pointer_from_objref(ctx_outer)
  Treal = real(T)
  # this this pre-compilable?
  fptr = cfunction(rhs_wrapper, PetscErrorCode, (C.TS{T}, Treal, C.Vec{T}, C.Vec{T}, Ptr{Void}))

  chk(C.TSSetRHSFunction(ts.p, r.p, fptr, ctx_ptr))
end

"""
  Wrapper for the right hand side function.  This function is always passed
  to PETSc as the right hand side function, and calls the user supplied
  function internally.  The user supplied function must be the first
  component of the ctx tuple

"""
function rhs_wrapper{T}(ts::C.TS{T}, t, u::C.Vec{T}, F::C.Vec{T}, ctx_ptr::Ptr{Void})

  Treal = real(T)
  # transform into high level objects
  bigts = TS{T}(ts, first_instance=false)

  tref = Array(C.VecType, 1)
#  tref = Ref{C.VecType}()
#  chk(C.VecGetType(u, tref))
  bigu = Vec{T}(u, first_instance=false)

#  tref2 = Array(C.VecType, 1)
#  tref2 = Ref{C.VecType}()
#  chk(C.VecGetType(F, tref2))
  bigF = Vec{T}(F, first_instance=false)

  ctx = unsafe_pointer_to_objref(ctx_ptr)
  func = ctx[1]
  ctx_inner = ctx[2]  # the user provided ctx

  ret_status = func(bigts, Treal(t), bigu, bigF, ctx_inner)
  return PetscErrorCode(ret_status)
end

# a PETSc provided rhs function for the linear, time invarient coefficient
# matrix case

function ComputeRHSFunctionLinear(ts::TS, t, u::Vec, F::Vec, ctx::Tuple)

  # this is amusing, wrapping things just to unwrap them again
  C.TSComputeRHSFunctionLinear(ts.p, t, u.p, F.p, C_NULL)
end

###############################################################################
# right hand side jacobian

function set_rhs_jac{T}(ts::TS{T}, A::Mat{T}, B::Mat{T}, f::Function, ctx=())

  ctx_outer = (f, ctx)
  push!(ts.data, ctx_outer)
  ctx_ptr = pointer_from_objref(ctx_outer)
  Treal = real(T)
  # this this pre-compilable?
  fptr = cfunction(rhs_jac_wrapper, PetscErrorCode, (C.TS{T}, Treal, C.Vec{T}, C.Mat{T}, C.Mat{T}, Ptr{Void}))

  chk(C.TSSetRHSJacobian(ts.p, A.p, B.p, fptr, ctx_ptr))
end



function rhs_jac_wrapper{T}(ts::C.TS{T}, t, u::C.Vec{T}, A::C.Mat{T}, B::C.Mat{T}, ctx_ptr::Ptr{Void})

  Treal = real(T) 
  bigts = TS{T}(ts, first_instance=false)
#  tref = Ref{C.VecType}()
#  chk(C.VecGetType(u, tref))
  bigu = Vec{T}(u, first_instance=false)

#  tref2 = Ref{C.MatType}()
#  chk(C.MatGetType(A, tref2))
  bigA = Mat{T}(A, first_instance=false)

  bigB = Mat{T}(B, first_instance=false)

  ctx = unsafe_pointer_to_objref(ctx_ptr)
  func = ctx[1]
  ctx_inner = ctx[2]

  ret_status = func(bigts, Treal(t), bigu, bigA, bigB, ctx_inner)

  return PetscErrorCode(ret_status)
end

# a PETSc provided rhs jacobian function for time invarient jacobians
function ComputeRHSJacobianConstant(ts::TS, t, u::Vec, A::Mat, B::Mat, ctx::Tuple)

  chk(C.TSComputeRHSJacobianConstant(ts.p, t, u.p, A.p, B.p, C_NULL))
end

###############################################################################
# left hand side function

function set_lhs_function{T}(ts::TS{T}, res::Vec{T}, f::Function, ctx::Tuple=())

  ctx_outer = (f, ctx)
  push!(ts.data, ctx_outer)
  ctx_ptr = pointer_from_objref(ctx_outer)
  Treal = real(T)

  fptr = cfunction(lhs_wrapper, PetscErrorCode, (C.TS{T}, Treal, C.Vec{T}, C.Vec{T}, C.Vec{T}, Ptr{Void}))

  chk(C.TSSetIFunction(ts.p, res.p, fptr, ctx_ptr))
  return nothing
end

function lhs_wrapper{T}(ts::C.TS{T}, t, u::C.Vec{T}, ut::C.Vec{T}, F::C.Vec{T}, ctx_ptr::Ptr{Void})

  Treal = real(T)
  # transform into high level objects
  bigts = TS{T}(ts, first_instance=false)

#  tref = Array(C.VecType, 1)
#  tref = Ref{C.VecType}()
#  chk(C.VecGetType(u, tref))
  bigu = Vec{T}(u, first_instance=false)

#  tref2 = Array(C.VecType, 1)
#  chk(C.VecGetType(ut, tref2))
  bigut = Vec{T}(ut, first_instance=false)

#  tref3 = Array(C.VecType, 1)
#  tref2 = Ref{C.VecType}()
#  chk(C.VecGetType(F, tref3))
  bigF = Vec{T}(F, first_instance=false)

  ctx = unsafe_pointer_to_objref(ctx_ptr)
  func = ctx[1]
  ctx_inner = ctx[2]  # the user provided ctx

  ret_status = func(bigts, Treal(t), bigu, bigut, bigF, ctx_inner)
  return PetscErrorCode(ret_status)
end

###############################################################################
# left hand side jacobian

function set_lhs_jac{T}(ts::TS{T}, A::Mat{T}, B::Mat{T}, f::Function, ctx::Tuple=())

  ctx_outer = (f, ctx)
  push!(ts.data, ctx_outer)
  ctx_ptr = pointer_from_objref(ctx_outer)
  Treal = real(T)

  fptr = cfunction(lhs_jac_wrapper, PetscErrorCode, (C.TS{T}, Treal, C.Vec{T}, C.Vec{T}, Treal, C.Mat{T}, C.Mat{T}, Ptr{Void}))

  chk(C.TSSetIJacobian(ts.p, A.p, B.p, fptr, ctx_ptr))
end


function lhs_jac_wrapper{T}(ts::C.TS{T}, t, u::C.Vec{T}, ut::C.Vec{T}, a, A::C.Mat{T}, B::C.Mat{T}, ctx_ptr::Ptr{Void})

  Treal = real(T) 
  bigts = TS{T}(ts, first_instance=false)

#  tref = Array(C.VecType, 1)
#  chk(C.VecGetType(u, tref))
  bigu = Vec{T}(u, first_instance=false)

#  tref2 = Array(C.VecType, 1)
#  chk(C.VecGetType(ut, tref2))
  bigut = Vec{T}(ut, first_instance=false)


#  tref3 = Array(C.VecType, 1)
#  chk(C.MatGetType(A, tref3))
  bigA = Mat{T}(A, first_instance=false)

#  tref4 = Array(C.VecType, 1)
#  chk(C.MatGetType(B, tref4))
  bigB = Mat{T}(B, first_instance=false)

  ctx = unsafe_pointer_to_objref(ctx_ptr)
  func = ctx[1]
  ctx_inner = ctx[2]  # the user provided ctx

  ret_status = func(bigts, bigu, bigut, Treal(a), bigA, bigB, ctx_inner)
  return PetscErrorCode(ret_status)
end

