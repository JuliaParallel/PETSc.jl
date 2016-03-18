# functions for PETSc TS (time stepping) algorithms


type TS{T}
  p::C.TS{T}
  data
  function TS(p, data=nothing, first_instance=false)
    ts = new(p, data)

    if first_instance
      finalizer(ts, PetscDestroy)
    end

    return ts
  end
end


function PetscDestroy{T}(ts::TS{T})

  if !PetscFinalized(T)  && !isfinalized(vec)
    ts_ref = Ref(ts)
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


"""
  Preferred constructor: set problem type explicitly, get method from
  options database
"""
function TS{T<:Scalar}(::Type{T}, tsptype::C.TSProblemType; comm=MPI.COMM_WORLD)

  ts = Ref{C.TS{T}}()
  chk(C.TSCreate(comm, ts))
  chk(C.TSSetProblemType(ts[], tsptype))

  return TS{T}(ts[])
end


"""
  More explicit constructor: set problem type, method directly
"""
function TS{T<:Scalar}(tsptype::C.TSProblemType, tstype::C.TSType; 
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
function set_times(ts::TS{T}, t0, dt0,  nsteps::Integer, tmax)
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


function set_rhs_function{T}(ts::TS{T}, r::Vec{T} f::Function, ctx=())

  ctx_outer = (f, ctx)
  ctx_ptr = pointer_from_objref(ctx_outer)
  treal = real(T)
  # this this pre-compilable?
  fptr = cfunction(rhs_wrapper, PetscErrorCode, (C.TS{T}, treal, C.Vec{T}, C.Vec{T}, Ptr{Void}))

  chk(C.TSSetRHSFunction(ts.p, r.p, fptr, ctx_ptr))
end

"""
  Wrapper for the right hand side function.  This function is always passed
  to PETSc as the right hand side function, and calls the user supplied
  function internally.  The user supplied function must be the first
  component of the ctx tuple
"""
function rhs_wrapper{T}(ts::C.TS{T}, t, u::C.Vec{T}, F::C.Vec{T}, ctx_ptr::Ptr{Void})

  # transform into high level objects
  bigts = TS{T}(ts, first_instance=false)

  tref = Ref{C.VecType}()
  chk(C.VecGetType(u, tref))
  bigu = Vec{T, tref[]}(u, first_instance=false)

  tref2 = Ref{C.VecType}()
  chk(C.VecGetType(F, tref2))
  bigF = Vec{T, tref2[]}(F, first_instance=false)

  ctx = unsafe_pointer_to_objref(ctx_ptr)
  func = ctx[1]
  ctx_inner = ctx[2]  # the user provided ctx

  ret_status = func(bigts, t, bigu, bigF, ctx_inner)

  return PetscErrorCode(ret_status)
end

# a PETSc provided rhs function for the linear, time invarient coefficient
# matrix case

function ComputeRHSFunctionLinear(ts::TS, t, u::Vec, F::vec, ctx)

  # this is amusing, wrapping things just to unwrap them again
  C.TSComputeRHSFunctionLinear(ts.p, t, u.p, F.p, C_NULL)
end
