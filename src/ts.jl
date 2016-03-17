# functions for PETSc TS (time stepping) algorithms


type TS{T}
  p::C.TS{T}
  data
  function TS(p, data=nothing)
    return new(p, data)
  end
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
