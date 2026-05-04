SciMLBase.u_modified!(i::PETScTSIntegrator, val::Bool) = (i.u_modified = val)

@static if isdefined(SciMLBase, :derivative_discontinuity!)
    SciMLBase.derivative_discontinuity!(i::PETScTSIntegrator, val::Bool) =
        (i.derivative_discontinuity = val)
end

function SciMLBase.savevalues!(integ::PETScTSIntegrator, force::Bool = false)
    integ.opts.save_on || return (false, false)
    saved = false
    saved_exactly = false

    # Drain due `saveat` times first so the trajectory stays sorted in the
    # integration direction even when `save_everystep` is also enabled.
    # `saveat` stores tdir*t in a descending-sorted Vector so `last` is the
    # next requested time in tdir order and `pop!` removes it in O(1).
    while !isempty(integ.opts.saveat) &&
          last(integ.opts.saveat) <= integ.tdir * integ.t
        t_save = pop!(integ.opts.saveat) / integ.tdir
        if t_save == integ.t
            # Saveat coincides with the step endpoint: just save the current
            # state below; do not interpolate or duplicate.
            saved_exactly = true
            continue
        end
        u_interp = similar(integ.u)
        v_interp = PETSc.VecSeq(integ.petsclib, length(integ.u))
        try
            PETSc.LibPETSc.TSInterpolate(
                integ.petsclib, integ.ts,
                integ.petsclib.PetscReal(t_save), v_interp,
            )
            PETSc.withlocalarray!(v_interp; read = true, write = false) do arr
                copyto!(u_interp, reshape(arr, integ.sizeu))
            end
        finally
            PETSc.destroy(v_interp)
        end
        push!(integ.sol.t, t_save)
        push!(integ.sol.u, u_interp)
        saved = true
    end

    if (integ.opts.save_everystep || force || saved_exactly) &&
       (isempty(integ.sol.t) || last(integ.sol.t) != integ.t)
        push!(integ.sol.t, integ.t)
        push!(integ.sol.u, copy(integ.u))
        saved = true
    end

    return (saved, !isempty(integ.sol.t) && last(integ.sol.t) == integ.t)
end

function SciMLBase.terminate!(
    i::PETScTSIntegrator,
    retcode = SciMLBase.ReturnCode.Terminated,
)
    i.retcode = retcode
    i.done = true
    while !isempty(i.opts.tstops)
        pop!(i.opts.tstops)
    end
    return nothing
end

# Local re-implementations of the DiffEqBase callback lifecycle helpers.
# All dependencies are in SciMLBase, so no DiffEqBase dependency is needed.
_initialize_callbacks!(cb::SciMLBase.CallbackSet{Tuple{}, Tuple{}}, u, t, integ) = false
function _initialize_callbacks!(cb::SciMLBase.CallbackSet, u, t, integ)
    _initialize_callbacks!(u, t, integ, false,
        cb.continuous_callbacks..., cb.discrete_callbacks...)
end
function _initialize_callbacks!(u, t, integ, any_modified,
    c::SciMLBase.DECallback, cs::SciMLBase.DECallback...,
)
    c.initialize(c, u, t, integ)
    _initialize_callbacks!(u, t, integ, any_modified || integ.u_modified, cs...)
end
function _initialize_callbacks!(u, t, integ, any_modified, c::SciMLBase.DECallback)
    c.initialize(c, u, t, integ)
    any_modified || integ.u_modified
end

_finalize_callbacks!(cb::SciMLBase.CallbackSet{Tuple{}, Tuple{}}, u, t, integ) = false
function _finalize_callbacks!(cb::SciMLBase.CallbackSet, u, t, integ)
    _finalize_callbacks!(u, t, integ, false,
        cb.continuous_callbacks..., cb.discrete_callbacks...)
end
function _finalize_callbacks!(u, t, integ, any_modified,
    c::SciMLBase.DECallback, cs::SciMLBase.DECallback...,
)
    c.finalize(c, u, t, integ)
    _finalize_callbacks!(u, t, integ, any_modified || integ.u_modified, cs...)
end
function _finalize_callbacks!(u, t, integ, any_modified, c::SciMLBase.DECallback)
    c.finalize(c, u, t, integ)
    any_modified || integ.u_modified
end

# Local re-implementation of the callback dispatch chain that DiffEqBase
# previously provided. All dependencies (condition, affect!, save_positions,
# initializealg, reeval_internals_due_to_modification!, savevalues!) are part
# of SciMLBase, so no DiffEqBase dependency is needed.
@inline function _apply_discrete_callback!(
    integrator::PETScTSIntegrator, callback::SciMLBase.DiscreteCallback,
)
    saved_in_cb = false
    if callback.condition(integrator.u, integrator.t, integrator)
        _, savedexactly = SciMLBase.savevalues!(integrator)
        saved_in_cb = true
        @inbounds if callback.save_positions[1]
            savedexactly || SciMLBase.savevalues!(integrator, true)
        end
        integrator.u_modified = true
        callback.affect!(integrator)
        if integrator.u_modified
            SciMLBase.reeval_internals_due_to_modification!(
                integrator, false;
                callback_initializealg = callback.initializealg,
            )
        end
        @inbounds if callback.save_positions[2]
            SciMLBase.savevalues!(integrator, true)
            saved_in_cb = true
        end
    end
    integrator.u_modified, saved_in_cb
end

@inline function _apply_discrete_callback!(
    integrator::PETScTSIntegrator,
    discrete_modified::Bool, saved_in_cb::Bool,
    callback::SciMLBase.DiscreteCallback,
)
    bool, saved2 = _apply_discrete_callback!(integrator, callback)
    discrete_modified || bool, saved_in_cb || saved2
end

@inline function _apply_discrete_callback!(
    integrator::PETScTSIntegrator,
    callback::SciMLBase.DiscreteCallback,
    rest::SciMLBase.DiscreteCallback...,
)
    _apply_discrete_callback!(
        integrator, _apply_discrete_callback!(integrator, callback)..., rest...,
    )
end

function handle_callbacks!(integ::PETScTSIntegrator)
    cbs = integ.opts.callback
    discrete = cbs.discrete_callbacks
    saved_in_cb = false
    if !isempty(discrete)
        modified, saved_in_cb = _apply_discrete_callback!(integ, discrete...)
        modified && (integ.u_modified = true)
    end
    saved_in_cb || SciMLBase.savevalues!(integ)
    return nothing
end
