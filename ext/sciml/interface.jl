DiffEqBase.u_modified!(i::PETScTSIntegrator, val::Bool) = (i.u_modified = val)

@static if isdefined(SciMLBase, :derivative_discontinuity!)
    SciMLBase.derivative_discontinuity!(i::PETScTSIntegrator, val::Bool) =
        (i.derivative_discontinuity = val)
end

function DiffEqBase.savevalues!(integ::PETScTSIntegrator, force::Bool = false)
    integ.opts.save_on || return (false, false)
    saved = false
    saved_exactly = false

    # Drain due `saveat` times first so the trajectory stays sorted in the
    # integration direction even when `save_everystep` is also enabled.
    # `saveat` stores tdir*t in a forward BinaryMinHeap, so first(...) is the
    # next requested time in tdir order.
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

function DiffEqBase.terminate!(
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

function handle_callbacks!(integ::PETScTSIntegrator)
    cbs = integ.opts.callback
    discrete = cbs.discrete_callbacks
    saved_in_cb = false
    if !isempty(discrete)
        modified, saved_in_cb =
            DiffEqBase.apply_discrete_callback!(integ, discrete...)
        modified && (integ.u_modified = true)
    end
    saved_in_cb || DiffEqBase.savevalues!(integ)
    return nothing
end
