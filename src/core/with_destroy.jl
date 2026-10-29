# ============================================================================
#   Destroying a group of objects together
# ============================================================================

# The objects `with_destroy` destroys when its block ends, in the order tracked
struct DestroyList
    objects::Vector{Any}
end

"""
    with_destroy(f)

Run `f(list)` and return its result, then call `destroy!` on every object
tracked in `list` with [`track!`](@ref), in reverse order of tracking. This
happens also when `f` throws, and the exception is rethrown afterwards.

On more than one process `destroy!` is collective and no finalizer frees an
object, so every rank has to destroy the same objects in the same order. Running
the same block on every rank does that.

```julia
iterations = PETSc.with_destroy() do list
    dm  = PETSc.track!(list, PETSc.DMStag(petsclib, comm, ...))
    ksp = PETSc.track!(list, PETSc.KSP(dm; ksp_type = "cg"))
    x   = PETSc.track!(list, PETSc.global_vec(dm))
    b   = PETSc.track!(list, PETSc.global_vec(dm))
    ...
    PETSc.solve!(x, ksp, b)
    PETSc.iteration_number(ksp)
end
```

An object that leaves the block, as part of the result or stored elsewhere,
must be taken out with [`untrack!`](@ref) first, or it is destroyed under its
new owner.

The cleanup after an exception is safe when every rank throws, as they do for a
[`ConvergenceError`](@ref), whose reason is the same on every rank. An exception
on some ranks only leaves the ranks at different points of the program, which
`with_destroy` cannot repair: the other ranks will wait in their next collective
call.
"""
function with_destroy(f)
    list = DestroyList(Any[])
    try
        return f(list)
    finally
        for obj in Iterators.reverse(list.objects)
            destroy!(obj)
        end
        empty!(list.objects)
    end
end

"""
    track!(list, obj)

Have `obj` destroyed when the [`with_destroy`](@ref) block that gave `list`
ends, and return `obj`. A borrowed handle (see [`owns`](@ref)) is returned
without being tracked, since `destroy!` on it does nothing.
"""
function track!(list::DestroyList, obj)
    applicable(destroy!, obj) ||
        throw(ArgumentError("`destroy!` has no method for a $(typeof(obj)), so it cannot be tracked"))
    any(o -> o === obj, list.objects) &&
        throw(ArgumentError("this $(nameof(typeof(obj))) is already tracked"))
    hasproperty(obj, :own) && !owns(obj) && return obj
    push!(list.objects, obj)
    return obj
end

"""
    untrack!(list, obj)

Stop tracking `obj` in `list`, so that the [`with_destroy`](@ref) block leaves
it alive, and return `obj`. Destroying it is then the caller's job. A borrowed
handle, which [`track!`](@ref) never tracks, is returned as it is.
"""
function untrack!(list::DestroyList, obj)
    hasproperty(obj, :own) && !owns(obj) && return obj
    i = findfirst(o -> o === obj, list.objects)
    isnothing(i) && throw(ArgumentError("this $(nameof(typeof(obj))) is not tracked"))
    deleteat!(list.objects, i)
    return obj
end
