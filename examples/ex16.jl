using PETSc, MPI, Printf

# Van der Pol example adapted from PETSc TS tutorial `ex16.c`, see
# https://petsc.org/main/src/ts/tutorials/ex16.c.html.
#
# The second-order ODE
#
#     y'' - mu * ((1 - y^2) * y' - y) = 0
#
# is rewritten as the first-order system
#
#     u1_t = u2
#     u2_t = mu * ((1 - u1^2) * u2 - u1)
#
# and then split into explicit and implicit pieces for IMEX timestepping:
#
#     G(u, t) = [u2, 0]
#     F(u_t, u, t) = [u1_t, u2_t - mu * ((1 - u1^2) * u2 - u1)]
#
# with the special case `imex = false` folding `u2` back into the implicit
# residual, matching the upstream example.

mutable struct Ex16Context{
    PetscLib <: PETSc.LibPETSc.PetscLibType,
    PetscReal <: Real,
}
    petsclib::PetscLib
    mu::PetscReal
    imex::Bool
    next_output::PetscReal
end

# In Julia it is common to load this file once and call `solve_ex16(...)`
# repeatedly from the same interactive session. The custom `myark2` ARKIMEX
# method therefore needs a small "register once per PETSc lifetime" guard so we
# do not try to re-register it on every solve. The upstream PETSc C example
# does not need this bookkeeping because it is a one-shot program: it starts,
# registers the method, runs once, and exits.
const EX16_REGISTERED_AGES = IdDict{Any, Int}()

# Split parsed command-line options into:
# - example-specific options handled directly in this file (`mu`, `imex`,
#   `monitor`)
# - all remaining PETSc solver options, which we still want to forward to
#   `TSSetFromOptions`
function ex16_solver_options(parsed_options::NamedTuple)
    filtered_pairs = Pair{Symbol, Union{String, Nothing}}[]
    for (key, value) in pairs(parsed_options)
        key in (:mu, :imex, :monitor) && continue
        push!(filtered_pairs, key => value)
    end
    return (; filtered_pairs...)
end

# Read the example-specific runtime options (`mu`, `imex`, `monitor`) using
# PETSc's own option-query routines so this part stays close to the upstream C
# example. This complements `ex16_solver_options(...)`: that function filters
# these keys out before the remaining options are passed on to
# `TSSetFromOptions`.
function ex16_runtime_options(
    petsclib,
    parsed_options::NamedTuple;
    mu::Real,
    imex::Bool,
    monitor::Bool,
)
    PetscReal = petsclib.PetscReal
    query_options = PETSc.Options(petsclib; parsed_options...)

    try
        mu_value, mu_set =
            PETSc.LibPETSc.PetscOptionsGetReal(petsclib, query_options, "", "-mu")
        imex_value, imex_set =
            PETSc.LibPETSc.PetscOptionsGetBool(petsclib, query_options, "", "-imex")
        monitor_value, monitor_set =
            PETSc.LibPETSc.PetscOptionsGetBool(petsclib, query_options, "", "-monitor")

        return (
            mu = Bool(mu_set) ? PetscReal(mu_value) : PetscReal(mu),
            imex = Bool(imex_set) ? Bool(imex_value) : imex,
            monitor = Bool(monitor_set) ? Bool(monitor_value) : monitor,
        )
    finally
        PETSc.destroy(query_options)
    end
end

function ex16_register_myark2!(petsclib)
    get(EX16_REGISTERED_AGES, petsclib, -1) == petsclib.age && return nothing

    # For the stage tables `At` and `A`, PETSc expects flat vectors in row-major
    # order, matching the layout used by C arrays.
    PetscReal = petsclib.PetscReal
    A = PetscReal[
        0.0, 0.0, 0.0,
        0.41421356237309504880, 0.0, 0.0,
        0.75, 0.25, 0.0,
    ]
    At = PetscReal[
        0.0, 0.0, 0.0,
        0.12132034355964257320, 0.29289321881345247560, 0.0,
        0.20710678118654752440, 0.5, 0.29289321881345247560,
    ]

    PETSc.LibPETSc.TSARKIMEXRegister(
        petsclib,
        "myark2",
        2,
        3,
        At,
        nothing,
        nothing,
        A,
        nothing,
        nothing,
        nothing,
        nothing,
        0,
        nothing,
        nothing,
    )
    EX16_REGISTERED_AGES[petsclib] = petsclib.age
    return nothing
end

"""
    ex16_rhs!

Right-hand side callback for the explicit part of the van der Pol split.
"""
function ex16_rhs!(
    ::PETSc.LibPETSc.CTS,
    ::PETSc.LibPETSc.PetscReal,
    x_ptr::PETSc.LibPETSc.CVec,
    f_ptr::PETSc.LibPETSc.CVec,
    ctx_ptr::Ptr{Cvoid},
)::PETSc.LibPETSc.PetscErrorCode
    ctx = unsafe_pointer_to_objref(ctx_ptr)::Ex16Context
    petsclib = ctx.petsclib
    x = PETSc.VecPtr(petsclib, x_ptr, false)
    f = PETSc.VecPtr(petsclib, f_ptr, false)

    PETSc.withlocalarray!(
        (x, f);
        read = (true, false),
        write = (false, true),
    ) do x_array, f_array
        f_array[1] = ctx.imex ? x_array[2] : 0.0
        f_array[2] = 0.0
    end

    return PETSc.LibPETSc.PetscErrorCode(0)
end

const EX16_RHS_FUNCTION_PTR = @cfunction(
    ex16_rhs!,
    PETSc.LibPETSc.PetscErrorCode,
    (
        PETSc.LibPETSc.CTS,
        PETSc.LibPETSc.PetscReal,
        PETSc.LibPETSc.CVec,
        PETSc.LibPETSc.CVec,
        Ptr{Cvoid},
    ),
)

"""
    ex16_ifunction!

Implicit residual callback for the van der Pol problem.
"""
function ex16_ifunction!(
    ::PETSc.LibPETSc.CTS,
    ::PETSc.LibPETSc.PetscReal,
    x_ptr::PETSc.LibPETSc.CVec,
    xdot_ptr::PETSc.LibPETSc.CVec,
    f_ptr::PETSc.LibPETSc.CVec,
    ctx_ptr::Ptr{Cvoid},
)::PETSc.LibPETSc.PetscErrorCode
    ctx = unsafe_pointer_to_objref(ctx_ptr)::Ex16Context
    petsclib = ctx.petsclib
    x = PETSc.VecPtr(petsclib, x_ptr, false)
    xdot = PETSc.VecPtr(petsclib, xdot_ptr, false)
    f = PETSc.VecPtr(petsclib, f_ptr, false)

    PETSc.withlocalarray!(
        (x, xdot, f);
        read = (true, true, false),
        write = (false, false, true),
    ) do x_array, xdot_array, f_array
        f_array[1] = xdot_array[1] + (ctx.imex ? 0.0 : x_array[2])
        f_array[2] =
            xdot_array[2] - ctx.mu * ((1.0 - x_array[1] * x_array[1]) * x_array[2] - x_array[1])
    end

    return PETSc.LibPETSc.PetscErrorCode(0)
end

const EX16_IFUNCTION_PTR = @cfunction(
    ex16_ifunction!,
    PETSc.LibPETSc.PetscErrorCode,
    (
        PETSc.LibPETSc.CTS,
        PETSc.LibPETSc.PetscReal,
        PETSc.LibPETSc.CVec,
        PETSc.LibPETSc.CVec,
        PETSc.LibPETSc.CVec,
        Ptr{Cvoid},
    ),
)

# PETSc's IMEX TS interface does not ask for the Jacobian of the full ODE
# right-hand side `u_t = f(u, t)`. It asks for the shifted Jacobian
# corresponding to the implicit residual supplied via `TSSetIFunction`, with
# the time-derivative contribution folded in through the scalar `a` passed to
# `IJacobian`. In this example
#   `G(u, t) = [u2, 0]`
#   `F(u_t, u, t) = [u1_t + (imex ? 0 : u2),
#                    u2_t - mu * ((1 - u1^2) * u2 - u1)]`
# so the 2x2 entries assembled below match that residual, not the explicit
# `RHSFunction` alone.
function ex16_fill_ijacobian!(
    A::PETSc.LibPETSc.PetscMat,
    B::PETSc.LibPETSc.PetscMat,
    a,
    x1,
    x2,
    ctx::Ex16Context,
)
    PetscScalar = ctx.petsclib.PetscScalar
    mu = ctx.mu

    A[1, 1] = PetscScalar(a)
    A[1, 2] = PetscScalar(ctx.imex ? 0.0 : 1.0)
    A[2, 1] = PetscScalar(mu * (2.0 * x1 * x2 + 1.0))
    A[2, 2] = PetscScalar(a - mu * (1.0 - x1 * x1))
    PETSc.assemble!(A)

    if B.ptr != A.ptr
        B[1, 1] = PetscScalar(a)
        B[1, 2] = PetscScalar(ctx.imex ? 0.0 : 1.0)
        B[2, 1] = PetscScalar(mu * (2.0 * x1 * x2 + 1.0))
        B[2, 2] = PetscScalar(a - mu * (1.0 - x1 * x1))
        PETSc.assemble!(B)
    end

    return nothing
end

"""
    ex16_ijacobian!

Implicit Jacobian callback corresponding to `ex16.c`.
"""
function ex16_ijacobian!(
    ::PETSc.LibPETSc.CTS,
    ::PETSc.LibPETSc.PetscReal,
    x_ptr::PETSc.LibPETSc.CVec,
    ::PETSc.LibPETSc.CVec,
    a::PETSc.LibPETSc.PetscReal,
    A_ptr::PETSc.LibPETSc.CMat,
    B_ptr::PETSc.LibPETSc.CMat,
    ctx_ptr::Ptr{Cvoid},
)::PETSc.LibPETSc.PetscErrorCode
    ctx = unsafe_pointer_to_objref(ctx_ptr)::Ex16Context
    petsclib = ctx.petsclib
    x = PETSc.VecPtr(petsclib, x_ptr, false)
    A = PETSc.LibPETSc.PetscMat(A_ptr, petsclib)
    B = PETSc.LibPETSc.PetscMat(B_ptr, petsclib)

    x1, x2 = PETSc.withlocalarray!(x; read = true, write = false) do x_array
        (x_array[1], x_array[2])
    end

    ex16_fill_ijacobian!(A, B, a, x1, x2, ctx)
    return PETSc.LibPETSc.PetscErrorCode(0)
end

const EX16_IJACOBIAN_PTR = @cfunction(
    ex16_ijacobian!,
    PETSc.LibPETSc.PetscErrorCode,
    (
        PETSc.LibPETSc.CTS,
        PETSc.LibPETSc.PetscReal,
        PETSc.LibPETSc.CVec,
        PETSc.LibPETSc.CVec,
        PETSc.LibPETSc.PetscReal,
        PETSc.LibPETSc.CMat,
        PETSc.LibPETSc.CMat,
        Ptr{Cvoid},
    ),
)

"""
    ex16_monitor!

Optional monitor callback that interpolates the TS solution at multiples of
`0.1`, mirroring the upstream example.
"""
function ex16_monitor!(
    ts_ptr::PETSc.LibPETSc.CTS,
    step::PETSc.LibPETSc.PetscInt,
    t::PETSc.LibPETSc.PetscReal,
    x_ptr::PETSc.LibPETSc.CVec,
    ctx_ptr::Ptr{Cvoid},
)::PETSc.LibPETSc.PetscErrorCode
    ctx = unsafe_pointer_to_objref(ctx_ptr)::Ex16Context
    petsclib = ctx.petsclib
    ts = PETSc.LibPETSc.TS(ts_ptr, petsclib)
    x = PETSc.VecPtr(petsclib, x_ptr, false)
    dt = PETSc.LibPETSc.TSGetTimeStep(petsclib, ts)
    tfinal = PETSc.LibPETSc.TSGetMaxTime(petsclib, ts)

    while ctx.next_output <= t && ctx.next_output <= tfinal
        interpolated_x = similar(x)
        try
            PETSc.LibPETSc.TSInterpolate(
                petsclib,
                ts,
                petsclib.PetscReal(ctx.next_output),
                interpolated_x,
            )
            PETSc.withlocalarray!(interpolated_x; read = true, write = false) do x_array
                @printf(
                    "[%.1f] %d TS %.6f (dt = %.6f) X % 12.6e % 12.6e\n",
                    ctx.next_output,
                    step,
                    t,
                    dt,
                    x_array[1],
                    x_array[2],
                )
            end
        finally
            PETSc.destroy(interpolated_x)
        end

        ctx.next_output += 0.1
    end

    return PETSc.LibPETSc.PetscErrorCode(0)
end

const EX16_MONITOR_PTR = @cfunction(
    ex16_monitor!,
    PETSc.LibPETSc.PetscErrorCode,
    (
        PETSc.LibPETSc.CTS,
        PETSc.LibPETSc.PetscInt,
        PETSc.LibPETSc.PetscReal,
        PETSc.LibPETSc.CVec,
        Ptr{Cvoid},
    ),
)

function ex16_create_jacobian_template(petsclib)
    PetscScalar = petsclib.PetscScalar
    jac = PETSc.MatSeqAIJ(petsclib, 2, 2, petsclib.PetscInt(2))
    jac[1, [1, 2]] = PetscScalar.([1.0, 1.0])
    jac[2, [1, 2]] = PetscScalar.([1.0, 1.0])
    PETSc.assemble!(jac)
    PETSc.LibPETSc.MatZeroEntries(petsclib, jac)
    return jac
end

function ex16_initial_condition!(u::PETSc.LibPETSc.PetscVec, mu::Real)
    PETSc.withlocalarray!(u; read = false, write = true) do u_array
        u_array[1] = 2.0
        u_array[2] = -2.0 / 3.0 + 10.0 / (81.0 * mu) - 292.0 / (2187.0 * mu * mu)
    end
    return nothing
end

"""
    solve_ex16(; kwargs...)

Solve the van der Pol IMEX example from PETSc TS tutorial `ex16.c`.

Keyword arguments:

- `petsclib`: PETSc library instance. Defaults to `Float64`.
- `mu`: van der Pol stiffness parameter. Default `1000.0`.
- `imex`: whether to keep the `u2` term in the explicit split. Default `true`.
- `monitor`: register the interpolation monitor from the upstream example.
  Default `false`.
- `final_time`: final integration time. Default `0.5`.
- `dt`: initial time step. Default `0.01`.
- `finalize_petsc`: whether to finalize PETSc at the end when this function had
  to initialize it. Default `false`.
- `options`: additional PETSc command-line options for this solve. By default
  the script uses `ARGS`.
- `verbose`: print a short run summary after the solve.

Returns a named tuple with the final time, step count, solution, and effective
problem parameters.
"""
function solve_ex16(;
    petsclib = PETSc.getlib(PetscScalar = Float64),
    mu::Real = 1000.0,
    imex::Bool = true,
    monitor::Bool = false,
    final_time::Real = 0.5,
    dt::Real = 0.01,
    finalize_petsc::Bool = false,
    options = nothing,
    verbose::Bool = true,
)
    options === nothing && (options = copy(String.(ARGS)))
    parsed_options = PETSc.parse_options(options)
    solver_options = ex16_solver_options(parsed_options)

    comm = MPI.COMM_WORLD
    PetscScalar = petsclib.PetscScalar

    did_initialize = !PETSc.initialized(petsclib)
    if did_initialize
        PETSc.initialize(petsclib)
    end

    MPI.Comm_size(comm) == 1 || error("This example only supports sequential runs.")

    runtime_options = ex16_runtime_options(
        petsclib,
        parsed_options;
        mu,
        imex,
        monitor,
    )

    ts = PETSc.LibPETSc.TS(petsclib)
    u = PETSc.LibPETSc.PetscVec(petsclib)
    jac = PETSc.LibPETSc.PetscMat(petsclib)
    current_time = petsclib.PetscReal(NaN)
    steps = petsclib.PetscInt(-1)
    solution = PetscScalar[]
    ctx = Ex16Context(
        petsclib,
        runtime_options.mu,
        runtime_options.imex,
        petsclib.PetscReal(0),
    )
    petsc_options = PETSc.Options(petsclib; solver_options...)
    pushed_options = false

    try
        ex16_register_myark2!(petsclib)

        ts = PETSc.LibPETSc.TSCreate(petsclib, comm)
        PETSc.LibPETSc.TSSetType(petsclib, ts, "beuler")
        PETSc.LibPETSc.TSSetProblemType(petsclib, ts, PETSc.LibPETSc.TS_NONLINEAR)

        u = PETSc.VecSeq(petsclib, 2)
        ex16_initial_condition!(u, runtime_options.mu)
        PETSc.LibPETSc.TSSetSolution(petsclib, ts, u)

        jac = ex16_create_jacobian_template(petsclib)

        callback_ctx_ptr = pointer_from_objref(ctx)
        GC.@preserve ctx begin
            PETSc.LibPETSc.TSSetRHSFunction(
                petsclib,
                ts,
                nothing,
                EX16_RHS_FUNCTION_PTR,
                callback_ctx_ptr,
            )
            PETSc.LibPETSc.TSSetIFunction(
                petsclib,
                ts,
                nothing,
                EX16_IFUNCTION_PTR,
                callback_ctx_ptr,
            )
            PETSc.LibPETSc.TSSetIJacobian(
                petsclib,
                ts,
                jac,
                jac,
                EX16_IJACOBIAN_PTR,
                callback_ctx_ptr,
            )
            PETSc.LibPETSc.TSSetMaxTime(petsclib, ts, petsclib.PetscReal(final_time))
            PETSc.LibPETSc.TSSetExactFinalTime(
                petsclib,
                ts,
                PETSc.LibPETSc.TS_EXACTFINALTIME_STEPOVER,
            )
            PETSc.LibPETSc.TSSetTimeStep(petsclib, ts, petsclib.PetscReal(dt))

            if runtime_options.monitor
                PETSc.LibPETSc.TSMonitorSet(
                    petsclib,
                    ts,
                    EX16_MONITOR_PTR,
                    callback_ctx_ptr,
                )
            end

            push!(petsc_options)
            pushed_options = true
            PETSc.LibPETSc.TSSetFromOptions(petsclib, ts)
            pop!(petsc_options)
            pushed_options = false

            PETSc.LibPETSc.TSSolve(petsclib, ts, u)
        end

        current_time = PETSc.LibPETSc.TSGetSolveTime(petsclib, ts)
        steps = PETSc.LibPETSc.TSGetStepNumber(petsclib, ts)
        solution = copy(u[:])

        if verbose
            if abs(current_time - final_time) > 100 * eps(petsclib.PetscReal)
                @printf(
                    "Note: prescribed final time %.16g differs from actual final time %.16g\n",
                    final_time,
                    current_time,
                )
            end

            ts_type = PETSc.LibPETSc.TSGetType(petsclib, ts)
            @printf("TS type: %s\n", ts_type)
            if ts_type == "arkimex"
                @printf("ARKIMEX type: %s\n", PETSc.LibPETSc.TSARKIMEXGetType(petsclib, ts))
            end
            @printf("mu: %.16g\n", runtime_options.mu)
            @printf("IMEX split: %s\n", runtime_options.imex ? "true" : "false")
            @printf("Final time: %.16g\n", current_time)
            @printf("Steps: %d\n", steps)
            @printf(
                "Final solution: [% .16e, % .16e]\n",
                solution[1],
                solution[2],
            )
        end

        return (
            final_time = current_time,
            steps = steps,
            solution = solution,
            mu = runtime_options.mu,
            imex = runtime_options.imex,
        )
    finally
        if pushed_options
            pop!(petsc_options)
        end
        if petsc_options.ptr != C_NULL
            PETSc.destroy(petsc_options)
        end
        if jac.ptr != C_NULL
            PETSc.destroy(jac)
        end
        if u.ptr != C_NULL
            PETSc.destroy(u)
        end
        if ts.ptr != C_NULL
            PETSc.LibPETSc.TSDestroy(petsclib, ts)
        end
        if did_initialize && finalize_petsc
            PETSc.finalize(petsclib)
        end
    end
end

if !isinteractive() && abspath(PROGRAM_FILE) == @__FILE__
    solve_ex16(finalize_petsc = true)
end
