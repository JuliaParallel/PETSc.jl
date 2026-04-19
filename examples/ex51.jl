using PETSc, MPI, Printf

# Small ODE to test TS accuracy, see
# https://petsc.org/main/src/ts/tutorials/ex51.c.html.
#
# The ODE
#                 u1_t = cos(t),
#                 u2_t = sin(u2)
# with analytical solution
#                 u1(t) = sin(t),
#                 u2(t) = 2 * atan(exp(t) * tan(0.5))
# is used to test the accuracy of TS schemes.

# PETSc callbacks take a raw `void*` context pointer. We keep the active
# `petsclib` inside a mutable Julia object so we can pass a stable reference
# through `pointer_from_objref()` and recover it inside the RHS routine
# to access the library-specific types, functions, and constants.
mutable struct Ex51Context{PetscLib <: PETSc.LibPETSc.PetscLibType}
    petsclib::PetscLib
end

"""
    ex51_rhs!

Right-hand side routine passed to `TSSetRHSFunction`, corresponding to
`RHSFunction` in PETSc's `ex51.c`.
"""
function ex51_rhs!(
    ::PETSc.LibPETSc.CTS, # optional argument `TS ts` not needed here
    t::Float64,
    u_ptr::PETSc.LibPETSc.CVec,
    f_ptr::PETSc.LibPETSc.CVec,
    ctx_ptr::Ptr{Cvoid},
)::PETSc.LibPETSc.PetscErrorCode
    # This function is called through PETSc's C callback interface, so its argument
    # types must match the low-level PETSc ABI. In particular, the state and RHS
    # vectors arrive as raw `CVec` pointers. We immediately wrap them using
    # PETSc.jl's higher-level `VecPtr(..., own = false)` helper so we can use
    # Julia-friendly helpers such as `withlocalarray!` without taking ownership
    # away from PETSc.
    ctx = unsafe_pointer_to_objref(ctx_ptr)::Ex51Context
    petsclib = ctx.petsclib
    # `own = false` since memory is managed by PETSc internally
    u = PETSc.VecPtr(petsclib, u_ptr, false)
    f = PETSc.VecPtr(petsclib, f_ptr, false)

    PETSc.withlocalarray!(
        (u, f);
        read = (true, false),
        write = (false, true),
    ) do u_array, f_array
        f_array[1] = cos(t)
        f_array[2] = sin(u_array[2])
    end

    return PETSc.LibPETSc.PetscErrorCode(0)
end

const EX51_RHS_FUNCTION_PTR = @cfunction(
    ex51_rhs!,
    PETSc.LibPETSc.PetscErrorCode,
    (
        PETSc.LibPETSc.CTS,
        Float64,
        PETSc.LibPETSc.CVec,
        PETSc.LibPETSc.CVec,
        Ptr{Cvoid},
    ),
)

function set_initial_condition!(u::PETSc.LibPETSc.PetscVec)
    PETSc.withlocalarray!(u; read = false, write = true) do u_array
        u_array[1] = 0.0
        u_array[2] = 1.0
    end
    return nothing
end

function exact_solution!(u::PETSc.LibPETSc.PetscVec, t::Real)
    PETSc.withlocalarray!(u; read = false, write = true) do u_array
        u_array[1] = sin(t)
        u_array[2] = 2 * atan(exp(t) * tan(0.5))
    end
    return nothing
end

"""
    solve_ex51(; kwargs...)

Solve the small ODE from PETSc TS tutorial `ex51.c` using PETSc.jl's low-level
TS wrappers.

Keyword arguments:

- `petsclib`: PETSc library instance. Defaults to `Float64`.
- `final_time`: final integration time. Default `1.0`.
- `dt`: initial time step. Default `0.25`.
- `save_trajectory`: whether to call `TSSetSaveTrajectory`. Default `true`.
- `finalize_petsc`: whether to finalize PETSc at the end when this function had
  to initialize it. Default `false` so repeated calls in one Julia session work
  reliably.
- `options`: additional PETSc command-line options for this solve. By default
  the script uses `ARGS`, so invocations like
  `julia --project=. examples/ex51.jl -ts_type rk -ts_rk_type 5dp` work, and
  programmatic calls can use values such as
  `options = ["-ts_type", "rk", "-ts_rk_type", "3bs"]`.
- `verbose`: print a short run summary.

Returns a named tuple with the final time, error norm, and numerical solution.
"""
function solve_ex51(;
    petsclib = PETSc.getlib(PetscScalar = Float64),
    final_time::Real = 1.0,
    dt::Real = 0.25,
    save_trajectory::Bool = true,
    finalize_petsc::Bool = false,
    options = nothing,
    verbose::Bool = true,
)
    options === nothing && (options = copy(String.(ARGS)))
    parsed_options = PETSc.parse_options(options)

    comm = MPI.COMM_WORLD
    PetscScalar = petsclib.PetscScalar

    # `did_initialize`: whether we initialized the library in this call
    did_initialize = !PETSc.initialized(petsclib)
    if did_initialize
        PETSc.initialize(petsclib)
    end

    MPI.Comm_size(comm) == 1 || error("This example only supports sequential runs.")

    # Keep these variables concretely typed across the whole function, even
    # though the actual PETSc objects are created later inside the `try` block.
    # The null-pointer placeholders are overwritten by `TSCreate`/`VecSeq`,
    # and the `finally` block checks `ptr != C_NULL` before destroying them.
    ts = PETSc.LibPETSc.TS(petsclib)
    u = PETSc.LibPETSc.PetscVec(petsclib)
    u_exact = PETSc.LibPETSc.PetscVec(petsclib)
    current_time = petsclib.PetscReal(NaN)
    error_norm = petsclib.PetscReal(NaN)
    solution = PetscScalar[]
    ctx = Ex51Context(petsclib)
    petsc_options = PETSc.Options(petsclib; parsed_options...)
    pushed_options = false

    try
        # Create timestepping solver context.
        ts = PETSc.LibPETSc.TSCreate(petsclib, comm)
        PETSc.LibPETSc.TSSetType(petsclib, ts, "rosw")
        PETSc.LibPETSc.TSSetProblemType(petsclib, ts, PETSc.LibPETSc.TS_NONLINEAR)

        # Set initial conditions.
        u = PETSc.VecSeq(petsclib, 2)
        set_initial_condition!(u)
        PETSc.LibPETSc.TSSetSolution(petsclib, ts, u)

        callback_ctx_ptr = pointer_from_objref(ctx)
        GC.@preserve ctx begin
            PETSc.LibPETSc.TSSetRHSFunction(
                petsclib,
                ts,
                nothing,
                EX51_RHS_FUNCTION_PTR,
                callback_ctx_ptr,
            )

            # Configure solver options. As in the C example, adaptivity is
            # forced to take constant time steps.
            save_trajectory && PETSc.LibPETSc.TSSetSaveTrajectory(petsclib, ts)
            PETSc.LibPETSc.TSSetMaxTime(petsclib, ts, petsclib.PetscReal(final_time))
            PETSc.LibPETSc.TSSetExactFinalTime(
                petsclib,
                ts,
                PETSc.LibPETSc.TS_EXACTFINALTIME_STEPOVER,
            )
            PETSc.LibPETSc.TSSetTimeStep(petsclib, ts, petsclib.PetscReal(dt))

            adapt = PETSc.LibPETSc.TSGetAdapt(petsclib, ts)
            PETSc.LibPETSc.TSAdaptSetType(petsclib, adapt, "none")

            push!(petsc_options)
            pushed_options = true
            PETSc.LibPETSc.TSSetFromOptions(petsclib, ts)
            pop!(petsc_options)
            pushed_options = false
            PETSc.LibPETSc.TSSolve(petsclib, ts, u)
        end

        # Compute the error against the analytical solution at the achieved
        # final time, matching the PETSc example's `VecAYPX` + `VecNorm` path.
        current_time = PETSc.LibPETSc.TSGetTime(petsclib, ts)
        u_exact = similar(u)
        exact_solution!(u_exact, current_time)

        PETSc.LibPETSc.VecAYPX(petsclib, u_exact, PetscScalar(-1), u)
        error_norm = PETSc.LibPETSc.VecNorm(petsclib, u_exact, PETSc.NORM_2)
        solution = copy(u[:])

        if verbose
            if abs(current_time - final_time) > 100 * eps(petsclib.PetscReal)
                @printf(
                    "Note: prescribed final time %.16g differs from actual final time %.16g\n",
                    final_time,
                    current_time,
                )
            end

            @printf("TS type: %s\n", PETSc.LibPETSc.TSGetType(petsclib, ts))
            @printf("Final time: %.16g\n", current_time)
            @printf("Error at final time: %.2E\n", error_norm)
        end

        return (final_time = current_time, error = error_norm, solution = solution)
    finally
        if pushed_options
            pop!(petsc_options)
        end
        if petsc_options.ptr != C_NULL
            PETSc.destroy(petsc_options)
        end
        if u_exact.ptr != C_NULL
            PETSc.destroy(u_exact)
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
    solve_ex51(finalize_petsc = true)
end
