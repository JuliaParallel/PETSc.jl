using PETSc, MPI, Printf

# Small ODE to test TS accuracy.
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
    ctx = unsafe_pointer_to_objref(ctx_ptr)::Ex51Context
    petsclib = ctx.petsclib
    age = petsclib.age
    u = PETSc.LibPETSc.PetscVec(u_ptr, petsclib, age)
    f = PETSc.LibPETSc.PetscVec(f_ptr, petsclib, age)

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
- `options`: additional PETSc command-line options. By default the script uses
  `ARGS`, so invocations like `julia --project=. examples/ex51.jl -ts_type rk
  -ts_rk_type 5dp` work.
- `verbose`: print a short run summary.

Returns a named tuple with the final time, error norm, and numerical solution.
"""
function solve_ex51(;
    petsclib = PETSc.getlib(PetscScalar = Float64),
    final_time::Real = 1.0,
    dt::Real = 0.25,
    save_trajectory::Bool = true,
    options = nothing,
    verbose::Bool = true,
)
    options === nothing && (options = copy(String.(ARGS)))

    comm = MPI.COMM_WORLD
    PetscInt = petsclib.PetscInt
    PetscScalar = petsclib.PetscScalar
    did_initialize = !PETSc.initialized(petsclib)

    if did_initialize
        PETSc.initialize(petsclib; options = copy(options))
    end

    MPI.Comm_size(comm) == 1 || error("This example only supports sequential runs.")

    # Keep these variables concretely typed across the whole function, even
    # though the actual PETSc objects are created later inside the `try` block.
    # The null-pointer placeholders are overwritten by `TSCreate`/`VecCreate`,
    # and the `finally` block checks `ptr != C_NULL` before destroying them.
    ts = PETSc.LibPETSc.TS(petsclib)
    u = PETSc.LibPETSc.PetscVec(petsclib)
    u_exact = PETSc.LibPETSc.PetscVec(petsclib)
    current_time = petsclib.PetscReal(NaN)
    error_norm = petsclib.PetscReal(NaN)
    solution = PetscScalar[]
    ctx = Ex51Context(petsclib)

    try
        # Create timestepping solver context.
        ts = PETSc.LibPETSc.TSCreate(petsclib, comm)
        PETSc.LibPETSc.TSSetType(petsclib, ts, "rosw")
        PETSc.LibPETSc.TSSetProblemType(petsclib, ts, PETSc.LibPETSc.TS_NONLINEAR)

        # Set initial conditions.
        u = PETSc.LibPETSc.VecCreate(petsclib, comm)
        PETSc.LibPETSc.VecSetSizes(petsclib, u, PetscInt(2), PetscInt(2))
        PETSc.LibPETSc.VecSetUp(petsclib, u)
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

            PETSc.LibPETSc.TSSetFromOptions(petsclib, ts)
            PETSc.LibPETSc.TSSolve(petsclib, ts, u)
        end

        # Compute the error against the analytical solution at the achieved
        # final time, matching the PETSc example's `VecAYPX` + `VecNorm` path.
        current_time = PETSc.LibPETSc.TSGetTime(petsclib, ts)
        u_exact = PETSc.LibPETSc.VecDuplicate(petsclib, u)
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
        if u_exact.ptr != C_NULL
            PETSc.destroy(u_exact)
        end
        if u.ptr != C_NULL
            PETSc.destroy(u)
        end
        if ts.ptr != C_NULL
            PETSc.LibPETSc.TSDestroy(petsclib, ts)
        end
        if did_initialize
            PETSc.finalize(petsclib)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    solve_ex51()
end
