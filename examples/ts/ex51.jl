using PETSc, MPI, Printf
using LinearAlgebra: norm

# Small ODE to test TS accuracy, see
# https://petsc.org/main/src/ts/tutorials/ex51.c.html.
#
# The ODE
#                 u1_t = cos(t),
#                 u2_t = sin(u2)
# with analytical solution
#                 u1(t) = sin(t),
#                 u2(t) = 2 * atan(exp(t) * tan(0.5))
# is a test of the accuracy of TS schemes.

"""
    ex51_rhs!(F, ts, t, u)

Right-hand side of the ODE, corresponding to `RHSFunction` in PETSc's `ex51.c`.
It fills the vector `F`.
"""
function ex51_rhs!(F, ts, t, u)
    PETSc.with_local_array!(
        (u, F);
        read = (true, false),
        write = (false, true),
    ) do u_array, f_array
        f_array[1] = cos(t)
        f_array[2] = sin(u_array[2])
    end
    return nothing
end

function set_initial_condition!(u::PETSc.LibPETSc.PetscVec)
    PETSc.with_local_array!(u; read = false, write = true) do u_array
        u_array[1] = 0.0
        u_array[2] = 1.0
    end
    return nothing
end

function exact_solution!(u::PETSc.LibPETSc.PetscVec, t::Real)
    PETSc.with_local_array!(u; read = false, write = true) do u_array
        u_array[1] = sin(t)
        u_array[2] = 2 * atan(exp(t) * tan(0.5))
    end
    return nothing
end

"""
    solve_ex51(; kwargs...)

Solve the small ODE from PETSc TS tutorial `ex51.c` with PETSc.jl's TS interface.

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
  `julia --project=. examples/ts/ex51.jl -ts_type rk -ts_rk_type 5dp` work, and
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
    did_initialize = !PETSc.isinitialized(petsclib)
    if did_initialize
        PETSc.initialize(petsclib)
    end

    MPI.Comm_size(comm) == 1 || error("This example only supports sequential runs.")

    try
        return PETSc.with_destroy() do list
            # The options are applied at the solve, after what is set below, so
            # `-ts_type` and the like on the command line win.
            ts = PETSc.track!(list, PETSc.TS(petsclib, comm;
                exact_final_time = PETSc.LibPETSc.TS_EXACTFINALTIME_STEPOVER,
                parsed_options...))
            PETSc.set_type!(ts, :rosw)
            PETSc.set_problem_type!(ts, PETSc.LibPETSc.TS_NONLINEAR)

            u = PETSc.track!(list, PETSc.PetscVec(petsclib, 2))
            set_initial_condition!(u)
            PETSc.set_solution!(ts, u)
            PETSc.set_rhs_function!(ex51_rhs!, ts)

            # As in the C example, adaptivity is forced to take constant time steps.
            save_trajectory && PETSc.LibPETSc.TSSetSaveTrajectory(petsclib, ts)
            PETSc.set_max_time!(ts, final_time)
            PETSc.set_timestep!(ts, dt)
            PETSc.set_adapt_type!(ts, :none)

            PETSc.solve!(u, ts)

            # Error against the analytical solution at the achieved final time
            current_time = PETSc.current_time(ts)
            u_exact = PETSc.track!(list, similar(u))
            exact_solution!(u_exact, current_time)
            err = PETSc.track!(list, u - u_exact)
            error_norm = norm(err)
            solution = copy(u[:])

            if verbose
                if abs(current_time - final_time) > 100 * eps(petsclib.PetscReal)
                    @printf(
                        "Note: prescribed final time %.16g differs from actual final time %.16g\n",
                        final_time,
                        current_time,
                    )
                end

                @printf("TS type: %s\n", PETSc.type_name(ts))
                @printf("Final time: %.16g\n", current_time)
                @printf("Error at final time: %.2E\n", error_norm)
            end

            return (final_time = current_time, error = error_norm, solution = solution)
        end
    finally
        if did_initialize && finalize_petsc
            PETSc.finalize(petsclib)
        end
    end
end

if !isinteractive() && abspath(PROGRAM_FILE) == @__FILE__
    solve_ex51(finalize_petsc = true)
end
