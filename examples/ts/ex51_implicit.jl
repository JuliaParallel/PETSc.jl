using PETSc, MPI, Printf
using LinearAlgebra: norm

# Small ODE to test implicit TS accuracy with Gauss/IRK schemes; adapted from
# https://petsc.org/main/src/ts/tutorials/ex51.c.html.
#
# The ODE
#                 u1_t = cos(t),
#                 u2_t = sin(u2)
# with analytical solution
#                 u1(t) = sin(t),
#                 u2(t) = 2 * atan(exp(t) * tan(0.5))
# is rewritten in implicit form as
#                 F(t, u, u_t) = u_t - f(t, u) = 0.
#
# PETSc 3.25 can solve this with `TSIRK` and `-ts_irk_type gauss`. The IRK setup
# builds a `MATKAIJ` stage operator on top of the Jacobian matrix, so that
# matrix must be AIJ. In practice this means:
#
# - `seqdense` Jacobians do not work with Gauss/IRK here
# - `-snes_mf_operator` does not work here, since PETSc still expects the stage
#   operator to be `MATKAIJ`
# - `-snes_mf` works: the whole stage operator becomes matrix-free
# - a sparse AIJ matrix is required even if PETSc computes the Jacobian values
#   for us via finite differences and coloring
#
# This example therefore creates a tiny `seqaij` Jacobian template and, by
# default, asks PETSc to fill in the values with `-snes_fd_color` so callers do
# not have to provide Jacobian entries manually.
#
# The number of Gauss stages is controlled by PETSc option
# `-ts_irk_nstages <s>`. For example:
#
#     julia --project=. examples/ts/ex51_implicit.jl -ts_irk_nstages 2
#
# or programmatically
#
#     solve_ex51_implicit(options = ["-ts_irk_nstages", "3"])

"""
    ex51_rhs_ifunction!(F, ts, t, u, u_t)

Implicit residual, corresponding to the DAE form `F(t, u, u_t) = 0` with
`F(t, u, u_t) = u_t - f(t, u)`. It fills the vector `F`.
"""
function ex51_rhs_ifunction!(F, ts, t, u, udot)
    PETSc.with_local_array!(
        (u, udot, F);
        read = (true, true, false),
        write = (false, false, true),
    ) do u_array, udot_array, f_array
        f_array[1] = udot_array[1] - cos(t)
        f_array[2] = udot_array[2] - sin(u_array[2])
    end
    return nothing
end

"""
    ex51_ijacobian(PetscScalar)

The optional analytic implicit Jacobian callback `(A, P, ts, t, u, u_t, shift)`
for the example, for a library with scalar type `PetscScalar`.

The default solve path uses PETSc's finite-difference coloring instead, but
keeping this explicit callback available is useful as a reference and a
regression-tested fallback. For the residual

    F(t, u, u_t) = [u1_t - cos(t), u2_t - sin(u2)],

the implicit Jacobian is

    dF/du + shift * dF/du_t = [shift   0;
                               0   shift - cos(u2)].
"""
function ex51_ijacobian(PetscScalar)
    function fill_jacobian!(M, shift, diag22)
        fill!(M, 0)
        M[1, [1, 2]] = PetscScalar.([shift, 0])
        M[2, [1, 2]] = PetscScalar.([0, diag22])
        PETSc.assemble!(M)
    end

    return function (A, P, ts, t, u, udot, shift)
        diag22 = PETSc.with_local_array!(u; read = true, write = false) do u_array
            shift - cos(u_array[2])
        end
        fill_jacobian!(A, shift, diag22)
        A.ptr == P.ptr || fill_jacobian!(P, shift, diag22)
        return nothing
    end
end

function ex51_initial_condition!(u::PETSc.LibPETSc.PetscVec)
    PETSc.with_local_array!(u; read = false, write = true) do u_array
        u_array[1] = 0.0
        u_array[2] = 1.0
    end
    return nothing
end

function ex51_exact_solution!(u::PETSc.LibPETSc.PetscVec, t::Real)
    PETSc.with_local_array!(u; read = false, write = true) do u_array
        u_array[1] = sin(t)
        u_array[2] = 2 * atan(exp(t) * tan(0.5))
    end
    return nothing
end

function ex51_implicit_default_options(parsed_options::NamedTuple, jacobian_mode::Symbol)
    # PETSc 3.25 rejects `-snes_mf_operator` inside the `TSIRK` Gauss setup
    # path for this problem, so fail early with a clear message instead of
    # letting PETSc error out later inside `TSSetUp_IRK`.
    if haskey(parsed_options, :snes_mf_operator)
        throw(
            ArgumentError(
                "PETSc 3.25 TSIRK/Gauss does not support `-snes_mf_operator` here. " *
                "Use `-snes_mf` for a matrix-free solve, or the AIJ-backed default path.",
            ),
        )
    end

    jacobian_mode in (:finite_difference_color, :analytic) ||
        throw(ArgumentError("Unsupported jacobian_mode: $(jacobian_mode)"))

    effective = merge(
        (
            ts_type = "irk",
            ts_irk_type = "gauss",
            ksp_type = "gmres",
            pc_type = "none",
        ),
        parsed_options,
    )

    if jacobian_mode == :finite_difference_color &&
       !haskey(effective, :snes_mf) &&
       !haskey(effective, :snes_fd) &&
       !haskey(effective, :snes_fd_color)
        effective = merge(effective, (snes_fd_color = nothing,))
    end

    return effective
end

function ex51_implicit_create_jacobian_template(petsclib)
    PetscScalar = petsclib.PetscScalar

    # `TSIRK/Gauss` needs an AIJ Jacobian matrix. We use a tiny full 2x2 AIJ
    # pattern so PETSc can later compute the values with finite differences and
    # coloring, without the caller having to provide sparse Jacobian entries.
    jac = PETSc.PetscMat(petsclib, 2, 2, petsclib.PetscInt(2))
    jac[1, [1, 2]] = PetscScalar.([1, 1])
    jac[2, [1, 2]] = PetscScalar.([1, 1])
    PETSc.assemble!(jac)
    fill!(jac, 0)
    return jac
end

"""
    solve_ex51_implicit(; kwargs...)

Solve the small ODE from PETSc TS tutorial `ex51.c`, but in implicit form so
that PETSc's Gauss/IRK timesteppers can be used.

Keyword arguments:

- `petsclib`: PETSc library instance. Defaults to `Float64`.
- `final_time`: final integration time. Default `1.0`.
- `dt`: initial time step. Default `0.25`.
- `save_trajectory`: whether to call `TSSetSaveTrajectory`. Default `true`.
- `finalize_petsc`: whether to finalize PETSc at the end when this function had
  to initialize it. Default `false`.
- `options`: additional PETSc command-line options for this solve. By default
  the script uses `ARGS`. This is also where PETSc time-stepper options such as
  `-ts_irk_nstages <s>` are passed, for example
  `options = ["-ts_irk_nstages", "3"]`.
- `jacobian_mode`: `:finite_difference_color` (default) asks PETSc to compute
  Jacobian values automatically on top of an AIJ sparsity template. `:analytic`
  uses the built-in fallback callback in this file.
- `verbose`: print a short run summary.

Notes:

- The default path is not matrix-free. PETSc 3.25 `TSIRK/Gauss` requires an
  AIJ sparse Jacobian matrix and rejects `seqdense` and `-snes_mf_operator` in
  this setup. `-snes_mf` makes the whole solve matrix-free and works.
- The default `:finite_difference_color` mode avoids requiring callers to
  provide Jacobian values manually.
"""
function solve_ex51_implicit(;
    petsclib = PETSc.getlib(PetscScalar = Float64),
    final_time::Real = 1.0,
    dt::Real = 0.25,
    save_trajectory::Bool = true,
    finalize_petsc::Bool = false,
    options = nothing,
    jacobian_mode::Symbol = :finite_difference_color,
    verbose::Bool = true,
)
    options === nothing && (options = copy(String.(ARGS)))
    parsed_options = ex51_implicit_default_options(PETSc.parse_options(options), jacobian_mode)

    comm = MPI.COMM_WORLD

    # `did_initialize`: whether we initialized the library in this call
    did_initialize = !PETSc.isinitialized(petsclib)
    if did_initialize
        PETSc.initialize(petsclib)
    end

    MPI.Comm_size(comm) == 1 || error("This example only supports sequential runs.")

    try
        return PETSc.with_destroy() do list
            # The options are applied at the solve, after what is set below.
            ts = PETSc.track!(list, PETSc.TS(petsclib, comm;
                exact_final_time = PETSc.LibPETSc.TS_EXACTFINALTIME_STEPOVER,
                parsed_options...))
            # Start from PETSc's implicit Runge-Kutta (`irk`) and let the
            # options choose the Gauss family member and the number of stages.
            PETSc.set_type!(ts, :irk)
            PETSc.set_problem_type!(ts, PETSc.LibPETSc.TS_NONLINEAR)

            u = PETSc.track!(list, PETSc.PetscVec(petsclib, 2))
            ex51_initial_condition!(u)
            PETSc.set_solution!(ts, u)

            # The sparse AIJ Jacobian template `TSIRK/Gauss` needs. With
            # `jacobian_mode = :finite_difference_color` PETSc fills in the
            # values from its sparsity pattern by finite differences and
            # coloring. The analytic callback is the reference.
            jac = PETSc.track!(list, ex51_implicit_create_jacobian_template(petsclib))
            PETSc.set_ifunction!(ex51_rhs_ifunction!, ts)
            PETSc.set_ijacobian!(ex51_ijacobian(petsclib.PetscScalar), ts, jac)

            # As in `ex51.jl`, adaptivity is forced to take constant time steps.
            # Users can change the Gauss stage count through `-ts_irk_nstages <s>`
            # in `options` or on the command line.
            save_trajectory && PETSc.LibPETSc.TSSetSaveTrajectory(petsclib, ts)
            PETSc.set_max_time!(ts, final_time)
            PETSc.set_timestep!(ts, dt)
            PETSc.set_adapt_type!(ts, :none)

            PETSc.solve!(u, ts)

            # Error against the analytical solution at the achieved final time
            current_time = PETSc.current_time(ts)
            u_exact = PETSc.track!(list, similar(u))
            ex51_exact_solution!(u_exact, current_time)
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

                ts_type = PETSc.type_name(ts)
                @printf("TS type: %s\n", ts_type)
                if ts_type == :irk
                    @printf("IRK type: %s\n", PETSc.LibPETSc.TSIRKGetType(petsclib, ts))
                    @printf("IRK stages: %d\n", PETSc.LibPETSc.TSIRKGetNumStages(petsclib, ts))
                end
                @printf("Jacobian mode: %s\n", String(Symbol(jacobian_mode)))
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
    solve_ex51_implicit(finalize_petsc = true)
end
