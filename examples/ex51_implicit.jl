using PETSc, MPI, Printf

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
# PETSc 3.22 can solve this with `TSIRK` and `-ts_irk_type gauss`, but the
# internal IRK setup currently requires an AIJ-style sparse Jacobian matrix so
# it can form a `MATKAIJ` stage operator. In practice this means:
#
# - `seqdense` Jacobians do not work with Gauss/IRK here
# - `-snes_mf` / `-snes_mf_operator` do not work here
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
#     julia --project=. examples/ex51_implicit.jl -ts_irk_nstages 2
#
# or programmatically
#
#     solve_ex51_implicit(options = ["-ts_irk_nstages", "3"])

mutable struct Ex51Context{PetscLib <: PETSc.LibPETSc.PetscLibType}
    petsclib::PetscLib
end

"""
    ex51_rhs_ifunction!

Implicit residual routine passed to `TSSetIFunction`, corresponding to the DAE
form `F(t, u, u_t) = 0` with `F(t, u, u_t) = u_t - f(t, u)`.
"""
function ex51_rhs_ifunction!(
    ::PETSc.LibPETSc.CTS,
    t::Float64,
    u_ptr::PETSc.LibPETSc.CVec,
    udot_ptr::PETSc.LibPETSc.CVec,
    f_ptr::PETSc.LibPETSc.CVec,
    ctx_ptr::Ptr{Cvoid},
)::PETSc.LibPETSc.PetscErrorCode
    # This function is called through PETSc's C callback interface, so its
    # argument types must match the low-level PETSc ABI. The current state,
    # time derivative, and residual vectors arrive as raw `CVec` pointers. We
    # wrap them with `VecPtr(..., own = false)` so we can use PETSc.jl's array
    # helpers without taking ownership away from PETSc.
    ctx = unsafe_pointer_to_objref(ctx_ptr)::Ex51Context
    petsclib = ctx.petsclib
    # `own = false` since memory is managed by PETSc internally
    u = PETSc.VecPtr(petsclib, u_ptr, false)
    udot = PETSc.VecPtr(petsclib, udot_ptr, false)
    f = PETSc.VecPtr(petsclib, f_ptr, false)

    PETSc.withlocalarray!(
        (u, udot, f);
        read = (true, true, false),
        write = (false, false, true),
    ) do u_array, udot_array, f_array
        f_array[1] = udot_array[1] - cos(t)
        f_array[2] = udot_array[2] - sin(u_array[2])
    end

    return PETSc.LibPETSc.PetscErrorCode(0)
end

const EX51_RHS_IFUNCTION_PTR = @cfunction(
    ex51_rhs_ifunction!,
    PETSc.LibPETSc.PetscErrorCode,
    (
        PETSc.LibPETSc.CTS,
        Float64,
        PETSc.LibPETSc.CVec,
        PETSc.LibPETSc.CVec,
        PETSc.LibPETSc.CVec,
        Ptr{Cvoid},
    ),
)

"""
    ex51_fill_ijacobian!

Fill the tiny 2x2 implicit Jacobian matrix used by this example.

For the residual

    F(t, u, u_t) = [u1_t - cos(t), u2_t - sin(u2)],

the implicit Jacobian is

    dF/du + shift * dF/du_t = [shift   0;
                               0   shift - cos(u2)].
"""
function ex51_fill_ijacobian!(
    A::PETSc.LibPETSc.PetscMat,
    B::PETSc.LibPETSc.PetscMat,
    shift,
    diag22,
    petsclib,
)
    PetscScalar = petsclib.PetscScalar

    PETSc.LibPETSc.MatZeroEntries(petsclib, A)
    A[1, [1, 2]] = PetscScalar.([shift, 0])
    A[2, [1, 2]] = PetscScalar.([0, diag22])
    PETSc.assemble!(A)

    if B.ptr != A.ptr
        PETSc.LibPETSc.MatZeroEntries(petsclib, B)
        B[1, [1, 2]] = PetscScalar.([shift, 0])
        B[2, [1, 2]] = PetscScalar.([0, diag22])
        PETSc.assemble!(B)
    end

    return nothing
end

"""
    ex51_ijacobian!

Optional analytic implicit Jacobian callback for the example.

The default solve path uses PETSc's finite-difference coloring instead, but
keeping this explicit callback available is useful as a reference and a
regression-tested fallback.
"""
function ex51_ijacobian!(
    ::PETSc.LibPETSc.CTS,
    ::Float64,
    u_ptr::PETSc.LibPETSc.CVec,
    ::PETSc.LibPETSc.CVec,
    shift::Float64,
    A_ptr::PETSc.LibPETSc.CMat,
    B_ptr::PETSc.LibPETSc.CMat,
    ctx_ptr::Ptr{Cvoid},
)::PETSc.LibPETSc.PetscErrorCode
    ctx = unsafe_pointer_to_objref(ctx_ptr)::Ex51Context
    petsclib = ctx.petsclib
    # `u` is borrowed from PETSc; do not take ownership.
    u = PETSc.VecPtr(petsclib, u_ptr, false)
    A = PETSc.LibPETSc.PetscMat(A_ptr, petsclib)
    B = PETSc.LibPETSc.PetscMat(B_ptr, petsclib)

    diag22 = PETSc.withlocalarray!(u; read = true, write = false) do u_array
        shift - cos(u_array[2])
    end

    ex51_fill_ijacobian!(A, B, shift, diag22, petsclib)
    return PETSc.LibPETSc.PetscErrorCode(0)
end

const EX51_IJACOBIAN_PTR = @cfunction(
    ex51_ijacobian!,
    PETSc.LibPETSc.PetscErrorCode,
    (
        PETSc.LibPETSc.CTS,
        Float64,
        PETSc.LibPETSc.CVec,
        PETSc.LibPETSc.CVec,
        Float64,
        PETSc.LibPETSc.CMat,
        PETSc.LibPETSc.CMat,
        Ptr{Cvoid},
    ),
)

function ex51_initial_condition!(u::PETSc.LibPETSc.PetscVec)
    PETSc.withlocalarray!(u; read = false, write = true) do u_array
        u_array[1] = 0.0
        u_array[2] = 1.0
    end
    return nothing
end

function ex51_exact_solution!(u::PETSc.LibPETSc.PetscVec, t::Real)
    PETSc.withlocalarray!(u; read = false, write = true) do u_array
        u_array[1] = sin(t)
        u_array[2] = 2 * atan(exp(t) * tan(0.5))
    end
    return nothing
end

function ex51_implicit_default_options(parsed_options::NamedTuple, jacobian_mode::Symbol)
    # PETSc 3.22 currently rejects matrix-free operators inside the `TSIRK`
    # Gauss setup path for this problem, so fail early with a clear message
    # instead of letting PETSc error out later inside `TSSetUp_IRK`.
    if haskey(parsed_options, :snes_mf) || haskey(parsed_options, :snes_mf_operator)
        throw(
            ArgumentError(
                "PETSc 3.22 TSIRK/Gauss does not support matrix-free Jacobians here. " *
                "Use the AIJ-backed default path or `jacobian_mode = :analytic` instead.",
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
    jac = PETSc.MatSeqAIJ(petsclib, 2, 2, petsclib.PetscInt(2))
    jac[1, [1, 2]] = PetscScalar.([1, 1])
    jac[2, [1, 2]] = PetscScalar.([1, 1])
    PETSc.assemble!(jac)
    PETSc.LibPETSc.MatZeroEntries(petsclib, jac)
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

- The default path is not matrix-free. PETSc 3.22 `TSIRK/Gauss` requires an
  AIJ sparse Jacobian matrix and currently rejects `seqdense` and `MATMFFD`
  operators in this setup.
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
    PetscScalar = petsclib.PetscScalar

    # `did_initialize`: whether we initialized the library in this call
    did_initialize = !PETSc.initialized(petsclib)
    if did_initialize
        PETSc.initialize(petsclib)
    end

    MPI.Comm_size(comm) == 1 || error("This example only supports sequential runs.")

    # Keep these variables concretely typed across the whole function, even
    # though the actual PETSc objects are created later inside the `try` block.
    # The null-pointer placeholders are overwritten by `TSCreate`/`VecSeq`/
    # `MatSeqAIJ`, and the `finally` block checks `ptr != C_NULL` before
    # destroying them.
    ts = PETSc.LibPETSc.TS(petsclib)
    u = PETSc.LibPETSc.PetscVec(petsclib)
    u_exact = PETSc.LibPETSc.PetscVec(petsclib)
    jac = PETSc.LibPETSc.PetscMat(petsclib)
    current_time = petsclib.PetscReal(NaN)
    error_norm = petsclib.PetscReal(NaN)
    solution = PetscScalar[]
    ctx = Ex51Context(petsclib)
    petsc_options = PETSc.Options(petsclib; parsed_options...)
    pushed_options = false

    try
        # Create timestepping solver context. We start from PETSc's implicit
        # Runge-Kutta (`irk`) implementation and let options choose the Gauss
        # family member and the number of stages. By default the options path
        # below keeps `-ts_irk_type gauss`.
        ts = PETSc.LibPETSc.TSCreate(petsclib, comm)
        PETSc.LibPETSc.TSSetType(petsclib, ts, "irk")
        PETSc.LibPETSc.TSSetProblemType(petsclib, ts, PETSc.LibPETSc.TS_NONLINEAR)

        # Set initial conditions.
        u = PETSc.VecSeq(petsclib, 2)
        ex51_initial_condition!(u)
        PETSc.LibPETSc.TSSetSolution(petsclib, ts, u)

        # Build the sparse AIJ Jacobian template required by PETSc's
        # `TSIRK/Gauss` implementation.
        jac = ex51_implicit_create_jacobian_template(petsclib)

        callback_ctx_ptr = pointer_from_objref(ctx)
        GC.@preserve ctx begin
            # Register the implicit residual and Jacobian callbacks. When
            # `jacobian_mode = :finite_difference_color`, PETSc uses the AIJ
            # sparsity pattern from `jac` and computes Jacobian values
            # automatically via finite differences and coloring. The analytic
            # callback remains available as a fallback and documentation aid.
            PETSc.LibPETSc.TSSetIFunction(
                petsclib,
                ts,
                nothing,
                EX51_RHS_IFUNCTION_PTR,
                callback_ctx_ptr,
            )
            PETSc.LibPETSc.TSSetIJacobian(
                petsclib,
                ts,
                jac,
                jac,
                EX51_IJACOBIAN_PTR,
                callback_ctx_ptr,
            )

            # Configure solver options. As in `ex51.jl`, adaptivity is forced
            # to take constant time steps. Users can change the Gauss stage
            # count through `-ts_irk_nstages <s>` in `options` or on the
            # command line.
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
        # final time, matching the `VecAYPX` + `VecNorm` pattern used in
        # `ex51.jl`.
        current_time = PETSc.LibPETSc.TSGetTime(petsclib, ts)
        u_exact = similar(u)
        ex51_exact_solution!(u_exact, current_time)

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

            ts_type = PETSc.LibPETSc.TSGetType(petsclib, ts)
            @printf("TS type: %s\n", ts_type)
            if ts_type == "irk"
                @printf("IRK type: %s\n", PETSc.LibPETSc.TSIRKGetType(petsclib, ts))
                @printf("IRK stages: %d\n", PETSc.LibPETSc.TSIRKGetNumStages(petsclib, ts))
            end
            @printf("Jacobian mode: %s\n", String(Symbol(jacobian_mode)))
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
        if jac.ptr != C_NULL
            PETSc.destroy(jac)
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
    solve_ex51_implicit(finalize_petsc = true)
end
