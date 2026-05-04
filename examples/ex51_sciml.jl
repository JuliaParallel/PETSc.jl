using PETSc, SciMLBase
using DiffEqBase # TODO: do we need this dependency?
using LinearAlgebra

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
#
# This version of the example demonstrates how to use the SciML integration
# of PETSc.jl, which allows to solve `SciMLBase.ODEProblem`s with PETSc's
# time-stepping solvers (TS).

function ex51_rhs!(du, u, p, t)
    du[1] = cos(t)
    du[2] = sin(u[2])
    return nothing
end

function exact_solution!(u, t::Real)
    u[1] = sin(t)
    u[2] = 2 * atan(exp(t) * tan(oftype(t, 0.5)))
    return nothing
end

function solve_ex51(;
    final_time::Real = 1.0,
    dt::Real = 0.25,
    alg = PETSc.TSRK("5dp"),
    verbose::Bool = true,
    kwargs... # further keyword arguements passed to the integrator/`solve`
)
    # Create the ODE problem
    RealT = typeof(final_time)
    tspan = (zero(RealT), final_time)
    u0 = zeros(RealT, 2)
    exact_solution!(u0, first(tspan))
    ode = ODEProblem(ex51_rhs!, u0, tspan)

    # Solve the ODE using an interface matching (mostly) the typical SciML interface
    sol = solve(ode, alg; dt, kwargs...)

    # Compute the error against the analytical solution at the achieved
    # final time.
    current_time = sol.t[end]
    u_exact = similar(u)
    exact_solution!(u_exact, current_time)

    error_norm = norm(sol.u[end] - u_exact)

    if verbose
        if !(current_time ≈ final_time)
            @warn "Note: prescribed `final_time` differs from `current_time`" final_time current_time
        end
        println("Error at final time: ", error_norm)
    end

    return (final_time = current_time, error = error_norm, solution = sol.u[end])
end
