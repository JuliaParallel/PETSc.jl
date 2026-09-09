using PETSc, SciMLBase

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
# In IMEX mode (`imex = true`) the right-hand side is split into a
# stiff/implicit part and a non-stiff/explicit part, matching the split used
# in examples/ex16.jl:
#
#     implicit (f1):  du[1] = 0,   du[2] = mu * ((1 - u1^2) * u2 - u1)
#     explicit (f2):  du[1] = u2,  du[2] = 0
#
# In non-IMEX mode (`imex = false`) the full right-hand side is treated as a
# single implicit term, matching the upstream `beuler` / `cn` / `theta` path.
#
# This version demonstrates how to use the SciML integration of PETSc.jl.
# Compare with examples/ex16.jl for the low-level PETSc TS interface.

# Stiff/implicit part of the van der Pol RHS (for SplitODEProblem).
# `p.mu` is the van der Pol parameter. The equation is autonomous, so `t`
# is unused.
function ex16_f1!(du, u, p, t)
    mu = p.mu
    du[1] = zero(u[1])
    du[2] = mu * ((1 - u[1]^2) * u[2] - u[1])
    return nothing
end

# Non-stiff/explicit part of the van der Pol RHS (for SplitODEProblem).
# This term has no parameter or time dependence.
function ex16_f2!(du, u, p, t)
    du[1] = u[2]
    du[2] = zero(u[2])
    return nothing
end

# Full right-hand side for the non-IMEX (plain ODEProblem) case.
function ex16_f!(du, u, p, t)
    mu = p.mu
    du[1] = u[2]
    du[2] = mu * ((1 - u[1]^2) * u[2] - u[1])
    return nothing
end

# Initial condition matching examples/ex16.jl (and upstream ex16.c).
function ex16_initial_condition(mu::Real)
    return Float64[
        2.0,
        -2.0 / 3.0 + 10.0 / (81.0 * mu) - 292.0 / (2187.0 * mu^2),
    ]
end

"""
    solve_ex16(; kwargs...)

Solve the van der Pol example from PETSc TS tutorial `ex16.c` using the SciML
interface of PETSc.jl.

In IMEX mode (`imex = true`, default) the problem is posed as a
`SplitODEProblem` and solved with `PETSc.TSARKIMEX`. In non-IMEX mode
(`imex = false`) a plain `ODEProblem` is solved with `PETSc.TSImplicit`.
Compare with `solve_ex16` in `examples/ex16.jl` for the equivalent low-level
PETSc TS interface.

# Keyword arguments
- `mu`: van der Pol stiffness parameter. Default `1000.0`.
- `imex`: split into stiff/non-stiff parts (IMEX) when `true`, fully implicit
  when `false`. Default `true`.
- `final_time`: final integration time. Default `0.5`.
- `dt`: initial time step. Default `0.01`.
- `alg`: PETSc TS algorithm. Defaults to `PETSc.TSARKIMEX("2e", ["-snes_fd"])`
  when `imex = true`, `PETSc.TSImplicit("beuler", ["-snes_fd"])` when
  `imex = false`. Any `PETSc.PETScTSAlgorithm` can be passed.
- `adaptive`: enable PETSc's adaptive time-step controller. Default `true`.
- `verbose`: print a summary after the solve. Default `true`.
- `kwargs...`: additional keyword arguments forwarded to `solve`.

# Returns
A named tuple `(final_time, steps, solution, mu, imex)` where `steps` is the
number of accepted time steps and `solution` is the state vector at the final
time.
"""
function solve_ex16(;
    mu::Real = 1000.0,
    imex::Bool = true,
    final_time::Real = 0.5,
    dt::Real = 0.01,
    alg = nothing,
    adaptive::Bool = true,
    verbose::Bool = true,
    kwargs...,
)
    u0 = ex16_initial_condition(mu)
    tspan = (0.0, Float64(final_time))
    p = (; mu = Float64(mu))

    if imex
        effective_alg = alg === nothing ? PETSc.TSARKIMEX("2e", ["-snes_fd"]) : alg
        prob = SplitODEProblem(ex16_f1!, ex16_f2!, u0, tspan, p)
    else
        effective_alg = alg === nothing ? PETSc.TSImplicit("beuler", ["-snes_fd"]) : alg
        prob = ODEProblem(ex16_f!, u0, tspan, p)
    end

    sol = solve(prob, effective_alg; dt, adaptive, kwargs...)

    current_time = sol.t[end]
    solution = sol.u[end]
    steps = sol.stats !== nothing ? sol.stats.naccept : length(sol.t) - 1

    if verbose
        if !(current_time ≈ final_time)
            @warn "prescribed `final_time` differs from `current_time`" final_time current_time
        end
        println("Algorithm:      ", effective_alg)
        println("mu:             ", mu)
        println("IMEX split:     ", imex)
        println("Final time:     ", current_time)
        println("Steps:          ", steps)
        println("Final solution: ", solution)
    end

    return (
        final_time = current_time,
        steps = steps,
        solution = solution,
        mu = mu,
        imex = imex,
    )
end

if !isinteractive() && abspath(PROGRAM_FILE) == @__FILE__
    solve_ex16()
end
