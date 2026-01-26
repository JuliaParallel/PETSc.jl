# INCLUDE IN MPI TEST
#
#!/usr/bin/env julia
using PETSc, MPI, Printf

# Convergence test usage (examples)
#
# Run convergence tests for `ex45` (3D) or `ex50` (2D) from the repository root.
# - Direct (LU) baseline for ex45 (3D):
#     julia --project=. examples/convergence_test.jl -Ns="65,129,257" -levels="1,1,1" -ksp_type preonly -pc_type lu -example ex45
# - Geometric multigrid for ex45:
#     julia --project=. examples/convergence_test.jl -Ns="65,129,257,513" -levels="3,4,5,6" -ksp_type cg -pc_type mg -example ex45
#
# - Direct (LU) baseline for ex50 (2D):
#     julia --project=. examples/convergence_test.jl -Ns="65,129,257" -levels="1,1,1" -ksp_type preonly -pc_type lu -example ex50
# - Geometric multigrid for ex50:
#     julia --project=. examples/convergence_test.jl -Ns="65,129,257,513" -levels="3,4,5,6" -ksp_type cg -pc_type mg -example ex50
#
# Notes:
# - `-levels` maps one MG level per entry in `-Ns` (comma-separated).
# - Use `-ksp_type preonly -pc_type lu` for a direct factor baseline.
# - Use `-ksp_type cg -pc_type mg` (and tune MG-specific flags) for geometric multigrid.
# - PETSc CLI flags (e.g. `-ksp_monitor`, `-pc_mg_type`, `-mg_levels_ksp_type`) are forwarded
#   to the solver; pass them on the command line to tune behavior at runtime.

# This script runs the selected example (`ex50` by default) for a sequence
# of grid resolutions and prints a convergence table similar to
# `ex50_run_convergence.jl`. Resolutions are shown as `NxNxN` (or `NxN` for 2D).

include(joinpath(@__DIR__, "ex45.jl"))
include(joinpath(@__DIR__, "ex50.jl"))

# Parse PETSc/CLI options and build a kwargs dict to forward to `solve_ex45`.
opts = PETSc.parse_options(ARGS)
ns_opt = get(opts, Symbol("Ns"), nothing)
N_start = parse(Int, get(opts, Symbol("N_start"), "9"))
reverse_order = get(opts, Symbol("reverse"), false)
example = get(opts, Symbol("example"), "ex50")
dim = parse(Int, get(opts, Symbol("dim"), example == "ex50" ? "2" : "3"))
levels_opt = get(opts, Symbol("levels"), nothing)

if ns_opt !== nothing
    Ns = [parse(Int, strip(s)) for s in split(string(ns_opt), ",") if strip(s) != ""]
else
    Ns = [N_start * (2^(i-1)) for i in 1:5]
end
if reverse_order
    Ns = reverse(Ns)
end

# Parse explicit per-resolution MG levels if provided (comma-separated)
if levels_opt !== nothing
    Levels = [parse(Int, strip(s)) for s in split(string(levels_opt), ",") if strip(s) != ""]
    if length(Levels) != length(Ns)
        error("Provided -levels must contain the same number of entries as Ns")
    end
else
    Levels = nothing
end

# Forwardable solver options: everything in opts except the control args above
forward_keys = Set([:Ns, :N_start, :reverse, :dim, :example])
solver_kwargs = Dict{Symbol,Any}()
for k in keys(opts)
    if !(k in forward_keys)
        solver_kwargs[k] = opts[k]
    end
end

petsclib = PETSc.getlib(; PetscScalar = Float64)
PETSc.initialize(petsclib)

results = []
for (i, N) in enumerate(Ns)
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        if example == "ex50"
            @printf("Running ex50 with N=%d\n", N)
        else
            @printf("Running ex45 with N=%d\n", N)
        end
    end
    # Build per-run options dict (copy solver_kwargs)
    opts_dict = Dict(solver_kwargs)
    # If user provided explicit per-resolution MG levels, use them
    if Levels !== nothing
        opts_dict[:pc_mg_levels] = string(Levels[i])
    else
        # keep any provided pc_mg_levels or leave absent
        if get(opts_dict, :pc_type, "") == "mg" && !haskey(opts_dict, :pc_mg_levels)
            # default geometric mapping: increase levels with resolution
            opts_dict[:pc_mg_levels] = string(1 + (i - 1))
        end
    end

    # For the smallest resolution run twice to exclude compilation time
    if i == 1
        if example == "ex50"
            solve_ex50(N; opts_dict...)  # warm-up / discard
        else
            solve_ex45(N; opts_dict...)  # warm-up / discard
        end
    end
    # call selected solver programmatically, forwarding any solver kwargs
    if example == "ex50"
        res = solve_ex50(N; opts_dict...)
    else
        res = solve_ex45(N; opts_dict...)
    end
    # res is a NamedTuple: (norm, final_grid, niter, solve_time, L2, max)
    final_grid = res.final_grid
    h = 1.0 / N
    push!(results, (N, h, res.L2, res.solve_time, res.niter))
end

# Compute orders and time-scaling (𝒪(Time))
orders = fill(NaN, length(results))
time_orders = fill(NaN, length(results))
for i in 1:length(results)-1
    e1 = results[i][3]
    e2 = results[i+1][3]
    h1 = results[i][2]
    h2 = results[i+1][2]
    if e1 > 0 && e2 > 0
        orders[i] = log(e1 / e2) / log(h1 / h2)
    end
    t1 = results[i][4]
    t2 = results[i+1][4]
    N1 = results[i][1]
    N2 = results[i+1][1]
    if t1 > 0 && t2 > 0 && N1 > 0 && N2 > 0
        time_orders[i] = log(t2 / t1) / (dim * log(N2 / N1))
    end
end

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    # Header matching ex50_convergence.jl but hide MG Lvl unless PC is multigrid
    pc_type_val = get(solver_kwargs, :pc_type, "gamg")
    show_mg = pc_type_val == "mg"
    @printf("\nConvergence Analysis Results (KSP: %s, PC: %s):\n", get(solver_kwargs, :ksp_type, "cg"), pc_type_val)
    if show_mg
        @printf("%-11s %-8s %-10s %-10s %-12s %-8s %-10s %-8s\n", "N", "MG Lvl", "h", "KSP Iters", "L2 Error", "𝒪(N)", "Time (s)", "𝒪(Time)")
        @printf("%-11s %-8s %-10s %-10s %-12s %-8s %-10s %-8s\n", repeat("-",11), repeat("-",8), repeat("-",10), repeat("-",10), repeat("-",12), repeat("-",8), repeat("-",10), repeat("-",8))
    else
        @printf("%-11s %-10s %-10s %-12s %-8s %-10s %-8s\n", "N", "h", "KSP Iters", "L2 Error", "𝒪(N)", "Time (s)", "𝒪(Time)")
        @printf("%-11s %-10s %-10s %-12s %-8s %-10s %-8s\n", repeat("-",11), repeat("-",10), repeat("-",10), repeat("-",12), repeat("-",8), repeat("-",10), repeat("-",8))
    end

    for i in 1:length(results)
        Nval, h, e, time, iters = results[i]
        ord = orders[i]
        tord = time_orders[i]
        # Format N as NxN or NxNxN depending on spatial dimension
        Nstr = dim == 3 ? string(Nval, "x", Nval, "x", Nval) : string(Nval, "x", Nval)
        if show_mg
            # display MG level for this row
            mglvl = get(solver_kwargs, :pc_mg_levels, "--")
            if Levels !== nothing
                mglvl = string(Levels[i])
            else
                if haskey(solver_kwargs, :pc_mg_levels)
                    mglvl = string(solver_kwargs[:pc_mg_levels])
                else
                    if pc_type_val == "mg"
                        mglvl = string(1 + (i - 1))
                    end
                end
            end

            if i == 1
                @printf("%-11s %-8s %-10.6f %-10d %-12.6e %-8s %-10.4f %-8s\n", Nstr, mglvl, h, iters, e, "-", time, "-")
            else
                s_ord = isnan(ord) ? "N/A" : @sprintf("%8.2f", ord)
                s_tord = isnan(tord) ? "N/A" : @sprintf("%8.2f", tord)
                if s_ord == "N/A" && s_tord == "N/A"
                    @printf("%-11s %-8s %-10.6f %-10d %-12.6e %-8s %-10.4f %-8s\n", Nstr, mglvl, h, iters, e, "N/A", time, "N/A")
                elseif s_ord == "N/A"
                    @printf("%-11s %-8s %-10.6f %-10d %-12.6e %-8s %-10.4f %-8.2f\n", Nstr, mglvl, h, iters, e, "N/A", time, tord)
                elseif s_tord == "N/A"
                    @printf("%-11s %-8s %-10.6f %-10d %-12.6e %-8.2f %-10.4f %-8s\n", Nstr, mglvl, h, iters, e, ord, time, "N/A")
                else
                    @printf("%-11s %-8s %-10.6f %-10d %-12.6e %-8.2f %-10.4f %-8.2f\n", Nstr, mglvl, h, iters, e, ord, time, tord)
                end
            end
        else
            if i == 1
                @printf("%-11s %-10.6f %-10d %-12.6e %-8s %-10.4f %-8s\n", Nstr, h, iters, e, "-", time, "-")
            else
                s_ord = isnan(ord) ? "N/A" : @sprintf("%8.2f", ord)
                s_tord = isnan(tord) ? "N/A" : @sprintf("%8.2f", tord)
                if s_ord == "N/A" && s_tord == "N/A"
                    @printf("%-11s %-10.6f %-10d %-12.6e %-8s %-10.4f %-8s\n", Nstr, h, iters, e, "N/A", time, "N/A")
                elseif s_ord == "N/A"
                    @printf("%-11s %-10.6f %-10d %-12.6e %-8s %-10.4f %-8.2f\n", Nstr, h, iters, e, "N/A", time, tord)
                elseif s_tord == "N/A"
                    @printf("%-11s %-10.6f %-10d %-12.6e %-8.2f %-10.4f %-8s\n", Nstr, h, iters, e, ord, time, "N/A")
                else
                    @printf("%-11s %-10.6f %-10d %-12.6e %-8.2f %-10.4f %-8.2f\n", Nstr, h, iters, e, ord, time, tord)
                end
            end
        end
    end
end

PETSc.finalize(petsclib)
