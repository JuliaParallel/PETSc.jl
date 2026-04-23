# EXCLUDE FROM TESTING
#!/usr/bin/env julia
# parse_scaling.jl
# Usage: julia parse_scaling.jl scaling_*.out

using Printf

struct ScalingResult
    jobid            :: String
    ntasks           :: Int
    nx               :: Int
    ny               :: Int
    nz               :: Int
    mg_levels        :: Int
    ncoarse          :: Int
    solve_time       :: Float64   # application "Solve time:" field
    ksp_solve_time   :: Float64   # PETSc KSPSolve event time (for weak scaling)
    ksp_iter         :: Int
    l2_error         :: Float64
    max_error        :: Float64
    residual         :: Float64
    total_gflops     :: Float64
    gflops_s         :: Float64
    converged_reason :: String
    backend          :: String    # "native" or "PETSc.jl"
end

function parse_file(filename::String)
    m = match(r"scaling([a-zA-Z]*)_(\d+)_n(\d+)_nx(\d+)_ny(\d+)_nz(\d+)_mg(\d+)_nc(\d+)\.out", basename(filename))
    if isnothing(m)
        @warn "Filename $filename does not match expected pattern, skipping"
        return nothing
    end
    variant   = m[1]   # e.g. "", "Galerkin", "NATIVE"
    jobid     = m[2]
    ntasks    = parse(Int, m[3])
    nx        = parse(Int, m[4])
    ny        = parse(Int, m[5])
    nz        = parse(Int, m[6])
    mg_levels = parse(Int, m[7])
    ncoarse   = parse(Int, m[8])
    backend = if uppercase(variant) == "NATIVE"
        "native"
    elseif uppercase(variant) == "LOCALLIB"
        "local lib"
    else
        "PETSc.jl"
    end

    solve_time       = NaN
    ksp_solve_time   = NaN
    ksp_iter         = -1
    l2_error         = NaN
    max_error        = NaN
    residual         = NaN
    total_gflops     = NaN
    gflops_s         = NaN
    converged_reason = "UNKNOWN"

    in_petsc_summary = false
    in_ksp_monitor   = false

    for line in eachline(filename)

        # Detect PETSc performance summary block
        if occursin("PETSc Performance Summary", line)
            in_petsc_summary = true
        end

        # Detect KSP monitor lines (leading spaces + digit + " KSP")
        if match(r"^\s+\d+ KSP", line) !== nothing
            in_ksp_monitor = true
        end

        # KSP monitor block ends at the "Linear solve" line
        if occursin("Linear solve", line) && occursin("due to", line)
            in_ksp_monitor = false
            m2 = match(r"Linear solve (\w+) due to (\S+)", line)
            isnothing(m2) || (converged_reason = string(m2[1] == "converged" ? "OK " : "FAIL ") * m2[2])
            continue
        end

        # Skip KSP monitor lines entirely
        in_ksp_monitor && continue

        # ---- application output (before PETSc summary) ----
        if !in_petsc_summary

            if occursin("Outer KSP iterations:", line)
                m2 = match(r"Outer KSP iterations:\s+(\d+)", line)
                isnothing(m2) || (ksp_iter = parse(Int, m2[1]))

            elseif occursin("Final residual norm", line)
                m2 = match(r"Final residual norm\s+([\d.e+\-]+)", line)
                isnothing(m2) || (residual = parse(Float64, m2[1]))

            elseif occursin("Solve time:", line)
                m2 = match(r"Solve time:\s+([\d.e+\-]+)", line)
                isnothing(m2) || isnan(solve_time) && (solve_time = parse(Float64, m2[1]))

            elseif occursin("L2 error:", line)
                m2 = match(r"L2 error:\s+([\d.e+\-]+)", line)
                isnothing(m2) || isnan(l2_error) && (l2_error = parse(Float64, m2[1]))

            elseif occursin("Max error:", line)
                m2 = match(r"Max error:\s+([\d.e+\-]+)", line)
                isnothing(m2) || isnan(max_error) && (max_error = parse(Float64, m2[1]))
            end

        # ---- PETSc performance summary block ----
        else
            # KSPSolve event line: extract the Max time (column 3 after the counts)
            # Format: "KSPSolve   count ratio  maxtime ratio  ..."
            if match(r"^\s*KSPSolve\s", line) !== nothing && isnan(ksp_solve_time)
                # fields: Name Count Ratio MaxTime Ratio ...
                vals = split(strip(line))
                # vals[1]=KSPSolve, vals[2]=count, vals[3]=ratio, vals[4]=maxtime
                if length(vals) >= 4
                    t = tryparse(Float64, vals[4])
                    isnothing(t) || (ksp_solve_time = t)
                end
            end

            if match(r"^\s*Flops:\s", line) !== nothing
                vals = split(strip(replace(line, r"^\s*Flops:\s*" => "")))
                length(vals) >= 4 && (total_gflops = parse(Float64, vals[4]) / 1e9)

            elseif match(r"^\s*Flops/sec:\s", line) !== nothing
                vals = split(strip(replace(line, r"^\s*Flops/sec:\s*" => "")))
                length(vals) >= 4 && (gflops_s = parse(Float64, vals[4]) / 1e9)
            end
        end
    end

    return ScalingResult(jobid, ntasks, nx, ny, nz, mg_levels, ncoarse,
                         solve_time, ksp_solve_time, ksp_iter,
                         l2_error, max_error, residual,
                         total_gflops, gflops_s,
                         converged_reason, backend)
end

function ndofs(r::ScalingResult)
    return r.nx * r.ny * r.nz
end

function print_table(results::Vector{ScalingResult})
    sort!(results, by = r -> (r.ntasks, r.backend))

    println()
    @printf("%-12s  %6s  %4s  %4s  %15s  %9s  %10s  %11s  %10s  %10s  %4s  %12s  %12s  %12s  %s\n",
            "JobID", "Ntasks", "MG", "NC", "Grid", "Backend",
            "SolveTime", "KSPSolveTime", "TotGFlops",
            "GFlops/s", "KSP", "L2_error", "Max_error", "Residual", "Converged")
    println(repeat("-", 190))

    for r in results
        grid  = @sprintf("%dx%dx%d", r.nx, r.ny, r.nz)
        ksp_t = isnan(r.ksp_solve_time) ? "        N/A" : @sprintf("%11.3f", r.ksp_solve_time)
        @printf("%-12s  %6d  %4d  %4d  %15s  %9s  %10.3f  %s  %10.2f  %10.2f  %4d  %12.4e  %12.4e  %12.4e  %s\n",
                r.jobid, r.ntasks, r.mg_levels, r.ncoarse, grid, r.backend,
                r.solve_time, ksp_t, r.total_gflops, r.gflops_s,
                r.ksp_iter, r.l2_error, r.max_error, r.residual,
                r.converged_reason)
    end
    println()
end

function print_convergence(results::Vector{ScalingResult})
    iso = filter(r -> r.nx == r.ny == r.nz, results)
    sort!(iso, by = r -> r.nx)

    if length(iso) >= 2
        println("Spatial convergence rate - isotropic runs (should be ~2.0 for second-order FD):")
        println(repeat("-", 65))
        for i in 2:length(iso)
            r_fine   = iso[i]
            r_coarse = iso[i-1]
            h_fine   = 1.0 / (r_fine.nx   - 1)
            h_coarse = 1.0 / (r_coarse.nx - 1)
            if !isnan(r_fine.l2_error) && !isnan(r_coarse.l2_error)
                rate = log(r_coarse.l2_error / r_fine.l2_error) / log(h_coarse / h_fine)
                @printf("  %dx%dx%d -> %dx%dx%d  L2: %.4e->%.4e  rate= %.2f  [%s]\n",
                        r_coarse.nx, r_coarse.ny, r_coarse.nz,
                        r_fine.nx,   r_fine.ny,   r_fine.nz,
                        r_coarse.l2_error, r_fine.l2_error, rate, r_fine.backend)
            end
        end
        println()
    end
end

function print_scaling(results::Vector{ScalingResult})
    # Print separate weak scaling tables per backend
    backends = unique(r.backend for r in results)
    for backend in sort(backends)
        subset = filter(r -> r.backend == backend, results)
        sorted = sort(subset, by = r -> r.ntasks)
        isempty(sorted) && continue

        ref    = sorted[1]
        t0_ksp   = ref.ksp_solve_time
        t0_solve = ref.solve_time
        use_ksp  = !isnan(t0_ksp)
        t0       = use_ksp ? t0_ksp : t0_solve

        label = use_ksp ? "KSPSolve time" : "Solve time (fallback)"
        println("Weak scaling efficiency [$backend] based on $label (relative to smallest ntasks):")
        println(repeat("-", 85))
        @printf("  %-6s  %-15s  %8s  %11s  %10s  %10s  %s\n",
                "Ntasks", "Grid", "DOFs/core", "KSPSolve(s)", "SolveTime", "Efficiency", "Converged")
        println(repeat("-", 85))
        for r in sorted
            t_ksp = isnan(r.ksp_solve_time) ? NaN : r.ksp_solve_time
            t_use = use_ksp ? t_ksp : r.solve_time
            eff   = isnan(t_use) ? NaN : t0 / t_use * 100
            dpc   = ndofs(r) / r.ntasks
            grid_str = @sprintf("%dx%dx%d", r.nx, r.ny, r.nz)
            ksp_str  = isnan(t_ksp) ? "        N/A" : @sprintf("%11.3f", t_ksp)
            eff_str  = isnan(eff)   ? "       N/A" : @sprintf("%9.1f%%", eff)
            @printf("  %-6d  %-15s  %8.0f  %s  %10.3f  %s  %s\n",
                    r.ntasks, grid_str, dpc, ksp_str, r.solve_time, eff_str, r.converged_reason)
        end
        println()
    end
end

function main()
    if isempty(ARGS)
        println("Usage: julia parse_scaling.jl scaling_*.out")
        exit(1)
    end

    results = ScalingResult[]
    for f in ARGS
        isfile(f) || (@warn "File not found: $f, skipping"; continue)
        r = parse_file(f)
        isnothing(r) || push!(results, r)
    end

    isempty(results) && (println("No valid files found."); exit(1))

    print_table(results)
    print_convergence(results)
    print_scaling(results)
end

main()