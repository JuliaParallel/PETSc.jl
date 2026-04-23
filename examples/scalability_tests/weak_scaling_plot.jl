using GLMakie
using Statistics

# All KSPSolve times per (ntasks, backend)
# Multiple runs where available

data = Dict(
    :native => Dict(
        64   => [26.917, 27.155, 26.854],
        512  => [31.492, 30.122, 30.607, 30.307],
        4096 => [30.694, 28.214, 28.498],
        32768 => [50.858],
    ),
    :petscjl => Dict(
        64   => [27.916, 27.108, 28.121],
        512  => [33.670, 32.482, 34.236],
        4096 => [34.332, 44.244],
        32768 => [63.300],
    ),
    :locallib => Dict(
        64   => [27.065, 27.006, 27.551],
        512  => [34.456, 32.033, 30.708],
        4096 => [38.414, 35.292, 34.220],
        32768 =>[73.108],   
    ),
)

ntasks = [64, 512, 4096, 32768]

function stats(d, backend)
    means = [isempty(d[backend][n]) ? NaN : mean(d[backend][n]) for n in ntasks]
    mins  = [isempty(d[backend][n]) ? NaN : minimum(d[backend][n]) for n in ntasks]
    maxs  = [isempty(d[backend][n]) ? NaN : maximum(d[backend][n]) for n in ntasks]
    return means, mins, maxs
end

mean_native,   min_native,   max_native   = stats(data, :native)
mean_petscjl,  min_petscjl,  max_petscjl  = stats(data, :petscjl)
mean_locallib, min_locallib, max_locallib = stats(data, :locallib)

fig = Figure(size = (900, 600), fontsize = 14)

ax = Axis(fig[1, 1],
    title  = "Weak Scaling on LUMI-C | KSPSolve Time (3D Poisson, 64 tasks/node), ex45.jl",
    xlabel = "Number of MPI tasks, [number of nodes], (global grid size)",
    ylabel = "Solution time (KSPSolve) (s)",
    xscale = log2,
    xticks = (ntasks, ["64\n[1]\n(513³)", "512\n[8]\n(1025³)", "4096\n[64]\n(2049³)", "32768\n[512]\n(4097³)"]),
    xminorticksvisible = false,
    yminorticksvisible = false,
    limits = (nothing, nothing, 0, 75),
)

col_native   = :black
col_petscjl  = :dodgerblue
col_locallib = :tomato

# Shaded min-max bands (only where data exists)
# native: all 4 points
band!(ax, ntasks, min_native,   max_native,   color = (col_native,   0.10))
# petscjl: all 4 points
band!(ax, ntasks, min_petscjl,  max_petscjl,  color = (col_petscjl,  0.15))
# locallib: all 4 points
band!(ax, ntasks[1:4], min_locallib[1:4], max_locallib[1:4], color = (col_locallib, 0.15))

# Mean lines
lines!(ax, ntasks,       mean_native,              color = col_native,   linewidth = 2.5, label = "Native C")
lines!(ax, ntasks,       mean_petscjl,             color = col_petscjl,  linewidth = 2.5, label = "PETSc.jl (_jll)")
lines!(ax, ntasks,  mean_locallib[1:4],       color = col_locallib, linewidth = 2.5, label = "Local lib")
# Dashed extension hint for local lib pending
scatter!(ax, [32768], [NaN], color = col_locallib, markersize = 14, marker = :diamond,
         strokewidth = 1.5, strokecolor = :white)  # placeholder, won't render NaN

# Individual run scatter points
for (backend, col, mkr) in [(:native, col_native, :circle), (:petscjl, col_petscjl, :rect), (:locallib, col_locallib, :diamond)]
    for n in ntasks
        vals = data[backend][n]
        isempty(vals) && continue
        scatter!(ax, fill(n, length(vals)), vals, color = col, markersize = 9, marker = mkr)
    end
end

# Mean markers (larger, white outline)
scatter!(ax, ntasks,      mean_native,        color = col_native,   markersize = 14, marker = :circle,  strokewidth = 1.5, strokecolor = :white)
scatter!(ax, ntasks,      mean_petscjl,       color = col_petscjl,  markersize = 14, marker = :rect,    strokewidth = 1.5, strokecolor = :white)
scatter!(ax, ntasks,  mean_locallib,       color = col_locallib, markersize = 14, marker = :diamond, strokewidth = 1.5, strokecolor = :white)

# Ideal flat line based on native baseline
hlines!(ax, [mean_native[1]], color = :gray, linewidth = 1.5, linestyle = :dash, label = "Ideal (flat)")

axislegend(ax, position = :lt, framevisible = true, labelsize = 13)

Label(fig[2, 1],
    "Shaded bands show min–max range across repeated runs. Large markers show mean.",
    fontsize = 11, color = :gray50, tellwidth = false)

save("weak_scaling.png", fig, px_per_unit = 2)
println("Saved weak_scaling.png")
