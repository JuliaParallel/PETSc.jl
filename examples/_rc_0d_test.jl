# 0D VEP test: mirrors rc_maxwell_stress_2d in ex62c.jl exactly.
# Checks whether RC's SeriesModel(LinearViscosity, IncompressibleElasticity, DruckerPrager)
# caps τ_II at the yield stress C = 0.5.
#
# Run: julia +1.12 --project=examples examples/rc_0d_test.jl

using RheologyCalculator, Printf, CairoMakie
include(joinpath(pkgdir(RheologyCalculator), "rheologies", "RheologyDefinitions.jl"))
include(joinpath(pkgdir(RheologyCalculator), "examples", "tensor_helpers.jl"))

# ── model — matches ex62c.jl defaults exactly ──────────────────────────────
η      = 1.0;  G = 5.0;  C = 0.5
η_vp   = 0.001   # ex62c.jl default eta_vp_bg; set 0.0 for pure DP
dt = 0.005; nsteps = 30

dp    = DruckerPrager(C, 0.0, 0.0)
plast = η_vp > 0 ? ParallelModel(dp, LinearViscosity(η_vp)) : dp
c     = SeriesModel(LinearViscosity(η), IncompressibleElasticity(G), plast)
println("Model: η=$η  G=$G  C=$C  η_vp=$η_vp  dt=$dt")

# Background pure-shear strain rate: exx=1, eyy=-1, exy=0 (Voigt 2D: xx,yy,xy)
εᵥ = (1.0, -1.0, 0.0)

# Timestepping — mirrors rc_maxwell_stress_2d call sequence

# Diagnostic: what does RC give when τ0 = 0.5 (plastic fixed point)?
let τ0_diag = (0.5, -0.5, 0.0)
    τ0_rc  = (τ0_diag,)
    others = (; dt = dt, P = 0.0, τ0 = τ0_rc, P0 = (0.0,))
    vars   = (; ε = εᵥ, θ = 0.0)
    x      = initial_guess_x(c, vars, (; τ = 0.0, λ = 0), others)
    τII    = solve(c, x, vars, others; verbose = false)[1]
    τdev   = elastic_stress_history_2D(c, τII, εᵥ, τ0_rc, others)[1]
    τ_II_d = sqrt(0.5 * (τdev[1]^2 + τdev[2]^2 + 2 * τdev[3]^2))
    @printf("Diagnostic: τ0=0.5 → RC τII=%.8f, stored τ_II=%.8f\n", τII, τ_II_d)
end

println("step   t        τ_II      τ_II_analytical(BE)")
τ_rel = η / G
relax = η / (η + G * dt)
τ0ᵥ_ref = Ref((0.0, 0.0, 0.0))

t_vec    = Float64[]
τ_num    = Float64[]
τ_an     = Float64[]

for step in 1:nsteps
    τ0ᵥ = τ0ᵥ_ref[]
    τ0_rc  = (τ0ᵥ,)
    others = (; dt = dt, P = 0.0, τ0 = τ0_rc, P0 = (0.0,))
    vars   = (; ε = εᵥ, θ = 0.0)

    x    = initial_guess_x(c, vars, (; τ = 0.0, λ = 0), others)
    τII  = solve(c, x, vars, others; verbose = false)[1]
    τdev = elastic_stress_history_2D(c, τII, εᵥ, τ0_rc, others)[1]

    τ_II_stored = sqrt(0.5 * (τdev[1]^2 + τdev[2]^2 + 2 * τdev[3]^2))
    t_now  = step * dt
    τ_disc = min(2η * abs(εᵥ[1]) * (1 - relax^step), C)

    @printf("  %2d   %.4f   %.5f   %.5f\n", step, t_now, τ_II_stored, τ_disc)
    push!(t_vec, t_now); push!(τ_num, τ_II_stored); push!(τ_an, τ_disc)

    τ0ᵥ_ref[] = τdev
end

# ── plot ──────────────────────────────────────────────────────────────────────
fig = Figure(size = (650, 380))
ax  = Axis(fig[1, 1]; xlabel = "time", ylabel = "τ_II",
           title = "0D Maxwell VEP  (η=$η, G=$G, C=$C, η_vp=$η_vp, dt=$dt)")
lines!(ax, t_vec, τ_num;  color = :steelblue, linewidth = 2, label = "RC numerical")
lines!(ax, t_vec, τ_an;   color = :black, linestyle = :dash, linewidth = 1.5,
       label = "analytical (discrete BE)")
hlines!(ax, [C]; color = :red, linestyle = :dot, linewidth = 1, label = "yield stress C=$C")
axislegend(ax; position = :rb)
save("rc_0d_plot.pdf", fig)
println("Plot saved to rc_0d_plot.pdf")
