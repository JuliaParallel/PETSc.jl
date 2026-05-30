using ForwardDiff, StaticArrays, RheologyCalculator
include(joinpath(pkgdir(RheologyCalculator), "rheologies", "RheologyDefinitions.jl"))
include(joinpath(pkgdir(RheologyCalculator), "examples", "tensor_helpers.jl"))

function rc_maxwell_stress_2d(c, εᵥ::NTuple{3,<:Real}, τ0ᵥ::NTuple{3,<:Real}, dt::Real)
    τ0_rc  = (τ0ᵥ,)
    others = (; dt = dt, P = 0.0, τ0 = τ0_rc, P0 = (0.0,))
    vars   = (; ε = εᵥ, θ = 0e0)
    x      = initial_guess_x(c, vars, (; τ = 0e0, λ = 0), others)
    τII    = solve(c, x, vars, others; verbose = false)[1]
    return elastic_stress_history_2D(c, τII, εᵥ, τ0_rc, others)[1]
end

function compute_stress_tensor_2d(ε::SVector{3,T}, c, τ0ᵥ, dt, index) where T
    τdev = rc_maxwell_stress_2d(c, (ε[1], ε[2], ε[3]), τ0ᵥ, dt)
    return τdev[index]
end

dt = 0.1
d  = ForwardDiff.derivative

println("=== Linear Maxwell (n=1) ===")
let c = SeriesModel(LinearViscosity(1.0), IncompressibleElasticity(1.0))
    eta_eff = 1.0 * 1.0 * dt / (1.0 + 1.0 * dt)
    println("  Analytical η_eff=$eta_eff → ∂τxx/∂εxx should be $(2*eta_eff)")
    for (label, εᵥ, τ0) in [
        ("ε=0     ", (0.0,  0.0,  0.0), (0.0, 0.0, 0.0)),
        ("ε=0.1   ", (0.1, -0.1,  0.0), (0.0, 0.0, 0.0)),
        ("ε=0.1+τ0", (0.1, -0.1,  0.0), (0.1,-0.1, 0.0)),
    ]
        ∂ = d(e -> compute_stress_tensor_2d(SA[e, εᵥ[2], εᵥ[3]], c, τ0, dt, 1), εᵥ[1])
        println("  $label: ∂τxx/∂εxx = $∂  $(isnan(∂) ? "NaN!" : "ok")")
    end
end

println("\n=== Power-law purely viscous (n=3, G=Inf) ===")
let c = SeriesModel(PowerLawViscosity(1.0, 3), IncompressibleElasticity(1e30))
    # Analytical: τ = (2η·ε)^(1/n), so ∂τxx/∂εxx = (1/n)*(2η)^(1/n)*εxx^(1/n-1)
    εxx = 0.1
    n = 3; η = 1.0
    tau_analytical = (2η * εxx)^(1/n)
    dtau_deps_analytical = (1/n) * (2η)^(1/n) * εxx^(1/n - 1)
    println("  τxx analytical = $tau_analytical")
    println("  ∂τxx/∂εxx analytical = $dtau_deps_analytical")
    for (label, εᵥ, τ0) in [
        ("ε=0     ", (0.0,  0.0,  0.0), (0.0, 0.0, 0.0)),
        ("ε=0.1   ", (0.1, -0.1,  0.0), (0.0, 0.0, 0.0)),
    ]
        δ = 1e-20
        ε_II_sq = 0.5 * (εᵥ[1]^2 + εᵥ[2]^2 + 2*εᵥ[3]^2)
        εᵥ_reg = ε_II_sq < δ^2 ? (δ, -δ, 0.0) : εᵥ
        ∂ = d(e -> compute_stress_tensor_2d(SA[e, εᵥ_reg[2], εᵥ_reg[3]], c, τ0, dt, 1), εᵥ_reg[1])
        println("  $label: ∂τxx/∂εxx = $∂  $(isnan(∂) ? "NaN!" : "ok")")
    end
end

println("\n=== Power-law Maxwell (n=3, G=1) ===")
let c = SeriesModel(PowerLawViscosity(1.0, 3), IncompressibleElasticity(1.0))
    for (label, εᵥ, τ0) in [
        ("ε=0.1   ", (0.1, -0.1,  0.0), (0.0, 0.0, 0.0)),
        ("ε=0.1+τ0", (0.1, -0.1,  0.0), (0.1,-0.1, 0.0)),
    ]
        ∂ = d(e -> compute_stress_tensor_2d(SA[e, εᵥ[2], εᵥ[3]], c, τ0, dt, 1), εᵥ[1])
        println("  $label: ∂τxx/∂εxx = $∂  $(isnan(∂) ? "NaN!" : "ok")")
    end
end
