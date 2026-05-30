using RheologyCalculator
include(joinpath(pkgdir(RheologyCalculator), "rheologies", "RheologyDefinitions.jl"))
include(joinpath(pkgdir(RheologyCalculator), "examples", "tensor_helpers.jl"))

function make_rc_visc(eta0, n, eps0, eta_min, eta_max)
    visc = if n ≈ 1.0
        LinearViscosity(eta0)
    else
        eta_RC = 2.0^(n-1) * eta0^n * eps0^(n-1)
        PowerLawViscosity(eta_RC, n)
    end
    v_upper = eta_max < 1e29 ? SeriesModel(visc, LinearViscosity(eta_max)) : visc
    v_full  = eta_min > 0.0  ? ParallelModel(v_upper, LinearViscosity(eta_min)) : v_upper
    return v_full
end

# Test 1: n=1 model type
c1 = SeriesModel(make_rc_visc(1.0, 1.0, 1.0, 0.0, 1e30), IncompressibleElasticity(1e30))
println("n=1 viscosity model: ", c1.leafs)

# Test 2: n=3 power-law stress vs analytical
# η₀=1, ε₀=1, n=3 → η_RC=1; τ_II = (2*ε_II)^(1/3)
c3 = SeriesModel(make_rc_visc(1.0, 3.0, 1.0, 0.0, 1e30), IncompressibleElasticity(1e30))
others = (; dt=1.0, P=0.0, τ0=((0.0,0.0,0.0),), P0=(0.0,))
args0  = (; τ=0.0, λ=0)
vars   = (; ε=(0.1,-0.1,0.0), θ=0.0)
x0 = initial_guess_x(c3, vars, args0, others)
tauII = solve(c3, x0, vars, others; verbose=false)[1]
eps_II = sqrt(0.5*(0.1^2+0.1^2))
tau_analytic = (2*eps_II)^(1/3)
println("n=3: tauII RC=", tauII, "  analytical=", tau_analytic, "  relerr=", abs(tauII-tau_analytic)/tau_analytic)
println("Test 2 passed: ", abs(tauII-tau_analytic)/tau_analytic < 1e-6)

# Test 3: n=3 with ε₀=0.5 (should still match Ostwald-de Waele)
# η₀=1, ε₀=0.5, n=3 → η_RC = 2^2 * 1 * 0.5^2 = 1.0
# τ_II = (2*η_RC*ε_II)^(1/3) = (2*ε_II)^(1/3) same as test 2 since η_RC=1.0
# Ostwald: τ = 2*η₀*(ε_II/ε₀)^(1/3-1)*ε_II = 2*1*(0.1/0.5)^(-2/3)*0.1 = 0.2*(0.2)^(-2/3)
c3b = SeriesModel(make_rc_visc(1.0, 3.0, 0.5, 0.0, 1e30), IncompressibleElasticity(1e30))
x0b = initial_guess_x(c3b, vars, args0, others)
tauII_b = solve(c3b, x0b, vars, others; verbose=false)[1]
tau_analytic_b = 2*1*(eps_II/0.5)^(1/3-1)*eps_II
println("n=3 ε₀=0.5: tauII RC=", tauII_b, "  analytical=", tau_analytic_b, "  relerr=", abs(tauII_b-tau_analytic_b)/tau_analytic_b)
println("Test 3 passed: ", abs(tauII_b-tau_analytic_b)/tau_analytic_b < 1e-6)

# Test 4: cutoffs don't error
vc = make_rc_visc(1.0, 3.0, 1.0, 0.1, 10.0)
println("n=3 with cutoffs: ", vc)
