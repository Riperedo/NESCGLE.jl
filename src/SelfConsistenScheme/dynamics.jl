using ModeCouplingTheory
include("postprocess.jl")
include("asymptotic.jl")


"""
    dyn_params(Δt::Float64, t_max::Float64, N::Integer, tolerance::Float64, verbose::Bool)

Dynamic parameters configuration for the Self consistent scheme.

More info:
https://github.com/IlianPihlajamaa/ModeCouplingTheory.jl
"""
struct dyn_params 
    Δt::Float64
    t_max::Float64
    N::Integer
    tolerance::Float64
    verbose::Bool
end

function dyn_params()
    Δt = 10^-5
    t_max = 10.0^12
    N = 8
    tolerance = 10^-8
    verbose = false
    return dyn_params(Δt, t_max, N, tolerance, verbose)
end


"""
    SCGLE(ϕ::Float64, k_array::Vector{Float64}, S_array::Vector{Float64}; dp = dyn_params(), k_max = 7.2)

Self consistent scheme of the NESCGLE theory. 

# Arguments
- `ϕ::Float64` Volume fraction
- `k_array::Vector{Float64}` wave vector
- `S_array::Vector{Float64}` static structure_factor

# Keywords
- `dp = dyn_params()` dynamic parameters configuration
- `k_max = 7.2` wave vector to save relaxation times

ϕ::Float64, k_array::Vector{Float64}, S_array::Vector{Float64}; dp = dyn_params(), k_max = 7.2
More info:
https://github.com/IlianPihlajamaa/ModeCouplingTheory.jl
"""
function SCGLE(ϕ::Float64, k_array::Vector{Float64}, S_array::Vector{Float64}; dp = dyn_params(), k_max = 7.2)
    Nk = length(k_array)
    # Structure factor
    s = ones(length(k_array))
    # update the long arrays
    k = [k_array; k_array]
    S = [s; S_array]

    # Initial config
    ∂F0 = zeros(2*Nk); α = 0.0; β = 1.0; γ = @. k^2/S; δ = 0.0

    kernel = SCGLEKernel(ϕ, k, S)
    equation = MemoryEquation(α, β, γ, δ, S, ∂F0, kernel)
    solver = TimeDoublingSolver(Δt=dp.Δt, t_max=dp.t_max, 
        N = dp.N, tolerance=dp.tolerance, verbose=dp.verbose)
    sol = solve(equation, solver)

    # post process
    idx = sum(k_array.<k_max)
    τ = sol.t
    Fs = get_F(sol, 1:length(τ), idx)
    F = get_F(sol, 1:length(τ), Nk + idx)
    t, Δr², Δζ = MSD(k, sol, solver)
    ΔG = get_ΔG(sol, k, S)
    return τ, Fs, F./S[Nk+idx], Δζ, ΔG, Δr²
end
