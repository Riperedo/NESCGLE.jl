include("src/API.jl")

ϕ = 0.2
T = 1.0
k = unique(vcat(ChebyshevGrids(0.01,2*π,200)[1], collect(2*π:0.1:15*π)))#
#k = collect(0.0:0.1:15*π)
λ = 1.5
I = Input_SW(ϕ, T, λ, k)
ϕ2 = 0.2
T2 = 0.7
I2 = Input_SW(ϕ2, T2, λ, k)
InstantaneousProcess(I, I2, mute = false)
#α = 0.001
#cooling_rate(α, I, I2, mute = false)
#hysteresis(α, I, I2, mute = false)
