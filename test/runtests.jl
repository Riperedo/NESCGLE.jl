using NESCGLE
using Test

@testset "NESCGLE.jl" begin
	# Write your tests here.
	ϕ = 0.61
	T = 1.0
	#k = unique(vcat(ChebyshevGrids(0.01,2*π,200)[1], collect(2*π:0.1:15*π)))#
	k = collect(0.0:0.1:15*π)
	#λ = 1.5
	I = Input_WCA(ϕ, T, k)
	#S = structure_factor(I)
	@test typeof(I) == Input_SCGLE
	#@test typeof(S) == typeof(k)
	#save_data("test.dat", [k S])
	τ, Fs, F, Δζ, Δη, D, W = SCGLE(I)
	@test length(τ) == length(W)
	ϕ2 = 0.61
	T2 = 0.0001
	I2 = Input_WCA(ϕ2, T2, k)
	InstantaneousProcess(I, I2, mute = false)
	#α = 10.0
	#cooling_rate(α, I, I2, mute = false)
	#hysteresis(α, I, I2, mute = false)
end
