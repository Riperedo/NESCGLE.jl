include("TransportProperties/TP.jl")
include("Steps/Steps.jl")

"""`dynamics(ϕ :: Float64, k :: Array{Float64}, S :: Array{Float64}, dt :: Float64, nT :: Integer, decimations :: Integer)`
Function to evaluate the dynamics of a colloidal system using the Self-Consistent Generalized Langevin Equation theory[1]. This function returns a set of arrays that contains the correlation time, the self and collective Intermediate Scatering Function, the friction memory function, the shear viscosisty relaxation, the diffusion coefficient and the mean squared displacement.
# Arguments
- `ϕ :: Float64`: volume fraction.
- `k :: Array{Float64}`: wave vector.
- `S :: Array{Float64}`: static structure factor.
- `dt :: Float64`: initial time step.
- `nT :: Integer`: Auxiliar number to set the number of intemediate steps.
- `decimations :: Integer`: Number of decades to perform.
# References
[1] Laura Yeomans-Reyna and Magdaleno Medina-Noyola. Overdamped van hove function of colloidal suspensions. Physical Review E - Statistical Physics, Plasmas, Fluids, and Related Interdisciplinary Topics, pages 3382–3399, September 2000.
"""
function dynamics(ϕ :: Float64, k :: Array{Float64}, S :: Array{Float64}, dt :: Float64, 
	nT :: Integer, decimations :: Integer)
	# temporal grid
	N = 2<<nT
	F = zeros(length(k), N)
	Fs = zeros(length(k), N)
	Δζ = zeros(N)
	T = zeros(N)
	D = zeros(N)
	W = zeros(N)

	# preparing output
	index = findfirst(item-> item == maximum(S), S)
	T_save = []
	Δζ_save = []
	Fs_save = []
	F_save = []
	Δη_save = []
	D_save = []
	W_save = []
	#initial conditions
	ONE = ones(length(k))
	kc = 2*π*1.305
	F[:, 1] = S
	Fs[:, 1] = ONE
	Δζ[1] = Δζ_(ϕ, k, S, F[:, 1], Fs[:, 1])
	# first steps
	# free diffusion
	for n in 1:N-1
		T, F, Fs, Δζ, D, W = step_free(ϕ, k, S, F, Fs, Δζ, T, D, W, dt, n)
	end
	# Intermediate steps
	for d in 1:decimations
		# decimation
		T, Δζ, D, W, Fs, F, T_save, Δζ_save, D_save, W_save, Fs_save, F_save, Δη_save = 
		decimation(k, S, T, Δζ, D, W, Fs, F, T_save, Δζ_save, D_save, W_save, Fs_save, F_save, Δη_save, index)
		dt *= 2
		for n in Int(N//2):N
			#T, F, Fs, Δζ, D, W = step_free(L, G, S, F, Fs, Δζ, T, D, W, dt, n-1)
			T, F, Fs, Δζ, D, W = step(ϕ, k, S, F, Fs, Δζ, T, D, W, dt, n-1)
		end
		if Δζ[end] < 1e-30 break end
	end

	# saving final steps
	for n in (Int(N//2) + 1):N
		append!(T_save, T[n]/2)
		append!(Δζ_save, Δζ[n])
		append!(D_save, D[n])
		append!(W_save, W[n])
		append!(Fs_save, Fs[index, n])
		append!(F_save, F[index, n])
		append!(Δη_save, Δη_(k, F[:,n], S))
	end

	return T_save, Fs_save, F_save, Δζ_save, Δη_save, D_save, W_save
end
