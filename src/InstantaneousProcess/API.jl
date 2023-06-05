# API application programming interface for call the function developed for the 
# Non-Equilibrium Self-Consistent Generalized Langevin Equation theory.

include("../SCGLE/API.jl")
include("../Utils/API.jl")
include("NETP/API.jl")
include("NETP/long_times_TP.jl")
include("AdaptiveGrid/API.jl")


"""`Asymptotic_NE(I :: Input_SCGLE, F :: Input_SCGLE, u :: Float64; flag= false)`
Compute the asymptotic dynaimics for a given non-equilibrium structural input. This function finds an arrested state if there is one.
# Arguments
- `I :: Input_SCGLE`: Initial Structural input.
- `F :: Input_SCGLE`: Final Structural input.
- `u :: Float64`: Inner clock.
# Keywords
- `flag = false`: flag to print internal procedure.
"""
function Asymptotic_NE(I :: Input_SCGLE, F :: Input_SCGLE, u :: Float64; flag= false)
	ϕ = volume_fraction(F)
	function condition(U)
		if flag println("Computing u = $U") end
		Input = Input_NESCGLE(I, F, U)
		k = wave_vector(Input)
		S = structure_factor(Input)
		system = "Dump"
		#println(minimum)
		if minimum(S) <= 0.0
			system = "Zero"
		else
			iterations, gammas, system = Asymptotic(ϕ, k, S, flag = flag)
		end
		return system != "Fluid"
	end
	return bisection(condition, u, 0.0, 1e-10, flag = false)
end

"""`Asymptotic_NE(I :: Input_SCGLE, F :: Input_SCGLE, u :: Float64; flag= false)`
Compute the asymptotic dynaimics for a given non-equilibrium structural input. This function finds an ergodic state.
# Arguments
- `I :: Input_SCGLE`: Initial Structural input.
- `F :: Input_SCGLE`: Final Structural input.
- `u :: Float64`: Inner clock.
# Keywords
- `flag = false`: flag to print internal procedure.
"""
function Asymptotic_NE_2(I :: Input_SCGLE, F :: Input_SCGLE, u :: Float64; flag= false)
	ϕ = volume_fraction(F)
	function condition(U)
		if flag println("Computing u = $U") end
		Input = Input_NESCGLE(I, F, U)
		k = wave_vector(Input)
		S = structure_factor(Input)
		system = "Dump"
		#println(minimum)
		if minimum(S) <= 0.0
			system = "Zero"
		else
			iterations, gammas, system = Asymptotic(ϕ, k, S, flag = flag)
		end
		return system == "Fluid"
	end
	return bisection(condition, 0.0, u, 1e-10, flag = false)
end

"""`NE_structure_factor(I :: Input_SCGLE, F :: Input_SCGLE, u :: Float64)`
Returns a function to evaluate the Non-Equilibrium Static structure factor[1]
# Arguments
- `I :: Input_SCGLE`: Initial Structural input.
- `F :: Input_SCGLE`: Final Structural input.
- `u :: Float64`: Inner clock.
# References
[1] Pedro Ramírez-González and Magdaleno Medina-Noyola. General nonequilibrium theory of colloid dynamics. Phys. Rev. E, 82:061503, Dec 2010.
"""
function NE_structure_factor(I :: Input_SCGLE, F :: Input_SCGLE, u :: Float64)
	#k = wave_vector(I)
	#@assert k == wave_vector(F) "The grid vector must be the same in both SCGLE inputs."
	Si = structure_factor_function(I)
	Sf = structure_factor_function(F)
	f(x) = Si(x)*exp(-2*u*(x^2)/Sf(x)) + Sf(x)*(1-exp(-2*u*(x^2)/Sf(x)))
	return f
end

"""`Input_NESCGLE(I :: Input_SCGLE, F :: Input_SCGLE, u :: Float64)`
Construct a Non Equilibroum structural input.
# Arguments
- `I :: Input_SCGLE`: Initial Structural input.
- `F :: Input_SCGLE`: Final Structural input.
- `u :: Float64`: Inner clock.
"""
function Input_NESCGLE(I :: Input_SCGLE, F :: Input_SCGLE, u :: Float64; update = false)
	k = wave_vector(F)
	S = NE_structure_factor(I, F, u)
	# addaptative k
	k = update ? adapt(k, S) : k
	return Input_SCGLE(parameters(F), S, k, system(F))
end

"""`InstantaneousProcess(I :: Input_SCGLE, F :: Input_SCGLE; np = 10 :: Int64, decades = 120 :: Int64, mute = false, save_path = "", label = "T" :: String)`
Perform an instantaneous process that starts at the point defined by the paramteres of I and ends in the parameters of F.
# Arguments
- `I :: Input_SCGLE`: Initial Structural input.
- `F :: Input_SCGLE`: Final Structural input.
# Keywords
- `np = 10 :: Int64`: Number of points beteeen decades.
- `decades = 120 :: Int64`: Maximum number of decades.
- `mute = false`: flag to print internal procedure.
- `save_path = ""`: base directory to save data.
- `label = "T" :: String`: label to save the second parameter in the inputs.
- `update_k = false :: Bool`: enable updating the grid in the wave vector.
"""
function InstantaneousProcess(I :: Input_SCGLE, F :: Input_SCGLE; np = 10 :: Int64,
	decades = 120 :: Int64, mute = false, save_path = "", label = "T" :: String, update_k = false :: Bool)
	
	# Save directories
	path = make_directory(save_path*system(F))
	params_i = parameters(I)
	params_f = parameters(F)
	ϕi = params_i[1]
	ϕf = params_f[1]
	phi_save = ""
	if ϕi == ϕf
		phi_save *= "vf"*num2text(ϕi)
	else
		phi_save *= "phii"*num2text(ϕi)*"phif"*num2text(ϕf)
	end
	path = make_directory(path*phi_save)

	if length(params_i) > 1 # for a system different than HS
		Ti = params_i[2]
		Tf = params_f[2]
		T_save = ""
		if Ti == Tf # this imples a compression
			T_save = label*num2text(Ti)
		else
			T_save = label*"i"*num2text(Ti)*label*"f"*num2text(Tf)
		end
		path = make_directory(path*T_save)
	end
		
	# Initial configurations
	if !mute println("Start") end
	# Time parameters
	tw = ExponentialGrid(np, decades, 1e-6)
	t₀ = 0.0
	# Inner clock
	u = 0.0
	# Saving index
	index = 0

	# Initial dynamics
	k = wave_vector(I)
	Si = structure_factor(I)
	τ, f, fs, Δζ, Δη, D, W = SCGLE(I; dt = 1e-10, nT = 5, decimations = 100)
	save_data(path*"/input"*string(index)*".dat", [k Si], header = "1 k\t2 S(k)", flag = false)
	save_data(path*"/fs"*string(index)*".dat", [τ fs f Δζ Δη D W], header = "1 τ\t2 Fs\t3 F\t4 Δζ\t5 Δη\t6 D\t7 W", flag = false)
	# Computing initial inverse mobility
	bI = b⁻¹(τ, Δζ)
	# Some arrays to save data
	index_save = [index]
	t_save = [0.0]
	u_save = [0.0]
	bi_save = [bI]
	tau_save = [τα(τ, fs)]
	deta_save = [dη(τ, Δη)]
	S_save = [maximum(Si)]
	k_save = [k[findfirst(item-> item == maximum(Si), Si)]]
	save_data(path*"/instantaneous_process.dat", [index_save t_save u_save bi_save tau_save deta_save S_save k_save], header = "1 index\t2 t\t3 u\t4 1/b\t5 τα\t6 Δη\t7 S_max\t8 k_max", flag = false)

	if !mute println("index = $index, tw = ", t₀) end

	# Final structure
	Sf = structure_factor(F)
	save_data(path*"/inputF.dat", [k Sf], header = "1 k\t2 S(k)", flag = false)
	# Looking for thermodynamic instabilities
	instability = minimum(Sf) < 0.0
	if !mute if instability println("Warning!\nThemodynamic Instability") end end
	# Computing the max value of u
	u_max = !instability ? Inf : Asymptotic_NE_2(I, F, 1000.0, flag = false)
	if !mute if instability println("u_max = $u_max") end end
	arrest = false
	# Dynamic evolution
	for (index, t) in enumerate(tw)
		# Step
		Δt = t - t₀
		# Inner clock
		u  = minimum([u+(Δt/bI), u_max])
		# Dynamic at time u
		I_NE = Input_NESCGLE(I, F, u, update=update_k)
		S = structure_factor(I_NE)
		k = wave_vector(I_NE)
		τ, f, fs, Δζ, Δη, D, W = SCGLE(I_NE; dt = 1e-10, nT = 5, decimations = 100)
		save_data(path*"/input"*string(index)*".dat", [k S], header = "1 k\t2 S(k)", flag = false)
		save_data(path*"/fs"*string(index)*".dat", [τ fs f Δζ Δη D W], header = "1 τ\t2 Fs\t3 F\t4 Δζ\t5 Δη\t6 D\t7 W", flag = false)
		# Update bI
		bI = b⁻¹(τ, Δζ)
		# Saving data
		append!(index_save, index)
		append!(t_save, t)
		append!(u_save, u)
		append!(bi_save, bI)
		append!(tau_save, τα(τ, fs))
		append!(deta_save, dη(τ, Δη))
		append!(S_save, maximum(S))
		append!(k_save, k[findfirst(item-> item == maximum(S), S)])
		save_data(path*"/instantaneous_process.dat", [index_save t_save u_save bi_save tau_save deta_save S_save k_save], header = "1 index\t2 t\t3 u\t4 1/b\t5 τα\t6 Δη\t7 S_max\t8 k_max", flag = false)
		if !mute println("index = $index, tw = $t") end

		# Update t₀
		t₀ = t

		# Interrupting cases
		if t > 1e15 break end
		if bI > 1e15 || u == u_max
			arrest = true
			break 
		end
		if fs[end] > exp(-1)
			pop!(index_save)
			pop!(t_save)
			pop!(u_save)
			pop!(bi_save)
			pop!(tau_save)
			pop!(deta_save)
			pop!(S_save)
			pop!(k_save)
			arrest = true
			break
		end
	end
	if arrest
		# Computing uₐ
		if !mute println("Arrested state, computing uₐ") end
		uₐ = Asymptotic_NE(I, F, 2*u_save[end], flag = !mute)
		if !mute println("Computing u grid without dynamics.") end
		uⁿ⁺¹ = uₐ - u_save[end]
		uⁿ = uₐ - u_save[end-1]
		tⁿ⁺¹ = t_save[end]
		tⁿ = t_save[end-1]
		μ = (log(uⁿ⁺¹)-log(uⁿ))/(log(tⁿ⁺¹)-log(tⁿ))
		t_long = tw[tw .> tⁿ⁺¹]
		U = -(uⁿ/(tⁿ^μ))*(t_long.^μ) .+ uₐ

		# long times funtions
		bi_save = longtimes(bi_save, t_save, t_long)
		tau_save = longtimes(tau_save, t_save, t_long)
		deta_save = longtimes(deta_save, t_save, t_long)
		# deleting last point
		for i in 1:length(t_long)
			I_NE = Input_NESCGLE(I, F, U[i])
			S = structure_factor(I_NE)
			# Saving data
			append!(index_save, index_save[end]+1)
			append!(t_save, t_long[i])
			append!(u_save, U[i])
			append!(S_save, maximum(S))
			append!(k_save, k[findfirst(item-> item == maximum(S), S)])
		end
		save_data(path*"/instantaneous_process.dat", [index_save t_save u_save bi_save tau_save deta_save S_save k_save], header = "1 index\t2 t\t3 u\t4 1/b\t5 τα\t6 Δη\t7 S_max\t8 k_max", flag = false)
		I_NE = Input_NESCGLE(I, F, uₐ)
		S = structure_factor(I_NE)
		k = wave_vector(I_NE)
		save_data(path*"/inputA.dat", [k S], header = "1 k\t2 S(k)", flag = false)
	end
end
