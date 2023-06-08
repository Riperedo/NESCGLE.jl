include("../../InstantaneousProcess/API.jl")

"""`hysteresis(α :: Float64, I :: Input_SCGLE, F :: Input_SCGLE, k :: Array{Float64}; np = 10 :: Int64, decades = 120 :: Int64, mute = false, save_path = "cooling_rate", label = "T" :: String)`
Perform a cooling down following a cooling rate `α` using equaly spaced `M` steps, and then a warming up back to the initial temperature.
# Arguments
- `α :: Float64`: cooling rate.
- `I :: Input_SCGLE`: Initial Structural input.
- `F :: Input_SCGLE`: Final Structural input.
# Keywords
- `np = 10 :: Int64`: Number of points beteeen decades.
- `decades = 120 :: Int64`: Maximum number of decades.
- `mute = false`: flag to print internal procedure.
- `save_path = "cooling_rate"`: base directory to save data.
- `label = "T" :: String`: label to save the second parameter in the inputs.
- `M = 30 :: Int64`: Number of steps.
- `update_k = false :: Bool`: enable updating the grid in the wave vector.
"""
function hysteresis(α :: Float64, I :: Input_SCGLE, F :: Input_SCGLE; 
	np = 10 :: Int64, decades = 120 :: Int64, mute = false, save_path = "hysteresis", 
	label = "T" :: String, M = 30 :: Int64, ID = "" :: String, update_k = false :: Bool)

	# Saving directory
	path = make_directory(save_path)
	path = make_directory(path*system(F)*ID)
	path = make_directory(path*"rate"*num2text(α))
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

	Ti = params_i[2]
	Tf = params_f[2]
	tp = (Ti-Tf)/α
	
	if !mute println("Current Process:\nHysteresis T_min = $Tf, T_max = $Ti, α = $α") end

	Δt_c = tp/float(M-1)
	t_m = collect(0.0:Δt_c:tp)
	if t_m[end] != tp
		append!(t_m, tp)
	end
	T_m = α*t_m .+ Ti

	tw_aux = ExponentialGrid(np, decades, 1e-6)
	tw_aux = unique(vcat(tw_aux[tw_aux .< Δt_c], [Δt_c]))
	
	T_save = ""
	if Ti == Tf # this imples a compression
		T_save = label*num2text(Ti)
	else
		T_save = label*"i"*num2text(Ti)*label*"f"*num2text(Tf)
	end
	path = make_directory(path*T_save)
	
	# Initial configurations
	if !mute println("Start") end
	if !mute println("Initial Evolution\nStarting cooling down.") end

	# Initial dynamics
	k = wave_vector(I)
	Si = structure_factor(I)
	τ, fs, f, Δζ, Δη, D, W = SCGLE(I; dt = 1e-10, nT = 5, decimations = 100)
	save_data(path*"inputI.dat", [k Si], header = "1 k\t2 S(k)", flag = false)
	save_data(path*"fsI.dat", [τ fs f Δζ Δη D W], header = "1 τ\t2 Fs\t3 F\t4 Δζ\t5 Δη\t6 D\t7 W", flag = false)
	bI = b⁻¹(τ, Δζ)
	# Some arrays to save data
	index_save = [0]
	t_save = [0.0]
	u_save = [0.0]
	bi_save = [bI]
	tau_save = [τα(τ, fs)]
	deta_save = [dη(τ, Δη)]
	S_save = [maximum(Si)]
	k_save = [k[findfirst(item-> item == maximum(Si), Si)]]
	T_save = [Ti]
	save_data(path*"non_instantaneous_process.dat", [index_save t_save u_save bi_save tau_save deta_save S_save k_save T_save], header = "1 index\t2 t\t3 u\t4 1/b\t5 τα\t6 Δη\t7 S_max\t8 k_max\t9 T", flag = false)

	# initial ramp to the final temperature
	for m in 1:length(T_m)
		if !mute println("Step $m of $(2*M)\nCurrent Temperature $(T_m[m])") end
		t₀ = 0.0
		u = 0.0
		# In each step we define the Final structure
		# Final structure
		params_f[2] = T_m[m]
		Sf = structure_factor(params_f, system(F))
		F = Input_SCGLE(params_f, Sf, k, system(F))
		Sf = structure_factor(F)
		save_data(path*"/target"*string(m)*".dat", [k Sf], header = "1 k\t2 S(k)", flag = false)
		# Looking for thermodynamic instabilities
		instability = minimum(Sf) < 0.0
		if !mute if instability println("Warning!\nThemodynamic Instability") end end
		# Computing the max value of u
		u_max = !instability ? Inf : Asymptotic_NE_2(I, F, 1000.0, flag = false)
		if !mute if instability println("u_max = $u_max") end end
		arrest = false
		# Dynamic evolution
		for (index, t) in enumerate(tw_aux)
			# Step
			Δt = t - t₀
			# Inner clock
			u  = minimum([u+(Δt/bI), u_max])
			# Dynamic at time u
			I_NE = Input_NESCGLE(I, F, u, update = update_k)
			S = structure_factor(I_NE)
			k = wave_vector(I_NE)
			τ, f, fs, Δζ, Δη, D, W = SCGLE(I_NE; dt = 1e-10, nT = 5, decimations = 100)
			# Update bI
			bI = b⁻¹(τ, Δζ)
			# Saving data
			if t == tw_aux[end] #|| m == M
				save_data(path*"/input"*string(index_save[end]+1)*".dat", [k S], header = "1 k\t2 S(k)", flag = false)
				save_data(path*"/fs"*string(index_save[end]+1)*".dat", [τ fs f Δζ Δη D W], header = "1 τ\t2 Fs\t3 F\t4 Δζ\t5 Δη\t6 D\t7 W", flag = false)
				append!(index_save, index_save[end]+1)
				append!(t_save, t_m[m]+t)
				#append!(u_save, u)
				append!(u_save, m == 1 ? u : u_save[end] + u)
				append!(bi_save, bI)
				append!(tau_save, τα(τ, fs))
				append!(deta_save, dη(τ, Δη))
				append!(S_save, maximum(S))
				append!(k_save, k[findfirst(item-> item == maximum(S), S)])
				append!(T_save, T_m[m])
				save_data(path*"non_instantaneous_process.dat", [index_save t_save u_save bi_save tau_save deta_save S_save k_save T_save], header = "1 index\t2 t\t3 u\t4 1/b\t5 τα\t6 Δη\t7 S_max\t8 k_max\t9 T", flag = false)
				# Updating I for the next step
				I = Input_NESCGLE(I, F, u, update = update_k)
			end
			if !mute println("index = $index, tw = $(t_m[m]+t)") end

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
				pop!(T_save)
				arrest = true
				break
			end
		end
	end
	
	if !mute println("Final temperature reached!\nStarting heating process") end

	# final ramp to the initial temperature
	reverse!(T_m)
	append!(T_m, Ti)
	for m in 2:length(T_m)
		if !mute println("Step $(M+m-1) of $(2*M)\nCurrent Temperature $(T_m[m])") end
		t₀ = 0.0
		u = 0.0
		# In each step we define the Final structure
		# Final structure
		params_f[2] = T_m[m]
		Sf = structure_factor(params_f, system(F))
		F = Input_SCGLE(params_f, Sf, k, system(F))
		Sf = structure_factor(F)
		save_data(path*"/target"*string(M+m-1)*".dat", [k Sf], header = "1 k\t2 S(k)", flag = false)
		# Looking for thermodynamic instabilities
		instability = minimum(Sf) < 0.0
		if !mute if instability println("Warning!\nThemodynamic Instability") end end
		# Computing the max value of u
		u_max = !instability ? Inf : Asymptotic_NE_2(I, F, 1000.0, flag = false)
		if !mute if instability println("u_max = $u_max") end end
		arrest = false
		# Dynamic evolution
		for (index, t) in enumerate(tw_aux)
			# Step
			Δt = t - t₀
			# Inner clock
			u  = minimum([u+(Δt/bI), u_max])
			# Dynamic at time u
			I_NE = Input_NESCGLE(I, F, u, update = update_k)
			S = structure_factor(I_NE)
			k = wave_vector(I_NE)
			τ, f, fs, Δζ, Δη, D, W = SCGLE(I_NE; dt = 1e-10, nT = 5, decimations = 100)
			# Update bI
			bI = b⁻¹(τ, Δζ)
			# Saving data
			if t == tw_aux[end] #|| m == 2# || m == M
				save_data(path*"/input"*string(index_save[end]+1)*".dat", [k S], header = "1 k\t2 S(k)", flag = false)
				save_data(path*"/fs"*string(index_save[end]+1)*".dat", [τ fs f Δζ Δη D W], header = "1 τ\t2 Fs\t3 F\t4 Δζ\t5 Δη\t6 D\t7 W", flag = false)
				append!(index_save, index_save[end]+1)
				append!(t_save, t_m[m-1]+t)
				#append!(u_save, u)
				append!(u_save, m == 1 ? u : u_save[end] + u)
				append!(bi_save, bI)
				append!(tau_save, τα(τ, fs))
				append!(deta_save, dη(τ, Δη))
				append!(S_save, maximum(S))
				append!(k_save, k[findfirst(item-> item == maximum(S), S)])
				append!(T_save, T_m[m])
				save_data(path*"non_instantaneous_process.dat", [index_save t_save u_save bi_save tau_save deta_save S_save k_save T_save], header = "1 index\t2 t\t3 u\t4 1/b\t5 τα\t6 Δη\t7 S_max\t8 k_max\t9 T", flag = false)
				# Updating I for the next step
				I = Input_NESCGLE(I, F, u, update = update_k)
			end
			if !mute println("index = $index, tw = $(t_m[m-1]+t)") end

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
				pop!(T_save)
				arrest = true
				break
			end
		end
	end
	
end
